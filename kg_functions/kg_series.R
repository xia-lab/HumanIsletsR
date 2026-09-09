# ==============================================================================
# kg_functions/kg_series.R  --  kgSeries: DONOR-LEVEL VALUES for plotting.
#
# ⚠ THIS COMPUTES NOTHING. It assembles values that already exist -- the feature's row from
# HI_omics_v2.sqlite and the phenotype's column from metadata_sum_norm -- so it is safe to
# fire on the fast (auto) path alongside retrieval, before any analysis runs.
#
# WHY IT EXISTS. The plot spec routes the STATS views (volcano / rank_bar / manhattan /
# locus / network / matrix) to neo4j and the COMPUTED views (forest / enrichment_bar) to an
# engine result -- but the DONOR-LEVEL views have no source at all:
#     scatter   feature value vs a continuous phenotype
#     box       feature value grouped by a categorical phenotype
#     trajectory  the same across an ORDERED set of levels
#     before_after  raw vs covariate-adjusted values, side by side
# One row shape serves all four.
#
# ------------------------------------------------------------------------------
# PARAMETERS
# ------------------------------------------------------------------------------
# features    NAMED features only -- "INS" | "INS,GCG" | "INS|GCG". ANY entity the resolver
#             knows (gene | protein | metabolite | contaminant | reaction).
#             ⚠ "all" is REFUSED (RES-NO-FEATURE). A plot is per-feature: the caller decides
#             WHICH features are worth drawing (the ranking step does that), and this
#             function will not silently emit a whole layer.
#             ⚠ NO CAP is applied here. How many panels to draw is a DISPLAY decision that
#             belongs to the caller, not a threshold invented in the data layer.
#
# phenotypes  one or more metadata columns, or a PARENT concept that expands to its children
#             (same expansion as kgAssoc). EMPTY = feature values only, with no phenotype
#             columns -- which is what a distribution/box-by-nothing view needs.
#             ⚠ "all" is not supported: a donor x feature x every-phenotype product is not a
#             plot, it is a table dump.
#
# layers      "all" (the feature's OWN layers, resolver-driven) or a comma list of layer
#             keys. Entity-driven exactly as in kgAssoc -- the layer map is SHARED, not
#             copied, so a series can never disagree with the analysis about which table a
#             feature came from.
#
# subset      "all" | "col=value" -- the same donor subset kgAssoc applies, so a plot can be
#             drawn over exactly the donors an analysis was run on.
#
# covariates  "none" (DEFAULT) -> `value_adjusted` is empty.
#             a comma list      -> `value_adjusted` is the feature value RESIDUALISED on
#                                  those covariates, which is the "after" panel of a
#                                  before/after view.
#             ⚠ This is a DONOR-LEVEL adjustment for display. The before/after of an EFFECT
#             (coefficient +- CI, per F26) comes from the two A(1) runs, not from here.
#
# use_raw     "false" (DEFAULT) -> metadata_sum_norm | "true" -> metadata_sum_raw.
#
# minN        minimum donors for a series to be emitted at all (default 10). A series below
#             it is dropped rather than drawn, because a scatter of 4 points invites a
#             conclusion the data cannot carry.
#
# ------------------------------------------------------------------------------
# OUTPUT -- kg_series.csv, one row per donor x feature x layer x phenotype
# ------------------------------------------------------------------------------
#   donor_id  feature  symbol  feature_type  layer
#   value  value_adjusted
#   phenotype  pheno_value  pheno_type      ("cont" | "disc")
#   n_donors  subset
#
# `pheno_value` carries the level itself for a categorical phenotype, so a separate
# group column would only duplicate it; `trajectory` gets its level ORDER from the parse
# (F29 fixes it there), not from the data.
#
# Returns "RES-OK:<file>" | "RES-NO" | "RES-NO-DB" | "RES-NO-META" | "RES-NO-FEATURE"
#       | "RES-NO-RESOLVE" | "RES-NO-RESOLVER" | "RES-NO-PHENO".
#
# Sourced AFTER kg_common.R and kg_assoc.R: it reuses the resolver + layer vocabulary from
# those (.kgSplitFeatures / .kgLayerKey / .kgAllLayerKeys) rather than restating them.
# ==============================================================================

kgSeries <- function(features = "", phenotypes = "", layers = "all",
                     subset = "all", covariates = "none", use_raw = "false",
                     minN = "10", mode = "tool"){
  library(RSQLite); library(DBI)
  .kgSetPaths()
  minN <- suppressWarnings(as.integer(minN)); if(is.na(minN)) minN <- 10L
  db.path <- paste0(sqlite.path, "HI_omics_v2.sqlite"); if(!file.exists(db.path)) return("RES-NO-DB")

  feat_raw <- trimws(features)
  if(!nzchar(feat_raw)) return("RES-NO-FEATURE")
  # a plot is per-feature: refuse a whole-layer dump rather than emit one
  if(identical(tolower(feat_raw), "all")) return("RES-NO-FEATURE")

  m <- .kgLoadMetaFrame(use_raw = tolower(trimws(use_raw)) %in% c("true", "1", "yes"))
  if(is.null(m)) return("RES-NO-META")
  meta_names <- names(m); rec <- as.character(m[["record_id"]])

  # ---- phenotypes (OPTIONAL; "all" deliberately unsupported) ----
  ptoks  <- trimws(strsplit(phenotypes, "[,;]")[[1]]); ptoks <- ptoks[nzchar(ptoks)]
  if(length(ptoks) == 1 && tolower(ptoks) == "all") return("RES-NO-PHENO")
  phenos <- if(length(ptoks) == 0) character(0) else .kgExpandPhenotypes(ptoks, meta_names)
  phenos <- phenos[phenos %in% meta_names]
  if(length(ptoks) > 0 && length(phenos) == 0) return("RES-NO-PHENO")

  # ---- donor subset (identical to kgAssoc, so a plot matches its analysis) ----
  keep_ids <- rec; subset_note <- "all donors"
  # -- MUST STAY IDENTICAL TO kgAssoc (the note above), so this moved to the shared
  #    `.kgDonorSubset` in the same change. Left on the old equality-only parser, a plot of
  #    "donors over 60" would have been drawn from EVERY donor while the analysis beside it was
  #    correctly restricted -- the exact mismatch that note exists to prevent.
  if(!identical(subset, "all") && nzchar(subset)){
    msub <- .kgDonorSubset(m, subset)
    if(is.null(msub)) return("RES-NO-SUBSET")
    keep_ids <- as.character(msub[["record_id"]]); subset_note <- subset
  }
  keep_ids <- keep_ids[!is.na(keep_ids)]

  # ---- covariates -> the residualisation basis for `value_adjusted` ----
  covuniv <- if(covariates %in% c("none", "NA", "")) character(0)
             else trimws(strsplit(covariates, "[,;]")[[1]])
  covuniv <- unique(covuniv[nzchar(covuniv) & covuniv %in% meta_names])
  covbase <- if(length(covuniv)){ cb <- .kgCoerce(m, covuniv); rownames(cb) <- rec; cb } else NULL

  # ---- phenotype values (coerced once) ----
  phframe <- if(length(phenos)){ pf <- .kgCoerce(m, phenos); rownames(pf) <- rec; pf } else NULL
  phtypes <- if(length(phenos)) vapply(phenos, function(p) .kgPhenoType(m, p), character(1)) else character(0)

  # ---- resolve the named features (ANY entity) -- same path as kgAssoc ----
  res <- .kgLoadFeatureResolver(); if(is.null(res)) return("RES-NO-RESOLVER")
  ftoks <- .kgSplitFeatures(feat_raw, res); ftoks <- ftoks[nzchar(ftoks)]
  if(length(ftoks) == 0) return("RES-NO-FEATURE")
  lay_filter <- if(identical(layers, "all")) NULL else tolower(trimws(strsplit(layers, "[,;]")[[1]]))

  rmap <- .kgResolveFeatures(res, ftoks, "auto")
  feats <- list()
  for(tk in ftoks){
    r <- rmap[[.kgNorm(tk)]]
    if(is.null(r) || length(r$layers) == 0 || is.na(r$type)) next
    for(la in r$layers){
      lk <- .kgLayerKey(la$display)
      if(!(lk %in% .kgAllLayerKeys)) next
      if(!is.null(lay_filter) && !(lk %in% lay_filter || tolower(la$display) %in% lay_filter)) next
      feats[[length(feats)+1]] <- list(token = tk, etype = r$type, layer = la)
    }
  }
  if(length(feats) == 0) return("RES-NO-RESOLVE")

  con <- dbConnect(SQLite(), db.path); on.exit(dbDisconnect(con), add = TRUE)
  rows <- list()

  # group by layer table so each layer is read ONCE (one IN-query over its features)
  key <- vapply(feats, function(f) paste(f$layer$display, f$layer$table, f$layer$read_col, sep = "|"),
                character(1))
  for(k in unique(key)){
    idx <- which(key == k); la <- feats[[idx[1]]]$layer
    rvals <- vapply(feats[idx], function(f) as.character(f$layer$read_val), character(1))
    FR <- .kgFeatureRows(con, la$table, la$read_col, rvals); if(is.null(FR)) next
    dcommon <- intersect(colnames(FR$mat), keep_ids); if(length(dcommon) < minN) next
    rpos <- match(rvals, FR$id)

    for(j in seq_along(idx)){
      f <- feats[[idx[j]]]; ri <- rpos[j]; if(is.na(ri)) next
      y <- FR$mat[ri, dcommon]; names(y) <- dcommon
      sym <- FR$symbol[ri]

      # residualised values -- computed ONCE per feature, then looked up by donor
      yadj <- if(!is.null(covbase)){
                cb <- covbase[intersect(dcommon, rownames(covbase)), , drop = FALSE]
                .kgResidualize(y, cb)
              } else NULL

      if(length(phenos) == 0){
        don <- dcommon[is.finite(y[dcommon])]
        if(length(don) < minN) next
        rows[[length(rows)+1]] <- data.frame(
          donor_id = don, feature = f$token, symbol = sym, feature_type = f$etype,
          layer = la$display,
          value = as.numeric(y[don]),
          value_adjusted = if(is.null(yadj)) NA_real_ else as.numeric(yadj[don]),
          phenotype = NA_character_, pheno_value = NA_character_, pheno_type = NA_character_,
          n_donors = length(don), subset = subset_note, stringsAsFactors = FALSE)
        next
      }

      for(pi in seq_along(phenos)){
        ph <- phenos[pi]
        pv <- phframe[dcommon, ph]; names(pv) <- dcommon
        ok  <- is.finite(y[dcommon]) & !is.na(pv) & !(as.character(pv) %in% c("", "NA"))
        don <- dcommon[ok]
        if(length(don) < minN) next
        rows[[length(rows)+1]] <- data.frame(
          donor_id = don, feature = f$token, symbol = sym, feature_type = f$etype,
          layer = la$display,
          value = as.numeric(y[don]),
          value_adjusted = if(is.null(yadj)) NA_real_ else as.numeric(yadj[don]),
          phenotype = ph, pheno_value = as.character(pv[don]), pheno_type = unname(phtypes[pi]),
          n_donors = length(don), subset = subset_note, stringsAsFactors = FALSE)
      }
    }
  }

  if(length(rows) == 0) return("RES-NO")
  out <- do.call(rbind, rows)

  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  safe <- function(s) gsub("[^A-Za-z0-9]", "_", substr(s, 1, 30))
  # ⚠ Return the UNIQUE (timestamped) filename in the "RES-OK:<file>" form the
  # frontend read-back expects (humanislets_utils.R:2000). A fixed name would let
  # two rapid selections race — the second write landing between the first's write
  # and its download — so the caller must read the exact file it was handed.
  fname <- paste0("kg_series_", safe(features), "_", safe(phenotypes), "_", ts, ".csv")
  utils::write.csv(out, fname, row.names = FALSE)
  utils::write.csv(out, "kg_series.csv", row.names = FALSE)   # stable copy, harmless
  paste0("RES-OK:", fname)
}
