# ==============================================================================
# kg_functions/kg_assoc.R  --  kgAssoc: the unified feature <-> phenotype ASSOCIATION
# entry point (op#1). ONE interface over the full input grid; TWO internal engines chosen
# by the FEATURE axis; orthogonal modifiers (covariates / subset / layers / type) applied
# the same way in both.
#
#   features    named set ("INS" | "INS,GCG,PDX1")  -> TARGETED engine (lm)
#               "all"                                -> SCREEN engine (limma, all features)
#               ANY entity: gene | metabolite | contaminant | reaction(flux) -- one engine,
#               same model; the resolver supplies each one's table + read key.
#   phenotypes  a leaf metadata column | a PARENT concept (e.g. "insulin secretion" /
#               "gsis") that expands to its child columns | a comma list | "all"
#   layers      "all" | comma list of layer keys. LAYERS ARE ENTITY-DRIVEN -- a gene covers
#               only the gene layers, a metabolite only the metabolite layers, etc.:
#                 gene        rnaseq | protein | nanostring | pbrna_alpha | pbrna_beta
#                 metabolite  metabolite_hg | metabolite_lg | metabolite_ratio
#                 contaminant contaminant
#                 flux        flux_hg | flux_lg | flux_ratio   (supported, NOT default)
#               "all" on a TARGETED query = that feature's own layers (resolver-driven);
#               "all" on a SCREEN = the GENE layers -- name a layer to screen a metabolite /
#               contaminant / flux. OUT OF SCOPE: pbrna_delta (n=3), pbrna_pp (n=7).
#   covariates  "default" (age,sex,BMI,culturetime; primary dropped) | "none" | list
#               ⚠ purity removed 2026-09-04 to match the precompute default (dropped there
#               2026-08); the two must name the same set -- see kg_association.R's header.
#   subset      "all" | "col=value" (donor subset before fitting)
#   contrast    categorical only: "A-B" (A vs B) | "anova" | "" (auto: 2-level = A-B,
#               >=3 levels = omnibus F)
#   fdr_family  "query" (BH across exactly the tests asked -> P_family; NEVER genome-wide)
#               | "none" | "genomewide" (screen only: copies limma's genome-wide adj-p)
#   method      "auto" (DEFAULT) | "lm" | "limma" | "kendall"
#               auto = KENDALL when the layer is `contaminant`, otherwise the normal split
#               (targeted -> lm, screen -> limma). This mirrors the DEPLOYED trigger, which
#               is layer-based rather than a flag: DonorRegression switches to Kendall on
#               omicsType == "proc_contaminants" (humanislets_statistics.R:1878), because
#               small n and non-normal abundances make limma's moderated-t inappropriate
#               there. Naming a method overrides the choice.
#               ⚠ Kendall's effect is TAU, so EffectType is "tau" (or "kruskalH" for a >=3
#               level omnibus, which has no per-feature effect size) -- never "log2FC".
#
# TARGETED = lm(value ~ phenotype + covariates): cont -> slope; disc 2-level/named -> mean
# difference; disc >=3 -> covariate-adjusted omnibus ANOVA F. Reads only the queried
# feature rows (one IN-query per layer). This is the fast path for one/few features.
# SCREEN  = limma (trend+robust eBayes) over ALL features, one fit per phenotype, reading
# each layer once. This is the genome-wide differential (feeds enrichment etc.).
#
# Retrieval note: kg_functions are on-demand COMPUTE. The "is this in the precompute ->
# retrieve" decision lives in the caller (it reads a different store, HI_precomputed);
# kgAssoc always COMPUTES the requested spec. all-features x all-phenotypes recompute is
# refused (RES-USE-PRECOMPUTE) unless narrowed, because that IS the precompute.
#
# Writes kg_assoc.csv (one row per feature x phenotype x layer). Returns "RES-OK;<n_rows>"
#   | "RES-USE-PRECOMPUTE" | "RES-NO" | "RES-NO-DB" | "RES-NO-META" | "RES-NO-PHENO"
#   | "RES-NO-FEATURE" | "RES-NO-RESOLVE" | "RES-NO-RESOLVER".
# ==============================================================================

.kgAssocDefaultCovs <- c("donorage", "donorsex", "bodymassindex", "predistributionculturetime")   # purity dropped 2026-09-04 to match the precompute default (4 covariates)

# ENTITY SCOPE: gene, metabolite, contaminant and reaction(flux) all run through the SAME
# engine -- the model is identical (lm on the feature's donor vector); only the source table
# differs, and the resolver already carries it. LAYERS ARE ENTITY-DRIVEN: a gene only ever
# covers the gene layers, a metabolite only the metabolite layers, a contaminant/reaction its
# own -- enforced because the resolver returns a feature's OWN layers, filtered to the
# supported keys below. FLUX is SUPPORTED but NOT a default layer (user, 2026-07-16): query a
# reaction token or name layers="flux_hg" to reach it; it never enters a gene/screen default.
# OUT OF SCOPE (user, 2026-07-16): pbrna_delta (n=3) and pbrna_pp (n=7) -- below minN, so
# always empty.
.kgScreenLayer <- function(name){
  switch(tolower(name),
    "rnaseq"           = c("proc_rnaseq",                  "accession"),
    "protein"          = c("proc_prot_combat",             "id"),
    "nanostring"       = c("proc_nanostring_merge",        "gene_id"),
    "pbrna_alpha"      = c("proc_pbrna_Alpha",             "symbol"),
    "pbrna_beta"       = c("proc_pbrna_Beta",              "symbol"),
    "metabolite_hg"    = c("proc_metabolite_combat_HG",    "compound"),
    "metabolite_lg"    = c("proc_metabolite_combat_LG",    "compound"),
    "metabolite_ratio" = c("proc_metabolite_combat_ratio", "compound"),
    "contaminant"      = c("proc_contaminants",            "Compound"),
    "flux_hg"          = c("proc_flux_HG",                 "rxn"),
    "flux_lg"          = c("proc_flux_LG",                 "rxn"),
    "flux_ratio"       = c("proc_flux_ratio",              "rxn"),
    NULL)
}

# SUPPORTED layer keys, grouped by entity. Note the split between SUPPORTED and DEFAULT:
# a SCREEN (features="all") defaults to the GENE layers only (what "run the differential
# across the omics" means, and what feeds enrichment); metabolite / contaminant / flux are
# SUPPORTED but never default -- reached by querying that entity's token (layers are
# entity-driven) or by naming the layer (e.g. layers="flux_hg").
.kgGeneLayerKeys        <- c("rnaseq","protein","nanostring","pbrna_alpha","pbrna_beta")
.kgMetaboliteLayerKeys  <- c("metabolite_hg","metabolite_lg","metabolite_ratio")
.kgContaminantLayerKeys <- c("contaminant")
.kgFluxLayerKeys        <- c("flux_hg","flux_lg","flux_ratio")
.kgAllLayerKeys <- c(.kgGeneLayerKeys, .kgMetaboliteLayerKeys,
                     .kgContaminantLayerKeys, .kgFluxLayerKeys)

# short layer key from a resolver display name, so `layers` takes short codes in BOTH engines
# (e.g. "Bulk gene expression (RNA-seq)" -> "rnaseq"; "Metabolite (LG)" -> "metabolite_lg").
.kgLayerKey <- function(display){
  d <- tolower(display)
  if(grepl("contaminant", d)) return("contaminant")
  if(grepl("metabolic flux", d))
    return(if(grepl("ratio", d)) "flux_ratio" else if(grepl("\\(lg\\)", d)) "flux_lg" else "flux_hg")
  if(grepl("metabolite", d))
    return(if(grepl("ratio", d)) "metabolite_ratio" else if(grepl("\\(lg\\)", d)) "metabolite_lg" else "metabolite_hg")
  if(grepl("pseudobulk", d)){
    if(grepl("alpha", d)) return("pbrna_alpha")
    if(grepl("beta",  d)) return("pbrna_beta")
    if(grepl("delta", d)) return("pbrna_delta")
    if(grepl("pp",    d)) return("pbrna_pp")
    return("pbrna")
  }
  if(grepl("protein", d))        return("protein")
  if(grepl("nanostring", d))     return("nanostring")
  if(grepl("rna-seq|rnaseq", d)) return("rnaseq")
  gsub("[^a-z0-9]+", "_", d)
}

# Split a feature list. Metabolite/contaminant NAMES legitimately contain commas
# (58 tokens, e.g. "1,2,3,4,6,7,8-HPCDD", "2,4,7,9-tetramethyl-5-decyne-4,7-diol"), while NO
# token contains ";" or "|". So: ";"/"|" are ALWAYS separators; "," is a separator only when
# the whole string is not itself one resolvable token. Keeps "INS,GCG,PDX1" working while
# never shredding a comma-bearing compound name.
.kgSplitFeatures <- function(s, res){
  if(grepl("[;|]", s)) return(trimws(strsplit(s, "[;|]")[[1]]))
  if(nrow(res[res$token == .kgNorm(s), , drop = FALSE]) > 0) return(s)
  trimws(strsplit(s, ",", fixed = TRUE)[[1]])
}

.kgLayerEntity <- function(key){
  if(key %in% c("metabolite_hg","metabolite_lg","metabolite_ratio")) "metabolite"
  else if(identical(key, "contaminant")) "contaminant"
  else if(key %in% c("flux_hg","flux_lg","flux_ratio")) "reaction"
  else "gene"
}

# Entity types whose donor n is too small to spend the default covariate DoF (metabolite ~37,
# contaminants ~31) -> the DEFAULT covariate set is dropped for them (an explicit user
# covariate list is still honored). Flux is NOT here: it has 371 donors, a full-size layer.
.kgSmallNTypes <- c("metabolite", "contaminant")

# lm -> limma switch for a NAMED feature set (user-set 2026-07-16). Past this many named
# features the query is a screen in all but name, so we fit the FULL layer and extract the
# named rows: moderated stats (prior from ALL features, not the named subset) instead of k
# ordinary-t lm's. NOT a statistical cliff -- nothing changes at exactly 500. Measured on
# rnaseq (1 thread): a full-layer limma is a FLAT ~3.3s at any k, while k lm's cost ~0.0036s
# each, so the two cross at k ~= 870. 500 sits below that crossover on purpose: the switch
# costs ~1.4s there but buys the better statistic, and it stays far above any real named set
# (1-50), so interactive queries keep the ~16x faster lm path.
.kgLmMaxFeatures <- 500L

# resolve the fitting method for ONE layer.
# ⚠ `auto` reproduces the DEPLOYED trigger, which is LAYER-based, not a flag: DonorRegression
# switches to Kendall on `omicsType == "proc_contaminants"`
# (humanislets_statistics.R:1878). So contaminant x phenotype gets Kendall and every other
# layer keeps the existing targeted-lm / screen-limma split -- no other frame changes.
# Scope is exactly that pairing: a collaborator decision that is owned, never derived from
# the data (deriving it as "small-n or tie-heavy" would silently capture metabolomics and
# methylation, which the decision never covered).
.kgAssocMethod <- function(method, layer_key, is_screen){
  mm <- tolower(trimws(if(is.null(method) || !nzchar(method)) "auto" else method))
  if(mm %in% c("lm", "limma", "kendall")) return(mm)
  if(identical(tolower(layer_key), "contaminant")) return("kendall")
  if(is_screen) "limma" else "lm"
}

.kgAssocImpl <- function(features = "", phenotypes = "",
                    layers = "all", covariates = "default",
                    subset = "all", contrast = "", ref = "",
                    fdr_family = "query", minN = "10", method = "auto",
                    mode = "tool"){
  library(RSQLite); library(DBI)
  .kgSetPaths()
  minN <- suppressWarnings(as.integer(minN)); if(is.na(minN)) minN <- 10L
  db.path <- paste0(sqlite.path, "HI_omics_v2.sqlite"); if(!file.exists(db.path)) return("RES-NO-DB")

  # parse a group contrast "A-B" once (targeted path passes level A + ref B to .kgLmCell;
  # "" = auto, "anova" = omnibus). "A-B" means A vs B.
  ct_lvl <- contrast; ct_ref <- ref
  if(nzchar(contrast) && !identical(contrast, "anova")){
    pp <- trimws(strsplit(contrast, "-", fixed = TRUE)[[1]])
    if(length(pp) == 2){ ct_lvl <- pp[1]; ct_ref <- pp[2] }
  }

  m <- .kgLoadMetaFrame(); if(is.null(m)) return("RES-NO-META")
  meta_names <- names(m); rec <- as.character(m[["record_id"]])

  # ---- phenotypes: expand parents -> leaves; "all" -> all catalogue leaves ----
  ptoks   <- trimws(strsplit(phenotypes, "[,;]")[[1]]); ptoks <- ptoks[nzchar(ptoks)]
  if(length(ptoks) == 0) return("RES-NO-PHENO")
  all_pheno <- (length(ptoks) == 1 && tolower(ptoks) == "all")
  phenos <- if(all_pheno) .kgAllPhenotypes(meta_names) else .kgExpandPhenotypes(ptoks, meta_names)
  phenos <- phenos[phenos %in% meta_names]; if(length(phenos) == 0) return("RES-NO-PHENO")

  # ---- feature axis -> engine (the list is split later, once the resolver is loaded) ----
  feat_raw  <- trimws(features); if(!nzchar(feat_raw)) return("RES-NO-FEATURE")
  is_screen <- identical(tolower(feat_raw), "all")

  # ---- guard (default #3): all features x all phenotypes recompute -> use precompute ----
  if(is_screen && all_pheno && identical(covariates, "default") && identical(subset, "all"))
    return("RES-USE-PRECOMPUTE")

  # ---- donor subset (applied in both engines) ----
  keep_ids <- rec; subset_note <- "all donors"
  # -- COMPARISONS, NOT JUST EQUALITY (2026-09-04), via the ONE parser in kg_common.R. This block
  #    tested grepl("=", subset) and split on "=", so a filter like "bodymassindex>30" fell straight
  #    through: `subset_note` stayed "all donors" and the analysis ran on EVERY donor while the
  #    caller believed the restriction had been applied. The thresholds come from the question
  #    ("donors over 60", "a BMI over 30"); none is invented here.
  # -- AN UNAPPLIABLE FILTER NOW REFUSES. A column absent from the frame used to skip the block and
  #    return the unrestricted analysis as the answer to a restricted question.
  if(!identical(subset, "all") && nzchar(subset)){
    msub <- .kgDonorSubset(m, subset)
    if(is.null(msub)) return("RES-NO-SUBSET")
    keep_ids <- as.character(msub[["record_id"]]); subset_note <- subset
  }
  keep_ids <- keep_ids[!is.na(keep_ids)]

  # ---- covariate universe (coerced once), and a per-phenotype selector (drop primary) ----
  covuniv <- if(covariates %in% c("none","NA","")) character(0)
             else if(identical(covariates,"default")) .kgAssocDefaultCovs
             else trimws(strsplit(covariates, "[,;]")[[1]])
  covuniv <- unique(covuniv[nzchar(covuniv) & covuniv %in% meta_names])
  covbase <- if(length(covuniv)){ cb <- .kgCoerce(m, covuniv); rownames(cb) <- rec; cb } else NULL
  covcols_for <- function(ph) setdiff(covuniv, ph)

  # ---- phenotype values (coerced once) ----
  phframe <- .kgCoerce(m, phenos); rownames(phframe) <- rec

  con <- dbConnect(SQLite(), db.path); on.exit(dbDisconnect(con), add = TRUE)
  rows <- list()

  parse_contrast <- function(){
    if(!nzchar(contrast) || identical(contrast, "anova")) return(list(c = "anova", r = ref))
    p <- trimws(strsplit(contrast, "-", fixed = TRUE)[[1]])
    if(length(p) == 2) list(c = p[1], r = p[2]) else list(c = "anova", r = ref)
  }

  # ---- resolve the named features (ANY entity) BEFORE choosing the engine ----
  feats <- list(); ftoks <- character(0)
  if(!is_screen){
    res <- .kgLoadFeatureResolver(); if(is.null(res)) return("RES-NO-RESOLVER")
    ftoks <- .kgSplitFeatures(feat_raw, res); ftoks <- ftoks[nzchar(ftoks)]
    if(length(ftoks) == 0) return("RES-NO-FEATURE")
    lay_filter <- if(identical(layers, "all")) NULL else tolower(trimws(strsplit(layers, "[,;]")[[1]]))

    # resolve each feature -> its layers. ANY entity type (gene | metabolite | contaminant |
    # reaction) -- the resolver carries the source table, and the model is the same.
    # ONE vectorised pass for the whole set (a per-token scan costs ~9ms each).
    rmap <- .kgResolveFeatures(res, ftoks, "auto")
    for(tk in ftoks){
      r <- rmap[[.kgNorm(tk)]]
      if(is.null(r) || length(r$layers) == 0 || is.na(r$type)) next
      for(la in r$layers){
        lk <- .kgLayerKey(la$display)
        # entity-driven + supported-only: the resolver hands back this feature's OWN layers,
        # and anything out of scope (pbrna_delta/pbrna_pp/flux) is dropped here.
        if(!(lk %in% .kgAllLayerKeys)) next
        if(!is.null(lay_filter) && !(lk %in% lay_filter || tolower(la$display) %in% lay_filter)) next
        feats[[length(feats)+1]] <- list(token = tk, etype = r$type, layer = la)
      }
    }
    if(length(feats) == 0) return("RES-NO-RESOLVE")
  }

  # ---- ENGINE CHOICE: named set -> lm, until it is big enough to be a screen in all but
  # name (.kgLmMaxFeatures), past which the full-layer limma fit is both better and faster.
  use_limma <- is_screen || (length(ftoks) >= .kgLmMaxFeatures)

  if(!use_limma){
    # =========================== TARGETED engine (lm) ===========================
    # group by layer table so each layer is read ONCE (IN-query over its features)
    key <- vapply(feats, function(f) paste(f$layer$display, f$layer$table, f$layer$read_col, sep = "@@"), character(1))
    for(k in unique(key)){
      idx <- which(key == k); la <- feats[[idx[1]]]$layer
      read_vals <- unique(vapply(feats[idx], function(f) as.character(f$layer$read_val), character(1)))
      FR <- .kgFeatureRows(con, la$table, la$read_col, read_vals); if(is.null(FR)) next
      dcommon <- intersect(colnames(FR$mat), keep_ids); if(length(dcommon) < minN) next

      # small-n entities (metabolite/contaminant) drop the DEFAULT covariates; an explicit
      # user list is always honored. Entity is a property of the layer, so this is per-layer.
      small <- identical(covariates, "default") && (feats[[idx[1]]]$etype %in% .kgSmallNTypes)
      # The design pieces depend only on (layer, phenotype) -- NOT on the feature -- so build
      # them ONCE per layer instead of re-subsetting per feature (that churn dominated a
      # many-feature query far more than the lm's themselves).
      P <- list()
      for(ph in phenos){
        cvn <- if(small) character(0) else covcols_for(ph)
        # ⚠ A DONOR SUBSET CAN LEAVE A COVARIATE CONSTANT ON THESE DONORS ("within males, adjusted
        # for sex"), and the fit drops it -- see `.kgVaryingCovs` in kg_common.R. Narrow it HERE
        # too, because `cvn` is also what the `Covariates` column reports: otherwise the row claims
        # an adjustment the model never made.
        if(length(cvn) && !is.null(covbase)) cvn <- .kgVaryingCovs(covbase[dcommon, , drop = FALSE], cvn)
        phv <- phframe[dcommon, ph]; names(phv) <- dcommon
        P[[ph]] <- list(pt = .kgPhenoType(m, ph), cvn = cvn, phv = phv,
                        cvd = if(length(cvn) && !is.null(covbase)) covbase[dcommon, cvn, drop = FALSE] else NULL)
      }
      rpos <- match(vapply(feats[idx], function(f) as.character(f$layer$read_val), character(1)), FR$id)
      # method depends only on the LAYER, so resolve it once here. `la$display` is the display
      # name, so map it back to the short key the resolver understands.
      meth    <- .kgAssocMethod(method, .kgLayerKey(la$display), is_screen)
      cellfun <- if(identical(meth, "kendall")) .kgKendallCell else .kgLmCell

      for(j in seq_along(idx)){
        f  <- feats[[idx[j]]]; ri <- rpos[j]; if(is.na(ri)) next
        y  <- FR$mat[ri, dcommon]; names(y) <- dcommon; sym <- FR$symbol[ri]
        for(ph in phenos){
          pp <- P[[ph]]; pt <- pp$pt; cvn <- pp$cvn
          cell <- cellfun(y, pp$phv, pp$cvd, pt, ct_lvl, ct_ref, minN); if(is.null(cell)) next
          rows[[length(rows)+1]] <- data.frame(
            Feature = f$token, Symbol = sym, Entity = f$etype, Layer = la$display, Phenotype = ph,
            N = cell$N, Effect = signif(cell$Effect, 3), EffectType = cell$EffectType,
            Statistic = signif(cell$t, 3), P_value = signif(cell$P, 4),
            Direction = cell$Direction,
            Covariates = if(length(cvn)) paste(cvn, collapse = ",") else "none",
            Comparison = if(pt == "disc") (if(nzchar(contrast)) contrast else paste0("auto(", cell$EffectType, ")")) else "continuous",
            Subset = subset_note, Method = meth, stringsAsFactors = FALSE)
        }
      }
    }

  } else {
    # ===================== LIMMA engine (ALWAYS a FULL-LAYER fit) =====================
    # Two callers: (a) features="all" -> report every feature; (b) a NAMED set >= the
    # threshold -> fit the FULL layer, then EXTRACT the named rows. Never subset-then-fit:
    # eBayes estimates its variance prior + mean-variance trend from the genes in the fit, so
    # fitting only a named (usually biased) subset would make each p-value depend on which
    # other features happened to be named. Fit all -> subset the RESULTS.
    work <- list()
    if(is_screen){
      # "all" on a SCREEN = the GENE layers (entity-coherent); a metabolite/contaminant/flux
      # screen is opt-in by naming its layer (e.g. layers="metabolite_hg").
      lays <- if(identical(layers, "all")) .kgGeneLayerKeys else tolower(trimws(strsplit(layers, "[,;]")[[1]]))
      for(ln in lays){
        lp <- .kgScreenLayer(ln); if(is.null(lp)) next
        work[[length(work)+1]] <- list(disp = ln, tbl = lp[1], id = lp[2], keep = NULL,
                                       ent = .kgLayerEntity(tolower(ln)))
      }
    } else {
      key <- vapply(feats, function(f) paste(f$layer$display, f$layer$table, f$layer$read_col, sep = ""), character(1))
      for(k in unique(key)){
        idx <- which(key == k); la <- feats[[idx[1]]]$layer
        work[[length(work)+1]] <- list(
          disp = la$display, tbl = la$table, id = la$read_col,
          keep = unique(vapply(feats[idx], function(f) as.character(f$layer$read_val), character(1))),
          ent  = feats[[idx[1]]]$etype)
      }
    }
    for(w in work){
      L <- .kgLoadLayerMatrix(con, w$tbl, w$id); if(is.null(L)) next
      dcommon <- intersect(colnames(L$mat), keep_ids); if(length(dcommon) < minN) next
      Lmat <- L$mat[, dcommon, drop = FALSE]
      ent   <- w$ent
      small <- identical(covariates, "default") && (ent %in% .kgSmallNTypes)
      meth <- .kgAssocMethod(method, w$disp, is_screen)
      for(ph in phenos){
        pt  <- .kgPhenoType(m, ph)
        cvn <- if(small) character(0) else covcols_for(ph)
        md  <- .kgCoerce(m, unique(c(ph, cvn))); rownames(md) <- rec
        md  <- md[intersect(rownames(md), dcommon), , drop = FALSE]
        # narrowed on the fitted rows, for the same reason as the screen path above: `cvn` is what
        # the `Covariates` column reports, and the fit has already dropped anything constant here.
        cvn <- .kgVaryingCovs(md, cvn)
        # ONE choice of fit; the two functions return the same columns (kg_common.R).
        fitfun <- if(identical(meth, "kendall")) .kgKendallFit else .kgLimmaFit
        if(pt == "disc"){
          pc <- parse_contrast()
          # omnibus needs a real reference level; default to the first level present
          rr <- pc$r
          if(identical(pc$c, "anova") && !nzchar(rr)){
            lv <- levels(factor(as.character(md[[ph]]))); if(length(lv) < 2) next; rr <- lv[1]
          }
          de <- if(identical(meth, "kendall"))
                  .kgKendallFit(Lmat, L$info, md, ph, pt, ref = rr, contrast = pc$c, minN = minN)
                else fitfun(Lmat, L$info, md, ph, pt, ref = rr, contrast = pc$c)
        } else {
          de <- if(identical(meth, "kendall"))
                  .kgKendallFit(Lmat, L$info, md, ph, pt, "", "", minN = minN)
                else fitfun(Lmat, L$info, md, ph, pt, "", "")
        }
        if(is.null(de) || nrow(de) == 0) next
        # EXTRACT the named rows AFTER the full-layer fit (prior/trend stay correct).
        if(!is.null(w$keep)){
          de <- de[as.character(de$FeatureId) %in% w$keep, , drop = FALSE]
          if(nrow(de) == 0) next
        }
        de <- data.frame(
          Feature = de$FeatureId, Symbol = de$Symbol, Entity = ent, Layer = w$disp, Phenotype = ph,
          N = de$N, Effect = signif(de$Effect, 3), EffectType = de$EffectType,
          Statistic = signif(de$T_statistic, 3), P_value = signif(de$P_value, 4),
          Direction = ifelse(is.na(de$Effect), "NA", ifelse(de$Effect > 0, "up", "down")),
          Covariates = if(length(cvn)) paste(cvn, collapse = ",") else "none",
          Comparison = if(pt == "disc") (if(nzchar(contrast)) contrast else "anova") else "continuous",
          Subset = subset_note, Method = meth,
          Adjusted_p_value = signif(de$Adjusted_p_value, 4), stringsAsFactors = FALSE)
        rows[[length(rows)+1]] <- de
      }
    }
  }

  if(length(rows) == 0) return("RES-NO")
  out <- do.call(rbind, rows)

  # ---- family FDR (default #1: BH across exactly the tests asked; never genome-wide) ----
  if(identical(fdr_family, "query")){
    out$P_family <- signif(stats::p.adjust(out$P_value, method = "BH"), 4)
  } else if(identical(fdr_family, "genomewide") && "Adjusted_p_value" %in% names(out)){
    out$P_family <- out$Adjusted_p_value
  }

  out <- out[order(out$P_value), ]
  ts  <- format(Sys.time(), "%Y%m%d_%H%M%S"); safe <- function(s) gsub("[^A-Za-z0-9]", "_", substr(s, 1, 30))
  utils::write.csv(out, paste0("kg_assoc_", safe(features), "_", safe(phenotypes), "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_assoc.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out))
}

# ---- crash guard (.kgGuard, kg_common.R) -------------------------------------
# The servlet calls the PUBLIC name; the body now lives in the .*Impl above. A thrown
# R error comes back as "RES-ERR;<fn>;<message>" instead of a null -> empty HTTP body.
# Every RES-OK / RES-NO* / RES-TOO-* RETURN passes through untouched: tryCatch sees
# conditions, not return values, so "no data" and "crashed" stay distinct answers.
kgAssoc <- function(...) .kgGuard("kgAssoc", .kgAssocImpl, ...)
