# ==============================================================================
# kg_functions/kg_corr.R  --  kgCorr: the unified CORRELATION entry point.
# Consolidates F31 (feature<->feature), F32 (phenotype<->phenotype), F19 (cross-omics
# concordance = the SAME token across its own layers) and the previously-missing
# feature<->phenotype correlation, and adds ONE-VS-ALL. One interface, two sides.
#
#   a, b        each side is: a feature token | a phenotype column | a PARENT concept
#               (e.g. "gsis" -> its children) | a "|"-separated list | "all"
#   a_type/b_type  "auto" | "feature" | "phenotype" (also gene|metabolite|contaminant|
#               reaction to force an entity). "auto" = a phenotype column/group if the
#               token is one, else a feature. "all" REQUIRES the side's type.
#   layers      "all" | layer keys -- applies to the FEATURE side(s), entity-driven.
#   covariates  "none" (default -- a raw correlation is not covariate-adjusted unless
#               asked) | "default" | a list  -> PARTIAL correlation (both sides
#               residualised on the covariates, then correlated).
#   subset      "all" | "col=value" (donor subset before correlating)
#   method      "both" | "pearson" | "spearman"
#   fdr_family  "query" (BH over exactly the pairs asked -> P_family) | "none"
#
# CORRELATION IS FOR CONTINUOUS DATA: categorical phenotypes are skipped (use kgAssoc's
# group contrast for those). Self-pairs (same token AND same layer) are dropped; the same
# token across DIFFERENT layers is kept -- that is F19 concordance.
#
# "all" is allowed on AT MOST ONE side (a structural guard, not an invented cap): all x all
# would be ~16.7k x 16.7k pairs.
#
# Stats: Pearson p is the exact t-test (identical to cor.test). Spearman p is the standard
# asymptotic t-approximation on the rank correlation (cor.test uses an exact method only for
# small n) -- chosen so one-vs-all stays vectorised instead of 16.7k cor.test calls.
#
# Writes kg_corr.csv; returns "RES-OK;<n_pairs>" | "RES-NO" | "RES-TOO-BROAD"
#   | "RES-NO-DB" | "RES-NO-META" | "RES-NO-RESOLVER" | "RES-NO-SIDE-A/B".
# ==============================================================================

.kgCorrDefaultCovs <- c("donorage", "donorsex", "bodymassindex", "predistributionculturetime")   # purity dropped 2026-09-04 to match the precompute default (4 covariates)

# residualise an (items x donors) matrix on a covariate frame -> partial correlation input.
# Rows with NO NAs share ONE projection (fast matrix fit); rows with NAs fall back to a
# per-row fit (their missing pattern differs, so a shared projection would be wrong).
.kgResidMatrix <- function(M, covdf, minN = 10){
  don <- rownames(covdf)[stats::complete.cases(covdf)]
  don <- intersect(colnames(M), don)
  if(length(don) < minN) return(NULL)
  M <- M[, don, drop = FALSE]
  # a covariate left CONSTANT by a donor subset ("within males, adjusted for sex") makes `~ .` an
  # UNCAUGHT model.matrix error here -- see `.kgVaryingCovs` in kg_common.R. Residualising on the
  # covariates that still vary is the same fit: a constant column is collinear with the intercept.
  cdf <- covdf[don, , drop = FALSE]
  cdf <- cdf[, .kgVaryingCovs(cdf, names(cdf)), drop = FALSE]
  X <- if(ncol(cdf)) stats::model.matrix(stats::as.formula("~ ."), data = cdf)
       else matrix(1, nrow = length(don), ncol = 1, dimnames = list(don, "(Intercept)"))
  if(nrow(X) < ncol(X) + 2) return(NULL)
  R  <- M
  ok <- rowSums(is.na(M)) == 0
  if(any(ok)){
    fit <- stats::lm.fit(X, t(M[ok, , drop = FALSE]))
    R[ok, ] <- t(fit$residuals)
  }
  for(i in which(!ok)){
    y <- M[i, ]; cc <- !is.na(y)
    if(sum(cc) < ncol(X) + 2){ R[i, ] <- NA_real_; next }
    f <- try(stats::lm.fit(X[cc, , drop = FALSE], y[cc]), silent = TRUE)
    R[i, ] <- NA_real_
    if(!inherits(f, "try-error")) R[i, cc] <- f$residuals
  }
  R
}

# Stack per-layer matrices on the UNION of donors, padding absent donors with NA. Layers do
# NOT share a donor set (rnaseq 371, nanostring 190, protein 234, pbrna_beta 25), so
# intersecting across them would collapse to almost nothing; NA + pairwise-complete is right.
.kgCorrBind <- function(mats){
  don <- sort(unique(unlist(lapply(mats, colnames))))
  do.call(rbind, lapply(mats, function(M){
    o <- matrix(NA_real_, nrow(M), length(don), dimnames = list(rownames(M), don))
    o[, colnames(M)] <- M
    o
  }))
}

# Build ONE side -> list(mat = items x donors, meta = data.frame(Token,Symbol,Id,Kind,Entity,Layer)).
# Id + Layer are the CANONICAL keys used to drop self-pairs (Token is only for display).
.kgCorrSide <- function(con, res, m, rec, spec, tp, lay_filter, phenos_all, minN){
  spec <- trimws(spec); if(!nzchar(spec)) return(NULL)
  is_all <- identical(tolower(spec), "all")
  tp     <- tolower(tp)
  want_pheno <- tp %in% c("phenotype", "pheno")
  want_feat  <- tp %in% c("feature", "gene", "metabolite", "contaminant", "reaction")
  if(is_all && !want_pheno && !want_feat) return(NULL)   # "all" needs an explicit side type

  # ---------------- phenotype side ----------------
  build_pheno <- function(cols){
    cols <- cols[cols %in% names(m)]
    keep <- cols[vapply(cols, function(c) identical(.kgPhenoType(m, c), "cont"), logical(1))]
    if(length(keep) == 0) return(NULL)
    d <- .kgCoerce(m, keep); rownames(d) <- rec
    mat <- t(as.matrix(d[, keep, drop = FALSE])); rownames(mat) <- keep
    list(mat = mat, meta = data.frame(Token = keep, Symbol = keep, Id = keep, Kind = "phenotype",
                                      Entity = "phenotype", Layer = "metadata", stringsAsFactors = FALSE))
  }
  if(is_all && want_pheno) return(build_pheno(phenos_all))

  toks <- .kgSplitFeatures(spec, res); toks <- toks[nzchar(toks)]
  if(!is_all && !want_feat){
    ph <- .kgExpandPhenotypes(toks, names(m))
    if(length(ph) > 0 && !want_feat) { s <- build_pheno(ph); if(!is.null(s)) return(s) }
    if(want_pheno) return(NULL)
  }

  # ---------------- feature side ----------------
  et_hint <- if(tp %in% c("gene","metabolite","contaminant","reaction")) tp else "auto"
  if(is_all){
    # every feature of the requested layer(s); "all" layers on a feature screen = gene layers
    lays <- if(is.null(lay_filter)) .kgGeneLayerKeys else lay_filter
    out <- list()
    for(ln in lays){
      lp <- .kgScreenLayer(ln); if(is.null(lp)) next
      L  <- .kgLoadLayerMatrix(con, lp[1], lp[2]); if(is.null(L)) next
      rownames(L$mat) <- paste0(L$info$FeatureId, "@@", ln)
      out[[length(out)+1]] <- list(mat = L$mat, meta = data.frame(
        Token = L$info$FeatureId, Symbol = L$info$Symbol, Id = L$info$FeatureId, Kind = "feature",
        Entity = .kgLayerEntity(tolower(ln)), Layer = ln, stringsAsFactors = FALSE))
    }
    if(length(out) == 0) return(NULL)
    return(list(mat  = .kgCorrBind(lapply(out, function(o) o$mat)),
                meta = do.call(rbind, lapply(out, function(o) o$meta))))
  }

  rmap <- .kgResolveFeatures(res, toks, et_hint)
  items <- list()
  for(tk in toks){
    r <- rmap[[.kgNorm(tk)]]; if(is.null(r) || is.na(r$type)) next
    for(la in r$layers){
      lk <- .kgLayerKey(la$display)
      if(!(lk %in% .kgAllLayerKeys)) next
      if(!is.null(lay_filter) && !(lk %in% lay_filter || tolower(la$display) %in% lay_filter)) next
      items[[length(items)+1]] <- list(token = tk, etype = r$type, layer = la, key = lk)
    }
  }
  if(length(items) == 0) return(NULL)
  grp <- vapply(items, function(f) paste(f$layer$table, f$layer$read_col, sep = "@@"), character(1))
  out <- list()
  for(g in unique(grp)){
    idx <- which(grp == g); la <- items[[idx[1]]]$layer
    rv  <- unique(vapply(items[idx], function(f) as.character(f$layer$read_val), character(1)))
    FR  <- .kgFeatureRows(con, la$table, la$read_col, rv); if(is.null(FR)) next
    pos <- match(vapply(items[idx], function(f) as.character(f$layer$read_val), character(1)), FR$id)
    keep <- !is.na(pos); if(!any(keep)) next
    mm  <- FR$mat[pos[keep], , drop = FALSE]
    rownames(mm) <- paste0(vapply(items[idx][keep], function(f) f$token, character(1)), "@@",
                           vapply(items[idx][keep], function(f) f$key, character(1)))
    out[[length(out)+1]] <- list(mat = mm, meta = data.frame(
      Token  = vapply(items[idx][keep], function(f) f$token, character(1)),
      Symbol = FR$symbol[pos[keep]],
      Id     = vapply(items[idx][keep], function(f) as.character(f$layer$read_val), character(1)),
      Kind   = "feature",
      Entity = vapply(items[idx][keep], function(f) f$etype, character(1)),
      Layer  = vapply(items[idx][keep], function(f) f$key, character(1)),
      stringsAsFactors = FALSE))
  }
  if(length(out) == 0) return(NULL)
  list(mat  = .kgCorrBind(lapply(out, function(o) o$mat)),
       meta = do.call(rbind, lapply(out, function(o) o$meta)))
}

# ⚠ AT MOST THIS MANY LAYERS ON AN "all" SIDE (user ruling 2026-09-09: "we can only do two layer
# all feature at most, and tell user to click for other layers"). An "all" side loads EVERY feature
# of each requested layer and correlates it against the other side, so the work grows with the
# number of layers, not with the question. MEASURED before this cap: TCF7L2 vs "all" over seven
# layers returned 42,302 rows in 20.3 s — and most of them answered a question nobody asked.
# ⚠ IT REFUSES, IT DOES NOT SILENTLY TRUNCATE. A truncated layer list is indistinguishable from a
# deliberate narrowing — the exact defect that made the dropped-half bug invisible for so long — so
# this returns a code naming the count and the cap, and the caller offers the rest as clicks.
# ⚠ THE CAP FOLLOWS THE MODALITY (user ruling 2026-09-09: "layer b can be 3 for other while first
# 2 for gene"). The cost of an "all" side is FEATURES x DONORS, and the two differ by two orders of
# magnitude: a gene layer loads ~20k features against up to 371 donors, a metabolite layer ~215
# against 37. Two gene layers and three non-gene layers are comparable work — and three is exactly
# the whole metabolomics modality (hg / lg / ratio), so a metabolite screen runs complete in one
# call and needs no click at all.
# ⚠ A MIXED REQUEST TAKES THE STRICTER CAP. If any requested layer is a gene layer the gene cost
# dominates, so 2 applies; nothing here can widen a gene screen by pairing it with a cheap layer.
.kgCorrAllSideLayerCap      <- 2L   # gene layers  (.kgGeneLayerKeys)
.kgCorrAllSideLayerCapOther <- 3L   # every other modality

.kgCorrImpl <- function(a = "", b = "", a_type = "auto", b_type = "auto",
                   layers = "all", covariates = "none", subset = "all",
                   method = "both", fdr_family = "query", minN = "10", mode = "tool",
                   layers_a = "", layers_b = ""){
  # ⚠ PER-SIDE LAYERS (user ruling 2026-09-09). `layers` filters BOTH sides, and both branches of
  # `.kgCorrSide` apply it — the named side drops a feature whose layers are all filtered out
  # (-> RES-NO-SIDE-A) and the "all" side iterates exactly the layers it is given. So a question
  # like "which METABOLITES correlate with TCF7L2" could not be expressed: any list narrow enough
  # to give metabolite partners also strips TCF7L2, and any list wide enough to keep TCF7L2 lets
  # side B return genes. MEASURED on the live backend: metabolite-only -> RES-NO-SIDE-A;
  # gene+metabolite layers -> RES-OK;42302 whose top rows are all gene-gene.
  # ⚠ EMPTY MEANS "USE `layers`", so every existing caller — Java, R, or a direct call — behaves
  # byte-identically and no deployed client has to change.
  library(RSQLite); library(DBI)
  .kgSetPaths()
  minN <- suppressWarnings(as.integer(minN)); if(is.na(minN)) minN <- 10L
  method <- tolower(method)
  db.path <- paste0(sqlite.path, "HI_omics_v2.sqlite"); if(!file.exists(db.path)) return("RES-NO-DB")

  # structural guard: "all" on at most ONE side
  if(identical(tolower(trimws(a)), "all") && identical(tolower(trimws(b)), "all")) return("RES-TOO-BROAD")

  m <- .kgLoadMetaFrame(); if(is.null(m)) return("RES-NO-META")
  rec <- as.character(m[["record_id"]])
  res <- .kgLoadFeatureResolver(); if(is.null(res)) return("RES-NO-RESOLVER")
  phenos_all <- .kgAllPhenotypes(names(m))
  lay_filter <- if(identical(layers, "all")) NULL else tolower(trimws(strsplit(layers, "[,;]")[[1]]))
  # ⚠ EACH SIDE TAKES ITS OWN FILTER, FALLING BACK TO THE SHARED ONE. `layers_a`/`layers_b` are the
  # whole point of the change; an empty one is exactly the old behaviour for that side.
  .side_lay <- function(spec){
    spec <- trimws(as.character(spec))
    if(!nzchar(spec)) return(lay_filter)
    if(identical(tolower(spec), "all")) return(NULL)
    tolower(trimws(strsplit(spec, "[,;]")[[1]]))
  }
  lay_a <- .side_lay(layers_a); lay_b <- .side_lay(layers_b)

  # ⚠ THE CAP APPLIES TO WHICHEVER SIDE IS "all" — see `.kgCorrAllSideLayerCap`. Only that side
  # loads every feature per layer; a NAMED side reads one feature per layer and is cheap however
  # many layers it spans, so capping it would refuse questions that cost nothing.
  .cap_check <- function(spec, lays, side){
    if(!identical(tolower(trimws(as.character(spec))), "all")) return(NULL)
    if(is.null(lays)) return(NULL)
    cap <- if(any(lays %in% .kgGeneLayerKeys)) .kgCorrAllSideLayerCap
           else .kgCorrAllSideLayerCapOther
    if(length(lays) <= cap) return(NULL)
    paste0("RES-TOO-MANY-LAYERS;", side, ";", length(lays), ";", cap,
           ";", paste(lays, collapse = "|"))
  }
  .capA <- .cap_check(a, lay_a, "A"); if(!is.null(.capA)) return(.capA)
  .capB <- .cap_check(b, lay_b, "B"); if(!is.null(.capB)) return(.capB)

  con <- dbConnect(SQLite(), db.path); on.exit(dbDisconnect(con), add = TRUE)
  A <- .kgCorrSide(con, res, m, rec, a, a_type, lay_a, phenos_all, minN); if(is.null(A)) return("RES-NO-SIDE-A")
  B <- .kgCorrSide(con, res, m, rec, b, b_type, lay_b, phenos_all, minN); if(is.null(B)) return("RES-NO-SIDE-B")

  # ---- donors: shared, optionally subset ----
  don <- intersect(colnames(A$mat), colnames(B$mat))
  subset_note <- "all donors"
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
    don <- intersect(don, as.character(msub[["record_id"]])); subset_note <- subset
  }
  don <- don[!is.na(don)]; if(length(don) < minN) return("RES-NO")
  Am <- A$mat[, don, drop = FALSE]; Bm <- B$mat[, don, drop = FALSE]

  # ---- optional PARTIAL correlation (residualise both sides on the covariates) ----
  cov_note <- "none"
  if(!(covariates %in% c("none", "NA", ""))){
    cvn <- if(identical(covariates, "default")) .kgCorrDefaultCovs else trimws(strsplit(covariates, "[,;]")[[1]])
    cvn <- unique(cvn[nzchar(cvn) & cvn %in% names(m)])
    if(length(cvn)){
      cd <- .kgCoerce(m, cvn); rownames(cd) <- rec; cd <- cd[don, , drop = FALSE]
      Ar <- .kgResidMatrix(Am, cd, minN); Br <- .kgResidMatrix(Bm, cd, minN)
      if(is.null(Ar) || is.null(Br)) return("RES-NO")
      dd <- intersect(colnames(Ar), colnames(Br)); if(length(dd) < minN) return("RES-NO")
      Am <- Ar[, dd, drop = FALSE]; Bm <- Br[, dd, drop = FALSE]
      # ⚠ REPORT WHAT WAS ACTUALLY RESIDUALISED ON. `.kgResidMatrix` drops a covariate left
      # constant by the donor subset (see `.kgVaryingCovs`), so naming the requested list here
      # would claim an adjustment that never happened. `cd` itself is left INTACT on purpose --
      # it still decides which donors are complete, so the donor set does not move.
      cov_use  <- .kgVaryingCovs(cd, cvn)
      cov_note <- if(length(cov_use)) paste(cov_use, collapse = ",") else "none"
    }
  }

  At <- t(Am); Bt <- t(Bm)                      # donors x items
  N  <- crossprod(!is.na(At), !is.na(Bt))       # pairwise complete n
  pfun <- function(R){ Tt <- R * sqrt((N - 2) / (1 - R^2)); p <- 2 * stats::pt(-abs(Tt), N - 2); p[!is.finite(p)] <- NA_real_; p }
  Rp <- Rs <- Pp <- Ps <- NULL
  if(method %in% c("both", "pearson")){
    Rp <- suppressWarnings(stats::cor(At, Bt, use = "pairwise.complete.obs")); Pp <- pfun(Rp)
  }
  if(method %in% c("both", "spearman")){
    Rs <- suppressWarnings(stats::cor(At, Bt, method = "spearman", use = "pairwise.complete.obs")); Ps <- pfun(Rs)
  }
  RR <- if(!is.null(Rp)) Rp else Rs
  if(is.null(RR)) return("RES-NO")

  ia <- rep(seq_len(nrow(A$meta)), times = nrow(B$meta))
  ib <- rep(seq_len(nrow(B$meta)), each  = nrow(A$meta))
  out <- data.frame(
    Feature_A = A$meta$Token[ia], Symbol_A = A$meta$Symbol[ia], Kind_A = A$meta$Kind[ia],
    Entity_A = A$meta$Entity[ia], Layer_A = A$meta$Layer[ia],
    Feature_B = B$meta$Token[ib], Symbol_B = B$meta$Symbol[ib], Kind_B = B$meta$Kind[ib],
    Entity_B = B$meta$Entity[ib], Layer_B = B$meta$Layer[ib],
    .IdA = A$meta$Id[ia], .IdB = B$meta$Id[ib],
    N = as.vector(N),
    Pearson_r    = if(!is.null(Rp)) signif(as.vector(Rp), 3) else NA_real_,
    Pearson_p    = if(!is.null(Pp)) signif(as.vector(Pp), 4) else NA_real_,
    Spearman_rho = if(!is.null(Rs)) signif(as.vector(Rs), 3) else NA_real_,
    Spearman_p   = if(!is.null(Ps)) signif(as.vector(Ps), 4) else NA_real_,
    Covariates = cov_note, Subset = subset_note, stringsAsFactors = FALSE)

  # drop self-pairs on the CANONICAL id + layer (a named token and an "all"-side row can be
  # the same feature under different display names).
  out <- out[!(out$.IdA == out$.IdB & out$Layer_A == out$Layer_B), , drop = FALSE]

  # F31 project rule (kept verbatim from the original kgFeatureCorr): two DIFFERENT features
  # of the SAME entity type are comparable only WITHIN one layer -- never across platforms,
  # which would conflate assay effects. Different entity types (gene vs metabolite) and
  # feature-vs-phenotype keep every layer combo. The SAME feature across layers is KEPT:
  # that is exactly F19 cross-omics concordance.
  sk <- out$Kind_A == "feature" & out$Kind_B == "feature" & out$Entity_A == out$Entity_B
  sf <- toupper(out$Symbol_A) == toupper(out$Symbol_B); sf[is.na(sf)] <- FALSE
  out <- out[!(sk & !sf & out$Layer_A != out$Layer_B), , drop = FALSE]

  # correlation is SYMMETRIC -> keep one row per unordered pair (a==b would give both ways)
  ka <- paste0(out$.IdA, "@@", out$Layer_A); kb <- paste0(out$.IdB, "@@", out$Layer_B)
  out <- out[!duplicated(ifelse(ka < kb, paste0(ka, "||", kb), paste0(kb, "||", ka))), , drop = FALSE]
  out$.IdA <- NULL; out$.IdB <- NULL
  out <- out[out$N >= minN, , drop = FALSE]
  prim <- if(!is.null(Pp)) out$Pearson_p else out$Spearman_p
  out  <- out[!is.na(prim), , drop = FALSE]; prim <- prim[!is.na(prim)]
  if(nrow(out) == 0) return("RES-NO")

  if(identical(fdr_family, "query")) out$P_family <- signif(stats::p.adjust(prim, method = "BH"), 4)
  out <- out[order(prim), , drop = FALSE]

  ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); safe <- function(s) gsub("[^A-Za-z0-9]", "_", substr(s, 1, 24))
  utils::write.csv(out, paste0("kg_corr_", safe(a), "_", safe(b), "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_corr.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out))
}

# ---- crash guard (.kgGuard, kg_common.R) -------------------------------------
# The servlet calls the PUBLIC name; the body now lives in the .*Impl above. A thrown
# R error comes back as "RES-ERR;<fn>;<message>" instead of a null -> empty HTTP body.
# Every RES-OK / RES-NO* / RES-TOO-* RETURN passes through untouched: tryCatch sees
# conditions, not return values, so "no data" and "crashed" stay distinct answers.
kgCorr <- function(...) .kgGuard("kgCorr", .kgCorrImpl, ...)
