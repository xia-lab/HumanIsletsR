# ==============================================================================
# kg_functions/kg_roc.R  --  kgPanelROC: ONE cross-validated ROC for a set of features
# taken TOGETHER (op#14).
#
#   features    the RANKED, SIGNIFICANT features the answer is showing -- a comma list.
#               NEVER "all": this panel describes the set the reader is looking at.
#   outcome     a metadata column holding the groups (e.g. "final_cluster", "diagnosis").
#   levels      "A,B" -- the TWO classes. A is the POSITIVE class.
#   layer       the omics layer the features were measured in (resolver display name or a
#               key). ONE layer only -- see the complete-case note below.
#   covariates  "default" (the 5) | "none" | list. Applied by RESIDUALISING each feature
#               on them (.kgResidualize, the same adjustment the correlation ops use),
#               NOT by adding them as predictors: a covariate in the design would put its
#               own separating power into the AUC and the curve would stop being about the
#               features.
#   folds/repeats/seed  the PRECOMPUTE'S OWN SCHEME: RepeatedStratifiedKFold, 5 folds,
#               5 repeats (3 when n is small), seed 42 -- so this AUC is comparable with
#               `model_quality` already stored on COHORT_PREDICTS instead of contradicting it.
#   boot        bootstrap resamples for the CI (0 = point estimate only).
#
# THE MODEL: cv.glmnet logistic elastic-net -- the R equivalent of the `logistic_en`
# family the precompute uses for every binary target. Lambda is chosen by an INNER
# cross-validation on the TRAINING folds only.
#
# !! THE CURVE IS BUILT FROM OUT-OF-FOLD PROBABILITIES, NEVER IN-SAMPLE. A panel of ~25
# features on ~94 donors fitted and scored on the same rows draws a near-perfect curve
# that means nothing. Every probability here was produced by a model that had not seen
# that donor.
#
# !! AND THE HONEST LIMIT THE CALLER MUST PRINT: the features were SELECTED for separating
# these same donors (they are the significant ones), so even out-of-fold this AUC is
# optimistic as an estimate for NEW donors. It is a fair description of the set on screen;
# it is not a validation. The selection is not repeated inside the folds -- doing so would
# describe a different set than the table shows, which is not the question asked.
#
# !! ONE LAYER ONLY. Feature coverage differs by layer (protein ~234 donors, RNA-seq ~371),
# so a joint fit across layers silently drops every donor missing any one feature. The
# caller passes the layer the ranked features came from.
#
# Writes kg_roc.csv -- one row per curve point (fpr, tpr, tpr_lo, tpr_hi) with the summary
# repeated on every row. Returns "RES-OK;<n_points>" | "RES-NO" | "RES-NO-DB"
#   | "RES-NO-META" | "RES-NO-OUTCOME" | "RES-NO-LEVELS" | "RES-NO-FEATURE"
#   | "RES-NO-RESOLVE" | "RES-NO-RESOLVER" | "RES-TOO-FEW;<n>;<minN>" | "RES-ONE-CLASS".
# ==============================================================================

.kgRocDefaultCovs <- c("donorage", "donorsex", "bodymassindex",
                       "predistributionculturetime")   # purity dropped 2026-09-04 to match the precompute default (4 covariates)

# Stratified fold assignment -- each class split across the folds in equal shares, so a
# fold can never end up single-class (which would make its AUC undefined).
.kgStratFolds <- function(y, k){
  f <- integer(length(y))
  for(cl in unique(y)){
    idx <- which(y == cl)
    f[idx[sample.int(length(idx))]] <- rep_len(seq_len(k), length(idx))
  }
  f
}

# AUC from ranks (Mann-Whitney) -- identical to the ROC AUC, and the SAME estimator
# kg_prediction.R's .kgAUC already uses, so two panels never disagree on a number.
.kgRocAUC <- function(p, y){
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if(n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# The ROC on a FIXED grid of false-positive rates, so bootstrap replicates are averageable
# (curves from different resamples have different numbers of points otherwise).
.kgRocCurve <- function(p, y, grid){
  ord <- order(p, decreasing = TRUE)
  y <- y[ord]
  tp <- cumsum(y == 1); fp <- cumsum(y == 0)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if(n1 == 0 || n0 == 0) return(rep(NA_real_, length(grid)))
  tpr <- c(0, tp / n1); fpr <- c(0, fp / n0)
  stats::approx(fpr, tpr, xout = grid, ties = "ordered", rule = 2)$y
}

.kgPanelROCImpl <- function(features, outcome, levels, layer = "all",
                       covariates = "default", use_raw = FALSE,
                       folds = 5, repeats = 5, seed = 42, boot = 1000,
                       minN = 30, alpha = 0.5, mode = "tool",
                       subset = "all"){   # `subset` LAST: positional callers unaffected
  .kgSetPaths()
  suppressPackageStartupMessages({ library(DBI); library(RSQLite); library(glmnet) })

  # every numeric arrives as a STRING over Rserve (KgRCenter assigns character values)
  .int <- function(x, d){ v <- suppressWarnings(as.integer(x)); if(is.na(v)) d else v }
  .num <- function(x, d){ v <- suppressWarnings(as.numeric(x)); if(is.na(v)) d else v }
  folds <- .int(folds, 5L); repeats <- .int(repeats, 5L); seed <- .int(seed, 42L)
  boot  <- .int(boot, 1000L); minN <- .int(minN, 30L); alpha <- .num(alpha, 0.5)
  use_raw <- isTRUE(use_raw) || identical(tolower(as.character(use_raw)), "true")

  db.path <- paste0(sqlite.path, "HI_omics_v2.sqlite")
  if(!file.exists(db.path)) return("RES-NO-DB")
  m <- .kgLoadMetaFrame(use_raw); if(is.null(m) || !nrow(m)) return("RES-NO-META")
  # the donor row filter -- ONE implementation, `.kgDonorSubset` in kg_common.R. Applied right after
  # the frame loads, so every vector, every N and every model below describes the SAME donor set.
  # NULL means the filter could not be applied, and that REFUSES: returning the unrestricted
  # analysis as the answer to a restricted question is a wider answer wearing the asked one's face.
  m <- .kgDonorSubset(m, subset); if(is.null(m)) return("RES-NO-SUBSET")
  meta_names <- names(m); rec <- as.character(m[["record_id"]])

  if(!nzchar(outcome) || !(outcome %in% meta_names)) return("RES-NO-OUTCOME")
  lv <- trimws(strsplit(levels, "[,;]")[[1]]); lv <- lv[nzchar(lv)]
  if(length(lv) != 2) return("RES-NO-LEVELS")
  ovals <- as.character(m[[outcome]]); names(ovals) <- rec
  keep <- rec[!is.na(ovals) & ovals %in% lv]
  if(length(keep) < minN) return(sprintf("RES-TOO-FEW;%d;%d", length(keep), minN))
  y_all <- stats::setNames(as.integer(ovals[keep] == lv[1]), keep)   # lv[1] = positive

  ftoks <- trimws(strsplit(features, "[,;]")[[1]]); ftoks <- ftoks[nzchar(ftoks)]
  ftoks <- unique(ftoks[!(tolower(ftoks) %in% c("all", ""))])
  if(length(ftoks) == 0) return("RES-NO-FEATURE")
  res <- .kgLoadFeatureResolver(); if(is.null(res)) return("RES-NO-RESOLVER")
  feats <- .kgResolveFeatures(res, ftoks, "auto")
  if(length(feats) == 0) return("RES-NO-RESOLVE")
  lay_filter <- if(identical(layer, "all")) NULL
                else tolower(trimws(strsplit(layer, "[,;]")[[1]]))

  # ---- ONE READ PER LAYER, not one per feature (the N+1 already fixed elsewhere) -------
  con <- dbConnect(SQLite(), db.path); on.exit(dbDisconnect(con), add = TRUE)
  groups <- list()
  for(tk in names(feats)){
    ls <- feats[[tk]]$layers
    if(!is.null(lay_filter)) ls <- Filter(function(l) tolower(l$display) %in% lay_filter, ls)
    if(length(ls)){
      la <- ls[[1]]                       # ONE layer per feature: the panel is layer-scoped
      key <- paste(la$table, la$read_col, sep = "\r")
      if(is.null(groups[[key]]))
        groups[[key]] <- list(table = la$table, read_col = la$read_col,
                              toks = character(0), vals = character(0))
      groups[[key]]$toks <- c(groups[[key]]$toks, tk)
      groups[[key]]$vals <- c(groups[[key]]$vals, la$read_val)
    }
  }
  if(length(groups) == 0) return("RES-NO-RESOLVE")
  cols <- list()
  for(g in groups){
    fr <- .kgFeatureRows(con, g$table, g$read_col, g$vals)
    if(is.null(fr)) next
    for(i in seq_along(g$toks)){
      ri <- match(g$vals[i], rownames(fr$mat)); if(is.na(ri)) next
      v <- fr$mat[ri, ]; names(v) <- colnames(fr$mat)
      cols[[g$toks[i]]] <- v
    }
  }
  if(length(cols) == 0) return("RES-NO-RESOLVE")

  # ---- ONE complete-case donor set across EVERY feature + covariate + the outcome ------
  don <- keep
  for(v in cols) don <- intersect(don, names(v)[!is.na(v)])
  covuniv <- if(covariates %in% c("none", "NA", "")) character(0)
             else if(identical(covariates, "default")) .kgRocDefaultCovs
             else trimws(strsplit(covariates, "[,;]")[[1]])
  covuniv <- unique(covuniv[nzchar(covuniv) & covuniv %in% meta_names])
  covdf <- NULL
  if(length(covuniv)){
    cb <- .kgCoerce(m, covuniv); rownames(cb) <- rec
    cb <- cb[don, , drop = FALSE]
    don <- don[stats::complete.cases(cb)]
    covdf <- cb[don, , drop = FALSE]
  }
  if(length(don) < minN) return(sprintf("RES-TOO-FEW;%d;%d", length(don), minN))
  y <- y_all[don]
  if(length(unique(y)) < 2) return("RES-ONE-CLASS")

  X <- do.call(cbind, lapply(cols, function(v) as.numeric(v[don])))
  colnames(X) <- names(cols)
  # !! ADJUSTMENT BY RESIDUALISATION, the same .kgResidualize the correlation ops use -- the
  # covariates shape the FEATURES, they never become predictors of their own.
  if(!is.null(covdf) && ncol(covdf)){
    for(j in seq_len(ncol(X))){
      rj <- try(.kgResidualize(stats::setNames(X[, j], don), covdf), silent = TRUE)
      if(!inherits(rj, "try-error") && length(rj) == length(don)) X[, j] <- as.numeric(rj)
    }
  }
  keep_col <- apply(X, 2, function(z) all(is.finite(z)) && stats::sd(z) > 0)
  X <- X[, keep_col, drop = FALSE]
  if(ncol(X) == 0) return("RES-NO")

  # ---- OUT-OF-FOLD PROBABILITIES: repeated stratified k-fold, the precompute's scheme ---
  set.seed(seed)
  n <- length(don); oof <- matrix(NA_real_, n, repeats)
  for(rp in seq_len(repeats)){
    fold <- .kgStratFolds(y, folds)
    for(k in seq_len(folds)){
      te <- which(fold == k); tr <- which(fold != k)
      if(length(unique(y[tr])) < 2 || !length(te)) next
      fit <- try(glmnet::cv.glmnet(X[tr, , drop = FALSE], y[tr], family = "binomial",
                                   alpha = alpha, nfolds = min(5, length(tr) - 1)),
                 silent = TRUE)
      if(inherits(fit, "try-error")) next
      pr <- try(as.numeric(stats::predict(fit, newx = X[te, , drop = FALSE],
                                          s = "lambda.min", type = "response")), silent = TRUE)
      if(!inherits(pr, "try-error") && length(pr) == length(te)) oof[te, rp] <- pr
    }
  }
  # average the repeats per donor -- each repeat is an independent out-of-fold estimate
  p_oof <- rowMeans(oof, na.rm = TRUE)
  ok <- is.finite(p_oof)
  if(sum(ok) < minN || length(unique(y[ok])) < 2) return("RES-NO")
  p_oof <- p_oof[ok]; y_ok <- y[ok]

  grid <- seq(0, 1, by = 0.01)
  tpr  <- .kgRocCurve(p_oof, y_ok, grid)
  auc  <- .kgRocAUC(p_oof, y_ok)

  # ---- CI: bootstrap over DONORS on the out-of-fold probabilities ----------------------
  auc_lo <- NA_real_; auc_hi <- NA_real_
  tpr_lo <- rep(NA_real_, length(grid)); tpr_hi <- tpr_lo
  if(boot > 0){
    nn <- length(y_ok); aucs <- rep(NA_real_, boot)
    band <- matrix(NA_real_, boot, length(grid))
    for(b in seq_len(boot)){
      i <- sample.int(nn, nn, TRUE)
      if(length(unique(y_ok[i])) < 2) next
      aucs[b]  <- .kgRocAUC(p_oof[i], y_ok[i])
      band[b, ] <- .kgRocCurve(p_oof[i], y_ok[i], grid)
    }
    aucs <- aucs[is.finite(aucs)]
    if(length(aucs) >= 100){
      qq <- stats::quantile(aucs, c(0.025, 0.975), names = FALSE)
      auc_lo <- qq[1]; auc_hi <- qq[2]
      tpr_lo <- apply(band, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
      tpr_hi <- apply(band, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
    }
  }

  out <- data.frame(
    fpr = grid, tpr = tpr, tpr_lo = tpr_lo, tpr_hi = tpr_hi,
    auc = auc, auc_lo = auc_lo, auc_hi = auc_hi,
    n_donors = length(y_ok), n_pos = sum(y_ok == 1), n_neg = sum(y_ok == 0),
    n_features = ncol(X), outcome = outcome, level_pos = lv[1], level_neg = lv[2],
    layer = layer, covariates = if(length(covuniv)) paste(covuniv, collapse = ";") else "none",
    folds = folds, repeats = repeats, alpha = alpha, boot = boot,
    model = "logistic elastic-net (cv.glmnet), out-of-fold probabilities",
    caveat = paste("features were SELECTED as significant on these same donors, so the AUC is",
                   "optimistic for new donors; this describes the displayed set, it is not a",
                   "validation"),
    stringsAsFactors = FALSE)

  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  safe <- function(s) gsub("[^A-Za-z0-9]", "_", substr(s, 1, 30))
  utils::write.csv(out, paste0("kg_roc_", safe(outcome), "_", safe(levels), "_", ts, ".csv"),
                   row.names = FALSE)
  utils::write.csv(out, "kg_roc.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out))
}

# ---- crash guard (.kgGuard, kg_common.R) -------------------------------------
# The servlet calls the PUBLIC name; the body now lives in the .*Impl above. A thrown
# R error comes back as "RES-ERR;<fn>;<message>" instead of a null -> empty HTTP body.
# Every RES-OK / RES-NO* / RES-TOO-* RETURN passes through untouched: tryCatch sees
# conditions, not return values, so "no data" and "crashed" stay distinct answers.
kgPanelROC <- function(...) .kgGuard("kgPanelROC", .kgPanelROCImpl, ...)
