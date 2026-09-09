# ==============================================================================
# kg_functions/kg_prediction.R
# ------------------------------------------------------------------------------
# Knowledge Search analysis engine -- Tier 2: SINGLE-FEATURE PREDICTION / BIOMARKER
# (F5 "is gene X a predictor" + single-feature biomarker, design sec.3/sec.5/sec.8).
#
# The STANDALONE predictive value of ONE feature for an outcome:
#   * binary outcome  -> ROC AUC (of the feature as a classifier) + Mann-Whitney p
#   * continuous out. -> linear R^2 + coefficient + p
# This is the "single predictor, no penalization" case (a multi-feature panel is
# elastic-net / DIABLO = RECOMMEND, not here). It is deliberately RAW (no covariate
# adjustment): a biomarker question asks how well the feature ALONE separates/predicts.
#
# PURPOSE-BUILT for a query (an entity + an outcome). Reuses the qs resolver +
# .kgFeatureVector from kg_correlation.R; ANY entity type may be a biomarker
# (gene, metabolite, contaminant, flux), so etype is not restricted.
#
# Feature   <- HI_omics_v2 proc_* (R### donors, via the resolver)
# Outcome   <- display_data/metadata_sum_norm.csv (record_id x phenotype)
#
# Returns "RES-OK;<n_rows>" (kg_prediction.csv) | "RES-NO" | "RES-NO-RESOLVE"
#   | "RES-NO-RESOLVER" | "RES-NO-META" | "RES-NO-OUTCOME" (outcome not a column)
#   | "RES-NO-CONTRAST" (a >2-level categorical outcome needs which two classes).
# ==============================================================================

# AUC of a single continuous score x against a binary label g (1=positive, 0=negative),
# computed from ranks (Mann-Whitney) -- identical to the ROC AUC of the feature and to a
# single-predictor logistic model's AUC, with no pROC dependency. Returns NA if degenerate.
.kgAUC <- function(x, g){
  n1 <- sum(g == 1); n0 <- sum(g == 0)
  if(n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(x)
  (sum(r[g == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# ==============================================================================
# kgSingleFeaturePredict
# ------------------------------------------------------------------------------
# var_id     : the predictor feature (gene | metabolite | contaminant | rxn | ids)
# outcome    : a metadata_sum_norm column (continuous or categorical)
# var_type   : "auto" | gene | ensembl | entrez | metabolite | contaminant | reaction
# contrast   : for a categorical outcome with >2 levels -> "A-B" (A = positive class);
#              ignored for continuous or already-binary outcomes
# omics_filter: "all" | pipe-separated display names
# use_raw    : "false" (metadata_sum_norm) | "true" (metadata_sum_raw)
# minN       : minimum complete-case donors (default 10)
# ==============================================================================
kgSingleFeaturePredict <- function(var_id, outcome, var_type = "auto",
                                   contrast = "", omics_filter = "all",
                                   use_raw = "false", minN = "10", mode = "tool",
                                   subset = "all"){
  # `subset` is LAST so every existing positional caller is unaffected.
  library(RSQLite); library(DBI)
  .kgSetPaths()
  minN <- suppressWarnings(as.integer(minN)); if(is.na(minN)) minN <- 10L

  res <- .kgLoadFeatureResolver(); if(is.null(res)) return("RES-NO-RESOLVER")
  r <- .kgResolveFeature(res, var_id, var_type)
  if(length(r$layers) == 0) return("RES-NO-RESOLVE")

  # ---- outcome from metadata_sum_norm ----
  m <- .kgLoadMetaFrame(tolower(use_raw) %in% c("true","1")); if(is.null(m)) return("RES-NO-META")
  # the donor row filter, applied BEFORE the outcome vector is built, so the model is fitted on the
  # donors the question named and `rid` cannot carry anyone outside them. NULL = refuse (the
  # unrestricted prediction is a different question); see .kgDonorSubset in kg_common.R.
  m <- .kgDonorSubset(m, subset); if(is.null(m)) return("RES-NO-SUBSET")
  if(nrow(m) < minN) return("RES-NO")
  if(!(outcome %in% names(m))){
    j <- match(.kgNorm(outcome), tolower(names(m)))
    if(is.na(j)) return("RES-NO-OUTCOME")
    outcome <- names(m)[j]
  }
  rid  <- as.character(m[["record_id"]])
  yraw <- m[[outcome]]
  ynum <- suppressWarnings(as.numeric(yraw))
  n_distinct <- length(unique(yraw[!is.na(yraw) & yraw != "" & yraw != "NA"]))
  is_cont <- (sum(!is.na(ynum)) / max(1, sum(!is.na(yraw) & yraw != "" & yraw != "NA")) > 0.95) && n_distinct > 6

  # binary label vector g (1/0) for a categorical outcome; else NULL
  glab <- NULL; pos_lab <- NULL; neg_lab <- NULL
  if(!is_cont){
    lv <- unique(yraw[!is.na(yraw) & yraw != "" & yraw != "NA"])
    if(nzchar(contrast) && grepl("-", contrast)){
      pr <- trimws(strsplit(contrast, "-", fixed = TRUE)[[1]])
      if(length(pr) != 2) return("RES-NO-CONTRAST")
      pos_lab <- pr[1]; neg_lab <- pr[2]
    } else if(length(lv) == 2){
      pos_lab <- lv[1]; neg_lab <- lv[2]
    } else return("RES-NO-CONTRAST")
    # ⚠ ONE-VS-REST — `neg_lab == "rest"` IS A RESERVED TOKEN, NOT A LEVEL (2026-08-31). Same
    # fix as `.kgLimmaFit` in kg_common.R, and needed for the same reason: F14's one-vs-rest
    # handoff fires op2 as well as op1, and no donor carries a `final_cluster` value of "rest",
    # so the 0-class matched NOTHING — every non-C0 donor stayed NA, the label had one class and
    # the AUC was undefined. The negative class is therefore every MEASURED donor outside the
    # named level, which is the atlas's own construction.
    # ⚠ UNMEASURED IS NOT "REST": a blank / NA / "NA" outcome stays NA and is dropped, exactly as
    # before. ⚠ An ordinary two-level contrast is untouched — the else branch is the old line.
    .yc <- as.character(yraw)
    .meas <- !is.na(.yc) & .yc != "" & .yc != "NA"
    glab <- rep(NA_integer_, length(yraw))
    glab[.meas & .yc == pos_lab] <- 1L
    if(identical(neg_lab, "rest")) glab[.meas & .yc != pos_lab] <- 0L
    else                           glab[.meas & .yc == neg_lab] <- 0L
    names(glab) <- rid
  } else {
    names(ynum) <- rid
  }

  allow <- if(identical(omics_filter,"all")) NULL else trimws(strsplit(omics_filter, "\\|")[[1]])
  con <- dbConnect(SQLite(), paste0(sqlite.path, "HI_omics_v2.sqlite")); on.exit(dbDisconnect(con), add = TRUE)

  rows <- list()
  for(la in r$layers){
    if(!is.null(allow) && !(la$display %in% allow)) next
    fv <- .kgFeatureVector(con, la); if(is.null(fv)) next

    if(is_cont){
      common <- intersect(names(fv), names(ynum))
      x <- fv[common]; y <- ynum[common]
      ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
      if(length(x) < minN) next
      fit <- try(stats::lm(y ~ x), silent = TRUE); if(inherits(fit,"try-error")) next
      sm  <- summary(fit)
      coefp <- if(nrow(sm$coefficients) >= 2) sm$coefficients[2, ] else c(NA,NA,NA,NA)
      rows[[length(rows)+1]] <- data.frame(
        Feature=var_id, Layer=la$display, N=length(x),
        OutcomeType="continuous", Metric="R2", MetricValue=signif(sm$r.squared,3),
        Coefficient=signif(unname(coefp[1]),3), P_value=signif(unname(coefp[4]),3),
        Direction=ifelse(is.na(coefp[1]),"NA",ifelse(coefp[1]>0,"up","down")),
        Outcome=outcome, stringsAsFactors=FALSE)
    } else {
      common <- intersect(names(fv), names(glab))
      x <- fv[common]; g <- glab[common]
      ok <- is.finite(x) & !is.na(g); x <- x[ok]; g <- g[ok]
      if(length(x) < minN || sum(g==1) < 3 || sum(g==0) < 3) next
      auc <- .kgAUC(x, g)
      pv  <- try(suppressWarnings(stats::wilcox.test(x[g==1], x[g==0])$p.value), silent = TRUE)
      if(inherits(pv,"try-error")) pv <- NA_real_
      rows[[length(rows)+1]] <- data.frame(
        Feature=var_id, Layer=la$display, N=length(x),
        OutcomeType="binary", Metric="AUC", MetricValue=signif(auc,3),
        Coefficient=NA_real_, P_value=signif(pv,3),
        Direction=ifelse(is.na(auc),"NA",ifelse(auc>=0.5,paste0("higher in ",pos_lab),paste0("higher in ",neg_lab))),
        Outcome=paste0(outcome," (",pos_lab," vs ",neg_lab,")"), stringsAsFactors=FALSE)
    }
  }
  if(length(rows) == 0) return("RES-NO")

  out <- do.call(rbind, rows)
  # rank by predictive strength: AUC distance from 0.5, or R^2
  ord <- ifelse(out$Metric=="AUC", abs(out$MetricValue-0.5), out$MetricValue)
  out <- out[order(-ord), ]
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  safe <- function(s) gsub("[^A-Za-z0-9]", "_", s)
  utils::write.csv(out, paste0("kg_prediction_", safe(var_id), "_", safe(outcome), "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_prediction.csv", row.names = FALSE)

  if(exists("rcmd")){
    if(!file.exists("savedAnalysis")) dir.create("savedAnalysis")
    write(gsub("tool","local", rcmd), file = "savedAnalysis/Rhistory.R", append = TRUE)
  }
  paste0("RES-OK;", nrow(out))
}
