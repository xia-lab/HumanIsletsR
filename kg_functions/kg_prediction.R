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
# outcome    : a metadata_sum_norm column (continuous or categorical), or several joined by "|" (one call, 2026-10-09)
# var_type   : "auto" | gene | ensembl | entrez | metabolite | contaminant | reaction
# contrast   : for a categorical outcome with >2 levels -> "A-B" (A = positive class);
#              ignored for continuous or already-binary outcomes
# omics_filter: "all" | pipe-separated display names or short keys | a staged spec ("<primary>/<fallback>" groups,
#              2026-10-09) -> the main layers run, the rest written to kg_prediction_layers.csv as `not_calculated`
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
  rid  <- as.character(m[["record_id"]])

  # ⚠ ONE CALL CAN CARRY SEVERAL OUTCOMES, "|"-separated (2026-10-09, user: "yes prepare 3"). A concept such as "beta cell
  # function" has 68 members the donor table holds, and one call per member paid the whole setup (connection, scripts,
  # resolver, donor table) every time — MEASURED through Payara: 1.63 s a call whether 2 layers or 5 were fitted, so
  # "GCK x beta cell function" spent 115 s in op2 and ran past the 120 s build limit. Now the setup is paid once and each
  # member is one more fit per layer. ⚠ A SINGLE OUTCOME IS UNCHANGED: the same checks in the same order, the same codes
  # (RES-NO-OUTCOME, RES-NO-CONTRAST) and the same rows. In a list, a member that is no column or needs a contrast is
  # skipped, and RES-NO-OUTCOME is answered only when no member is a column.
  outs <- trimws(strsplit(as.character(outcome), "|", fixed = TRUE)[[1]]); outs <- outs[nzchar(outs)]
  multi <- length(outs) > 1
  .prep <- function(outcome){
    yraw <- m[[outcome]]
    ynum <- suppressWarnings(as.numeric(yraw))
    n_distinct <- length(unique(yraw[!is.na(yraw) & yraw != "" & yraw != "NA"]))
    is_cont <- (sum(!is.na(ynum)) / max(1, sum(!is.na(yraw) & yraw != "" & yraw != "NA")) > 0.95) && n_distinct > 6
    if(is_cont){ names(ynum) <- rid; return(list(is_cont = TRUE, ynum = ynum)) }
    # binary label vector g (1/0) for a categorical outcome
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
    list(is_cont = FALSE, glab = glab, pos_lab = pos_lab, neg_lab = neg_lab)
  }
  spec <- list()
  for(oc in outs){
    if(!(oc %in% names(m))){
      j <- match(.kgNorm(oc), tolower(names(m)))
      if(is.na(j)){ if(multi) next; return("RES-NO-OUTCOME") }
      oc <- names(m)[j]
    }
    p <- .prep(oc)
    if(is.character(p)){ if(multi) next; return(p) }
    spec[[oc]] <- p
  }
  if(!length(spec)) return("RES-NO-OUTCOME")

  # ⚠ THE FILTER TAKES THE LAYER NAME OR ITS SHORT KEY (2026-10-08). The page sends full names joined by "|"; the
  # agent sends the short keys every other kg function takes ("rnaseq", "pbrna_alpha", "metabolite_hg,metabolite_lg")
  # — MEASURED on F5: all 180 op2 calls that named a layer matched no layer and returned RES-NO. The key is read with
  # kg_assoc.R's `.kgLayerKey` (sourced by KgSingleFeaturePredict.java, as KgSeries/KgCorr do), not a second map.
  allow <- if(identical(omics_filter,"all")) NULL else trimws(strsplit(omics_filter, "[|,]")[[1]])
  con <- dbConnect(SQLite(), paste0(sqlite.path, "HI_omics_v2.sqlite")); on.exit(dbDisconnect(con), add = TRUE)
  # one layer per contaminant TISSUE (`.kgTissueLayers`, kg_common.R), each read through its own row filter and
  # reported under its tissue label; a filter naming the resolver's own "Environmental contaminants" keeps both
  lays <- unlist(lapply(r$layers, .kgTissueLayers), recursive = FALSE)

  # ⚠ A STAGED SPEC RUNS THE MAIN LAYERS AND OFFERS THE REST (user 2026-10-04 — "in default we do only at most three layer
  # Bulk RNA-seq, Protein, Metabolomics HG / LG / ratio, and tell user what is provide but not calculate yet" — the rule a
  # concept's op1 screen follows, applied to op2 on 2026-10-09). The spec is `generate_analysis_v3.CONCEPT_STAGED_LAYERS`
  # ("<primary>/<fallback>" groups split by "|") and kg_assoc.R's own `.kgStagedLayers` decides — on THIS feature's own
  # layers only, so a fallback runs when no main layer the feature is measured in has donors with the outcome. Each
  # offered layer is written to kg_prediction_layers.csv as `not_calculated` with its donor counts, as kg_assoc_layers.csv
  # does, so the answer can say it can still be calculated. A filter without "/" is read exactly as before.
  staged <- grepl("/", omics_filter, fixed = TRUE)
  offered <- character(0)
  if(staged){
    own <- unique(vapply(lays, function(l) .kgLayerKey(l$display), character(1)))
    .own_only <- function(s) paste(intersect(tolower(trimws(strsplit(s, "[,;]")[[1]])), own), collapse = ",")
    gspec <- vapply(strsplit(omics_filter, "|", fixed = TRUE)[[1]], function(g){
      parts <- strsplit(g, "/", fixed = TRUE)[[1]]
      paste0(.own_only(parts[1]), if(length(parts) > 1) paste0("/", .own_only(parts[2])) else "")
    }, character(1))
    phframe <- as.data.frame(lapply(names(spec), function(oc) m[[oc]]), col.names = names(spec),
                             check.names = FALSE, stringsAsFactors = FALSE)
    rownames(phframe) <- rid
    st <- .kgStagedLayers(con, paste(gspec, collapse = "|"), rid, phframe, names(spec), minN)
    allow <- trimws(strsplit(st$run, ",", fixed = TRUE)[[1]]); allow <- allow[nzchar(allow)]
    if(!length(allow)) allow <- "__no_layer__"
    offered <- st$offered
    lst <- list()
    for(k in offered){
      dc <- intersect(.kgLayerDonors(con, k), rid)
      labs <- unique(unlist(lapply(lays, function(l) if(identical(.kgLayerKey(l$display), k)) l$display)))
      if(!length(labs)) labs <- k
      for(lb in labs) for(oc in names(spec)){
        v <- m[[oc]][match(dc, rid)]
        lst[[length(lst)+1]] <- data.frame(
          Layer = lb, Phenotype = oc, N_layer = length(dc), N_usable = sum(!is.na(v) & !(as.character(v) %in% c("", "NA"))),
          MinN = minN, Status = "not_calculated", Subset = if(nzchar(subset) && subset != "all") subset else "all donors",
          stringsAsFactors = FALSE)
      }
    }
    utils::write.csv(if(length(lst)) do.call(rbind, lst) else
                       data.frame(Layer = character(0), Phenotype = character(0), N_layer = integer(0),
                                  N_usable = integer(0), MinN = integer(0), Status = character(0), Subset = character(0)),
                     "kg_prediction_layers.csv", row.names = FALSE)
  }

  rows <- list()
  for(la in lays){
    if(!is.null(allow) && !(la$display %in% allow) &&
       !(!is.null(la$base_display) && la$base_display %in% allow) &&
       !(.kgLayerKey(la$display) %in% tolower(allow))) next
    fv <- .kgFeatureVector(con, la); if(is.null(fv)) next      # read once per layer, used by every outcome

    for(outcome in names(spec)){
      p <- spec[[outcome]]
      if(p$is_cont){
        ynum <- p$ynum
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
        glab <- p$glab; pos_lab <- p$pos_lab; neg_lab <- p$neg_lab
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
  }
  if(length(rows) == 0) return("RES-NO")

  out <- do.call(rbind, rows)
  # rank by predictive strength: AUC distance from 0.5, or R^2
  ord <- ifelse(out$Metric=="AUC", abs(out$MetricValue-0.5), out$MetricValue)
  out <- out[order(-ord), ]
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  safe <- function(s) gsub("[^A-Za-z0-9]", "_", s)
  # a list of outcomes is named by its count — the joined ids would make a file name no file system accepts
  oname <- if(multi) paste0(length(spec), "_outcomes") else names(spec)[1]
  utils::write.csv(out, paste0("kg_prediction_", safe(var_id), "_", safe(oname), "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_prediction.csv", row.names = FALSE)

  if(exists("rcmd")){
    if(!file.exists("savedAnalysis")) dir.create("savedAnalysis")
    write(gsub("tool","local", rcmd), file = "savedAnalysis/Rhistory.R", append = TRUE)
  }
  paste0("RES-OK;", nrow(out))
}
