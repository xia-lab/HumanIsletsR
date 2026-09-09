# ==============================================================================
# kg_functions/kg_correlation.R  --  computed-only correlations + cross-omics concordance.
#   * F31  feature <-> feature correlation   ("does A track B?" / "what co-varies with X?")
#   * F32  phenotype <-> phenotype correlation ("how does secretion relate to BMI?")
#   * F19  cross-omics concordance (op#3)     (same gene across its gene modalities)
# Method = Pearson AND Spearman over the donors that have both; optional partial
# correlation when the user asks to control for covariates.
#
# Self-contained: all shared helpers (paths, resolver, correlation, metadata) live in
# kg_common.R (sourced first by the servlet). Nothing here calls another folder's code.
# Platform contract: string args; read data via the kg globals; write a results CSV;
# return "RES-OK;<n_rows>" | "RES-NO" | "RES-NO-RESOLVE" | "RES-NO-RESOLVER"
#   | "RES-NO-DB" | "RES-NO-META".
# ==============================================================================


# ==============================================================================
# kgFeatureCorr - feature <-> feature correlation (F31)
# ==============================================================================
kgFeatureCorr <- function(var_id_a, var_id_b,
                          var_type_a = "auto", var_type_b = "auto",
                          omics_filter = "all", covariates = "NA",
                          minN = "10", method = "both", mode = "tool"){
  library(RSQLite); library(DBI)
  minN <- suppressWarnings(as.integer(minN)); if(is.na(minN)) minN <- 10L
  method <- tolower(method)

  res <- .kgLoadFeatureResolver(); if(is.null(res)) return("RES-NO-RESOLVER")
  .kgSetPaths()
  db.path <- paste0(sqlite.path, "HI_omics_v2.sqlite")
  if(!file.exists(db.path)) return("RES-NO-DB")

  ra <- .kgResolveFeature(res, var_id_a, var_type_a)
  rb <- .kgResolveFeature(res, var_id_b, var_type_b)
  if(length(ra$layers) == 0 || length(rb$layers) == 0) return("RES-NO-RESOLVE")

  cov_names <- if(identical(covariates,"NA") || covariates == "") character(0) else trimws(strsplit(covariates, ",")[[1]])
  cov.df <- .kgLoadCovariates(cov_names)
  allow  <- if(identical(omics_filter,"all")) NULL else trimws(strsplit(omics_filter, "\\|")[[1]])

  con <- dbConnect(SQLite(), db.path); on.exit(dbDisconnect(con), add = TRUE)

  # same-modality pairing: same base type (gene<->gene) -> same table only (same platform,
  # same donors); different-type features (gene<->metabolite) -> all layer combinations.
  base_type <- function(t) if(t %in% c("gene","ensembl","entrez")) "gene" else t
  same_base <- base_type(ra$type) == base_type(rb$type)

  rows <- list()
  for(la in ra$layers){
    for(lb in rb$layers){
      if(same_base && !identical(la$table, lb$table)) next
      if(!is.null(allow) && !(la$display %in% allow) && !(lb$display %in% allow)) next
      va <- .kgFeatureVector(con, la); vb <- .kgFeatureVector(con, lb)
      if(is.null(va) || is.null(vb)) next
      used_cov <- ""
      if(!is.null(cov.df)){
        va <- .kgResidualize(va, cov.df); vb <- .kgResidualize(vb, cov.df)
        used_cov <- paste(intersect(names(cov.df), cov_names), collapse = ",")
      }
      cp <- .kgCorPair(va, vb, method)
      if(cp$n < minN) next
      rows[[length(rows)+1]] <- data.frame(
        FeatureA = var_id_a, LayerA = la$display, FeatureB = var_id_b, LayerB = lb$display, N = cp$n,
        Pearson_r = signif(cp$pear_r,3), Pearson_p = signif(cp$pear_p,3),
        Spearman_rho = signif(cp$spear_rho,3), Spearman_p = signif(cp$spear_p,3),
        Covariates = ifelse(used_cov=="", "none", used_cov), stringsAsFactors = FALSE)
    }
  }
  if(length(rows) == 0) return("RES-NO")

  out <- do.call(rbind, rows)
  ord <- ifelse(is.na(out$Spearman_rho), abs(out$Pearson_r), abs(out$Spearman_rho))
  out <- out[order(-ord), ]
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); safe <- function(s) gsub("[^A-Za-z0-9]", "_", s)
  utils::write.csv(out, paste0("kg_feature_corr_", safe(var_id_a), "_", safe(var_id_b), "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_feature_corr.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out))
}


# ==============================================================================
# kgPhenoCorr - phenotype <-> phenotype correlation (F32)
# ==============================================================================
# ⚠ THE MANY-PAIR PATH, kept beside the pair path rather than replacing it (2026-09-01). It reuses
# the SAME loader, the SAME residualizer and the SAME `.kgCorPair`, so a row produced here is
# computed identically to a row produced one pair at a time — only the number of round trips
# differs. Written as a helper so the single-pair branch above stays visibly untouched.
# ⚠ EVERY PAIR IS REPORTED, INCLUDING THE ONES BELOW `minN`. A member the cohort measures on too
# few shared donors is a real part of the answer to a family question ("13 of the 67 could not be
# tested"), and dropping it silently would let the reader read the survivors as the whole family.
# `N` carries the shared-donor count and the correlation columns are NA for those rows.
# ⚠ RES-NO ONLY WHEN NOTHING WAS TESTABLE — otherwise a single un-testable member would refuse a
# question the rest of the family answers.
.kgPhenoCorrMany <- function(m, cols_a, cols_b, covariates, minN, method){
  rid <- as.character(m[["record_id"]])
  cov_names <- if(identical(covariates, "NA") || covariates == "") character(0)
               else trimws(strsplit(covariates, ",")[[1]])
  cov_names <- cov_names[cov_names %in% names(m)]
  cov.df <- NULL
  if(length(cov_names) > 0){ cov.df <- .kgCoerce(m, cov_names); rownames(cov.df) <- rid }
  used_cov <- if(length(cov_names)) paste(cov_names, collapse = ",") else ""
  vec <- function(col){
    v <- suppressWarnings(as.numeric(m[[col]])); names(v) <- rid
    if(!is.null(cov.df)) v <- .kgResidualize(v, cov.df)
    v
  }
  cache <- list()
  getv <- function(col){
    if(is.null(cache[[col]])) cache[[col]] <<- vec(col)
    cache[[col]]
  }
  rows <- list(); n_ok <- 0L
  for(ca in cols_a) for(cb in cols_b){
    if(identical(ca, cb)) next                       # a phenotype against itself is not a finding
    cp <- .kgCorPair(getv(ca), getv(cb), method)
    ok <- isTRUE(cp$n >= minN)
    if(ok) n_ok <- n_ok + 1L
    rows[[length(rows) + 1L]] <- data.frame(
      PhenotypeA = ca, PhenotypeB = cb, N = cp$n,
      Pearson_r    = if(ok) signif(cp$pear_r, 3)    else NA_real_,
      Pearson_p    = if(ok) signif(cp$pear_p, 3)    else NA_real_,
      Spearman_rho = if(ok) signif(cp$spear_rho, 3) else NA_real_,
      Spearman_p   = if(ok) signif(cp$spear_p, 3)   else NA_real_,
      Covariates   = ifelse(used_cov == "", "none", used_cov),
      Tested       = ok, stringsAsFactors = FALSE)
  }
  if(!length(rows)) return("RES-NO-RESOLVE")
  if(n_ok == 0L)    return("RES-NO")
  out <- do.call(rbind, rows)
  # ⚠ THE MULTIPLE-TESTING CORRECTION BELONGS HERE, BESIDE THE TESTS IT CORRECTS (user ruling
  # 2026-09-01: "the correction should be in R"). This function issues ONE test per pair and is the
  # only place that knows how many were asked, so `P_family` is BH over exactly those pairs —
  # identical to `kg_corr.R`:255-259, the same rule and the same column name.
  #   prim <- Pearson_p when present, else Spearman_p        (kg_corr.R:255)
  #   P_family <- signif(p.adjust(prim, "BH"), 4)            (kg_corr.R:259)
  # ⚠ NO NEW @QueryParam, SO NO JAVA REBUILD. `KgPhenoCorr.java` passes a fixed argument list and
  # exposes no `fdrFamily`; emitting the column UNCONDITIONALLY adds it for every caller without
  # touching the servlet. kgCorr gates on `fdr_family` because its servlet offers the switch — this
  # one has no switch to read, and a single pair is its own family (BH over n=1 is the raw p), so
  # the column is correct in both shapes.
  # ⚠ WHY IT MATTERS: §L's `correlation` family counts a SIGNIFICANT result, and the caller was
  # left adjusting raw p-values itself because this column did not exist — the correction sat one
  # process away from the tests it belonged to.
  .prim <- if("Pearson_p" %in% names(out) && any(!is.na(out$Pearson_p))) out$Pearson_p
           else out$Spearman_p
  out$P_family <- signif(stats::p.adjust(.prim, method = "BH"), 4)
  out <- out[order(out$Pearson_p, na.last = TRUE), , drop = FALSE]
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); safe <- function(s) gsub("[^A-Za-z0-9]", "_", s)
  utils::write.csv(out, paste0("kg_pheno_corr_", safe(cols_a[1]), "_",
                               safe(cols_b[1]), "_many_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_pheno_corr.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out))
}


.kgPhenoCorrImpl <- function(pheno_a, pheno_b,
                        covariates = "NA", use_raw = "false",
                        minN = "10", method = "both", mode = "tool",
                        subset = "all"){
  # ⚠ `subset` IS NEW (2026-09-04) AND IS LAST ON PURPOSE, so every existing positional caller is
  # unaffected. op10 is the corpus's heaviest user of donor filters (F32 alone carries 404 filtered
  # questions) and accepted none of them, so every one had its computed op DROPPED upstream.
  minN <- suppressWarnings(as.integer(minN)); if(is.na(minN)) minN <- 10L
  method <- tolower(method)

  m <- .kgLoadMetaFrame(tolower(use_raw) %in% c("true","1")); if(is.null(m)) return("RES-NO-META")
  # the donor row filter, applied BEFORE any vector is built from the frame -- so x, y, the
  # covariates and n all describe the same restricted set. NULL means it could not be applied, and
  # that REFUSES: the unrestricted correlation is a different question from the one asked.
  m <- .kgDonorSubset(m, subset); if(is.null(m)) return("RES-NO-SUBSET")
  if(nrow(m) < minN) return("RES-NO")

  # resolve each token to a metadata column (exact or case-insensitive)
  col_lc <- tolower(names(m))
  resolve_pheno <- function(tok){
    if(tok %in% names(m)) return(tok)
    j <- match(.kgNorm(tok), col_lc); if(!is.na(j)) return(names(m)[j]); NA_character_
  }
  # ⚠ ONE SIDE MAY NAME MANY PHENOTYPES (2026-09-01) — the contract `kgCorr` has always had
  # ("a '|'-separated list | 'all' ... so one-vs-all stays vectorised instead of 16.7k cor.test
  # calls"), which this function never received. Without it a CONCEPT question has to be answered
  # by one HTTP call per member: MEASURED, `beta cell function` x HbA1c is 67 pairs and took 49.4 s
  # through the auto-fire path — 67 Rserve evals and 67 CSV round trips for arithmetic that is
  # milliseconds. As ONE call it is a single round trip.
  # ⚠ A LIST IS ALLOWED ON AT MOST ONE SIDE, the same structural guard `kgCorr` states: all-x-all
  # is a different question (a matrix), not a bigger version of this one.
  # ⚠ THE SINGLE-PAIR PATH IS UNCHANGED. One token on each side takes exactly the code below, so
  # every existing caller returns byte-identical rows; the loop simply runs once.
  # ⚠ AN UNRESOLVABLE MEMBER IS DROPPED, NOT FATAL. A concept's members are heterogeneous — some
  # are measured columns and some are not — and refusing the whole question because one member is
  # unmeasured would lose the answer for the rest. If NONE resolve, the old RES-NO-RESOLVE stands.
  .split_side <- function(tok){
    p <- trimws(strsplit(as.character(tok), "|", fixed = TRUE)[[1]])
    p[nzchar(p)]
  }
  toks_a <- .split_side(pheno_a); toks_b <- .split_side(pheno_b)
  if(length(toks_a) > 1 && length(toks_b) > 1) return("RES-TOO-BROAD")
  cols_a <- unique(stats::na.omit(vapply(toks_a, resolve_pheno, character(1))))
  cols_b <- unique(stats::na.omit(vapply(toks_b, resolve_pheno, character(1))))
  if(length(cols_a) == 0 || length(cols_b) == 0) return("RES-NO-RESOLVE")
  if(length(cols_a) > 1 || length(cols_b) > 1)
    return(.kgPhenoCorrMany(m, cols_a, cols_b, covariates, minN, method))

  ca <- cols_a[1]; cb <- cols_b[1]
  if(is.na(ca) || is.na(cb)) return("RES-NO-RESOLVE")

  cov_names <- if(identical(covariates,"NA") || covariates == "") character(0) else trimws(strsplit(covariates, ",")[[1]])
  cov_names <- cov_names[cov_names %in% names(m)]

  rid <- as.character(m[["record_id"]])
  x <- suppressWarnings(as.numeric(m[[ca]])); names(x) <- rid
  y <- suppressWarnings(as.numeric(m[[cb]])); names(y) <- rid

  used_cov <- ""
  if(length(cov_names) > 0){
    cov.df <- .kgCoerce(m, cov_names); rownames(cov.df) <- rid
    x <- .kgResidualize(x, cov.df); y <- .kgResidualize(y, cov.df)
    # report what was actually residualised on: a covariate left constant by the donor subset is
    # dropped inside `.kgResidualize` (see `.kgVaryingCovs`), so naming the requested list here
    # would claim an adjustment that never happened.
    cov_use  <- .kgVaryingCovs(cov.df, cov_names)
    used_cov <- paste(cov_use, collapse = ",")
  }

  cp <- .kgCorPair(x, y, method)
  if(cp$n < minN) return("RES-NO")

  out <- data.frame(PhenotypeA = ca, PhenotypeB = cb, N = cp$n,
    Pearson_r = signif(cp$pear_r,3), Pearson_p = signif(cp$pear_p,3),
    Spearman_rho = signif(cp$spear_rho,3), Spearman_p = signif(cp$spear_p,3),
    Covariates = ifelse(used_cov=="", "none", used_cov), stringsAsFactors = FALSE)
  # ⚠ THE SAME COLUMN ON THE SINGLE-PAIR PATH, so both shapes of this function answer with the same
  # schema and a caller never has to ask which one produced the file. One pair IS its own family:
  # BH over a single test returns the raw p unchanged, so this adds a column, never a correction
  # that was not asked for. Same rule as `.kgPhenoCorrMany` above and `kg_corr.R`:255-259.
  .prim <- if(!is.na(out$Pearson_p[1])) out$Pearson_p else out$Spearman_p
  out$P_family <- signif(stats::p.adjust(.prim, method = "BH"), 4)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); safe <- function(s) gsub("[^A-Za-z0-9]", "_", s)
  utils::write.csv(out, paste0("kg_pheno_corr_", safe(ca), "_", safe(cb), "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_pheno_corr.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out))
}


# ==============================================================================
# kgPhenoGroup - phenotype x GROUPING screen (frame F35)
# ------------------------------------------------------------------------------
# Which phenotypes differ across the LEVELS of a grouping variable (donor STATE
# final_cluster / diagnosis / donorsex)? feature x state IS precomputed, but
# phenotype x grouping is NOT, so this computes it live over metadata_sum_norm,
# REUSING kgPhenoCorr's loader (.kgLoadMetaFrame) and residualizer (.kgResidualize).
# Non-parametric by default: Wilcoxon (2 levels) / Kruskal-Wallis (>2); method='t'
# or 'anova' switches to parametric. CONTINUOUS phenotype columns only (.kgInferType).
# Covariates optional (residualized, default none). BH-FDR across the phenotypes.
# Writes kg_pheno_group.csv into the user folder and returns "RES-OK;<nrows>".
#   grouping   : column to group by ("final_cluster" | "diagnosis" | "donorsex")
#   levels     : "C0,C2" -> pairwise (2 groups) ; "" -> ALL non-blank levels (omnibus)
#   phenotypes : "all" -> every continuous non-grouping column ; else a named list
.kgPhenoGroupImpl <- function(grouping, levels = "", phenotypes = "all",
                         covariates = "NA", use_raw = "false",
                         minN = "10", method = "auto", mode = "tool",
                         subset = "all"){   # `subset` LAST: positional callers unaffected
  minN <- suppressWarnings(as.integer(minN)); if(is.na(minN)) minN <- 10L
  method <- tolower(method)

  m <- .kgLoadMetaFrame(tolower(use_raw) %in% c("true","1")); if(is.null(m)) return("RES-NO-META")

  # the donor row filter, applied BEFORE the grouping is read off the frame -- so the levels, the
  # per-group n and every test below describe the SAME restricted donor set. NULL means the filter
  # could not be applied, and that REFUSES: the unrestricted group comparison is a different
  # question from the one asked, and returning it would be a wider answer wearing the asked one's face.
  m <- .kgDonorSubset(m, subset); if(is.null(m)) return("RES-NO-SUBSET")

  col_lc <- tolower(names(m))
  resolve_col <- function(tok){
    if(tok %in% names(m)) return(tok)
    j <- match(.kgNorm(tok), col_lc); if(!is.na(j)) return(names(m)[j]); NA_character_
  }

  gcol <- resolve_col(grouping); if(is.na(gcol)) return("RES-NO-RESOLVE")
  rid  <- as.character(m[["record_id"]])
  grp  <- as.character(m[[gcol]])

  # restrict to requested levels (pairwise); empty = ALL non-blank levels (omnibus)
  want <- if(identical(levels,"") || is.null(levels)) character(0) else trimws(strsplit(levels, ",")[[1]])
  keep <- !is.na(grp) & !(grp %in% c("", "NA"))
  if(length(want) > 0) keep <- keep & (grp %in% want)
  grp <- grp[keep]; rid <- rid[keep]
  glev <- sort(unique(grp)); if(length(glev) < 2) return("RES-NO-GROUPS")
  if(length(grp) < minN) return("RES-NO")

  cov_names <- if(identical(covariates,"NA") || covariates == "") character(0) else trimws(strsplit(covariates, ",")[[1]])
  cov_names <- cov_names[cov_names %in% names(m) & cov_names != gcol]
  cov.df <- NULL
  if(length(cov_names) > 0){ cov.df <- .kgCoerce(m[keep, , drop = FALSE], cov_names); rownames(cov.df) <- rid }

  if(identical(phenotypes,"all") || phenotypes == ""){
    cand <- setdiff(names(m), c("record_id", gcol, cov_names))
    phen_cols <- cand[vapply(cand, function(cc) identical(.kgInferType(m[[cc]]), "cont"), logical(1))]
  } else {
    toks <- trimws(strsplit(phenotypes, ",")[[1]])
    phen_cols <- unique(stats::na.omit(vapply(toks, resolve_col, character(1))))
  }
  if(length(phen_cols) == 0) return("RES-NO-PHENO")

  two <- length(glev) == 2
  rows <- vector("list", 0)
  for(pc in phen_cols){
    y <- suppressWarnings(as.numeric(m[[pc]][keep])); names(y) <- rid
    if(!is.null(cov.df)) y <- .kgResidualize(y, cov.df)
    yy <- y[is.finite(y)]; gg <- grp[match(names(yy), rid)]
    tab <- table(gg)
    if(length(yy) < minN || sum(tab >= 2) < 2) next   # need >= 2 groups each with >= 2 obs
    p <- NA_real_; stat <- NA_real_; dir <- ""
    if(two){
      st <- if(method %in% c("t","ttest","parametric")) try(stats::t.test(yy ~ gg), silent = TRUE)
            else try(stats::wilcox.test(yy ~ gg, exact = FALSE), silent = TRUE)
      if(inherits(st, "try-error")) next
      p <- st$p.value; stat <- unname(st$statistic)
      med <- tapply(yy, gg, stats::median, na.rm = TRUE)
      dir <- if(med[[glev[2]]] >= med[[glev[1]]]) paste0("higher in ", glev[2]) else paste0("higher in ", glev[1])
    } else {
      st <- if(method %in% c("anova","parametric")) try(summary(stats::aov(yy ~ factor(gg)))[[1]][["Pr(>F)"]][1], silent = TRUE)
            else try(stats::kruskal.test(yy ~ factor(gg))$p.value, silent = TRUE)
      if(inherits(st, "try-error") || is.null(st) || is.na(st)) next
      p <- as.numeric(st); dir <- "differs across groups (non-directional)"
    }
    rows[[length(rows)+1]] <- data.frame(
      Phenotype = pc, N = length(yy), Groups = paste(glev, collapse = "/"),
      Statistic = signif(stat, 3), P_value = signif(p, 3), Direction = dir,
      stringsAsFactors = FALSE)
  }
  if(length(rows) == 0) return("RES-NO")

  out <- do.call(rbind, rows)
  out$Adjusted_p_value <- signif(stats::p.adjust(out$P_value, method = "BH"), 3)
  out <- out[order(out$P_value), ]
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); safe <- function(s) gsub("[^A-Za-z0-9]", "_", s)
  utils::write.csv(out, paste0("kg_pheno_group_", safe(gcol), "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_pheno_group.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out))
}


# ==============================================================================
# kgCrossOmicsConcordance - cross-omics concordance (F19 / op#3)
# ------------------------------------------------------------------------------
# Within-gene agreement across its gene-measuring modalities: correlate the SAME gene's
# per-donor values between every pair of its donor-level layers (RNA-seq, Nanostring,
# protein, pseudobulk). Only defined for a gene present in >=2 gene modalities.
# ==============================================================================
kgCrossOmicsConcordance <- function(var_id, var_type = "auto",
                                    covariates = "NA", minN = "10",
                                    method = "both", mode = "tool",
                                    subset = "all"){
  # `subset` is LAST so every existing positional caller is unaffected.
  library(RSQLite); library(DBI)
  minN <- suppressWarnings(as.integer(minN)); if(is.na(minN)) minN <- 10L
  method <- tolower(method)

  res <- .kgLoadFeatureResolver(); if(is.null(res)) return("RES-NO-RESOLVER")
  .kgSetPaths()
  db.path <- paste0(sqlite.path, "HI_omics_v2.sqlite")
  if(!file.exists(db.path)) return("RES-NO-DB")

  r <- .kgResolveFeature(res, var_id, var_type)
  if(!identical(r$type, "gene") || length(r$layers) < 2) return("RES-NO")

  cov_names <- if(identical(covariates,"NA") || covariates == "") character(0) else trimws(strsplit(covariates, ",")[[1]])
  cov.df <- .kgLoadCovariates(cov_names)

  con <- dbConnect(SQLite(), db.path); on.exit(dbDisconnect(con), add = TRUE)

  # the donor row filter. This op holds no metadata frame of its own -- it works on per-donor
  # vectors keyed by record_id -- so the filter is resolved once here and each layer vector is
  # narrowed to those donors. Restricting the VECTORS (not just the pairing) keeps every N, every
  # residualisation and every correlation on the same donor set.
  keep_ids <- NULL
  if(!identical(subset, "all") && nzchar(subset)){
    mfr <- .kgLoadMetaFrame(FALSE); if(is.null(mfr)) return("RES-NO-META")
    msub <- .kgDonorSubset(mfr, subset); if(is.null(msub)) return("RES-NO-SUBSET")
    keep_ids <- as.character(msub[["record_id"]])
    if(length(keep_ids) < minN) return("RES-NO")
  }

  vecs <- list()
  for(la in r$layers){
    v <- .kgFeatureVector(con, la)
    if(!is.null(v) && !is.null(keep_ids)) v <- v[names(v) %in% keep_ids]
    if(!is.null(v) && length(v) >= minN) vecs[[length(vecs)+1]] <- list(display = la$display, v = v)
  }
  if(length(vecs) < 2) return("RES-NO")

  rows <- list()
  for(i in 1:(length(vecs)-1)) for(j in (i+1):length(vecs)){
    va <- vecs[[i]]$v; vb <- vecs[[j]]$v; used_cov <- ""
    if(!is.null(cov.df)){
      va <- .kgResidualize(va, cov.df); vb <- .kgResidualize(vb, cov.df)
      used_cov <- paste(intersect(names(cov.df), cov_names), collapse = ",")
    }
    cp <- .kgCorPair(va, vb, method)
    if(cp$n < minN) next
    rows[[length(rows)+1]] <- data.frame(
      Feature = var_id, LayerA = vecs[[i]]$display, LayerB = vecs[[j]]$display, N = cp$n,
      Pearson_r = signif(cp$pear_r,3), Pearson_p = signif(cp$pear_p,3),
      Spearman_rho = signif(cp$spear_rho,3), Spearman_p = signif(cp$spear_p,3),
      Covariates = ifelse(used_cov=="", "none", used_cov), stringsAsFactors = FALSE)
  }
  if(length(rows) == 0) return("RES-NO")

  out <- do.call(rbind, rows)
  ord <- ifelse(is.na(out$Spearman_rho), abs(out$Pearson_r), abs(out$Spearman_rho))
  out <- out[order(-ord), ]
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); safe <- function(s) gsub("[^A-Za-z0-9]", "_", s)
  utils::write.csv(out, paste0("kg_crossomics_", safe(var_id), "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_crossomics.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out))
}

# ---- crash guard (.kgGuard, kg_common.R) -------------------------------------
# The servlet calls the PUBLIC name; the body now lives in the .*Impl above. A thrown
# R error comes back as "RES-ERR;<fn>;<message>" instead of a null -> empty HTTP body.
# Every RES-OK / RES-NO* / RES-TOO-* RETURN passes through untouched: tryCatch sees
# conditions, not return values, so "no data" and "crashed" stay distinct answers.
kgPhenoCorr <- function(...) .kgGuard("kgPhenoCorr", .kgPhenoCorrImpl, ...)
kgPhenoGroup <- function(...) .kgGuard("kgPhenoGroup", .kgPhenoGroupImpl, ...)
