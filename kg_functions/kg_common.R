# ==============================================================================
# kg_functions/kg_common.R  --  SELF-CONTAINED shared engine for the KG analysis funcs.
# ------------------------------------------------------------------------------
# RULE: the kg_functions folder is PURELY INDEPENDENT for ANALYSIS. Nothing here (or in
# any kg_* file) calls an analysis function from web_tool_functions (no DonorRegression,
# PatchseqSpearman, performGSEA/ORA). Every method is implemented here directly. kg always
# reads the current data; there is no version switch.
# THE ONE REUSE = filepath_utils.R (setPaths): the shared, relative, machine-detected path
# config -- kg reuses it rather than duplicating the paths (the servlet sources it first).
#
# Data (read directly, never through another folder's function):
#   omics       HI_omics_v2.sqlite  (proc_* only; feature x donor; joinable donors = R###)
#   phenotypes  display_data/metadata_sum_norm.csv (record_id x phenotype; values match
#               the precompute's proc_metadata, so an on-demand result is comparable)
#   resolver    display_kg/kg_feature_resolver.qs (token -> proc_* table + read key)
#   single-cell hdf5_V2/sc_<cell>_<glucose>.h5 + HI_tables.sqlite ephys_cell (op#12)
#   libraries   libraries/<lib>.rds ($sets = pathway -> entrez members, $term = names)
# ==============================================================================

# ---- paths: REUSE the shared filepath_utils.R setPaths() (the one reused piece) --------
# The servlet sources filepath_utils.R first, so setPaths() is defined; we never hardcode
# or duplicate the machine paths here.
.kgSetPaths <- function(){
  if(!exists("sqlite.path") || !exists("other.tables.path") || !exists("h5.v2.path")){
    if(exists("setPaths")) setPaths()
  }
  invisible()
}

.kgNorm <- function(x) tolower(trimws(as.character(x)))

# ---- crash guard: AN ERROR IS NOT "NO DATA" ----------------------------------
# Every kg entry function ANSWERS with a code: "RES-OK;<n>" when it computed something,
# and one of the RES-NO* / RES-TOO-* / RES-ONE-CLASS / RES-USE-PRECOMPUTE codes when the
# data cannot support the question. Those are RETURN VALUES, and `tryCatch` intercepts
# CONDITIONS only -- so every no-data answer travels through this guard untouched. The
# guard cannot turn a no-data answer into an error, or an error into a no-data answer.
#
# WHY IT EXISTS: `KgRCenter`'s twenty op wrappers all end
#     } catch (Exception rse) { System.out.println(rse); } return null;
# so an uncaught R error becomes a null -> a zero-length HTTP body -> a caller that reads
# "" with no error, which is indistinguishable from an honest empty result. That is how
# the F-statistic crash in `.kgLimmaFit` stayed invisible on every >=3-level omnibus: the
# analysis never ran, and the answer read as "the analysis found nothing". The guard gives
# a crash its OWN code -- "RES-ERR;<fn>;<message>" -- so the two can never be confused.
#
# It closes the OTHER silent-null route too: `RC.eval(...).asString()` throws whenever R
# hands back something that is not a string, and that lands in the same Java catch. A
# return that cannot be delivered as a character answer is reported as an error rather
# than vanishing. Anything that works today -- a character result of length >= 1 with a
# non-empty first element, which is what `asString()` reads -- is passed straight through.
#
# `;` is the field separator in every RES code, so it is stripped from the message.
.kgGuard <- function(.fn, .impl, ...){
  .clean <- function(s){
    s <- gsub("[\r\n\t;]+", " ", as.character(s)[1])
    s <- gsub("[[:space:]]+", " ", s)
    s <- gsub("^[[:space:]]+|[[:space:]]+$", "", s)
    if(is.na(s) || !nzchar(s)) s <- "unknown R error"
    if(nchar(s) > 300) s <- paste0(substr(s, 1, 300), " ...")
    s
  }
  .msg <- NULL
  out  <- tryCatch(.impl(...), error = function(e){ .msg <<- .clean(conditionMessage(e)); NULL })
  if(!is.null(.msg)) return(paste0("RES-ERR;", .fn, ";", .msg))
  if(is.character(out) && length(out) >= 1L && !is.na(out[1L]) && nzchar(out[1L])) return(out)
  paste0("RES-ERR;", .fn, ";returned no deliverable answer (",
         paste(class(out), collapse = "/"), ", length ", length(out), ")")
}

# ---- feature resolver (qs) ---------------------------------------------------
.kgLoadFeatureResolver <- function(){
  library(qs)
  .kgSetPaths()
  qs.path <- paste0(other.tables.path, "display_kg/kg_feature_resolver.qs")
  if(!file.exists(qs.path)) return(NULL)
  mt <- file.info(qs.path)$mtime
  if(is.null(.GlobalEnv$.kg_feat_cache) || !identical(.GlobalEnv$.kg_feat_mtime, mt)){
    .GlobalEnv$.kg_feat_cache <- qs::qread(qs.path)
    .GlobalEnv$.kg_feat_mtime <- mt
  }
  .GlobalEnv$.kg_feat_cache
}

# resolve a token -> list(type, layers[]{display, table, read_col, read_val}); type
# disambiguation (auto): gene > reaction > contaminant > metabolite.
.kgResolveFeature <- function(res, var_id, var_type = "auto"){
  q    <- .kgNorm(var_id)
  hits <- res[res$token == q, , drop = FALSE]
  if(nrow(hits) == 0) return(list(type = NA, layers = list()))
  want <- if(var_type %in% c("gene","ensembl","entrez")) "gene"
          else if(var_type %in% c("metabolite","contaminant","reaction")) var_type
          else NA_character_
  if(!is.na(want)){
    hits <- hits[hits$etype == want, , drop = FALSE]
  } else {
    for(et in c("gene","reaction","contaminant","metabolite"))
      if(any(hits$etype == et)){ hits <- hits[hits$etype == et, , drop = FALSE]; break }
  }
  if(nrow(hits) == 0) return(list(type = NA, layers = list()))
  layers <- lapply(seq_len(nrow(hits)), function(i)
    list(display = hits$layer[i], table = hits$tbl[i], read_col = hits$read_col[i], read_val = hits$read_val[i]))
  list(type = hits$etype[1], layers = layers)
}

# resolve MANY tokens in ONE pass -- same contract/priority as .kgResolveFeature, but for a
# vector. .kgResolveFeature scans+subsets the whole 290k-row resolver PER token (~9ms each,
# i.e. ~4.5s for a 500-feature query); this does a single %in% pass + one split, so the cost
# is flat in the number of tokens. (A persistent hashed index is NOT used on purpose: Rserve
# forks per request and re-sources, so .GlobalEnv would not survive and the index would be
# rebuilt every call -- a net loss for the common single-feature query.)
# Returns a named list keyed by NORMALISED token: token -> list(type, layers[]).
.kgResolveFeatures <- function(res, var_ids, var_type = "auto"){
  q <- unique(.kgNorm(var_ids)); if(length(q) == 0) return(list())
  hits <- res[res$token %in% q, , drop = FALSE]
  if(nrow(hits) == 0) return(list())
  by_tok <- split(seq_len(nrow(hits)), hits$token)
  want <- if(var_type %in% c("gene","ensembl","entrez")) "gene"
          else if(var_type %in% c("metabolite","contaminant","reaction")) var_type
          else NA_character_
  out <- list()
  for(tk in names(by_tok)){
    h <- hits[by_tok[[tk]], , drop = FALSE]
    if(!is.na(want)){
      h <- h[h$etype == want, , drop = FALSE]
    } else {
      for(et in c("gene","reaction","contaminant","metabolite"))
        if(any(h$etype == et)){ h <- h[h$etype == et, , drop = FALSE]; break }
    }
    if(nrow(h) == 0) next
    out[[tk]] <- list(type = h$etype[1], layers = lapply(seq_len(nrow(h)), function(i)
      list(display = h$layer[i], table = h$tbl[i], read_col = h$read_col[i], read_val = h$read_val[i])))
  }
  out
}

# a single feature's per-donor values (R### donors; bulk = all, pbrna = its R### subset).
.kgFeatureVector <- function(con, tuple){
  sql <- sprintf('SELECT * FROM %s WHERE "%s" = ?', tuple$table, tuple$read_col)
  row <- try(DBI::dbGetQuery(con, sql, params = list(tuple$read_val)), silent = TRUE)
  if(inherits(row, "try-error") || nrow(row) == 0) return(NULL)
  row <- row[1, , drop = FALSE]
  donor_cols <- grep("^R[0-9]+$", names(row), value = TRUE)
  if(length(donor_cols) == 0) return(NULL)
  v <- suppressWarnings(as.numeric(as.character(unlist(row[, donor_cols], use.names = FALSE))))
  names(v) <- donor_cols
  v
}

# a gene's symbol (the sc hdf5 / display key) from any of its resolver gene layers.
.kgGeneSymbol <- function(con, layers){
  for(la in layers){
    q <- sprintf('SELECT symbol FROM %s WHERE "%s" = ? LIMIT 1', la$table, la$read_col)
    s <- try(DBI::dbGetQuery(con, q, params = list(la$read_val)), silent = TRUE)
    if(!inherits(s, "try-error") && nrow(s) && !is.na(s$symbol[1]) && nzchar(s$symbol[1]))
      return(as.character(s$symbol[1]))
  }
  NA_character_
}

# ---- correlation helpers -----------------------------------------------------
.kgCorPair <- function(x, y, method){
  common <- intersect(names(x), names(y)); x <- x[common]; y <- y[common]
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  n <- length(x)
  res <- list(n = n, pear_r = NA, pear_p = NA, spear_rho = NA, spear_p = NA)
  if(n < 3) return(res)
  if(method %in% c("both","pearson")){
    p <- try(stats::cor.test(x, y, method = "pearson"), silent = TRUE)
    if(!inherits(p, "try-error")){ res$pear_r <- unname(p$estimate); res$pear_p <- p$p.value }
  }
  if(method %in% c("both","spearman")){
    s <- try(suppressWarnings(stats::cor.test(x, y, method = "spearman")), silent = TRUE)
    if(!inherits(s, "try-error")){ res$spear_rho <- unname(s$estimate); res$spear_p <- s$p.value }
  }
  res
}

.kgResidualize <- function(v, cov.df){
  if(is.null(cov.df) || ncol(cov.df) == 0) return(v)
  d <- data.frame(.y = v[rownames(cov.df)], cov.df, check.names = FALSE)
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if(nrow(d) < 3) return(setNames(numeric(0), character(0)))
  # ⚠ A COVARIATE LEFT CONSTANT BY A DONOR SUBSET makes `lm(.y ~ .)` a try-error, and the line
  # below turns that into an EMPTY vector -- so the pair silently disappears from the result
  # instead of being residualised. Drop the constant instead: it is collinear with the intercept,
  # so the residuals are identical. See `.kgVaryingCovs`. Rows are already complete-cased, so the
  # donor set does not move either way.
  keep <- .kgVaryingCovs(d, setdiff(names(d), ".y"))
  if(length(keep) == 0)                       # nothing left to adjust for: intercept-only
    return(setNames(d$.y - mean(d$.y), rownames(d)))   # == residuals(lm(.y ~ 1))
  d   <- d[, c(".y", keep), drop = FALSE]
  fit <- try(stats::lm(.y ~ ., data = d), silent = TRUE)
  if(inherits(fit, "try-error")) return(setNames(numeric(0), character(0)))
  setNames(stats::residuals(fit), rownames(d))
}

# ---- phenotypes (metadata_sum_norm) ------------------------------------------
# full metadata frame (character), keyed by record_id; cached per session by mtime.
.kgLoadMetaFrame <- function(use_raw = FALSE){
  .kgSetPaths()
  f <- paste0(other.tables.path, "display_data/", if(use_raw) "metadata_sum_raw.csv" else "metadata_sum_norm.csv")
  if(!file.exists(f)) return(NULL)
  utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
}

# coerce named metadata columns to numeric (where numeric) else factor.
.kgCoerce <- function(m, cols){
  d <- m[, cols, drop = FALSE]
  for(c in cols){
    num <- suppressWarnings(as.numeric(d[[c]]))
    nonblank <- !is.na(d[[c]]) & !(d[[c]] %in% c("", "NA"))
    d[[c]] <- if(all(is.na(num[nonblank])) && any(nonblank)) as.factor(d[[c]]) else num
  }
  d
}

# covariate frame (rownames = record_id) for partial correlation.
.kgLoadCovariates <- function(cov_names){
  if(length(cov_names) == 0) return(NULL)
  m <- .kgLoadMetaFrame(); if(is.null(m)) return(NULL)
  cov_names <- cov_names[cov_names %in% names(m)]; if(length(cov_names) == 0) return(NULL)
  d <- .kgCoerce(m, cov_names); rownames(d) <- as.character(m[["record_id"]]); d
}

# infer "cont" | "disc" from a metadata column's raw values.
.kgInferType <- function(vals){
  nb <- vals[!is.na(vals) & !(vals %in% c("", "NA"))]
  if(length(nb) == 0) return(NA_character_)
  num <- suppressWarnings(as.numeric(nb))
  if(sum(!is.na(num)) / length(nb) > 0.95 && length(unique(nb)) > 6) "cont" else "disc"
}

# ---- the OMICS DIFFERENTIAL (limma over all features of one proc_* table) -----
# Self-contained limma matching the project's model: covariate-adjusted, trend+robust
# eBayes, BH-adjusted. cat phenotype -> makeContrasts(contrast - ref) or omnibus "anova";
# cont phenotype -> the phenotype coefficient. Returns a data.frame for ALL features:
#   FeatureId (the table's id_col value), Symbol, Gene_ID (entrez if present),
#   Effect, EffectType (log2FC|coefficient), T_statistic, P_value, Adjusted_p_value.
# meta = coerced design frame (rownames record_id; cols = phenotype + covariates).
# load a proc_* layer ONCE into a numeric matrix (rows = features, cols = R### donors).
# Donor columns are stored REAL, so numeric columns are kept as-is; only genuinely
# non-numeric columns are coerced. (The old as.numeric(as.character()) round-trip on
# already-numeric columns cost ~20s per table -- pure waste.) Returns list(mat, info).
.kgLoadLayerMatrix <- function(con, tbl, id_col){
  ft <- try(DBI::dbReadTable(con, tbl), silent = TRUE)
  if(inherits(ft, "try-error")) return(NULL)
  ids <- as.character(ft[[id_col]])
  keep <- !is.na(ids) & nzchar(ids) & !duplicated(ids)
  ft <- ft[keep, , drop = FALSE]; ids <- ids[keep]
  sym    <- if("symbol" %in% names(ft)) as.character(ft$symbol) else ids
  entrez <- if("gene_id" %in% names(ft)) as.character(ft$gene_id) else rep(NA_character_, length(ids))
  donor_cols <- grep("^R[0-9]+$", names(ft), value = TRUE)
  if(length(donor_cols) == 0) return(NULL)
  sub <- ft[, donor_cols, drop = FALSE]
  notnum <- !vapply(sub, is.numeric, logical(1))
  if(any(notnum)) for(j in which(notnum)) sub[[j]] <- suppressWarnings(as.numeric(as.character(sub[[j]])))
  mat <- as.matrix(sub); rownames(mat) <- ids
  list(mat = mat, info = data.frame(FeatureId = ids, Symbol = sym, Gene_ID = entrez, stringsAsFactors = FALSE))
}

# limma fit for ONE phenotype on a PRELOADED layer matrix (from .kgLoadLayerMatrix). Split
# out so a screen over many phenotypes reads the layer once and only refits. Same model as
# before: trend+robust eBayes; disc -> makeContrasts(contrast - ref) or omnibus "anova";
# cont -> the phenotype coefficient. Returns the all-feature data.frame (see .kgOmicsDE).
.kgLimmaFit <- function(mat, info, meta, phenotype, ptype, ref = "", contrast = ""){
  library(limma)
  common <- intersect(colnames(mat), rownames(meta))
  if(length(common) < 10) return(NULL)
  mat <- mat[, common, drop = FALSE]; md <- meta[common, , drop = FALSE]

  # ⚠ ONE-VS-REST — `ref = "rest"` IS A RESERVED TOKEN, NOT A LEVEL (2026-08-31). F14's
  # one-vs-rest mode ("which features distinguish cluster C0?") has no precomputed contrast, and
  # no donor carries a `final_cluster` value of "rest", so the subset below kept ONLY the C0
  # donors, the factor collapsed to a single level and the call returned RES-NO. Verified live
  # before this: /kgassoc contrast=C0-rest -> RES-NO, while C1-C0 -> RES-OK;1.
  # ⚠ THE CONSTRUCTION IS THE ATLAS'S OWN, not a new one: `humanislets_atlas.R` builds exactly
  # this binary grouping (`ifelse(all_donors %in% cluster_donors, "cluster", "rest")`, then
  # `makeContrasts(cluster - rest)`). Recoding the column here means the existing design,
  # contrast and covariate machinery below runs UNCHANGED — `grp` becomes c(contrast, "rest")
  # and `makeContrasts("C0 - rest")` is an ordinary two-level fit.
  # ⚠ NA STAYS NA. A donor with no value for this variable is not part of "the rest" — it is
  # unmeasured, and `complete.cases` below must still drop it. Folding NA into the reference
  # group would silently inflate n and bias every effect.
  # ⚠ THE FLOOR IS THE ATLAS'S (>= 3 donors in the named group), reused rather than re-chosen.
  # ⚠ INERT FOR EVERY OTHER CALL: a pairwise contrast names two real levels and the omnibus
  # passes "anova", so `ref` is never "rest" and nothing here fires.
  if(ptype == "disc" && identical(ref, "rest")){
    .v <- as.character(md[[phenotype]])
    .measured <- !is.na(.v) & nzchar(.v) & .v != "NA"
    md[[phenotype]] <- ifelse(.measured, ifelse(.v == contrast, contrast, "rest"), NA_character_)
    if(sum(md[[phenotype]] == contrast, na.rm = TRUE) < 3) return(NULL)
  }

  if(ptype == "disc" && !identical(contrast, "anova")){
    kd <- rownames(md)[as.character(md[[phenotype]]) %in% c(ref, contrast)]
    md <- md[kd, , drop = FALSE]; mat <- mat[, kd, drop = FALSE]
  }
  cc <- stats::complete.cases(md); md <- md[cc, , drop = FALSE]; mat <- mat[, rownames(md), drop = FALSE]
  if(ncol(mat) < 10) return(NULL)

  fk <- rowSums(!is.na(mat)) >= 9
  mat <- mat[fk, , drop = FALSE]; info <- info[fk, , drop = FALSE]
  if(nrow(mat) == 0) return(NULL)

  # judged on `md` AFTER complete.cases and the contrast subset above -- the rows the fit sees.
  covs <- .kgVaryingCovs(md, setdiff(colnames(md), phenotype))
  if(ptype == "disc"){
    md[[phenotype]] <- factor(md[[phenotype]])
    grp <- levels(md[[phenotype]])
    if(length(covs)) md <- droplevels(md)
    form   <- if(length(covs)) paste0("~ 0 + ", phenotype, " + ", paste(covs, collapse = " + ")) else paste0("~ 0 + ", phenotype)
    design <- stats::model.matrix(stats::as.formula(form), data = md)
    colnames(design)[seq_along(grp)] <- grp
    if(identical(contrast, "anova")){
      cs <- grp[grp != ref]; args <- as.list(paste0(cs, " - ", ref))
    } else {
      args <- list(paste0(contrast, " - ", ref))
    }
    args[["levels"]] <- design
    cmat <- do.call(limma::makeContrasts, args)
    fit <- limma::lmFit(mat, design, trend = TRUE, robust = TRUE)
    fit <- limma::contrasts.fit(fit, cmat); fit <- limma::eBayes(fit)
    tt  <- limma::topTable(fit, number = Inf, sort.by = "none")
    eff <- if(identical(contrast, "anova")) rep(NA_real_, nrow(tt)) else tt$logFC
    eff.type <- "log2FC"
  } else {
    form   <- if(length(covs)) paste0("~ ", phenotype, " + ", paste(covs, collapse = " + ")) else paste0("~ ", phenotype)
    design <- stats::model.matrix(stats::as.formula(form), data = md)
    fit <- limma::lmFit(mat, design, trend = TRUE, robust = TRUE); fit <- limma::eBayes(fit)
    coef.nm <- if(phenotype %in% colnames(design)) phenotype else colnames(design)[2]
    tt  <- limma::topTable(fit, coef = coef.nm, number = Inf, sort.by = "none")
    eff <- tt$logFC; eff.type <- "coefficient"
  }
  # ⚠ AN OMNIBUS HAS NO t, AND READING ONE CRASHED THE WHOLE CALL (2026-08-31). With TWO OR MORE
  # contrasts `topTable` returns limma's F-test table, whose columns are the per-contrast log2FCs
  # plus AveExpr / F / P.Value / adj.P.Val — there is no `t`. So `tt$t` was NULL, `data.frame()`
  # saw a zero-length column against 7,822 features and threw
  # "arguments imply differing number of rows: 7822, 1, 0". The R error was uncaught, so Rserve
  # returned error code 127 and the servlet answered with an EMPTY BODY — indistinguishable from a
  # successful call that found nothing. TRACED on `collagenase` (levels Roche|Serva|VitaCyte ->
  # contrasts "Serva - Roche" and "VitaCyte - Roche").
  # ⚠ WHY IT LOOKED INTERMITTENT: the count that matters is the levels PRESENT among the donors in
  # that layer, not the levels the variable declares. `donorsex` (2 levels -> 1 contrast) always
  # worked, and `diagnosis` worked on a layer where only two of its levels are represented.
  # ⚠ `Effect` ALREADY HANDLED THIS — `rep(NA_real_, nrow(tt))` for the anova, because an omnibus
  # has no single signed effect either. `T_statistic` was the one column that still assumed a
  # pairwise fit. F is the statistic the omnibus actually produces, so it is what is reported.
  .tstat <- if(!is.null(tt$t)) tt$t else if(!is.null(tt$F)) tt$F else rep(NA_real_, nrow(tt))
  data.frame(FeatureId = rownames(mat), Symbol = info$Symbol, Gene_ID = info$Gene_ID,
             Effect = eff, EffectType = eff.type, T_statistic = .tstat,
             P_value = tt$P.Value, Adjusted_p_value = tt$adj.P.Val,
             N = ncol(mat), stringsAsFactors = FALSE)
}

# original signature (kgEnrichment / kgMultiOmics / kgAssociation): load the layer, fit one.
.kgOmicsDE <- function(con, tbl, id_col, meta, phenotype, ptype, ref = "", contrast = ""){
  L <- .kgLoadLayerMatrix(con, tbl, id_col); if(is.null(L)) return(NULL)
  .kgLimmaFit(L$mat, L$info, meta, phenotype, ptype, ref, contrast)
}

# ---- pathway library loader --------------------------------------------------
# Normalises every library to list(sets = <id -> members>, term = <display names>).
# GENE libraries are <name>.rds and already have that shape. METABOLITE libraries are .qs
# and do not, so they are converted here (structure taken from the deployed performGSEA,
# humanislets_statistics.R:33-44).
# ⚠ THE LIBRARY NAME IS NOT ALWAYS THE FILE NAME. `hsa_kegg` (metabolite KEGG) lives in
# kegg_hsa_met.qs -- there is ALSO a hsa_kegg.qs in the same folder, but that is the
# MUMMICHOG library (cpd.lib / cpd.tree / pathways, for untargeted LC-MS peaks) and loading
# it here would silently produce a meaningless test. The mapping below is explicit for
# exactly that reason.
.kgLoadLibrary <- function(library){
  .kgSetPaths()
  lp  <- paste0(other.tables.path, "libraries/")
  lib.l <- tolower(library)

  if(lib.l == "hsa_kegg"){                      # metabolite KEGG -- keyed by KEGG compound id
    f <- paste0(lp, "kegg_hsa_met.qs")          # NOT hsa_kegg.qs (that is mummichog)
    if(!file.exists(f)) return(NULL)
    l <- try(qs::qread(f), silent = TRUE); if(inherits(l, "try-error")) return(NULL)
    sets <- l$mset.list; if(is.null(sets)) return(NULL)
    nm   <- setNames(names(l$path.ids), l$path.ids)     # hsa00010 -> "Glycolysis..."
    tm   <- nm[names(sets)]; tm[is.na(tm)] <- names(sets)[is.na(tm)]
    return(list(sets = sets, term = unname(tm)))
  }

  # ⚠ kegg_hsa_met.qs is the ONLY correct metabolite library (user ruling 2026-07-20), which
  # is also what §P states -- metabolite sets are KEGG-only. hsaGem_metabolite.* is present in
  # the same folder and is deliberately NOT reachable from here.

  f <- paste0(lp, library, ".rds")              # gene libraries: kegg / reactome / go_bp / go_mf
  if(!file.exists(f)) return(NULL)
  l <- try(readRDS(f), silent = TRUE); if(inherits(l, "try-error") || is.null(l$sets)) return(NULL)
  list(sets = l$sets, term = l$term)
}

# ---- enrichment (fgsea) on a ranked differential -----------------------------
# ranks: named numeric; the NAME is whatever identifier the library is keyed on -- entrez
# for the gene libraries, kegg_id / gem_id for the metabolite ones (the caller supplies it).
# minimum set size, matching the deployed performGSEA (humanislets_statistics.R:274-276):
# 15 for the gene libraries, 3 for the metabolite one (KEGG metabolite sets are much smaller).
# ⚠ This is NOT cosmetic filtering -- fgsea's default is 1, so without it 1-2 member "sets"
# are tested, and because they enter the multiple-testing pool they shift `padj` for EVERY
# pathway, including the real ones.
.kgLibMinSize <- function(library) if(identical(tolower(library), "hsa_kegg")) 3L else 15L

# ⚠ `sym` MAPS THE LIBRARY KEY BACK TO A READABLE SYMBOL. fgsea's leadingEdge is expressed in
# whatever the sets are keyed on (entrez for the gene libraries), which is unreadable in a tooltip.
# The caller holds the Key -> Symbol map, so it is passed in rather than re-derived here.
.kgFgsea <- function(ranks, library, fdr = 0.05, sym = NULL){
  library(fgsea); .kgSetPaths()
  # ⚠ fgsea's p-values come from a permutation procedure, so WITHOUT A SEED the same question
  # on the same data returns slightly different numbers each run. Every deployed GSEA path
  # seeds (performGSEA :20, endotypeOmicsPathway :2933). Same value, so kg and the web tool
  # agree run for run.
  set.seed(42)
  lib <- .kgLoadLibrary(library)
  if(is.null(lib)) return(NULL)
  ranks <- ranks[!is.na(ranks) & is.finite(ranks)]
  if(length(ranks) < 10) return(NULL)
  res <- try(fgsea::fgsea(pathways = lib$sets, stats = ranks,
                          minSize = .kgLibMinSize(library)), silent = TRUE)
  if(inherits(res, "try-error") || is.null(res) || nrow(res) == 0) return(NULL)
  res <- as.data.frame(res)
  nm  <- if(!is.null(lib$term)) setNames(lib$term, names(lib$sets)) else NULL
  # ⚠ THE LEADING EDGE IS THE ANSWER TO "WHICH GENES DROVE THIS PATHWAY" and fgsea already
  # computes it — it was simply never carried out of this function, so every consumer could name a
  # pathway but never the features behind it. Collapsed to a string because a list-column cannot be
  # written to CSV; the FULL size is reported alongside so a truncated list can never read as complete.
  le <- res$leadingEdge
  le.txt <- vapply(seq_len(nrow(res)), function(i){
    g <- le[[i]]
    if(is.null(g) || !length(g)) return("")
    if(!is.null(sym)) { mapped <- unname(sym[as.character(g)]); g <- ifelse(is.na(mapped), g, mapped) }
    paste(utils::head(g, 25), collapse = ";")
  }, character(1))
  le.n <- vapply(le, function(g) if(is.null(g)) 0L else length(g), integer(1))
  data.frame(
    Set_Name = if(!is.null(nm)) ifelse(is.na(nm[res$pathway]), res$pathway, nm[res$pathway]) else res$pathway,
    Set_ID = res$pathway, P_value = signif(res$pval, 3), Adjusted_p_value = signif(res$padj, 3),
    Normalized_ES = signif(res$NES, 3), Set_Size = res$size,
    Leading_Edge = le.txt, Leading_Edge_N = le.n, stringsAsFactors = FALSE)
}

# over-representation (fisher) of a significant entrez set against a library.
.kgOra <- function(sig_entrez, universe_entrez, library, fdr = 0.05){
  library(fgsea); .kgSetPaths()
  lib <- .kgLoadLibrary(library)
  if(is.null(lib)) return(NULL)
  sig_entrez <- unique(as.character(sig_entrez)); universe_entrez <- unique(as.character(universe_entrez))
  if(length(sig_entrez) < 1) return(NULL)
  # no seed needed here: fora is a hypergeometric test and is deterministic. minSize applies
  # for the same reason it does in GSEA -- it changes the multiple-testing pool.
  res <- try(fgsea::fora(pathways = lib$sets, genes = sig_entrez, universe = universe_entrez,
                         minSize = .kgLibMinSize(library)), silent = TRUE)
  if(inherits(res, "try-error") || is.null(res) || nrow(res) == 0) return(NULL)
  res <- as.data.frame(res)
  nm  <- if(!is.null(lib$term)) setNames(lib$term, names(lib$sets)) else NULL
  data.frame(
    Set_Name = if(!is.null(nm)) ifelse(is.na(nm[res$pathway]), res$pathway, nm[res$pathway]) else res$pathway,
    Set_ID = res$pathway, P_value = signif(res$pval, 3), Adjusted_p_value = signif(res$padj, 3),
    Overlap = res$overlap, Set_Size = res$size, stringsAsFactors = FALSE)
}

# ==============================================================================
# kgAssoc support: phenotype catalogue, targeted multi-feature read, per-cell lm.
# ==============================================================================

# phenotype catalogue (display_interface/proc_variable_summary_v2.csv): each row is a leaf
# column with group_id / group / type. Lets a phenotype token be a LEAF or a GROUP (which
# expands to its leaves). Cached per session by mtime.
.kgLoadVarSummary <- function(){
  .kgSetPaths()
  f <- paste0(other.tables.path, "display_interface/proc_variable_summary_v2.csv")
  if(!file.exists(f)) return(NULL)
  mt <- file.info(f)$mtime
  if(is.null(.GlobalEnv$.kg_vsum_cache) || !identical(.GlobalEnv$.kg_vsum_mtime, mt)){
    .GlobalEnv$.kg_vsum_cache <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    .GlobalEnv$.kg_vsum_mtime <- mt
  }
  .GlobalEnv$.kg_vsum_cache
}

# expand phenotype tokens -> leaf metadata columns present in meta. Each token matches, in
# order: exact meta column; var-summary group_id; var-summary group display; normalized meta
# column. Returns unique columns in request order (a parent like "insulin secretion" ->
# all its children).
.kgExpandPhenotypes <- function(tokens, meta_names){
  vs <- .kgLoadVarSummary(); out <- character(0)
  for(tk in tokens){
    tk <- trimws(tk); if(!nzchar(tk)) next
    if(tk %in% meta_names){ out <- c(out, tk); next }
    matched <- character(0)
    if(!is.null(vs)){
      gi <- vs$column[.kgNorm(vs$group_id) == .kgNorm(tk)]
      gp <- vs$column[.kgNorm(vs$group)    == .kgNorm(tk)]
      matched <- unique(c(gi, gp))
    }
    if(length(matched) == 0){
      j <- which(.kgNorm(meta_names) == .kgNorm(tk)); if(length(j)) matched <- meta_names[j[1]]
    }
    out <- c(out, matched[matched %in% meta_names])
  }
  unique(out)
}

# "all" phenotypes = every catalogue leaf present in meta (fallback: all meta cols minus id).
.kgAllPhenotypes <- function(meta_names){
  vs <- .kgLoadVarSummary()
  if(!is.null(vs)) return(intersect(vs$column, meta_names))
  setdiff(meta_names, "record_id")
}

# a phenotype column's type: catalogue "type" (cont|disc) first, else inferred from values.
.kgPhenoType <- function(m, col){
  vs <- .kgLoadVarSummary()
  if(!is.null(vs)){ t <- vs$type[.kgNorm(vs$column) == .kgNorm(col)]; if(length(t) && t[1] %in% c("cont","disc")) return(t[1]) }
  .kgInferType(m[[col]])
}

# read NAMED features' rows from a layer in ONE query (features x donor matrix, coercion-
# correct). read_vals = the resolver read_val per feature. Returns list(mat, id, symbol).
.kgFeatureRows <- function(con, tbl, id_col, read_vals){
  read_vals <- unique(as.character(read_vals)); if(length(read_vals) == 0) return(NULL)
  ph  <- paste(rep("?", length(read_vals)), collapse = ",")
  sql <- sprintf('SELECT * FROM %s WHERE "%s" IN (%s)', tbl, id_col, ph)
  ft  <- try(DBI::dbGetQuery(con, sql, params = as.list(read_vals)), silent = TRUE)
  if(inherits(ft, "try-error") || nrow(ft) == 0) return(NULL)
  donor_cols <- grep("^R[0-9]+$", names(ft), value = TRUE); if(length(donor_cols) == 0) return(NULL)
  sub <- ft[, donor_cols, drop = FALSE]
  notnum <- !vapply(sub, is.numeric, logical(1))
  if(any(notnum)) for(j in which(notnum)) sub[[j]] <- suppressWarnings(as.numeric(as.character(sub[[j]])))
  mat <- as.matrix(sub); rownames(mat) <- as.character(ft[[id_col]])
  list(mat = mat, id = as.character(ft[[id_col]]),
       symbol = if("symbol" %in% names(ft)) as.character(ft$symbol) else as.character(ft[[id_col]]))
}

# ---- KENDALL (contaminant x phenotype) ---------------------------------------
# METHOD OF RECORD: DonorRegression's proc_contaminants branch
# (web_tool_functions/humanislets_statistics.R:1878-1960). Reproduced here rather than
# called, per this folder's independence rule. Its rationale, verbatim from :1871-1874 --
# "Small sample sizes and non-normal abundance distributions make limma's moderated-t
# inappropriate for proc_contaminants."
#
#   continuous primary      -> Kendall's tau on the numeric phenotype
#   discrete 2-level/named  -> 0/1 encoding (ref = 0, contrast = 1), then tau
#   discrete >=3 + omnibus  -> Kruskal-Wallis; NO per-feature effect size exists
#   covariates present      -> PARTIAL Kendall via ppcor::pcor.test(method = "kendall")
#   covariates absent       -> cor.test(method = "kendall", exact = FALSE)
#
# ⚠ THE EFFECT IS TAU, NEVER A LOG FOLD CHANGE (:1930-1933), so EffectType is "tau" and the
# omnibus is "kruskalH". Reporting a tau under an EffectType that means something else is the
# specific mistake this naming prevents.
# ⚠ SCOPE IS EXACTLY contaminant x phenotype -- a collaborator decision that is OWNED, not a
# rule derived from the data. Do NOT generalise it to "small-n or tie-heavy layers": that
# would silently swallow metabolomics and methylation, which the decision never covered.
# ⚠ n floor: the caller's `minN` (default 10). The deployed branch hardcodes 5; kgAssoc has
# always exposed minN, so it is respected here rather than overridden. Flagged, not decided.
.kgKendallCell <- function(y, ph, covdf, ptype, contrast = "", ref = "", minN = 10){
  don <- intersect(names(y), names(ph))
  if(!is.null(covdf) && ncol(covdf)) don <- intersect(don, rownames(covdf))
  if(length(don) < minN) return(NULL)
  yv <- suppressWarnings(as.numeric(y[don]))

  if(ptype == "cont"){
    x <- suppressWarnings(as.numeric(as.character(ph[don])))
  } else {
    pc    <- as.character(ph[don])
    valid <- !is.na(pc) & pc != "" & pc != "NA"
    lv    <- unique(pc[valid])
    named <- nzchar(contrast) && !identical(contrast, "anova")
    if(named || length(lv) == 2){
      if(named){ a <- contrast; b <- ref } else { a <- lv[1]; b <- lv[2] }
      x <- ifelse(pc == a, 1L, ifelse(pc == b, 0L, NA_integer_))
    } else {
      # >=3 levels, no named contrast -> Kruskal-Wallis. It yields a chi-squared statistic
      # and a p-value only; there is NO per-feature coefficient, so Effect stays NA.
      d <- data.frame(y = yv, g = pc, stringsAsFactors = FALSE)
      d <- d[is.finite(d$y) & valid, , drop = FALSE]
      if(nrow(d) < minN || length(unique(d$g)) < 2) return(NULL)
      kw <- try(stats::kruskal.test(d$y ~ as.factor(d$g)), silent = TRUE)
      if(inherits(kw, "try-error")) return(NULL)
      return(list(N = nrow(d), Effect = NA_real_, EffectType = "kruskalH",
                  t = unname(kw$statistic), P = kw$p.value,
                  Direction = "NA", n_levels = length(unique(d$g))))
    }
  }

  ok <- is.finite(yv) & !is.na(x)
  cm <- NULL
  if(!is.null(covdf) && ncol(covdf)){
    cm <- covdf[don, , drop = FALSE]
    ok <- ok & stats::complete.cases(cm)
  }
  if(sum(ok) < minN || length(unique(x[ok])) < 2) return(NULL)

  use.pcor <- !is.null(cm) && requireNamespace("ppcor", quietly = TRUE)
  r <- try(if(use.pcor) ppcor::pcor.test(yv[ok], x[ok], cm[ok, , drop = FALSE], method = "kendall")
           else         stats::cor.test(yv[ok], x[ok], method = "kendall", exact = FALSE),
           silent = TRUE)
  if(inherits(r, "try-error")) return(NULL)
  tau <- unname(r$estimate); pv <- unname(r$p.value)
  list(N = sum(ok), Effect = tau, EffectType = "tau", t = tau, P = pv,
       Direction = if(is.na(tau)) "NA" else ifelse(tau > 0, "up", "down"),
       n_levels = if(ptype == "cont") NA_integer_ else 2L)
}

# Kendall across a WHOLE preloaded layer -- the screen counterpart of .kgLimmaFit, and it
# returns the SAME columns so the caller only chooses which fit to run. BH across the
# features tested here (as the deployed branch does at :1937).
.kgKendallFit <- function(mat, info, meta, phenotype, ptype, ref = "", contrast = "", minN = 10){
  common <- intersect(colnames(mat), rownames(meta))
  if(length(common) < minN) return(NULL)
  mat <- mat[, common, drop = FALSE]; md <- meta[common, , drop = FALSE]
  ph  <- setNames(md[[phenotype]], rownames(md))
  covn  <- setdiff(colnames(md), phenotype)
  covdf <- if(length(covn)) md[, covn, drop = FALSE] else NULL

  res <- lapply(seq_len(nrow(mat)), function(i){
    y <- mat[i, ]; names(y) <- colnames(mat)
    cell <- .kgKendallCell(y, ph, covdf, ptype, contrast, ref, minN)
    if(is.null(cell)) return(NULL)
    data.frame(FeatureId = rownames(mat)[i], Symbol = info$Symbol[i], Gene_ID = info$Gene_ID[i],
               Effect = cell$Effect, EffectType = cell$EffectType, T_statistic = cell$t,
               P_value = cell$P, N = cell$N, stringsAsFactors = FALSE)
  })
  res <- res[!vapply(res, is.null, logical(1))]
  if(length(res) == 0) return(NULL)
  out <- do.call(rbind, res)
  out$Adjusted_p_value <- stats::p.adjust(out$P_value, method = "BH")
  out[, c("FeatureId","Symbol","Gene_ID","Effect","EffectType","T_statistic",
          "P_value","Adjusted_p_value","N")]
}

# TARGETED single-feature association: lm(value ~ phenotype + covariates).
#  cont            -> phenotype slope (coefficient); EffectType "slope".
#  disc 2-level or named "A-B" -> covariate-adjusted mean difference; EffectType "meanDiff".
#  disc >=3, no contrast       -> covariate-adjusted omnibus ANOVA F; EffectType "anovaF".
# y = named feature vector (donor -> value); ph = named phenotype (donor -> value, coerced);
# covdf = donor x covariate frame or NULL. Returns a one-row list or NULL.
.kgLmCell <- function(y, ph, covdf, ptype, contrast = "", ref = "", minN = 10){
  don <- names(y); don <- intersect(don, names(ph))
  if(!is.null(covdf) && ncol(covdf)) don <- intersect(don, rownames(covdf))
  if(length(don) < minN) return(NULL)
  d <- data.frame(.y = as.numeric(y[don]), .ph = ph[don], stringsAsFactors = FALSE)
  if(!is.null(covdf) && ncol(covdf)) d <- cbind(d, covdf[don, , drop = FALSE])
  covn   <- setdiff(names(d), c(".y", ".ph"))
  # ⚠ BUILT AT EACH FIT, NOT ONCE HERE. Every branch below narrows `d` again (complete.cases, and
  # the two-level contrast subset), and a donor subset can leave a covariate constant on those rows
  # -- see `.kgVaryingCovs`. A constant covariate here does not crash, it makes lm() a try-error and
  # returns NULL for the feature, so the whole screen comes back RES-NO with nothing to explain it.
  .rhs <- function(dd){
    cn <- .kgVaryingCovs(dd, covn)
    if(length(cn)) paste0(" + ", paste(sprintf("`%s`", cn), collapse = " + ")) else ""
  }

  if(ptype == "cont"){
    d$.ph <- suppressWarnings(as.numeric(as.character(d$.ph)))
    d <- d[stats::complete.cases(d), , drop = FALSE]; if(nrow(d) < minN) return(NULL)
    fit <- try(stats::lm(stats::as.formula(paste0(".y ~ .ph", .rhs(d))), data = d), silent = TRUE)
    if(inherits(fit, "try-error")) return(NULL)
    co <- summary(fit)$coefficients; if(!(".ph" %in% rownames(co))) return(NULL)
    eff <- co[".ph","Estimate"]
    return(list(N = nrow(d), Effect = eff, EffectType = "slope",
                t = co[".ph","t value"], P = co[".ph","Pr(>|t|)"],
                Direction = ifelse(eff > 0, "up", "down"), n_levels = NA_integer_))
  }
  # discrete
  d$.ph <- as.character(d$.ph)
  d <- d[!is.na(d$.ph) & d$.ph != "" & d$.ph != "NA", , drop = FALSE]
  # ⚠ ONE-VS-REST — `ref = "rest"` IS A RESERVED TOKEN, NOT A LEVEL. Same rule and same reason as
  # `.kgLimmaFit` above; BOTH are needed because `kg_assoc.R` has two paths and they use different
  # fitters — a NAMED feature takes the targeted `.kgLmCell` path (this one) while a screen takes
  # `.kgLimmaFit`. Fixing only the screen left /kgassoc?features=TCF7L2&contrast=C0-rest returning
  # RES-NO, which is exactly how this was found.
  # ⚠ RECODED BEFORE THE SUBSET BELOW, so `%in% c(a, b)` keeps every donor and the ordinary
  # two-level fit runs unchanged; without it the subset kept only the named cluster and
  # `nlevels(droplevels(...)) < 2` returned NULL.
  # ⚠ NA IS ALREADY GONE (filtered on the line above), so "the rest" can only be measured donors.
  # ⚠ THE FLOOR IS THE ATLAS'S (>= 3 in the named group), the same constant used in `.kgLimmaFit`.
  # ⚠ INERT OTHERWISE: fires only when the caller passes the literal ref "rest", which no pairwise
  # or omnibus call does.
  if(identical(ref, "rest") && nzchar(contrast) && !identical(contrast, "anova")){
    d$.ph <- ifelse(d$.ph == contrast, contrast, "rest")
    if(sum(d$.ph == contrast) < 3) return(NULL)
  }
  lv    <- unique(d$.ph)
  named <- nzchar(contrast) && !identical(contrast, "anova")
  if(named || length(lv) == 2){
    if(named){ a <- contrast; b <- ref } else { a <- lv[1]; b <- lv[2] }
    d <- d[d$.ph %in% c(a, b), , drop = FALSE]
    d$.ph <- factor(d$.ph, levels = c(b, a))
    d <- d[stats::complete.cases(d), , drop = FALSE]
    if(nrow(d) < minN || nlevels(droplevels(d$.ph)) < 2) return(NULL)
    fit <- try(stats::lm(stats::as.formula(paste0(".y ~ .ph", .rhs(d))), data = d), silent = TRUE)
    if(inherits(fit, "try-error")) return(NULL)
    co <- summary(fit)$coefficients; rn <- paste0(".ph", a)
    if(!(rn %in% rownames(co))) return(NULL)
    eff <- co[rn,"Estimate"]
    return(list(N = nrow(d), Effect = eff, EffectType = "meanDiff",
                t = co[rn,"t value"], P = co[rn,"Pr(>|t|)"],
                Direction = ifelse(eff > 0, "up", "down"), n_levels = 2L))
  }
  d$.ph <- factor(d$.ph); d <- d[stats::complete.cases(d), , drop = FALSE]
  d <- droplevels(d); if(nrow(d) < minN || nlevels(d$.ph) < 2) return(NULL)
  # ONE covariate set for both, so anova(red, full) stays nested.
  .rc  <- .rhs(d)
  full <- try(stats::lm(stats::as.formula(paste0(".y ~ .ph", .rc)), data = d), silent = TRUE)
  red  <- try(stats::lm(stats::as.formula(paste0(".y ~ 1", .rc)),   data = d), silent = TRUE)
  if(inherits(full, "try-error") || inherits(red, "try-error")) return(NULL)
  an <- try(stats::anova(red, full), silent = TRUE); if(inherits(an, "try-error")) return(NULL)
  list(N = nrow(d), Effect = NA_real_, EffectType = "anovaF",
       t = an$F[2], P = an$`Pr(>F)`[2], Direction = "NA", n_levels = nlevels(d$.ph))
}

# ------------------------------------------------------------------------------------------------
# .kgVaryingCovs -- the covariates that STILL VARY in the frame actually being modelled.
# ------------------------------------------------------------------------------------------------
# ⚠ A DONOR SUBSET MAKES ITS OWN COLUMN CONSTANT. "Within male donors, adjusted for sex" leaves
# `donorsex` with one value, and `donorsex` is in `.kgMoDefaultCovs` -- so this is the NORMAL case
# for a subsetted call, not an edge case. Keeping a constant covariate is not conservative, it is
# UNDEFINED: model.matrix() errors on a single-level factor ("contrasts can be applied only to
# factors with 2 or more levels"), and at kg_common.R:303/319 and kg_corr.R:45 that error is
# UNCAUGHT -- the caller gets an empty body, which is how this was found (op7 returned "" for
# donorsex=Male while a bad column correctly returned RES-NO-SUBSET). In `.kgLmCell` the same
# error is caught by try() and the feature is dropped instead, so a subsetted screen returns
# RES-NO rather than crashing -- silent, and worse.
# ⚠ DROPPING IS AN IDENTITY, NOT A MODELLING CHOICE. A column with one distinct value is collinear
# with the intercept (`~ .`) or with the phenotype dummies (`~ 0 + pheno`); it carries no
# information and estimates no coefficient. Removing it moves NO effect, t or p below. This is
# why it is safe to apply unconditionally rather than only when a subset was passed.
# ⚠ NO CHANGE WITHOUT A SUBSET. A default covariate can only be constant across the whole cohort
# if the data itself is constant in it -- which crashes or silently drops TODAY. Nothing that
# currently returns a correct answer takes a different path.
# ⚠ CALLED AT THE FIT, NOT AT SELECTION: complete.cases() and a contrast subset both narrow the
# frame further, so constancy is judged on the rows the model actually sees.
.kgVaryingCovs <- function(d, covs){
  covs <- covs[covs %in% names(d)]
  if(!length(covs)) return(covs)
  keep <- vapply(covs, function(cn){
    v <- d[[cn]]
    if(is.factor(v)) v <- as.character(v)
    v <- v[!is.na(v)]
    v <- v[!(as.character(v) %in% c("", "NA"))]
    length(unique(v)) >= 2
  }, logical(1))
  covs[keep]
}

# ------------------------------------------------------------------------------------------------
# .kgDonorSubset -- THE donor row filter, used by every op that accepts one.
# ------------------------------------------------------------------------------------------------
# `subset` is `column<op>value`, op in = > < >= <=, e.g. "donorsex=Male", "bodymassindex>30".
# Returns the filtered frame, or NULL meaning "this filter could not be applied".
#
# -- ONE IMPLEMENTATION ON PURPOSE. kg_association.R had its own inline parser and kg_correlation.R
#    had none; adding a second copy is how the two drift into disagreeing about what a filter means.
# -- NULL MEANS REFUSE, NEVER "RUN UNRESTRICTED". The old inline version fell through when the
#    column was absent, returning the whole-cohort analysis as the answer to a restricted question --
#    a wider answer presented as the asked one, undetectable downstream. Callers must return
#    RES-NO-SUBSET on NULL.
# -- COMPARISONS CARRY THE CALLER'S OWN THRESHOLD ("donors over 60", "a BMI over 30"); no cut-off is
#    invented here.
.kgDonorSubset <- function(m, subset){
  if(is.null(subset) || identical(subset, "all") || !nzchar(subset)) return(m)
  mm <- regmatches(subset, regexec("^\\s*([^><=]+?)\\s*(>=|<=|=|>|<)\\s*(.+?)\\s*$", subset))[[1]]
  if(length(mm) != 4) return(NULL)
  scol <- trimws(mm[2]); sop <- mm[3]; sval <- trimws(mm[4])
  if(!(scol %in% names(m))) return(NULL)
  if(identical(sop, "=")){
    keep <- as.character(m[[scol]]) == sval
  } else {
    lhs <- suppressWarnings(as.numeric(as.character(m[[scol]])))
    rhs <- suppressWarnings(as.numeric(sval))
    if(all(is.na(lhs)) || is.na(rhs)) return(NULL)
    keep <- switch(sop, ">" = lhs > rhs, "<" = lhs < rhs,
                        ">=" = lhs >= rhs, "<=" = lhs <= rhs, rep(FALSE, length(lhs)))
  }
  keep[is.na(keep)] <- FALSE
  m[keep, , drop = FALSE]
}
