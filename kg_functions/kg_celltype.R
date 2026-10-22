# ==============================================================================
# kg_functions/kg_celltype.R  --  which islet cell types express a gene (or each gene of a set)
#
# Source: the pooled islet single-cell atlas the website's gene page already reads -- hdf5/ma_singlecell.h5,
# web_tool_functions/humanislets_figures.R `plotExpressionByCell` (same file, same fields: meta/genes/symbol with the
# Entrez fallback, meta/cells/celltype, data/norm_expression). One HDF5 row per gene, so each gene is one small read.
#
# ⚠ NOTHING HERE IS FIXED TO A QUESTION (user 2026-10-19: "need not hardcoded for the question, they can ask any other
# question which need cell type"). The genes come from the caller; the cell types are the atlas's OWN labels -- every
# type it holds by default, or the subset the caller names. A type with fewer than `minCells` cells is not summarised,
# and `kg_celltype_types.csv` lists every type with its cell count so the page can say which were skipped.
#
# features   comma list of gene symbols (or Entrez ids), any number
# celltypes  "all" (DEFAULT) | comma list of the atlas's own cell-type labels (case-insensitive)
# minCells   cells a type needs to be summarised (DEFAULT 20)
#
# Writes kg_celltype_expr.csv -- one row per (gene x cell type):
#   feature, symbol, celltype, n_cells, mean_expression, pct_expressing, mean_in_expressing
#   (mean_expression over all cells of the type, as `plotExpressionByCell`; pct_expressing = % cells with value > 0;
#    mean_in_expressing = mean over those cells only)
# and kg_celltype_types.csv -- celltype, n_cells, summarised.
# Returns "RES-OK;<n_rows>" | "RES-NO-DB" | "RES-NO-FEATURE" | "RES-NO-RESOLVE" | "RES-NO-CELLTYPE".
# ==============================================================================

kgCellTypeExpr <- function(...) .kgGuard("kgCellTypeExpr", .kgCellTypeExprImpl, ...)

.kgCellTypeExprImpl <- function(features = "", celltypes = "all", minCells = "20", mode = "tool"){
  suppressPackageStartupMessages(library(rhdf5))
  .kgSetPaths()
  h5f <- paste0(h5.path, "ma_singlecell.h5")
  if(!file.exists(h5f)) return("RES-NO-DB")
  minCells <- suppressWarnings(as.integer(minCells)); if(is.na(minCells) || minCells < 1) minCells <- 20L
  toks <- unique(trimws(strsplit(as.character(features), "[,;]")[[1]])); toks <- toks[nzchar(toks)]
  if(length(toks) == 0) return("RES-NO-FEATURE")

  # the atlas stores fixed-width strings; strip the padding exactly as the case-study reader did
  clean <- function(x) sub("\\s+$", "", gsub("\\\\000|\\x00", "", as.character(x)))
  on.exit(H5close(), add = TRUE)
  sym <- clean(h5read(h5f, "meta/genes/symbol"))
  ent <- tryCatch(clean(h5read(h5f, "meta/genes/entrez")), error = function(e) rep(NA_character_, length(sym)))
  ct  <- clean(h5read(h5f, "meta/cells/celltype"))

  # ---- the cell types: the atlas's own labels, all by default ----
  tab <- table(ct)
  types <- names(tab)
  if(!identical(tolower(trimws(as.character(celltypes))), "all")){
    want <- tolower(trimws(strsplit(as.character(celltypes), "[,;]")[[1]]))
    types <- types[tolower(types) %in% want]
    if(length(types) == 0) return("RES-NO-CELLTYPE")
  }
  types <- types[as.integer(tab[types]) >= minCells]
  if(length(types) == 0) return("RES-NO-CELLTYPE")
  write.csv(data.frame(celltype = names(tab), n_cells = as.integer(tab), summarised = names(tab) %in% types,
                       stringsAsFactors = FALSE),
            "kg_celltype_types.csv", row.names = FALSE)

  # ---- one HDF5 row per gene, summarised per cell type ----
  idx <- split(seq_along(ct), ct)
  rows <- list()
  for(tk in toks){
    j <- which(toupper(sym) == toupper(tk))[1]
    if(is.na(j)) j <- which(ent == tk)[1]
    if(is.na(j)) next
    v <- as.numeric(h5read(h5f, "data/norm_expression", index = list(j, NULL)))
    for(t in types){
      x <- v[idx[[t]]]; x <- x[!is.na(x)]
      if(length(x) == 0) next
      pos <- x > 0
      rows[[length(rows) + 1]] <- data.frame(
        feature = tk, symbol = sym[j], celltype = t, n_cells = length(x),
        mean_expression = signif(mean(x), 4), pct_expressing = round(100 * mean(pos), 2),
        mean_in_expressing = if(any(pos)) signif(mean(x[pos]), 4) else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  if(length(rows) == 0) return("RES-NO-RESOLVE")
  out <- do.call(rbind, rows)
  write.csv(out, "kg_celltype_expr.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out))
}
