# ==============================================================================
# kg_functions/kg_singlecell.R  --  single-cell patch-seq gene <-> electrophysiology
# (F33 / op#12). Patch-seq measures transcriptome AND ephys on the SAME cell, so the
# correlation is across cells (donor aggregation would destroy it).
#
# Self-contained: reads the gene's expression row directly from the hdf5 (a one-row
# hyperslab -> fast, no full-matrix load) and correlates it (Spearman) with the ephys
# feature across cells, per cell type x glucose. Nothing here calls another folder's
# function. Expression <- hdf5_V2/sc_<cell>_<glucose>.h5 (symbol-keyed); ephys <-
# HI_tables.sqlite ephys_cell (per cell_id). NO donor covariate adjustment (per-cell).
#
# emit_cells  "false" (DEFAULT) | "true" -> ALSO write kg_sc_cells.csv, the PER-CELL values
#             behind the grid: cell_id, cell_type, glucose, feature, expression, ephys,
#             ephys_value, n_cells. The correlation already reads exactly these and then
#             discards them, so this costs one extra write and no extra I/O.
#             ⚠ It is what F33's per-cell SCATTER needs -- "ONE POINT PER CELL, the unit of
#             analysis". The grid (kg_sc_ephys.csv) serves the MATRIX view; the matrix
#             cannot be un-summarised back into cells, so without this the scatter has no
#             source at all.
#             ⚠ Cells are collected BEFORE the >=10 correlation floor, so an under-powered
#             grid cell can still be DRAWN and labelled under-powered rather than vanishing
#             (F33: "a grid cell under 10 cells is UNDER-POWERED, never a null").
#
# Returns "RES-OK;<n_rows>" (kg_sc_ephys.csv) | "RES-NO" | "RES-NO-RESOLVE"
#   | "RES-NO-RESOLVER" | "RES-NO-EPHYS".
#
# ⚠ KNOWN SPEC GAP, NOT ADDRESSED HERE (flagged 2026-07-20): the GRID still SKIPS a cell
# with <10 cells silently (`if(length(x) < 10) next`), so kg_sc_ephys.csv contains no row
# for it. F33 requires it to be reported as under-powered, distinct from "not measured".
# Changing that alters the existing grid output, so it needs a ruling first. `emit_cells`
# partly mitigates it: the cells are in kg_sc_cells.csv even when the grid row is absent.
# ==============================================================================

.kgEphysCols <- c("cell_size_pF_cell", "total_exocytosis_fF_pF_cell",
  "early_exocytosis_fF_pF_cell", "late_exocytosis_fF_pF_cell", "calcium_entry_pC_pF_cell",
  "na_current_amp_pA_pF_cell", "na_half_inactivation_mV_cell",
  "early_ca_current_pA_pF_cell", "late_ca_current_pA_pF_cell")

# read one gene's per-cell expression from a patch-seq hdf5 (symbol-keyed, one-row slab).
.kgScGeneExpr <- function(h5f, symbol){
  genes  <- rhdf5::h5read(h5f, "meta/genes")   # a LIST: ensembl/entrez/name/symbol
  gs     <- if(!is.null(genes$symbol)) as.character(genes$symbol) else as.character(unlist(genes))
  idx <- which(toupper(gs) == toupper(symbol))
  if(length(idx) == 0) return(NULL)
  idx <- idx[1]
  expr <- try(rhdf5::h5read(h5f, "data/norm_expression", index = list(idx, NULL)), silent = TRUE)
  if(inherits(expr, "try-error")) return(NULL)
  cellid <- as.character(rhdf5::h5read(h5f, "meta/cells/cellid"))
  rhdf5::H5close()
  v <- as.numeric(expr)
  if(length(v) != length(cellid)) return(NULL)
  setNames(v, cellid)
}

kgSingleCellEphys <- function(var_id, ephys_var, var_type = "auto",
                              cell = "Alpha,Beta", glucose = "1,5,10",
                              emit_cells = "false", mode = "tool"){
  library(RSQLite); library(DBI); library(rhdf5)
  .kgSetPaths()
  if(!(ephys_var %in% .kgEphysCols)) return("RES-NO-EPHYS")
  emit <- tolower(trimws(as.character(emit_cells))) %in% c("true", "1", "yes")

  res <- .kgLoadFeatureResolver(); if(is.null(res)) return("RES-NO-RESOLVER")
  r <- .kgResolveFeature(res, var_id, var_type)
  if(length(r$layers) == 0 || !identical(r$type, "gene")) return("RES-NO-RESOLVE")

  con <- dbConnect(SQLite(), paste0(sqlite.path, "HI_omics_v2.sqlite"))
  sym <- .kgGeneSymbol(con, r$layers); dbDisconnect(con)
  if(is.na(sym)) return("RES-NO-RESOLVE")

  # ephys per cell_id (once)
  ht <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
  ephys <- dbReadTable(ht, "ephys_cell"); dbDisconnect(ht)
  ecol  <- setNames(suppressWarnings(as.numeric(ephys[[ephys_var]])), as.character(ephys$cell_id))

  cells <- trimws(strsplit(cell,    "[,;]")[[1]]); cells <- cells[nzchar(cells)]
  glucs <- trimws(strsplit(glucose, "[,;]")[[1]]); glucs <- glucs[nzchar(glucs)]

  rows <- list(); cellrows <- list()
  for(cl in cells){
    for(gl in glucs){
      h5f <- paste0(h5.v2.path, "sc_", cl, "_", gl, ".h5")
      if(!file.exists(h5f)) next
      ev <- .kgScGeneExpr(h5f, sym); if(is.null(ev)) next
      common <- intersect(names(ev), names(ecol))
      x <- ev[common]; y <- ecol[common]
      ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]

      # per-cell values for the F33 scatter -- ONE POINT PER CELL, which is the unit of
      # analysis. They are collected BEFORE the >=10 guard so an under-powered grid cell
      # can still be drawn and labelled as such rather than vanishing.
      if(emit && length(x) > 0){
        cellrows[[length(cellrows)+1]] <- data.frame(
          cell_id = names(x), cell_type = cl, glucose = gl, feature = sym,
          expression = as.numeric(x), ephys = ephys_var, ephys_value = as.numeric(y),
          n_cells = length(x), stringsAsFactors = FALSE)
      }

      if(length(x) < 10) next
      ct <- try(suppressWarnings(stats::cor.test(x, y, method = "spearman")), silent = TRUE)
      if(inherits(ct, "try-error")) next
      rr <- unname(ct$estimate)
      rows[[length(rows)+1]] <- data.frame(
        Feature = sym, Cell = cl, Glucose = gl, N_cells = length(x),
        Spearman_r = signif(rr, 3), P_value = signif(ct$p.value, 3),
        Direction = ifelse(is.na(rr), "NA", ifelse(rr > 0, "up", "down")),
        Ephys = ephys_var, stringsAsFactors = FALSE)
    }
  }
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); safe <- function(s) gsub("[^A-Za-z0-9]", "_", s)

  # the per-cell file is written whenever cells were collected -- INCLUDING when no grid
  # cell reached the >=10 correlation floor. That case is a real, drawable answer
  # ("we have cells, too few to correlate"), and it is exactly what must not be
  # indistinguishable from "we have nothing".
  if(emit && length(cellrows) > 0){
    cw <- do.call(rbind, cellrows)
    utils::write.csv(cw, paste0("kg_sc_cells_", safe(var_id), "_", safe(ephys_var), "_", ts, ".csv"), row.names = FALSE)
    utils::write.csv(cw, "kg_sc_cells.csv", row.names = FALSE)
  }

  if(length(rows) == 0) return("RES-NO")

  out <- do.call(rbind, rows); out <- out[order(out$P_value), ]
  utils::write.csv(out, paste0("kg_sc_ephys_", safe(var_id), "_", safe(ephys_var), "_", ts, ".csv"), row.names = FALSE)
  utils::write.csv(out, "kg_sc_ephys.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out))
}
