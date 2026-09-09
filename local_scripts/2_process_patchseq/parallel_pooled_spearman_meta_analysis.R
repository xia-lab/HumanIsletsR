## Parallel pooled (Rubin + study-meta) for PatchSeq MI correlations
## - Parallelizes over celltypes and optionally over imputations
## - Shows progress bars
## - Produces pooled_by_ct

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(future)
  library(future.apply)
  library(progressr)
})

## ----------------------------
## Setup (same as your script)
## ----------------------------
pclamp_patched_all <- readRDS("20251029_patchseq_symbs.rds") #replace with the other object
imp_bundle <- readRDS("20251210_imp_bundle.rds")

id_col    <- "CellID"
study_col <- "Study"

.to_num <- function(x) {
  if (is.ordered(x) || is.factor(x)) suppressWarnings(as.numeric(as.character(x)))
  else if (is.character(x))          suppressWarnings(as.numeric(x))
  else if (is.logical(x))            as.numeric(x)
  else                               x
}

.align_df_to_cells <- function(df, cell_order, id_col) {
  df <- df[df[[id_col]] %in% cell_order, , drop = FALSE]
  df <- df[match(cell_order, df[[id_col]]), , drop = FALSE]
  stopifnot(identical(df[[id_col]], cell_order))
  df
}

DefaultAssay(pclamp_patched_all) <- "RNA"
pclamp_patched_all <- JoinLayers(pclamp_patched_all) #if necessary

pclamp_patched_all <- NormalizeData(pclamp_patched_all)

# expression (cells x genes)
mRNA_data <- pclamp_patched_all[["RNA"]]$data
mRNA_data <- as.matrix(mRNA_data)
mRNA_data <- t(mRNA_data)
storage.mode(mRNA_data) <- "double"

# global gene filter (>=5% cells): can adjust depending on run time
#keep <- Matrix::colMeans(mRNA_data > 0) >= 0.05
#mRNA_data <- mRNA_data[, keep, drop = FALSE]

# metadata used for celltype labels (aligned to mRNA_data)
mi_mean <- pclamp_patched_all@meta.data %>%
  dplyr::select(CellID, SingleR.celltype, CellSize_pF, contains("mi_mean"))

cell_order <- rownames(mRNA_data)

mi_mean2 <- .align_df_to_cells(mi_mean, cell_order, id_col)

meta_df <- pclamp_patched_all@meta.data
meta_df <- meta_df[cell_order, , drop = FALSE]
stopifnot(identical(rownames(meta_df), cell_order))

# MI completed datasets aligned to cell_order
completed_list <- imp_bundle$completed_list
completed_list_filt <- lapply(completed_list, function(d) .align_df_to_cells(d, cell_order, id_col))

stopifnot(
  is.list(completed_list_filt),
  length(completed_list_filt) >= 2,
  all(vapply(completed_list_filt, nrow, integer(1)) == nrow(mRNA_data)),
  all(vapply(completed_list_filt, function(d) identical(d[[id_col]], cell_order), logical(1)))
)

## ----------------------------
## Ephys vars (same)
## ----------------------------
ephys_vars_norm <- c(
  "NormalizedTotalCapacitance_fF.pF",
  "NormalizedFirstDepolarizationCapacitance_fF.pF",
  "NormalizedLateDepolarizationCapacitance",
  "CalciumIntegralNormalizedtoCellSize_pC.pF",
  "CapacitanceNormalizedtoCalcium_fF.pC",
  "NormalizedEarlyPeakCalciumCurrentAmplitudeat.10mV_pA.pF",
  "NormalizedLateCalciumCurrentAmplitudeat.10mV_pA.pF",
  "NormalizedLVA_.60mV_pA.pF",
  "NormalizedHVA_.20mV_pA.pF",
  "NormalizedPeakSodiumCurrentAmplitudeat.10mV_pA.pF",
  "HalfInactivationofSodiumCurrent_mV",
  "VoltageforSodiumPeakCurrent_mV",
  "ReversalPotentialbyramp_mV",
  "NormalizedHyperpolarizationactivatedcurrent.at.140mV_pA.pF"
)

vars_to_test <- c("CellSize_pF", ephys_vars_norm)

celltypes <- sort(unique(completed_list_filt[[1]]$SingleR.celltype))

## ----------------------------
## Core functions (same math)
## ----------------------------
.meta_by_study <- function(y, groups, gene_mat) {
  y <- .to_num(y)
  groups <- factor(groups)
  levs <- levels(groups)
  
  rs <- vector("list", length(levs))
  ns <- integer(length(levs))
  
  for (i in seq_along(levs)) {
    ix <- which(groups == levs[i] & is.finite(y))
    if (length(ix) >= 4L) {
      rs[[i]] <- suppressWarnings(cor(
        y[ix],
        gene_mat[ix, , drop = FALSE],
        method = "spearman",
        use = "pairwise.complete.obs"
      ))
      ns[i] <- sum(is.finite(y[ix]))
    } else {
      rs[[i]] <- rep(NA_real_, ncol(gene_mat))
      ns[i] <- 0L
    }
  }
  
  R <- do.call(rbind, rs)                         # S x G
  R <- pmax(pmin(R, 0.999999), -0.999999)
  Z <- atanh(R)
  w <- pmax(ns - 3L, 0L); w[w == 0] <- NA_real_
  
  Wmat <- matrix(w, nrow = length(w), ncol = ncol(Z), byrow = TRUE)
  num  <- colSums(Z * Wmat, na.rm = TRUE)
  den  <- colSums(Wmat,     na.rm = TRUE)
  
  zbar <- num / den
  se_z <- 1 / sqrt(den)
  rbar <- tanh(zbar)
  
  list(r = rbar, z = zbar, se_z = se_z)
}

pool_rubin_one <- function(meta_by_imp_for_var, gene_names) {
  m <- length(meta_by_imp_for_var)
  if (m < 2) stop("Need >=2 imputations for Rubin pooling")
  
  Z   <- sapply(meta_by_imp_for_var, `[[`, "z")    # G x m
  SEz <- sapply(meta_by_imp_for_var, `[[`, "se_z") # G x m
  
  W  <- rowMeans(SEz^2, na.rm = TRUE)              # within-imp var
  B  <- apply(Z, 1L, var, na.rm = TRUE)            # between-imp var
  Tz <- W + (1 + 1/m) * B                          # total var
  
  zbar  <- rowMeans(Z, na.rm = TRUE)
  se    <- sqrt(Tz)
  zstat <- zbar / se
  
  data.frame(
    gene  = gene_names,
    z     = zbar,
    se_z  = se,
    r     = tanh(zbar),
    zstat = zstat,
    p     = 2 * pnorm(-abs(zstat)),
    stringsAsFactors = FALSE
  )
}

## ----------------------------
## Parallel worker for ONE celltype
## ----------------------------
.pool_one_celltype <- function(ct,
                               completed_list_filt,
                               mRNA_data,
                               vars_to_test,
                               study_col,
                               log_file,
                               min_cells = 10L) {
  
  ix <- which(completed_list_filt[[1]]$SingleR.celltype == ct)
  if (length(ix) < min_cells) return(NULL)
  
  # expression subset + celltype-specific gene filter
  gene_sub <- mRNA_data[ix, , drop = FALSE]
  keep_genes <- Matrix::colMeans(gene_sub > 0) >= 0.05
  gene_sub <- gene_sub[, keep_genes, drop = FALSE]
  gene_names <- colnames(gene_sub)
  
  if (ncol(gene_sub) < 10L) return(NULL)
  
  cat(sprintf("[%s] ct=%s: computing meta for %d imputations\n",
              Sys.time(), ct, length(completed_list_filt)),
      file = log_file, append = TRUE)
  
  # meta per imputation (serial, with per-imp heartbeat)
  meta_list_ct <- lapply(seq_along(completed_list_filt), function(k) {
    cat(sprintf("[%s] ct=%s: imp %d/%d\n",
                Sys.time(), ct, k, length(completed_list_filt)),
        file = log_file, append = TRUE)
    
    df_k <- completed_list_filt[[k]]
    df_sub <- df_k[ix, , drop = FALSE]
    grp    <- df_sub[[study_col]]
    
    out <- lapply(vars_to_test, function(v) {
      .meta_by_study(
        y        = df_sub[[v]],
        groups   = grp,
        gene_mat = gene_sub
      )
    })
    names(out) <- vars_to_test
    out
  })
  
  pooled_vars <- lapply(vars_to_test, function(v) {
    df <- pool_rubin_one(lapply(meta_list_ct, `[[`, v), gene_names = gene_names)
    df$celltype  <- ct
    df$ephys_var <- v
    df
  })
  names(pooled_vars) <- vars_to_test
  
  pooled_vars
}
## ----------------------------
## Run in parallel + progress
## ----------------------------
# Choose workers
n_workers <- max(1L, future::availableCores()/2)
options(future.globals.maxSize = 3e+09)
plan(multisession, workers = n_workers)


progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")  # <- key change for terminal
options(progressr.enable = TRUE)

log_file <- "pooled_progress.log"
cat(sprintf("[%s] Starting run\n", Sys.time()), file = log_file, append = TRUE)

## Choose workers conservatively for local machine
n_workers <- max(1L, future::availableCores()/2)
future::plan(multisession, workers = n_workers)

with_progress({
  p <- progressor(along = celltypes)
  
  pooled_list <- future_lapply(celltypes, function(ct) {
    # heartbeat at start (goes to log file)
    cat(sprintf("[%s] START ct=%s\n", Sys.time(), ct), file = log_file, append = TRUE)
    
    p(sprintf("Pooling: %s", ct))
    
    t0 <- proc.time()[["elapsed"]]
    res <- .pool_one_celltype(
      ct = ct,
      completed_list_filt = completed_list_filt,
      mRNA_data = mRNA_data,
      vars_to_test = vars_to_test,
      study_col = study_col,
      log_file = log_file
    )
    dt <- proc.time()[["elapsed"]] - t0
    
    # heartbeat at end
    cat(sprintf("[%s] END   ct=%s elapsed=%.1fs\n", Sys.time(), ct, dt),
        file = log_file, append = TRUE)
    
    res
  }, future.seed = TRUE)
})
qs::qsave(pooled_list,"pooled_list.qs")
pooled_by_ct <- setNames(pooled_list, celltypes)
pooled_by_ct <- pooled_by_ct[!vapply(pooled_by_ct, is.null, logical(1))]

cat(sprintf("[%s] Finished run\n", Sys.time()), file = log_file, append = TRUE)
saveRDS(pooled_by_ct, file = "pooled_by_ct_parallel.rds")
