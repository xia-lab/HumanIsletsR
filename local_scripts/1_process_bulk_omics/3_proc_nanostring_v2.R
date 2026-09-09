# Process Nanostring data for webtool -- v2
#
# v2 of 3_proc_nanostring.R (original author: Jessica Ewald). Same method:
# per-sample median scaling -> log2 -> ComBat.
#
# What changed from v1:
#   * input is the cleaned REDCap export (3b_extract_nanostring_xlsx.py), not
#     the two raw codeset .txt files, so the 69 new donors are included
#   * THREE batches instead of two: C6555 (82), C8898 (106), NEW (69).
#     The new batch is donors R511+, confirmed three ways -- absent from both
#     raw codeset exports, absent from HI_omics_v2.sqlite (which stops at R510),
#     and matching the boundary the collaborator reported.
#   * NEG_*/POS_* spike-in controls are dropped before normalization. They are
#     technical controls, not expression, and are only present for some donors.
#     v1 never saw them (they are absent from the raw .txt exports).
#   * ComBat runs per gene group rather than once. v1 corrected only the genes
#     complete in every codeset and passed the rest through untouched; here a
#     gene is corrected across whichever batches measured it, so the 24 module-X
#     genes shared by C8898+NEW get corrected too. Only genes complete in fewer
#     than two batches remain uncorrected -- those are confounded with batch by
#     construction and cannot be corrected by any method.
#   * placeholder zeros are converted to NA before log2 (see below). v1 never
#     hit this because the raw .txt exports omit the probe entirely where the
#     REDCap export writes 0.0.
#   * PCA is written out before AND after ComBat, with the batch effect
#     quantified as the variance in each PC explained by batch. Gene groups
#     corrected across a subset of batches get their own before/after pair.
#   * a feature annotation table records, per gene, which batches measured it
#     and whether it was corrected, so single-batch genes cannot be picked up
#     for cross-batch analysis by accident.
#
# Like v1, this drops NO donors.
#
# Run:
#   cd local_scripts/1_process_bulk_omics && Rscript 3_proc_nanostring_v2.R

library(ggplot2)
library(dplyr)
library(sva)

source("../set_paths.R")
setPaths()

clean.path <- paste0(other.tables.path, "updated_data/")
proc.omics.path <- paste0(other.tables.path, "omics_processing_input/proc/")
raw.omics.path  <- paste0(other.tables.path, "omics_processing_input/unproc/")
if (!dir.exists(proc.omics.path)) { dir.create(proc.omics.path, recursive = TRUE) }
if (!dir.exists(raw.omics.path))  { dir.create(raw.omics.path, recursive = TRUE) }

fig.path <- paste0(proc.omics.path, "nanostring_v2_figures/")
if (!dir.exists(fig.path)) { dir.create(fig.path, recursive = TRUE) }

# ---- locate the newest cleaned export -------------------------------------

mat.file <- paste0(clean.path, "nanostring_clean_v2.csv")
qc.file  <- paste0(clean.path, "nanostring_clean_donor_qc_v2.csv")
stopifnot(file.exists(mat.file), file.exists(qc.file))
cat("matrix:", basename(mat.file), "\nqc    :", basename(qc.file), "\n")

# Write a CSV, moving any existing file aside to a timestamped .bak first.
writeWithBackup <- function(x, path, ...) {
  if (file.exists(path)) {
    backup <- paste0(path, ".bak_", format(Sys.time(), "%Y%m%d_%H%M"))
    file.rename(path, backup)
    cat("  backed up existing ->", basename(backup), "\n")
  }
  write.csv(x, path, ...)
}

cs <- read.csv(mat.file, row.names = 1, check.names = FALSE)
cs <- as.matrix(cs)
mode(cs) <- "numeric"

qc <- read.csv(qc.file, row.names = 1, check.names = FALSE)
stopifnot(all(colnames(cs) %in% rownames(qc)))
batch <- qc[colnames(cs), "batch"]
names(batch) <- colnames(cs)

cat("\nloaded", nrow(cs), "probes x", ncol(cs), "donors\n")
print(table(batch))

# ---- drop technical control probes ----------------------------------------

is.ctrl <- grepl("^(NEG|POS)_", rownames(cs))
cat("\ndropping", sum(is.ctrl), "NEG_/POS_ control probes ->", sum(!is.ctrl), "features\n")
cs <- cs[!is.ctrl, ]

# ---- normalize -------------------------------------------------------------
# Unchanged from v1: scale each sample by its own median, then log2. Note this
# makes the result invariant to any per-donor constant factor, which is why the
# differing normalization constants between the xlsx and the raw C6555 export
# (0.951720 / 0.945915) have no effect on the output.

normNanostring <- function(csMat) {
  csMat <- apply(csMat, 2, function(x) { x / median(x, na.rm = TRUE) })
  csMat[!is.na(csMat)] <- log2(csMat[!is.na(csMat)])
  return(csMat)
}

# Placeholder zeros -> NA, BEFORE log2. The REDCap export records probes that
# were not run as an exact 0.0 rather than leaving them blank: the 27 C8898
# donors carrying module-Y genes (CCL22, DPPA4, DPY30, KDR, KIT, KMT2A, KMT2B,
# KMT2D, WT1) are exactly 0.0 in 243/243 cells, while C6555 has 738/738 nonzero
# real values for the same genes and exact zeros appear nowhere else in the
# matrix. Left alone these become log2(0) = -Inf and poison every downstream
# step. Treating them as missing is also the only defined option.
n.zero <- sum(cs == 0, na.rm = TRUE)
if (n.zero > 0) {
  zero.genes <- rownames(cs)[rowSums(cs == 0, na.rm = TRUE) > 0]
  cat("\nplaceholder zeros -> NA:", n.zero, "cells across", length(zero.genes),
      "genes:", paste(zero.genes, collapse = ", "), "\n")
  cs[!is.na(cs) & cs == 0] <- NA
}

# Save the merged matrix at the pre-normalization stage: all three codesets
# together, control probes dropped and placeholder zeros converted to NA, but
# still on the original nCounter normalized-count scale -- no median scaling,
# no log2, no ComBat. This is the v2 counterpart of the unproc_nanostring_*
# codeset exports, which are likewise per-codeset and untransformed.
unproc.out <- paste0(raw.omics.path, "unproc_nanostring_merge_v2.csv")
writeWithBackup(cs, unproc.out)
cat("\nwrote", basename(unproc.out), " (", nrow(cs), "features x", ncol(cs),
    "donors, pre-normalization )\n")

cs.norm <- normNanostring(cs)
stopifnot(all(is.finite(cs.norm) | is.na(cs.norm)))

# ---- PCA helper ------------------------------------------------------------
# ComBat and prcomp both need complete cases, so both run on the genes measured
# in every batch. Genes unique to a subset of codesets are carried through
# uncorrected, exactly as in v1.

complete.genes <- rownames(cs.norm)[complete.cases(cs.norm)]
cat("\ngenes complete across all donors (ComBat set):", length(complete.genes), "\n")
cat("genes with gaps (passthrough, UNCORRECTED)   :",
    nrow(cs.norm) - length(complete.genes), "\n")

batchR2 <- function(scores, batch, npc = 5) {
  sapply(seq_len(min(npc, ncol(scores))), function(i) {
    summary(lm(scores[, i] ~ factor(batch)))$r.squared
  })
}

plotPCA <- function(mat, batch, title, file) {
  mat <- mat[complete.cases(mat), , drop = FALSE]
  # prcomp(scale. = TRUE) errors on a constant gene; drop any before scaling.
  keep <- apply(mat, 1, sd) > 0
  if (any(!keep)) {
    cat("  note: dropping", sum(!keep), "zero-variance genes from PCA:",
        paste(rownames(mat)[!keep], collapse = ", "), "\n")
  }
  pr <- prcomp(t(mat[keep, , drop = FALSE]), scale. = TRUE)
  x <- as.data.frame(pr$x)
  x$batch <- batch[rownames(x)]
  pct <- round((pr$sdev^2) / sum(pr$sdev^2) * 100, 1)
  r2 <- batchR2(pr$x, x$batch)

  p <- ggplot(x, aes(x = PC1, y = PC2, color = batch)) +
    geom_point(size = 2, alpha = 0.8) +
    xlab(paste0("PC1 (", pct[1], "%)")) +
    ylab(paste0("PC2 (", pct[2], "%)")) +
    theme_bw() +
    ggtitle(paste0(title, "\nvariance explained by batch: PC1 ",
                   round(r2[1] * 100, 1), "%, PC2 ", round(r2[2] * 100, 1), "%"))
  ggsave(file, p, width = 7, height = 5, dpi = 150)

  cat("\n", title, "\n", sep = "")
  cat("  PC variance %       :", paste(pct[1:5], collapse = ", "), "\n")
  cat("  batch R2 per PC (%) :", paste(round(r2 * 100, 1), collapse = ", "), "\n")
  invisible(r2)
}

r2.before <- plotPCA(cs.norm, batch, "Merged codesets, normalized (BEFORE ComBat)",
                     paste0(fig.path, "pca_before_combat_v2.png"))

# ---- batch correction ------------------------------------------------------
# mod = ~1 matches v1: no biological covariate is protected. See the note at the
# bottom of this script before changing that.

# A gene can only be corrected across the batches that actually measured it.
# Group genes by which batches have complete data, then run ComBat within each
# group over just those batches. This corrects the 123 genes shared by all three
# AND the 24 module-X genes shared by C8898+NEW, instead of passing the latter
# through uncorrected as v1 would. Genes complete in fewer than two batches
# cannot be corrected at all and are carried through as-is.

batches <- sort(unique(batch))
completeIn <- sapply(batches, function(b) {
  sub <- cs.norm[, batch == b, drop = FALSE]
  rowSums(is.na(sub)) == 0
})
rownames(completeIn) <- rownames(cs.norm)

group.key <- apply(completeIn, 1, function(r) paste(batches[r], collapse = "+"))
n.batches <- rowSums(completeIn)

cat("\ngene groups by batches with complete data:\n")
print(table(ifelse(n.batches >= 2, group.key, paste0("<uncorrectable> ", group.key))))

cs.combat <- cs.norm
corrected.in <- setNames(rep(NA_character_, nrow(cs.norm)), rownames(cs.norm))

for (key in unique(group.key[n.batches >= 2])) {
  genes <- rownames(cs.norm)[group.key == key]
  keep.batches <- strsplit(key, "+", fixed = TRUE)[[1]]
  donors <- colnames(cs.norm)[batch %in% keep.batches]

  sub <- cs.norm[genes, donors, drop = FALSE]
  sub.batch <- batch[donors]
  stopifnot(!any(is.na(sub)))

  # ComBat divides by the within-batch variance, so a gene that is constant
  # inside any batch would come back NaN. Exclude those rather than corrupt them.
  within.sd <- sapply(unique(sub.batch), function(b) {
    apply(sub[, sub.batch == b, drop = FALSE], 1, sd)
  })
  const <- rownames(sub)[apply(as.matrix(within.sd), 1, function(v) any(v == 0))]
  if (length(const) > 0) {
    cat("  WARNING: zero within-batch variance, left UNCORRECTED:",
        paste(const, collapse = ", "), "\n")
    genes <- setdiff(genes, const)
    sub <- sub[genes, , drop = FALSE]
  }
  if (length(genes) == 0) next

  cat("\nComBat on", length(genes), "genes x", length(donors), "donors  [", key, "]\n")
  fixed <- ComBat(dat = sub, batch = sub.batch,
                  mod = model.matrix(~1, data = data.frame(b = sub.batch)))
  stopifnot(all(is.finite(fixed)))
  cs.combat[genes, donors] <- fixed
  corrected.in[genes] <- key

  # The all-donor PCA cannot show a correction confined to a subset of batches
  # (those genes are NA elsewhere, so they never enter it). Plot the subset too.
  if (length(keep.batches) < length(batches)) {
    tag <- gsub("[^A-Za-z0-9]+", "_", key)
    plotPCA(cs.norm[genes, donors, drop = FALSE], batch,
            paste0("Subset ", key, ", normalized (BEFORE ComBat)"),
            paste0(fig.path, "pca_before_combat_", tag, "_v2.png"))
    plotPCA(cs.combat[genes, donors, drop = FALSE], batch,
            paste0("Subset ", key, ", batch-corrected (AFTER ComBat)"),
            paste0(fig.path, "pca_after_combat_", tag, "_v2.png"))
  }
}

stopifnot(ncol(cs.combat) == ncol(cs), nrow(cs.combat) == nrow(cs.norm))
stopifnot(identical(is.na(cs.combat), is.na(cs.norm)))

r2.after <- plotPCA(cs.combat, batch, "Merged codesets, batch-corrected (AFTER ComBat)",
                    paste0(fig.path, "pca_after_combat_v2.png"))

# ---- write out -------------------------------------------------------------

# Fixed names mirroring v1 (proc_nanostring_merge.csv), so downstream scripts
# can read them by a stable path. NOTE: these are overwritten on every run.
norm.out   <- paste0(proc.omics.path, "proc_nanostring_normalized_v2.csv")
combat.out <- paste0(proc.omics.path, "proc_nanostring_merge_v2.csv")
pca.out    <- paste0(proc.omics.path, "proc_nanostring_batch_r2_v2.csv")
feat.out   <- paste0(proc.omics.path, "proc_nanostring_feature_info_v2.csv")

# Feature annotation: which batches measured each gene, and whether it was
# batch corrected. A gene complete in only one batch is fully confounded with
# batch -- any cross-batch comparison on it measures the batch, not biology.
feat <- data.frame(
  feature          = rownames(cs.norm),
  batches_complete = group.key,
  n_batches        = n.batches,
  combat_status    = ifelse(n.batches >= 2,
                            paste0("corrected (", corrected.in, ")"),
                            ifelse(n.batches == 1,
                                   "NOT corrected - single batch, confounded",
                                   "NOT corrected - incomplete in every batch")),
  cross_batch_safe = n.batches >= 2,
  n_donors_with_data = rowSums(!is.na(cs.norm)),
  row.names = NULL
)

writeWithBackup(cs.norm, norm.out)
writeWithBackup(cs.combat, combat.out)
writeWithBackup(feat, feat.out, row.names = FALSE)
writeWithBackup(data.frame(PC = paste0("PC", seq_along(r2.before)),
                     batch_r2_before = r2.before,
                     batch_r2_after  = r2.after),
          pca.out, row.names = FALSE)

cat("\nwrote", basename(norm.out), "\n")
cat("wrote", basename(combat.out), "  (", nrow(cs.combat), "features x",
    ncol(cs.combat), "donors )\n")
cat("wrote", basename(pca.out), "\n")
cat("wrote", basename(feat.out), "\n")
cat("figures in", fig.path, "\n")

cat("\nfeature summary:\n")
print(table(feat$combat_status))

cat("\ndonors in  :", ncol(cs), "\ndonors out :", ncol(cs.combat),
    "\ndonors dropped:", ncol(cs) - ncol(cs.combat), "\n")

# ---- note ------------------------------------------------------------------
# Two open methodological points, unchanged from v1 and NOT decided here:
#
# 1. mod = ~1 protects no biological covariate. T2D is unevenly distributed
#    across codesets (11/82 in C6555 vs 5/106 in C8898), so a null model can
#    shrink real disease signal along with batch. Protecting diagnosis would
#    need donor metadata for the new batch, which donor_info_processed.csv does
#    not yet have for R511+.
#
# 2. The per-sample median is taken over whichever genes a donor carries, and
#    that gene set differs by codeset. A median over the common gene set would
#    make the scaling factor comparable across batches. Kept as-is to match v1.
