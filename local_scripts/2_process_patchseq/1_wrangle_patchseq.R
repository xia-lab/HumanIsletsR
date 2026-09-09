# Wrangle patchSeq data
# Author: Jessica Ewald
# [CHANGED] Updated by LZY 2026-03 to use RNA assay + NormalizeData for consistency
#           with parallel_pooled_spearman_meta_analysis.R (precomputed Feature View)

## Set your working directory to the "2_process_patchseq" directory

library(Seurat)
library(dplyr)
library(data.table)

# single-cell filtering thresholds taken from here: https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html

source("../set_paths.R")
setPaths()

proc.omics.path <- paste0(other.tables.path, "omics_processing_input/proc/patchseq_proc/")
if(!dir.exists(proc.omics.path)){ dir.create(proc.omics.path) }

pclamp_patched_all <- readRDS(paste0(other.tables.path, "omics_processing_input/unproc/20251029_patchseq_symbs.rds"))

# extract metadata
metadata <- pclamp_patched_all@meta.data

## [CHANGED] Use RNA assay instead of SCT
## Previously: counts <- pclamp_patched_all@assays[["SCT"]]@counts %>% as.matrix() %>% as.data.frame()
## Now: extract raw RNA counts (for pseudobulk aggregation) AND normalized RNA data (for single-cell)
DefaultAssay(pclamp_patched_all) <- "RNA"
pclamp_patched_all <- JoinLayers(pclamp_patched_all)

# Raw counts for pseudobulk (sum raw counts -> TMM -> lcpm)
raw_counts <- pclamp_patched_all[["RNA"]]$counts %>% as.matrix() %>% as.data.frame()

# Normalized data for single-cell (consistent with parallel_pooled_spearman)
## [CHANGED] Reset future plan to sequential to avoid globals size limit error
## (multisession from parallel_pooled_spearman causes NormalizeData to fail
##  when serializing the large Seurat object to workers)
if (requireNamespace("future", quietly = TRUE)) future::plan(sequential)
pclamp_patched_all <- NormalizeData(pclamp_patched_all)
norm_expr <- pclamp_patched_all[["RNA"]]$data %>% as.matrix() %>% as.data.frame()

rm(pclamp_patched_all) # remove large object from memory

## [CHANGED] Use SingleR.celltype instead of celltype for consistency with precomputed Spearman
## Previously: metadata$celltype (lowercase: alpha, beta, delta, PP)
## Now: metadata$SingleR.celltype (check values — may already be capitalized)
## NOTE: verify that SingleR.celltype values match expected labels in your data

# look at total #s before filtering
meta <- metadata[,c("Donor", "SingleR.celltype", "Glucose_mM")]
meta <- as.data.table(meta)
meta <- meta[,.N, by = .(Donor, SingleR.celltype, Glucose_mM)]
print(unique(meta$SingleR.celltype)) ## [CHANGED] print to verify cell type labels

## [CHANGED] Filter using SingleR.celltype — adjust labels if needed
## If SingleR.celltype uses capitalized names (Alpha, Beta, etc.), update the filter accordingly
ct_keep <- c("alpha", "beta", "delta", "PP", "Alpha", "Beta", "Delta") # accept both cases
cells.keep <- metadata$CellID[metadata$SingleR.celltype %in% ct_keep &
                                (metadata$nCount_RNA > 500) &
                                (metadata$nFeature_RNA > 300) &
                                (metadata$percent.mt < 20)
                              ]
raw_counts <- raw_counts[, which(colnames(raw_counts) %in% cells.keep)]
norm_expr <- norm_expr[, which(colnames(norm_expr) %in% cells.keep)]
metadata <- metadata[metadata$CellID %in% cells.keep, ]

# look at total #s after filtering
meta <- metadata[,c("Donor", "SingleR.celltype", "Glucose_mM")]
meta <- as.data.table(meta)
meta <- meta[,.N, by = .(Donor, SingleR.celltype, Glucose_mM)]
summary.filt <- meta[, .(cells = sum(N), donors = .N), by = .(SingleR.celltype, Glucose_mM)]

# extract patchseq metadata v1
# vars <- c("CellID", "Donor", "celltype", "predicted.celltype.score", "Glucose_mM", "CellSize_pF", "NormalizedTotalCapacitance_fF.pF",
# "NormalizedFirstDepolarizationCapacitance_fF.pF", "NormalizedLateDepolarizationCapacitance",
# "CalciumIntegralNormalizedtoCellSize_pC.pF", "NormalizedPeakSodiumCurrentAmplitudeat.10mV_pA.pF",
# "HalfInactivationofSodiumCurrent_mV", "NormalizedEarlyPeakCalciumCurrentAmplitudeat.10mV_pA.pF",
# "NormalizedLateCalciumCurrentAmplitudeat.10mV_pA.pF",
# "Sex", "Age", "BMI", "HbA1c", "DiabetesStatus")

## [CHANGED] Use SingleR.celltype instead of celltype
vars <-c(
  "CellID",
  "Donor",
  "SingleR.celltype",  ## [CHANGED] was "celltype"
  "SingleR.score", # "predicted.celltype.score"
  "Glucose_mM",

  "CellSize_pF",
  "mi_mean_NormalizedTotalCapacitance_fF.pF",
  "mi_mean_NormalizedFirstDepolarizationCapacitance_fF.pF",
  "mi_mean_NormalizedLateDepolarizationCapacitance",
  "mi_mean_CalciumIntegralNormalizedtoCellSize_pC.pF",
  "mi_mean_NormalizedPeakSodiumCurrentAmplitudeat.10mV_pA.pF",
  "mi_mean_HalfInactivationofSodiumCurrent_mV",
  "mi_mean_NormalizedEarlyPeakCalciumCurrentAmplitudeat.10mV_pA.pF",
  "mi_mean_NormalizedLateCalciumCurrentAmplitudeat.10mV_pA.pF",

  "Sex",
  "Age",
  "BMI",
  "HbA1c",
  "DiabetesStatus"
)

patch.meta <- metadata[,vars]

# harmonize missing values
patch.meta[patch.meta == "n/a"] <- NA
patch.meta[patch.meta == ""] <- NA
patch.meta[patch.meta == "--"] <- NA
patch.meta$DiabetesStatus[patch.meta$DiabetesStatus == "Non-Diabetic"] <- "None"
patch.meta$DiabetesStatus[patch.meta$DiabetesStatus == "T2DM"] <- "Type2"
patch.meta$DiabetesStatus[patch.meta$DiabetesStatus == "T1DM"] <- "Type1"
patch.meta$DiabetesStatus[patch.meta$DiabetesStatus == "T1DM*"] <- "Type1"

## [CHANGED] Capitalize cell types from SingleR.celltype if needed
patch.meta$SingleR.celltype <- as.character(patch.meta$SingleR.celltype)
patch.meta$SingleR.celltype[patch.meta$SingleR.celltype == "alpha"] <- "Alpha"
patch.meta$SingleR.celltype[patch.meta$SingleR.celltype == "beta"] <- "Beta"
patch.meta$SingleR.celltype[patch.meta$SingleR.celltype == "delta"] <- "Delta"
## NOTE: if SingleR.celltype is already capitalized, these lines are harmless no-ops

colnames(patch.meta) <- c("cell_id", "record_id", "cell_type", "cell_score", "glucose_mM", "cell_size_pF", "total_exocytosis_fF_pF",
                           "early_exocytosis_fF_pF", "late_exocytosis_fF_pF", "calcium_entry_pC_pF", "na_current_amp_pA_pF",
                           "na_half_inactivation_mV", "early_ca_current_pA_pF", "late_ca_current_pA_pF",
                           "donorsex", "donorage", "bodymassindex", "hba1c", "diagnosis")

# Convert to numeric
patch.meta[,c(4,6:14,16:18)] <- apply(patch.meta[,c(4,6:14,16:18)], 2, as.numeric)

# Flip signs of some variables so that greater values correspond to higher insulin secretion
patch.meta$calcium_entry_pC_pF <- -1*patch.meta$calcium_entry_pC_pF
patch.meta$na_current_amp_pA_pF <- -1*patch.meta$na_current_amp_pA_pF
patch.meta$early_ca_current_pA_pF <- -1*patch.meta$early_ca_current_pA_pF
patch.meta$late_ca_current_pA_pF <- -1*patch.meta$late_ca_current_pA_pF

# Collapse negative exocytosis values to zero
patch.meta$total_exocytosis_fF_pF[patch.meta$total_exocytosis_fF_pF < 0] <- 0
patch.meta$early_exocytosis_fF_pF[patch.meta$early_exocytosis_fF_pF < 0] <- 0
patch.meta$late_exocytosis_fF_pF[patch.meta$late_exocytosis_fF_pF < 0] <- 0

# separate types of metadata
filter.meta <- patch.meta[,c(1:5, 15:19)]
patch.meta <- patch.meta[,-c(15:19)]

## [CHANGED] Save both raw counts (for pseudobulk) and normalized expression (for single-cell)
## Previously: only saved SCT counts as counts_v2.rds
saveRDS(raw_counts, paste0(proc.omics.path, "counts_v2.rds"))        # raw RNA counts for pseudobulk
saveRDS(norm_expr, paste0(proc.omics.path, "norm_expr_v2.rds"))      # [NEW] NormalizeData output for single-cell
saveRDS(patch.meta, paste0(proc.omics.path, "metadata_v2.rds"))
saveRDS(patch.meta, paste0(other.tables.path, "outcomes_processing_input/patchseq_metadata_v2.rds")) ### this is for merging with the other ephys metadata
saveRDS(filter.meta, paste0(other.tables.path, "analysis_input/patchseq_donors_v2.rds")) ### this is for performing donor filtering



