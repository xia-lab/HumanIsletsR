# Wrangle patchSeq data
# Author: Jessica Ewald

## Set your working directory to the "2_process_patchseq" directory

library(Seurat)
library(dplyr)
library(data.table)

# single-cell filtering thresholds taken from here: https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html

source("../set_paths.R")
setPaths()

proc.omics.path <- paste0(other.tables.path, "omics_processing_input/proc/patchseq_proc/")
if(!dir.exists(proc.omics.path)){ dir.create(proc.omics.path) }

load(paste0(other.tables.path, "omics_processing_input/raw/20231128_pclamp_integrated_patched-only_joined.RData"))

# extract data
metadata <- pclamp_patched_all@meta.data
counts <- pclamp_patched_all@assays[["SCT"]]@counts %>% as.matrix() %>% as.data.frame()
rm(pclamp_patched_all) # remove large object from memory

# look at total #s before filtering
meta <- metadata[,c("Donor", "celltype", "Glucose_mM")]
meta <- as.data.table(meta)
meta <- meta[,.N, by = .(Donor, celltype, Glucose_mM)]
meta <- meta[celltype %in% c("alpha", "beta")]
summary <- meta[, .(cells = sum(N), donors = .N), by = .(celltype, Glucose_mM)]

# remove cells other than alpha and beta, plus QA/QC filters, plus cryo donors, plus R296
cells.keep <- metadata$CellID[(metadata$celltype == "alpha" | metadata$celltype == "beta") &
                                (metadata$nCount_RNA > 500) &
                                (metadata$nFeature_RNA > 300) &
                                (metadata$percent.mt < 20) &
                                (!(metadata$Donor %in% c("R119", "R132", "R079", "R296")))]
counts <- counts[, which(colnames(counts) %in% cells.keep)]
metadata <- metadata[metadata$CellID %in% cells.keep, ]

# look at total #s after filtering
meta <- metadata[,c("Donor", "celltype", "Glucose_mM")]
meta <- as.data.table(meta)
meta <- meta[,.N, by = .(Donor, celltype, Glucose_mM)]
meta <- meta[celltype %in% c("alpha", "beta")]
summary.filt <- meta[, .(cells = sum(N), donors = .N), by = .(celltype, Glucose_mM)]

# extract patchseq metadata
vars <- c("CellID", "Donor", "celltype", "predicted.celltype.score", "Glucose_mM", "CellSize_pF", "NormalizedTotalCapacitance_fF.pF", 
"NormalizedFirstDepolarizationCapacitance_fF.pF", "NormalizedLateDepolarizationCapacitance", 
"CalciumIntegralNormalizedtoCellSize_pC.pF", "NormalizedPeakSodiumCurrentAmplitudeat.10mV_pA.pF", 
"HalfInactivationofSodiumCurrent_mV", "NormalizedEarlyPeakCalciumCurrentAmplitudeat.10mV_pA.pF", 
"NormalizedLateCalciumCurrentAmplitudeat.10mV_pA.pF",
"Sex", "Age", "BMI", "HbA1c", "DiabetesStatus")

patch.meta <- metadata[,vars]

# harmonize missing values
patch.meta[patch.meta == "n/a"] <- NA
patch.meta[patch.meta == ""] <- NA
patch.meta[patch.meta == "--"] <- NA
patch.meta$DiabetesStatus[patch.meta$DiabetesStatus == "Non-Diabetic"] <- "None"
patch.meta$DiabetesStatus[patch.meta$DiabetesStatus == "T2DM"] <- "Type2"
patch.meta$DiabetesStatus[patch.meta$DiabetesStatus == "T1DM"] <- "Type1"
patch.meta$DiabetesStatus[patch.meta$DiabetesStatus == "T1DM*"] <- "Type1"
patch.meta$celltype <- as.character(patch.meta$celltype)
patch.meta$celltype[patch.meta$celltype == "alpha"] <- "Alpha"
patch.meta$celltype[patch.meta$celltype == "beta"] <- "Beta"

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

# save as an object
saveRDS(counts, paste0(proc.omics.path, "counts.rds"))
saveRDS(patch.meta, paste0(proc.omics.path, "metadata.rds"))
saveRDS(patch.meta, paste0(other.tables.path, "processing_input/patchseq_metadata.rds")) ### this is for merging with the other ephys metadata
saveRDS(filter.meta, paste0(other.tables.path, "analysis_input/patchseq_donors.rds")) ### this is for performing donor filtering

