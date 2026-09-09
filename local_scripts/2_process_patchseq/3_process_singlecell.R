# Process single-cell expression for Omics View
# Author: Jessica Ewald
# [CHANGED] Updated by LZY 2026-03: use NormalizeData log-normalized RNA values
#           instead of SCT counts -> TMM -> lcpm, for consistency with
#           parallel_pooled_spearman_meta_analysis.R (precomputed Feature View).
#           Gene filter changed from <=80% zeros to >=5% expressed (matching spearman).
#           Zero-to-NA replacement removed (spearman handles zeros in log-norm data).

## Set your working directory to the "2_process_patchseq" directory

library(data.table)
library(dplyr)
library(ggplot2)
library(RSQLite)
library(rhdf5)
library(HDF5Array)

source("../set_paths.R")
setPaths()

proc.omics.path <- paste0(other.tables.path, "omics_processing_input/proc/patchseq_proc/")
if(!dir.exists(proc.omics.path)){ dir.create(proc.omics.path) }

gene_anno <- qs::qread(paste0(other.tables.path, "libraries/gene_anno_full.qs"))

## [CHANGED] Read NormalizeData output instead of raw counts
## Previously: expr <- readRDS(paste0(proc.omics.path, "counts_v2.rds"))
## Now: use log-normalized RNA data (consistent with parallel_pooled_spearman)
expr <- readRDS(paste0(proc.omics.path, "norm_expr_v2.rds"))
meta <- readRDS(paste0(proc.omics.path, "metadata_v2.rds"))

# remove cells with weird glucose
cells.keep <- meta$cell_id[meta$glucose_mM %in% c("1", "5", "10")]
expr <- expr[,colnames(expr) %in% cells.keep]

# separate by cell type and glucose conc
types <- c("Alpha", "Beta")
gluc.conc <- c("1", "5", "10")
sc.norm <- list()  ## [CHANGED] renamed from sc.lcpm to sc.norm
sc.names <- c()
counter <- 1
for(i in c(1:length(types))){
  for(k in c(1:length(gluc.conc))){
    type <- types[i]
    conc <- gluc.conc[k]

    meta.temp <- meta[meta$cell_type == type & meta$glucose_mM == conc, ]
    expr.temp <- expr[,colnames(expr) %in% meta.temp$cell_id]

    remove_patterns <- c(
      "^RNA5S",        # 5S ribosomal RNA (Ensembl naming)
      "^RNA5-8S",      # 5.8S ribosomal RNA (Ensembl naming)
      "^RNA18S",       # 18S ribosomal RNA (Ensembl naming)
      "^RNA28S",       # 28S ribosomal RNA (Ensembl naming)
      "^RNA45S",       # 45S pre-ribosomal RNA (Ensembl naming)
      "5S-rRNA",       # 5S ribosomal RNA (alternative naming)
      "5-8S-rRNA",     # 5.8S ribosomal RNA (alternative naming)
      "^Y-RNA",        # Y RNA (non-coding, Ro particle component)
      "^MTRNR",        # mitochondrial ribosomal RNA
      "^RN7SK",        # 7SK small nuclear RNA
      "^RN7SL",        # 7SL signal recognition particle RNA
      "^RNU",          # spliceosomal small nuclear RNAs
      "^Metazoa-SRP$"  # metazoan signal recognition particle RNA
    )
    pattern <- paste(remove_patterns, collapse = "|")
    genes_to_remove <- grep(pattern, rownames(expr.temp), value = TRUE)
    expr.temp <- expr.temp[!rownames(expr.temp) %in% genes_to_remove, ]
    cat("Removed", length(genes_to_remove), "genes\n")

    ## [CHANGED] Gene filter: >=5% cells expressing (matching parallel_pooled_spearman)
    ## Previously: genes.keep <- apply(expr.temp, 1, function(x){sum(x == 0)/length(x) < 0.8})
    genes.keep <- rowMeans(expr.temp > 0) >= 0.05
    expr.temp <- expr.temp[genes.keep, ]

    ## [CHANGED] No TMM/lcpm normalization — data is already NormalizeData log-normalized
    ## Previously:
    ##   nf <- calcNormFactors(expr.temp, method="TMM")
    ##   lcpm.temp <- cpm(expr.temp, lib.size=colSums(expr.temp)*nf, log = TRUE, prior.count = 1)
    ##   lcpm.temp[expr.temp == 0] <- NA
    ## Now: use the log-normalized values directly

    sc.norm[[counter]] <- expr.temp
    sc.names <- c(sc.names, paste0(type, "_", conc))
    counter <- counter + 1
  }
}

names(sc.norm) <- sc.names
saveRDS(sc.norm, paste0(proc.omics.path, "sc_lcpm_v2.rds"))  # keep filename for backward compatibility

## Write out data to HDF5

# annotate to Entrez
# myDb <- dbConnect(SQLite(), paste0(other.tables.path, "libraries/hsa_genes.sqlite"))
# entrez <- dbReadTable(myDb, "entrez")
# dbDisconnect(myDb)
meta <- readRDS(paste0(proc.omics.path, "metadata_v2.rds"))

# add single-cell gene expression tables to hdf5 files
# build metadata of interest directly into here (?)
for(i in c(7:length(sc.names))){
  type <- sc.names[i]
  temp <- sc.norm[[type]]  ## [CHANGED] was sc.lcpm, renamed to sc.norm
  
  # get gene info
  feature.vec <- rownames(temp)
  #hit.inx <- match(feature.vec, entrez[, "gene_id"])
  entrez.temp <- gene_anno[match(feature.vec,gene_anno$input_symbol),c("input_symbol","ensembl_gene_id","entrez_selected","description")]
  
  # get cell info
  meta.temp <- meta[meta$cell_id %in% colnames(temp), ]
  meta.temp <- meta.temp[match(colnames(temp), meta.temp$cell_id), ]
  
  # write out an hdf5 file here
  file.nm <- paste0("/Users/lzy/humanislet/hdf5_v2/sc_", type, ".h5")
  h5createFile(file.nm)
  h5createGroup(file.nm, "meta") # group is like a folder 
  h5createGroup(file.nm, "data")
  
  cd <- c(min(50, nrow(temp)), min(50, ncol(temp)))
  writeHDF5Array(temp, file.nm, "data/norm_expression", chunkdim = cd)
  H5close()
  
  h5createGroup(file.nm, "meta/genes")
  h5write(entrez.temp$entrez_selected, file.nm, "meta/genes/entrez")
  h5write(entrez.temp$input_symbol, file.nm, "meta/genes/symbol")
  h5write(entrez.temp$description, file.nm, "meta/genes/name")
  h5write(entrez.temp$ensembl_gene_id, file.nm, "meta/genes/ensembl")
  H5close()
  
  h5createGroup(file.nm, "meta/cells")
  h5write(meta.temp$cell_id, file.nm, "meta/cells/cellid")
  h5write(meta.temp$cell_score, file.nm, "meta/cells/cellscore")
  H5close()
}


