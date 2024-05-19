# Process pseudobulk
# Author: Jessica Ewald

## Set your working directory to the "2_process_patchseq" directory

library(data.table)
library(dplyr)
library(edgeR)
library(ggplot2)
library(RSQLite)
library(rhdf5)
library(HDF5Array)

source("../set_paths.R")
setPaths()

proc.omics.path <- paste0(other.tables.path, "omics_processing_input/proc/patchseq_proc/")
if(!dir.exists(proc.omics.path)){ dir.create(proc.omics.path) }

# annotate to Entrez
myDb <- dbConnect(SQLite(), paste0(other.tables.path, "libraries/hsa_genes.sqlite"))
entrez <- dbReadTable(myDb, "entrez")
ensembl <- dbReadTable(myDb, "entrez_embl_gene")
dbDisconnect(myDb)

expr <- readRDS(paste0(proc.omics.path, "counts.rds"))
meta <- readRDS(paste0(proc.omics.path, "metadata.rds"))

ensG <- rownames(expr)

hit.inx <- match(ensG, ensembl[, "accession"])
gene.ids <- ensembl$gene_id[hit.inx]
expr <- expr[!is.na(gene.ids), ]
gene.ids <- gene.ids[!is.na(gene.ids)]

# sum duplicated gene ids
expr <- cbind(gene_id = gene.ids, expr)
expr <- aggregate(.~ gene_id, data = expr, FUN = sum, na.rm = TRUE)
rownames(expr) <- expr$gene_id
expr <- expr[,-1]

# remove cells with weird glucose
cells.keep <- meta$cell_id[meta$glucose_mM %in% c("1", "5", "10")]
expr <- expr[,colnames(expr) %in% cells.keep]

# separate by cell type and glucose conc
types <- c("Alpha", "Beta")
gluc.conc <- c("1", "5", "10")
sc.lcpm <- list()
sc.names <- c()
counter <- 1
for(i in c(1:length(types))){
  for(k in c(1:length(gluc.conc))){
    type <- types[i]
    conc <- gluc.conc[k]
    
    meta.temp <- meta[meta$cell_type == type & meta$glucose_mM == conc, ]
    
    expr.temp <- expr[,colnames(expr) %in% meta.temp$cell_id]
    genes.keep <- apply(expr.temp, 1, function(x){sum(x == 0)/length(x) < 0.8}) # can have maximum 80% zeros
    expr.temp <- expr.temp[genes.keep, ]
    
    # convert to lcpm
    # this link suggested TMM is better than RLE for single-cell: https://bioinformatics-core-shared-training.github.io/SingleCell_RNASeq_June22/UnivCambridge_ScRnaSeqIntro_Base/Markdowns/05_Normalisation.html#cpm 
    nf <- calcNormFactors(expr.temp, method="TMM")
    lcpm.temp <- cpm(expr.temp, lib.size=colSums(expr.temp)*nf, log = TRUE, prior.count = 1)
    
    # replace zeros with NA so they won't be considered in the correlation analysis
    lcpm.temp[expr.temp == 0] <- NA
    
    sc.lcpm[[counter]] <- lcpm.temp
    sc.names <- c(sc.names, paste0(type, "_", conc))
    counter <- counter + 1
  }
}

names(sc.lcpm) <- sc.names
saveRDS(sc.lcpm, paste0(proc.omics.path, "sc_lcpm.rds"))

## Write out data to HDF5

# annotate to Entrez
myDb <- dbConnect(SQLite(), paste0(other.tables.path, "libraries/hsa_genes.sqlite"))
entrez <- dbReadTable(myDb, "entrez")
dbDisconnect(myDb)
meta <- readRDS(paste0(proc.omics.path, "metadata.rds"))

# add single-cell gene expression tables to hdf5 files
# build metadata of interest directly into here (?)
for(i in c(1:length(sc.names))){
  type <- sc.names[i]
  temp <- sc.lcpm[[type]]
  
  # get gene info
  feature.vec <- rownames(temp)
  hit.inx <- match(feature.vec, entrez[, "gene_id"])
  entrez.temp <- entrez[hit.inx, ]
  
  # get cell info
  meta.temp <- meta[meta$cell_id %in% colnames(temp), ]
  meta.temp <- meta.temp[match(colnames(temp), meta.temp$cell_id), ]
  
  # write out an hdf5 file here
  file.nm <- paste0("/Users/jessicaewald/hdf5/sc_", type, ".h5")
  h5createFile(file.nm)
  h5createGroup(file.nm, "meta") # group is like a folder 
  h5createGroup(file.nm, "data")
  
  writeHDF5Array(temp, file.nm, "data/norm_expression", chunkdim = c(100,100))
  H5close()
  
  h5createGroup(file.nm, "meta/genes")
  h5write(entrez.temp$gene_id, file.nm, "meta/genes/entrez")
  h5write(entrez.temp$symbol, file.nm, "meta/genes/symbol")
  h5write(entrez.temp$name, file.nm, "meta/genes/name")
  H5close()
  
  h5createGroup(file.nm, "meta/cells")
  h5write(meta.temp$cell_id, file.nm, "meta/cells/cellid")
  h5write(meta.temp$cell_score, file.nm, "meta/cells/cellscore")
  H5close()
}


