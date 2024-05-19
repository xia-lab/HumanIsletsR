# Proteomics processing
# Author: Jessica Ewald

## Set your working directory to the "1_process_bulk_omics" directory

library(WGCNA)
library(dplyr)
library(ggplot2)
library(imp4p)
library(RSQLite)

source("../set_paths.R")
setPaths()

raw.omics.path <- paste0(other.tables.path, "omics_processing_input/unproc/")
proc.omics.path <- paste0(other.tables.path, "omics_processing_input/proc/")

if(!dir.exists(proc.omics.path)){ dir.create(proc.omics.path) }

# RF works well for proteomics data: https://www.nature.com/articles/s41598-021-81279-4#:~:text=Missing%20value%20imputation%20is%20a,and%20disadvantages%20of%20each%20method
# imp4p requires log-normalized data (takes a long time to run)
# imp4p RF just a wrapper for 'MissForest' package: https://academic.oup.com/bioinformatics/article/28/1/112/219101

# get proteomics data
prot <- read.csv(paste0(raw.omics.path, "unproc_prot.csv"))
prot.ids <- prot$X
rownames(prot) <- prot$X
prot <- prot[,-1]

# remove proteins with >50% missing values
prot <- prot[goodSamplesGenes(t(prot))$goodGenes, ]

# normalize
prot <- as.matrix(prot)
rand.inds <- sample(1:dim(prot)[2], 25)
boxplot(prot[,rand.inds], main = "Raw")

prot <- apply(prot, 2, function(x){x/median(x, na.rm=T)});
boxplot(prot[,rand.inds], main = "Sample Median")

prot[!is.na(prot)] <- log2(prot[!is.na(prot)])
boxplot(prot[,rand.inds], main = "Log transformed")

prot <- impute.RF(prot, conditions = factor(rep("1", dim(prot)[2])), verbose = TRUE)
boxplot(prot[,rand.inds], main = "Imputed")

# annotate to Entrez
myDb <- dbConnect(SQLite(), paste0(other.tables.path, "libraries/hsa_genes.sqlite"))
uniprot <- dbReadTable(myDb, "entrez_uniprot")
entrez <- dbReadTable(myDb, "entrez")
dbDisconnect(myDb)

feature.vec <- rownames(prot)
hit.inx <- match(feature.vec, uniprot[, "accession"])
gene.ids <- uniprot[hit.inx, ]
prot <- as.data.frame(prot)
prot$gene.id <- gene.ids$gene_id
prot <- prot[!is.na(prot$gene.id), ]

# average duplicates
prot <- aggregate(prot[which(colnames(prot) != "gene.id")], prot[which(colnames(prot) == "gene.id")], mean)
rownames(prot) <- prot$gene.id
prot <- prot[,-which(colnames(prot) == "gene.id")]

# filter based on variance
prot.var <- apply(prot, 1, var)
var.thresh <- quantile(prot.var, 0.20)
prot.keep <- which(prot.var > var.thresh)
prot <- prot[prot.keep, ]

# write out results
write.csv(prot, paste0(proc.omics.path, "proc_prot.csv"))
