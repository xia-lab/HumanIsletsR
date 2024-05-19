# Proportion stats for paper
# Author: Jessica Ewald

## Set your working directory to the "6_for_manuscript" directory

library(dplyr)
library(RSQLite)
library(ggplot2)
library(data.table)
library(pheatmap)
library(UpSetR)

source("../set_paths.R")
setPaths()

# read in data
proportions <- read.csv(paste0(other.tables.path, "outcomes_processing_input/composition.csv"))
mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
meta <- dbReadTable(mydb, "proc_metadata")
dbDisconnect(mydb)

#### Calculate correlation between purity and trapped
cor.test(meta$puritypercentage, meta$trappedpercentage, use = "complete")

#### Calculate correlation between percent exocrine and percent purity
cor.test(meta$puritypercentage, meta$exo_per, use = "complete")

#### Calculate t-stat for % purity in T2D vs. None and T1D vs. None
# t-test purity, diabetes status
res.t2d <- t.test(puritypercentage ~ diagnosis, data = meta[meta$diagnosis %in% c("Type2", "None"),], var.equal = TRUE)
res.t1d <- t.test(puritypercentage ~ diagnosis, data = meta[meta$diagnosis %in% c("Type1", "None"),], var.equal = TRUE)

### Count sample size for each omics layer
# read in processed data
omics.path <- paste0(other.tables.path, "omics_processing_input/proc/")
dat.p <- read.csv(paste0(omics.path, "proc_prot.csv"), row.names = 1)
dat.r <- read.csv(paste0(omics.path, "proc_rnaseq.csv"), row.names = 1)
dat.n <- read.csv(paste0(omics.path, "proc_nanostring_merge.csv"), row.names = 1)

sum(colnames(dat.p) %in% proportions$record_id)
sum(colnames(dat.r) %in% proportions$record_id)
sum(colnames(dat.n) %in% proportions$record_id)

### Calculate overlap between features associated with cell type across omics layers
results.path <- "../4_deconvolution_analysis/cell_proportion_associations/"
cell.types <- c("nonendo", "alpha", "beta", "delta", "gamma")

# read in association results
prot.files <- list.files(paste0(results.path, "proc_prot"),
                         full.names = TRUE)
prot <- lapply(prot.files, read.csv)
file.names <- gsub(paste0(results.path, "proc_prot\\/"), "", prot.files)
file.names <- gsub("\\.csv", "", file.names)
file.names <- gsub("_.*", "", file.names)
file.names[file.names == "exo"] <- "nonendo"
names(prot) <- file.names

rna.files <- list.files(paste0(results.path, "proc_rnaseq"),
                        full.names = TRUE)
rna <- lapply(rna.files, read.csv)
file.names <- gsub(paste0(results.path, "proc_rnaseq\\/"), "", rna.files)
file.names <- gsub("\\.csv", "", file.names)
file.names <- gsub("_.*", "", file.names)
file.names[file.names == "exo"] <- "nonendo"
names(rna) <- file.names

nano.files <- list.files(paste0(results.path, "proc_nanostring_merge"),
                         full.names = TRUE)
nano <- lapply(nano.files, read.csv)
file.names <- gsub(paste0(results.path, "proc_nanostring_merge\\/"), "", nano.files)
file.names <- gsub("\\.csv", "", file.names)
file.names <- gsub("_.*", "", file.names)
file.names[file.names == "exo"] <- "nonendo"
names(nano) <- file.names

# read in markers 
markers <- readRDS("../4_deconvolution_analysis/input_data/markers_herrera.rds")
markers <- distinct(markers[,c(1:3)])

# append marker labels to DEA results
for(i in c(1:length(cell.types))){
  cell <- cell.types[i]
  marker.temp <- markers[markers$cell_type == cell, ]
  
  # flag markers
  if (cell != "nonendo") {
    nano[[cell]]$marker = nano[[cell]]$Gene_ID %in% marker.temp$gene_id
    rna[[cell]]$marker = rna[[cell]]$Gene_ID %in% marker.temp$gene_id
    prot[[cell]]$marker = prot[[cell]]$Gene_ID %in% marker.temp$gene_id
  } else {
    nano$nonendo$marker = FALSE
    rna$nonendo$marker = FALSE
    prot$nonendo$marker = FALSE
  }
}

# Count number of features associated with each cell proportion
prot.feats <- list()
rna.feats <- list()
nano.feats <- list()
for(i in c(1:length(cell.types))){
  cell <- cell.types[i]
  
  rna.IDs <- rna[[cell]]$Gene_ID[rna[[cell]]$sig == "up"]
  prot.IDs <- prot[[cell]]$Gene_ID[prot[[cell]]$sig == "up"]
  nano.IDs <- nano[[cell]]$Gene_ID[nano[[cell]]$sig == "up"]
  
  rna.marker <- length(markers$gene_id[(markers$cell_type) == cell & (markers$gene_id %in% rna.IDs)])
  prot.marker <- length(markers$gene_id[(markers$cell_type) == cell & (markers$gene_id %in% prot.IDs)])
  nano.marker <- length(markers$gene_id[(markers$cell_type) == cell & (markers$gene_id %in% nano.IDs)])
  
  print(cell)
  
  print(length(prot.IDs))
  print(prot.marker/length(prot.IDs))
  
  print(length(rna.IDs))
  print(rna.marker/length(rna.IDs))
  
  print(length(nano.IDs))
  print(nano.marker/length(nano.IDs))
  
  prot.feats[[cell]] <- prot.IDs
  rna.feats[[cell]] <- rna.IDs
  nano.feats[[cell]] <- nano.IDs
}

# look at intersect of lists
upset(fromList(prot.feats), order.by = "freq")
upset(fromList(rna.feats), order.by = "freq")
upset(fromList(nano.feats), order.by = "freq")

#### look at correlations across omics layers
markers.both <- markers[markers$gene_id %in% prot$alpha$Gene_ID, ]
markers.both <- markers.both[markers.both$gene_id %in% rna$alpha$Gene_ID, ]

d.p <- dat.p[,colnames(dat.p) %in% colnames(dat.r)]
d.r <- dat.r[,colnames(dat.r) %in% colnames(dat.p)]
d.p <- d.p[,match(colnames(d.r), colnames(d.p))]
identical(colnames(d.p), colnames(d.r))

markers.both$cor <- NA
markers.both$p.value <- NA
for(i in c(1:dim(markers.both)[1])){
  l1 <- unlist(d.p[markers.both$gene_id[i],])
  l2 <- unlist(d.r[markers.both$gene_id[i],])
  res <- cor.test(l1, l2)
  markers.both$cor[i] <- res$estimate
  markers.both$p.value[i] <- res$p.value
}
markers.both$fdr <- p.adjust(markers.both$p.value, method = "fdr")

sum(markers.both$fdr < 0.05)
dim(markers.both)
# 122 / 195 are significantly correlated across omics layers

# compare coefficients
prot.beta <- prot$beta[,c("Gene_ID", "Coefficient", "negLogPval", "marker", "Feature", "Description")]
colnames(prot.beta) <- c("Gene_ID", "coef.prot", "nlp.prot", "marker", "Feature", "Description")
rna.beta <- rna$beta[,c("Gene_ID", "Coefficient", "negLogPval")]
colnames(rna.beta) <- c("Gene_ID", "coef.rna", "nlp.rna")
beta.signal <- merge(rna.beta, prot.beta, by = "Gene_ID")

plot(beta.signal$coef.prot, beta.signal$coef.rna)
plot(beta.signal$nlp.prot, beta.signal$nlp.rna) # nlp = neg log10 pval

cor.test(beta.signal$coef.prot, beta.signal$coef.rna)


