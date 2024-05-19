# Process Raw RNA-seq counts for web-tool
# Author: Jessica Ewald

## Set your working directory to the "1_process_bulk_omics" directory

library(dplyr)
library(ggplot2)
library(RSQLite)
library(edgeR)
library(sva)

source("../set_paths.R")
setPaths()

raw.omics.path <- paste0(other.tables.path, "omics_processing_input/unproc/")
proc.omics.path <- paste0(other.tables.path, "omics_processing_input/proc/")
if(!dir.exists(proc.omics.path)){ dir.create(proc.omics.path) }

# read in counts
fileName <- paste0(raw.omics.path, "unproc_rnaseq_ensembl.txt")
orig.counts <- data.table::fread(fileName, header=TRUE, check.names=FALSE, data.table=FALSE)

# annotate to Entrez
myDb <- dbConnect(SQLite(), paste0(other.tables.path, "libraries/hsa_genes.sqlite"))
ensg <- dbReadTable(myDb, "entrez_embl_gene")
entrez <- dbReadTable(myDb, "entrez")
dbDisconnect(myDb)

feature.vec <- orig.counts$V1
hit.inx <- match(feature.vec, ensg[, "accession"])
gene.ids <- ensg[hit.inx, ]
counts <- orig.counts
counts$V1 <- gene.ids$gene_id
counts <- counts[!is.na(counts$V1), ]

# remove duplicates
counts <- aggregate(counts[-1], counts[1], sum)
rownames(counts) <- counts$V1
counts <- counts[,-1]
write.csv(counts, paste0(raw.omics.path, "unproc_rnaseq.csv"))

# remove 80% zeros
num.zeros <- apply(counts, 1, function(x){sum(x == 0)})
inds.remove <- which(num.zeros > 0.8*dim(counts)[2])
counts <- counts[-inds.remove, ]

# normalize for sequencing depth
nf <- calcNormFactors(counts, method="RLE")
lcpm <- cpm(counts, lib.size=colSums(counts)*nf, log = TRUE)

# PCA before batch effect
batch <- read.table(paste0(raw.omics.path, "rnaseq_batch_info.txt"), sep="\t")
colnames(batch) <- c("record_id", "batch")

pr <- prcomp(t(lcpm), scale = TRUE)
pr.x <- pr$x
pr.x <- merge(pr.x, batch, by.x = "row.names", by.y = "record_id")

ggplot(pr.x, aes(x = PC1, y = PC2, color = batch)) +
  geom_point() +
  xlab(paste0("PC1 (", round((pr$sdev[1]^2)/sum(pr$sdev^2)*100,1), "%)")) +
  ylab(paste0("PC2 (", round((pr$sdev[2]^2)/sum(pr$sdev^2)*100,1), "%)")) +
  theme_bw() +
  ggtitle("LCPM")

# one batch (Batch2_Oxford) is EXTREMELY different. remove and re-normalize
batch <- batch[batch$batch != "Batch2_Oxford",]
counts <- counts[,colnames(counts) %in% batch$record_id]
nf <- calcNormFactors(counts, method="RLE")
lcpm <- cpm(counts, lib.size=colSums(counts)*nf, log = TRUE)

pr <- prcomp(t(lcpm), scale = TRUE)
pr.x <- pr$x
pr.x <- merge(pr.x, batch, by.x = "row.names", by.y = "record_id", all = F)

ggplot(pr.x, aes(x = PC1, y = PC2, color = batch)) +
  geom_point() +
  xlab(paste0("PC1 (", round((pr$sdev[1]^2)/sum(pr$sdev^2)*100,1), "%)")) +
  ylab(paste0("PC2 (", round((pr$sdev[2]^2)/sum(pr$sdev^2)*100,1), "%)")) +
  theme_bw() +
  ggtitle("LCPM (Oxford2 batch removed)")

# adjust for batch effect
identical(batch$record_id, colnames(counts)) # check to make sure order preserved
counts <- as.matrix(counts)
adjusted <- ComBat_seq(counts, batch=batch$batch, group=NULL)

nf <- calcNormFactors(adjusted, method="RLE")
lcpm <- cpm(adjusted, lib.size=colSums(adjusted)*nf, log = TRUE)

pr <- prcomp(t(lcpm), scale = TRUE)
pr.x <- pr$x
pr.x <- merge(pr.x, batch, by.x = "row.names", by.y = "record_id", all = F)

ggplot(pr.x, aes(x = PC1, y = PC2, color = batch)) +
  geom_point() +
  xlab(paste0("PC1 (", round((pr$sdev[1]^2)/sum(pr$sdev^2)*100,1), "%)")) +
  ylab(paste0("PC2 (", round((pr$sdev[2]^2)/sum(pr$sdev^2)*100,1), "%)")) +
  theme_bw() +
  ggtitle("LCPM (adjusted for batch effect)")

lcpm[lcpm == min(lcpm)] <- NA
lcpm[lcpm == min(lcpm, na.rm = TRUE)] <- NA

# filter based on variance
lcpm.var <- apply(lcpm, 1, var, na.rm = TRUE)
lcpm.mean <- apply(lcpm, 1, mean, na.rm = TRUE)
plot(lcpm.mean, lcpm.var)
var.thresh <- quantile(lcpm.var, 0.20)
lcpm.keep <- which(lcpm.var > var.thresh)
lcpm <- lcpm[lcpm.keep, ]

# write out files
write.csv(adjusted, paste0(raw.omics.path, "unproc_rnaseq_batch.csv"))
write.csv(lcpm, paste0(proc.omics.path, "proc_rnaseq.csv"))


