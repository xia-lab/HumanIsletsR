# Process pseudobulk
# Author: Jessica Ewald

## Set your working directory to the "2_process_patchseq" directory

library(data.table)
library(dplyr)
library(edgeR)
library(ggplot2)
library(RSQLite)

source("../set_paths.R")
setPaths()

raw.omics.path <- paste0(other.tables.path, "omics_processing_input/unproc/")
proc.omics.path <- paste0(other.tables.path, "omics_processing_input/proc/patchseq_proc/")
if(!dir.exists(proc.omics.path)){ dir.create(proc.omics.path) }

# read in data
counts <- readRDS(paste0(proc.omics.path, "counts.rds"))
meta <- readRDS(paste0(proc.omics.path, "metadata.rds"))

# require cell score to be >0.8
identical(meta$cell_id, colnames(counts))
inds <- which(meta$cell_score > 0.8)
meta <- meta[inds, ]
counts <- counts[, inds]

# create pseudobulk profiles
# sum by cell type, donor
pb.counts <- list()
types <- unique(meta$cell_type)
for(i in c(1:length(types))){
    type <- types[i]
    print(type)
    
    meta.temp <- meta[meta$cell_type == type, ]
    cells <- meta.temp$cell_id
    donor.temp <- meta.temp$record_id
    
    # filter donors to include on those with at least 5 cells
    donor.freq <- table(donor.temp) %>% as.data.frame()
    donor.keep <- donor.freq$donor.temp[donor.freq$Freq >= 5] %>% as.character()
    
    meta.temp <- meta.temp[meta.temp$record_id %in% donor.keep, ]
    cells <- meta.temp$cell_id
    counts.temp <- counts[,colnames(counts) %in% cells]
    donor.temp <- donor.temp[donor.temp %in% donor.keep]
    
    counts.temp <- t(counts.temp)
    mode(counts.temp) <- "numeric"
    counts.temp <- as.data.frame(counts.temp)
    meta.temp <- meta.temp[meta.temp$cell_id %in% rownames(counts.temp),]
    meta.temp <- meta.temp[match(rownames(counts.temp), meta.temp$cell_id), ]
    
    counts.temp <- cbind(data.frame(record_id = meta.temp$record_id), counts.temp)
    
    dt <- as.data.table(counts.temp)
    dt <- dt[,lapply(.SD, sum), by = .(record_id)]
    
    counts.temp <- dt %>% as.data.frame()
    rownames(counts.temp) <- counts.temp$record_id
    counts.temp <- counts.temp[,-1]
    counts.temp <- t(counts.temp) %>% as.data.frame()
    
    pb.counts[[i]] <- counts.temp
}
names(pb.counts) <- types

# write out 'raw' files to omics folder
write.csv(pb.counts$Alpha, paste0(raw.omics.path, "unproc_pbrna_Alpha.csv"))
write.csv(pb.counts$Beta, paste0(raw.omics.path, "unproc_pbrna_Beta.csv"))

## convert pseudobulk data to lcpm
# annotate to Entrez
myDb <- dbConnect(SQLite(), paste0(other.tables.path, "libraries/hsa_genes.sqlite"))
entrez <- dbReadTable(myDb, "entrez")
ens <- dbReadTable(myDb, "entrez_embl_gene")
dbDisconnect(myDb)

entrez <- merge(entrez, ens, by = "gene_id")

# create lcpm matrices
pb.lcpm <- list()
types <- names(pb.counts)
for(i in c(1:length(types))){
  temp <- pb.counts[[types[i]]]
  
  # convert to Entrez
  feature.vec <- rownames(temp)
  hit.inx <- match(feature.vec, ens[, "accession"])
  entrez.temp <- ens[hit.inx, ]
  gene.ids <- ens$gene_id[hit.inx]
  
  temp <- temp[!is.na(gene.ids), ]
  entrez.temp <- entrez.temp[!is.na(gene.ids), ]
  gene.ids <- gene.ids[!is.na(gene.ids)]
  temp$gene_id <- gene.ids
  rownames(temp) <- NULL
  
  # sum duplicated IDs
  id.ind <- which(colnames(temp) == "gene_id")
  temp <- aggregate(temp[-id.ind], temp[id.ind], sum)
  rownames(temp) <- temp$gene_id
  temp <- temp[,-1]
  
  # filter by abundance
  # can have maximum 80% zeros
  genes.keep <- apply(temp, 1, function(x){sum(x == 0)/length(x) < 0.8})
  temp <- temp[genes.keep, ]
  
  # convert to lcpm
  nf <- calcNormFactors(temp, method="TMM")
  lcpm.temp <- cpm(temp, lib.size=colSums(temp)*nf, log = TRUE, prior.count = 1)
  
  # filter out expression equal to zero or 1 counts (high dropout even in pseudobulk)
  lcpm.temp[lcpm.temp == min(lcpm.temp)] <- NA
  lcpm.temp[lcpm.temp == min(lcpm.temp, na.rm = TRUE)] <- NA
  
  pb.lcpm[[i]] <- lcpm.temp
}
names(pb.lcpm) <- types

saveRDS(pb.lcpm, paste0(proc.omics.path, "pb_lcpm.rds"))
