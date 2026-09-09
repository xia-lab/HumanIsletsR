# Process pseudobulk
# Author: Jessica Ewald
# [CHANGED] Updated by LZY 2026-03: input counts_v2.rds now contains RNA raw counts
#           (previously SCT counts). Pseudobulk pipeline (sum -> TMM -> lcpm) is unchanged
#           since raw counts are the correct input for pseudobulk aggregation.

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
## [CHANGED] counts_v2.rds now contains RNA raw counts (was SCT counts)
## This is correct for pseudobulk: raw counts -> sum by donor -> TMM -> lcpm
counts <- readRDS(paste0(proc.omics.path, "counts_v2.rds"))
meta <- readRDS(paste0(proc.omics.path, "metadata_v2.rds"))

# require cell score to be >0.8
identical(meta$cell_id, colnames(counts))
#inds <- which(meta$cell_score > 0.8)  skip in v2
#meta <- meta[inds, ]
#counts <- counts[, inds]

# create pseudobulk profiles
# sum by cell type, donor
pb.counts <- list()
types <- unique(meta$cell_type)
for(i in c(2:length(types))){
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
write.csv(pb.counts$Alpha, paste0(raw.omics.path, "unproc_pbrna_Alpha_v2.csv"))
write.csv(pb.counts$Beta, paste0(raw.omics.path, "unproc_pbrna_Beta_v2.csv"))

## convert pseudobulk data to lcpm
# annotate to Entrez
# myDb <- dbConnect(SQLite(), paste0(other.tables.path, "libraries/hsa_genes.sqlite"))
# entrez <- dbReadTable(myDb, "entrez")
# ens <- dbReadTable(myDb, "entrez_embl_gene")
# dbDisconnect(myDb)
# 
# entrez <- merge(entrez, ens, by = "gene_id")

 
# create lcpm matrices
pb.lcpm <- list()
types <- names(pb.counts)
for(i in c(2:length(types))){
  temp <- pb.counts[[types[i]]]
  
  # convert to Entrez
   #feature.vec <- rownames(temp)
  # hit.inx <- match(feature.vec, out_anno$hgnc_corrected_symbol)
  # entrez.temp <- ens[hit.inx, ]
 #gene.ids <- ens$gene_id[hit.inx]
  
  # temp <- temp[!is.na(gene.ids), ]
  # entrez.temp <- entrez.temp[!is.na(gene.ids), ]
  # gene.ids <- gene.ids[!is.na(gene.ids)]
#  temp$gene_id <- gene.ids
   #rownames(temp) <- NULL
  
  # sum duplicated IDs
  # id.ind <- which(colnames(temp) == "gene_id")
  # temp <- aggregate(temp[-id.ind], temp[id.ind], sum)
  # rownames(temp) <- temp$gene_id
  # temp <- temp[,-1]
  
  ## remove informative genes
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
  genes_to_remove <- grep(pattern, rownames(temp), value = TRUE)
  temp <- temp[!rownames(temp) %in% genes_to_remove, ]
  cat("Removed", length(genes_to_remove), "genes\n")
  grep("rRNA|RNA|snoRNA|snRNA|scRNA|Y-RNA|SRP", rownames(temp), value = TRUE)
  
  
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

saveRDS(pb.lcpm, paste0(proc.omics.path, "pb_lcpm_v2.rds"))

write.csv(pb.lcpm$Alpha, paste0(proc.omics.path, "proc_pbrna_Alpha_v2.csv"))
write.csv(pb.lcpm$Beta, paste0(proc.omics.path, "proc_pbrna_Beta_v2.csv"))

write.csv(pb.lcpm$Alpha, paste0(proc.omics.path, "proc_pbrna_Delta_v2.csv"))
write.csv(pb.lcpm$Beta, paste0(proc.omics.path, "proc_pbrna_PP_v2.csv"))

myDb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_omics_v2.sqlite"))

for(i in c(1:length(types))){
  temp <- data.frame(pb.lcpm[[types[i]]])
  anno <- gene_anno[match(rownames(temp) ,gene_anno$input_symbol),c("input_symbol","entrez_selected","ensembl_gene_id","description")]
 names(anno) <- c("symbol","gene_id","ensembl","name")
 temp <- cbind(anno,temp)
 rownames(temp) <- NULL
 dbWriteTable(myDb, paste0("proc_pbrna_",types[i]),temp, overwrite = TRUE)
 
}

dbDisconnect(myDb)


