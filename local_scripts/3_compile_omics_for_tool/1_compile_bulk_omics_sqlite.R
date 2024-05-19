# Compile processed omics data into SQLite tables
# Author: Jessica Ewald

## Set your working directory to the "3_compile_omics_for_tool" directory
# All previous omics processing scripts must have been run, or the processed files will not exist for compilation

library(RSQLite)

source("../set_paths.R")
setPaths()

proc.filepath <- paste0(other.tables.path, "omics_processing_input/proc/")
raw.filepath <- paste0(other.tables.path, "omics_processing_input/unproc/")

proc.names <- c("proc_nanostring_merge.csv", "proc_pbrna_Alpha.csv", 
                "proc_pbrna_Beta.csv", "proc_prot.csv", "proc_rnaseq.csv")

raw.names <- c("unproc_nanostring_C6555.txt", "unproc_nanostring_C8898.txt", 
               "unproc_prot.csv", "unproc_rnaseq.csv", "unproc_pbrna_Alpha.csv", 
               "unproc_pbrna_Beta.csv")

myDb <- dbConnect(SQLite(), paste0(other.tables.path, "libraries/hsa_genes.sqlite"))
entrez <- dbReadTable(myDb, "entrez")
dbDisconnect(myDb)
entrez$name <- gsub('\\"', "", entrez$name)


myDb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_omics.sqlite"))
for(i in c(1:length(proc.names))){
  if(grepl("csv", proc.names[i], fixed = TRUE)){
    dat <- read.csv(paste0(proc.filepath, proc.names[i]))
    colnames(dat)[1] <- "feature_id"
    dat.name <- gsub("\\.csv", "", proc.names[i])
  } else {
    dat <- read.table(paste0(proc.filepath, proc.names[i]), sep="\t", header = TRUE)
    colnames(dat)[1] <- "feature_id"
    dat.name <- gsub("\\.txt", "", proc.names[i])
  }
  
  # get symbols
  feature.vec <- dat$feature_id
  if(i == 1) {
    hit.inx <- match(feature.vec, entrez[, "symbol"])
    gene.ids <- entrez[hit.inx, ]
    na.inx <- which(is.na(gene.ids$gene_id))
    gene.ids$gene_id[na.inx] <- feature.vec[na.inx]
    gene.ids$symbol[na.inx] <- feature.vec[na.inx]
    gene.ids$name[na.inx] <- "--"
    dat <- dat[,-1]
    dat <- cbind(gene.ids, dat)
    rownames(dat) <- NULL
  } else {
    hit.inx <- match(feature.vec, entrez[, "gene_id"])
    gene.ids <- entrez[hit.inx, ]
    na.inx <- which(is.na(gene.ids$gene_id))
    gene.ids$gene_id[na.inx] <- feature.vec[na.inx]
    gene.ids$symbol[na.inx] <- feature.vec[na.inx]
    gene.ids$name[na.inx] <- "--"
    dat <- dat[,-1]
    dat <- cbind(gene.ids, dat)
    rownames(dat) <- NULL
  }
  
  dbWriteTable(myDb, dat.name, dat, overwrite = TRUE)
}

for(i in c(1:length(raw.names))){
  if(grepl("csv", raw.names[i], fixed = TRUE)){
    dat <- read.csv(paste0(raw.filepath, raw.names[i]))
    colnames(dat)[1] <- "feature_id"
    dat.name <- gsub("\\.csv", "", raw.names[i])
  } else {
    dat <- read.table(paste0(raw.filepath, raw.names[i]), sep="\t", header = TRUE)
    colnames(dat)[1] <- "feature_id"
    dat.name <- gsub("\\.txt", "", raw.names[i])
  }
  dbWriteTable(myDb, dat.name, dat, overwrite = TRUE)
}
dbDisconnect(myDb)
