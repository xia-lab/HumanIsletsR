# Process raw RNA-seq counts for web-tool (v2)
# Based on 1_proc_rnaseq.R (Jessica Ewald)
#
# Changes from v1:
#  - no Entrez ID transfer: Ensembl ID (accession) is the only key. Entrez is not used as a key
#    or intermediate at any step; genes are never mapped, merged or dropped by Entrez ID.
#    Entrez ID (gene_id), symbol, name and biotype are annotation columns matched directly by
#    Ensembl ID from libraries/gene_anno_full.qs, and may be empty
#  - expression filter (edgeR filterByExpr: >= 10 reads, scaled to library size,
#    in >= 20% of donors) replaces the zero filter (>= 1 read in >= 20% of donors)
#  - no variance filter: it removed stably expressed genes such as TCF7L2.
#    Variance-based selection belongs inside analyses that need it (PCA, MOFA, WGCNA, DIABLO)
#  - all batches kept; batch removed with ComBat on RLE log-CPM, protecting diagnosis and sex.
#    Both are unevenly spread across batches, and ComBat_seq without protection shrank their
#    real effects (T2D ~8%, sex ~12%); the protected log-CPM ComBat matched a batch-covariate
#    model without inflating false positives
#  - output is in the same format as the proc_rnaseq / unproc_rnaseq tables in HI_omics_v2.sqlite;
#    unproc keeps the raw (uncorrected) counts
#
# Order: raw counts -> expression filter -> RLE + log-CPM -> ComBat (diagnosis + sex protected) -> annotation

## Set your working directory to the "1_process_bulk_omics" directory

library(data.table)
library(edgeR)
library(sva)
library(qs)
library(RSQLite)

source("../set_paths.R")
setPaths()

raw.omics.path <- paste0(other.tables.path, "omics_processing_input/unproc/")
proc.omics.path <- paste0(other.tables.path, "omics_processing_input/proc/")

# back up an existing output file before it is overwritten
backupIfExists <- function(f){
  if(file.exists(f)){
    bak <- paste0(f, ".bak_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    stopifnot(file.copy(f, bak))
    message("Backed up ", f, " -> ", bak)
  }
}

# read in raw counts (Ensembl IDs, one column per donor)
fileName <- paste0(other.tables.path, "updated_data/From_Anna_group/HumanIslets_OxfordStanfordCombined_Edmonton_GeneCounts.txt")
orig.counts <- fread(fileName, header = TRUE, data.table = FALSE)
counts <- as.matrix(orig.counts[, -(1:2)])
rownames(counts) <- orig.counts$GeneID
stopifnot(!anyDuplicated(rownames(counts)), !anyNA(counts))

# batch info (same content as RNAseqEdmontonSampleBatch.xlsx)
batch <- read.table(paste0(raw.omics.path, "rnaseq_batch_info_v2.txt"), sep = "\t", header = TRUE)
batch <- batch[match(colnames(counts), batch$record_id), ]
stopifnot(identical(batch$record_id, colnames(counts)))

# donor diagnosis and sex (protected during batch correction)
myDb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"), flags = SQLITE_RO)
meta <- dbGetQuery(myDb, "SELECT record_id, diagnosis, donorsex FROM proc_metadata")
dbDisconnect(myDb)
meta <- meta[match(colnames(counts), meta$record_id), ]
stopifnot(identical(meta$record_id, colnames(counts)), !anyNA(meta$diagnosis), !anyNA(meta$donorsex))

# expression filter: >= 10 reads (CPM at the median library size) in >= 20% of donors
keep <- filterByExpr(counts, min.count = 10, large.n = 0, min.prop = 0.2)
counts <- counts[keep, ]

# normalize for sequencing depth
nf <- calcNormFactors(counts, method = "RLE")
lcpm <- cpm(counts, lib.size = colSums(counts)*nf, log = TRUE)

# adjust for batch effect on the log-CPM scale, protecting diagnosis and sex
mod <- model.matrix(~ diagnosis + donorsex, data = meta)
lcpm <- ComBat(lcpm, batch = batch$batch, mod = mod, par.prior = TRUE)

# annotate from gene_anno_full.qs (annotation only, no Entrez ID transfer)
# the primary Ensembl-based entry is listed first; later duplicates of an Ensembl ID are
# alias look-ups (old symbols), so the first entry per Ensembl ID is kept
anno <- qread(paste0(other.tables.path, "libraries/gene_anno_full.qs"))
anno <- anno[!is.na(anno$ensembl_gene_id) & !duplicated(anno$ensembl_gene_id), ]
anno <- anno[match(rownames(lcpm), anno$ensembl_gene_id), ]
feature.info <- data.frame(accession = rownames(lcpm),
                           gene_id = as.character(anno$entrez_selected),
                           symbol = anno$hgnc_symbol,
                           name = trimws(sub("\\s*\\[Source:.*\\]$", "", anno$description)),
                           biotype = anno$gene_biotype)

# keep protein-coding and lncRNA genes, as in the web tool
keep.bt <- feature.info$biotype %in% c("protein_coding", "lncRNA")
proc <- cbind(feature.info[keep.bt, ], as.data.frame(lcpm[keep.bt, ]))
unproc <- cbind(data.frame(feature_id = rownames(counts)[keep.bt]), as.data.frame(counts[keep.bt, ]))
rownames(proc) <- rownames(unproc) <- NULL

# write out files (same columns as proc_rnaseq / unproc_rnaseq in HI_omics_v2.sqlite; NA written as empty)
proc.file <- paste0(proc.omics.path, "proc_rnaseq_v2.csv")
unproc.file <- paste0(raw.omics.path, "unproc_rnaseq_v2.csv")
backupIfExists(proc.file)
backupIfExists(unproc.file)
fwrite(proc, proc.file)
fwrite(unproc, unproc.file)
