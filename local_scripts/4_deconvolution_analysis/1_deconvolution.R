# Deconvolution version 2
# Author: Jessica Ewald

## Set your working directory to the "4_deconvolution_analysis" directory

library(dplyr)
library(stringr)
library(RSQLite)
library(ggplot2)
library(pheatmap)
library(ggpubr)

source("../set_paths.R")
setPaths()

raw.omics.path <- paste0(other.tables.path, "omics_processing_input/raw/")

# read in proteomics data
prot <- read.csv(paste0(raw.omics.path, "raw_prot.csv"),
                 header = TRUE, row.names = 1)

# process uniprot annotation data
uniprot.dat <- read.table("./input_data/uniprot_id_mass.tsv",
                          header = TRUE, sep = "\t", fill = TRUE, quote="")
IDs <- str_split(uniprot.dat$Gene.Names, " ")
names(IDs) <- uniprot.dat$Entry
IDs <- data.table::melt(IDs)

# normalize by sample volume
sample.mass <- apply(prot, 2, sum, na.rm = TRUE)
prot <- sweep(prot, 2, sample.mass, "/")
prot[is.na(prot)] <- min(prot, na.rm = T)/5
  
########## COMPUTE CELL TYPE PROPORTIONS ###########

# purity & cell type prior
max.purity <- 0.98
proB.mean <- 0.55
proA.mean <- 0.30
proD.mean <- 0.10
proG.mean <- 0.05

# Read in and process marker genes
markers <- readRDS("./input_data/markers_herrera.rds")

# get marker genes
bcm.dat <- prot[rownames(prot) %in% markers$uniprot[markers$cell_type == "beta"], ]
acm.dat <- prot[rownames(prot) %in% markers$uniprot[markers$cell_type == "alpha"], ]
dcm.dat <- prot[rownames(prot) %in% markers$uniprot[markers$cell_type == "delta"], ]
gcm.dat <- prot[rownames(prot) %in% markers$uniprot[markers$cell_type == "gamma"], ]

# compute median of each marker
bcell <- apply(bcm.dat, 1, median, na.rm = TRUE) %>% unlist()
bcell <- data.frame(uniprot = names(bcell), median.bulk = unname(bcell))

acell <- apply(acm.dat, 1, median, na.rm = TRUE) %>% unlist()
acell <- data.frame(uniprot = names(acell), median.bulk = unname(acell))

dcell <- apply(dcm.dat, 1, median, na.rm = TRUE) %>% unlist()
dcell <- data.frame(uniprot = names(dcell), median.bulk = unname(dcell))

gcell <- apply(gcm.dat, 1, median, na.rm = TRUE) %>% unlist()
gcell <- data.frame(uniprot = names(gcell), median.bulk = unname(gcell))


# compute estimated median level of marker within cell type
bcell$median.beta <- bcell$median.bulk/proB.mean
acell$median.alpha <- acell$median.bulk/proA.mean
dcell$median.delta <- dcell$median.bulk/proD.mean
gcell$median.gamma <- gcell$median.bulk/proG.mean

# re-order
bcm.dat <- bcm.dat[match(bcell$uniprot, rownames(bcm.dat)), ]
acm.dat <- acm.dat[match(acell$uniprot, rownames(acm.dat)), ]
dcm.dat <- dcm.dat[match(dcell$uniprot, rownames(dcm.dat)), ]
gcm.dat <- gcm.dat[match(gcell$uniprot, rownames(gcm.dat)), ]

# check
identical(bcell$uniprot, rownames(bcm.dat)) & 
  identical(acell$uniprot, rownames(acm.dat)) & 
  identical(dcell$uniprot, rownames(dcm.dat)) & 
  identical(gcell$uniprot, rownames(gcm.dat))

# compute proportion of each cell type in each sample for each marker gene
bcell.pro <- sweep(bcm.dat, 1, bcell$median.beta, "/")
acell.pro <- sweep(acm.dat, 1, acell$median.alpha, "/")
dcell.pro <- sweep(dcm.dat, 1, dcell$median.delta, "/")
gcell.pro <- sweep(gcm.dat, 1, gcell$median.gamma, "/")

# compute median proportion of each cell type in each sample, across all marker gene estimates
bcell.pro.med <- apply(bcell.pro, 2, median, na.rm = TRUE)
acell.pro.med <- apply(acell.pro, 2, median, na.rm = TRUE)
dcell.pro.med <- apply(dcell.pro, 2, median, na.rm = TRUE)
gcell.pro.med <- apply(gcell.pro, 2, median, na.rm = TRUE)

# check
identical(names(bcell.pro.med), names(acell.pro.med)) &
  identical(names(bcell.pro.med), names(dcell.pro.med)) &
  identical(names(bcell.pro.med), names(gcell.pro.med))

# compute proportions
cell.pro <- data.frame(record_id = names(acell.pro.med), beta = bcell.pro.med, alpha = acell.pro.med, 
                       delta = dcell.pro.med, gamma = gcell.pro.med)
cell.pro$end_total <- apply(cell.pro[,-1], 1, sum)

# Assume max possible purity is 98%. Calculate the final "total" from which to subtract everything else.
max.tot.end <- max(cell.pro$end_total)
mass.total <- max.tot.end/max.purity
cell.pro$exo_total <- mass.total - cell.pro$end_total
cell.pro$total <- cell.pro$end_total + cell.pro$exo_total

cell.pro$beta_per <- cell.pro$beta/cell.pro$total
cell.pro$alpha_per <- cell.pro$alpha/cell.pro$total
cell.pro$delta_per <- cell.pro$delta/cell.pro$total
cell.pro$gamma_per <- cell.pro$gamma/cell.pro$total
cell.pro$exo_per <- cell.pro$exo_total/cell.pro$total

cell.pro$beta_end <- cell.pro$beta/cell.pro$end_total
cell.pro$alpha_end <- cell.pro$alpha/cell.pro$end_total
cell.pro$delta_end <- cell.pro$delta/cell.pro$end_total
cell.pro$gamma_end <- cell.pro$gamma/cell.pro$end_total

# Write out results
write.csv(cell.pro[,c("record_id", "beta_end", "alpha_end", "delta_end", "gamma_end", "exo_per")],
          paste0(other.tables.path, "processing_input/composition.csv"),
          row.names = FALSE)

