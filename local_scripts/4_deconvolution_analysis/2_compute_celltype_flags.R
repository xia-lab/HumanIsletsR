# Compute the endocrine and non-endocrine cell type flags
# Jessica Ewald

# Set working directory to "4_deconvolution_analysis"

library(dplyr)
library(RSQLite)

source("../set_paths.R")
source("/Users/jessicaewald/NetbeansProjects/omicsquareapi/src/main/webapp/resources/rscripts/humanislets_statistics.R")
setPaths()


# Make results directory
results.path <- "./cell_proportion_associations/"
dir.create("./cell_proportion_associations/")

cell.vars <- c("exo_per", "alpha_end", "beta_end", "delta_end", "gamma_end")
omics.types <- c("proc_prot", "proc_rnaseq", "proc_nanostring_merge")

# perform dea
for(i in c(1:length(omics.types))){
  new.dir <- paste0(results.path, omics.types[i])
  dir.create(new.dir)
  
  for(k in c(1:length(cell.vars))){
    DonorRegression("cell_pro", cell.vars[k], "NA", "NA", "NA", omics.types[i], 0.05, "all", "NA", "NA", "true")
    file.copy("dea_results.csv", paste0(new.dir, "/", cell.vars[k], ".csv"))
    file.remove("dea_results.csv")
  }
}

cell.types <- c("nonendo", "alpha", "beta", "delta", "gamma")

# read in association results
rna.files <- list.files(paste0(results.path, "proc_rnaseq"),
                        full.names = TRUE)
rna <- lapply(rna.files, read.csv)
file.names <- gsub(paste0(results.path, "proc_rnaseq\\/"), "", rna.files)
file.names <- gsub("\\.csv", "", file.names)
file.names <- gsub("_.*", "", file.names)
file.names[file.names == "exo"] <- "nonendo"
names(rna) <- file.names

# write out non-endocrine signal strength
mydb <- dbConnect(SQLite(), paste0(other.tables.path, "libraries/hsa_genes.sqlite"))
gene.info <- dbReadTable(mydb, "entrez")
dbDisconnect(mydb)
gene.info <- data.frame(Gene_ID = unique(gene.info$gene_id))

# make non-endocrine signal summary based on RNA-seq data (most genes)
nonendo <- rna$nonendo
cellpro_assns <- gene.info
cellpro_assns$`Associated proportions` <- "--"
cellpro_assns$`Associated proportions`[cellpro_assns$Gene_ID %in% 
                                         nonendo$Gene_ID[nonendo$Adjusted.p_value < 0.05 & nonendo$Coefficient > 0]] <- "Non-endocrine"
cellpro_assns$`Associated proportions`[cellpro_assns$Gene_ID %in% 
                                         nonendo$Gene_ID[nonendo$Adjusted.p_value < 0.05 & nonendo$Coefficient < 0]] <- "Endocrine"
cellpro_assns$`Associated proportions`[cellpro_assns$Gene_ID %in% 
                                         nonendo$Gene_ID[nonendo$Adjusted.p_value >= 0.05]] <- "None"
saveRDS(cellpro_assns, paste0(other.tables.path, "analysis_input/cell_signal.rds"))
