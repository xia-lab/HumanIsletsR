# Manuscript Figure 2
# Author: Jessica Ewald

## Set your working directory to the "6_for_manuscript" directory

library(dplyr)
library(RSQLite)
library(pheatmap)
library(ggplot2)

source("../set_paths.R")
setPaths()


###### Figure 2A #######
# This figure shows the distribution of the main metadata across all donors

# get data
mydb <- dbConnect(RSQLite::SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
donor <- dbReadTable(mydb, "donor")
dbDisconnect(mydb)

donor$diagnosis_computed[donor$diagnosis_computed == "Type1"] <- "Type 1"
donor$diagnosis_computed[donor$diagnosis_computed == "Type2"] <- "Type 2"
donor$diagnosis_computed[donor$diagnosis_computed == "Pre.T2D"] <- "Pre-diabetes"

# sex
pdf(file="./figures/fig2A.pdf", width=3, height=3.8)
ggplot(donor, aes(x = donorsex)) +
  geom_bar(stat="count") +
  theme_classic() +
  xlab("Sex") +
  ylab("Count") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=1))
dev.off()

# age
pdf(file="./figures/fig2B.pdf", width=4, height=4)
ggplot(donor, aes(x = donorage)) +
  geom_histogram() +
  theme_classic() +
  xlab("Age") +
  ylab("Count")
dev.off()

# diagnosis
pdf(file="./figures/fig2C.pdf", width=5.5, height=4)
ggplot(donor, aes(x = diagnosis_computed)) +
  geom_bar(stat="count") +
  theme_classic() +
  xlab("Diabetes diagnosis") +
  ylab("Count") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

# hba1c
pdf(file="./figures/fig2D.pdf", width=4, height=4)
ggplot(donor, aes(x = hba1c)) +
  geom_histogram() +
  theme_classic() +
  xlab("HbA1c (%)") +
  ylab("Count")
dev.off()

## heatmap of data availability
mydb <- dbConnect(RSQLite::SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
donors <- dbReadTable(mydb, "donor")
donors <- donors$record_id
dbDisconnect(mydb)

avail_df <- data.frame(donor = rep(1, length(donors)),
                       tech = rep(1, length(donors)),
                       gsis = rep(0, length(donors)),
                       peri = rep(0, length(donors)),
                       seahorse = rep(0, length(donors)),
                       ephys = rep(0, length(donors)),
                       prot = rep(0, length(donors)),
                       nano = rep(0, length(donors)),
                       bulkrna = rep(0, length(donors)),
                       pbrna = rep(0, length(donors)),
                       scrna = rep(0, length(donors)))

# get lists of record_ids in each outcomes table
mydb <- dbConnect(RSQLite::SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
gsis <- dbGetQuery(mydb, "SELECT record_id FROM gsis")[,1] %>% unique() 
peri <- dbGetQuery(mydb, "SELECT record_id FROM function_summary")[,1] %>% unique()
seahorse <- dbGetQuery(mydb, "SELECT record_id FROM seahorse")[,1] %>% unique()
ephys <- dbGetQuery(mydb, "SELECT record_id FROM ephys_donor")[,1] %>% unique()
scrna <- dbGetQuery(mydb, "SELECT record_id FROM ephys_cell")[,1] %>% unique()
dbDisconnect(mydb)

# get lists of record_ids in each bulk omics table
mydb <- dbConnect(RSQLite::SQLite(), paste0(sqlite.path, "HI_omics.sqlite"))
nano <- dbListFields(mydb, "proc_nanostring_merge")
bulkrna <- dbListFields(mydb, "proc_rnaseq")
pbrna <- c(dbListFields(mydb, "proc_pbrna_Alpha"), dbListFields(mydb, "proc_pbrna_Beta")) %>% unique()
prot <- dbListFields(mydb, "proc_prot")
dbDisconnect(mydb)

# update data avail matrix
avail_df$gsis[donors %in% gsis] <- 1
avail_df$peri[donors %in% peri] <- 1
avail_df$seahorse[donors %in% seahorse] <- 1
avail_df$ephys[donors %in% ephys] <- 1
avail_df$prot[donors %in% prot] <- 1
avail_df$nano[donors %in% nano] <- 1
avail_df$bulkrna[donors %in% bulkrna] <- 1
avail_df$pbrna[donors %in% pbrna] <- 1
avail_df$scrna[donors %in% scrna] <- 1

disp.colnames <- c("Clinical metadata", "Technical metadata", "Static insulin secretion", "Dynamic insulin secretion",
                   "Mitochondrial function", "Electrophysiology", "Proteomics (bulk)", "Nanostring (bulk)", "RNAseq (bulk)",
                   "RNAseq (pseudobulk)", "RNAseq (single-cell)")

pheatmap(avail_df, legend = FALSE, show_rownames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE,
         color = c("white", "#43A047"), angle_col=90, labels_col = disp.colnames, border_color = NA,
         filename="./figures/fig2E.pdf", width=4.5, height=9)





