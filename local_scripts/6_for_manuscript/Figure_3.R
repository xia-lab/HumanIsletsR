# Manuscript Figure 3
# Author: Jessica Ewald

## Set your working directory to the "6_for_manuscript" directory

library(dplyr)
library(RSQLite)
library(ggplot2)

source("../set_paths.R")
source("/Users/jessicaewald/NetbeansProjects/restxialab/src/main/webapp/resources/rscripts/humanislets_statistics.R")
setPaths()

# get data
mydb <- dbConnect(RSQLite::SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
donor <- dbReadTable(mydb, "donor")
isolation <- dbReadTable(mydb, "isolation")
dbDisconnect(mydb)


###### Culture time #######
# culture time
pdf(file="./figures/fig3A.pdf", width=6, height=4)
ggplot(donor, aes(x = predistributionculturetime)) +
  geom_histogram() +
  theme_classic() +
  xlab("Culture time (h)") +
  ylab("Count")
dev.off()

# % purity
pdf(file="./figures/fig3B.pdf", width=6, height=4)
ggplot(isolation, aes(x = puritypercentage)) +
  geom_histogram() +
  theme_classic() +
  xlab("Islet purity (%)") +
  ylab("Count")
dev.off()

###### Figure 2A #######
# This figure shows the distribution of the main metadata across all donors

mydb <- dbConnect(RSQLite::SQLite(), paste0(sqlite.path, "HI_omics.sqlite"))
rnaseq <- dbReadTable(mydb, "proc_rnaseq")
dbDisconnect(mydb)

# perform analysis
t2d_nocovar = DonorRegression(varGroup = "donor", 
                              analysisVar = "diagnosis",
                              ref = "None",
                              contrast = "Type2",
                              fixedEffects = "NA",
                              omicsType = "proc_rnaseq",
                              pvalThresh = 0.05,
                              mode = "tool")
t2d_nocovar <- read.csv("dea_results.csv")


t2d_culturetime = DonorRegression(varGroup = "donor", 
                                  analysisVar = "diagnosis",
                                  ref = "None",
                                  contrast = "Type2",
                                  fixedEffects = "predistributionculturetime",
                                  omicsType = "proc_rnaseq",
                                  pvalThresh = 0.05,
                                  mode = "tool")
t2d_culturetime <- read.csv("dea_results.csv")


t2d_purity = DonorRegression(varGroup = "donor", 
                             analysisVar = "diagnosis",
                             ref = "None",
                             contrast = "Type2",
                             fixedEffects = "puritypercentage",
                             omicsType = "proc_rnaseq",
                             pvalThresh = 0.05,
                             mode = "tool")
t2d_purity <- read.csv("dea_results.csv")

t2d_nonendo = DonorRegression(varGroup = "donor", 
                              analysisVar = "diagnosis",
                              ref = "None",
                              contrast = "Type2",
                              fixedEffects = "exo_per",
                              omicsType = "proc_rnaseq",
                              pvalThresh = 0.05,
                              mode = "tool")
t2d_nonendo <- read.csv("dea_results.csv")


# Extract columns of interest
res_nocovar <- t2d_nocovar[,c("Gene_ID", "negLogPval", "P_value", "Log2FC", "Adjusted.p_value")]
colnames(res_nocovar)[2:5] <- paste0("None_", colnames(res_nocovar)[2:5])

res_culturetime <- t2d_culturetime[,c("Gene_ID", "negLogPval", "P_value", "Log2FC", "Adjusted.p_value")]
colnames(res_culturetime)[2:5] <- paste0("CT_", colnames(res_culturetime)[2:5])

res_purity <- t2d_purity[,c("Gene_ID", "negLogPval", "P_value", "Log2FC", "Adjusted.p_value")]
colnames(res_purity)[2:5] <- paste0("Purity_", colnames(res_purity)[2:5])

res_nonendo <- t2d_nonendo[,c("Gene_ID", "negLogPval", "P_value", "Log2FC", "Adjusted.p_value")]
colnames(res_nonendo)[2:5] <- paste0("Nonendo_", colnames(res_nonendo)[2:5])


# Merge results
dea_res <- merge(res_nocovar, res_culturetime, by = "Gene_ID", all = TRUE)
dea_res <- merge(dea_res, res_purity, by = "Gene_ID", all = TRUE)
dea_res <- merge(dea_res, res_nonendo, by = "Gene_ID", all = TRUE)

# Get consensus significance label
dea_res$Sig_None_CT <- NA
dea_res$Sig_None_CT[dea_res$None_Adjusted.p_value < 0.05 & dea_res$CT_Adjusted.p_value < 0.05] <- "Both"
dea_res$Sig_None_CT[dea_res$None_Adjusted.p_value > 0.05 & dea_res$CT_Adjusted.p_value > 0.05] <- "Neither"
dea_res$Sig_None_CT[dea_res$None_Adjusted.p_value < 0.05 & dea_res$CT_Adjusted.p_value > 0.05] <- "Without covariate"
dea_res$Sig_None_CT[dea_res$None_Adjusted.p_value > 0.05 & dea_res$CT_Adjusted.p_value < 0.05] <- "With covariate"

dea_res$Sig_None_Purity <- NA
dea_res$Sig_None_Purity[dea_res$None_Adjusted.p_value < 0.05 & dea_res$Purity_Adjusted.p_value < 0.05] <- "Both"
dea_res$Sig_None_Purity[dea_res$None_Adjusted.p_value > 0.05 & dea_res$Purity_Adjusted.p_value > 0.05] <- "Neither"
dea_res$Sig_None_Purity[dea_res$None_Adjusted.p_value < 0.05 & dea_res$Purity_Adjusted.p_value > 0.05] <- "Without covariates"
dea_res$Sig_None_Purity[dea_res$None_Adjusted.p_value > 0.05 & dea_res$Purity_Adjusted.p_value < 0.05] <- "With % purity"

dea_res$Sig_None_Nonendo <- NA
dea_res$Sig_None_Nonendo[dea_res$None_Adjusted.p_value < 0.05 & dea_res$Nonendo_Adjusted.p_value < 0.05] <- "Both"
dea_res$Sig_None_Nonendo[dea_res$None_Adjusted.p_value > 0.05 & dea_res$Nonendo_Adjusted.p_value > 0.05] <- "Neither"
dea_res$Sig_None_Nonendo[dea_res$None_Adjusted.p_value < 0.05 & dea_res$Nonendo_Adjusted.p_value > 0.05] <- "Without covariates"
dea_res$Sig_None_Nonendo[dea_res$None_Adjusted.p_value > 0.05 & dea_res$Nonendo_Adjusted.p_value < 0.05] <- "With % purity"

# Plot results
pdf(file="./figures/fig3F_legend.pdf", width=5, height=5)
ggplot(dea_res, aes(x = None_negLogPval, y = CT_negLogPval, color = Sig_None_CT)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "black", lwd = 1) +
  theme_bw() +
  scale_color_manual(values = c("#636363", "#cccccc", "red", "blue")) + 
  xlab("-log10(p-value): no covariate adjustment") +
  ylab("-log10(p-value): adjusting for culture time")
dev.off()

pdf(file="./figures/fig3F.pdf", width=5, height=5)
ggplot(dea_res, aes(x = None_negLogPval, y = CT_negLogPval, color = Sig_None_CT)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "black", lwd = 1) +
  theme_bw() +
  scale_color_manual(values = c("#636363", "#cccccc", "red", "blue")) + 
  xlab("-log10(p-value): no covariate adjustment") +
  ylab("-log10(p-value): adjusting for culture time") +
  theme(legend.position = "None")
dev.off()

ggplot(dea_res, aes(x = None_negLogPval, y = Purity_negLogPval, color = Sig_None_Purity)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "black", lwd = 1) +
  theme_bw() +
  scale_color_manual(values = c("#636363", "#cccccc", "blue")) + 
  xlab("-log10(p-value): no covariate adjustment") +
  ylab("-log10(p-value): adjusting for % purity") +
  theme(legend.position = "None")


pdf(file="./figures/fig3G.pdf", width=5, height=5)
ggplot(dea_res, aes(x = None_negLogPval, y = Nonendo_negLogPval, color = Sig_None_Nonendo)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "black", lwd = 1) +
  theme_bw() +
  scale_color_manual(values = c("#636363", "#cccccc", "red", "blue")) + 
  xlab("-log10(p-value): no covariate adjustment") +
  ylab("-log10(p-value): adjusting for % non-endocrine") +
  theme(legend.position = "None")
dev.off()
  