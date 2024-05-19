# Process Nanostring data for webtool
# Author: Jessica Ewald

## Set your working directory to the "1_process_bulk_omics" directory

library(ggplot2)
library(dplyr)
library(sva)

source("../set_paths.R")
setPaths()

raw.omics.path <- paste0(other.tables.path, "omics_processing_input/unproc/")
proc.omics.path <- paste0(other.tables.path, "omics_processing_input/proc/")
if(!dir.exists(proc.omics.path)){ dir.create(proc.omics.path) }

# read in two codesets
cs1 <- read.table(paste0(raw.omics.path, "unproc_nanostring_C6555.txt"),
                  header = T, sep = "\t")
cs1 <- cs1[cs1$X != '', ]
rownames(cs1) <- toupper(cs1$X)
cs1 <- cs1[,-1]
cs1 <- as.matrix(cs1)
mode(cs1) <- "numeric"

cs2 <- read.table(paste0(raw.omics.path, "unproc_nanostring_C8898.txt"),
                  header = T, sep = "\t")
cs2 <- cs2[cs2$X != '', ]
rownames(cs2) <- toupper(cs2$X)
cs2 <- cs2[,-1]
cs2 <- as.matrix(cs2)
mode(cs2) <- "numeric"

# merge together and see
cs.merge <- merge(cs1, cs2, by = "row.names", all = TRUE)
rownames(cs.merge) <- toupper(cs.merge$Row.names)
cs.merge <- cs.merge[,-1]
cs.merge <- as.matrix(cs.merge)
mode(cs.merge) <- "numeric"


# normalize
normNanostring <- function(csMat, nBox){
  rand.inds <- sample(1:dim(csMat)[2], nBox)
  boxplot(csMat[,rand.inds], main = "Raw")
  
  csMat<- apply(csMat, 2, function(x){x/median(x, na.rm=T)});
  boxplot(csMat[,rand.inds], main = "Sample Median")
  
  csMat[!is.na(csMat)] <- log2(csMat[!is.na(csMat)])
  boxplot(csMat[,rand.inds], main = "Log transformed")
  
  return(csMat)
}

cs1 <- normNanostring(cs1, 25)
cs2 <- normNanostring(cs2, 25)
cs.merge <- normNanostring(cs.merge, 25)

# get metadata


# look at PCA
cs1.pr <- prcomp(t(cs1), scale = TRUE)
cs1.x <- cs1.pr$x %>% as.data.frame()

cs2.pr <- prcomp(t(cs2), scale = TRUE)
cs2.x <- cs2.pr$x %>% as.data.frame()

cs.merge.pr <- prcomp(t(na.omit(cs.merge)), scale = TRUE)
cs.merge.x <- cs.merge.pr$x %>% as.data.frame()
cs.merge.x$batch <- c(rep(NA, dim(cs.merge.x)[1]))
cs.merge.x$batch[rownames(cs.merge.x) %in% colnames(cs1)] <- "CS1"
cs.merge.x$batch[rownames(cs.merge.x) %in% colnames(cs2)] <- "CS2"

ggplot(cs1.x, aes(x = PC1, y = PC2)) +
  geom_point() +
  xlab(paste0("PC1 (", round((cs1.pr$sdev[1]^2)/sum(cs1.pr$sdev^2)*100,1), "%)")) +
  ylab(paste0("PC2 (", round((cs1.pr$sdev[2]^2)/sum(cs1.pr$sdev^2)*100,1), "%)")) +
  theme_bw() +
  ggtitle("CS1")

ggplot(cs2.x, aes(x = PC1, y = PC2)) +
  geom_point() +
  xlab(paste0("PC1 (", round((cs2.pr$sdev[1]^2)/sum(cs2.pr$sdev^2)*100,1), "%)")) +
  ylab(paste0("PC2 (", round((cs2.pr$sdev[2]^2)/sum(cs2.pr$sdev^2)*100,1), "%)")) +
  theme_bw() +
  ggtitle("CS2")

ggplot(cs.merge.x, aes(x = PC1, y = PC2, color = batch)) +
  geom_point() +
  xlab(paste0("PC1 (", round((cs.merge.pr$sdev[1]^2)/sum(cs.merge.pr$sdev^2)*100,1), "%)")) +
  ylab(paste0("PC2 (", round((cs.merge.pr$sdev[2]^2)/sum(cs.merge.pr$sdev^2)*100,1), "%)")) +
  theme_bw() +
  ggtitle("Merged CS")

# perform batch effect correction of merged data
batch <- data.frame(batch = cs.merge.x$batch)
mod.combat <- model.matrix(~1, data = batch)
# ComBat only works if there are values for both batches. Only batch correct merged data and then recombine.
cs.combat <- ComBat(dat = na.omit(cs.merge), batch = batch$batch, mod = mod.combat) 
cs.combat <- rbind(cs.combat, cs.merge[!(rownames(cs.merge) %in% rownames(cs.combat)), ])

cs.combat.pr <- prcomp(t(na.omit(cs.combat)), scale = TRUE)
cs.combat.x <- cs.combat.pr$x %>% as.data.frame()
cs.combat.x$batch <- c(rep(NA, dim(cs.combat.x)[1]))
cs.combat.x$batch[rownames(cs.combat.x) %in% colnames(cs1)] <- "CS1"
cs.combat.x$batch[rownames(cs.combat.x) %in% colnames(cs2)] <- "CS2"

ggplot(cs.combat.x, aes(x = PC1, y = PC2, color = batch)) +
  geom_point() +
  xlab(paste0("PC1 (", round((cs.combat.pr$sdev[1]^2)/sum(cs.combat.pr$sdev^2)*100,1), "%)")) +
  ylab(paste0("PC2 (", round((cs.combat.pr$sdev[2]^2)/sum(cs.combat.pr$sdev^2)*100,1), "%)")) +
  theme_bw() +
  ggtitle("Merged & Batch-corrected CS")


# write out results
write.csv(cs1, paste0(proc.omics.path, "proc_nanostring_C6555.csv"))
write.csv(cs2, paste0(proc.omics.path, "proc_nanostring_C8898.csv"))
write.csv(cs.combat, paste0(proc.omics.path, "proc_nanostring_merge.csv"))



