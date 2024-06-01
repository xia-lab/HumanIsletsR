# Count sig features for each omics ~ outcome relationship (SI Table 2)
# Author: Jessica Ewald

## Set your working directory to the "6_for_manuscript" directory

library(dplyr)
library(RSQLite)
library(data.table)

source("../set_paths.R")
setPaths()

##### First, re-run scripts from "1_omic_meta_associations.R" but without covariate adjustment

raw.omics.path <- paste0(other.tables.path, "omics_processing_input/raw/")
proc.omics.path <- paste0(other.tables.path, "omics_processing_input/proc/")

# use same functions as used by web-tool
source(paste0(r.path, "/humanislets_statistics.R")) # put path to your local HumanIsletsR directory

# read in metadata
mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
meta <- dbReadTable(mydb, "proc_metadata")
meta.info <- dbReadTable(mydb, "proc_variable_summary")
dbDisconnect(mydb)

# split by variable type, also, all meta.ephys variables are continuous
meta.cont <- meta.info[meta.info$type == "cont", ]
meta.cont <- meta.cont[meta.cont$group_id != "ephys_cell" & meta.cont$group_id != "ephys_donor", ]
meta.disc <- meta.info[meta.info$type == "disc", ]
meta.ephys.donor <- meta.info[meta.info$group_id == "ephys_donor", ]
meta.ephys.cell <- meta.info[meta.info$group_id == "ephys_cell", ]

# Will need to consider discrete metadata differently - must specify contrasts. 
# For some, it makes more sense to do ANOVA because there is no obvious contrast.

# donorsex: F v M
# diagnosis: T2D v ND, T1D v ND
# diagnosis_computed: PD vs ND, T2D vs ND, T1D vs ND
# hla_a2: Positive vs Negative
# fatinfiltration: Yes vs No

# All others ANOVA: donationtype, pancreasconsistency, collagenase, collagenasetype, super_pop, pop

# read in all datasets
omics <- c("proc_nanostring_merge", "proc_pbrna_Alpha", "proc_pbrna_Beta", "proc_prot", "proc_rnaseq")

# Continuous metadata
dea.results <- data.frame()
for(i in c(1:length(omics))){
  for(k in c(1:dim(meta.cont)[1])){
    
    meta.temp <- meta.cont[k, ]
    print(paste0(omics[i], ", ", meta.temp$column))
    res <- DonorRegression(varGroup = meta.temp$group_id, analysisVar = meta.temp$column, omicsType = omics[i],
                           fixedEffects="NA", mode="tool")
    
    if(grepl("RES-OK", res)){
      res <- read.csv("dea_results.csv")
      res$variable <- meta.temp$column
      res$var_group <- meta.temp$group_id
      res$omics_type <- omics[i]
      colnames(res)[2] <- "Coefficient"
      
      dea.results <- rbind(dea.results, res)
    }
  }
}

dea.results$contrast <- NA
saveRDS(dea.results, "./deg_count_results/dea_results.rds")


## Discrete metadata with obvious contrast choices 
varGroups <- c("donor", "donor", "donor", "donor", "donor", "donor", "donor", "isolation")
varNames <- c("donorsex", "diagnosis", "diagnosis", "diagnosis_computed", "diagnosis_computed", "diagnosis_computed",
              "hla_a2", "fatinfiltration")
refs <- c("Female", "None", "None", "None", "None", "None", "Negative", "No")
contrasts <- c("Male", "Type2", "Type1", "Pre.T2D", "Type2", "Type1", "Positive", "Yes")

disc.contrast.results <- data.frame()
for(i in c(1:length(omics))){
  for(k in c(1:length(varNames))){
    print(paste0(omics[i], ", ", varNames[k]))
    
    res <- DonorRegression(varGroup = varGroups[k], analysisVar = varNames[k], 
                           ref = refs[k], contrast = contrasts[k], omicsType = omics[i], 
                           fixedEffects="NA", mode="tool")
    
    if(grepl("RES-OK", res)){
      res <- read.csv("dea_results.csv")
      res$variable <- varNames[k]
      res$var_group <- varGroups[k]
      res$omics_type <- omics[i]
      res$contrast <- paste0(contrasts[k], "-", refs[k])
      colnames(res)[2] <- "Coefficient"
      
      disc.contrast.results <- rbind(disc.contrast.results, res)
    }
    
  }
}
saveRDS(disc.contrast.results, "./deg_count_results/disc_contrast_results.rds")


## Discrete metadata with ANOVA stats
meta.anova <- meta.disc[!(meta.disc$column %in% varNames), ]

anova.refs <- c()
for(i in c(1:dim(meta.anova)[1])){
  meta.var <- meta.anova$column[i]
  freq <- table(meta[,meta.var]) %>% as.data.frame()
  freq <- freq[order(freq$Freq, decreasing = TRUE), ]
  anova.refs <- c(anova.refs, as.character(freq$Var1[1]))
}

disc.anova.results <- data.frame()
for(i in c(1:length(omics))){
  for(k in c(1:dim(meta.anova)[1])){
    
    meta.temp <- meta.anova[k, ]
    print(paste0(omics[i], ", ", anova.refs[k]))
    
    res <- DonorRegression(varGroup = meta.temp$group_id, analysisVar = meta.temp$column, 
                           ref = anova.refs[k], contrast = "anova", omicsType = omics[i], 
                           fixedEffects="NA", mode="tool")
    
    if(grepl("RES-OK", res)){
      res <- read.csv("dea_results.csv")
      res$variable <- meta.temp$column
      res$var_group <- meta.temp$group_id
      res$omics_type <- omics[i]
      res$contrast <- "ANOVA"
      colnames(res)[2] <- "Coefficient"
      
      disc.anova.results <- rbind(disc.anova.results, res)
    }
    
  }
}
disc.anova.results$Coefficient <- NA
disc.anova.results$sig[disc.anova.results$sig == "up"] <- "sig"
saveRDS(disc.anova.results, "./deg_count_results/disc_anova_results.rds")

# Ephys (donor-level) continuous outcomes
cells <- c("Alpha", "Beta")
gluc.concs <- c("1", "5", "10")
ephys.donor.results <- data.frame()

for(i in c(1:length(omics))){
  for(k in c(1:dim(meta.ephys.donor)[1])){
    for(j in c(1:length(cells))){
      for(l in c(1:length(gluc.concs))){
        print(paste0(omics[i], ", ", k, ", ", j, ", ", l))
        
        meta.temp <- meta.ephys.donor[k, ]
        res <- DonorRegression(varGroup = meta.temp$group_id, analysisVar = meta.temp$column, fixedEffects = "NA",
                               cell = cells[j], glucose = gluc.concs[l], omicsType = omics[i], mode="tool")
        
        if(grepl("RES-OK", res)){
          res <- read.csv("dea_results.csv")
          res$variable <- paste0(meta.temp$column, "_", cells[j], "_", gluc.concs[l], "mM")
          res$var_group <- meta.temp$group_id
          res$omics_type <- paste0(omics[i])
          colnames(res)[2] <- "Coefficient"
          
          ephys.donor.results <- rbind(ephys.donor.results, res)
        }
      }
    }
  }
}

ephys.donor.results$contrast <- NA
saveRDS(ephys.donor.results, "./deg_count_results/ephys_donor_results.rds")

# Ephys outcomes - single cell
ephys.cell.results <- data.frame()
for(k in c(7:dim(meta.ephys.cell)[1])){
  for(j in c(1:length(cells))){
    for(l in c(1:length(gluc.concs))){
      print(paste0(k, ", ", j, ", ", l))
      
      meta.temp <- meta.ephys.cell[k, ]
      
      tryCatch({
        res <- PatchseqSpearman(analysisVar = meta.temp$column, cell = cells[j], glucose = gluc.concs[l], mode="tool")
        
        if(grepl("RES-OK", res)){
          res <- read.csv("dea_results.csv")
          res$variable <- paste0(meta.temp$column, "_", cells[j], "_", gluc.concs[l], "mM")
          res$var_group <- meta.temp$group_id
          res$omics_type <- paste0("proc_scrna_", cells[j], "_", gluc.concs[l], "mM")
          colnames(res)[2] <- "Coefficient"
          
          ephys.cell.results <- rbind(ephys.cell.results, res)
          print(paste0(meta.temp$column, "_", cells[j], "_", gluc.concs[l], "mM"))
        }
      }, error = function(e) {
        print("An error occurred")
      })
    }
  }
}

ephys.cell.results$contrast <- NA
ephys.cell.results$T_statistic <- NA
ephys.cell.results$Associated.proportions <- NA
ephys.cell.results <- ephys.cell.results[,c(1:3, 14, 4:5, 15, 6:13)]
saveRDS(ephys.cell.results, "./deg_count_results/ephys_cell_results.rds")

### Combine everything
res.objs <- list(dea.results, disc.contrast.results, disc.anova.results, ephys.donor.results, ephys.cell.results)
all.results <- rbindlist(res.objs)

# Only keep some columns
all.results <- all.results[,c("Feature", "Gene_ID", "omics_type", "variable", "contrast", "Coefficient", "P_value", "Adjusted.p_value",
                              "Description", "var_group", "Average.level", "T_statistic")]

# Give better column names
colnames(all.results) <- c("Feature", "Gene_ID", "Omics_type", "Metadata", "Contrast", "Coefficient", "P__value", "Adjusted_p__value",
                           "Description", "Metadata_group", "Average_level", "Test_statistic")

# add display omics types
all.results$Omics_ID <- all.results$Omics_type
all.results$Omics_type[all.results$Omics_type == "proc_prot"] <- "Bulk protein expression"
all.results$Omics_type[all.results$Omics_type == "proc_rnaseq"] <- "Bulk gene expression (RNA-seq)"
all.results$Omics_type[all.results$Omics_type == "proc_nanostring_merge"] <- "Bulk gene expression (Nanostring)"
all.results$Omics_type[all.results$Omics_type == "proc_pbrna_Alpha"] <- "Pseudobulk gene expression (alpha cells)"
all.results$Omics_type[all.results$Omics_type == "proc_pbrna_Beta"] <- "Pseudobulk gene expression (beta cells)"
all.results$Omics_type[all.results$Omics_type == "proc_scrna_Alpha_1mM"] <- "Single-cell gene expression"
all.results$Omics_type[all.results$Omics_type == "proc_scrna_Alpha_5mM"] <- "Single-cell gene expression"
all.results$Omics_type[all.results$Omics_type == "proc_scrna_Alpha_10mM"] <- "Single-cell gene expression"
all.results$Omics_type[all.results$Omics_type == "proc_scrna_Beta_1mM"] <- "Single-cell gene expression"
all.results$Omics_type[all.results$Omics_type == "proc_scrna_Beta_5mM"] <- "Single-cell gene expression"
all.results$Omics_type[all.results$Omics_type == "proc_scrna_Beta_10mM"] <- "Single-cell gene expression"

# add display variable groups
all.results$Metadata_Group_ID <- all.results$Metadata_group
group.map <- meta.info[,c("group_id", "group")] %>% distinct()
all.results$Metadata_group[all.results$Metadata_group == "donor"] <- "Donor Characteristics"
all.results$Metadata_group[all.results$Metadata_group == "organ"] <- "Organ Characteristics and Processing" 
all.results$Metadata_group[all.results$Metadata_group == "isolation"] <- "Isolation Outcomes"
all.results$Metadata_group[all.results$Metadata_group == "cell_pro"] <- "Cell Type Proportions"
all.results$Metadata_group[all.results$Metadata_group == "cell_culture"] <- "Cell Culture Outcomes"
all.results$Metadata_group[all.results$Metadata_group == "gsis"] <- "Static Insulin Secretion"
all.results$Metadata_group[all.results$Metadata_group == "seahorse"] <- "Islet Oxygen Consumption (Seahorse Assay)"
all.results$Metadata_group[all.results$Metadata_group == "perifusion"] <- "Dynamic Insulin Responses to Macronutrients"
all.results$Metadata_group[all.results$Metadata_group == "ephys_donor"] <- "Single-cell Function"
all.results$Metadata_group[all.results$Metadata_group == "ephys_cell"] <- "Single-cell Function"

# add display variable names
all.results$Metadata_ID <- all.results$Metadata
all.results$Metadata <- NA

# make for loop
meta.vars <- unique(meta.info$column)
for(i in c(1:length(meta.vars))){
  all.results$Metadata[all.results$Metadata_ID == meta.vars[i]] <- meta.info$display[meta.info$column == meta.vars[i]]
}

other.vars <- unique(all.results$Metadata_ID[is.na(all.results$Metadata)])
for(i in c(1:length(other.vars))){
  temp <- other.vars[i]
  if(grepl("Alpha", temp)){
    temp <- strsplit(temp, "_Alpha_")[[1]]
    label <- meta.info$display[meta.info$column == temp[1]]
    label <- paste0(label, ", alpha cells, ", temp[2], " glucose")
  } else {
    temp <- strsplit(temp, "_Beta_")[[1]]
    label <- meta.info$display[meta.info$column == temp[1]]
    label <- paste0(label, ", beta cells, ", temp[2], " glucose")
  }
  all.results$Metadata[all.results$Metadata_ID == other.vars[i]] <- label
}

# add contrast info
all.results$Metadata[!is.na(all.results$Contrast)] <- paste0(all.results$Metadata[!is.na(all.results$Contrast)], ", ", all.results$Contrast[!is.na(all.results$Contrast)])

# remove ephys ~ pseudobulk associations (don't make sense)
library(stringr)
inds1 <- grep("ephys", all.results$Metadata_Group_ID)
inds2 <- grep("pbrna", all.results$Omics_ID)
inds.remove <- intersect(inds1, inds2)

res <- all.results[-c(inds.remove), ]


#### Count DEGs for each association

meta <- unique(res$Metadata)
ephys.vars <- c(meta[grep("alpha cells", meta)], meta[grep("beta cells", meta)])
res = as.data.table(res)

# summarize number of features with FDR < 0.05 for each omics-outcome combination
all.combinations <- distinct(res[, .(Omics_type, Metadata)])
all.combinations$N <- 0
#all.combinations$N[(all.combinations$Omics_type == "Single-cell gene expression") & (all.combinations$Metadata %in% ephys.vars)] <- NA
all.combinations$unique <- paste0(all.combinations$Omics_type, "_", all.combinations$Metadata)

res <- res[Adjusted_p__value < 0.05]
res <- res[, .N, by = .(Omics_type, Metadata)]

# remove all.combinations not in res
unique.res <- paste0(res$Omics_type, "_", res$Metadata) %>% unique()
all.combinations <- all.combinations[!(unique %in% unique.res)]

# concatenate results
all.combinations <- all.combinations[,.(Omics_type, Metadata, N)]
res <- rbind(res, all.combinations)

# summarize deg counts into wide format
deg.summary <- dcast(res, Metadata ~ Omics_type, value.var = "N")
deg.summary <- as.data.frame(deg.summary)
deg.summary <- deg.summary[match(meta, deg.summary$Metadata), ]

write.csv(deg.summary, "./si_tables/SI_Table_2.csv", row.names = F)
