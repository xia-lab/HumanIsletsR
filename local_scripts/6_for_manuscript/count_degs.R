# Count sig features for each omics ~ outcome relationship (SI Table 2)
# Author: Jessica Ewald

## Set your working directory to the "6_for_manuscript" directory

### NEED TO RE-DO PRECOMPUTED ANALYSIS AND SAVE IN DIFFERENT TABLE WITH NO COVARIATE ADJUSTMENT - WILL INCREASE N AND BE CONSISTENT WITH RESULTS TABLE IN PAPER

library(dplyr)
library(RSQLite)
library(data.table)

source("../set_paths.R")
setPaths()

# read in data
mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_precomputed.sqlite"))
res <- dbReadTable(mydb, "omics_outcomes")
dbDisconnect(mydb)

meta <- unique(res$Metadata)
ephys.vars <- c(meta[grep("alpha cells", meta)], meta[grep("beta cells", meta)])
res = as.data.table(res)

# summarize number of features with FDR < 0.05 for each omics-outcome combination
all.combinations <- distinct(res[, .(Omics_type, Metadata)])
all.combinations$N <- 0
all.combinations$N[(all.combinations$Omics_type == "Single-cell gene expression") & (all.combinations$Metadata %in% ephys.vars)] <- NA
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
deg.summary <- dcast(res, Metadata ~ Omics_type) %>% as.data.frame()
deg.summary <- deg.summary[match(meta, deg.summary$Metadata), ]

write.csv(deg.summary, "./si_tables/SI_Table_2.csv", row.names = F)
