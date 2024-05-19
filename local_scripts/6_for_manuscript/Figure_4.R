# Manuscript Figure 4
# Author: Jessica Ewald

## Set your working directory to the "6_for_manuscript" directory

library(dplyr)
library(RSQLite)
library(ggplot2)

source("../set_paths.R")
setPaths()

# get data
mydb <- dbConnect(RSQLite::SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
composition <- dbReadTable(mydb, "computed")
dbDisconnect(mydb)

composition$beta <- (1-composition$exo_per)*composition$beta_end
composition$alpha <- (1-composition$exo_per)*composition$alpha_end
composition$delta <- (1-composition$exo_per)*composition$delta_end
composition$gamma <- (1-composition$exo_per)*composition$gamma_end

composition <- composition[,c(1,6:10)]

# make cell composition plot

stacked.bar <- reshape2::melt(composition)
colnames(stacked.bar) <- c("record_id", "cell_type", "proportion")
stacked.bar$proportion = stacked.bar$proportion*100
stacked.bar$cell_type <- as.character(stacked.bar$cell_type)
stacked.bar$cell_type[stacked.bar$cell_type == "exo_per"] <- "% non-endocrine"
stacked.bar$cell_type[stacked.bar$cell_type == "beta"] <- "% beta"
stacked.bar$cell_type[stacked.bar$cell_type == "alpha"] <- "% alpha"
stacked.bar$cell_type[stacked.bar$cell_type == "delta"] <- "% delta"
stacked.bar$cell_type[stacked.bar$cell_type == "gamma"] <- "% gamma"
stacked.bar$record_id <- factor(stacked.bar$record_id, levels = composition$record_id[order(composition$exo_per, composition$beta)])

pdf(file="./figures/fig4A.pdf", width=12, height=5)
ggplot(stacked.bar, aes(fill=cell_type, y=proportion, x=record_id)) + 
  geom_bar(position='stack', stat='identity', width=1) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(),
        legend.position = "top") +
  xlab("Donor ID") +
  ylab("Cell type %")
dev.off()
