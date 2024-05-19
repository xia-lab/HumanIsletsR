# R Functions for HumanIslets web tool
# Author: Jessica Ewald

################################################################################

# visually summarize the sex, diagnosis, age, bmi, and HbA1c of donors in the downloaded data.

datasetSummary_fun <- function(tables){

  require(RSQLite)
  require(stringr)
  require(ggplot2)
  require(ggpubr)
  require(RColorBrewer)
  
  # set sqlite file path
  outcomes.path <- paste0(sqlite.path, "HI_tables.sqlite")
  omics.path <- paste0(sqlite.path, "HI_omics.sqlite")
  
  donors <- readRDS("donors.rds")
  tables <- str_split(tables, ",")[[1]]
  
  outcome.tables <- tables[tables %in% c('donor', 'isolation', 'gsis', 'peri_gluc', 'peri_leu', 'peri_olp', 'ephys_cell', 'seahorse')]
  omics.tables <- tables[tables %in% c('proc_nanostring_merge', 'proc_prot', 'proc_rnaseq', 'proc_pbrna_Alpha', 'proc_pbrna_Beta')]
  
  # get donors from outcome tables
  con <- dbConnect(RSQLite::SQLite(), outcomes.path)
  dat <- dbReadTable(con, "donor")
  for(i in c(1:length(outcome.tables))){
    temp.nm <- outcome.tables[i]
    temp <- dbReadTable(con, temp.nm)
    temp <- temp[temp$record_id %in% donors, ]
  }
  dbDisconnect(con)
  
  # get donors from omics tables
  if(length(omics.tables) > 0){
    con <- dbConnect(RSQLite::SQLite(), omics.path)
    for(i in c(1:length(omics.tables))){
      temp.nm <- omics.tables[i]
      temp <- dbReadTable(con, temp.nm)
      donors <- intersect(donors, colnames(temp))
    }
    dbDisconnect(con)
  }
  
  # get donors in filtered list
  dat$selected <- dat$record_id %in% donors
  
  # sex in dataset
  df <- data.frame(
    Sex = c("Male","Female"),
    Percent = c(sum(dat[dat$selected,c("donorsex")] == "Male"), sum(dat[dat$selected,c("donorsex")] == "Female"))
  )
  df$Percent <- round((df$Percent/length(donors))*100)
  df$Donor <- "In dataset"
  
  # sex all
  df2 <- data.frame(
    Sex = c("Male","Female"),
    Percent = c(sum(dat$donorsex == "Male", na.rm = TRUE), sum(dat$donorsex == "Female", na.rm = TRUE))
  )
  df2$Percent <- round((df2$Percent/sum(df2$Percent))*100)
  df2$Donor <- "All"
  
  df <- rbind(df, df2)
  
  
  bp.sex <- ggplot(df, aes(x=Donor, y=Percent, fill=Sex)) +
    geom_bar(width=0.9, stat="identity") +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    ggtitle("Sex") + 
    theme(legend.position = c(0.5, 0.8), legend.title = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(size=0.25, linetype = "solid", colour="black")) 
  
  # diabetes type in dataset
  df <- data.frame(
    Type = c("Type 1", "Type 2", "None"),
    Percent = c(sum(dat[dat$selected,c("T1_diabetes")] == "1"), sum(dat[dat$selected,c("T2_diabetes")] == "1"), 
                sum(dat[dat$selected,c("T1_diabetes")] == "0" & dat[dat$selected,c("T2_diabetes")] == "0"))
  )
  df$Percent <- round((df$Percent/length(donors))*100)
  df$Donor <- "In dataset"
  
  # diabetes type all
  df2 <- data.frame(
    Type = c("Type 1", "Type 2", "None"),
    Percent = c(sum(dat$T1_diabetes == "1", na.rm = TRUE), sum(dat$T2_diabetes == "1", na.rm = TRUE), 
                sum(dat$T1_diabetes == "0" & dat$T2_diabetes == "0", na.rm = TRUE))
  )
  df2$Percent <- round((df2$Percent/sum(df2$Percent))*100)
  df2$Donor <- "All"
  df <- rbind(df, df2)
  
  bp.type <- ggplot(df, aes(x=Donor, y=Percent, fill=Type)) +
    geom_bar(width=0.9, stat="identity") +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    ggtitle("Diabetes Type") +
    theme(legend.position = c(0.5, 0.8), legend.title = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(size=0.25, linetype = "solid", colour="black"))
  
  # hist
  df <- data.frame(
    Age = dat$donorage,
    BMI = dat$bodymassindex,
    HbA1c = dat$hba1c,
    Donor = "All"
  )
  df <- rbind(df, data.frame(
    Age = dat$donorage[dat$selected],
    BMI = dat$bodymassindex[dat$selected],
    HbA1c = dat$hba1c[dat$selected],
    Donor = "In dataset"
  ))
  
  hist.age <- ggplot(df, aes(x = Age, fill = Donor)) +
    geom_histogram(alpha=0.5, position = "identity") +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    theme(legend.position = c(0.2, 0.8), legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(size=0.25, linetype = "solid", colour="black")) +
    ylab("Count")
  
  hist.bmi <- ggplot(df, aes(x = BMI, fill = Donor)) +
    geom_histogram(alpha=0.5, position = "identity") +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    theme(legend.position = c(0.8, 0.8), legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(size=0.25, linetype = "solid", colour="black")) +
    ylab("Count")
  
  hist.hba1c <- ggplot(df, aes(x = HbA1c, fill = Donor)) +
    geom_histogram(alpha=0.5, position = "identity") +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    theme(legend.position = c(0.8, 0.8), legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(size=0.25, linetype = "solid", colour="black")) +
    ylab("Count")
  
  # arrange plots
  col1 <- ggarrange(bp.sex, bp.type, ncol = 1)
  col2 <- ggarrange(hist.age, hist.bmi, hist.hba1c, ncol = 1)
  fig <- ggarrange(col1, col2, nrow = 1, widths = c(1,2))
  
  # plot data 
  Cairo::Cairo(file = "data_summary.png", unit="in", res=300, width=5.5, 
               height= 9.5, type="PNG", bg="white");
  print(fig)
  dev.off()
  
 ggsave("data_summary.svg", plot = fig, width = 5.5, height = 9.5, dpi = 300)
  ggsave("data_summary.pdf", plot = fig, width = 5.5, height = 9.5, dpi = 300)

  # string to return
  donor.num <- length(donors)
  
  per.missing <- NA
  
  res <- paste(c(donor.num, per.missing, "RES-OK"), collapse = ";")
  return(res)

}


################################################################################

metaplot_fun = function(meta, x.label, imgNm) {
  
  require(RColorBrewer)
  require(RSQLite)
  require(dplyr)
  require(ggplot2)
  require(Cairo)
  
  # set file paths
  sqlite.path <- paste0(sqlite.path, "HI_tables.sqlite");
  
  donor.meta <- c("donorage", "donorsex", "donationtype", "bodymassindex", "hba1c", "hla_a2", 
                  "diagnosis", "other_condition", "percentieqrecoverypostculture", "pdisletparticleindex",
                  "pdinsulinperieq", "pdinsulindnaratio", "predistributionculturetime")
  
  if (meta %in% donor.meta) {
    table.name <- "donor"
  } else {
    table.name <- "isolation"
  }
  
  query <- paste0('SELECT ', meta, ' FROM ', table.name)
  
  # get data
  mydb <- dbConnect(RSQLite::SQLite(), sqlite.path)
  dat <- dbGetQuery(mydb, query)
  dbDisconnect(mydb)
  if(colnames(dat) %in% c("cryotubesremaining","sftubesremaining")){
      dat =data.frame(meta=as.numeric(dat[!is.na(dat[,1]),]),stringsAsFactors = F) 
   }else{
     colnames(dat) <- "meta"
   }
  x.label <- gsub("\\.", " ", x.label)
  
  # create ggplot plot object
  if( class(dat$meta) == "numeric" ) {
    a <- ggplot(dat, aes(x = meta)) +
      geom_histogram(fill="black") +
      theme_bw() +
      xlab(x.label) +
      ylab("Count") +
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    a <- ggplot(dat, aes(x = meta)) +
      geom_bar(fill="black") +
      theme_bw() +
      xlab(x.label) +
      ylab("Count") +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  # plot data 
  imgNm <- paste0(imgNm, ".png");
  Cairo::Cairo(file = imgNm, unit="in", res=300, width=6, 
               height= 5, type="PNG", bg="white");
  print(a)
  dev.off()
  
  # add status code
  res <- paste0(imgNm, "RES-OK");
  return(res);

}

################################################################################

plotOmicsFeature <- function(gene.id, meta.var, omics.type, donors = "all", cell, glucose, mode="tool"){

  library(RSQLite)
  library(dplyr)
  library(ggplot2)
  library(see)
  library(Cairo)
  library(rhdf5)

  # set omics category
  if(omics.type == "proc_scrna"){
    table.path <- paste0(h5.path, "sc_", cell, "_", glucose, ".h5")
    meta.table <- "ephys_cell"
  } else if(omics.type == "proc_pbrna"){
    omics.category <- "bulk"
    omics.type <- paste0(omics.type, "_", cell)
    meta.table <- "proc_metadata"
  } else {
    omics.category <- "bulk"
    if(grepl("_donor", meta.var)){
        meta.table <- "ephys_donor"
    } else {
        meta.table <- "proc_metadata"
    }
  }

  if(mode == "tool"){
    # get data
    if(omics.type == "proc_scrna"){
      genes <- h5read(table.path, "meta/genes/entrez")
      gene.ind <- which(genes == gene.id)

      feat.info <- h5read(table.path, "meta/genes") %>% as.data.frame()
      feat.info <- feat.info[gene.ind, ]
      cells <- h5read(table.path, "meta/cells/cellid")

      feat.dat <- h5read(table.path, "data/norm_expression", index = list(gene.ind, NULL)) %>% t() %>% as.data.frame()
      rownames(feat.dat) <- cells

    } else {
      mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_omics.sqlite"))
        query <- paste0("SELECT * FROM ", omics.type, " WHERE gene_id=", gene.id)
        feat.dat <- dbGetQuery(mydb, query)
      dbDisconnect(mydb)

      feat.info <- feat.dat[,c(1:3)]
      feat.dat <- t(feat.dat[,-c(1:3)])
    }

    # get metadata
    mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
    query <- paste0("SELECT * FROM ", meta.table)
    metadata <- dbGetQuery(mydb, query)

    # filter by cell and glucose if relevant
    if(grepl("_donor", meta.var)){ 
        metadata <- metadata[metadata$cell_type == cell, ]
    } else if(omics.type == "proc_scrna"){
        metadata <- metadata[metadata$cell_type == cell & metadata$glucose_mM == glucose, ]
    }

    query <- paste0("SELECT * FROM proc_variable_summary WHERE column='", meta.var, "'")
    meta.info <- dbGetQuery(mydb, query)
    dbDisconnect(mydb)
     if(!file.exists("savedAnalysis")){
    dir.create("savedAnalysis")
    cat("## R history", file = "savedAnalysis/Rhistory.R", sep = "\n")
  }

    saveRDS(feat.info,"savedAnalysis/feat.info.rds")
    saveRDS(feat.dat,"savedAnalysis/feat.dat.rds")
    saveRDS(meta.info,"savedAnalysis/meta.info.rds")
    rcmd <- gsub("tool","local",rcmd)
    write(rcmd, file = "savedAnalysis/Rhistory.R", append = TRUE);

    if(!file.exists("savedAnalysis/plotOmicsFeature.R")){
       dump("plotOmicsFeature", file = "savedAnalysis/plotOmicsFeature.R",append=T)
    }  
     
  } else {
      feat.info <- readRDS("feat.info.rds")
      feat.dat <- readRDS("feat.dat.rds")
      meta.info <- readRDS("meta.info.rds")
  }

  meta.type <- meta.info$type
  
  if(omics.type != "proc_scrna"){
    df <- merge(metadata[,c("record_id", meta.var)], feat.dat, by.x = "record_id", by.y = "row.names", all = FALSE)
    df <- na.omit(df)
    colnames(df) <- c("record_id", "meta", "feature")

    # filter by donors
    if(donors == "subset"){
      donor.list <- readRDS("donors.rds")
      df <- df[df$record_id %in% donor.list,]
    }
      
    if(meta.type == "disc"){
       p <- ggplot2::ggplot(df, aes(x = meta, y = feature, fill = meta)) +
       geom_violin(trim = FALSE, aes(color = meta), show.legend = FALSE) + 
       geom_jitter(height = 0, width = 0.05, show.legend = FALSE) +
       theme(legend.position = "none") +  xlab(meta.info$axis_title) + ylab(feat.info$symbol) +
       stat_summary(fun=mean, colour="yellow", geom="point", shape=18, size=3, show.legend = FALSE) +
       scale_fill_okabeito() + 
       scale_color_okabeito() + 
       theme(axis.text.x = element_text(angle=90, hjust=1)) +
       theme_bw() +
       theme(plot.margin = margin(t=0.35, r=0.25, b=0.15, l=0.25, "cm"), axis.text = element_text(size=10))
    } else {
      p <- ggplot2::ggplot(df, aes(x=meta, y=feature))+
      geom_point(size=2) + theme_bw()  + geom_smooth(method=lm,se=T)+
      xlab(meta.info$axis_title) + ylab(feat.info$symbol) +
      theme(axis.text.x = element_text(angle=90, hjust=1)) + guides(size="none") + theme_bw() +
      theme(plot.margin = margin(t=0.35, r=0.25, b=0.15, l=0.25, "cm"), axis.text = element_text(size=10))
    }

  } else {
    df <- merge(metadata[,c("cell_id", meta.var)], feat.dat, by.x = "cell_id", by.y = "row.names", all = FALSE)
    df <- na.omit(df)
    colnames(df) <- c("cell_id", "meta", "feature")
    df$meta <- as.numeric(df$meta)
    df$feature <- as.numeric(df$feature)

    # filter by donors
    if(donors == "subset"){
      donor.list <- readRDS("donors.rds")
      cells.keep <- metadata$cell_id[metadata$record_id %in% donor.list]
      df <- df[df$cell_id %in% cells.keep,]
    }
    
    p <- ggplot2::ggplot(df, aes(x=meta, y=feature))+
        geom_point(size=2, alpha = 0.25) + theme_bw() + 
        xlab(meta.info$display) + ylab(feat.info$symbol) +
        scale_x_continuous(trans='log10')+
        theme(axis.text.x = element_text(angle=90, hjust=1)) + guides(size="none") + theme_bw() +
        theme(plot.margin = margin(t=0.35, r=0.25, b=0.15, l=0.25, "cm"), axis.text = element_text(size=10))
  }
  
  write.csv(df,"df.csv",row.names=F)
  # print out plot
  imgNm <- paste0(gene.id, ".png")
  Cairo(file = paste0(gene.id, ".png"), width=4, height=4, type="png", bg="white", unit="in", dpi=150)
    print(p)
  dev.off()
   ggplot2::ggsave(paste0(gene.id, ".svg"), plot = p, width = 4, height = 4, dpi = 300)  
 ggplot2::ggsave(paste0(gene.id, ".pdf"), plot = p, width = 4, height = 4, dpi = 300)  

  return(paste0("RES-OK",imgNm))
}


################################################################################

plotExpressionByCell <- function(gene.id, display = FALSE){
  library(rhdf5)
  library(ggplot2)
  library(ggbeeswarm)
  library(dplyr)
  library(data.table)
  library(Cairo)

  # set hdf5 file path
  hdf5.path <- paste0(h5.path, "ma_singlecell.h5");

  genes <- h5read(hdf5.path, "meta/genes/entrez")

  gene.ind <- which(genes == gene.id)
  if(length(gene.ind) == 0){return("RES-NO")}

  symbol <- h5read(hdf5.path, "meta/genes/symbol")[gene.ind]
  
  df <- data.frame(celltypes = h5read(hdf5.path, "meta/cells/celltype"),
                   value = c(h5read(hdf5.path, "data/norm_expression", index = list(gene.ind, 1:143181))))
  H5close()
  
  dt <- as.data.table(df)
  summary <- dt[,.(mean = mean(value), num_obs = .N, num_zero = sum(value==0),
                   pro_zero = sum(value==0)/.N, pro_zero = sum(value==0)/.N), by = celltypes]
  
  df <- df[df$value > 0, ]
  
  df <- merge(df, summary, by = "celltypes")
  df$pro_zero <- round(df$pro_zero*100,2)
  
  # plotting takes a long time if so many cells
  max.num <- 1500
  cell.freq <- table(df$celltypes) %>% as.data.frame()
  abun.cells <- cell.freq$Var1[cell.freq$Freq > max.num] %>% as.character()
  num.cells <- cell.freq$Freq[cell.freq$Freq > max.num] %>% as.numeric()
  if(length(abun.cells) > 0){
    for(i in c(1:length(abun.cells))){
      inds <- which(df$celltypes == abun.cells[i])
      rem.inds <- sample(inds, num.cells[i]-max.num, replace = FALSE)
      df <- df[-rem.inds, ]
    }
  }
  
  df$celltypes <- factor(df$celltypes, 
                         levels = c("beta", "alpha", "delta", "PP", "epsilon",
                                    "acinar", "ductal", "endothelial",
                                    "Schwann",
                                    "activated stellate", "quiescent stellate", "mast", "macrophage", "Cytotoxic T"))
  p <- ggplot(df, aes(x = celltypes, y = value, fill = pro_zero, color = pro_zero)) +
    geom_quasirandom(varwidth = TRUE, alpha = 0.5) +
    theme_bw() +
    coord_flip() +
    ggtitle(paste0("Non-zero expression of ", symbol, " in non-diabetic donors"))+
    ylab("Normalized expression values") +
    theme(legend.justification = c(1, 1), legend.position = c(1, 1), plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(), 
          axis.text.y=element_text(size=12, color = "black"),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA)) +
    labs(fill = "% Zero Counts") +
    guides(color = "none") +
    scale_x_discrete(drop=FALSE)
 
  if(display){
    return(p)
  } else {
    # print out plot
    imgNm <- paste0(gene.id, "_singlecell.png")
    Cairo(file = imgNm, width=7, height=9, type="png", bg="white", unit="in", dpi=150)
    print(p)
    dev.off()
    ggplot2::ggsave(paste0(gene.id, "_singlecell.svg"), plot = p, width = 7, height = 9, dpi = 300)  
    ggplot2::ggsave(paste0(gene.id, "_singlecell.pdf"), plot = p, width = 7, height = 9, dpi = 300)  

    return(paste0("RES-OK",imgNm))
  }
}


################################################################################

plotPathwayHeatmap <- function(pathName, funcLib, analysisVar, omicsType, varGroup, cell, glucose, donors = "all",mode="local"){
   
    library(dplyr)
    library(pheatmap)
    library(RSQLite)
    library(rhdf5)
    library(smoother)
    library(RColorBrewer)
    
    # set omics category
    if(omicsType == "proc_scrna"){
      table.path <- paste0(h5.path, "sc_", cell, "_", glucose, ".h5")
      meta.table <- "ephys_cell"
    } else if(omicsType == "proc_pbrna"){
      omics.category <- "bulk"
      omics.type <- paste0(omicsType, "_", cell)
      meta.table <- "proc_metadata"
    } else {
      omics.category <- "bulk"
      if(grepl("_donor", analysisVar)){
        meta.table <- "ephys_donor"
      } else {
        meta.table <- "proc_metadata"
      }
    }
    if(mode=="tool"){
     # set library file path
    lib.path <- paste0(other.tables.path, "libraries/");
    
    # get omics data
    if(omicsType == "proc_scrna"){

      feature_info <- h5read(table.path, "meta/genes") %>% as.data.frame()
      H5close()
      
      feature_info <- feature_info[,c(1,3,2)]
      colnames(feature_info) <- c("gene_id", "symbol", "name")
      
      cells <- h5read(table.path, "meta/cells/cellid")
      H5close()
      
      feature_table <- h5read(table.path, "data/norm_expression") %>% as.data.frame()
      H5close()
      
      colnames(feature_table) <- cells
      
      genes.keep <- !is.na(feature_info$gene_id)
      feature_table <- feature_table[genes.keep, ]
      feature_info <- feature_info[genes.keep, ]
      rownames(feature_table) <- feature_info$gene_id
      feature_table[feature_table < 0.0001] <- NA
    } else {
      mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_omics.sqlite"))
      if(omicsType == "proc_pbrna"){
        table.nm <- paste0(omicsType, "_", cell)
        feature_table <- dbReadTable(mydb, table.nm)
      } else {
        feature_table <- dbReadTable(mydb, omicsType)
      }
      dbDisconnect(mydb)
      rownames(feature_table) <- feature_table$gene_id
      feature_info <- feature_table[,c(1:3)]
      feature_table <- feature_table[,-c(1:3)]
    }
    
    # get metadata
    mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
    query <- paste0("SELECT * FROM ", meta.table)
    metadata <- dbGetQuery(mydb, query)
    
    # filter by cell and glucose if relevant
    if(grepl("_donor", analysisVar)){ 
      metadata <- metadata[metadata$cell_type == cell, ]
    } else if(omicsType == "proc_scrna"){
      metadata <- metadata[metadata$cell_type == cell & metadata$glucose_mM == glucose, ]
    }
    
    query <- paste0("SELECT * FROM proc_variable_summary WHERE column='", analysisVar, "'")
    meta.info <- dbGetQuery(mydb, query)
    dbDisconnect(mydb)
 
    # get pathway genes
    libraryRDS <- readRDS(paste0(lib.path, funcLib, ".rds"))
    pathGenes <- libraryRDS$sets[[which(libraryRDS$term == pathName)]]

    saveRDS(feature_info,"savedAnalysis/feature_info.rds")
    saveRDS(feature_table,"savedAnalysis/feature_table.rds")
    saveRDS(meta.info,"savedAnalysis/meta.info.rds")
    saveRDS(pathGenes,"savedAnalysis/pathGenes.rds")
    rcmd <- gsub("tool","local",rcmd)
    write(rcmd, file = "savedAnalysis/Rhistory.R", append = TRUE);
     if(!file.exists("savedAnalysis/plotPathwayHeatmap.R")){
       dump("plotPathwayHeatmap", file = "savedAnalysis/plotPathwayHeatmap.R",append=T)
     }     
    
   }else{
       feature_info <- readRDS("feature_info.rds")
     feature_table <- readRDS("feature_table.rds")
     meta.info <- readRDS("meta.info.rds")
      pathGenes <- readRDS("pathGenes.rds")
   }
    meta.type <- meta.info$type

    # get only relevant metadata
    if(omicsType == "proc_scrna"){
      # filter by donors
      if(donors == "subset"){
        donor.list <- readRDS("donors.rds")
        cells.keep <- metadata$cell_id[metadata$record_id %in% donor.list]
        metadata <- metadata[metadata$cell_id %in% cells.keep,]
        if(dim(metadata)[1] < 10){return("RES-NO")}
      }        
      metadata <- metadata[,c("cell_id", analysisVar)]
      colnames(metadata)[1] <- "record_id"
    } else {
      # filter by donors
      if(donors == "subset"){
        donor.list <- readRDS("donors.rds")
        metadata <- metadata[metadata$record_id %in% donor.list,]
        if(dim(metadata)[1] < 10){return("RES-NO")}
      }

      metadata <- metadata[,c("record_id", analysisVar)]
    }

    if(meta.type == "cont"){metadata[,analysisVar] <- as.numeric(metadata[,analysisVar])}
    metadata <- metadata[metadata$record_id %in% colnames(feature_table), ]
    metadata <- metadata[order(metadata[,analysisVar]), ]
    metadata <- metadata[!is.na(metadata[,analysisVar]), ]

    # make feature table match metadata table
    feature_table <- feature_table[rownames(feature_table) %in% pathGenes, colnames(feature_table) %in% metadata$record_id]
    feature_table <- feature_table[,match(metadata$record_id, colnames(feature_table))]

    # convert Entrez to gene symbol
    feat.vec <- rownames(feature_table)
    hit.inx <- match(feat.vec, feature_info$gene_id)
    feature_info <- feature_info[hit.inx, ]

    # convert Entrez to gene symbol
    feat.vec <- rownames(feature_table)
    hit.inx <- match(feat.vec, feature_info$gene_id)
    feature_info <- feature_info[hit.inx, ]
    feature_info$symbol[is.na(feature_info$symbol)] <- feature_info$gene_id[is.na(feature_info$symbol)]
    rownames(feature_table) <- feature_info$symbol

    # highlight DEGs
    dea.res <- read.csv("dea_results.csv")
    deg.ind <- which(rownames(feature_table) %in% dea.res$Feature[dea.res$sig != "NS"])
    if(length(deg.ind) > 0){
      rownames(feature_table)[deg.ind] <- paste0(rownames(feature_table)[deg.ind], "***")
    }

    # remove features with too many NAs
    num.vals <- apply(feature_table, 1, function(x){sum(!is.na(x))})
    feature_table <- feature_table[num.vals > 7, ]
    
    # create heatmap column annotation
    colAnn <- data.frame(Metadata = metadata[,2])
    rownames(colAnn) <- metadata$record_id

    # smooth values
    smooth.df <- feature_table
    smooth.df[is.na(feature_table)] <- min(smooth.df, na.rm = T) # must replace NA with low value for smoothing
    smooth.df <- apply(smooth.df, 1, function(x){
      smth.gaussian(as.numeric(x), window = 0.03, alpha = 2.5, tails = TRUE, na.rm = TRUE)
    }) %>% t() %>% as.data.frame()
    colnames(smooth.df) <- colnames(feature_table)
    smooth.df[is.na(feature_table)] <- NA # add NAs back
    
    # create heatmap
    if(dim(feature_table)[2] < 50){ show.donorID = TRUE } else { show.donorID = FALSE }

    hm <- pheatmap(smooth.df, show_colnames = show.donorID, annotation_col = colAnn, cluster_cols = FALSE,
                   legend = TRUE, annotation_names_col = FALSE, scale = "row", silent = TRUE, border_color = NA,
                   main = pathName, na_col = "black", color=colorRampPalette(c("#000080", "#FFF300"))(30))

    # plot heatmap
    hm.height <- dim(feature_table)[1]*0.2 + 0.75
    Cairo::Cairo(file = "pathway_heatmap.png", unit="in", res=300, width=10, 
                 height= hm.height, type="PNG", bg="white");
        print(hm)
    dev.off()

      ggplot2::ggsave("pathway_heatmap.svg", plot = hm, width = 10, height =  hm.height, dpi = 300)
    ggplot2::ggsave("pathway_heatmap.pdf", plot = hm, width = 10, height =  hm.height, dpi = 300)
    

    return("RES-OK")
}


################################################################################
plotDataAvailHeatmap <- function(){
  
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
  
  avail_df <- t(avail_df)
  
  disp.rownames <- c("Clinical metadata", "Technical metadata", "Static insulin secretion", "Dynamic insulin secretion",
                     "Oxygen consumption", "Electrophysiology", "Proteomics (bulk)", "Nanostring (bulk)", "RNAseq (bulk)",
                     "RNAseq (pseudobulk)", "RNAseq (patchSeq)")
  
  pheatmap(avail_df, legend = FALSE, show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE,
           color = c("white", "#43A047"), angle_col=90, labels_row = disp.rownames, border_color = NA,
           filename=paste0(other.tables.path, "/display_data/avail_heatmap.png"), width=10, height=5)
}
