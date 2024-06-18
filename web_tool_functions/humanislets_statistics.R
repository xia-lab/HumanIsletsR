# R Functions for HumanIslets web tool
# Author: Jessica Ewald

################################################################################


performGSEA <- function(funcLib = "kegg", rank.stat = "coef", fdr = 0.05, collapse = "true",mode="local"){

  library(fgsea)
  library(data.table)
  
  set.seed(42)
  
  fdr <- as.numeric(fdr)
  
  if(mode=="tool"){
    # set library filepath
    lib.path <- paste0(other.tables.path, "libraries/");    
    
    # process libraries
    libraryRDS <- readRDS(paste0(lib.path, funcLib, ".rds"))
    libraryList <- libraryRDS$sets
    saveRDS(libraryRDS,"savedAnalysis/libraryRDS.rds")
    saveRDS(libraryList,"savedAnalysis/libraryList.rds")
    rcmd <- gsub("tool","local", rcmd)
    write(rcmd, file = "savedAnalysis/Rhistory.R", append = TRUE);
   if(!file.exists("savedAnalysis/performGSEA.R")){
     dump("performGSEA", file = "savedAnalysis/performGSEA.R",append=T)
   }
       
    }else{
       libraryRDS <- readRDS("libraryRDS.rds")
       libraryList <- readRDS("libraryList.rds")
    }
  
  # read in DEGs
  degs <- read.csv("dea_results.csv")
  
  if(rank.stat == "coef"){
    ranks <- as.numeric(degs[,2])
  } else {
    ranks <- as.numeric(degs$T_statistic)
  }
  names(ranks) <- degs$Gene_ID
  set.names <- data.frame(names = libraryRDS$term, IDs = names(libraryRDS$sets))
  
  # perform analysis
  res <- fgsea(pathways = libraryList,
               stats = ranks,
               minSize = 15,
               maxSize = 500,
               )
  
  # find independent pathways if collapse set to true
  if(collapse == "true"){
    ind.pathways <- collapsePathways(fgseaRes = res,
                                     pathways = libraryList,
                                     stats = ranks)
    res <- res[res$pathway %in% ind.pathways$mainPathways]
  }

  res <- as.data.frame(res)

  # format results for download
  if(dim(res)[1] > 0){
    res <- res[,-c(4,5,8)] # remove unneeded columns
    res <- merge(res, set.names, by.x = "pathway", by.y = "IDs")
    res <- res[order(res$pval),]
    res <- res[,c(6,2:5,1)]
    res[,2:3] <- signif(res[,2:3], digits = 3)
    res[,4] <- round(res[,4], digits = 3)
    colnames(res) <- c("Set Name", "P_value", "Adj P_value", "Normalized ES", "Set Size", "Set ID")
    write.csv(res, "functional_results.csv", row.names = FALSE)
    num.sig <- sum(res$`Adj P_value` < fdr, na.rm = TRUE)

    # format table for ridgeline plot
    genesPW <- reshape::melt(libraryList, level = 1)
    colnames(genesPW) <- c("GeneID", "Set ID")
    res$pathNegLogP <- -log10(res$P_value)
    res$RankSig <- rank(-res$pathNegLogP)
    genesPW <- merge(genesPW, res, by = "Set ID", all.y = TRUE, all.x = FALSE)
    genesPW$Sig <- as.character(genesPW$`Adj P_value` < fdr)
    genesPW <- genesPW[,c("Set ID", "GeneID", "Set Name", "Set Size", "Adj P_value", "Sig", "pathNegLogP", "RankSig")]
    genesPW <- merge(genesPW, data.frame(GeneID = names(ranks), GeneStat = unname(ranks)), by = "GeneID", all.x = TRUE, all.y = FALSE)
    genesPW <- na.omit(genesPW)
    genesPW <- as.data.table(genesPW)
    genesPW <- genesPW[,.(GeneID, `Set Name`, `Set Size`, pathFDR = `Adj P_value`, Sig, pathNegLogP, GeneStat, MeanGeneStat = mean(GeneStat), RankSig), by = .(`Set ID`)]
    genesPW <- genesPW[order(MeanGeneStat)]
    write.csv(genesPW, "ridgeline.csv", row.names = FALSE)

    return(paste0("RES-OK;", num.sig))
  } else {
    return("RES-NO")
  }
}

################################################################################


performORA <- function(funcLib = "kegg", fdr = 0.05, collapse = "true",mode="local"){

  library(fgsea)
  library(data.table)

  fdr <- as.numeric(fdr)
  
  # set library filepath
 
  # read in DEGs
  degs <- read.csv("dea_results.csv")
  universe <- as.character(degs$Gene_ID)
  hits <- as.character(degs$Gene_ID[degs$sig != "NS"])
  
  if(length(hits) > 10){
    # process libraries
    if(mode=="tool"){
        lib.path <- paste0(other.tables.path, "libraries/")
       libraryRDS <- readRDS(paste0(lib.path, funcLib, ".rds"))
       libraryList <- libraryRDS$sets
      
       saveRDS(libraryRDS,"savedAnalysis/libraryRDS.rds")
       saveRDS(libraryList,"savedAnalysis/libraryList.rds")
        rcmd <- gsub("tool","local",rcmd)
      write(rcmd, file = "savedAnalysis/Rhistory.R", append = TRUE);
     if(!file.exists("savedAnalysis/performORA.R")){
       dump("performORA", file = "savedAnalysis/performORA.R",append=T)
     }
    }else{
       libraryRDS <- readRDS("libraryRDS.rds")
       libraryList <- readRDS("libraryList.rds")
    }


    set.names <- data.frame(names = libraryRDS$term, IDs = names(libraryRDS$sets))
    
    # perform analysis
    res <- fora(pathways = libraryList, 
                genes = hits,
                universe = universe,
                minSize = 15, maxSize = 500)

    # find independent pathways if collapse set to true
    if(collapse == "true"){
      ind.pathways <- collapsePathwaysORA(foraRes = res,
                                       pathways = libraryList,
                                       genes = hits,
                                       universe = universe)
      res <- res[res$pathway %in% ind.pathways$mainPathways]
    }
    res <- as.data.frame(res)
    
    # format results
    res <- merge(res, set.names, by.x = "pathway", by.y = "IDs")
    res <- res[order(res$pval),]
    res <- res[,-6] # remove list of genes
    res <- res[,c(6,2:5,1)]
    res[,2:3] <- signif(res[,2:3], digits = 3)
    colnames(res) <- c("Set Name", "P_value", "Adj P_value", "Hits", "Set Size", "Set ID")
    write.csv(res, "functional_results.csv", row.names = FALSE)
    num.sig <- sum(res$`Adj P_value` < fdr, na.rm = TRUE)

    # format results for ridgeline
    genesPW <- reshape::melt(libraryList, level = 1)
    colnames(genesPW) <- c("GeneID", "Set ID")
    res$pathNegLogP <- -log10(res$P_value)
    res$RankSig <- rank(res$P_value)
    genesPW <- merge(genesPW, res, by = "Set ID", all.y = TRUE, all.x = FALSE)
    genesPW$Sig <- as.character(genesPW$`Adj P_value` < fdr)
    genesPW <- genesPW[,c("Set ID", "GeneID", "Set Name", "Set Size", "Adj P_value", "Sig", "pathNegLogP", "RankSig")]
    genesPW <- merge(genesPW, data.frame(GeneID = degs$Gene_ID, GeneStat = degs[,2]), by = "GeneID", all.x = TRUE, all.y = FALSE)
    genesPW <- na.omit(genesPW)
    genesPW <- as.data.table(genesPW)
    genesPW <- genesPW[,.(GeneID, `Set Name`, `Set Size`, pathFDR = `Adj P_value`, Sig, pathNegLogP, GeneStat, MeanGeneStat = mean(GeneStat), RankSig), by = .(`Set ID`)]
    genesPW <- genesPW[order(MeanGeneStat)]
    write.csv(genesPW, "ridgeline.csv", row.names = FALSE)

    return(paste0("RES-OK;", num.sig))
  } else {
    return("RES-NO; fewer than 10 DEGs")
  }
  
}


################################################################################


## perform comparison analysis based on  multi-variate linear regression 
DonorRegression <- function(
    varGroup, # what table it comes from
    analysisVar, # metadata variable name
    ref = NULL, # reference class from analysis.var metadata (only if categorical)
    contrast = NULL,  # comparison class from analysis.var (only if categorical)
    fixedEffects = NULL,  # metadata variables to adjust for
    omicsType = NULL,
    pvalThresh = "0.05",
    donors = "all",
    cell = "alpha",
    glucose = "1",
    fdr = "true",
   mode = "local"){ 
 
  #mode = "tool";
  if(length(contrast) == 0){contrast = 'NULL'}
  
  # load libraries
  library(limma)
  library(dplyr)
  library(RSQLite)

  # process inputs
  pvalThresh <- as.numeric(pvalThresh)
  if(fixedEffects == "NA"){
    fixedEffects <- NULL
  } else {
    fixedEffects <- strsplit(fixedEffects, ",")[[1]]
  }
  
if(mode == "tool"){
  # get data
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
  
  
  # get metadata
  mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
  if(varGroup == 'ephys_donor'){
    metadata <- dbReadTable(mydb, "ephys_donor")
    metadata <- metadata[metadata$cell_type == cell & metadata$glucose_mM == glucose, ]
    
    if(length(fixedEffects) > 0){
        other.metadata <- dbReadTable(mydb, "proc_metadata")
        metadata <- merge(metadata, other.metadata, by = "record_id")
    }
  } else {
    metadata <- dbReadTable(mydb, "proc_metadata")
  }
  meta.info <- dbReadTable(mydb, "proc_variable_summary")
  dbDisconnect(mydb)

  # get meta info and modify analysisVar
  analysis.type <- meta.info$type[meta.info$column == analysisVar]
  nonendo <- readRDS(paste0(other.tables.path, "analysis_input/cell_signal.rds"))
  # make dir: savedAnalysis
  if(!file.exists("savedAnalysis")){
    dir.create("savedAnalysis")
    cat("## R history", file = "savedAnalysis/Rhistory.R", sep = "\n")
  }
   saveRDS(feature_table, "savedAnalysis/features.rds")
   saveRDS(metadata, "savedAnalysis/meta.rds")
   saveRDS(meta.info, "savedAnalysis/meta_info.rds")
   saveRDS(feature_info, "savedAnalysis/feature_info.rds")
   saveRDS(nonendo, "savedAnalysis/nonendo.rds")
   
   if( exists("rcmd") ){
     rcmd <- gsub("tool","local",rcmd)
     write(rcmd, file = "savedAnalysis/Rhistory.R", append = TRUE);
     if(!file.exists("savedAnalysis/DonorRegression.R")){
       dump( "DonorRegression", file = "savedAnalysis/DonorRegression.R",append=T)
       
     }
   }
 
} else { # mode = "local"
    feature_table <- readRDS("features.rds")
    metadata <- readRDS("meta.rds")
    meta.info <- readRDS("meta_info.rds")
    feature_info <- readRDS("feature_info.rds")
    nonendo <- readRDS("nonendo.rds")
    analysis.type <- meta.info$type[meta.info$column == analysisVar]
}
  all.vars <- c(analysisVar, fixedEffects)

  # filter to keep only data with complete metadata & omics data
  metadata <- metadata[metadata$record_id %in% colnames(feature_table), c("record_id", all.vars)]
  metadata <- na.omit(metadata)
  if(dim(metadata)[1] < 10){return("RES-NO")}
  rownames(metadata) <- metadata$record_id
  feature_table <- feature_table[,colnames(feature_table) %in% metadata$record_id]
  
  # make sure metadata and omics data donors are in the same order
  feature_table <- feature_table[,match(metadata$record_id, colnames(feature_table))]

  # eliminate any feature with fewer than 10 observations after complete case analysis with metadata variables
  feature.keep <- apply(feature_table, 1, function(x){sum(!is.na(x)) >= 9})
  feature_table <- feature_table[feature.keep, ]
  feature_info <- feature_info[feature.keep, ]
 
  # filter by donors
  if(donors == "subset"){
    donor.list <- readRDS("donors.rds")
    feature_table <- feature_table[,colnames(feature_table) %in% donor.list]
    metadata <- metadata[metadata$record_id %in% donor.list, ]
  }

  # check if there are enough samples
  n.samps <- dim(feature_table)[2]
  if(analysis.type == "disc" && contrast != "anova"){
    if(n.samps < 10){
      return("RES-NO")
    }
    
    if(sum(metadata[,analysisVar] == ref) < 5){
      return("RES-NO")
    }
    
    if(sum(metadata[,analysisVar] == contrast) < 5){
      return("RES-NO")
    }
    
  } else {
    if(n.samps < 10){
      return("RES-NO")
    }
  }
  
  # perform analysis
  if(analysis.type == "disc"){
    
    # make design matrix
    grp.nms <- sort(unique(metadata[,analysisVar]))
    if(length(all.vars) == 1){
      design <- model.matrix(formula(paste0("~ 0 + ", all.vars)), data = metadata)
    } else {
      design <- model.matrix(formula(paste0("~ 0 + ", all.vars[1], paste0(" + ", all.vars[2:length(all.vars)], collapse = ""))), data = metadata)
    }
    colnames(design)[1:length(grp.nms)] <- grp.nms
    
    # make contrast matrix
    myargs <- list();
    if(contrast == "anova"){ 
      contrasts <- grp.nms[grp.nms != ref];
      myargs <- as.list(paste(contrasts, "-", ref, sep = "")); 
    } else {
      myargs <- as.list(paste(contrast, "-", ref, sep = ""));
    }
    myargs[["levels"]] <- design
    contrast.matrix <- do.call(makeContrasts, myargs)
    
    # get results
    fit <- lmFit(feature_table, design, trend = TRUE, robust = TRUE)
    fit <- contrasts.fit(fit, contrast.matrix)
    fit <- eBayes(fit)
    res.table <- topTable(fit, number = Inf)
    
    # Remove coefficients for ANOVA contrasts
    if(contrast == "anova"){
      if(length(myargs) > 2){
        res.table <- res.table[,-c(2:(length(myargs)-1))]
      }
      res.table[,1] <- "N/A"
    }
    
    colnames(res.table)[1] <- "Log2FC"

  } else { 
    
    # make design matrix
    if(length(all.vars) == 1){
      design <- model.matrix(formula(paste0("~ ", all.vars)), data = metadata)
    } else {
      design <- model.matrix(formula(paste0("~ ", all.vars[1], paste0(" + ", all.vars[2:length(all.vars)], collapse = ""))), data = metadata)
    }
    
    # get results
    fit <- lmFit(feature_table, design, trend = TRUE, robust = TRUE)
    fit <- eBayes(fit)
    res.table <- topTable(fit, number = Inf, coef = analysisVar)
    colnames(res.table)[1] <- "Coefficient"
    
  }

  # Remove results rows with NAs
  res.table <- res.table[!is.na(res.table$P.Value), ]

  # Process results for output
  res.table <- merge(res.table, feature_info, by.x = "row.names", by.y = "gene_id")
  if(contrast != "anova"){
    res.table <- res.table[,-7]
  } else if (contrast == "anova" & length(myargs) == 2){
    res.table <- res.table[,-7]
  }

  # add endocrine signal

  colnames(res.table)[1] <- "Gene_ID"
  res.table <- merge(res.table, nonendo, by = "Gene_ID")
  
  
  res.table <- res.table[,c(7,2:6,9,1,8)]
  colnames(res.table)[1] <- "Feature"
  colnames(res.table)[3] <- "Average level"
  colnames(res.table)[4] <- "T_statistic"
  colnames(res.table)[5] <- "P_value"
  colnames(res.table)[6] <- "Adjusted p_value"
  colnames(res.table)[9] <- "Description"

  if(contrast != "anova"){
    res.table[,2:6] <- signif(res.table[,2:6], digits = 3)
  } else {
    res.table[,3:6] <- signif(res.table[,3:6], digits = 3)
  }
  res.table <- res.table[order(res.table$P_value), ]
  res.table$sig <- rep("NS", dim(res.table)[1])

  if(fdr == "true"){
    res.table$sig[res.table$`Adjusted p_value` < pvalThresh & res.table[,2] > 0] <- "up"
    res.table$sig[res.table$`Adjusted p_value` < pvalThresh & res.table[,2] < 0] <- "down"
    res.table$negLogPval <- -log10(res.table[,"P_value"])

    sig.num <- sum(res.table$`Adjusted p_value` < pvalThresh)
    sig.up <- sum(res.table$`Adjusted p_value` < pvalThresh & res.table[,2] > 0)
    sig.down <- sum(res.table$`Adjusted p_value` < pvalThresh & res.table[,2] < 0)
  } else {
    res.table$sig[res.table$P_value < pvalThresh & res.table[,2] > 0] <- "up"
    res.table$sig[res.table$P_value < pvalThresh & res.table[,2] < 0] <- "down"
    res.table$negLogPval <- -log10(res.table[,"P_value"])

    sig.num <- sum(res.table$P_value < pvalThresh)
    sig.up <- sum(res.table$P_value < pvalThresh & res.table[,2] > 0)
    sig.down <- sum(res.table$P_value < pvalThresh & res.table[,2] < 0)
  }

  write.csv(res.table, "dea_results.csv", row.names = FALSE)

  return(paste0("RES-OK;", sig.num, ";", sig.up, ";", sig.down, ";", n.samps))
}

################################################################################


## perform Spearman ranked correlation for patch-seq data

PatchseqSpearman <- function(
    analysisVar, # metadata variable name
    pvalThresh = "0.05",
    donors = "all",
    cell = "Alpha",
    glucose = "1",
    fdr = "true",
    mode="local"
){ 
  
  # load libraries
  library(dplyr)
  library(RSQLite)
  library(rhdf5)

  # process inputs
  pvalThresh <- as.numeric(pvalThresh)

 if(mode == "tool"){
    # get data
  table.path <- paste0(h5.path, "sc_", cell, "_", glucose, ".h5")

  feature_table <- h5read(table.path, "data/norm_expression") %>% as.data.frame()
  cells <- h5read(table.path, "meta/cells/cellid")
  genes <- h5read(table.path, "meta/genes") %>% as.data.frame()
  H5close()

  feature_info <- genes[,c(1,3,2)]
  colnames(feature_info) <- c("gene_id", "symbol", "name")
  colnames(feature_table) <- cells

  genes.keep <- !is.na(feature_info$gene_id)
  feature_table <- feature_table[genes.keep, ]
  feature_info <- feature_info[genes.keep, ]
  rownames(feature_table) <- feature_info$gene_id

  # get metadata
  mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
  metadata <- dbReadTable(mydb, "ephys_cell")
  dbDisconnect(mydb)
  
  # filter to keep only data with complete metadata & omics data
  metadata <- metadata[metadata$cell_id %in% colnames(feature_table), c("cell_id", "record_id", analysisVar)]
  metadata <- na.omit(metadata)
  rownames(metadata) <- metadata$cell_id
  feature_table <- feature_table[,colnames(feature_table) %in% metadata$cell_id]
   if(dim(metadata)[1] < 10){return("RES-NO")}
  
  # make sure metadata and omics data donors are in the same order
  feature_table <- feature_table[,match(metadata$cell_id, colnames(feature_table))]
   # filter by donors
  if(donors == "subset"){
    donor.list <- readRDS("donors.rds")
    cells.keep <- metadata$cell_id[metadata$record_id %in% donor.list]
    feature_table <- feature_table[,colnames(feature_table) %in% cells.keep]
    metadata <- metadata[metadata$cell_id %in% cells.keep, ]
  }

  # check if there are enough samples
  n.samps <- length(unique(metadata$record_id))
  n.cells <- dim(feature_table)[2]
  if(n.cells < 10){
    return("RES-NO")
   }
  
  # make dir: savedAnalysis
  if(!file.exists("savedAnalysis")){
    dir.create("savedAnalysis")
    cat("## R history", file = "savedAnalysis/Rhistory.R", sep = "\n")
  }
   nm = paste(analysisVar, pvalThresh, donors, cell, glucose,fdr,sep="_")
   saveRDS(feature_table, paste0("savedAnalysis/features_",nm,".rds"))
   saveRDS(metadata, paste0("savedAnalysis/meta_",nm,".rds"))
   if( exists("rcmd") ){
     rcmd <- gsub("tool","local",rcmd)
     write(rcmd, file = "savedAnalysis/Rhistory.R", append = TRUE);
     if(!file.exists("savedAnalysis/PatchseqSpearman.R")){
       dump( "PatchseqSpearman", file = "savedAnalysis/PatchseqSpearman.R",append=T)
       
     }
   }
 
} else { # mode = "local"
   nm = paste(analysisVar, pvalThresh, donors, cell, glucose,fdr,sep="_")
    feature_table <- readRDS(paste0("features_",nm,".rds"))
    metadata <- readRDS(paste0("meta_",nm,".rds"))
}
  
 
  # perform analysis
  # spearman correlation
  outcome <- metadata[,analysisVar] %>% as.numeric()
  res <- try(apply(feature_table, 1, function(x){
    temp <- cor.test(x, outcome, use = "complete.obs", method = "spearman")
    return(c(temp$estimate, temp$p.value))
  }), silent = TRUE)
  
  if(!inherits(res, "try-error")){
    res <- t(res) %>% as.data.frame()
    colnames(res) <- c("r", "pval")
    res$fdr <- p.adjust(res$pval, method = "fdr")
    res$avg.level <- apply(feature_table, 1, mean, na.rm = TRUE)
    res <- na.omit(res)
  } else {
    return("RES-NO")
  }
    

  # Process results for output
  res.table <- merge(res, feature_info, by.x = "row.names", by.y = "gene_id")
  res.table <- res.table[,c(6,2,5,3:4,1,7)]
  colnames(res.table)[1] <- "Feature"
  colnames(res.table)[2] <- "Coefficient"
  colnames(res.table)[3] <- "Average level"
  colnames(res.table)[4] <- "P_value"
  colnames(res.table)[5] <- "Adjusted p_value"
  colnames(res.table)[6] <- "Gene_ID"
  colnames(res.table)[7] <- "Description"
  res.table[,2:5] <- signif(res.table[,2:5], digits = 3)
  res.table <- res.table[order(res.table$P_value), ]
  res.table$sig <- rep("NS", dim(res.table)[1])

  if(fdr == "true"){
    res.table$sig[res.table$`Adjusted p_value` < pvalThresh & res.table[,2] > 0] <- "up"
    res.table$sig[res.table$`Adjusted p_value` < pvalThresh & res.table[,2] < 0] <- "down"
    res.table$negLogPval <- -log10(res.table[,"P_value"])

    sig.num <- sum(res.table$`Adjusted p_value` < pvalThresh)
    sig.up <- sum(res.table$`Adjusted p_value` < pvalThresh & res.table[,2] > 0)
    sig.down <- sum(res.table$`Adjusted p_value` < pvalThresh & res.table[,2] < 0)
  } else {
    res.table$sig[res.table$P_value < pvalThresh & res.table[,2] > 0] <- "up"
    res.table$sig[res.table$P_value < pvalThresh & res.table[,2] < 0] <- "down"
    res.table$negLogPval <- -log10(res.table[,"P_value"])

    sig.num <- sum(res.table$P_value < pvalThresh)
    sig.up <- sum(res.table$P_value < pvalThresh & res.table[,2] > 0)
    sig.down <- sum(res.table$P_value < pvalThresh & res.table[,2] < 0)
  }

  write.csv(res.table, "dea_results.csv", row.names = FALSE)
  
  return(paste0("RES-OK;", sig.num, ";", sig.up, ";", sig.down, ";", n.samps, ";", n.cells))
}