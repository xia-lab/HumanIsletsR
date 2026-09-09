# R Functions for HumanIslets web tool
# Author: Jessica Ewald, Yao Lu

################################################################################


performGSEA <- function(funcLib = "kegg", rank.stat = "coef", fdr = 0.05, collapse = "true",mode="local"){

  tryCatch({
    library(fgsea)
  }, error = function(e){
    return(paste0("RES-NO; Failed to load fgsea package: ", e$message))
  })
  if(!requireNamespace("fgsea", quietly = TRUE)){
    return("RES-NO; fgsea package is not installed")
  }
  library(fgsea)
  library(data.table)

  set.seed(42)

  fdr <- as.numeric(fdr)

  # Methylation pathway analysis (GSEA): dea_results.csv is CpG-level (has
  # gene + pos_hg38, no Gene_ID). Branch to the methylation engine (methylGSA
  # methylRRA, GSEA mode) and return. The gene-level path below is untouched
  # for every other omic.
  .mhdr <- tryCatch(colnames(data.table::fread("dea_results.csv", nrows = 0)),
                    error = function(e) character(0))
  if(all(c("gene", "pos_hg38") %in% .mhdr) && !("Gene_ID" %in% .mhdr)){
    return(.runMethylationGSEA(funcLib, fdr, collapse, mode))
  }

  if(mode=="tool"){
    # set library filepath
    lib.path <- paste0(other.tables.path, "libraries/");

    # process libraries
    if(funcLib == "hsa_kegg"){
      lib <- qs::qread(paste0(lib.path, "kegg_hsa_met.qs"))
      libraryList <- lib$mset.list                          # keyed by KEGG pathway IDs (hsa00010, ...)
      name_lookup  <- setNames(names(lib$path.ids), lib$path.ids)  # hsa00010 → "Glycolysis..."
      readable_names <- name_lookup[names(libraryList)]
      readable_names[is.na(readable_names)] <- names(libraryList)[is.na(readable_names)]
      libraryRDS <- list(term = unname(readable_names), sets = libraryList)

    } else if(funcLib == "gem_metabolite"){
      lib <- qs::qread(paste0(lib.path, "hsaGem_metabolite.qs"))
      libraryList <- lib$sets
      libraryRDS <- lib

    } else if(funcLib == "gem_react"){
      # Load GEM reaction library
      lib = qs::qread(paste0(lib.path,"hsaGem_reaction.qs"))
      rxn_tab <- lib$reacTab
      rules <- lib$rules
 
      # Create library list from subsystems
      libraryList <- split(rxn_tab$humanGEM, rxn_tab$sysID)

      # Get unique pathway names for each subsystem ID
      subsystem_names <- unique(data.frame(
        sysID = rxn_tab$sysID,
        pathwayVis = rxn_tab$pathwayVis,
        stringsAsFactors = FALSE
      ))

      libraryRDS <- list()
      libraryRDS$term <- subsystem_names$pathwayVis[match(names(libraryList), subsystem_names$sysID)]
      libraryRDS$sets <- libraryList
     } else {
      lib_file <- paste0(lib.path, funcLib, ".rds")
      if(!file.exists(lib_file)){
        return(paste0("RES-NO; Library file not found: ", lib_file))
      }
      libraryRDS <- readRDS(lib_file)
      libraryList <- libraryRDS$sets
    }

     if(!file.exists("savedAnalysis")){
    dir.create("savedAnalysis")
    cat("## R history", file = "savedAnalysis/Rhistory.R", sep = "\n")
  }


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
  if(!file.exists("dea_results.csv")){
    return("RES-NO; dea_results.csv not found. Please run differential analysis first.")
  }
  degs <- read.csv("dea_results.csv")
  if(nrow(degs) == 0){
    return("RES-NO; dea_results.csv is empty")
  }
  print(paste("GSEA: Loaded", nrow(degs), "features. Columns:", paste(colnames(degs), collapse=", ")))

  # Check if this is flux data (has Feature column but no Gene_ID)
  is_flux_data <- "Feature" %in% colnames(degs) && !("Gene_ID" %in% colnames(degs))

  if(funcLib == "hsa_kegg"){
    if(!"kegg_id" %in% colnames(degs)){
      return("RES-NO; kegg_id column not found in DEA results. Run metabolite regression first.")
    }
    tmp <- degs[degs$kegg_id != "" & !is.na(degs$kegg_id), ]
    tmp <- tmp[!duplicated(tmp$kegg_id), ]
    if(rank.stat == "coef"){
      ranks <- setNames(as.numeric(tmp[, 2]), tmp$kegg_id)
    } else {
      ranks <- setNames(as.numeric(tmp$T_statistic), tmp$kegg_id)
    }
    ranks <- ranks[!is.na(ranks)]

  } else if(funcLib == "gem_metabolite"){
    if(!"gem_id" %in% colnames(degs)){
      return("RES-NO; gem_id column not found in DEA results. Run metabolite regression first.")
    }
    tmp <- degs[degs$gem_id != "" & !is.na(degs$gem_id), ]
    tmp <- tmp[!duplicated(tmp$gem_id), ]
    if(rank.stat == "coef"){
      ranks <- setNames(as.numeric(tmp[, 2]), tmp$gem_id)
    } else {
      ranks <- setNames(as.numeric(tmp$T_statistic), tmp$gem_id)
    }
    ranks <- ranks[!is.na(ranks)]

  } else if(is_flux_data){
    # Flux data - use Feature column (reaction IDs)
    degs <- degs[!duplicated(degs$Feature)&!is.na(degs$Feature),]
    if(rank.stat == "coef"){
      ranks <- as.numeric(degs[,2])
    } else {
      ranks <- as.numeric(degs$T_statistic)
    }
    names(ranks) <- degs$Feature

  } else if(funcLib == "gem_react"){
    # Gene-based data with gem_react library - need to convert genes to reactions
    # GPR rules use ENSG IDs, so we need ENSG as the gene identifier

    if(all(grepl("ENSG", head(degs$Feature)))){
      # Feature already has ENSG IDs (e.g. RNA-seq) - use directly, no conversion
      degs <- degs[!duplicated(degs$Feature)&!is.na(degs$Feature),]
      if(rank.stat == "coef"){
        gene_ranks <- setNames(as.numeric(degs[,2]), degs$Feature)
      } else {
        gene_ranks <- setNames(as.numeric(degs$T_statistic), degs$Feature)
      }
    } else {
      # Feature has symbols (proteomics, pbrna) - convert symbol to ENSG directly
      degs <- degs[!duplicated(degs$Feature)&!is.na(degs$Feature),]
      out_anno <- qs::qread(paste0(other.tables.path,"libraries/gene_anno_full.qs"))
      degs$accession <- out_anno$ensembl_gene_id[match(degs$Feature, out_anno$hgnc_symbol)]
      degs <- degs[!is.na(degs$accession),]
      if(rank.stat == "coef"){
        gene_ranks <- setNames(as.numeric(degs[,2]), degs$accession)
      } else {
        gene_ranks <- setNames(as.numeric(degs$T_statistic), degs$accession)
      }
    }

    # Map genes to reactions using GPR rules
    library(stringr)

    print(paste("Converting", length(gene_ranks), "genes to reactions for GSEA"))
    print(paste("Total rules available:", length(rules)))

    # Map all genes to reactions using GPR rules (similar to ORA approach)
    universe_map <- lapply(rules, function(r) eval_gpr(r, names(gene_ranks)))
     # Filter to reactions that map to at least one gene
    rxn_mapped <- unlist(lapply(universe_map, function(x) x[["map"]]))
    universe_map <- universe_map[rxn_mapped]


    # Initialize list to store matched genes for each reaction
    matched_genes <- list()

    # For each mapped reaction, aggregate gene ranks
    # Handle compound gene IDs (with &) in GPR rules, same as ORA
    rxn_ranks <- sapply(names(universe_map), function(rxn_id){
      gpr_result <- universe_map[[rxn_id]]
      map <- gpr_result$expr
      map <- map[map$map, ]  # filter to only mapped genes
      if(nrow(map) == 0) return(NA)

      # Resolve compound gene IDs (e.g. "ENSG001 & ENSG002") to individual genes
      for(i in 1:nrow(map)){
        if(grepl("\\&", map$id[i])){
          ind_genes <- str_split(map$id[i], "&")[[1]] %>% str_trim()
          ind_vals <- gene_ranks[ind_genes]
          kp <- which.max(abs(ind_vals))
          if(length(kp) > 0){
            map$id[i] <- names(ind_vals)[kp]
          }
        }
      }

      gene_vals <- gene_ranks[map$id]
      gene_vals <- gene_vals[!is.na(gene_vals)]

      if(length(gene_vals) > 0){
        # Use the gene with maximum absolute value to represent the reaction
        max_gene <- names(gene_vals)[which.max(abs(gene_vals))]
        matched_genes[[rxn_id]] <<- max_gene
        return(gene_vals[which.max(abs(gene_vals))])
      }
      return(NA)
    })

    # Remove reactions with no valid gene values
    rxn_ranks <- rxn_ranks[!is.na(rxn_ranks)]

    # Strip gene suffix from compound names (e.g. "MAR00001.ENSG00000163686" -> "MAR00001")
    rxn_ids <- sub("\\..*$", "", names(rxn_ranks))
    # For duplicate reaction IDs, keep the one with max absolute rank
    ranks_df <- data.frame(rxn = rxn_ids, val = unname(rxn_ranks), stringsAsFactors = FALSE)
    ranks_df <- ranks_df[order(-abs(ranks_df$val)), ]
    ranks_df <- ranks_df[!duplicated(ranks_df$rxn), ]
    ranks <- setNames(ranks_df$val, ranks_df$rxn)

    print(paste("Successfully mapped to", length(ranks), "reactions"))

    # Verify ranks has names
    if(is.null(names(ranks)) || length(names(ranks)) != length(ranks)){
      print("ERROR: ranks vector does not have proper names")
      return("RES-NO")
    }

    if(length(ranks) == 0){
      print("ERROR: No reactions could be mapped from genes")
      return("RES-NO")
    }

  } else {
    # Gene-based data with non-gem libraries
    if(!"Gene_ID" %in% colnames(degs)){
      return(paste0("RES-NO; Gene_ID column not found in DEA results. Available columns: ", paste(colnames(degs), collapse=", ")))
    }
    degs <- degs[!duplicated(degs$Gene_ID)&!is.na(degs$Gene_ID),]
    if(nrow(degs) == 0){
      return("RES-NO; No valid Gene_ID entries found in DEA results")
    }
    if(rank.stat == "coef"){
      ranks <- as.numeric(degs[,2])
    } else {
      if(!"T_statistic" %in% colnames(degs)){
        return(paste0("RES-NO; T_statistic column not found. Available columns: ", paste(colnames(degs), collapse=", ")))
      }
      ranks <- as.numeric(degs$T_statistic)
    }
    names(ranks) <- degs$Gene_ID
    ranks <- ranks[!is.na(ranks)]
    if(length(ranks) == 0){
      return("RES-NO; No valid ranks computed from DEA results")
    }
    print(paste("GSEA: Using", length(ranks), "ranked genes for", funcLib, "library"))
  }

  set.names <- data.frame(names = libraryRDS$term, IDs = names(libraryRDS$sets))

  # Set minSize based on library type
  if(funcLib %in% c("gem_react", "hsa_kegg", "gem_metabolite")){
    minSize <- 3
  } else {
    minSize <- 15
  }

  # perform analysis
  print(paste("Running GSEA with", length(libraryList), "pathways and", length(ranks), "ranked features"))
 
  all_lib_rxns <- unique(unlist(libraryList))
  overlap <- intersect(all_lib_rxns, names(ranks))
 
  # Check pathway sizes after intersection
  pw_sizes <- sapply(libraryList, function(pw) length(intersect(pw, names(ranks))))
  
  res <- tryCatch({
    fgsea(pathways = libraryList,
          stats = ranks,
          minSize = minSize,
          maxSize = 500)
  }, error = function(e){
    print(paste("Error in fgsea:", e$message))
    return(NULL)
  })

  if(is.null(res)){
    return("RES-NO")
  }

  print(paste("GSEA completed with", nrow(res), "results"))

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

  if(!requireNamespace("fgsea", quietly = TRUE)){
    return("RES-NO; fgsea package is not installed")
  }
  library(fgsea)
  library(data.table)
print(c(funcLib,fdr,collapse,mode))
  fdr <- as.numeric(fdr)

  # Methylation pathway analysis (ORA): CpG-level dea_results.csv -> gsameth
  # (Wallenius CpG-count bias correction). Self-contained early return; the
  # gene-level ORA path below is untouched for every other omic.
  .mhdr <- tryCatch(colnames(data.table::fread("dea_results.csv", nrows = 0)),
                    error = function(e) character(0))
  if(all(c("gene", "pos_hg38") %in% .mhdr) && !("Gene_ID" %in% .mhdr)){
    return(.runMethylationORA(funcLib, fdr, mode))
  }

  # set library filepath

  # read in DEGs
  if(!file.exists("dea_results.csv")){
    return("RES-NO; dea_results.csv not found. Please run differential analysis first.")
  }
  degs <- read.csv("dea_results.csv")
  if(nrow(degs) == 0){
    return("RES-NO; dea_results.csv is empty")
  }
  print(paste("ORA: Loaded", nrow(degs), "features. Columns:", paste(colnames(degs), collapse=", ")))

  # Determine universe and hits based on library type
  if(funcLib == "hsa_kegg"){
    if(!"kegg_id" %in% colnames(degs)){
      return("RES-NO; kegg_id column not found in DEA results. Run metabolite regression first.")
    }
    universe <- unique(as.character(degs$kegg_id[degs$kegg_id != "" & !is.na(degs$kegg_id)]))
    hits <- unique(as.character(degs$kegg_id[degs$sig != "NS" & degs$kegg_id != "" & !is.na(degs$kegg_id)]))
  } else if(funcLib == "gem_metabolite"){
    if(!"gem_id" %in% colnames(degs)){
      return("RES-NO; gem_id column not found in DEA results. Run metabolite regression first.")
    }
    universe <- unique(as.character(degs$gem_id[degs$gem_id != "" & !is.na(degs$gem_id)]))
    hits <- unique(as.character(degs$gem_id[degs$sig != "NS" & degs$gem_id != "" & !is.na(degs$gem_id)]))
  } else if(funcLib == "gem_react" && "Feature" %in% colnames(degs) && !("Gene_ID" %in% colnames(degs))){
    universe <- unique(as.character(degs$Feature))
    hits <- unique(as.character(degs$Feature[degs$sig != "NS"]))
  } else {
    if(!"Gene_ID" %in% colnames(degs)){
      return(paste0("RES-NO; Gene_ID column not found in DEA results. Available columns: ", paste(colnames(degs), collapse=", ")))
    }
    # Drop features whose Gene_ID (Entrez) is missing, for EVERY omics -- newer fgsea's fora()
    # errors with "NAs in universe are not allowed". This used to be gated on Protein_Group, i.e.
    # proteomics v2 only, which is why patch-seq broke: PatchseqSpearman's v2 branch keeps genes
    # on a non-empty SYMBOL, so the 1,813 patch-seq genes that have a symbol but no Entrez
    # (lncRNAs etc.) carried Gene_ID = NA straight into the universe. performGSEA already does
    # exactly this for all omics (`degs[!duplicated(degs$Gene_ID) & !is.na(degs$Gene_ID),]`);
    # ORA was the outlier.
    # The gene-set libraries are Entrez-keyed (go_bp/kegg sets are Entrez integers), so Entrez
    # remains the only possible pathway key. Only the Entrez already present in the data is used
    # -- no external annotation is consulted to back-fill missing IDs. A feature with no Entrez
    # cannot belong to any set, so it belongs in neither the background nor the hits.
    universe <- unique(as.character(degs$Gene_ID[degs$Gene_ID != "" & !is.na(degs$Gene_ID)]))
    hits <- unique(as.character(degs$Gene_ID[degs$sig != "NS" & degs$Gene_ID != "" & !is.na(degs$Gene_ID)]))
  }
  
 
  if(length(hits) > 10){
    # process libraries
    if(mode=="tool"){
 lib.path <- paste0(other.tables.path, "libraries/")

if(funcLib == "hsa_kegg"){
  lib <- qs::qread(paste0(lib.path, "kegg_hsa_met.qs"))
  libraryList <- lib$mset.list                          # keyed by KEGG pathway IDs (hsa00010, ...)
  name_lookup  <- setNames(names(lib$path.ids), lib$path.ids)  # hsa00010 → "Glycolysis..."
  readable_names <- name_lookup[names(libraryList)]
  readable_names[is.na(readable_names)] <- names(libraryList)[is.na(readable_names)]
  libraryRDS <- list(term = unname(readable_names), sets = libraryList)
  minS <- 3
  set.names <- data.frame(names = unname(readable_names), IDs = names(libraryList), stringsAsFactors = FALSE)
  tmp <- degs[degs$kegg_id != "" & !is.na(degs$kegg_id), ]
  tmp <- tmp[!duplicated(tmp$kegg_id), ]
  degs2 <- data.frame(Gene_ID = tmp$kegg_id, Log2FC = tmp[, 2], stringsAsFactors = FALSE)

}else if(funcLib == "gem_metabolite"){
  lib <- qs::qread(paste0(lib.path, "hsaGem_metabolite.qs"))
  libraryList <- lib$sets
  libraryRDS <- lib
  minS <- 3
  set.names <- data.frame(names = lib$term, IDs = names(lib$sets), stringsAsFactors = FALSE)
  tmp <- degs[degs$gem_id != "" & !is.na(degs$gem_id), ]
  tmp <- tmp[!duplicated(tmp$gem_id), ]
  degs2 <- data.frame(Gene_ID = tmp$gem_id, Log2FC = tmp[, 2], stringsAsFactors = FALSE)

}else if(funcLib=="gem_react"){

# Check if universe contains reaction IDs (flux data) or gene IDs
is_flux_data <- "Feature" %in% colnames(degs) && !("Gene_ID" %in% colnames(degs))

if(is_flux_data){
  # For flux data, universe and hits are already reaction IDs
  # No gene conversion needed - use reactions directly for pathway analysis
  lib = qs::qread(paste0(lib.path,"hsaGem_reaction.qs"))
  rxn_tab <- lib$reacTab

  # Filter rxn_tab to only include reactions in our universe
  rxn_tab <- rxn_tab[rxn_tab$humanGEM %in% universe, ]
  rxn_tab$universe <- TRUE
  rxn_tab$hits <- rxn_tab$humanGEM %in% hits

  libraryList <- split(rxn_tab$humanGEM, rxn_tab$sysID)
  universe <- rxn_tab$humanGEM
  hits <- rxn_tab$humanGEM[rxn_tab$hits]
  minS <- 3
  libraryRDS <- rxn_tab
  set.names <- unique(data.frame(names = rxn_tab$pathwayVis, IDs = rxn_tab$sysID))

  # Use Feature column for fold changes
  gene_fc <- setNames(degs$Log2FC, degs$Feature)

  # Create degs2 structure for flux data - use reactions directly
  degs2 <- data.frame(
    Gene_ID = degs$Feature,
    Log2FC = degs$Log2FC,
    stringsAsFactors = FALSE
  )

} else {
  # Gene-based data: convert Feature to ENSG for GPR rule matching
  out_anno <- qs::qread(paste0(other.tables.path,"libraries/gene_anno_full.qs"))

  if(all(grepl("ENSG", head(degs$Feature)))){
    # Feature already has ENSG (RNA-seq) - use directly
    universe <- unique(as.character(degs$Feature))
    hits <- unique(as.character(degs$Feature[degs$sig != "NS"]))
  } else {
    # Feature has symbols (proteomics, pbrna) - convert symbol to ENSG
    universe <- out_anno$ensembl_gene_id[match(unique(as.character(degs$Feature)), out_anno$hgnc_symbol)]
    hits <- out_anno$ensembl_gene_id[match(unique(as.character(degs$Feature[degs$sig != "NS"])), out_anno$hgnc_symbol)]
  }
  universe <- universe[!is.na(universe)]
  hits <- hits[!is.na(hits)]

  lib = qs::qread(paste0(lib.path,"hsaGem_reaction.qs"))

  rxn_tab <- lib$reacTab
  rules<- lib$rules

  universe_map <- lapply(rules,function(r) eval_gpr(r, universe))
   rxn_ids <- sub("\\..*$", "", names(universe_map))
   rxn_tab <- rxn_tab[match(rxn_ids,rxn_tab$humanGEM),]
   rxn_tab$universe <- unlist(lapply(universe_map,function(x) x[["map"]]))
  kp = rxn_tab$universe
   rxn_tab <- rxn_tab[kp,]
   universe_map <- universe_map[kp]

   hit_rxn_names <- names(universe_map)
   hit_rules <- rules[hit_rxn_names]
   hits_map <- lapply(hit_rules,function(r) eval_gpr(r, hits))
   rxn_tab$hits <- unlist(lapply(hits_map,function(x) x[["map"]]))

   libraryList <- split(rxn_tab$humanGEM, rxn_tab$sysID)
  universe     <- rxn_tab$humanGEM
  hits          <- rxn_tab$humanGEM[rxn_tab$hits]
  minS <- 3
  libraryRDS <- rxn_tab
      set.names <- unique(data.frame(names = rxn_tab$pathwayVis, IDs =  rxn_tab$sysID))

  # Map degs to ENSG for fold change lookup
  if(all(grepl("ENSG", head(degs$Feature)))){
    degs <- degs[!duplicated(degs$Feature)&!is.na(degs$Feature),]
    gene_fc <- setNames(degs[,2], degs$Feature)
  } else {
    degs <- degs[!duplicated(degs$Feature)&!is.na(degs$Feature),]
    degs$accesion <- out_anno$ensembl_gene_id[match(degs$Feature, out_anno$hgnc_symbol)]
    degs <- degs[!is.na(degs$accesion),]
    gene_fc <- setNames(degs[,2], degs$accesion)
  }

   library(stringr)
   degs2 = lapply(universe_map,function(x){
     map = x[[2]]
     map =map[map$map,]
     for(i in 1:nrow(map)){
       if(grepl("\\&",map$id[i])){
         fc=  str_split(map$id[i], "&")[[1]] %>%
           str_trim()
         fc_vals <- gene_fc[fc]
         kp <- which.min(abs(fc_vals))
         map$id[i] <- names(fc_vals)[kp]

       }

     }
     map$Log2FC  = gene_fc[map$id]
     map[which.max(abs(map$Log2FC)),c(1,4)]
    })
    degs2 <- do.call(
     rbind,
     lapply(names(degs2), function(nm) {
       df <- degs2[[nm]]
       cbind(rxn = nm, df)
     })
   )
   degs2$id = NULL
   names(degs2) = c("Gene_ID","Log2FC")
   # Strip gene suffix from compound reaction names
   degs2$Gene_ID <- sub("\\..*$", "", degs2$Gene_ID)
} 

}else{

  lib_file <- paste0(lib.path, funcLib, ".rds")
  if(!file.exists(lib_file)){
    return(paste0("RES-NO; Library file not found: ", lib_file))
  }
       libraryRDS <- readRDS(lib_file)
       libraryList <- libraryRDS$sets

minS <- 15
    set.names <- data.frame(names = libraryRDS$term, IDs = names(libraryRDS$sets))
degs2 =degs;
}
 
   if(!file.exists("savedAnalysis")){
    dir.create("savedAnalysis")
    cat("## R history", file = "savedAnalysis/Rhistory.R", sep = "\n")
  }


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


    # perform analysis
    print(paste("ORA: Running fora with", length(hits), "hits,", length(universe), "universe genes,", length(libraryList), "pathways"))
    res <- tryCatch({
      fora(pathways = libraryList,
                genes = hits,
                universe = universe,
                minSize = minS, maxSize = 500)
    }, error = function(e){
      return(paste0("ERROR:", e$message))
    })
    if(is.character(res) && grepl("^ERROR:", res)){
      return(paste0("RES-NO; fora() failed: ", gsub("^ERROR:", "", res)))
    }

    # find independent pathways if collapse set to true
    if(collapse == "true"){
      res_precollapse <- res
      ind.pathways <- tryCatch({
        collapsePathwaysORA(foraRes = res,
                                       pathways = libraryList,
                                       genes = hits,
                                       universe = universe)
      }, error = function(e){
        print(paste("collapsePathwaysORA failed:", e$message))
        return(NULL)
      })
      if(!is.null(ind.pathways)){
        res <- res[res$pathway %in% ind.pathways$mainPathways]
      }
      # Never blank the results table: if collapsing removed every row (id-space
      # mismatch / edge case), fall back to the full un-collapsed set.
      if(nrow(as.data.frame(res)) == 0){ res <- res_precollapse }
    }
    res <- as.data.frame(res)

    # format results
    res <- tryCatch({
      res <- merge(res, set.names, by.x = "pathway", by.y = "IDs")
      res <- res[order(res$pval),]
      res$overlapGenes <- NULL # remove list of genes
      res <- res[,c("names" , "pval" , "padj" , "overlap" ,"size" ,"pathway")]
      res[,2:3] <- signif(res[,2:3], digits = 3)
      colnames(res) <- c("Set Name", "P_value", "Adj P_value", "Hits", "Set Size", "Set ID")
      res
    }, error = function(e){
      return(paste0("RES-NO; ORA result formatting failed: ", e$message))
    })
    if(is.character(res) && grepl("^RES-NO", res)){
      return(res)
    }

   write.csv(res, "functional_results.csv", row.names = FALSE)
    num.sig <- sum(res$`Adj P_value` < fdr, na.rm = TRUE)

    # format results for ridgeline
    genesPW <- tryCatch({
      genesPW <- reshape::melt(libraryList, level = 1)
      colnames(genesPW) <- c("GeneID", "Set ID")

      res$pathNegLogP <- -log10(res$P_value)
      res$RankSig <- rank(res$P_value)
      genesPW <- merge(genesPW, res, by = "Set ID", all.y = TRUE, all.x = FALSE)
      genesPW$Sig <- as.character(genesPW$`Adj P_value` < fdr)
      genesPW <- genesPW[,c("Set ID", "GeneID", "Set Name", "Set Size", "Adj P_value", "Sig", "pathNegLogP", "RankSig")]
      genesPW <- merge(genesPW, data.frame(GeneID = degs2$Gene_ID, GeneStat = degs2[,2]), by = "GeneID", all.x = TRUE, all.y = FALSE)
      genesPW <- na.omit(genesPW)
      genesPW <- as.data.table(genesPW)
      genesPW <- genesPW[,.(GeneID, `Set Name`, `Set Size`, pathFDR = `Adj P_value`, Sig, pathNegLogP, GeneStat, MeanGeneStat = mean(GeneStat), RankSig), by = .(`Set ID`)]
      genesPW <- genesPW[order(MeanGeneStat)]
      genesPW
    }, error = function(e){
      print(paste("Ridgeline formatting failed:", e$message))
      return(NULL)
    })
    if(!is.null(genesPW)){
      write.csv(genesPW, "ridgeline.csv", row.names = FALSE)
    }

    return(paste0("RES-OK;", num.sig))
  } else {
    return("RES-NO; fewer than 10 DEGs")
  }
  
}


################################################################################

## CAMERA-based pathway enrichment for metabolic flux data
## Uses saved limma fit from DEA step; accounts for intra-subsystem reaction correlation.
## Only called when omicsType == "proc_flux" + library == "gem_react".

performCAMERA <- function(funcLib = "gem_react", fdr = 0.05, mode = "local") {
  library(limma)
  library(data.table)

  fdr <- as.numeric(fdr)

  # ── Load saved limma objects from the DEA step ──────────────────────────────
  for(f in c("camera_fit.qs", "camera_features.qs", "camera_type.qs")){
    if(!file.exists(f)) return(paste0("RES-NO; CAMERA requires running differential analysis first (missing: ", f, ")"))
  }
  fit           <- qs::qread("camera_fit.qs")
  feature_table <- qs::qread("camera_features.qs")
  analysis_type <- qs::qread("camera_type.qs")

  if(analysis_type == "disc"){
    if(!file.exists("camera_contrast.qs")) return("RES-NO; camera_contrast.qs not found.")
    contrast_or_coef <- qs::qread("camera_contrast.qs")
    is_coef <- FALSE
  } else {
    if(!file.exists("camera_coef.qs")) return("RES-NO; camera_coef.qs not found.")
    contrast_or_coef <- qs::qread("camera_coef.qs")  # character: column name in design
    is_coef <- TRUE
  }

  # ── Load GEM reaction library ────────────────────────────────────────────────
  if(mode == "tool"){
    lib.path <- paste0(other.tables.path, "libraries/")
  } else {
    lib.path <- paste0(other.tables.path, "libraries/")
  }
  lib     <- qs::qread(paste0(lib.path, "hsaGem_reaction.qs"))
  rxn_tab <- lib$reacTab

  # ── Build pathway index list (reactions present in our data only) ────────────
  all_features  <- rownames(feature_table)
  rxn_filtered  <- rxn_tab[rxn_tab$humanGEM %in% all_features, ]
  pathway_list  <- split(rxn_filtered$humanGEM, rxn_filtered$sysID)
  subsys_names  <- unique(data.frame(sysID = rxn_filtered$sysID,
                                     pathwayVis = rxn_filtered$pathwayVis,
                                     stringsAsFactors = FALSE))

  pathway_indices <- lapply(pathway_list, function(ids){
    idx <- match(ids, all_features); idx[!is.na(idx)]
  })
  pathway_indices <- pathway_indices[sapply(pathway_indices, length) >= 3]

  if(length(pathway_indices) == 0)
    return("RES-NO; No subsystems with >= 3 reactions found in the data.")

  # ── Remove non-finite rows; re-fit if any were removed ──────────────────────
  row_ok <- apply(feature_table, 1, function(x) all(is.finite(x)))
  if(any(!row_ok)){
    feature_table   <- feature_table[row_ok, , drop = FALSE]
    all_features    <- rownames(feature_table)
    pathway_indices <- lapply(pathway_list, function(ids){
      idx <- match(ids, all_features); idx[!is.na(idx)]
    })
    pathway_indices <- pathway_indices[sapply(pathway_indices, length) >= 3]
    fit <- lmFit(feature_table, fit$design)
  }

  # ── Run CAMERA ───────────────────────────────────────────────────────────────
  cam_res <- tryCatch({
    if(is_coef){
      coef_idx <- which(colnames(fit$design) == contrast_or_coef)
      if(length(coef_idx) == 0) coef_idx <- 2  # fallback to 2nd column
      camera(feature_table, index = pathway_indices,
             design = fit$design, coef = coef_idx)
    } else {
      camera(feature_table, index = pathway_indices,
             design = fit$design, contrast = contrast_or_coef)
    }
  }, error = function(e){
    return(paste0("RES-NO; CAMERA error: ", e$message))
  })
  if(is.character(cam_res)) return(cam_res)
  if(is.null(cam_res) || nrow(cam_res) == 0) return("RES-NO; CAMERA returned no results.")

  # ── Per-reaction effect size for mean score and ridgeline ───────────────────
  if(is_coef){
    feature_stats <- fit$coefficients[, contrast_or_coef, drop = TRUE]
  } else {
    fit2          <- contrasts.fit(fit, contrast_or_coef)
    feature_stats <- fit2$coefficients[, 1]
  }
  names(feature_stats) <- rownames(feature_table)

  # ── Compute direction-signed mean score per pathway ─────────────────────────
  matched_sysIDs <- intersect(rownames(cam_res), names(pathway_list))
  pw_means <- sapply(matched_sysIDs, function(pw){
    members <- intersect(pathway_list[[pw]], names(feature_stats))
    if(length(members) == 0) return(0)
    m <- mean(feature_stats[members], na.rm = TRUE)
    if(cam_res[pw, "Direction"] == "Down") -abs(m) else abs(m)
  })

  # ── Format functional_results.csv (6 cols, same structure as GSEA/ORA) ──────
  res <- data.frame(
    `Set Name`    = subsys_names$pathwayVis[match(matched_sysIDs, subsys_names$sysID)],
    P_value       = cam_res[matched_sysIDs, "PValue"],
    `Adj P_value` = cam_res[matched_sysIDs, "FDR"],
    Mean_score    = round(pw_means, 4),
    Set_Size      = cam_res[matched_sysIDs, "NGenes"],
    Set_ID        = matched_sysIDs,
    stringsAsFactors = FALSE,
    check.names   = FALSE
  )
  res <- res[!is.na(res[["Set Name"]]), ]
  res <- res[order(res$P_value), ]
  res[, c("P_value","Adj P_value")] <- signif(res[, c("P_value","Adj P_value")], digits = 3)
  write.csv(res, "functional_results.csv", row.names = FALSE)
  num_sig <- sum(res[["Adj P_value"]] < fdr, na.rm = TRUE)

  # ── Format ridgeline.csv (same columns as GSEA/ORA ridgeline) ───────────────
  res$pathNegLogP <- -log10(res$P_value)
  res$RankSig     <- rank(res$P_value)

  genesPW <- tryCatch({
    genesPW <- reshape::melt(pathway_list[matched_sysIDs], level = 1)
    colnames(genesPW) <- c("GeneID", "Set_ID")        # underscore — matches res column name
    genesPW <- merge(genesPW, res[, c("Set_ID","Set Name","Set_Size","Adj P_value","pathNegLogP","RankSig")],
                     by = "Set_ID", all.x = TRUE)
    genesPW$Sig <- as.character(genesPW[["Adj P_value"]] < fdr)
    genesPW <- merge(genesPW,
                     data.frame(GeneID = names(feature_stats), GeneStat = unname(feature_stats)),
                     by = "GeneID", all.x = TRUE)
    genesPW <- na.omit(genesPW)
    genesPW <- as.data.table(genesPW)
    # Group by Set_ID (kept in ridgeline.csv for internal deduplication); output Set Size with space
    genesPW <- genesPW[, .(GeneID, `Set Name`, `Set Size` = Set_Size, pathFDR = `Adj P_value`,
                            Sig, pathNegLogP, GeneStat,
                            MeanGeneStat = mean(GeneStat), RankSig), by = .(`Set ID` = Set_ID)]
    genesPW[order(MeanGeneStat)]
  }, error = function(e){
    print(paste("CAMERA ridgeline formatting failed:", e$message)); NULL
  })
  if(!is.null(genesPW)) write.csv(genesPW, "ridgeline.csv", row.names = FALSE)

  if(mode == "tool" && !file.exists("savedAnalysis")){
    dir.create("savedAnalysis")
  }
  if(mode == "tool"){
    cat("## R history", file = "savedAnalysis/Rhistory.R", sep = "\n", append = TRUE)
    dump("performCAMERA", file = "savedAnalysis/performCAMERA.R")
  }

  return(paste0("RES-OK;", num_sig))
}

################################################################################

## Mummichog-based pathway enrichment for untargeted metabolomics
## Adapted from MetaboAnalyst's mummichog algorithm
## Uses m/z to compound matching via adduct mass lookup, then Fisher's exact test

performMummichog <- function(funcLib = "hsa_kegg", fdr = 0.05, permNum = 100,
                             instrument = 10, mode = "local"){

  library(qs)
  library(data.table)

  fdr <- as.numeric(fdr)
  permNum <- as.numeric(permNum)
  instrument <- as.numeric(instrument) # mass accuracy in ppm

  # Read DEA results (already saved by DonorRegression)
  degs <- read.csv("dea_results.csv")

  # Check if this is metabolomics data (has mode column)
  if(!"mode" %in% colnames(degs)){
    return("RES-NO; Not metabolomics data")
  }

  # Extract m/z and retention time from peak names (format: mz__rt)
  peak_parts <- strsplit(as.character(degs$Feature), "__", fixed = TRUE)
  degs$mz <- as.numeric(sapply(peak_parts, `[`, 1))
  degs$rt <- as.numeric(sapply(peak_parts, `[`, 2))
  degs <- degs[!is.na(degs$mz),]

  # Get p-values and determine significance
  if("P_value" %in% colnames(degs)){
    degs$pval <- degs$P_value
  } else {
    return("RES-NO; No P_value column found")
  }

  # Define significant peaks (p < 0.05 by default from DEA cutoff)
  sig_peaks <- degs$pval < 0.05

  if(sum(sig_peaks) < 3){
    return("RES-NO; Too few significant peaks for pathway analysis")
  }

  # Load mummichog KEGG library
  if(mode == "tool"){
    lib.path <- paste0(other.tables.path, "libraries/")
  } else {
    lib.path <- "libraries/"
  }

  lib_file <- paste0(lib.path, funcLib, ".qs")
  if(!file.exists(lib_file)){
    return(paste0("RES-NO; Library file not found: ", lib_file))
  }

  mum.lib <- qs::qread(lib_file)

  # Extract library components
  cpd.lib <- mum.lib$cpd.lib
  cpd.treep <- mum.lib$cpd.tree[["positive"]]
  cpd.treen <- mum.lib$cpd.tree[["negative"]]
  pathways <- mum.lib$pathways

  # Adduct mass matrices
  mz.matp <- cpd.lib$adducts[["positive"]]
  mz.matn <- cpd.lib$adducts[["negative"]]

  # Determine ion mode from data
  peak_modes <- unique(degs$mode)
  has_positive <- "positive" %in% peak_modes
  has_negative <- "negative" %in% peak_modes

  if(has_positive && has_negative){
    data_mode <- "mixed"
  } else if(has_positive){
    data_mode <- "positive"
  } else {
    data_mode <- "negative"
  }

  # --- Step 1: Match m/z to compounds using adduct mass lookup ---
  ref_mzlist <- degs$mz
  pos_inx <- degs$mode == "positive"
  ref_mzlistp <- ref_mzlist[pos_inx]
  ref_mzlistn <- ref_mzlist[!pos_inx]

  # m/z tolerance function
  mz_tol <- function(mz, ppm) ppm * 1e-06 * mz

  # Match positive mode peaks
  matched_list <- list()

  if(data_mode != "negative" && length(ref_mzlistp) > 0){
    my.tolsp <- mz_tol(ref_mzlistp, instrument)
    self.mzsp <- floor(ref_mzlistp)
    all.mzsp <- cbind(self.mzsp - 1, self.mzsp, self.mzsp + 1)
    modified.statesp <- colnames(mz.matp)
    tree_len_p <- length(cpd.treep)

    for(i in seq_along(ref_mzlistp)){
      mz <- ref_mzlistp[i]
      my.tol <- my.tolsp[i]
      all.mz <- all.mzsp[i,]
      all.mz <- all.mz[all.mz >= 1 & all.mz <= tree_len_p]
      if(length(all.mz) == 0) next
      pos.all <- as.numeric(unique(unlist(cpd.treep[all.mz])))
      pos.all <- pos.all[!is.na(pos.all) & pos.all > 0]

      if(length(pos.all) > 0){
        for(pos in pos.all){
          id <- cpd.lib$id[pos]
          mw.all <- mz.matp[pos,]
          diffs <- abs(mw.all - mz)
          hit.inx <- which(diffs < my.tol)
          if(length(hit.inx) > 0){
            for(spot in seq_along(hit.inx)){
              hit.pos <- hit.inx[spot]
              matched_list[[length(matched_list) + 1]] <- c(
                Query.Mass = mz,
                Matched.Compound = id,
                Matched.Form = modified.statesp[hit.pos],
                Mass.Diff = diffs[hit.pos]
              )
            }
          }
        }
      }
    }
  }

  # Match negative mode peaks
  if(data_mode != "positive" && length(ref_mzlistn) > 0){
    my.tolsn <- mz_tol(ref_mzlistn, instrument)
    self.mzsn <- floor(ref_mzlistn)
    all.mzsn <- cbind(self.mzsn - 1, self.mzsn, self.mzsn + 1)
    modified.statesn <- colnames(mz.matn)
    tree_len_n <- length(cpd.treen)

    for(i in seq_along(ref_mzlistn)){
      mz <- ref_mzlistn[i]
      my.tol <- my.tolsn[i]
      all.mz <- all.mzsn[i,]
      all.mz <- all.mz[all.mz >= 1 & all.mz <= tree_len_n]
      if(length(all.mz) == 0) next
      pos.all <- as.numeric(unique(unlist(cpd.treen[all.mz])))
      pos.all <- pos.all[!is.na(pos.all) & pos.all > 0]

      if(length(pos.all) > 0){
        for(pos in pos.all){
          id <- cpd.lib$id[pos]
          mw.all <- mz.matn[pos,]
          diffs <- abs(mw.all - mz)
          hit.inx <- which(diffs < my.tol)
          if(length(hit.inx) > 0){
            for(spot in seq_along(hit.inx)){
              hit.pos <- hit.inx[spot]
              matched_list[[length(matched_list) + 1]] <- c(
                Query.Mass = mz,
                Matched.Compound = id,
                Matched.Form = modified.statesn[hit.pos],
                Mass.Diff = diffs[hit.pos]
              )
            }
          }
        }
      }
    }
  }

  if(length(matched_list) == 0){
    return("RES-NO; No compound matches found from peak list")
  }

  matched_res <- data.frame(do.call(rbind, matched_list), stringsAsFactors = FALSE)
  matched_res$Query.Mass <- as.numeric(matched_res$Query.Mass)
  matched_res$Mass.Diff <- as.numeric(matched_res$Mass.Diff)

  print(paste0("Matched ", length(unique(matched_res$Matched.Compound)), " unique compounds from ", nrow(matched_res), " m/z-compound pairs"))

  # --- Step 2: Build dictionaries ---
  # mz -> compound dictionary
  mz2cpd <- split(matched_res$Matched.Compound, as.character(matched_res$Query.Mass))
  mz2cpd <- lapply(mz2cpd, unique)

  # compound -> mz dictionary
  cpd2mz <- split(as.character(matched_res$Query.Mass), matched_res$Matched.Compound)
  cpd2mz <- lapply(cpd2mz, unique)

  # All matched compounds
  total_matched_cpds <- unique(matched_res$Matched.Compound)

  # --- Step 3: Identify significant compound list ---
  sig_mz <- degs$mz[sig_peaks]

  # Map significant m/z to compounds using the same character keys as dictionaries
  all_matched_mz <- unique(matched_res$Query.Mass)
  sig_matched_mz <- all_matched_mz[all_matched_mz %in% sig_mz]
  sig_matched_mz_char <- as.character(sig_matched_mz)

  input_cpdlist <- unique(unlist(mz2cpd[sig_matched_mz_char]))
  input_cpdlist <- input_cpdlist[!is.null(input_cpdlist)]

  # Also get the significant m/z list for permutation (as character to match dict keys)
  input_mzlist <- sig_matched_mz_char
  N <- length(input_mzlist)

  print(paste0("Significant compounds: ", length(input_cpdlist), " from ", N, " significant peaks"))

  if(length(input_cpdlist) < 3){
    return("RES-NO; Too few compound matches from significant peaks")
  }

  # --- Step 4: Filter pathways by minimum size ---
  minLib <- 3
  path.length <- sapply(pathways$cpds, length)
  min.inx <- which(path.length >= minLib)

  cleaned.pathways <- list()
  cleaned.pathways$cpds <- pathways$cpds[min.inx]
  cleaned.pathways$name <- pathways$name[min.inx]

  # --- Step 5: Compute observed pathway enrichment (Fisher's exact test) ---
  qset <- unique(input_cpdlist)
  query_set_size <- length(qset)
  total_feature_num <- length(total_matched_cpds)

  current.mset <- cleaned.pathways$cpds
  path.num <- unlist(lapply(current.mset, length))

  # Compounds in each pathway that overlap with all matched compounds
  cpds <- lapply(current.mset, function(x) intersect(x, total_matched_cpds))
  set.num <- unlist(lapply(cpds, length))

  # Compounds in each pathway that overlap with significant compounds
  feats <- lapply(current.mset, function(x) intersect(x, qset))
  feat_len <- unlist(lapply(feats, length))

  # Adjust for multiple m/z per compound
  count_cpd2mz_local <- function(cpd.ids, inputmzlist){
    if(length(cpd.ids) == 0) return(0)
    mzs <- as.numeric(unique(unlist(cpd2mz[cpd.ids])))
    if(length(mzs) == 0) return(0)
    return(length(intersect(as.character(mzs), inputmzlist)))
  }

  sizes <- negneg <- vector("list", length(current.mset))
  for(i in seq_along(current.mset)){
    sizes[[i]] <- min(feat_len[i], count_cpd2mz_local(unlist(feats[i]), input_mzlist))
    negneg[[i]] <- total_feature_num + sizes[[i]] - set.num[i] - query_set_size
  }
  negneg <- rapply(negneg, function(x) ifelse(x < 0, 0, x), how = "replace")

  unsize <- as.integer(unlist(sizes))

  # Fisher's exact test
  fishermatrix <- cbind(unsize - 1, set.num, (query_set_size + unlist(negneg) - unsize), query_set_size)
  fisher.p <- apply(fishermatrix, 1, function(x) phyper(x[1], x[2], x[3], x[4], lower.tail = FALSE))

  # EASE score (modified Fisher)
  first <- unlist(lapply(sizes, function(x) max(0, x - 1)))
  easematrix <- cbind(first, (set.num - unsize + 1), (query_set_size - unsize), unlist(negneg))
  ease.p <- apply(easematrix, 1, function(x) {
    tryCatch(fisher.test(matrix(x, nrow = 2), alternative = "greater")$p.value, error = function(e) 1)
  })

  # --- Step 6: Permutation for background distribution ---
  print(paste0("Running ", permNum, " permutations..."))
  set.seed(123)

  perm_record <- vector("list", permNum)
  perm_hits <- vector("list", permNum)
  # Sample from all m/z features in the dataset (most won't match compounds - that's expected)
  all_ref_mz <- as.character(degs$mz)

  for(p in 1:permNum){
    set.seed(123 + p)
    perm_mz <- sample(all_ref_mz, N)
    perm_cpds <- unique(unlist(mz2cpd[perm_mz]))
    perm_cpds <- perm_cpds[!is.null(perm_cpds)]

    perm_feats <- lapply(current.mset, function(x) intersect(x, perm_cpds))
    perm_feat_len <- unlist(lapply(perm_feats, length))
    perm_query_size <- length(perm_cpds)

    perm_sizes <- perm_negneg <- vector("list", length(current.mset))
    for(i in seq_along(current.mset)){
      perm_sizes[[i]] <- min(perm_feat_len[i], count_cpd2mz_local(unlist(perm_feats[i]), perm_mz))
      perm_negneg[[i]] <- total_feature_num + perm_sizes[[i]] - set.num[i] - perm_query_size
    }
    perm_negneg <- rapply(perm_negneg, function(x) ifelse(x < 0, 0, x), how = "replace")

    perm_unsize <- as.integer(unlist(perm_sizes))
    perm_fisher <- cbind(perm_unsize - 1, set.num, (perm_query_size + unlist(perm_negneg) - perm_unsize), perm_query_size)
    perm_p <- apply(perm_fisher, 1, function(x) phyper(x[1], x[2], x[3], x[4], lower.tail = FALSE))

    perm_record[[p]] <- perm_p
    perm_hits[[p]] <- perm_unsize
  }

  # --- Step 7: Compute empirical p-values ---
  record_matrix <- do.call(cbind, perm_record)
  emp.p <- sapply(seq_along(fisher.p), function(i) sum(record_matrix[i,] <= fisher.p[i]) / permNum)

  # Gamma-adjusted p-values
  perm_minus <- abs(0.9999999999 - unlist(perm_record))
  gamma.p <- tryCatch({
    fit.gamma <- fitdistrplus::fitdist(perm_minus, distr = "gamma", method = "mle",
                                        lower = c(0, 0), start = list(scale = 1, shape = 1))
    1 - pgamma(1 - ease.p, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["scale"])
  }, error = function(e){
    rep(NA, length(ease.p))
  })

  # --- Step 8: Build result matrix ---
  feat_vec <- sapply(feats, function(x) paste(x, collapse = ";"))

  res.mat <- data.frame(
    Pathway = cleaned.pathways$name,
    Pathway_Total = path.num,
    Hits_Total = set.num,
    Hits_Sig = unsize,
    Expected = round(query_set_size * (path.num / length(unique(unlist(current.mset)))), 2),
    FET = signif(fisher.p, 5),
    EASE = signif(ease.p, 5),
    Gamma = signif(gamma.p, 5),
    Empirical = signif(emp.p, 5),
    stringsAsFactors = FALSE
  )

  # Filter to pathways with at least 1 significant hit
  hit.inx <- res.mat$Hits_Sig > 0
  if(sum(hit.inx) < 1){
    return("RES-NO; No pathways with significant compound hits")
  }
  res.mat <- res.mat[hit.inx,]

  # Order by Gamma p-value (or Fisher if gamma failed)
  if(all(is.na(res.mat$Gamma))){
    res.mat <- res.mat[order(res.mat$FET),]
  } else {
    res.mat <- res.mat[order(res.mat$Gamma),]
  }

  # Format output: Set Name, P(Fisher), P(Gamma), Hits, Set Size, Set ID
  output <- data.frame(
    `Set Name` = res.mat$Pathway,
    `P(Fisher)` = signif(res.mat$FET, 5),
    `P(Gamma)` = signif(res.mat$Gamma, 5),
    Hits = res.mat$Hits_Sig,
    `Set Size` = res.mat$Hits_Total,
    `Set ID` = res.mat$Pathway,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  write.csv(output, "functional_results.csv", row.names = FALSE)
  num.sig <- sum(output$`P(Gamma)` < fdr, na.rm = TRUE)

  # --- Step 9: Generate ridgeline data ---
  # For metabolomics, use compound-level stats mapped through pathways
  # Map compounds to t-statistics via m/z
  cpd_stats <- sapply(total_matched_cpds, function(cpd){
    mzs <- cpd2mz[[cpd]]
    if(length(mzs) == 0) return(NA)
    matched_degs <- degs[as.character(degs$mz) %in% mzs,]
    if(nrow(matched_degs) == 0) return(NA)
    # Use the most significant m/z for this compound
    matched_degs$T_statistic[which.min(matched_degs$pval)]
  })
  cpd_stats <- cpd_stats[!is.na(cpd_stats)]

  # Build ridgeline data
  pathway_cpd_list <- cleaned.pathways$cpds[hit.inx]
  pathway_names_list <- cleaned.pathways$name[hit.inx]

  ridge_rows <- list()
  for(j in seq_along(pathway_cpd_list)){
    pw_cpds <- intersect(unlist(pathway_cpd_list[j]), names(cpd_stats))
    if(length(pw_cpds) > 0){
      pw_name <- pathway_names_list[j]
      pw_idx <- match(pw_name, output$`Set Name`)
      if(!is.na(pw_idx)){
        for(cpd in pw_cpds){
          ridge_rows[[length(ridge_rows) + 1]] <- data.frame(
            `Set ID` = pw_name,
            GeneID = cpd,
            `Set Name` = pw_name,
            `Set Size` = output$`Set Size`[pw_idx],
            pathFDR = output$`P(Gamma)`[pw_idx],
            Sig = as.character(output$`P(Gamma)`[pw_idx] < fdr),
            pathNegLogP = -log10(output$`P(Fisher)`[pw_idx]),
            GeneStat = cpd_stats[cpd],
            MeanGeneStat = 0,
            RankSig = pw_idx,
            check.names = FALSE,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if(length(ridge_rows) > 0){
    ridgeData <- do.call(rbind, ridge_rows)
    # Calculate MeanGeneStat per pathway
    ridgeData <- as.data.table(ridgeData)
    ridgeData[, MeanGeneStat := mean(GeneStat, na.rm = TRUE), by = `Set ID`]
    ridgeData <- ridgeData[order(MeanGeneStat)]
    write.csv(ridgeData, "ridgeline.csv", row.names = FALSE)
  } else {
    # Create empty ridgeline file
    write.csv(data.frame(), "ridgeline.csv", row.names = FALSE)
  }

  return(paste0("RES-OK;", num.sig))
}


################################################################################

## Region-level (DMR) analysis for methylation. Reads the per-CpG limma result that the
## feature-level analysis ALREADY computed + saved (dea_results.csv: T_statistic,
## Coefficient, P_value, Adjusted p_value, chr, pos_hg38) and clusters it into regions
## with DMRcate. It does NOT reload the M matrix or re-fit limma. Our own coordinates
## avoid any EPICv2 probe-ID matching. Writes dmr_results.csv, returns "RES-OK;<n>".
## Validated standalone in humanislets/methylation/dmr_methylation.R.
.runMethylationDMRfromSaved <- function(fdr.thresh, lambda = 1000, C = 2){
  if(!file.exists("dea_results.csv")){
    return("RES-NO; please run the per-CpG (feature-level) analysis first, then region analysis")
  }
  # DMRcate -> minfi -> HDF5Array, whose .onLoad creates global-counter files in
  # tempdir() and FAILS if they already exist. Rserve forks share the daemon's tempdir,
  # so the 2nd+ request finds the previous fork's counters and DMRcate "fails to load".
  # Remove the stale counters so .onLoad can recreate them; surface the real error if any.
  unlink(list.files(tempdir(), pattern = "HDF5Array.*global_counter", full.names = TRUE), force = TRUE)
  dmr_load <- tryCatch({ suppressMessages({ library(DMRcate); library(GenomicRanges) }); TRUE },
                       error = function(e) conditionMessage(e))
  if(!isTRUE(dmr_load)) return(paste0("RES-NO; DMRcate could not be loaded: ", dmr_load))

  # Per-CpG result saved by the feature-level run. Fixed column layout:
  # 1=Feature 2=Coefficient 3=Average level 4=T_statistic 5=P_value 6=Adjusted p_value 7=chr 8=pos_hg38
  d    <- data.table::fread("dea_results.csv", data.table = FALSE)
  feat <- d[[1]]; coef <- d[[2]]; stat <- d[[4]]; rawp <- d[[5]]; adjp <- d[[6]]; chrom <- d[[7]]
  pos  <- suppressWarnings(as.integer(d[[8]]))
  gene <- if(ncol(d) >= 10) as.character(d[[10]]) else rep("", nrow(d))   # per-CpG gene (col 10)
  ok   <- is.finite(stat) & is.finite(adjp) & !is.na(chrom) & !is.na(pos)
  if(sum(ok) < 10) return("RES-NO; too few annotated CpGs for region analysis")

  gr <- GRanges(seqnames = chrom[ok],
                ranges   = IRanges(start = pos[ok], width = 1),
                stat     = stat[ok],
                diff     = coef[ok],                 # M-scale effect (the DEA "Coefficient")
                rawpval  = rawp[ok],
                ind.fdr  = adjp[ok],
                is.sig   = adjp[ok] < fdr.thresh,
                gene     = gene[ok])
  names(gr) <- feat[ok]
  gr <- gr[order(as.character(GenomicRanges::seqnames(gr)), GenomicRanges::start(gr))]

  # The CpGannotated object requires a betas slot, but DMRcate only uses it for plotting
  # and the betacutoff filter (neither of which we use), so a placeholder avoids loading
  # the 252MB beta matrix. Do NOT pass betacutoff -- it would touch the (placeholder) betas.
  betas <- matrix(0, nrow = length(gr), ncol = 1); rownames(betas) <- names(gr)
  annot <- new("CpGannotated", ranges = gr, betas = betas)

  dmrc <- tryCatch(dmrcate(annot, lambda = lambda, C = C), error = function(e) NULL)
  if(is.null(dmrc) || length(dmrc@coord) == 0) return("RES-NO; no differentially methylated regions found")

  # Build the region table from the DMResults slots directly. We deliberately do NOT call
  # extractRanges(): it fetches gene annotation from the DMRcatedata package and, when that
  # is absent, drops into an interactive "Install DMRcatedata [yes/no]" prompt that hangs
  # forever under Rserve. Instead we parse region coordinates from dmrc@coord and annotate
  # genes from the per-CpG `gene` column we already read.
  coord   <- dmrc@coord
  chr_r   <- sub(":.*", "", coord)
  se      <- sub(".*:", "", coord)
  start_r <- suppressWarnings(as.integer(sub("-.*", "", se)))
  end_r   <- suppressWarnings(as.integer(sub(".*-", "", se)))

  # genes per region = union of the array-annotated genes of the CpGs inside it
  reg     <- GRanges(chr_r, IRanges(start_r, end_r))
  hits    <- GenomicRanges::findOverlaps(reg, gr)
  gmap    <- tapply(gr$gene[S4Vectors::subjectHits(hits)], S4Vectors::queryHits(hits),
                    function(g){ g <- unlist(strsplit(g, "[;,]")); g <- unique(trimws(g))
                                 g <- g[!is.na(g) & nzchar(g)]; paste(sort(g), collapse = ", ") })
  genes_r <- rep("", length(reg)); if(length(gmap)) genes_r[as.integer(names(gmap))] <- as.character(gmap)

  # exact CpG ids covered by each region (for the on-screen tooltip + download)
  cmap    <- tapply(names(gr)[S4Vectors::subjectHits(hits)], S4Vectors::queryHits(hits),
                    function(c) paste(c, collapse = ", "))
  cpgs_r  <- rep("", length(reg)); if(length(cmap)) cpgs_r[as.integer(names(cmap))] <- as.character(cmap)

  out <- data.frame(
    Feature      = coord,
    chr          = chr_r,
    start        = start_r,
    end          = end_r,
    width        = end_r - start_r + 1,
    n_cpgs       = dmrc@no.cpgs,
    mean_diff_M  = signif(dmrc@meandiff, 3),
    max_diff_M   = signif(dmrc@maxdiff, 3),
    Stouffer_FDR = signif(dmrc@Stouffer, 3),
    HMFDR        = signif(dmrc@HMFDR, 3),
    genes        = genes_r,
    cpgs         = cpgs_r,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$Stouffer_FDR), ]
  data.table::fwrite(out, "dmr_results.csv")
  return(paste0("RES-OK;", nrow(out)))
}


################################################################################
## Methylation PATHWAY analysis (missMethyl gsameth = ORA; methylGSA methylRRA =
## GSEA). Dispatched from performORA/performGSEA when dea_results.csv is CpG-level
## (has gene + pos_hg38, no Gene_ID). Reuses the project's Entrez gene-set
## libraries; maps CpGs->genes WITH CpG-count bias correction (the methylation-
## specific step). Writes the SAME functional_results.csv + ridgeline.csv as the
## gene-level path, so the frontend table/ridgeline/dialog are unchanged.
## CpG->Entrez comes from the precomputed libraries/methylation_cpg_anno.qs
## (built offline by methylation/build_methylation_cpg_anno.R) -> no manifest load.
################################################################################

# load a gene-based Entrez library (.rds with $sets + $term)
.loadMethylLibrary <- function(funcLib){
  lib_file <- paste0(other.tables.path, "libraries/", funcLib, ".rds")
  if(!file.exists(lib_file)) return(NULL)
  rds <- readRDS(lib_file)
  list(sets = lapply(rds$sets, as.character),
       set.names = data.frame(names = rds$term, IDs = names(rds$sets),
                              stringsAsFactors = FALSE))
}

# per-gene(Entrez) statistic for the ridgeline = mean T-stat across the gene's CpGs
.methylGeneStat <- function(degs, map){
  m <- merge(data.frame(cpg = degs$Feature, stat = as.numeric(degs$T_statistic),
                        stringsAsFactors = FALSE),
             map[, c("cpg", "entrez")], by = "cpg")
  m <- m[!is.na(m$stat) & !is.na(m$entrez), ]
  agg <- tapply(m$stat, m$entrez, mean)
  data.frame(GeneID = as.character(names(agg)), GeneStat = as.numeric(agg),
             stringsAsFactors = FALSE)
}

# write ridgeline.csv (same columns as the gene-level path) from a per-gene stat
.writeMethylRidgeline <- function(geneStat, libraryList, resOut, fdr){
  geneStat$GeneID <- as.character(geneStat$GeneID)
  genesPW <- reshape::melt(libraryList, level = 1)
  colnames(genesPW) <- c("GeneID", "Set ID")
  genesPW$GeneID <- as.character(genesPW$GeneID)
  rr <- data.frame(`Set ID`      = resOut$`Set ID`,
                   `Set Name`    = resOut$`Set Name`,
                   `Set Size`    = resOut$`Set Size`,
                   `Adj P_value` = resOut$`Adj P_value`,
                   pathNegLogP   = -log10(resOut$P_value),
                   RankSig       = rank(resOut$P_value),
                   check.names = FALSE, stringsAsFactors = FALSE)
  genesPW <- merge(genesPW, rr, by = "Set ID", all.y = TRUE, all.x = FALSE)
  genesPW$Sig <- as.character(genesPW$`Adj P_value` < fdr)
  genesPW <- merge(genesPW, geneStat, by = "GeneID", all.x = TRUE, all.y = FALSE)
  genesPW <- na.omit(genesPW)
  genesPW <- data.table::as.data.table(genesPW)
  genesPW <- genesPW[, .(GeneID, `Set Name`, `Set Size`, pathFDR = `Adj P_value`,
                         Sig, pathNegLogP, GeneStat,
                         MeanGeneStat = mean(GeneStat), RankSig), by = .(`Set ID`)]
  genesPW <- genesPW[order(MeanGeneStat)]
  write.csv(genesPW, "ridgeline.csv", row.names = FALSE)
}

## --- GSEA (RRA-aggregated preranked GSEA via fgsea; no clusterProfiler) -----
## Same statistics as methylGSA methylRRA(method="GSEA"): each gene's Robust Rank
## Aggregation score (RobustRankAggreg::rhoScores over its CpG p-values) becomes a
## gene-level GSEA statistic, run through the tool's existing fgsea on the Entrez
## libraries directly. Advantages over calling methylRRA: (a) keeps MULTIPLE genes
## per CpG (a CpG feeds every gene it is annotated to), and (b) avoids
## clusterProfiler's heavy stack -> sub-minute + a few hundred MB (not ~3 GB).
.runMethylationGSEA <- function(funcLib, fdr, collapse = "true", mode = "tool"){
  suppressMessages({ library(fgsea); library(RobustRankAggreg); library(data.table) })
  set.seed(42)
  fdr <- as.numeric(fdr)
  lib <- .loadMethylLibrary(funcLib)
  if(is.null(lib)) return(paste0("RES-NO; Library file not found for methylation GSEA: ", funcLib, ".rds"))
  libraryList <- lib$sets

  degs <- data.table::fread("dea_results.csv", select = c("Feature", "P_value", "T_statistic"))
  cpg.pval <- setNames(as.numeric(degs$P_value), degs$Feature)
  cpg.pval <- cpg.pval[!is.na(cpg.pval) & cpg.pval > 0]
  if(length(cpg.pval) < 100) return("RES-NO; too few CpG p-values for methylation GSEA")

  # CpG -> Entrez, ALL annotated genes: a CpG shared by N genes feeds each of them.
  map <- qs::qread(paste0(other.tables.path, "libraries/methylation_cpg_anno.qs"))
  mp  <- map[map$cpg %in% names(cpg.pval), c("cpg", "entrez")]
  mp$pval <- cpg.pval[mp$cpg]

  # per-gene Robust Rank Aggregation score -> gene-level GSEA statistic (as methylRRA)
  pv.by.gene <- split(mp$pval, as.character(mp$entrez))
  rho <- vapply(pv.by.gene, function(p) RobustRankAggreg::rhoScores(p), numeric(1))
  z <- qnorm(rho / 2, lower.tail = FALSE)
  if(any(is.infinite(z))) z[is.infinite(z)] <- max(z[is.finite(z)], na.rm = TRUE)
  z <- sort(z, decreasing = TRUE)
  if(length(z) < 10) return("RES-NO; too few genes after CpG->gene aggregation")

  # preranked GSEA via fgsea (scoreType='pos': the RRA statistic is non-negative)
  res <- tryCatch(
    fgsea::fgsea(pathways = libraryList, stats = z, minSize = 15, maxSize = 500, scoreType = "pos"),
    error = function(e) paste0("ERROR:", conditionMessage(e)))
  if(is.character(res) && grepl("^ERROR:", res)) return(paste0("RES-NO; fgsea (methylation GSEA) failed: ", sub("^ERROR:", "", res)))
  res <- as.data.frame(res)
  if(nrow(res) == 0) return("RES-NO; methylation GSEA returned no gene sets")

  out <- data.frame(
    `Set Name`      = lib$set.names$names[match(as.character(res$pathway), lib$set.names$IDs)],
    P_value         = signif(res$pval, 3),
    `Adj P_value`   = signif(res$padj, 3),
    `Normalized ES` = round(res$NES, 3),
    `Set Size`      = res$size,
    `Set ID`        = res$pathway,
    check.names = FALSE, stringsAsFactors = FALSE)
  out$`Set Name`[is.na(out$`Set Name`)] <- as.character(out$`Set ID`)[is.na(out$`Set Name`)]
  out <- out[order(out$P_value), ]
  write.csv(out, "functional_results.csv", row.names = FALSE)
  num.sig <- sum(out$`Adj P_value` < fdr, na.rm = TRUE)

  geneStat <- tryCatch(.methylGeneStat(degs, map), error = function(e) NULL)
  if(!is.null(geneStat)) tryCatch(.writeMethylRidgeline(geneStat, libraryList, out, fdr),
                                  error = function(e) print(paste("methyl GSEA ridgeline failed:", e$message)))
  return(paste0("RES-OK;", num.sig))
}

## --- ORA (Wallenius bias-corrected over-representation; no missMethyl/HDF5Array)
## Reproduces missMethyl::gsameth's bias-corrected ORA, computed directly from our
## precomputed CpG->Entrez map so we never load missMethyl -> minfi -> HDF5Array
## (HDF5Array's .onLoad collides across forked Rserve requests). Pieces, identical
## to gsameth: per-gene fractional CpG count (bias) + DE indicator -> probability
## weighting function (limma::tricubeMovingAverage == gsameth's .estimatePWF) ->
## Wallenius noncentral hypergeometric test (BiasedUrn). All genes per CpG are kept,
## with fractional counts for CpGs shared by multiple genes.
.runMethylationORA <- function(funcLib, fdr, mode = "tool"){
  suppressMessages({ library(data.table) })   # limma + BiasedUrn via :: (no minfi/HDF5Array)
  fdr <- as.numeric(fdr)
  lib <- .loadMethylLibrary(funcLib)
  if(is.null(lib)) return(paste0("RES-NO; Library file not found for methylation ORA: ", funcLib, ".rds"))
  libraryList <- lib$sets

  degs <- data.table::fread("dea_results.csv", select = c("Feature", "sig"))
  sig.cpg <- unique(degs$Feature[degs$sig != "NS"])
  all.cpg <- unique(degs$Feature)
  if(length(sig.cpg) < 1)
    return("RES-NO; no significant CpGs (sig != NS) for ORA. For continuous phenotypes with high FDR, use GSEA instead.")

  # flat CpG->Entrez annotation (ALL annotated genes per CpG), restricted to tested CpGs
  map  <- qs::qread(paste0(other.tables.path, "libraries/methylation_cpg_anno.qs"))
  flat <- map[map$cpg %in% all.cpg, c("cpg", "entrez")]
  flat$entrez <- as.character(flat$entrez)
  if(nrow(flat) == 0) return("RES-NO; no CpG->gene mappings available for ORA")

  # multimap = # genes per CpG; inv = 1/multimap (fractional weight for shared CpGs)
  mmc <- table(flat$cpg)
  flat$inv <- 1 / as.numeric(mmc[flat$cpg])

  # universe genes + fractional CpG count per gene (bias) + DE indicator + fractional DE weight
  eg.universe   <- sort(unique(flat$entrez))
  equivN        <- as.numeric(tapply(flat$inv, flat$entrez, sum)[eg.universe])   # bias = frac #CpGs/gene
  sigflat       <- flat[flat$cpg %in% sig.cpg, ]
  eg.sig        <- unique(sigflat$entrez)
  test.de       <- as.integer(eg.universe %in% eg.sig)                            # 1 if gene has a sig CpG
  sorted.eg.sig <- eg.universe[test.de == 1]
  if(length(sorted.eg.sig) < 1) return("RES-NO; significant CpGs map to no annotated genes")
  fracw         <- pmin(tapply(sigflat$inv, sigflat$entrez, sum), 1)
  frac          <- as.numeric(fracw[sorted.eg.sig])

  # probability weighting function = P(DE) smoothed along the bias (== gsameth's .estimatePWF)
  o   <- order(equivN)
  pwf <- numeric(length(equivN)); pwf[o] <- limma::tricubeMovingAverage(test.de[o], span = 0.5)

  Nuniverse <- length(eg.universe); m <- length(sorted.eg.sig)
  sets <- lapply(libraryList, function(x) intersect(as.character(x), eg.universe))
  sets <- sets[vapply(sets, length, 0L) > 0L]
  if(length(sets) == 0) return("RES-NO; no gene sets overlap the methylation universe")

  Nset <- vapply(sets, length, 0L)
  DE   <- vapply(sets, function(s) sum(frac[sorted.eg.sig %in% s]), 0)
  P.DE <- vapply(seq_along(sets), function(i){
    inset <- eg.universe %in% sets[[i]]
    odds  <- (sum(pwf[inset]) / Nset[i]) / (sum(pwf[!inset]) / (Nuniverse - Nset[i]))
    p <- BiasedUrn::pWNCHypergeo(DE[i], Nset[i], Nuniverse - Nset[i], m, odds, lower.tail = FALSE) +
         BiasedUrn::dWNCHypergeo(DE[i], Nset[i], Nuniverse - Nset[i], m, odds)
    if(is.na(p) || p == 0)
      p <- stats::phyper(q = DE[i] - 0.5, m = m, n = Nuniverse - m, k = Nset[i], lower.tail = FALSE)
    p
  }, 0)
  FDR <- stats::p.adjust(P.DE, method = "BH")

  ids <- names(sets)
  out <- data.frame(
    `Set Name`    = lib$set.names$names[match(ids, lib$set.names$IDs)],
    P_value       = signif(P.DE, 3),
    `Adj P_value` = signif(FDR, 3),
    Hits          = floor(DE),
    `Set Size`    = Nset,
    `Set ID`      = ids,
    check.names = FALSE, stringsAsFactors = FALSE)
  out$`Set Name`[is.na(out$`Set Name`)] <- as.character(out$`Set ID`)[is.na(out$`Set Name`)]
  out <- out[order(out$P_value), ]
  write.csv(out, "functional_results.csv", row.names = FALSE)
  num.sig <- sum(out$`Adj P_value` < fdr, na.rm = TRUE)

  # ridgeline: per-gene mean T across the gene's CpGs (uses the full multi-gene map)
  degs2 <- data.table::fread("dea_results.csv", select = c("Feature", "T_statistic"))
  geneStat <- tryCatch(.methylGeneStat(degs2, map), error = function(e) NULL)
  if(!is.null(geneStat)) tryCatch(.writeMethylRidgeline(geneStat, libraryList, out, fdr),
                                  error = function(e) print(paste("methyl ORA ridgeline failed:", e$message)))
  return(paste0("RES-OK;", num.sig))
}


## Explain WHY a complete-case filter left too few donors.
## DonorRegression drops any donor missing ANY selected variable, so one sparse
## control variable can wipe out the analysis. It used to return a bare "RES-NO",
## leaving the UI to say only "not enough donors" -- this names the variable(s)
## responsible and how many donors each one costs.
## Reporting only: called on the failure path, it never changes which donors are
## analysed, and the string still starts with "RES-NO" so existing callers are
## unaffected.
##   meta.cand  donors that have omics data, BEFORE na.omit; record_id + one
##              column per selected variable
##   vars       selected variables (analysisVar first, then the covariates)
##   meta.info  proc_variable_summary, used for the human-readable labels
##   headline   sentence stating which requirement was not met
##   min.n      donors the failed check required; a "remove this one" hint is only
##              offered when dropping that variable actually reaches it
.explainMissingVars <- function(meta.cand, vars, analysisVar, meta.info, headline, min.n){
  detail <- tryCatch({
    vars <- intersect(vars, colnames(meta.cand))
    label.of <- function(v){
      lab <- meta.info$display[meta.info$column == v]
      lab <- lab[!is.na(lab) & nzchar(lab)]
      if(length(lab) > 0) lab[1]
      else if(identical(v, "batch")) "Batch (added automatically for this dataset)"
      else v
    }
    n.missing <- vapply(vars, function(v) sum(is.na(meta.cand[[v]])), integer(1))
    culprits  <- vars[n.missing > 0]
    if(length(culprits) == 0){
      ""
    } else {
      culprits <- culprits[order(n.missing[culprits], decreasing = TRUE)]
      lst <- paste0(vapply(culprits, function(v){
               paste0("'", label.of(v), "' is missing for ", n.missing[[v]], " of ",
                      nrow(meta.cand), " donors",
                      if(identical(v, analysisVar)) " (this is your variable of interest)" else "")
             }, character(1)), collapse = "; ")
      # Of the control variables, which single one gives back the most donors?
      # Only suggest it when dropping it actually clears the threshold -- a hint
      # that still leaves the analysis impossible would just mislead.
      covs <- setdiff(culprits, analysisVar)
      gain <- ""
      if(length(covs) > 0){
        n.without <- vapply(covs, function(v){
          keep <- setdiff(vars, v)
          if(length(keep) == 0) nrow(meta.cand)
          else sum(stats::complete.cases(meta.cand[, keep, drop = FALSE]))
        }, integer(1))
        best <- covs[which.max(n.without)]
        if(n.without[[best]] >= min.n){
          gain <- paste0(" Removing the control variable '", label.of(best), "' would leave ",
                         n.without[[best]], " donors.")
        }
      }
      paste0(" ", lst, ".", gain)
    }
  }, error = function(e) "")
  paste0("RES-NO; ", headline, detail)
}

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
    batch="b1",
    peakmode="positive",
    metaboGluc="LG",
    LOD_corrected = "true",
    tissueType = "Adipose",
    contamClass = "Dioxin",
     version="v1",
    mode = "local"){ 
print(c("version",version))
print(c( LOD_corrected,tissueType,contamClass))
  #  print(c(metaboGluc,"metaboGluc"))
#  print(c(varGroup,analysisVar,omicsType,batch,peakmode,metaboGluc,mode))
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
    # Accept either "," or ";" as the covariate separator. Historically this
    # was comma-only, but the Java Regression endpoint's isValidList
    # validator only allows semicolons in the whitelist (it rejects commas
    # as illegal characters). The frontend now sends ";" to satisfy the
    # validator — splitting on both keeps backward compatibility with any
    # caller still sending commas.
    fixedEffects <- strsplit(fixedEffects, "[,;]")[[1]]
  }
  
  # Region-level (DMR) analysis (omicsType "proc_methylation_dmr"): the per-CpG limma
  # result was already computed + saved by the feature-level run, so do NOT recompute.
  # Read dea_results.csv and cluster it into regions. Early return -> never loads the
  # M matrix or re-fits limma. No-op for any other omicsType.
  if(identical(omicsType, "proc_methylation_dmr")){
    return(.runMethylationDMRfromSaved(pvalThresh))
  }

if(mode == "tool"){
  # get data
if(omicsType == "proc_methylation"){
  # Methylation lives in its own database; everything else is handled identically.
  mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_methylation.sqlite"))
 }else if(version=="v2"){
  mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_omics_v2.sqlite"))

 }else{
  mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_omics.sqlite"))
}
 
    if(omicsType == "proc_pbrna"){
        table.nm <- paste0(omicsType, "_", cell)
        feature_table <- dbReadTable(mydb, table.nm)
        if(version=="v2"){
          # v2: use symbol as row ID, 4 info columns (symbol, gene_id, ensembl, name)
          feature_table <- feature_table[!is.na(feature_table$symbol) & feature_table$symbol != "", ]
          feature_table <- feature_table[!duplicated(feature_table$symbol), ]
          rownames(feature_table) <- feature_table$symbol
          feature_info <- feature_table[,c(1:4)]
          feature_table <- feature_table[,-c(1:4)]
        } else {
          rownames(feature_table) <- feature_table$gene_id
          feature_info <- feature_table[,c(1:3)]
          feature_table <- feature_table[,-c(1:3)]
        }
      }else if(omicsType == "proc_metabolite"){
        # batch param: "combat" -> ComBat-corrected table (no batch covariate);
        # anything else -> raw table, batch forced as covariate.
        # DB stores pre-split tables: proc_metabolite_{HG|LG|ratio}
        #                         and proc_metabolite_combat_{HG|LG|ratio}
        use_combat <- (batch == "combat")
        base_nm <- if(use_combat) "proc_metabolite_combat" else "proc_metabolite"
        # Map frontend metaboGluc value to the table suffix used in the database
        gluc_suffix <- switch(metaboGluc,
          "HG"          = "HG",
          "LG"          = "LG",
          "HG_LG_ratio" = "ratio",
          "LG"   # default fallback
        )
        table.nm <- paste0(base_nm, "_", gluc_suffix)
        feature_table <- dbReadTable(mydb, table.nm)

        # De-duplicate compound rows so they can be used as rownames
        feature_table <- feature_table[!is.na(feature_table$compound) & feature_table$compound != "", ]
        feature_table <- feature_table[!duplicated(feature_table$compound), ]
        rownames(feature_table) <- feature_table$compound
        feature_info <- feature_table[, 1:9]
        feature_table <- feature_table[, -(1:9)]
        # Column names are already clean record_ids — no _HG/_LG suffix stripping needed

        # Read batch info only when batch will be used as a covariate
        if(!use_combat){
          batch_info <- read.csv(paste0(other.tables.path, "display_interface/batch_info.csv"))
          # For ratio, one batch per record; pick the LG sample row as representative
          gluc_tag <- if(metaboGluc == "HG_LG_ratio") "LG" else metaboGluc
          batch_info <- batch_info[grepl(gluc_tag, batch_info$sample), ]
        }

      }else if(omicsType == "proc_metabo_optilcms"| omicsType == "proc_metabo_asari"){
        table.nm <- paste0(omicsType, "_", batch)

        feature_table <- dbReadTable(mydb, table.nm)
        if(peakmode!="mix"){
           feature_table =  feature_table[feature_table$mode==peakmode,]
        }

        rownames(feature_table) <- feature_table$peak
        feature_info <- feature_table[,c(1:7)]
       feature_table <- feature_table[,-c(1:7)]
       feature_table <- feature_table[,grepl(metaboGluc,colnames(feature_table))]
       colnames(feature_table) <- gsub(paste0("_",metaboGluc),"",   colnames(feature_table) )

      }else if(omicsType == "proc_flux"){
        # Handle flux simulation data with glucose treatment
        if(metaboGluc == "HG_LG_ratio"){
          table.nm <- paste0(omicsType, "_ratio")
        } else {
          table.nm <- paste0(omicsType, "_", metaboGluc)
        }
        feature_table <- dbReadTable(mydb, table.nm)
        rownames(feature_table) <- feature_table$rxn
        feature_info <- feature_table[,c(1,3,2,4)]
        feature_table <- feature_table[,-c(1:4)]

      }else if(omicsType == "proc_contaminants"){
       feature_table <- dbReadTable(mydb, omicsType)
       feature_table <- feature_table[feature_table$Tissue ==tissueType & feature_table$LOG_corrected==LOD_corrected,]
      if(contamClass !="all"){
     feature_table <- feature_table[feature_table$Class==contamClass,]
    }
       feature_info <- feature_table[,1:4]
        rownames(feature_table) <- feature_table$Compound

       feature_table <- feature_table[,-c(1:4)]
        
      }else if(omicsType == "proc_methylation"){
        # DEA runs on M-values. First 8 cols are annotation
        # (feature_id, chr, pos_hg38, strand, gene, gene_region, cgi_relation, cgi_name).
        # Fast-load .qs (~1.3s) instead of sqlite (~27s); built once by
        # scripts/convert_methylation_to_qs.R from HI_methylation.sqlite. Path is
        # derived from other.tables.path (same base used for cell_signal.rds below).
        feature_table <- qs::qread(paste0(other.tables.path, "omics_processing_input/proc/proc_methylation_M.qs"))
        rownames(feature_table) <- feature_table$feature_id
        feature_info <- feature_table[,c(1:8)]
        feature_table <- feature_table[,-c(1:8)]
      }else if(omicsType=="proc_rnaseq" & version=="v2"){
       feature_table <- dbReadTable(mydb, omicsType)
  rownames(feature_table) <- feature_table$accession
      feature_info <- feature_table[,c(1:5)]
  feature_table <- feature_table[,-c(1:5)]
      }else if(omicsType %in% c("proc_prot_b1","proc_prot_b2","proc_prot_v2")){
  # Proteomics v2: id = unique feature key (cleaned Protein_Group), symbol = display,
  # gene_id = pathway. proc_prot_v2 = combined batches: batch=="combat" -> the ComBat
  # table; anything else -> the merged (no-ComBat) table + batch forced as a covariate
  # (read below). proc_prot_b1/b2 load their own table directly.
  prot.tab <- if(omicsType == "proc_prot_v2"){ if(batch == "combat") "proc_prot_combat" else "proc_prot_combine" } else omicsType
  feature_table <- dbReadTable(mydb, prot.tab)
  rownames(feature_table) <- feature_table$id
  prot.info.cols <- c("id","gene_id","symbol","name","Protein_Group")
  feature_info  <- feature_table[, prot.info.cols]
  feature_table <- feature_table[, !(colnames(feature_table) %in% prot.info.cols)]
  if(omicsType == "proc_prot_v2" && batch != "combat"){
    batch_info <- read.csv(paste0(other.tables.path, "display_interface/batch_info_prot.csv"))
  }
      }else if(omicsType == "proc_nanostring_merge" & version=="v2"){
  # v2: KEY = symbol (convention: non-rnaseq genes calculate on symbol).
  # proc_nanostring_merge cols: symbol, gene_id, ensembl, name
  feature_table <- dbReadTable(mydb, omicsType)
  feature_table <- feature_table[!is.na(feature_table$symbol) & feature_table$symbol != "", ]
  feature_table <- feature_table[!duplicated(feature_table$symbol), ]
  rownames(feature_table) <- feature_table$symbol
  feature_info <- feature_table[,c(1:4)]
  feature_table <- feature_table[,-c(1:4)]
      }else{
  feature_table <- dbReadTable(mydb, omicsType)
  rownames(feature_table) <- feature_table$gene_id
  feature_info <- feature_table[,c(1:3)]
  feature_table <- feature_table[,-c(1:3)]

      }
  dbDisconnect(mydb)
 
  # get metadata
  mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
  if(varGroup == 'ephys_donor'){
    metadata <- dbReadTable(mydb, "ephys_donor")
    metadata <- metadata[metadata$cell_type == cell & metadata$glucose_mM == glucose, ]
    
    if(length(fixedEffects) > 0){
        other.metadata <- dbReadTable(mydb, "proc_metadata")
        metadata <- merge(metadata, other.metadata, by = "record_id")
    }
  }else if(varGroup=='prohormone'){
   metadata <- dbReadTable(mydb, "prohormone")

   # prohormone lives in its own table (no covariate columns). When covariates
   # are requested, merge them in from proc_metadata by record_id — but ONLY the
   # covariate columns, because proc_metadata also holds the prohormone columns
   # and a full merge would collide. Mirrors the ephys_donor branch above.
   if(length(fixedEffects) > 0){
       other.metadata <- dbReadTable(mydb, "proc_metadata")[, c("record_id", fixedEffects)]
       metadata <- merge(metadata, other.metadata, by = "record_id")
   }
   }else {
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
   # Methylation's feature table is huge (729k rows); default gzip saveRDS costs
   # ~20s. qsave is ~2s (and qread ~1.3s for the local "reproduce in R" path
   # below). Other omics keep .rds, so they are completely unaffected.
   if(omicsType == "proc_methylation"){
     qs::qsave(feature_table, "savedAnalysis/features.qs", preset = "high")
   } else {
     saveRDS(feature_table, "savedAnalysis/features.rds")
   }
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
    feature_table <- if(omicsType == "proc_methylation") qs::qread("features.qs") else readRDS("features.rds")
    metadata <- readRDS("meta.rds")
    meta.info <- readRDS("meta_info.rds")
    feature_info <- readRDS("feature_info.rds")
    nonendo <- readRDS("nonendo.rds")
    analysis.type <- meta.info$type[meta.info$column == analysisVar]
}
  # Fix: SQLite stores all columns as TEXT; convert continuous columns to numeric
  cont_cols <- meta.info$column[meta.info$type == "cont"]
  for(col in intersect(cont_cols, colnames(metadata))) {
    if(is.character(metadata[[col]])) {
      metadata[[col]][metadata[[col]] == "NA" | metadata[[col]] == ""] <- NA
      metadata[[col]] <- as.numeric(metadata[[col]])
    }
  }

  # Guard: non-finite values (e.g. Inf from a ratio with a zero denominator) -> NA,
  # so na.omit drops that donor instead of lmFit failing on "NA/NaN/Inf".
  for(col in intersect(cont_cols, colnames(metadata))) {
    if(is.numeric(metadata[[col]])) metadata[[col]][!is.finite(metadata[[col]])] <- NA
  }

  all.vars <- c(analysisVar, fixedEffects)

  # For proc_metabolite (raw, non-combat), force batch as a covariate for LG / HG / ratio.
  # For proc_metabolite_combat (batch == "combat") the batch effect is already removed,
  # so we never add it as a covariate.
  if((omicsType == "proc_metabolite" || omicsType == "proc_prot_v2") && batch != "combat"){
    if(exists("batch_info")){
      batch_map <- unique(batch_info[, c("record_id", "batch")])
      metadata <- merge(metadata, batch_map, by = "record_id", all.x = FALSE)
      if(!"batch" %in% all.vars){
        all.vars <- c(all.vars, "batch")
      }
    }
  }

  # filter to keep only data with complete metadata & omics data
  metadata <- metadata[metadata$record_id %in% colnames(feature_table), c("record_id", all.vars)]
  # Keep the table as it stands BEFORE the complete-case filter: it is the only
  # thing that still knows which variable each dropped donor was missing, which is
  # what the "insufficient data" message reports back. Never used to pick donors --
  # it is kept in step with `metadata` through the filters below purely so the
  # message counts match the analysis that was attempted.
  meta.cand <- metadata
  metadata <- na.omit(metadata)
  if(dim(metadata)[1] < 5){
    return(.explainMissingVars(meta.cand, all.vars, analysisVar, meta.info,
      paste0("Only ", dim(metadata)[1], " of ", nrow(meta.cand),
             " donors with this omics data have a value for every selected variable (at least 5 are needed)."), 5))
  }
  rownames(metadata) <- metadata$record_id
  feature_table <- feature_table[,colnames(feature_table) %in% metadata$record_id]

  # make sure metadata and omics data donors are in the same order
  feature_table <- feature_table[,match(metadata$record_id, colnames(feature_table))]

  # eliminate any feature with fewer than 10 observations after complete case analysis with metadata variables
 if(omicsType != "proc_metabolite"){
  # Vectorized: identical result to the per-row apply() but ~100x faster on the
  # 729k-row methylation matrix (apply() over rows was the dominant cost).
  feature.keep <- rowSums(!is.na(feature_table)) >= 9
  feature_table <- feature_table[feature.keep, ]
  feature_info <- feature_info[feature.keep, ]
 }
  # filter by donors
  if(donors == "subset"){
    donor.list <- readRDS("donors.rds")
    feature_table <- feature_table[,colnames(feature_table) %in% donor.list]
    metadata <- metadata[metadata$record_id %in% donor.list, ]
    meta.cand <- meta.cand[meta.cand$record_id %in% donor.list, ]
  }

  # for pairwise discrete comparisons, keep only the two selected groups
  if(analysis.type == "disc" && contrast != "anova"){
    metadata <- metadata[metadata[,analysisVar] %in% c(ref, contrast), ]
    feature_table <- feature_table[,colnames(feature_table) %in% metadata$record_id]
    meta.cand <- meta.cand[meta.cand[,analysisVar] %in% c(ref, contrast), ]
  }

  # check if there are enough samples
  n.samps <- dim(feature_table)[2]
  if(analysis.type == "disc" && contrast != "anova"){
    if(n.samps < 10){
      return(.explainMissingVars(meta.cand, all.vars, analysisVar, meta.info,
        paste0("Only ", n.samps, " of ", nrow(meta.cand),
               " donors in the two selected groups have a value for every selected variable (at least 10 are needed for a two-group comparison)."), 10))
    }

    if(sum(metadata[,analysisVar] == ref) < 5){
      return(.explainMissingVars(meta.cand[meta.cand[,analysisVar] == ref, , drop = FALSE],
        all.vars, analysisVar, meta.info,
        paste0("Only ", sum(metadata[,analysisVar] == ref), " of ",
               sum(meta.cand[,analysisVar] == ref), " donors in group '", ref,
               "' have a value for every selected variable (at least 5 per group are needed)."), 5))
    }

    if(sum(metadata[,analysisVar] == contrast) < 5){
      return(.explainMissingVars(meta.cand[meta.cand[,analysisVar] == contrast, , drop = FALSE],
        all.vars, analysisVar, meta.info,
        paste0("Only ", sum(metadata[,analysisVar] == contrast), " of ",
               sum(meta.cand[,analysisVar] == contrast), " donors in group '", contrast,
               "' have a value for every selected variable (at least 5 per group are needed)."), 5))
    }

  } else {
    if(n.samps < 5){
      return(.explainMissingVars(meta.cand, all.vars, analysisVar, meta.info,
        paste0("Only ", n.samps, " of ", nrow(meta.cand),
               " donors have a value for every selected variable (at least 5 are needed)."), 5))
    }
  }
  
  # perform analysis
  # ---- Kendall rank correlation branch for environmental contaminants -----
  # Small sample sizes and non-normal abundance distributions make limma's
  # moderated-t inappropriate for proc_contaminants. We use Kendall's tau
  # (partial Kendall via ppcor when covariates are present) for continuous
  # primary metadata and 2-group discrete contrasts (0/1 encoding).
  # For >=3-group ANOVA we fall back to Kruskal-Wallis with N/A coefficients.
  # res.table is built with the same column layout topTable() returns, so the
  # downstream merge/rename/post-processing block applies unchanged.
  if(omicsType == "proc_contaminants"){
    feature_mat <- as.matrix(feature_table)  # rows = compounds, cols = donors

    if(analysis.type == "disc" && contrast == "anova" && length(unique(metadata[,analysisVar])) >= 3){
      # Kruskal-Wallis for >=3 groups; no per-feature coefficient available
      grp <- metadata[,analysisVar]
      kw <- apply(feature_mat, 1, function(x){
        ok <- !is.na(x)
        if(sum(ok) < 5 || length(unique(grp[ok])) < 2) return(c(NA, NA))
        tryCatch({
          res <- kruskal.test(x[ok] ~ as.factor(grp[ok]))
          c(NA, unname(res$p.value))
        }, error = function(e) c(NA, NA))
      })
      coef.vec <- kw[1, ]
      pval.vec <- kw[2, ]
      # Kruskal-Wallis has no per-feature effect-size, only a chi-squared p-value.
      # Use "Coefficient" header to match the contaminants Kendall path; downstream
      # will overwrite the column with "N/A" for ANOVA contrasts.
      coef.col.name <- "Coefficient"
    } else {
      # Encode 2-group discrete as 0/1 (ref=0, contrast=1); otherwise use as-is for continuous.
      if(analysis.type == "disc"){
        x.var <- ifelse(metadata[,analysisVar] == contrast, 1L,
                 ifelse(metadata[,analysisVar] == ref, 0L, NA_integer_))
      } else {
        x.var <- as.numeric(metadata[,analysisVar])
      }
      # Covariate matrix for partial Kendall (excludes analysisVar; lipid 'batch' covariate
      # never reaches contaminants — proc_contaminants is handled in its own branch above).
      cov.names <- setdiff(all.vars, analysisVar)
      cov.mat   <- if(length(cov.names) > 0) as.data.frame(metadata[, cov.names, drop = FALSE]) else NULL
      use.pcor  <- !is.null(cov.mat) && requireNamespace("ppcor", quietly = TRUE)

      tau.pval <- apply(feature_mat, 1, function(y){
        ok <- !is.na(y) & !is.na(x.var)
        if(!is.null(cov.mat)){
          ok <- ok & complete.cases(cov.mat)
        }
        if(sum(ok) < 5 || length(unique(x.var[ok])) < 2) return(c(NA, NA))
        tryCatch({
          if(use.pcor){
            r <- ppcor::pcor.test(y[ok], x.var[ok], cov.mat[ok, , drop = FALSE], method = "kendall")
            c(unname(r$estimate), unname(r$p.value))
          } else {
            r <- cor.test(y[ok], x.var[ok], method = "kendall", exact = FALSE)
            c(unname(r$estimate), unname(r$p.value))
          }
        }, error = function(e) c(NA, NA))
      })
      coef.vec <- tau.pval[1, ]
      pval.vec <- tau.pval[2, ]
      # Always use "Coefficient" (= Kendall's tau) — never "Log2FC", since tau is
      # not a log fold change. The downstream merge/post-processing path treats
      # column 1's name as opaque, so this is safe.
      coef.col.name <- "Coefficient"
    }

    avg.vec <- rowMeans(feature_mat, na.rm = TRUE)
    adj.vec <- p.adjust(pval.vec, method = "fdr")

    # Match topTable's column layout: <coef>, AveExpr, t, P.Value, adj.P.Val, B
    # (B is a placeholder that is dropped at the existing line `res.table <- res.table[,-7]`.)
    res.table <- data.frame(
      Coefficient = coef.vec,
      AveExpr     = avg.vec,
      t           = coef.vec,   # no t-statistic — reuse tau so downstream rename to T_statistic stays meaningful
      P.Value     = pval.vec,
      adj.P.Val   = adj.vec,
      B           = NA_real_,
      stringsAsFactors = FALSE
    )
    rownames(res.table) <- rownames(feature_mat)
    res.table <- res.table[!is.na(res.table$P.Value), ]
    res.table <- res.table[order(res.table$P.Value), ]
    colnames(res.table)[1] <- coef.col.name
    if(analysis.type == "disc" && contrast == "anova"){
      res.table[, 1] <- "N/A"
    }
    # Downstream column-drop step (line ~1696) references `length(myargs)`,
    # which is only defined inside the limma disc branch. Set a 2-element
    # placeholder so the ANOVA path correctly drops the B placeholder column.
    myargs <- list("kendall", "levels")

  } else if(analysis.type == "disc"){

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
    # Save limma objects for CAMERA (flux only)
    if(omicsType == "proc_flux"){
      qs::qsave(fit,             "camera_fit.qs")
      qs::qsave(contrast.matrix, "camera_contrast.qs")
      qs::qsave(feature_table,   "camera_features.qs")
      qs::qsave("disc",          "camera_type.qs")
    }
    fit <- contrasts.fit(fit, contrast.matrix)
    fit <- eBayes(fit)
    res.table <- topTable(fit, number = Inf)
    
    # Remove coefficients for ANOVA contrasts
    if(contrast == "anova"){
      if(length(myargs) > 2){
        res.table <- res.table[,-c(2:(length(myargs)-1))]
      }
      res.table[,1] <- "N/A"
      # Multi-group ANOVA topTable has NO 'B' column (F-test), so it is one column
      # short of the pairwise layout [<coef>,AveExpr,<stat>,P.Value,adj.P.Val,B].
      # That off-by-one misaligns the symbol/entrez/name columns in the annotation
      # block below (symbol -> Description, entrez -> Feature/ID). Add a B placeholder
      # and force myargs to length 2 so the layout + the downstream col-7 (B) drop
      # match the pairwise / KW paths exactly.
      # NOTE: unreachable from Omics View (which never sends contrast="anova"); this
      # path is used by the precompute pipeline, which sources THIS file.
      if(!"B" %in% colnames(res.table)) res.table$B <- NA_real_
      myargs <- list("kendall", "levels")
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
    # Save limma objects for CAMERA (flux only)
    if(omicsType == "proc_flux"){
      qs::qsave(fit,          "camera_fit.qs")
      qs::qsave(analysisVar,  "camera_coef.qs")
      qs::qsave(feature_table,"camera_features.qs")
      qs::qsave("cont",       "camera_type.qs")
    }
    fit <- eBayes(fit)
    res.table <- topTable(fit, number = Inf, coef = analysisVar)
    colnames(res.table)[1] <- "Coefficient"
    
  }

  # Remove results rows with NAs
  res.table <- res.table[!is.na(res.table$P.Value), ]  
  # Process results for output
 if(omicsType == "proc_metabolite"){
  res.table <- merge(
    res.table,
    feature_info[, c("compound","inchikey","hmdb_id","kegg_id","gem_id","super_class","main_class","sub_class")],
    by.x = "row.names", by.y = "compound"
  )
  for(col in c("inchikey","hmdb_id","kegg_id","gem_id","super_class","main_class","sub_class")){
    res.table[[col]][is.na(res.table[[col]])] <- ""
  }
  res.table$Description <- apply(res.table, 1, function(row){
    parts <- c()
    if(row["main_class"] != "") parts <- c(parts, row["main_class"])
    if(row["sub_class"] != "")  parts <- c(parts, paste0("(", row["sub_class"], ")"))
    paste(parts, collapse = " ")
  })
  }else if(omicsType=="proc_contaminants"){
  res.table <- merge(res.table, feature_info, by.x = "row.names", by.y = "Compound")
  }else if(omicsType=="proc_methylation"){
  # Attach the 7 remaining annotation cols (chr, pos_hg38, strand, gene,
  # gene_region, cgi_relation, cgi_name) by CpG id. No Gene_ID -> excluded
  # from the endocrine-signal merge below.
  res.table <- merge(res.table, feature_info, by.x = "row.names", by.y = "feature_id")
  }else if(omicsType=="proc_rnaseq" & version=="v2"){
  res.table <- merge(res.table, feature_info, by.x = "row.names", by.y = "accession")
   colnames(res.table)[8] <- "Gene_ID"
  }else if(omicsType=="proc_flux" & version=="v2"){
  res.table <- merge(res.table, feature_info, by.x = "row.names", by.y = "rxn")
  }else if(omicsType %in% c("proc_pbrna","proc_nanostring_merge") & version=="v2"){
  # v2 proc_pbrna: rownames are symbols, feature_info has (symbol, gene_id, ensembl, name)
  # Merge by symbol, only keep gene_id and name (skip ensembl) for standard 9-col format
  res.table <- merge(res.table, feature_info[,c("symbol","gene_id","name")], by.x = "row.names", by.y = "symbol")
  # After merge: col1=symbol, cols2-7=topTable, col8=gene_id(entrez), col9=name
  # Swap col1(symbol) and col8(gene_id) so Gene_ID=entrez for nonendo merge
  tmp_symbols <- res.table[,1]
  res.table[,1] <- res.table[,8]
  res.table[,8] <- tmp_symbols
  colnames(res.table)[1] <- "Gene_ID"
  colnames(res.table)[8] <- "symbol"
  }else if(omicsType %in% c("proc_prot_b1","proc_prot_b2","proc_prot_v2")){
  # Proteomics v2: key on id; expose Gene_ID (entrez -> pathway) + keep symbol (display).
  res.table <- merge(res.table, feature_info, by.x = "row.names", by.y = "id")
  colnames(res.table)[colnames(res.table) == "gene_id"] <- "Gene_ID"
  }else{
  res.table <- merge(res.table, feature_info, by.x = "row.names", by.y = "gene_id")
   colnames(res.table)[1] <- "Gene_ID"
  }
   
  if(contrast != "anova"){
    res.table <- res.table[,-7]
  } else if (contrast == "anova" & length(myargs) == 2){
    res.table <- res.table[,-7]
  }
 
  # add endocrine signal
 
  
   if(!omicsType %in% c("proc_metabolite","proc_contaminants","proc_flux","proc_methylation","proc_prot_b1","proc_prot_b2","proc_prot_v2")){
  res.table <- merge(res.table, nonendo, by = "Gene_ID", all.x = TRUE)
   if(omicsType=="proc_rnaseq" & version=="v2"){
     res.table <- res.table[,c(2,3:7,11,1,9,8)]
  colnames(res.table)[10] <- "Symbol" 
  }else{
  res.table <- res.table[,c(7,2:6,9,1,8)]

  }
  colnames(res.table)[9] <- "Description" 
   
 }

  colnames(res.table)[1] <- "Feature"
  colnames(res.table)[3] <- "Average level"
  colnames(res.table)[4] <- "T_statistic"
  colnames(res.table)[5] <- "P_value"
  colnames(res.table)[6] <- "Adjusted p_value"
  
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
 
  # Methylation's full table is ~90MB; fwrite is ~18x faster than write.csv and
  # produces a standard CSV. No methylation backend path re-reads this file
  # (only the user downloads it). Other omics keep write.csv -> byte-identical.
  if(omicsType == "proc_methylation"){
    data.table::fwrite(res.table, "dea_results.csv")
  } else {
    write.csv(res.table, "dea_results.csv", row.names = FALSE)
  }

  # Methylation returns ~729k probes -- too many for the browser to parse/render.
  # The on-screen table + volcano read the top 10000 by p-value (res.table is
  # already ordered by P_value above); the full result stays in dea_results.csv
  # for download. Other omics are unaffected (no second file written).
  if(omicsType == "proc_methylation"){
    write.csv(head(res.table, 10000), "dea_results_top.csv", row.names = FALSE)
  }

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
    version = "v2",
    mode="local"
){
  
  # load libraries
  library(dplyr)
  library(RSQLite)
  library(rhdf5)

  # process inputs
  pvalThresh <- as.numeric(pvalThresh)

 if(mode == "tool"){
    # get data (v2 uses hdf5_V2 directory; v1 uses hdf5)
  sc.h5.path <- ifelse(version=="v2", h5.v2.path, h5.path)
  table.path <- paste0(sc.h5.path, "sc_", cell, "_", glucose, ".h5")

  feature_table <- h5read(table.path, "data/norm_expression") %>% as.data.frame()
  cells <- h5read(table.path, "meta/cells/cellid")
  genes <- h5read(table.path, "meta/genes") %>% as.data.frame()
  H5close()

  if(version=="v2"){
    # v2: genes has columns: ensembl, entrez, name, symbol
    feature_info <- data.frame(gene_id = genes$entrez, symbol = genes$symbol, name = genes$name, stringsAsFactors = FALSE)
  } else {
    feature_info <- genes[,c(1,3,2)]
    colnames(feature_info) <- c("gene_id", "symbol", "name")
  }
  colnames(feature_table) <- cells

  if(version=="v2"){
    # v2: use symbol as primary identifier
    genes.keep <- !is.na(feature_info$symbol) & feature_info$symbol != ""
    feature_table <- feature_table[genes.keep, ]
    feature_info <- feature_info[genes.keep, ]
    dup.idx <- duplicated(feature_info$symbol)
    feature_table <- feature_table[!dup.idx, ]
    feature_info <- feature_info[!dup.idx, ]
    rownames(feature_table) <- feature_info$symbol
  } else {
    # v1: use gene_id (entrez) as primary identifier
    genes.keep <- !is.na(feature_info$gene_id)
    feature_table <- feature_table[genes.keep, ]
    feature_info <- feature_info[genes.keep, ]
    rownames(feature_table) <- feature_info$gene_id
  }

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
  if(version=="v2"){
    # v2: row.names are symbols; merge by symbol
    res.table <- merge(res, feature_info, by.x = "row.names", by.y = "symbol")
    # After merge: col1=Row.names(symbol), col2=r, col3=pval, col4=fdr, col5=avg.level, col6=gene_id(entrez), col7=name
    res.table <- res.table[,c(1,2,5,3:4,6,7)]
  } else {
    # v1: row.names are gene_ids (entrez); merge by gene_id
    res.table <- merge(res, feature_info, by.x = "row.names", by.y = "gene_id")
    # After merge: col1=Row.names(gene_id), col2=r, col3=pval, col4=fdr, col5=avg.level, col6=symbol, col7=name
    res.table <- res.table[,c(6,2,5,3:4,1,7)]
  }
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

####helpers
convert_gpr_to_r <- function(rule, present_genes) {
  if (is.na(rule) || rule == "" || rule == "NA") return(NA)
  
  r <- rule
  
  # Normalize case
  r <- gsub("AND", "and", r)
  r <- gsub("OR",  "or",  r)
  
  # Convert "and/or" to &, |
  r <- gsub("\\sand\\s", " & ", r)
  r <- gsub("\\sor\\s",  " | ", r)
  
  # Wrap gene IDs in membership test
  r <- gsub("(ENSG[0-9]+)", "\"\\1\" %in% present_genes", r)
  
  r
}
 eval_gpr <- function(expr, present_genes) {
   expr$map= sapply(expr$expr,function(x) eval(parse(text = x)))
   list(map=any(expr$map),expr=expr)
 }

## ============================================================================
## compareEndotypeStates -- Donor State Explorer "advanced comparison" web tool function.
## Uses the SAME default phenotype statistic as the in-page enrichment table:
##   numeric phenotype, NO covariate    -> Welch two-sample t-test (group A vs group B)
##   numeric phenotype, WITH covariate  -> lm(value ~ group + covariates); report group term
##   categorical phenotype              -> chi-square (Cramer's V effect size)
## BH-FDR across phenotypes. Inputs read by RELATIVE path via setPaths() (other.tables.path).
## DISTINCT file names (endotype_compare_*) + qs so it never clashes with the omics DEA
## (donors.rds / dea_results.csv). Result written to the user folder; returns RES-OK;sig;up;down;n.
##   groupA, groupB : semicolon-separated state numbers, e.g. "4" vs "1;2;3;5"
##   covariates     : semicolon-separated phenotype column names ("" = none -> Welch)
##   subsetKey      : "" = all donors; else reads donors_<subsetKey>.rds written by the app-filter phenotype filter
##   sourceFilter   : "all" | "omics" | "predicted"
## ============================================================================
compareEndotypeStates <- function(groupA, groupB, covariates = "", subsetKey = "",
                                  sourceFilter = "all", fdr = "true", pvalThresh = "0.05") {
  pvalThresh <- suppressWarnings(as.numeric(pvalThresh)); if(is.na(pvalThresh)) pvalThresh <- 0.05
  deliv <- paste0(other.tables.path, "donor_clustering_pipeline_protv2/endotype_deliverable_FINAL/")
  pheno_file <- paste0(other.tables.path, "processed_comprehensive/phenotype/all_features.csv")
  if(!file.exists(paste0(deliv, "donor_states.csv")) || !file.exists(pheno_file)){
    return("RES-NO; required input tables not found")
  }

  states  <- read.csv(paste0(deliv, "donor_states.csv"), stringsAsFactors = FALSE)
  featdef <- read.csv(paste0(deliv, "compare_features.csv"), stringsAsFactors = FALSE)
  pheno   <- read.csv(pheno_file, stringsAsFactors = FALSE, check.names = FALSE)
  states$donor_id <- as.character(states$donor_id)
  pheno$sample_id <- as.character(pheno$sample_id)
  meta <- merge(states, pheno, by.x = "donor_id", by.y = "sample_id")

  if(sourceFilter == "omics")         meta <- meta[meta$source == "omics", ]
  if(sourceFilter == "predicted")     meta <- meta[meta$source == "predicted", ]
  if(sourceFilter == "pred_high")     meta <- meta[meta$source == "omics" | (meta$source == "predicted" & meta$reliability == "high"), ]
  if(sourceFilter == "pred_high_med") meta <- meta[meta$source == "omics" | (meta$source == "predicted" & meta$reliability %in% c("high","medium")), ]
  # "all" / "pred_all" -> no filter (every donor)

  # donor subset from the reusable phenotype filter (app-filter writes donors_<subsetKey>.rds;
  # a DISTINCT name from the omics page's donors.rds, so the two never clash).
  if(nzchar(subsetKey)){
    sf <- paste0("donors_", subsetKey, ".rds")
    if(file.exists(sf)){
      keep <- as.character(readRDS(sf))
      meta <- meta[meta$donor_id %in% keep, ]
    }
  }

  gA <- suppressWarnings(as.integer(strsplit(groupA, ";")[[1]]))
  gB <- suppressWarnings(as.integer(strsplit(groupB, ";")[[1]]))
  grp <- ifelse(meta$state_num %in% gA, "A", ifelse(meta$state_num %in% gB, "B", NA))
  meta <- meta[!is.na(grp), , drop = FALSE]; grp <- grp[!is.na(grp)]
  if(sum(grp == "A") < 3 || sum(grp == "B") < 3) return("RES-NO; fewer than 3 donors per group")
  meta$.grp <- factor(grp, levels = c("B", "A"))   # B = reference; coefficient / diff = A - B

  covs <- if(nzchar(covariates)) strsplit(covariates, ";")[[1]] else character(0)
  covs <- intersect(covs, colnames(meta))
  use.lm <- length(covs) > 0
  cat.names <- featdef$name[featdef$type == "categorical"]

  # typed covariate frame (numeric covariate -> numeric; categorical covariate -> factor)
  covdf <- NULL
  if(use.lm){
    covdf <- data.frame(row.names = seq_len(nrow(meta)))
    for(cv in covs){
      num <- suppressWarnings(as.numeric(meta[[cv]]))
      covdf[[cv]] <- if(cv %in% cat.names || all(is.na(num))) factor(as.character(meta[[cv]])) else num
    }
  }

  rows <- list()
  for(i in seq_len(nrow(featdef))){
    nm <- featdef$name[i]; ty <- featdef$type[i]
    if(!(nm %in% colnames(meta)) || nm %in% covs) next

    if(ty == "categorical"){
      v <- as.character(meta[[nm]]); ok <- !is.na(v) & v != "" & v != "NA"
      if(sum(ok) < 10) next
      tab <- table(meta$.grp[ok], v[ok])
      if(nrow(tab) < 2 || ncol(tab) < 2 || any(rowSums(tab) == 0)) next
      ct <- suppressWarnings(tryCatch(chisq.test(tab), error = function(e) NULL))
      if(is.null(ct)) next
      nn <- sum(tab); cv <- sqrt(as.numeric(ct$statistic) / (nn * (min(dim(tab)) - 1)))
      rows[[length(rows) + 1]] <- data.frame(
        Feature = nm, Type = "categorical", Method = "chi-square",
        Effect = round(cv, 4), Mean_A = NA_real_, Mean_B = NA_real_,
        Statistic = round(as.numeric(ct$statistic), 3), P_value = as.numeric(ct$p.value),
        n_A = sum(meta$.grp[ok] == "A"), n_B = sum(meta$.grp[ok] == "B"), stringsAsFactors = FALSE)
    } else {
      y <- suppressWarnings(as.numeric(meta[[nm]]))
      ok <- !is.na(y)
      yA <- y[ok & meta$.grp == "A"]; yB <- y[ok & meta$.grp == "B"]
      mA <- if(length(yA)) mean(yA) else NA_real_; mB <- if(length(yB)) mean(yB) else NA_real_

      if(!use.lm){
        if(length(yA) < 3 || length(yB) < 3) next
        tt <- suppressWarnings(tryCatch(t.test(yA, yB), error = function(e) NULL))   # Welch (var.equal = FALSE)
        if(is.null(tt)) next
        sp <- sqrt(((length(yA) - 1) * var(yA) + (length(yB) - 1) * var(yB)) / (length(yA) + length(yB) - 2))
        d  <- if(is.finite(sp) && sp > 0) (mA - mB) / sp else NA_real_
        rows[[length(rows) + 1]] <- data.frame(
          Feature = nm, Type = "numeric", Method = "Welch t-test",
          Effect = round(d, 4), Mean_A = round(mA, 4), Mean_B = round(mB, 4),
          Statistic = round(as.numeric(tt$statistic), 3), P_value = as.numeric(tt$p.value),
          n_A = length(yA), n_B = length(yB), stringsAsFactors = FALSE)
      } else {
        df <- data.frame(.y = y, .grp = meta$.grp, stringsAsFactors = FALSE)
        df <- cbind(df, covdf)
        df <- df[complete.cases(df), , drop = FALSE]
        if(sum(df$.grp == "A") < 3 || sum(df$.grp == "B") < 3) next
        form <- as.formula(paste(".y ~ .grp +", paste(sprintf("`%s`", covs), collapse = " + ")))
        fit <- suppressWarnings(tryCatch(lm(form, data = df), error = function(e) NULL))
        if(is.null(fit)) next
        co <- summary(fit)$coefficients
        if(!(".grpA" %in% rownames(co))) next
        rows[[length(rows) + 1]] <- data.frame(
          Feature = nm, Type = "numeric", Method = "linear model (covariate-adjusted)",
          Effect = round(co[".grpA", "Estimate"], 4),
          Mean_A = round(mean(df$.y[df$.grp == "A"]), 4), Mean_B = round(mean(df$.y[df$.grp == "B"]), 4),
          Statistic = round(co[".grpA", "t value"], 3), P_value = as.numeric(co[".grpA", "Pr(>|t|)"]),
          n_A = sum(df$.grp == "A"), n_B = sum(df$.grp == "B"), stringsAsFactors = FALSE)
      }
    }
  }

  if(length(rows) == 0) return("RES-NO; no testable features for this comparison")
  res <- do.call(rbind, rows)
  res$Adjusted_p <- p.adjust(res$P_value, method = "fdr")
  thr <- if(fdr == "true") res$Adjusted_p else res$P_value
  res$sig <- "NS"
  res$sig[res$Type == "numeric" & thr < pvalThresh & res$Effect > 0] <- "up"
  res$sig[res$Type == "numeric" & thr < pvalThresh & res$Effect < 0] <- "down"
  res$sig[res$Type == "categorical" & thr < pvalThresh] <- "diff"
  sig.num  <- sum(thr < pvalThresh, na.rm = TRUE)
  sig.up   <- sum(res$sig == "up"); sig.down <- sum(res$sig == "down")
  res$P_value    <- signif(res$P_value, 3)
  res$Adjusted_p <- signif(res$Adjusted_p, 3)
  res <- res[order(res$P_value), ]
  res <- res[, c("Feature","Type","Method","Effect","Mean_A","Mean_B","Statistic","P_value","Adjusted_p","n_A","n_B","sig")]

  write.csv(res, "endotype_compare_result.csv", row.names = FALSE)
  paste0("RES-OK;", sig.num, ";", sig.up, ";", sig.down, ";", nrow(res))
}

## ============================================================================
## compareEndotypeOmics -- Donor State Explorer omics differential expression (limma).
## Reads the normalized omics matrix directly (processed_comprehensive/omics/), groups
## donors by endotype state (Group A vs Group B), runs limma with optional covariate
## adjustment, writes endotype_omics_de.csv to the user folder. DISTINCT names from the
## omics view's DonorRegression (dea_results.csv) -- no clash. Inputs read by relative path.
##   omicsType  : proc_rnaseq | proc_prot_v2 | proc_nanostring | proc_methylation |
##                proc_pbrna_alpha | proc_pbrna_beta | proc_metabolite | proc_flux
##   groupA/B   : ";"-joined state numbers ; covariates : ";"-joined phenotype names
##   subsetKey  : "" or donors_<subsetKey>.rds (app-filter) ; sourceFilter : page donor-set tier
## ============================================================================
## Map an endotype omicsType to its source SQLite DB + table. ALL omics are read from SQLite
## (the canonical source, matching the Omics page): proteomics -> ComBat table proc_prot_combat,
## metabolites -> ComBat tables proc_metabolite_combat_{LG|HG|ratio}, methylation -> the methylation DB.
.endotype_sqlite_tbl <- function(omicsType) switch(omicsType,
    proc_rnaseq           = list(db = "HI_omics_v2.sqlite",    tbl = "proc_rnaseq"),
    proc_prot_v2          = list(db = "HI_omics_v2.sqlite",    tbl = "proc_prot_combat"),
    proc_nanostring       = list(db = "HI_omics_v2.sqlite",    tbl = "proc_nanostring_merge"),
    proc_pbrna_alpha      = list(db = "HI_omics_v2.sqlite",    tbl = "proc_pbrna_Alpha"),
    proc_pbrna_beta       = list(db = "HI_omics_v2.sqlite",    tbl = "proc_pbrna_Beta"),
    proc_methylation      = list(db = "HI_methylation.sqlite", tbl = "proc_methylation_M"),
    proc_metabolite       = list(db = "HI_omics_v2.sqlite",    tbl = "proc_metabolite_combat_LG"),
    proc_metabolite_lg    = list(db = "HI_omics_v2.sqlite",    tbl = "proc_metabolite_combat_LG"),
    proc_metabolite_hg    = list(db = "HI_omics_v2.sqlite",    tbl = "proc_metabolite_combat_HG"),
    proc_metabolite_ratio = list(db = "HI_omics_v2.sqlite",    tbl = "proc_metabolite_combat_ratio"),
    NULL)

## Read the endotype omics feature matrix (info cols + donor cols) as a data.frame, OR a
## "RES-NO; ..." string on failure. Every omics is read from its SQLite table (proteomics +
## metabolites = ComBat-corrected); donor columns are record_ids that match donor_states.
.endotype_read_omics <- function(omicsType) {
  m <- .endotype_sqlite_tbl(omicsType)
  if(is.null(m)) return("RES-NO; unknown omics type")
  if(!requireNamespace("RSQLite", quietly = TRUE)) return("RES-NO; RSQLite not installed")
  library(DBI)
  db <- dbConnect(RSQLite::SQLite(), paste0(sqlite.path, m$db)); on.exit(dbDisconnect(db))
  if(!m$tbl %in% dbListTables(db)) return(paste0("RES-NO; table not found: ", m$tbl))
  dbReadTable(db, m$tbl)
}

compareEndotypeOmics <- function(omicsType, groupA, groupB, covariates = "", subsetKey = "",
                                 sourceFilter = "all", fdr = "true", pvalThresh = "0.05") {
  if(!requireNamespace("limma", quietly = TRUE)) return("RES-NO; limma not installed")
  library(limma)
  pvalThresh <- suppressWarnings(as.numeric(pvalThresh)); if(is.na(pvalThresh)) pvalThresh <- 0.05
  deliv <- paste0(other.tables.path, "donor_clustering_pipeline_protv2/endotype_deliverable_FINAL/")
  if(!file.exists(paste0(deliv, "donor_states.csv"))) return("RES-NO; required files not found")

  states <- read.csv(paste0(deliv, "donor_states.csv"), stringsAsFactors = FALSE)
  states$donor_id <- as.character(states$donor_id)

  # omics matrix (features x donors); metabolites = ComBat SQLite, else CSV
  mat <- .endotype_read_omics(omicsType); if(is.character(mat)) return(mat)
  donor_cols <- intersect(colnames(mat), states$donor_id)
  if(length(donor_cols) < 10) return("RES-NO; fewer than 10 donors have this omics type")
  info_cols <- setdiff(colnames(mat), donor_cols)
  symcol <- intersect(c("symbol","Symbol","gene_name","hgnc_symbol","gene","name"), info_cols)
  fid_col <- if(length(symcol)) symcol[1] else if(length(info_cols)) info_cols[1] else NULL
  fid <- if(is.null(fid_col)) as.character(seq_len(nrow(mat))) else as.character(mat[[fid_col]])
  bad <- is.na(fid) | fid == ""; if(any(bad)) fid[bad] <- paste0("feature_", which(bad))  # no NA/empty rownames (topTable rejects them)
  expr <- as.matrix(mat[, donor_cols, drop = FALSE]); rownames(expr) <- make.unique(fid)
  storage.mode(expr) <- "double"

  # metadata: state + (optional) covariates
  meta <- merge(data.frame(donor_id = donor_cols, stringsAsFactors = FALSE),
                states[, c("donor_id","state_num","source","reliability")], by = "donor_id")
  covs <- if(nzchar(covariates)) strsplit(covariates, ";")[[1]] else character(0)
  if(length(covs)){
    pheno <- read.csv(paste0(other.tables.path, "processed_comprehensive/phenotype/all_features.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
    pheno$sample_id <- as.character(pheno$sample_id)
    covs <- intersect(covs, colnames(pheno))
    if(length(covs)) meta <- merge(meta, pheno[, c("sample_id", covs)], by.x = "donor_id", by.y = "sample_id")
  }

  if(sourceFilter == "omics")         meta <- meta[meta$source == "omics", ]
  if(sourceFilter == "predicted")     meta <- meta[meta$source == "predicted", ]
  if(sourceFilter == "pred_high")     meta <- meta[meta$source == "omics" | (meta$source == "predicted" & meta$reliability == "high"), ]
  if(sourceFilter == "pred_high_med") meta <- meta[meta$source == "omics" | (meta$source == "predicted" & meta$reliability %in% c("high","medium")), ]
  if(nzchar(subsetKey)){
    sf <- paste0("donors_", subsetKey, ".rds")
    if(file.exists(sf)) meta <- meta[meta$donor_id %in% as.character(readRDS(sf)), ]
  }

  gA <- suppressWarnings(as.integer(strsplit(groupA, ";")[[1]]))
  gB <- suppressWarnings(as.integer(strsplit(groupB, ";")[[1]]))
  grp <- ifelse(meta$state_num %in% gA, "A", ifelse(meta$state_num %in% gB, "B", NA))
  meta <- meta[!is.na(grp), , drop = FALSE]; grp <- grp[!is.na(grp)]
  meta$.grp <- factor(grp, levels = c("B","A"))   # B reference; coef .grpA = A - B (log2FC)

  if(length(covs)){
    for(cv in covs){
      num <- suppressWarnings(as.numeric(meta[[cv]]))
      meta[[cv]] <- if(all(is.na(num))) factor(as.character(meta[[cv]])) else num
    }
  }
  meta <- meta[complete.cases(meta[, c(".grp", covs), drop = FALSE]), , drop = FALSE]
  if(sum(meta$.grp == "A") < 3 || sum(meta$.grp == "B") < 3) return(sprintf("RES-NO; not enough donors with this omics for the contrast: %d in Group A, %d in Group B (need at least 3 per group)", sum(meta$.grp == "A"), sum(meta$.grp == "B")))

  expr <- expr[, meta$donor_id, drop = FALSE]
  keep <- rowSums(!is.na(expr)) >= max(5, 0.5 * ncol(expr))
  expr <- expr[keep, , drop = FALSE]
  if(nrow(expr) < 5) return("RES-NO; too few measurable features")

  form <- if(length(covs)) as.formula(paste("~ .grp +", paste(sprintf("`%s`", covs), collapse = " + "))) else as.formula("~ .grp")
  design <- model.matrix(form, data = meta)
  fit <- tryCatch(eBayes(lmFit(expr, design, trend = TRUE, robust = TRUE)), error = function(e) NULL)
  if(is.null(fit)) return("RES-NO; limma fit failed")
  if(!(".grpA" %in% colnames(fit$coefficients))) return("RES-NO; group coefficient missing")
  res <- topTable(fit, coef = ".grpA", number = Inf, sort.by = "P")

  out <- data.frame(Feature = rownames(res), log2FC = round(res$logFC, 4),
                    AveExpr = round(res$AveExpr, 4), t = round(res$t, 3),
                    P_value = signif(res$P.Value, 3), Adjusted_p = signif(res$adj.P.Val, 3),
                    stringsAsFactors = FALSE)
  thr <- if(fdr == "true") out$Adjusted_p else out$P_value
  out$sig <- ifelse(thr < pvalThresh & out$log2FC > 0, "up", ifelse(thr < pvalThresh & out$log2FC < 0, "down", "NS"))
  write.csv(out, "endotype_omics_de.csv", row.names = FALSE)
  paste0("RES-OK;", sum(thr < pvalThresh, na.rm = TRUE), ";", sum(out$sig == "up"), ";", sum(out$sig == "down"), ";", nrow(out))
}

## Within-cluster omics: split ONE state's donors by a phenotype variable, then limma.
##   splitMode = "num" -> design ~ continuous splitVar (slope per feature = linear association)
##   splitMode = "cat" -> design ~ group (groupA categories vs groupB categories)
## Writes endotype_omics_de.csv (same schema as compareEndotypeOmics) so the frontend loads it unchanged.
compareWithinOmics <- function(omicsType, clusterNum, splitVar, splitMode = "num",
                               groupA = "", groupB = "", sourceFilter = "all", fdr = "true", pvalThresh = "0.05",
                               covariates = "") {
  if(!requireNamespace("limma", quietly = TRUE)) return("RES-NO; limma not installed")
  library(limma)
  pvalThresh <- suppressWarnings(as.numeric(pvalThresh)); if(is.na(pvalThresh)) pvalThresh <- 0.05
  cl <- suppressWarnings(as.integer(clusterNum)); if(is.na(cl)) return("RES-NO; bad cluster")
  deliv <- paste0(other.tables.path, "donor_clustering_pipeline_protv2/endotype_deliverable_FINAL/")
  if(!file.exists(paste0(deliv, "donor_states.csv"))) return("RES-NO; required files not found")
  states <- read.csv(paste0(deliv, "donor_states.csv"), stringsAsFactors = FALSE); states$donor_id <- as.character(states$donor_id)

  mat <- .endotype_read_omics(omicsType); if(is.character(mat)) return(mat)
  donor_cols <- intersect(colnames(mat), states$donor_id)
  if(length(donor_cols) < 10) return("RES-NO; fewer than 10 donors have this omics type")
  info_cols <- setdiff(colnames(mat), donor_cols)
  symcol <- intersect(c("symbol","Symbol","gene_name","hgnc_symbol","gene","name"), info_cols)
  fid_col <- if(length(symcol)) symcol[1] else if(length(info_cols)) info_cols[1] else NULL
  fid <- if(is.null(fid_col)) as.character(seq_len(nrow(mat))) else as.character(mat[[fid_col]])
  bad <- is.na(fid) | fid == ""; if(any(bad)) fid[bad] <- paste0("feature_", which(bad))
  expr <- as.matrix(mat[, donor_cols, drop = FALSE]); rownames(expr) <- make.unique(fid); storage.mode(expr) <- "double"

  meta <- merge(data.frame(donor_id = donor_cols, stringsAsFactors = FALSE),
                states[, c("donor_id","state_num","source","reliability")], by = "donor_id")
  pheno <- read.csv(paste0(other.tables.path, "processed_comprehensive/phenotype/all_features.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE); pheno$sample_id <- as.character(pheno$sample_id)
  if(!(splitVar %in% colnames(pheno))) return(paste0("RES-NO; split variable not found: ", splitVar))
  meta <- merge(meta, pheno[, c("sample_id", splitVar)], by.x = "donor_id", by.y = "sample_id")
  if(splitVar %in% c("diagnosis","diagnosis_computed")){ vv <- as.character(meta[[splitVar]]); vv[is.na(vv) | vv == ""] <- "none"; meta[[splitVar]] <- vv }

  # optional numeric covariates: merge from all_features and stage as .cov1, .cov2, ...
  covs <- if(nzchar(covariates)) setdiff(strsplit(covariates, ";")[[1]], splitVar) else character(0)
  covs <- intersect(covs, colnames(pheno)); covCols <- character(0)
  if(length(covs)){
    meta <- merge(meta, pheno[, c("sample_id", covs)], by.x = "donor_id", by.y = "sample_id")
    for(i in seq_along(covs)){ cn <- paste0(".cov", i); meta[[cn]] <- suppressWarnings(as.numeric(meta[[covs[i]]])); covCols <- c(covCols, cn) }
  }

  if(sourceFilter == "omics")         meta <- meta[meta$source == "omics", ]
  if(sourceFilter == "predicted")     meta <- meta[meta$source == "predicted", ]
  if(sourceFilter == "pred_high")     meta <- meta[meta$source == "omics" | (meta$source == "predicted" & meta$reliability == "high"), ]
  if(sourceFilter == "pred_high_med") meta <- meta[meta$source == "omics" | (meta$source == "predicted" & meta$reliability %in% c("high","medium")), ]
  meta <- meta[meta$state_num == cl, , drop = FALSE]
  if(nrow(meta) < 6) return("RES-NO; too few donors in this state for the selected donor set")

  covTerm <- if(length(covCols)) paste0(" + ", paste(covCols, collapse = " + ")) else ""
  if(splitMode == "cat"){
    A <- strsplit(groupA, ";")[[1]]; B <- strsplit(groupB, ";")[[1]]
    g <- ifelse(as.character(meta[[splitVar]]) %in% A, "A", ifelse(as.character(meta[[splitVar]]) %in% B, "B", NA))
    meta <- meta[!is.na(g), , drop = FALSE]; g <- g[!is.na(g)]
    if(sum(g == "A") < 3 || sum(g == "B") < 3) return(sprintf("RES-NO; not enough donors per group: %d vs %d", sum(g == "A"), sum(g == "B")))
    meta$.grp <- factor(g, levels = c("B","A")); coefName <- ".grpA"; form <- as.formula(paste0("~ .grp", covTerm))
  } else {
    meta$.x <- suppressWarnings(as.numeric(meta[[splitVar]]))
    meta <- meta[!is.na(meta$.x), , drop = FALSE]
    if(nrow(meta) < 6) return("RES-NO; too few donors with this variable measured")
    coefName <- ".x"; form <- as.formula(paste0("~ .x", covTerm))
  }
  if(length(covCols)){
    ok <- stats::complete.cases(meta[, covCols, drop = FALSE]); meta <- meta[ok, , drop = FALSE]
    if(nrow(meta) < 6) return("RES-NO; too few donors with the covariate(s) measured")
    if(splitMode == "cat" && (sum(meta$.grp == "A") < 3 || sum(meta$.grp == "B") < 3)) return("RES-NO; not enough donors per group after covariate filtering")
  }
  expr <- expr[, meta$donor_id, drop = FALSE]
  keep <- rowSums(!is.na(expr)) >= max(5, 0.5 * ncol(expr)); expr <- expr[keep, , drop = FALSE]
  if(nrow(expr) < 5) return("RES-NO; too few measurable features")
  design <- model.matrix(form, data = meta)
  if(!(coefName %in% colnames(design))) return("RES-NO; design coefficient missing")
  fit <- tryCatch(eBayes(lmFit(expr, design, trend = TRUE, robust = TRUE)), error = function(e) NULL)
  if(is.null(fit)) return("RES-NO; limma fit failed")
  res <- topTable(fit, coef = coefName, number = Inf, sort.by = "P")
  out <- data.frame(Feature = rownames(res), log2FC = round(res$logFC, 4), AveExpr = round(res$AveExpr, 4),
                    t = round(res$t, 3), P_value = signif(res$P.Value, 3), Adjusted_p = signif(res$adj.P.Val, 3), stringsAsFactors = FALSE)
  thr <- if(fdr == "true") out$Adjusted_p else out$P_value
  out$sig <- ifelse(thr < pvalThresh & out$log2FC > 0, "up", ifelse(thr < pvalThresh & out$log2FC < 0, "down", "NS"))
  write.csv(out, "endotype_omics_de.csv", row.names = FALSE)
  paste0("RES-OK;", sum(thr < pvalThresh, na.rm = TRUE), ";", sum(out$sig == "up"), ";", sum(out$sig == "down"), ";", nrow(out))
}

## ============================================================================
## endotypeOmicsPathway -- pathway enrichment for the Donor State Explorer.
## Reuses HumanIsletsR performGSEA (fgsea preranked, ranked by the limma t-statistic)
## across one or more gene-set / metabolite-set libraries, on the SAME state-vs-group
## contrast as compareEndotypeOmics. Genes -> kegg / reactome / go_mf (entrez Gene_ID);
## metabolites -> hsa_kegg (KEGG compound sets). Writes endotype_omics_pathway.csv.
## No MetaboAnalystR -- HumanIsletsR machinery only.
##   libs : ";"-joined funcLib values (gene libraries; ignored for metabolites)
## ============================================================================

## shared helper: run the endotype limma contrast and return a ranked data.frame
## carrying entrez Gene_ID + kegg_id annotation (or a "RES-NO; ..." string on failure).
.endotypeOmicsRanked <- function(omicsType, groupA, groupB, covariates = "", subsetKey = "",
                                 sourceFilter = "all") {
  if(!requireNamespace("limma", quietly = TRUE)) return("RES-NO; limma not installed")
  library(limma)
  deliv <- paste0(other.tables.path, "donor_clustering_pipeline_protv2/endotype_deliverable_FINAL/")
  if(!file.exists(paste0(deliv, "donor_states.csv"))) return("RES-NO; required files not found")

  states <- read.csv(paste0(deliv, "donor_states.csv"), stringsAsFactors = FALSE)
  states$donor_id <- as.character(states$donor_id)

  mat <- .endotype_read_omics(omicsType); if(is.character(mat)) return(mat)
  donor_cols <- intersect(colnames(mat), states$donor_id)
  if(length(donor_cols) < 10) return("RES-NO; fewer than 10 donors have this omics type")
  info_cols <- setdiff(colnames(mat), donor_cols)
  symcol <- intersect(c("symbol","Symbol","gene_name","hgnc_symbol","gene","name"), info_cols)
  fid_col <- if(length(symcol)) symcol[1] else if(length(info_cols)) info_cols[1] else NULL
  fid <- if(is.null(fid_col)) as.character(seq_len(nrow(mat))) else as.character(mat[[fid_col]])
  bad <- is.na(fid) | fid == ""; if(any(bad)) fid[bad] <- paste0("feature_", which(bad))  # no NA/empty rownames (topTable rejects them)
  uid <- make.unique(fid)
  expr <- as.matrix(mat[, donor_cols, drop = FALSE]); rownames(expr) <- uid
  storage.mode(expr) <- "double"

  gid_col <- intersect(c("gene_id","Gene_ID","entrez","entrezgene"), info_cols)
  keg_col <- intersect(c("kegg_id"), info_cols)
  anno <- data.frame(Feature = uid,
                     Gene_ID = if(length(gid_col)) as.character(mat[[gid_col[1]]]) else NA_character_,
                     kegg_id = if(length(keg_col)) as.character(mat[[keg_col[1]]]) else NA_character_,
                     stringsAsFactors = FALSE)

  meta <- merge(data.frame(donor_id = donor_cols, stringsAsFactors = FALSE),
                states[, c("donor_id","state_num","source","reliability")], by = "donor_id")
  covs <- if(nzchar(covariates)) strsplit(covariates, ";")[[1]] else character(0)
  if(length(covs)){
    pheno <- read.csv(paste0(other.tables.path, "processed_comprehensive/phenotype/all_features.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
    pheno$sample_id <- as.character(pheno$sample_id)
    covs <- intersect(covs, colnames(pheno))
    if(length(covs)) meta <- merge(meta, pheno[, c("sample_id", covs)], by.x = "donor_id", by.y = "sample_id")
  }

  if(sourceFilter == "omics")         meta <- meta[meta$source == "omics", ]
  if(sourceFilter == "predicted")     meta <- meta[meta$source == "predicted", ]
  if(sourceFilter == "pred_high")     meta <- meta[meta$source == "omics" | (meta$source == "predicted" & meta$reliability == "high"), ]
  if(sourceFilter == "pred_high_med") meta <- meta[meta$source == "omics" | (meta$source == "predicted" & meta$reliability %in% c("high","medium")), ]
  if(nzchar(subsetKey)){
    sf <- paste0("donors_", subsetKey, ".rds")
    if(file.exists(sf)) meta <- meta[meta$donor_id %in% as.character(readRDS(sf)), ]
  }

  gA <- suppressWarnings(as.integer(strsplit(groupA, ";")[[1]]))
  gB <- suppressWarnings(as.integer(strsplit(groupB, ";")[[1]]))
  grp <- ifelse(meta$state_num %in% gA, "A", ifelse(meta$state_num %in% gB, "B", NA))
  meta <- meta[!is.na(grp), , drop = FALSE]; grp <- grp[!is.na(grp)]
  meta$.grp <- factor(grp, levels = c("B","A"))

  if(length(covs)){
    for(cv in covs){
      num <- suppressWarnings(as.numeric(meta[[cv]]))
      meta[[cv]] <- if(all(is.na(num))) factor(as.character(meta[[cv]])) else num
    }
  }
  meta <- meta[complete.cases(meta[, c(".grp", covs), drop = FALSE]), , drop = FALSE]
  if(sum(meta$.grp == "A") < 3 || sum(meta$.grp == "B") < 3) return(sprintf("RES-NO; not enough donors with this omics for the contrast: %d in Group A, %d in Group B (need at least 3 per group)", sum(meta$.grp == "A"), sum(meta$.grp == "B")))

  expr <- expr[, meta$donor_id, drop = FALSE]
  keep <- rowSums(!is.na(expr)) >= max(5, 0.5 * ncol(expr))
  expr <- expr[keep, , drop = FALSE]
  if(nrow(expr) < 5) return("RES-NO; too few measurable features")

  form <- if(length(covs)) as.formula(paste("~ .grp +", paste(sprintf("`%s`", covs), collapse = " + "))) else as.formula("~ .grp")
  design <- model.matrix(form, data = meta)
  fit <- tryCatch(eBayes(lmFit(expr, design, trend = TRUE, robust = TRUE)), error = function(e) NULL)
  if(is.null(fit)) return("RES-NO; limma fit failed")
  if(!(".grpA" %in% colnames(fit$coefficients))) return("RES-NO; group coefficient missing")
  res <- topTable(fit, coef = ".grpA", number = Inf, sort.by = "P")

  out <- data.frame(Feature = rownames(res), log2FC = round(res$logFC, 4),
                    T_statistic = round(res$t, 4), stringsAsFactors = FALSE)
  out <- merge(out, anno, by = "Feature", all.x = TRUE, sort = FALSE)
  out
}

## Load a HumanIslets gene-set / metabolite-set library as list(sets = named member-id vectors,
## term = pathwayID -> readable name). Genes: kegg/reactome/go_mf .rds ($sets entrez + $term);
## metabolites: kegg_hsa_met.qs ($mset.list keyed by hsa IDs, $path.ids name<->id). Same libraries
## the Omics page (performGSEA) uses.
.endotype_load_lib <- function(funcLib){
  lib.path <- paste0(other.tables.path, "libraries/")
  if(funcLib == "hsa_kegg"){
    lib  <- qs::qread(paste0(lib.path, "kegg_hsa_met.qs"))
    sets <- lib$mset.list
    nm   <- setNames(names(lib$path.ids), lib$path.ids)        # hsa00010 -> "Glycolysis ..."
    terms <- nm[names(sets)]; terms[is.na(terms)] <- names(sets)[is.na(terms)]
    return(list(sets = sets, term = setNames(unname(terms), names(sets))))
  }
  f <- paste0(lib.path, funcLib, ".rds")
  if(!file.exists(f)) return(NULL)
  rds <- readRDS(f)
  list(sets = rds$sets, term = setNames(rds$term, names(rds$sets)))
}

## Donor State Explorer: per-donor measured values for ONE omics feature (box plot).
## Writes endotype_omics_feature.csv (donor_id, value); the frontend groups the donors by the
## active contrast (states / within-state split) since it already holds the donor metadata.
endotypeOmicsFeature <- function(omicsType, feature) {
  mat <- .endotype_read_omics(omicsType); if(is.character(mat)) return(mat)
  deliv <- paste0(other.tables.path, "donor_clustering_pipeline_protv2/endotype_deliverable_FINAL/")
  if(!file.exists(paste0(deliv, "donor_states.csv"))) return("RES-NO; required files not found")
  states <- read.csv(paste0(deliv, "donor_states.csv"), stringsAsFactors = FALSE); states$donor_id <- as.character(states$donor_id)
  donor_cols <- intersect(colnames(mat), states$donor_id)
  if(!length(donor_cols)) return("RES-NO; no donors have this omics type")
  info_cols <- setdiff(colnames(mat), donor_cols)
  symcol <- intersect(c("symbol","Symbol","gene_name","hgnc_symbol","gene","name"), info_cols)
  fid_col <- if(length(symcol)) symcol[1] else if(length(info_cols)) info_cols[1] else NULL
  fid <- if(is.null(fid_col)) as.character(seq_len(nrow(mat))) else as.character(mat[[fid_col]])
  bad <- is.na(fid) | fid == ""; if(any(bad)) fid[bad] <- paste0("feature_", which(bad))
  uid <- make.unique(fid)                                   # SAME id construction as the DE tables
  idx <- match(feature, uid)
  if(is.na(idx)) return("RES-NO; feature not found in this omics type")
  vals <- suppressWarnings(as.numeric(as.matrix(mat[idx, donor_cols, drop = FALSE])[1, ]))
  out <- data.frame(donor_id = donor_cols, value = round(vals, 5), stringsAsFactors = FALSE)
  out <- out[!is.na(out$value), , drop = FALSE]
  if(!nrow(out)) return("RES-NO; no measured values for this feature")
  write.csv(out, "endotype_omics_feature.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(out))
}

## within-state analog of .endotypeOmicsRanked: split ONE state's donors by a phenotype
## variable (cat -> Group A vs B; num -> per-feature slope), optionally adjusting for numeric
## covariates, and return the SAME ranked data.frame (Feature, log2FC, T_statistic, Gene_ID,
## kegg_id) so the pathway machinery can rank by T_statistic exactly as for the cross-state case.
.endotypeWithinRanked <- function(omicsType, clusterNum, splitVar, splitMode = "num",
                                  groupA = "", groupB = "", sourceFilter = "all", covariates = "") {
  if(!requireNamespace("limma", quietly = TRUE)) return("RES-NO; limma not installed")
  library(limma)
  cl <- suppressWarnings(as.integer(clusterNum)); if(is.na(cl)) return("RES-NO; bad cluster")
  deliv <- paste0(other.tables.path, "donor_clustering_pipeline_protv2/endotype_deliverable_FINAL/")
  if(!file.exists(paste0(deliv, "donor_states.csv"))) return("RES-NO; required files not found")
  states <- read.csv(paste0(deliv, "donor_states.csv"), stringsAsFactors = FALSE); states$donor_id <- as.character(states$donor_id)

  mat <- .endotype_read_omics(omicsType); if(is.character(mat)) return(mat)
  donor_cols <- intersect(colnames(mat), states$donor_id)
  if(length(donor_cols) < 10) return("RES-NO; fewer than 10 donors have this omics type")
  info_cols <- setdiff(colnames(mat), donor_cols)
  symcol <- intersect(c("symbol","Symbol","gene_name","hgnc_symbol","gene","name"), info_cols)
  fid_col <- if(length(symcol)) symcol[1] else if(length(info_cols)) info_cols[1] else NULL
  fid <- if(is.null(fid_col)) as.character(seq_len(nrow(mat))) else as.character(mat[[fid_col]])
  bad <- is.na(fid) | fid == ""; if(any(bad)) fid[bad] <- paste0("feature_", which(bad))
  uid <- make.unique(fid)
  expr <- as.matrix(mat[, donor_cols, drop = FALSE]); rownames(expr) <- uid; storage.mode(expr) <- "double"

  gid_col <- intersect(c("gene_id","Gene_ID","entrez","entrezgene"), info_cols)
  keg_col <- intersect(c("kegg_id"), info_cols)
  anno <- data.frame(Feature = uid,
                     Gene_ID = if(length(gid_col)) as.character(mat[[gid_col[1]]]) else NA_character_,
                     kegg_id = if(length(keg_col)) as.character(mat[[keg_col[1]]]) else NA_character_,
                     stringsAsFactors = FALSE)

  meta <- merge(data.frame(donor_id = donor_cols, stringsAsFactors = FALSE),
                states[, c("donor_id","state_num","source","reliability")], by = "donor_id")
  pheno <- read.csv(paste0(other.tables.path, "processed_comprehensive/phenotype/all_features.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE); pheno$sample_id <- as.character(pheno$sample_id)
  if(!(splitVar %in% colnames(pheno))) return(paste0("RES-NO; split variable not found: ", splitVar))
  meta <- merge(meta, pheno[, c("sample_id", splitVar)], by.x = "donor_id", by.y = "sample_id")
  if(splitVar %in% c("diagnosis","diagnosis_computed")){ vv <- as.character(meta[[splitVar]]); vv[is.na(vv) | vv == ""] <- "none"; meta[[splitVar]] <- vv }

  covs <- if(nzchar(covariates)) setdiff(strsplit(covariates, ";")[[1]], splitVar) else character(0)
  covs <- intersect(covs, colnames(pheno)); covCols <- character(0)
  if(length(covs)){
    meta <- merge(meta, pheno[, c("sample_id", covs)], by.x = "donor_id", by.y = "sample_id")
    for(i in seq_along(covs)){ cn <- paste0(".cov", i); meta[[cn]] <- suppressWarnings(as.numeric(meta[[covs[i]]])); covCols <- c(covCols, cn) }
  }

  if(sourceFilter == "omics")         meta <- meta[meta$source == "omics", ]
  if(sourceFilter == "predicted")     meta <- meta[meta$source == "predicted", ]
  if(sourceFilter == "pred_high")     meta <- meta[meta$source == "omics" | (meta$source == "predicted" & meta$reliability == "high"), ]
  if(sourceFilter == "pred_high_med") meta <- meta[meta$source == "omics" | (meta$source == "predicted" & meta$reliability %in% c("high","medium")), ]
  meta <- meta[meta$state_num == cl, , drop = FALSE]
  if(nrow(meta) < 6) return("RES-NO; too few donors in this state for the selected donor set")

  covTerm <- if(length(covCols)) paste0(" + ", paste(covCols, collapse = " + ")) else ""
  if(splitMode == "cat"){
    A <- strsplit(groupA, ";")[[1]]; B <- strsplit(groupB, ";")[[1]]
    g <- ifelse(as.character(meta[[splitVar]]) %in% A, "A", ifelse(as.character(meta[[splitVar]]) %in% B, "B", NA))
    meta <- meta[!is.na(g), , drop = FALSE]; g <- g[!is.na(g)]
    if(sum(g == "A") < 3 || sum(g == "B") < 3) return(sprintf("RES-NO; not enough donors per group: %d vs %d", sum(g == "A"), sum(g == "B")))
    meta$.grp <- factor(g, levels = c("B","A")); coefName <- ".grpA"; form <- as.formula(paste0("~ .grp", covTerm))
  } else {
    meta$.x <- suppressWarnings(as.numeric(meta[[splitVar]]))
    meta <- meta[!is.na(meta$.x), , drop = FALSE]
    if(nrow(meta) < 6) return("RES-NO; too few donors with this variable measured")
    coefName <- ".x"; form <- as.formula(paste0("~ .x", covTerm))
  }
  if(length(covCols)){
    ok <- stats::complete.cases(meta[, covCols, drop = FALSE]); meta <- meta[ok, , drop = FALSE]
    if(nrow(meta) < 6) return("RES-NO; too few donors with the covariate(s) measured")
    if(splitMode == "cat" && (sum(meta$.grp == "A") < 3 || sum(meta$.grp == "B") < 3)) return("RES-NO; not enough donors per group after covariate filtering")
  }
  expr <- expr[, meta$donor_id, drop = FALSE]
  keep <- rowSums(!is.na(expr)) >= max(5, 0.5 * ncol(expr)); expr <- expr[keep, , drop = FALSE]
  if(nrow(expr) < 5) return("RES-NO; too few measurable features")
  design <- model.matrix(form, data = meta)
  if(!(coefName %in% colnames(design))) return("RES-NO; design coefficient missing")
  fit <- tryCatch(eBayes(lmFit(expr, design, trend = TRUE, robust = TRUE)), error = function(e) NULL)
  if(is.null(fit)) return("RES-NO; limma fit failed")
  res <- topTable(fit, coef = coefName, number = Inf, sort.by = "P")
  out <- data.frame(Feature = rownames(res), log2FC = round(res$logFC, 4),
                    T_statistic = round(res$t, 4), stringsAsFactors = FALSE)
  out <- merge(out, anno, by = "Feature", all.x = TRUE, sort = FALSE)
  out
}

## endotypeOmicsPathway: cross-state by default; when splitVar is supplied it instead enriches the
## WITHIN-state contrast (one state split by splitVar), so pathway follows the same DE the table shows.
endotypeOmicsPathway <- function(omicsType, groupA, groupB, covariates = "", subsetKey = "",
                                 sourceFilter = "all", libs = "kegg;reactome;go_mf",
                                 fdr = "0.05", collapse = "true",
                                 clusterNum = "", splitVar = "", splitMode = "num") {
  if(!requireNamespace("fgsea", quietly = TRUE)) return("RES-NO; fgsea not installed")
  library(fgsea); set.seed(42)
  is_metab   <- omicsType %in% c("proc_metabolite","proc_metabolite_lg","proc_metabolite_hg","proc_metabolite_ratio")
  gene_omics <- omicsType %in% c("proc_rnaseq","proc_prot_v2","proc_nanostring","proc_pbrna_alpha","proc_pbrna_beta")
  if(!is_metab && !gene_omics) return("RES-NO; pathway enrichment not available for this omics type")

  r <- if(nzchar(splitVar)) .endotypeWithinRanked(omicsType, clusterNum, splitVar, splitMode, groupA, groupB, sourceFilter, covariates)
       else .endotypeOmicsRanked(omicsType, groupA, groupB, covariates, subsetKey, sourceFilter)
  if(is.character(r)) return(r)

  # Build the preranked stats keyed by entrez (genes) / kegg_id (metabolites); map id -> display name.
  keycol <- if(is_metab) "kegg_id" else "Gene_ID"
  rr <- r[!is.na(r[[keycol]]) & r[[keycol]] != "" & !is.na(r$T_statistic), , drop = FALSE]
  rr <- rr[!duplicated(rr[[keycol]]), , drop = FALSE]
  if(nrow(rr) < 5) return("RES-NO; too few mapped features for pathway analysis")
  ranks  <- setNames(as.numeric(rr$T_statistic), as.character(rr[[keycol]]))
  id2name <- setNames(rr$Feature, as.character(rr[[keycol]]))   # entrez->symbol or kegg->compound

  lib_vec <- if(is_metab) "hsa_kegg" else strsplit(libs, ";")[[1]]
  lib_vec <- lib_vec[nzchar(lib_vec)]
  if(!length(lib_vec)) return("RES-NO; no libraries requested")
  fdrn <- suppressWarnings(as.numeric(fdr)); if(is.na(fdrn)) fdrn <- 0.05
  lib_label <- c(kegg = "KEGG", reactome = "Reactome", go_mf = "GO Molecular Function",
                 go_bp = "GO Biological Process", go_cc = "GO Cellular Component", hsa_kegg = "KEGG")

  acc <- list()
  for(lib in lib_vec){
    L <- tryCatch(.endotype_load_lib(lib), error = function(e) NULL)
    if(is.null(L) || !length(L$sets)) next
    minSize <- if(lib == "hsa_kegg") 3 else 15
    fres <- tryCatch(fgsea(pathways = L$sets, stats = ranks, minSize = minSize, maxSize = 500),
                     error = function(e) NULL)
    if(is.null(fres) || !nrow(fres)) next
    fres <- fres[!is.na(fres$pval) & !is.na(fres$NES) & !is.na(fres$padj), ]   # drop pathways with no computable stats
    if(!nrow(fres)) next
    fres <- fres[order(fres$pval), ]
    if(collapse == "true"){
      ind <- tryCatch(collapsePathways(fres, L$sets, ranks), error = function(e) NULL)
      if(!is.null(ind)) fres <- fres[fres$pathway %in% ind$mainPathways, ]
    }
    if(!nrow(fres)) next
    src  <- if(lib %in% names(lib_label)) unname(lib_label[lib]) else lib
    term <- L$term[as.character(fres$pathway)]; term[is.na(term)] <- as.character(fres$pathway)[is.na(term)]
    # leading-edge members -> display names (gene symbols / metabolite names)
    le_names <- vapply(fres$leadingEdge, function(ids){
      nm <- id2name[as.character(ids)]; nm <- nm[!is.na(nm) & nm != ""]
      if(!length(nm)) nm <- as.character(ids)
      if(length(nm) > 50) nm <- c(nm[1:50], paste0("(+", length(nm) - 50, " more)"))  # Hits keeps the true count
      paste(nm, collapse = ", ")
    }, character(1))
    acc[[length(acc) + 1]] <- data.frame(
      Term = unname(term), Source = src,
      Direction = ifelse(fres$NES > 0, "up", "down"), NES = round(fres$NES, 3),
      P_value = signif(fres$pval, 3), Adjusted_p = signif(fres$padj, 3),
      Hits = lengths(fres$leadingEdge), SetSize = as.integer(fres$size),
      Genes = le_names, stringsAsFactors = FALSE)
  }
  if(!length(acc)) return("RES-NO; no pathways enriched for this contrast")
  comb <- do.call(rbind, acc)
  comb <- comb[order(comb$P_value), ]
  write.csv(comb, "endotype_omics_pathway.csv", row.names = FALSE)
  paste0("RES-OK;", nrow(comb), ";", sum(comb$Adjusted_p < fdrn, na.rm = TRUE))
}
 
