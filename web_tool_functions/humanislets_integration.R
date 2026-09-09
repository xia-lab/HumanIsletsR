# R Functions for metadata analysis on HumanIslets web tool V2
# Author: Yao Lu

################################################################################
# ---------------------------------------------------------------------------
# load_integ_omics(): read ONE omics block for the integration from the single
# canonical SQLite source (HI_omics_v2.sqlite), via the already-open connection
# `mydb`. Returns a feature x donor data.frame (rownames = the block feature id)
# or NULL for an unrecognised type. Proteomics + metabolite use the ComBat
# batch-corrected tables; info columns are dropped and donor columns kept.
# ---------------------------------------------------------------------------
load_integ_omics <- function(mydb, omicsType, cell = "Alpha", metaboGluc = "LG",
                             tissueType = "", LOD_corrected = "", contamClass = "all") {
  drop_info <- function(ft, info_cols) ft[, !(colnames(ft) %in% info_cols), drop = FALSE]
  cond <- if (metaboGluc %in% c("HG_LG_ratio", "ratio")) "ratio" else if (metaboGluc %in% c("HG", "LG")) metaboGluc else "ratio"

  if (omicsType == "proc_rnaseq") {
    ft <- dbReadTable(mydb, "proc_rnaseq")
    rownames(ft) <- ft$accession
    # Canonical feature -> entrez(gene_id) + symbol mapping, straight from this table (attached as an
    # attribute so the pathway code uses YOUR IDs instead of re-matching an external annotation).
    an <- data.frame(feature = rownames(ft), entrez = as.character(ft$gene_id), symbol = as.character(ft$symbol), stringsAsFactors = FALSE)
    mat <- drop_info(ft, c("accession", "gene_id", "symbol", "name", "biotype"))
    attr(mat, "feat_anno") <- an; mat

  } else if (omicsType == "proc_prot_combat") {
    ft <- dbReadTable(mydb, "proc_prot_combat")
    sym <- as.character(ft$symbol)
    empty <- is.na(sym) | sym == ""
    sym[empty] <- as.character(ft$id[empty])
    rn <- make.unique(sym)
    # Canonical feature -> entrez(gene_id) + symbol mapping (captured before the columns are dropped).
    an <- data.frame(feature = rn, entrez = as.character(ft$gene_id), symbol = as.character(ft$symbol), stringsAsFactors = FALSE)
    mat <- drop_info(ft, c("id", "gene_id", "symbol", "name", "Protein_Group"))
    rownames(mat) <- rn
    attr(mat, "feat_anno") <- an; mat

  } else if (omicsType == "proc_nanostring_merge") {
    ft <- dbReadTable(mydb, "proc_nanostring_merge")
    rownames(ft) <- ft$gene_id
    an <- data.frame(feature = rownames(ft), entrez = as.character(ft$gene_id), symbol = as.character(ft$symbol), stringsAsFactors = FALSE)
    mat <- drop_info(ft, c("gene_id", "symbol", "ensembl", "name"))
    attr(mat, "feat_anno") <- an; mat

  } else if (omicsType == "proc_pbrna") {
    ft <- dbReadTable(mydb, paste0("proc_pbrna_", cell))
    rownames(ft) <- make.unique(as.character(ft$symbol))
    an <- data.frame(feature = rownames(ft), entrez = as.character(ft$gene_id), symbol = as.character(ft$symbol), stringsAsFactors = FALSE)
    mat <- drop_info(ft, c("symbol", "gene_id", "ensembl", "name"))
    attr(mat, "feat_anno") <- an; mat

  } else if (omicsType == "proc_metabolite") {
    ft <- dbReadTable(mydb, paste0("proc_metabolite_combat_", cond))
    rownames(ft) <- ft$compound
    drop_info(ft, c("compound", "inchikey", "ik_conn", "hmdb_id", "kegg_id",
                    "gem_id", "super_class", "main_class", "sub_class"))

  } else if (omicsType == "proc_flux") {
    ft <- dbReadTable(mydb, paste0("proc_flux_", cond))
    rownames(ft) <- ft$rxn
    drop_info(ft, c("rxn", "genes", "subsystem", "description"))

  } else if (omicsType == "proc_contaminants") {
    ft <- dbReadTable(mydb, "proc_contaminants")
    ft <- ft[ft$Tissue == tissueType & ft$LOG_corrected == LOD_corrected, ]
    if (contamClass != "all") ft <- ft[ft$Class == contamClass, ]
    rownames(ft) <- ft$Compound
    drop_info(ft, c("LOG_corrected", "Tissue", "Compound", "Class"))

  } else {
    NULL
  }
}


# Wrapper: return any R error as a readable string so it reaches the web interface
# (instead of throwing -> generic "Processing failed"). Real work is in GetIntegPlot_impl below.
GetIntegPlot <- function(...) {
  tryCatch(GetIntegPlot_impl(...), error = function(e) paste0("Error: DIABLO - ", conditionMessage(e)))
}

GetIntegPlot_impl <- function(
    omics1, omics2,omics3,
    varGroup,
    analysisVar,
    class1 = NULL,
    class2 = NULL,
    primaryType = "disc", 
    donors = "all",
    cell = "alpha",
    glucose = "1", 
    batch="b1",
    peakmode="positive",
    metaboGluc="LG",
    LOD_corrected = "true",
    tissueType = "Adipose",
    contamClass = "Dioxin",
    covariates = "",
    topFeat = 100,
    mode = "local"){
 
  # Headless-server graphics hardening: close any device left open by an earlier crashed run
  # (Rserve reuses one R session, so open devices accumulate) and force the Cairo PNG backend
  # so png()/ggsave don't need X11 -> fixes "unable to start device PNG".
  graphics.off()
  if (capabilities("cairo")) options(bitmapType = "cairo")

  if(length(class2) == 0){class2 = 'NULL'}
 
  if(!file.exists("fig_counts.qs")){
    fig_counts <- 0
  }else{
   fig_counts <- qs::qread("fig_counts.qs")
    fig_counts <- fig_counts+1
  }
  qs::qsave(fig_counts,"fig_counts.qs") 
  library(dplyr)

selected_omics <- unique(c(omics1, omics2, omics3))
  selected_omics <- selected_omics[selected_omics != "" & !is.na(selected_omics)]
   selected_omics <- selected_omics[!selected_omics %in% c("NA", "null", "NULL")]
  
  if(length(selected_omics) == 0) {
    stop("No valid omics selected. Please select at least one omics type.")
  }
  

  dataList <- list()

  # All omics blocks are read from the single canonical source (SQLite) over one
  # shared connection. Proteomics + metabolite use the ComBat batch-corrected tables.
  library(RSQLite)
  mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_omics_v2.sqlite"))
  feat_anno_all <- list()
  for(ot in selected_omics){
    blk <- load_integ_omics(mydb, ot, cell = cell, metaboGluc = metaboGluc,
                            tissueType = tissueType, LOD_corrected = LOD_corrected,
                            contamClass = contamClass)
    if(!is.null(blk)) {
      dataList[[ot]] <- blk
      if(!is.null(attr(blk, "feat_anno"))) feat_anno_all[[ot]] <- attr(blk, "feat_anno")
    }
  }
  dbDisconnect(mydb)
  # Persist the canonical feature -> entrez + symbol mapping (from the proc tables' own gene_id/symbol
  # columns) for the pathway ORA/GSEA steps -- so they use your IDs, not a re-matched annotation file.
  if(length(feat_anno_all) > 0) {
    feature_anno <- do.call(rbind, Map(function(df, nm) { df$omics <- nm; df }, feat_anno_all, names(feat_anno_all)))
    feature_anno <- feature_anno[!is.na(feature_anno$feature) & feature_anno$feature != "", ]
    qs::qsave(feature_anno, "feature_anno.qs")
  }

  
  dataList <- dataList[lengths(dataList) > 0] 
 

  if(length(dataList)<2){
    return("NO;Not enough omics data! Please make sure at least two omics datasets has been selected.")
  }

  idx <- match(selected_omics, names(dataList))
  idx <- idx[!is.na(idx)]  
  dataList <- dataList[idx]


meta <- read.csv(paste0(other.tables.path, "display_data/metadata_sum_norm.csv"))
  meta_full <- meta   # keep all columns for optional covariate adjustment
  # Inject cluster (or any proc_metadata-only) column if selected var is missing from the CSV
  if (!is.null(analysisVar) && !(analysisVar %in% colnames(meta))) {
    library(RSQLite)
    meta_db <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
    pheno <- dbReadTable(meta_db, "proc_metadata")
    dbDisconnect(meta_db)
    if (analysisVar %in% colnames(pheno)) {
      extra <- pheno[, c("record_id", analysisVar)]
      meta <- merge(meta, extra, by = "record_id", all.x = TRUE)
    } else {
      stop(paste0("Analysis failed: variable '", analysisVar,
                  "' not found in metadata or donor table."))
    }
  }
  shared <- Reduce(intersect, lapply(dataList, colnames))


  # Multi-omics page uses donors_multiomics.rds (isolated from omics/phenotype views' donors.rds).
  if(donors == "subset" && file.exists("donors_multiomics.rds")){
    donor.list <- readRDS("donors_multiomics.rds")
    shared <- intersect(donor.list, shared)
  }


if(primaryType == "disc"){
 
    meta <- meta[meta[[analysisVar]] %in% c(class1, class2), c("record_id", analysisVar)]
    shared <- shared[shared %in% meta$record_id]
 
    meta <- meta[match(shared, meta$record_id), ]
 
    groups <- factor(meta[[analysisVar]], levels = c(class1, class2))
    # Relabel group levels to their display names (disc_groups.csv) so plots/legends/messages
    # show e.g. "State 2 Stress" instead of the raw id "C0". Falls back to the raw id if unmapped.
    disc_map <- tryCatch(read.csv(paste0(other.tables.path, "display_interface/disc_groups.csv"),
                                  stringsAsFactors = FALSE), error = function(e) NULL)
    if(!is.null(disc_map)){
      sub <- disc_map[disc_map$column == analysisVar, ]
      if(nrow(sub) > 0){
        disp <- sub$display[match(levels(groups), sub$group)]
        disp[is.na(disp)] <- levels(groups)[is.na(disp)]
        levels(groups) <- disp
      }
    }
    counts <- table(groups)
    lev = levels(groups)
    
    # Check 1: Are there actually two groups?
    if(any(counts == 0)){
          missing_group <- names(counts)[counts == 0][1]
        stop(paste0("Analysis failed: Group '", missing_group, "' has zero common samples across omics. Two groups are required."))
    }
    
    # Check 2: Is any group size too small for DIABLO?
    min_found <- min(counts)
    if(min_found < 3){
        too_small_group <- names(counts)[which.min(counts)]
        stop(paste0("Analysis failed: Group '", too_small_group, "' only has ", min_found, 
                    " samples. A minimum of 3 is required for DIABLO."))
    }
    Y <- groups
names(Y) <- shared
} else {
 
    meta <- meta[match(shared, meta$record_id), c("record_id", analysisVar)]
     keep_idx <- !is.na(meta[[analysisVar]])
    meta <- meta[keep_idx, ]
    shared <- shared[keep_idx]
      Y <- matrix( as.numeric(as.character(meta[[analysisVar]])));
    rownames(Y) <- shared
}

 
  if (length(shared) < 10){
    stop("Too few common samples across omics and groups.")
    return("NO;Too few common samples across omics and groups.");
  }

# Final Safety Check
stopifnot(all(meta$record_id == shared))

  # ---- Optional covariate adjustment (limma::removeBatchEffect) -----------------
  # Regress the chosen covariates out of each omics block while preserving the
  # analysis variable (design = ~outcome). Numeric covariates and factor covariates
  # (via model.matrix dummies) are both supported. Covariate values come from the
  # summary metadata, falling back to the full donor table (proc_metadata) so any
  # platform variable the user picks can be resolved.
  cov_mm <- NULL; cov_des <- NULL
  cov_sel <- if (is.null(covariates) || covariates %in% c("", "NA", "null", "NULL")) character(0)
             else trimws(strsplit(covariates, ",")[[1]])
  cov_src <- meta_full
  if (length(cov_sel) > 0) {
    need <- setdiff(cov_sel, colnames(cov_src))
    if (length(need) > 0) {
      library(RSQLite)
      meta_db <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
      pheno <- dbReadTable(meta_db, "proc_metadata")
      dbDisconnect(meta_db)
      have <- intersect(need, colnames(pheno))
      if (length(have) > 0) cov_src <- merge(cov_src, pheno[, c("record_id", have), drop = FALSE], by = "record_id", all.x = TRUE)
    }
  }
  cov_sel <- intersect(cov_sel, colnames(cov_src))
  if (length(cov_sel) > 0) {
    library(limma)
    cov_df <- cov_src[match(shared, cov_src$record_id), cov_sel, drop = FALSE]
    for (cc in cov_sel) {
      v <- cov_df[[cc]]
      if (is.numeric(v)) {
        v[is.na(v)] <- median(v, na.rm = TRUE)
      } else {
        v <- as.character(v); v[is.na(v) | v == ""] <- "NA"
        v <- factor(v)
      }
      cov_df[[cc]] <- v
    }
    cform <- as.formula(paste("~", paste(sprintf("`%s`", cov_sel), collapse = " + ")))
    cov_mm <- model.matrix(cform, data = cov_df)
    cov_mm <- cov_mm[, colnames(cov_mm) != "(Intercept)", drop = FALSE]
    yvec   <- if (primaryType == "disc") factor(Y) else as.numeric(Y)
    cov_des <- model.matrix(~ yvec)
    print(paste("Adjusting for covariates:", paste(cov_sel, collapse = ", ")))
  }

  proc_variable <- read.csv(paste0(other.tables.path, "display_interface/proc_variable_multiomics.csv"))
  out_anno <- qs::qread(paste0(other.tables.path,"libraries/gene_anno_full.qs")) 
  
  # Outcome vector for the SUPERVISED feature screen inside the loop below.
  # Discrete -> factor (2-group t / multi-group F); continuous -> numeric (regression t = correlation).
  yvec_screen <- if (primaryType == "disc") droplevels(factor(Y)) else as.numeric(Y)

# Per-omics SIGNED limma t of the outcome contrast, over the FULL block (all features), for a
# proper GSEA ranking. Collected via <<- so kpList (the DIABLO input) is left byte-for-byte
# unchanged. Covariate-adjusted when covariates are selected, marginal otherwise -- so the GSEA
# ranking follows the same covariate choice as the rest of the run.
screen_tstat <- list()
kpList <- lapply(names(dataList), function(omics_name) {

  x <- dataList[[omics_name]][, shared]
 
  n_features <- nrow(x)
  
  # Adaptive threshold based on omics size
  cap <- if (grepl("rnaseq|pbrna", omics_name)) 5000 
         else if (grepl("prot", omics_name)) 3000 
         else if (grepl("metabo", omics_name)) 500 
         else if (grepl("contaminants", omics_name)) 200 
         else 1000
  
  top_k <- min(n_features, cap)

  # SUPERVISED screen: keep the top_k features most ASSOCIATED with the outcome.
  # limma moderated t covers every outcome type in one path: |t| of the group term
  # (2 groups), max |t| across group terms (>2 groups), or |t| of the slope for a
  # CONTINUOUS outcome (a robust correlation). This replaces the old unsupervised MAD
  # filter, which ranked by raw variability and so discarded low-variance but
  # discriminative features (the real biomarkers) before DIABLO ever saw them.
  xm <- as.matrix(x)
  score <- tryCatch({
    fit <- limma::eBayes(limma::lmFit(xm, model.matrix(~ yvec_screen)))
    yc  <- grep("yvec_screen", colnames(fit$t))
    if (length(yc) == 0) yc <- ncol(fit$t)
    sc  <- matrixStats::rowMaxs(abs(as.matrix(fit$t[, yc, drop = FALSE])))
    sc[is.na(sc)] <- 0
    sc
  }, error = function(e) {
    print(paste("supervised screen failed for", omics_name, "-", conditionMessage(e), "; MAD fallback"))
    matrixStats::rowMads(xm, na.rm = TRUE)
  })
  # SIGNED limma t of the outcome contrast over the FULL block (for GSEA ranking). Uses the
  # covariate-adjusted design when covariates are selected (cbind(outcome, covariates)), else the
  # marginal outcome design. Independent of the selection above, so DIABLO's feature set is unchanged.
  screen_tstat[[omics_name]] <<- tryCatch({
    des <- if (!is.null(cov_mm)) cbind(cov_des, cov_mm) else model.matrix(~ yvec_screen)
    ft  <- limma::eBayes(limma::lmFit(xm, des))
    yc2 <- grep("yvec", colnames(ft$t)); yc2 <- if (length(yc2) >= 1) yc2[1] else 2L
    st  <- as.numeric(ft$t[, yc2]); names(st) <- rownames(xm); st[is.na(st)] <- 0; st
  }, error = function(e) { st <- rep(0, nrow(xm)); names(st) <- rownames(xm); st })
  idx <- order(score, decreasing = TRUE)[seq_len(top_k)]
  x <- x[idx, ]
  # Regress out the chosen covariates (keeps the outcome) before scaling, if requested.
  # Guarded: a collinear/degenerate covariate just leaves that block unadjusted.
  if (!is.null(cov_mm)) {
    x <- tryCatch(
      limma::removeBatchEffect(as.matrix(x), covariates = cov_mm, design = cov_des),
      error = function(e) { print(paste("covariate adjustment skipped for", omics_name, ":", conditionMessage(e))); x }
    )
  }
  filtX <- as.matrix(t(scale(t(x), center=TRUE, scale=TRUE)))
  
  return(filtX)
})

names(kpList) <- names(dataList)

   qs::qsave(kpList,"omics_kpList.qs")
   qs::qsave(screen_tstat, "screen_tstat.qs")   # signed limma t per omics (full block) for GSEA ranking
 
  library(mixOmics)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  
  print(c(primaryType,"primaryType"))
  
  nm1 <- paste0("total_var_1_",fig_counts,".png")
  nm2 <- paste0("factor_var_1_",fig_counts,".png")
  nm3 <- paste0("sel_feature_1_",fig_counts,".png")
  nm4 <- paste0("sample_view_1_",fig_counts,".png")
  
  # Create display names for each omics (use view2omics for nice names, fallback to original)
  display_names <- sapply(selected_omics, function(x) {
    nice_name <- view2omics(x)  # Returns factor with nice names
    if(is.na(nice_name) || as.character(nice_name) == "other") {
      return(x)  # Fallback to original name if not in mapping
    } else {
      return(as.character(nice_name))
    }
  })
 
 

  # CRITICAL: mixOmics expects samples as ROWS, features as COLUMNS
  # Current format: features x samples, need to transpose to: samples x features
  X_diablo <- lapply(kpList, t)
  names(X_diablo) <- sapply(names(X_diablo), function(x) as.character(view2omics(x)))
 

  n_comp <- min(3, length(Y) - 1, min(sapply(kpList, nrow)))   
 
  # design = block<->block correlation weight (mixOmics auto-links each block to Y at 1).
  # 0.2 = light-to-moderate cross-omics integration: leans toward outcome discrimination but keeps a
  # cross-block link. Raise toward 1 for stronger integration, lower toward 0 to decouple the blocks.
  design <- matrix(0.2, nrow = length(X_diablo), ncol = length(X_diablo), dimnames = list(names(X_diablo), names(X_diablo)))
  diag(design) <- 0
  
  print("Running DIABLO with fixed-keepX sparse feature selection...")
 
 
    #test.keepX <- lapply(X_diablo, function(x) {
   #p <- ncol(x)
   # if (p < 200)    return(c(5, 10, 20))
   # if (p < 1000)   return(c(10, 30, 60))
  # if (p < 5000)   return(c(30, 60, 120))
  # return(c(50, 100, 200))})
 
  # Cross-validation to find optimal number of features
  #print("Running cross-validation (this may take a few minutes)...")
  
  # SET SEED for reproducibility
  set.seed(12345)
  

if (primaryType == "disc") {
    if(!all(sapply(kpList, function(x) identical(colnames(x), names(Y))))) {
    warning("Sample names don't match between Y and X_list, aligning now...")
    sample_order <- names(Y)
    kpList <- lapply(kpList, function(x) x[, sample_order, drop = FALSE])
  }
      
  #tune.result <- tune.block.splsda(  X = X_diablo, Y = Y, ncomp = n_comp, test.keepX = test.keepX, design = design, nrepeat = 2, folds = 3, progressBar = TRUE)
  
    #keepX_optimal <- tune.result$choice.keepX
  #deep_ncomp <- min(10, length(Y) - 1)
}  
 
  # ---- Sparse biomarker selection (fixed keepX, no cross-validation) ----------
  # DIABLO selects features via an L1 penalty: keep the top `keepX_per_comp`
  # features per omics block PER component (non-selected features get a zero
  # loading). Adjust this single number to make the signature larger/smaller.
  # (CV tuning via tune.block.splsda is intentionally left off for speed.)
  # keepX = how many features DIABLO selects per block per component. Kept small (sparse) per
  # mixOmics DIABLO practice: the canonical case-study tuning grids span ~5-25 per component.
  # 25 is a fixed global default (sparse, interpretable) — the site can't CV-tune keepX per request,
  # and a one-off tune.block.splsda for State 2 vs 3 put the comp-1 optimum at ~20 prot / 5 rna
  # (BER 0.23), so 25 is a sound one-size default. The "Top features" slider mirrors it.
  keepX_per_comp <- 25
  if(primaryType == "disc") {
    kx <- lapply(X_diablo, function(x) rep(min(keepX_per_comp, ncol(x)), 3))
    diablo_model <- block.splsda(X = X_diablo, Y = Y, ncomp = 3, keepX = kx, design = design)
  } else {
    kx <- lapply(X_diablo, function(x) rep(min(keepX_per_comp, ncol(x)), 5))
    diablo_model <- block.spls(X = X_diablo, Y = Y, ncomp = 5, keepX = kx, design = design, mode = "regression")
  }
  # Persist the fitted model so the circos similarity can be (re)built quickly
  # from it (e.g. when the user changes the "top features" slider) without
  # re-running the expensive integration.
  qs::qsave(diablo_model, "diablo_model.qs")
   # Extract variance from prop_expl_var which contains per-omics variance
  # Structure: prop_expl_var$proteomics, prop_expl_var$transcriptomics (NOT $X!)
  if(!is.null(diablo_model$prop_expl_var)) {
    # Get omics names (exclude 'Y' which is the response)
    omics_names <- setdiff(names(diablo_model$prop_expl_var), "Y")
    
    total_var <- sapply(omics_names, function(nm) {
      ev <- diablo_model$prop_expl_var[[nm]]
      if(is.null(ev)) return(0)
      sum(as.numeric(ev), na.rm = TRUE) * 100  # Convert to percentage
    })
  } else {
    total_var <- numeric(0)
  }
  
  print("Total variance extracted:")
  print(total_var)
  
  # Handle empty variance
  use_feature_counts <- FALSE
  if(length(total_var) == 0 || all(is.na(total_var)) || all(total_var == 0)) {
    warning("No valid variance extracted from DIABLO model. Using feature count as proxy.")
    use_feature_counts <- TRUE
    # Use number of selected features as proxy for "importance"
    total_var <- sapply(names(kpList), function(nm) {
      # Find matching block in diablo_model
      block_name <- names(diablo_model$loadings)[which(names(kpList) == nm)]
      if(length(block_name) == 0) return(10)  # Fallback
      sum(diablo_model$loadings[[block_name]] != 0)
    })
    names(total_var) <- names(kpList)
    print("Using feature counts as proxy:")
    print(total_var)
  }
  
  # Plot 1: Variance scree plot (line plot per component)
  if(!is.null(diablo_model$prop_expl_var)) {
    omics_names <- setdiff(names(diablo_model$prop_expl_var), "Y")
    
    scree_df <- data.frame()
    for(omics_name in omics_names) {
      var_vec <- diablo_model$prop_expl_var[[omics_name]]
      if(!is.null(var_vec) && length(var_vec) > 0) {
        for(i in 1:min(length(var_vec), n_comp)) {
          scree_df <- rbind(scree_df, data.frame(
            Component = i,
            omics = omics_name,
            variance = as.numeric(var_vec[i]) * 100,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    if(nrow(scree_df) > 0) {
      omics_colors <- sapply(unique(scree_df$omics), get_omics_color)
      names(omics_colors) <- unique(scree_df$omics)
      
      p1 <- ggplot(scree_df, aes(x = Component, y = variance, color = omics, group = omics)) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 3) +
        scale_color_manual(values = omics_colors, name = "Omics") +
        scale_x_continuous(breaks = 1:n_comp) +
        labs(x = "Component #", y = "Var. (%)") +
        theme_minimal(base_size = 10, base_family = "Arial") +
        theme(legend.position = "top", legend.title = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color = "black"),
              axis.title = element_text(color = "black"))
      ggsave(filename = nm1, plot = p1, width = 760, height = 820, units = "px", dpi = 150)
    } else {
      png(nm1, width = 1000, height = 600, units = "px", res = 150)
      plot(1, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), axes=FALSE)
      text(0.5, 0.5, "No variance data", cex=1.5)
      dev.off()
    }
  } else {
    png(nm1, width = 582, height = 580, units = "px", res = 200)
    plot(1, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), axes=FALSE)
    text(0.5, 0.5, "No variance data", cex=1.5)
    dev.off()
  }
  
  # Plot 2: AUROC — DISABLED. The AUROC is in-sample (fit + evaluated on the same donors,
  # on features DIABLO selected from those donors), so it is a training-set mirage and is
  # NOT displayed on the interface. Skip the expensive mixOmics auroc()/perf() call and just
  # write a tiny placeholder so the nm2 image slot in the RES-OK response stays valid.
  png(nm2, width = 8, height = 8, units = "px", res = 72); par(mar = c(0,0,0,0)); plot.new(); dev.off()

  # Extract selected features from loadings
  # ==================================================
  # IMPORTANT: Understanding DIABLO Feature Selection
  # ==================================================
  # 1. CV (tune.block.splsda) determines OPTIMAL NUMBER of features per component
  #    Example: Component1 = 10 features, Component2 = 50 features
  # 2. DIABLO (block.splsda) then selects WHICH specific features using keepX
  #    It selects features with highest discriminative power
  # 3. Loadings matrix: features × components
  #    Non-zero loading = feature was SELECTED by DIABLO for that component
  #    Zero loading = feature was NOT selected
  # 4. Loading magnitude = feature importance (higher = more important)
  # 5. We combine features from ALL components for network analysis
  # ==================================================
  
  print("Extracting sparse-selected features from DIABLO model...")
  print("Note: features with non-zero loadings are those DIABLO selected (fixed keepX per component)")
  feature_list <- list()
  
  for(omics_nm in display_names) {
    loadings_mat <- diablo_model$loadings[[omics_nm]]
    # Guard: a selected omics can be absent from the fitted model (dropped when it loaded empty,
    # lost all common donors, or is unwired like methylation) -> loadings_mat is NULL. Skip it
    # instead of crashing on loadings_mat[, comp] ("incorrect number of dimensions").
    if (is.null(loadings_mat)) { print(paste("Skipping", omics_nm, "- not in fitted model")); next }

    print(paste("Loadings for", omics_nm, ":", nrow(loadings_mat), "features x", ncol(loadings_mat), "components"))
    validated_ncomp <- min(3,ncol(loadings_mat)) 
    # Get all non-zero features across ALL components
    features_data <- data.frame()
    for(comp in 1:validated_ncomp) {
      comp_loadings <- loadings_mat[, comp]
      selected_idx <- comp_loadings != 0
    
      if(any(selected_idx)) {
        features_data <- rbind(features_data, data.frame(
          feature = rownames(loadings_mat)[selected_idx],
          loading = abs(comp_loadings[selected_idx]),
          loading_signed = as.numeric(comp_loadings[selected_idx]),  # keep the sign (direction)
          component = comp,
          view = omics_nm,
          stringsAsFactors = FALSE
        ))
      }
    }
 
    
    # Aggregate by feature: keep the component where the feature's |loading| is
    # largest, along with that loading. Keeping the component lets the UI show
    # which DIABLO component each feature's importance comes from (loadings are
    # only directly comparable within the same omics + component).
    if(nrow(features_data) > 0) {
      features_agg <- features_data %>%
        group_by(feature, view) %>%
        slice_max(order_by = loading, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        rename(mean_abs_loading = loading)
      feature_list[[omics_nm]] <- features_agg
    }
  }
  
  all_features <- bind_rows(feature_list)

  # Annotate ALL features for client-side JSON
  all_features$label <- as.character(all_features$feature)
  all_features$omics <- as.character(view2omics(all_features$view))
  all_features$label[all_features$label %in% out_anno$ensembl_gene_id] <- out_anno$hgnc_symbol[match(all_features$label[all_features$label %in% out_anno$ensembl_gene_id], out_anno$ensembl_gene_id)]
  all_features$label[all_features$label %in% out_anno$entrez_selected] <- out_anno$hgnc_symbol[match(all_features$label[all_features$label %in% out_anno$entrez_selected], out_anno$entrez_selected)]
  all_features$color <- sapply(all_features$omics, get_omics_color)
  # For a 2-class discriminant run, component 1 IS the group-vs-group separating axis;
  # components 2/3 capture variation orthogonal to that split. Rank component 1 first so the
  # real discriminants top the list instead of high-loading comp-2/3 features. Continuous
  # (block.spls) runs keep the plain loading-magnitude order (all components are meaningful).
  if (primaryType == "disc") {
    all_features <- all_features %>% arrange(component, desc(mean_abs_loading))
  } else {
    all_features <- all_features %>% arrange(desc(mean_abs_loading))
  }
  # Remove duplicates (keep highest loading per feature)
  all_features <- all_features[!duplicated(all_features$feature), ]

  qs::qsave(all_features, "all_features.qs")
  jsonlite::write_json(all_features[, c("feature", "label", "omics", "mean_abs_loading", "loading_signed", "component", "color")], "all_features.json")

  top_features <- all_features %>% slice_head(n = topFeat)
  top_features <- top_features %>% mutate(label = factor(label, levels = rev(unique(label))))

  # Create dynamic color palette
  feature_colors <- sapply(unique(top_features$omics), get_omics_color)
  names(feature_colors) <- unique(top_features$omics)
  
  p3 <- ggplot(top_features, aes(x = mean_abs_loading, y = label, colour = omics)) +
    geom_segment(aes(x = 0, xend = mean_abs_loading, yend = label), linewidth = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = feature_colors, name = "Omics") +
    labs(title = " ", x = "Loading", y = NULL) +
    coord_cartesian(xlim = c(0, max(top_features$mean_abs_loading, na.rm = TRUE) * 1.08)) +
    theme_minimal(base_size = 9, base_family = "Arial") +
    theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"), 
          axis.title.x = element_text(colour = "black"), plot.title = element_text(face = "bold"), legend.position = "right")
  n_features <- nrow(top_features)
  ggsave(nm3, p3, width = 6, height = max(4, n_features * 0.12), dpi = 300)
  
  # Plot 4: Arrow plot (ggplot2 version)
  # Extract data for arrow plot

  if(primaryType == "disc") {

 plot4 <- plotArrow(diablo_model,
                  block = 'average',
                  ind.names = FALSE,    # hide per-donor IDs (unreadable with hundreds of samples)
                  ellipse = TRUE,
                  size.legend = 0.6,    # Slightly larger than 0.5 but still small
                  legend = TRUE,
                  title = '',
                  abline = TRUE)        # Adds guide lines for better context
 
 # Apply custom metadata colors (distinct from omics colors)
 n_groups <- length(unique(Y))
 group_colors <- get_metadata_colors(n = n_groups, continuous = FALSE)
 names(group_colors) <- levels(Y)
 
 plot4 <- plot4 + 
   scale_fill_manual(values = group_colors) +
   scale_color_manual(values = group_colors)

}else{
# 1. Identify indices and names dynamically
y_idx <- diablo_model$indY 
all_block_names <- names(diablo_model$X)
predictor_blocks <- all_block_names[-y_idx] 

# 2. Extract numeric outcome for coloring
outcome_values <- as.numeric(diablo_model$X[[y_idx]])

# 3. Create a temporary model that ONLY contains predictor blocks
# This prevents plotArrow from ever "seeing" the Y node
temp_model <- diablo_model
temp_model$X <- diablo_model$X[-y_idx]
temp_model$variates <- diablo_model$variates[-y_idx]


plot4 <- plotArrow(temp_model,
                   block = predictor_blocks,
                   ind.names = FALSE,
                   group = outcome_values,
                   legend = F,
                   title = '',
                   abline = TRUE)
# 4. FIX THE CONTINUOUS ERROR (Update internal layers to numeric)
plot4$data$group <- as.numeric(as.character(plot4$data$group))
for(i in seq_along(plot4$layers)) {
  if(!is.null(plot4$layers[[i]]$data)) {
    plot4$layers[[i]]$data$group <- as.numeric(as.character(plot4$layers[[i]]$data$group))
  }
}

# 5. FULLY DYNAMIC SHAPE MAPPING (Zero Hard-Coding)
# Assigns 8 to 'centroid', then 1, 2, 3... to whatever blocks are in predictor_blocks
dynamic_shapes <- c(8, seq_along(predictor_blocks))
names(dynamic_shapes) <- c("centroid", predictor_blocks)

plot4 <- plot4 + 
  scale_shape_manual(
    name = "Omics Type",
    values = dynamic_shapes,
    breaks = c("centroid", predictor_blocks) # Programmatically removes 'Y'
  )



plot4 <- plot4 + 
  scale_color_gradientn(colors = get_metadata_colors(continuous = TRUE)) +
  labs(color = all_block_names[y_idx], x = "Dimension 1", y = "Dimension 2") + 
  theme_classic() + 
  theme(
    legend.position = "right",
    legend.box = "vertical",
    axis.line = element_line(color = "black") # Only X and Y axes
  )

 
}


# 2. Increase dimensions to accommodate 300 DPI
# For 300 DPI, a width of 1800-2400 is standard for a clear publication plot
Cairo::Cairo(file = nm4,
             unit = "px",
             res = 300,
             width = 1850,   # wider (landscape) so the sample view fills its wider card
             height = 1250,
             type = "png",
             bg = "white")

# 3. Print the object
print(plot4)

# 4. Close device
dev.off()

  msg <- ""
  print(c("RES-OK;",nm1,nm2,nm3,nm4,msg))
  qs::qsave(top_features,"top_features.qs")

  # Standard mixOmics DIABLO cross-omics similarity circos, built from the
  # fitted model in this same run (fast, no bootstrap). Guarded so a circos
  # failure can never break the main integration result.
  circos_nm <- tryCatch(
    build_diablo_circos(diablo_model, top_features, fig_counts),
    error = function(e) { print(paste("circos build failed:", conditionMessage(e))); "" }
  )

  # Per-component, per-omics signed loadings for the interactive loading plots.
  loadings_nm <- tryCatch(
    build_diablo_loadings(diablo_model, all_features, primaryType, Y, fig_counts),
    error = function(e) { print(paste("loadings build failed:", conditionMessage(e))); "" }
  )

  # Generate dynamic message for HTML display
  # Count features per omics
  feature_counts <- table(top_features$view)
  feature_msg <- paste(sapply(names(feature_counts), function(nm) {
    paste0(feature_counts[nm], " ", nm)
  }), collapse=" and ")
  
  # Build complete message
  if(primaryType == "disc") {
    method_msg <- sprintf(
      "Omics datasets were integrated using DIABLO (sparse block PLS-DA). %s features were selected across %d components (keepX = %d per block per component), optimised for discriminating %s.",
      feature_msg,
      n_comp,
      keepX_per_comp,
      analysisVar
    )
  } else {
    method_msg <- sprintf(
      "Omics datasets were integrated using DIABLO (sparse block PLS). %s features were selected across %d components (keepX = %d per block per component), optimised for predicting %s.",
      feature_msg,
      n_comp,
      keepX_per_comp,
      analysisVar
    )
  }
  
  gc()
  return(paste0("RES-OK;",nm1,";",nm2,";",nm3,";",nm4,";",method_msg,";",circos_nm,";",loadings_nm))
}

#########################
######################################################
##############################################################

# ------------------------------------------------------------------
# Standard mixOmics DIABLO cross-omics similarity for the circos plot.
#
# Similarity(feature i in block h, feature j in block k) =
#     sum_c  cor(X_h[, i], variate_h[, c]) * cor(X_k[, j], variate_k[, c])
#
# This is the construction mixOmics network()/circosPlot use: deterministic
# and fast (NO bootstrap / permutation). Cross-omics (between-block) pairs
# only -- the multi-omics integration story. One similarity per pair.
# Returns the JSON filename {nodes, edges}, or "" on failure.
# ------------------------------------------------------------------
build_diablo_circos <- function(model, top_features, fig_counts) {
  library(jsonlite)
  feats <- as.character(top_features$feature)

  # Combined matrix of the selected features across omics blocks (samples x features)
  # + a feature -> omics-block lookup. model$X is the centred/scaled input, so column
  # correlations are the Pearson correlations between features across donors.
  mats <- list(); ftype <- character(0)
  for (b in names(model$X)) {
    Xb <- model$X[[b]]
    if (is.null(Xb)) next
    sel <- intersect(colnames(Xb), feats)        # skips the Y block automatically
    if (length(sel) == 0) next
    mats[[b]] <- Xb[, sel, drop = FALSE]
    ftype[sel] <- b
  }

  nodes <- data.frame(
    id    = as.character(top_features$feature),
    type  = as.character(top_features$omics),
    label = as.character(top_features$label),
    freq  = as.numeric(top_features$mean_abs_loading),
    color = as.character(top_features$color),   # shared omics palette (get_omics_color)
    stringsAsFactors = FALSE
  )

  if (length(mats) == 0) {
    netnm <- paste0("circos_", fig_counts, ".json")
    write(toJSON(list(nodes = nodes, edges = data.frame(), n = 0), auto_unbox = TRUE, pretty = TRUE), file = netnm)
    return(netnm)
  }

  Xcomb <- do.call(cbind, mats)
  n <- nrow(Xcomb)

  # Pairwise Pearson correlation + two-sided p-value (t) + BH FDR for every unique
  # feature pair. Cross-omics (intra = FALSE) and within-omics (intra = TRUE).
  R <- suppressWarnings(stats::cor(Xcomb, use = "pairwise.complete.obs"))
  R[is.na(R)] <- 0
  idx <- which(upper.tri(R), arr.ind = TRUE)
  fa <- rownames(R)[idx[, 1]]; fb <- colnames(R)[idx[, 2]]
  r  <- R[idx]
  rc <- pmin(pmax(r, -0.999999), 0.999999)
  tval <- rc * sqrt((n - 2) / (1 - rc^2))
  pval <- 2 * stats::pt(-abs(tval), df = n - 2)
  padj <- stats::p.adjust(pval, method = "BH")

  edges <- data.frame(
    source = fa, target = fb,
    weight = round(r, 4),
    pval   = signif(pval, 3),
    padj   = signif(padj, 3),
    intra  = unname(ftype[fa] == ftype[fb]),
    sign   = ifelse(r >= 0, "pos", "neg"),
    stringsAsFactors = FALSE
  )

  # Keep edges_A_B.qs (cross-omics from/to) so the pathway "network features only" option still works.
  cross <- edges[!edges$intra, , drop = FALSE]
  qs::qsave(data.frame(from = cross$source, to = cross$target, stringsAsFactors = FALSE), "edges_A_B.qs")

  netnm <- paste0("circos_", fig_counts, ".json")
  write(toJSON(list(nodes = nodes, edges = edges, n = n), auto_unbox = TRUE, pretty = TRUE), file = netnm)
  netnm
}

# ------------------------------------------------------------------
# Per-component, per-omics signed loadings (long format) for the interactive
# loading plots. One row per (feature, component) the feature was selected on,
# with the signed loading and, for a discriminant outcome, the group in which
# that feature is highest (median). Returns the JSON filename, or "".
# ------------------------------------------------------------------
build_diablo_loadings <- function(model, all_features, primaryType, Y, fig_counts) {
  library(jsonlite)
  lab <- setNames(as.character(all_features$label), as.character(all_features$feature))
  grpLevels <- if (primaryType == "disc" && !is.null(Y)) levels(factor(Y)) else character(0)
  yf <- if (length(grpLevels) > 0) factor(Y) else NULL

  blocks <- setdiff(names(model$loadings), "Y")
  rows <- list()
  for (b in blocks) {
    L <- model$loadings[[b]]
    if (is.null(L) || is.null(nrow(L))) next
    Xb <- model$X[[b]]
    om <- as.character(view2omics(b))
    for (cmp in seq_len(ncol(L))) {
      v <- L[, cmp]
      sel <- which(v != 0)
      for (i in sel) {
        f <- rownames(L)[i]
        grp <- ""
        if (!is.null(yf) && !is.null(Xb) && f %in% colnames(Xb)) {
          meds <- tapply(Xb[, f], yf, median, na.rm = TRUE)
          if (length(meds)) grp <- names(which.max(meds))
        }
        rows[[length(rows) + 1]] <- data.frame(
          feature = f,
          label   = if (!is.na(lab[f])) lab[f] else f,
          omics   = om,
          component = cmp,
          loading = round(as.numeric(v[i]), 4),
          group   = grp,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0) return("")
  df <- do.call(rbind, rows)
  out <- list(primaryType = primaryType, groups = grpLevels,
              ncomp = max(df$component), loadings = df)
  netnm <- paste0("loadings_", fig_counts, ".json")
  write(toJSON(out, auto_unbox = TRUE, pretty = TRUE), file = netnm)
  netnm
}

# Rebuild the circos similarity from the saved DIABLO model -- used when the
# user changes the "top features" slider after a run. Fast (no bootstrap):
# the old slow differential-correlation test was replaced by the standard
# mixOmics model-based similarity in build_diablo_circos() above.
GetCorrNet <- function(
    omics1, omics2, omics3,
    varGroup,
    analysisVar,
    class1 = NULL,
    class2 = NULL,
    topFeat = 100 ){

  fig_counts <- qs::qread("fig_counts.qs")

  if (!file.exists("diablo_model.qs") || !file.exists("all_features.qs")) {
    return("NO;No integration model found. Please run the integration first.")
  }
  diablo_model <- qs::qread("diablo_model.qs")
  all_features <- qs::qread("all_features.qs")

  # Re-slice all_features to the requested top N and refresh top_features.qs
  topFeat <- max(10, min(as.numeric(topFeat), nrow(all_features)))
  all_features <- all_features[order(-all_features$mean_abs_loading), ]
  top_features <- head(all_features, topFeat)
  qs::qsave(top_features, "top_features.qs")

  netnm <- tryCatch(
    build_diablo_circos(diablo_model, top_features, fig_counts),
    error = function(e) { print(paste("circos build failed:", conditionMessage(e))); "" }
  )
  if (identical(netnm, "") || !file.exists(netnm)) {
    return("NO;Failed to build the cross-omics similarity.")
  }

  net <- jsonlite::fromJSON(netnm)
  n_correlations <- if (is.null(net$edges)) 0 else nrow(net$edges)
  n_features     <- if (is.null(net$nodes)) 0 else nrow(net$nodes)
  omics_counts   <- if (n_features > 0) table(net$nodes$type) else integer(0)
  omics_summary  <- paste(sapply(names(omics_counts), function(nm) paste0(omics_counts[nm], " ", nm)),
                          collapse = " and ")
  net_msg <- sprintf(
    "Cross-omics correlation network among the DIABLO-selected features: %d pairs among %d features (%s).",
    n_correlations, n_features, omics_summary
  )
  return(paste0("RES-OK;", netnm, ";", net_msg))
}
 

# Wrapper: return any R error as a readable string so it reaches the web interface.
# Without this, an R error propagated to RCenter.getMediation, which catches it and
# returns null; JAX-RS turns a null TEXT_PLAIN body into HTTP 204, and the frontend
# reported every distinct failure as the same opaque "Failed to process!". Mirrors
# the GetIntegPlot / GetIntegPlot_impl pattern above. Real work is in the _impl below.
GetMediation <- function(...) {
  tryCatch(GetMediation_impl(...), error = function(e) paste0("Error: Mediation - ", conditionMessage(e)))
}

GetMediation_impl <- function(
     predictor, mediator, predictorType, mediatorType,
    varGroup,
    analysisVar,
    class1 = NULL,
    class2 = NULL,primaryType ){

print(c( predictor, mediator, predictorType, mediatorType,
    varGroup,  
    analysisVar,class1,class2,primaryType))
 

kpList <- qs::qread( "omics_kpList.qs")


 meta<-read.csv(paste0(other.tables.path,"display_data/metadata_sum_norm.csv"))
  # Inject cluster (or any proc_metadata-only) column if selected var is missing from the CSV
  if (!is.null(analysisVar) && !(analysisVar %in% colnames(meta))) {
    library(RSQLite)
    meta_db <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
    pheno <- dbReadTable(meta_db, "proc_metadata")
    dbDisconnect(meta_db)
    if (analysisVar %in% colnames(pheno)) {
      extra <- pheno[, c("record_id", analysisVar)]
      meta <- merge(meta, extra, by = "record_id", all.x = TRUE)
    } else {
      stop(paste0("Analysis failed: variable '", analysisVar,
                  "' not found in metadata or donor table."))
    }
  }

   types <- omics2view(c(predictorType,mediatorType))
  
   outcome <- analysisVar
   predictor_matrix <- kpList[[types[1]]][rownames(kpList[[types[1]]]) %in% predictor, ,drop=F]
   mediator_matrix <- kpList[[types[2]]][rownames(kpList[[types[2]]]) %in% mediator, ,drop=F]
   
   input_data <- as.data.frame(t(rbind(as.matrix(predictor_matrix), as.matrix(mediator_matrix))))
   input_data$outcome <-  meta[[outcome]][match(rownames(input_data),meta$record_id)]
   input_data <- input_data[input_data$outcome %in% c(class1,class2),]
   
   if( primaryType=="disc"){
     input_data$outcome <- as.character(input_data$outcome)
     input_data$outcome[input_data$outcome==class1] <- 0
     input_data$outcome[input_data$outcome==class2] <- 1
     input_data$outcome<-as.factor(as.numeric(input_data$outcome))
   }else{
     input_data$outcome <- as.numeric(input_data$outcome)
   }
   result <- PerformCausalMediation(input_data, primaryType,)
   result$Effect <- gsub("control",class1, result$Effect)
   result$Effect <- gsub("treated",class2, result$Effect)
   
   effect_data <- result[grepl("ACME|ADE|Total", result$Effect), ]
   effect_data$Type <- factor(ifelse(grepl("ACME", effect_data$Effect), "Indirect\n(ACME)",
                                     ifelse(grepl("ADE", effect_data$Effect), "Direct\n(ADE)", "Total")),
                              levels = c("Total", "Direct\n(ADE)", "Indirect\n(ACME)"))
   effect_data$Group <- factor(ifelse(grepl(class1, effect_data$Effect),class1,
                                      ifelse(grepl(class2, effect_data$Effect), class2,
                                             ifelse(grepl("average", effect_data$Effect), "Average", "Combined"))),
                               levels = c(class1, class2, "Average", "Combined"))
   
   
   
   prop_data <- result[grepl("Prop. Mediated", result$Effect), ]
   prop_data$Group <- factor(ifelse(grepl(class1, prop_data$Effect), class1,
                                    ifelse(grepl(class2, prop_data$Effect), class2, "Average")),
                             levels = c(class1, class2, "Average"))
   
   # Color palette
   # Use custom metadata colors for mediation groups
   meta_colors <- get_metadata_colors(n = 4, continuous = FALSE)
   effect_colors <- c(class1 = meta_colors[1], 
                      class2 = meta_colors[2], 
                      "Average" = meta_colors[3],
                      "Combined" = meta_colors[4])
   names(effect_colors) = c(class1, class2, "Average", "Combined")
   # Horizontal Effect Plot
 
   
   library(ggplot2)  
   library(scales)
   p_effects <- ggplot(effect_data, aes(y = Type, x = Estimate, color = Group)) +
     geom_point(position = position_dodge(width = 0.6), size = 3) +
     geom_errorbarh(aes(xmin = CI.lower, xmax = CI.upper), 
                    height = 0.2, position = position_dodge(width = 0.6)) +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
     labs(y = "", x = "Effect Size") +
     scale_color_manual(values = effect_colors) +
     theme_minimal(base_size = 14) +
     theme(legend.position = "bottom",
           panel.grid.minor = element_blank(),
           axis.text.y = element_text(face = "bold"),
           legend.title = element_blank()) +
     geom_text(aes(label = ifelse(!is.na(p.value) & p.value < 0.05, "*", ""), 
                   group = Group),
               position = position_dodge(width = 0.6), 
               hjust = -0.3, size = 7, color = "black") +
     scale_x_continuous(limits = c(min(effect_data$CI.lower, na.rm = TRUE) - 0.001,
                                   max(effect_data$CI.upper, na.rm = TRUE) + 0.001))
   
   # Feature names go into the PLOT FILENAME, and the file is fetched back over
   # /download, where Download.java does `fileNm.replace("-", "/")` to allow
   # subdirectory navigation. A hyphen in a feature id therefore became a path
   # separator: "MT-CO1_SCG5_med.png" was looked up as "MT/CO1_SCG5_med.png",
   # 404'd, and the user saw a failure even though the analysis had SUCCEEDED.
   # Affects hyphenated symbols (MT-CO1, HLA-DRB1, NKX2-2) and metabolites
   # (L-Tyrosine, L-Carnitine), but not ENSG*/MAR* ids -- hence the apparent
   # randomness. Sanitise here rather than in Download.java, whose '-' rule other
   # pages rely on. Anything outside [A-Za-z0-9_.] becomes '_'; the sanitised name
   # is what we ggsave AND what we return, so the two always agree.
   safe_nm <- function(x) gsub("[^A-Za-z0-9_.]", "_", as.character(x))
   nm1 = paste0(safe_nm(predictor),"_",safe_nm(mediator),"_","med.png")

      ggsave(nm1, plot = p_effects, width = 5.5, height = 6, dpi = 300)

   print(paste0("RES-OK;",nm1))
  return(paste0("RES-OK;",nm1))


}

################################################################################
################################################################################
##########################HELPERS###############################################
################################################################################
################################################################################
 


align_factors <- function(Z_ref, Z_run) {
  # Match samples between models
  shared_samples <- intersect(rownames(Z_ref), rownames(Z_run))
  if (length(shared_samples) < 5)
    stop("Too few shared samples to align factors")
  
  Z_ref_sub <- Z_ref[shared_samples, , drop = FALSE]
  Z_run_sub <- Z_run[shared_samples, , drop = FALSE]
  
  # Compute correlation matrix between factors
  cor_mat <- cor(Z_ref_sub, Z_run_sub, use = "pairwise.complete.obs")
  
  # Use the Hungarian algorithm for one-to-one mapping
  assignment <- clue::solve_LSAP(abs(cor_mat), maximum = TRUE)
  mapping <- as.integer(assignment)
  sign_flip <- sign(cor_mat[cbind(seq_len(ncol(Z_ref_sub)), mapping)])
  
  list(mapping = mapping, sign_flip = sign_flip)
}

test_mofa_assoc <- function(df, analysisVar,meta) {
  
  if (!all(c("sample", "factor", "value") %in% names(df))) {
    stop("Data must contain 'sample', 'factor', and 'value' columns.")
  }
  if (!(analysisVar %in% names(df))) {
    stop(paste("Metadata variable", meta_var, "not found in df."))
  }
  
 
  results <- df %>%
    group_by(factor) %>%
    summarise({
      y <- get(analysisVar)
      x <- value
      ok <- complete.cases(x, y)
      if (sum(ok) < 3 || length(unique(y[ok])) < 2) {
        tibble(value = NA, pval = NA, test = NA)
      } else if (is.numeric(y)) {
        n=length(unique(df$sample))
        if (n < 50) {
          method <- "spearman" 
        } else if (n >= 50 && n <= 100 && shapiro.test(meta[[analysisVar]])$p.value >= 0.05) {
          method <- "pearson" 
        } else if (n >= 50 && n <= 100 && shapiro.test(meta[[analysisVar]])$p.value < 0.05) {
          method <- "spearman" 
        } else if (n > 100) {
          method <- "pearson" 
        } else {
          method <- NA 
        }
        
        test <- suppressWarnings(cor.test(x[ok], y[ok],method=method))
        tibble(value = unname(test$estimate), pval = test$p.value, test = "cor")
      } else {
        fit <- lm(x[ok] ~ y[ok])
        a <- anova(fit)
        tibble(value = sqrt(a["y[ok]", "F value"]),
               pval = a["y[ok]", "Pr(>F)"],
               test = "anova")
      }
    }, .groups = "drop") %>%
    mutate(
      covariate = analysisVar,
      padj = p.adjust(pval, "fdr")
    ) %>%
    arrange(padj)
  
  return(results)
}

 get_top_weighted_features <- function(m,var_exp, factors = c("Factor1"),
         n_total = 100,
         min_view_contrib = 0.05) {
 
  # --- 2. Extract loadings ---
  loadings <- get_weights(m, as.data.frame = TRUE)
  
  # --- 3. Container for results ---
  all_top <- list()
  
  for (f in factors) {
    # Get view contributions for this factor
    ve <- var_exp[f, , drop = TRUE]
    ve <- ve / sum(ve)
    
    # Filter views by min contribution
    ve <- ve[ve >= min_view_contrib]
    if (length(ve) <= 1) next
    
    # Compute number of features per view
    n_per_view <- round(ve * n_total)
    names(n_per_view) <- names(ve)
    
    # Subset loadings
    load_f <- subset(loadings, factor == f & view %in% names(ve))%>% mutate(view = as.character(view))
    
    # Select proportionally top |W| per view
    library(dplyr)
    top_feats <- load_f %>%
      mutate(view = as.character(view)) %>%
      group_split(view) %>%
      map_dfr(function(df) {
        view_name <- unique(df$view)
        n <- n_per_view[view_name]
        if (is.na(n) || length(n) == 0) n <- 50
        
        df %>%
          slice_max(order_by = abs(value), n = n, with_ties = FALSE) %>%
          mutate(
            selected_factor = f,
            rank_in_view = row_number()
          )
      }) %>%
      ungroup()
    
    
    
    attr(top_feats, "n_per_view") <- n_per_view
    attr(top_feats, "var_explained") <- ve
    
    all_top[[f]] <- top_feats
  }
  
  if(length(all_top)==0){ 
    return(NULL)
  }
  
  # --- 4. Combine results across factors ---
  combined <- dplyr::bind_rows(all_top, .id = "factor_id")
  attr(combined, "details") <- lapply(all_top, attributes)
  return(combined)
}


plot_mofa_pca <- function(factor_df, top_factors,varGroup,var_df) {
 

  plot_df  <- factor_df %>%
    select(any_of(c("sample", "factor", "value", varGroup))) %>%
    pivot_wider(names_from = factor, values_from = value) %>%
    distinct(sample, .keep_all = TRUE) %>%
    relocate(sample, all_of(varGroup))
  
  v <- plot_df[[varGroup]]
  is_continuous <- is.numeric(v)


  x_lab <- sprintf("%s (%.1f%% var.)", top_factors[1],
                   var_df[top_factors[1],"total"] )
  y_lab <- sprintf("%s (%.1f%% var.)", top_factors[2],
                   var_df[top_factors[2],"total"])
  
  if (is_continuous) {
    color_scale <- scale_color_gradientn(
       colours = c("#60a5fa", "#7cc9d3", "#a0d292"), 
      name = varGroup
    )
  } else {
    color_scale <- scale_color_manual(
    values = c("#60a5fa", "#7cc9d3", "#a0d292"), 
    name = varGroup
  )
  }
  
 
  p <- ggplot(plot_df,
              aes(x = .data[[top_factors[1]]],
                  y = .data[[top_factors[2]]],
                  color = .data[[varGroup]])) +
    geom_point(size = 2, alpha = 1) +
    color_scale +
    theme_minimal(base_size = 8,base_family = "Arial") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      legend.position = "right",
      legend.key.size = unit(0.7, "lines"),
      legend.title = element_text(size = 8),
      legend.text  = element_text(size = 7)
    )+
    labs(
     # title = paste("MOFA sample representation colored by", varGroup),
    #  subtitle = sprintf("Using %s and %s (top-correlated factors)", top_factors[1], top_factors[2]),
      x = x_lab,
      y = y_lab
    )
  
  return(p)
}

# Omics colors: Simpsons palette from ggsci
library(ggsci)
# Get Simpsons colors
simpsons_pal <- pal_simpsons("springfield")(16)
# Omics colours: single source of truth shared by every tab (static PNGs, the
# feature table swatches AND the d3 circos). Kept off blue so they never clash
# with the red/blue used for correlation sign in the circos.
pal <- c(
  transcriptomics = "#FF7F0E",  # orange
  proteomics      = "#2CA02C",  # green
  metabolomics    = "#9467BD",  # purple
  flux            = "#BCBD22",  # olive
  "single-cell RNA" = "#17BECF", # teal
  "pseudobulk RNA"  = "#E377C2"  # pink
)

# Edge colors - distinct from Simpsons node colors
edge_colors <- c(
  positive = "#E91E63",  # Magenta - no conflict
  negative = "#00BCD4"   # Cyan - no conflict
)

# Function to get color for any omics type (with fallback)
get_omics_color <- function(omics_name) {
  if(omics_name %in% names(pal)) {
    return(pal[omics_name])
  } else {
    return(simpsons_pal[6])  # Light green fallback
  }
}

# Function to get metadata colors using viridis
get_metadata_colors <- function(n = NULL, continuous = FALSE) {
  library(viridis)
  if(continuous) {
    return(viridis(100))
  } else {
    if(is.null(n)) n <- 8
    return(viridis(n, option = "viridis", begin = 0, end = 1))
  }
}


view2omics <- function(v) {
  x <- as.character(v)
  
  # Define standard mappings
  omics_map <- c(
    "proc_rnaseq" = "transcriptomics",
    "proc_prot_combat" = "proteomics",
    "proc_metabolite" = "metabolomics",
    "proc_flux"   = "flux",
    "proc_metab"  = "metabolics",
    "proc_scrna"  = "single-cell RNA",
    "proc_pbrna"  = "pseudobulk RNA"
  )
  
  # Apply mapping with fallback to original name
  out <- sapply(x, function(name) {
    if(name %in% names(omics_map)) {
      omics_map[name]
    } else {
      name  # Keep original if not in map
    }
  }, USE.NAMES = FALSE)
  
  return(out)
}

omics2view <- function(o) { 
  x <- as.character(o)
  
  # Define reverse mappings
  view_map <- c(
    "transcriptomics" = "proc_rnaseq",
    "proteomics"      = "proc_prot_combat",
    "metabolomics"    = "proc_metabolite",
    "flux"            = "proc_flux",
    "metabolics"      = "proc_metab",
    "single-cell RNA" = "proc_scrna",
    "pseudobulk RNA"  = "proc_pbrna"
  )
  
  # Apply mapping with fallback to original name  
  out <- sapply(x, function(name) {
    if(name %in% names(view_map)) {
      view_map[name]
    } else {
      name  # Keep original if not in map
    }
  }, USE.NAMES = FALSE)
  
  return(out)
}




PerformCausalMediation <- function(input,outcome_type  ,  sims = 1000, boot = TRUE){
  
  
  names(input) <-c( "X", "M","Y")
  if (outcome_type == "disc" && !all(input$Y %in% c(0, 1))) {
    stop("For a binary outcome, Y must be coded 0/1.")
  }
  
  ## ---- 1. Mediator model (always linear here) ----
  model.M <- lm(M ~ X, data = input)
  
  ## ---- 2. Outcome model, conditional on type ----
  if (outcome_type == "disc") {
    model.Y <- glm(Y ~ X + M, data = input, family = binomial())
  } else {
    model.Y <- lm(Y ~ X + M, data = input)
  }
  
  ## ---- 3. Mediation analysis ----
  library(mediation)
  med.out <- mediate(model.M, model.Y,
                     treat     = "X",
                     mediator  = "M",
                     boot      = boot,
                     sims      = sims)
  
  ## ---- 4. Tidy results table ----
  if(outcome_type == "disc"){
    out_tbl <- data.frame(
      Effect = c("ACME (control)", "ACME (treated)", 
                 "ADE (control)", "ADE (treated)", 
                 "Total Effect", 
                 "Prop. Mediated (control)", "Prop. Mediated (treated)", 
                 "ACME (average)", "ADE (average)", "Prop. Mediated (average)"),
      
      Estimate = c(med.out$d0, med.out$d1,
                   med.out$z0, med.out$z1,
                   med.out$tau.coef,
                   med.out$n0, med.out$n1,
                   med.out$d.avg, med.out$z.avg, med.out$n.avg),
      
      CI.lower = c(med.out$d0.ci[1], med.out$d1.ci[1],
                   med.out$z0.ci[1], med.out$z1.ci[1],
                   med.out$tau.ci[1],
                   med.out$n0.ci[1], med.out$n1.ci[1],
                   med.out$d.avg.ci[1], med.out$z.avg.ci[1], med.out$n.avg.ci[1]),
      
      CI.upper = c(med.out$d0.ci[2], med.out$d1.ci[2],
                   med.out$z0.ci[2], med.out$z1.ci[2],
                   med.out$tau.ci[2],
                   med.out$n0.ci[2], med.out$n1.ci[2],
                   med.out$d.avg.ci[2], med.out$z.avg.ci[2],  med.out$n.avg.ci[2]),
      
      p.value = c(med.out$d0.p, med.out$d1.p,
                  med.out$z0.p, med.out$z1.p,
                  med.out$tau.p,
                  med.out$n0.p, med.out$n1.p,
                  med.out$d.avg.p, med.out$z.avg.p,  med.out$n.avg.p),
      
      stringsAsFactors = FALSE
    )
  }else{
    
  }
  
  return(out_tbl)
}
################################################################################
impute_median <- function(M) {
  # M: features x samples; fill NAs per feature (row median)
  if (!is.matrix(M)) M <- as.matrix(M)
  apply(M, 1, function(x){ m <- stats::median(x, na.rm = TRUE); x[is.na(x)] <- m; x }) |>
    t()
}

npn_transform <- function(M) {
  # Nonparanormal (Gaussian copula) transform, robust to outliers/heavy tails.
  # huge.npn expects samples x variables; here M is features x samples.
  X <- t(M)
  X_npn <- huge::huge.npn(X, npn.func = "shrinkage")   # returns samples x variables
  t(X_npn)  # back to features x samples
}

cross_cor_shrink <- function(X, Y, use_npn = TRUE) {
  stopifnot(ncol(X) == ncol(Y))
  X <- impute_median(X); Y <- impute_median(Y)
  if (use_npn) { X <- npn_transform(X); Y <- npn_transform(Y) }
  
  # corpcor::cov.shrink needs samples x variables
  Z <- cbind(t(X), t(Y))           # samples x (pX + pY)
  S <- corpcor::cov.shrink(Z)      # shrinkage covariance
  R <- cov2cor(S)                  # correlation matrix
  
  pX <- nrow(X); pY <- nrow(Y)
  # extract cross-block correlation (X rows, Y cols)
  RXY <- R[seq_len(pX), pX + seq_len(pY), drop = FALSE]
  rownames(RXY) <- rownames(X); colnames(RXY) <- rownames(Y)
  RXY
}

.fz <- function(r) 0.5 * log((1 + pmin(pmax(r, -0.999999), 0.999999)) /
                               (1 - pmin(pmax(r, -0.999999), 0.999999)))



diff_cross_cor_continous  <- function(X, Y, covariate, B = 1000, df = 3,
                                     use_npn = TRUE, seed = 1) {
  stopifnot(ncol(X) == ncol(Y))
  set.seed(seed)
  library(dplyr)
  n <- ncol(X)
  pX <- nrow(X); pY <- nrow(Y)
  
  # impute + optional NPN
  X <- impute_median(X); Y <- impute_median(Y)
  if (use_npn) { X <- npn_transform(X); Y <- npn_transform(Y) }

  # observed: unadjusted
  S_unadj <- corpcor::cov.shrink(cbind(t(X), t(Y))) 

  R_unadj <- cov2cor(S_unadj) 
  RXY_unadj <- R_unadj[seq_len(pX), pX + seq_len(pY), drop = FALSE]  
  # observed: age-adjusted
  Xr <- resid_mat(X, covariate,df)
  Yr <- resid_mat(Y, covariate,df)
 
  S_adj <- corpcor::cov.shrink(cbind(t(Xr), t(Yr)))
 
  R_adj <- cov2cor(S_adj)
 
  RXY_adj <- R_adj[seq_len(pX), pX + seq_len(pY), drop = FALSE]
 
  # table of observed stats
  obs <- as.data.frame(as.table(RXY_unadj), stringsAsFactors = FALSE) |>
    dplyr::rename(x = Var1, y = Var2, r_unadj = Freq) |>
    dplyr::mutate(r_adj = RXY_adj[cbind(match(x, rownames(RXY_adj)),
                                        match(y, colnames(RXY_adj)))],
                  z_unadj = .fz(r_unadj),
                  z_adj   = .fz(r_adj),
                  dz      = z_unadj - z_adj)
  
  # --- permutations: shuffle covariate values
  perm_dz <- matrix(0, nrow = nrow(obs), ncol = B)
  for (b in seq_len(B)) {
 
    cov_perm <- sample(covariate)
    Xrp <- resid_mat(X, cov_perm,df)
    Yrp <- resid_mat(Y, cov_perm,df)
    R_perm <- cov2cor(corpcor::cov.shrink(cbind(t(Xrp), t(Yrp))))
    RXY_perm <- R_perm[seq_len(pX), pX + seq_len(pY), drop = FALSE]
    z_unadj <- .fz(RXY_unadj[cbind(match(obs$x, rownames(RXY_unadj)),
                                   match(obs$y, colnames(RXY_unadj)))])
    z_adj   <- .fz(RXY_perm[cbind(match(obs$x, rownames(RXY_perm)),
                                  match(obs$y, colnames(RXY_perm)))])
    perm_dz[, b] <- z_unadj - z_adj
  }
  
  # permutation p-values
  p_perm <- rowMeans(abs(perm_dz) >= abs(obs$dz), na.rm = TRUE)
  
  res <- obs |>
    dplyr::mutate(p_perm = (p_perm * B + 1) / (B + 1),
                  FDR = p.adjust(p_perm, "BH")) |>
    dplyr::arrange(FDR, desc(abs(dz)))
  
  res[,c(1,2,4,3,6,5,7:9)]
}


edge_gate <- function(dc, ptype, r_min = 0.3, fdr_max = 0.05) {
  # Step 1: Filter by significance
  if(ptype == "pval") {
    dc <- dc %>% filter(p_perm <= fdr_max) 
  } else {
    dc <- dc %>% filter(FDR <= fdr_max) 
  }
  
  # Step 2: Classify edge types based on correlation presence/absence
  dc <- dc %>%
    mutate(
      # Define "edge exists" as |r| > r_min
      edge_A = abs(rA) >= r_min,
      edge_B = abs(rB) >= r_min,
      
      # Classify edge type
      edge_type = case_when(
        # Gain/Loss edges: present in one group but not the other
        edge_A & !edge_B ~ "A-specific",
        !edge_A & edge_B ~ "B-specific",
        
        # Sign-flip edges: present in both but opposite signs
        edge_A & edge_B & sign(rA) != sign(rB) ~ "Sign-flip",
        
        # Magnitude-change edges: same sign but different strength
        edge_A & edge_B & sign(rA) == sign(rB) ~ "Magnitude-change",
        
        # Weak-weak but differential: neither passes threshold but significantly different
        !edge_A & !edge_B ~ "Weak-differential",
        
        TRUE ~ "Other"
      ),
      
      # Calculate useful metrics
      abs_r_max = pmax(abs(rA), abs(rB)),
      r_diff = abs(rA - rB),
      stronger_in = if_else(abs(rA) > abs(rB), "A", "B")
    )
  
  # Step 3: Keep edges that meet criteria
  # For gain/loss: require strong correlation in at least one group
  # For others: require significance of difference
  dc %>%
    filter(
      edge_type %in% c("A-specific", "B-specific", "Sign-flip", "Magnitude-change")
    ) %>%
    transmute(
      from = x, 
      to = y, 
      rA, 
      rB, 
      dz, 
      p_perm, 
      FDR,
      edge_type,
      abs_r_max,
      r_diff,
      stronger_in
    )
}

resid_mat <- function(M, cov,df) {
  Z <- model.matrix(~ splines::ns(cov, df = df))
  B <- solve(crossprod(Z), crossprod(Z, t(M)))
  R <- t(M) - Z %*% B
  t(R)
}


diff_cross_cor <- function(X=selist[[1]], Y=selist[[2]], group=groups, B = 1000, use_npn = TRUE,
                           prefilter = TRUE, top_per_feature = 50, seed = 1) {
  stopifnot(ncol(X) == ncol(Y))
  stopifnot(nlevels(group) == 2)
  set.seed(seed)
  library(dplyr)

  # Compute group-specific cross correlations (shrinkage)
  lev <- levels(group); l1 <- lev[1]; l2 <- lev[2]
  RXA <- cross_cor_shrink(X[, group == l1, drop = FALSE],
                          Y[, group == l1, drop = FALSE],
                          use_npn = use_npn)
  RXB <- cross_cor_shrink(X[, group == l2, drop = FALSE],
                          Y[, group == l2, drop = FALSE],
                          use_npn = use_npn)
  
  # Pre-screen to limit number of pairs (optional but recommended when p_x*p_y is large)
  if (prefilter) {
    # for each X feature, keep top |r| pairs in either group
    keep_pairs <- function(R) {
      as.data.frame(as.table(R), stringsAsFactors = FALSE) |>
        rename(x = Var1, y = Var2, r = Freq) |>
        group_by(x) |>
        slice_max(order_by = abs(r), n = top_per_feature, with_ties = FALSE) |>
        ungroup()
    }
    dfA <- keep_pairs(RXA); dfB <- keep_pairs(RXB)
    cand <- bind_rows(dfA, dfB) |>
      distinct(x, y)
  } else {
    cand <- as.data.frame(expand.grid(x = rownames(RXA), y = colnames(RXA), stringsAsFactors = FALSE))
  }

  
  # observed stats
  obs <- cand |>
         mutate(rA = RXA[cbind(match(x, rownames(RXA)), match(y, colnames(RXA)))],
           rB = RXB[cbind(match(x, rownames(RXB)), match(y, colnames(RXB)))],
           zA = .fz(rA), zB = .fz(rB),
           dz = zA - zB)
  
  # Permutation test: shuffle labels, recompute dz each time
  # (We recompute shrinkage correlations per perm on the WHOLE cross-block – heavy but robust.)
  # To save compute, we reuse the NPN transform once upfront.
  X_imp <- impute_median(X); Y_imp <- impute_median(Y)
  if (use_npn) { X_imp <- npn_transform(X_imp); Y_imp <- npn_transform(Y_imp) }
  
  n <- length(group)
  perm_dz <- matrix(0, nrow = nrow(obs), ncol = B)
  
  for (b in seq_len(B)) {
    gperm <- sample(group)  # permute labels (keeps group sizes)
    RXA_p <- cross_cor_shrink(X_imp[, gperm == l1, drop = FALSE],
                              Y_imp[, gperm == l1, drop = FALSE],
                              use_npn = FALSE)     # already npn'd
    RXB_p <- cross_cor_shrink(X_imp[, gperm == l2, drop = FALSE],
                              Y_imp[, gperm == l2, drop = FALSE],
                              use_npn = FALSE)
    
    zA_p <- .fz(RXA_p[cbind(match(obs$x, rownames(RXA_p)), match(obs$y, colnames(RXA_p)))])
    zB_p <- .fz(RXB_p[cbind(match(obs$x, rownames(RXB_p)), match(obs$y, colnames(RXB_p)))])
    perm_dz[, b] <- zA_p - zB_p
  }
  
  # two-sided p via permutation
  p_perm <- rowMeans(abs(perm_dz) >= abs(obs$dz), na.rm = TRUE)
  res <- obs |>
       mutate(p_perm = (p_perm * B + 1) / (B + 1)) |>
     mutate(FDR = p.adjust(p_perm, method = "BH")) |>
   arrange(FDR, desc(abs(dz)))
  
  res
}

auto_trim <- function(mat, cap) {
  mat <- as.matrix(mat)
  madv <- apply(mat, 1, mad, na.rm = TRUE)
  keep <- madv > 0
  mat <- mat[keep, , drop = FALSE]
  madv <- madv[keep]

  if (nrow(mat) > cap) {
    mat <- mat[order(madv, decreasing = TRUE)[seq_len(cap)], , drop = FALSE]
  }
  mat
}

plot_multi_auroc <- function(model, n_comp = 3) {
  library(ggplot2)
  
  # 1. Use a NULL device to prevent Quartz errors in Rserve
  pdf(NULL)
  
  block_names <- names(model$X)
  # Generate the master list with TRUE plots
  p_temp <- auroc(model, roc.block = block_names, roc.comp = n_comp, plot = TRUE)
  
  all_roc_data <- list()
  
  for (bn in block_names) {
    # 2. Extract the actual ggplot object
    graph_name <- paste0("graph.", bn)
    p_obj <- p_temp[[graph_name]][[1]] 
    
    # 3. Find the ROC curve layer coordinates
    built_layers <- ggplot_build(p_obj)$data
    coords <- NULL
    for(layer in built_layers) {
      if(!is.null(layer$x) && length(layer$x) > 1) {
        coords <- layer
        break
      }
    }
    
    # 4. Get exact cumulative AUC
    # Proteomics: 0.9849, Transcriptomics: 1
    auc_val <- round(as.numeric(p_temp[[bn]][[paste0("comp", n_comp)]][1]), 4)
    
    all_roc_data[[bn]] <- data.frame(
      x = coords$x,
      y = coords$y,
      Omics = paste0(bn, " (AUC: ", auc_val, ")")
    )
  }
  
  dev.off() # Close NULL device
  
  # 5. Rebuild the combined plot with custom colors
  final_df <- do.call(rbind, all_roc_data)
  
  # Extract omics names and create color mapping
  omics_colors <- sapply(block_names, get_omics_color)
  names(omics_colors) <- sapply(block_names, function(bn) {
    # Match the format "omics_name (AUC: value)"
    auc_val <- round(as.numeric(p_temp[[bn]][[paste0("comp", n_comp)]][1]), 4)
    paste0(bn, " (AUC: ", auc_val, ")")
  })
  
  ggplot(final_df, aes(x = x, y = y, color = Omics)) +
    geom_line(linewidth = 1.2) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    scale_color_manual(values = omics_colors) +
    scale_x_continuous(limits = c(0, 100)) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(title = "", x = "100 - Specificity (%)", y = "Sensitivity (%)") +
    theme_minimal() +
    theme(
      # Places legend inside the plot area
      legend.position = c(0.70, 0.15), 
      legend.title = element_blank(),
      # Shrinks the text size
      legend.text = element_text(size = 8),
      # Removes background box for a cleaner look
      legend.background = element_blank(),
      legend.key = element_blank()
    ) +
    # Shrinks the colored line keys in the legend
    guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5))
}



plot_continuous_performance <- function(model, n_comp = 3) {
  library(ggplot2)
  
  # 1. Headless safe device
  pdf(NULL)
  
  # 2. Dynamically identify the Y index (it is 3 in your image)
  y_index <- model$indY
  block_names <- names(model$X)
  
  # 3. Extract the actual Y values from the correct slot
  actual_values <- as.numeric(model$X[[y_index]])
  
  # 4. FIX THE ERROR: Subset model$X to ONLY include predictors
  # This removes the "Y" block so predict() doesn't fail
  predictor_data <- model$X[-y_index]
  
  # 5. Generate predictions
  pred_res <- predict(model, newdata = predictor_data)
  
  all_data <- list()
  predictor_names <- names(predictor_data)
  
  for (bn in predictor_names) {
    # Extract cumulative predicted values
    predicted_values <- as.numeric(pred_res$predict[[bn]][, 1, n_comp])
    
    # Calculate R-squared
    r2_val <- round(cor(actual_values, predicted_values, use = "complete.obs")^2, 4)
    
    all_data[[bn]] <- data.frame(
      Observed = actual_values,
      Predicted = predicted_values,
      Omics = paste0(bn, " (R2: ", r2_val, ")")
    )
  }
  
  dev.off() 
  
  # 6. Final Plot with custom colors
  df_final <- do.call(rbind, all_data)
  
  # Create color mapping
  omics_colors <- sapply(predictor_names, get_omics_color)
  names(omics_colors) <- sapply(predictor_names, function(bn) {
    r2_val <- round(cor(actual_values, as.numeric(pred_res$predict[[bn]][, 1, n_comp]), use = "complete.obs")^2, 4)
    paste0(bn, " (R2: ", r2_val, ")")
  })
  
  ggplot(df_final, aes(x = Observed, y = Predicted, color = Omics)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    scale_color_manual(values = omics_colors) +
    theme_minimal() +
    theme(
      legend.position = c(0.70, 0.15), 
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.background = element_blank()
    ) +
    guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
    labs(x = "Actual Value", y = paste("Predicted (Comp 1-", n_comp, ")"))
}

GetORAMultiOmics <- function(library, fdr, collapse, mode, omics_for_pathway = NULL, use_network_only = FALSE, circos_features = NULL) {
    # Precomputed WGCNA / MOFA pathway routing (same ora Java step as DIABLO).
    if (!is.null(omics_for_pathway) && grepl("^PC(WGCNA|MOFA),", omics_for_pathway)) {
        return(.pc_run_pathway(omics_for_pathway, library, fdr, collapse, mode, analysis = "ora"))
    }
    top_features <- qs::qread("top_features.qs")
    kpList <- qs::qread("omics_kpList.qs")
    
    # Use provided omics or default to transcriptomics
    if(!is.null(omics_for_pathway) && omics_for_pathway != "") {
        selected_views <- strsplit(omics_for_pathway, ",")[[1]]
    } else {
        selected_views <- c("proc_rnaseq", "proc_nanostring_merge", "proc_prot_combat")
    }
    
    # Get ALL features from selected views as background
    background_features <- NULL
    for(view in names(kpList)) {
        if(view %in% selected_views) {
            background_features <- c(background_features, rownames(kpList[[view]]))
        }
    }
    background_features <- unique(background_features)
    
    if(length(background_features) == 0) {
        return("NO;No features found in selected omics for pathway analysis")
    }
    
    # Filter by network features if requested
    if(use_network_only && file.exists("edges_A_B.qs")) {
        edges_A_B <- qs::qread("edges_A_B.qs")
        network_features <- unique(c(edges_A_B$from, edges_A_B$to))
        background_features <- intersect(background_features, network_features)
        
        if(length(background_features) == 0) {
            return("NO;No network features found in selected omics")
        }
    }
    
    # Foreground (hits) = the CIRCOS-selected features when provided ('|'-separated feature IDs
    # from the frontend circos selection), else the full DIABLO signature. This makes the pathway
    # test exactly the feature set shown in the circos.
    #
    # DELIMITER IS '|', NOT ','. Feature names legitimately contain commas -- 28 metabolites and
    # 114 contaminants ("1,2,3,4,6,7,8-HPCDD", "Chlordane, alpha (cis)") -- so the old comma split
    # shredded them into fragments matching nothing: no error, silently wrong pathway results.
    # ';' is unusable as well (multi-gene protein symbols such as GATD3A;GATD3B). '|' appears in
    # 0 of 94,168 distinct feature names across every omics table.
    # fixed = TRUE is REQUIRED: '|' is a regex metacharacter, and strsplit(x, "|") without it
    # splits on the empty string, i.e. into individual characters.
    hits <- if (!is.null(circos_features) && nchar(circos_features) > 0) {
        unique(trimws(strsplit(circos_features, "|", fixed = TRUE)[[1]]))
    } else {
        unique(top_features$feature)
    }
    
    # Filter hits by network if requested
    if(use_network_only && file.exists("edges_A_B.qs")) {
        edges_A_B <- qs::qread("edges_A_B.qs")
        network_features <- unique(c(edges_A_B$from, edges_A_B$to))
        hits <- intersect(hits, network_features)
    }
    
    # Convert to entrez if needed
    out_anno <- qs::qread(paste0(other.tables.path,"libraries/gene_anno_full.qs")) 
    
    # Map feature IDs -> entrez using YOUR canonical per-feature mapping (feature_anno.qs, built from
    # the proc tables' own gene_id). Exact, no symbol re-matching. Falls back to the annotation file
    # (ENSG only) if feature_anno.qs is absent (an integration run made before this feature existed).
    if (file.exists("feature_anno.qs")) {
        .fa   <- qs::qread("feature_anno.qs")
        .fmap <- setNames(as.character(.fa$entrez), as.character(.fa$feature))
        .toEntrez <- function(v) { e <- unname(.fmap[as.character(v)]); e[is.na(e) | e == ""] <- NA; unique(na.omit(e)) }
        background_features <- .toEntrez(background_features)
        hits                <- .toEntrez(hits)
    } else {
        background_features[grepl("ENSG", background_features)] <-
            out_anno$entrez_selected[match(background_features[grepl("ENSG", background_features)], out_anno$ensembl_gene_id)]
        background_features <- unique(na.omit(background_features))
        hits[grepl("ENSG", hits)] <-
            out_anno$entrez_selected[match(hits[grepl("ENSG", hits)], out_anno$ensembl_gene_id)]
        hits <- unique(na.omit(hits))
    }

    # ORA universe = the screened features (5000 RNA / 3000 protein) that are PRESENT IN the selected
    # pathway library (measured-in-analysis  ∩  library) -- not all measured genes, per your design.
    .lib_file <- paste0(other.tables.path, "libraries/", library, ".rds")
    if (file.exists(.lib_file)) {
        .lib_genes <- unique(as.character(unlist(readRDS(.lib_file)$sets)))
        background_features <- intersect(background_features, .lib_genes)
        hits                <- intersect(hits, .lib_genes)
    }

    # Create a lookup table for loading scores by entrez ID
    loading_lookup <- data.frame(
        entrez = character(),
        loading = numeric(),
        stringsAsFactors = FALSE
    )
    
    for(i in 1:nrow(top_features)) {
        gene_id <- top_features$feature[i]
        if(grepl("ENSG", gene_id)) {
            gene_id <- out_anno$entrez_selected[match(gene_id, out_anno$ensembl_gene_id)]
        }
        if(!is.na(gene_id)) {
            loading_lookup <- rbind(loading_lookup, data.frame(
                entrez = as.character(gene_id),
                loading = top_features$mean_abs_loading[i],
                stringsAsFactors = FALSE
            ))
        }
    }
    
    # Create dea_results.csv: background=NS, hits=sig
    dea_mat <- data.frame(
        Feature = out_anno$hgnc_symbol[match(background_features,out_anno$entrez_selected)],
        Gene_ID = background_features,
        Log2FC = 0,
        t.stat = 0,
        Description = out_anno$description[match(background_features,out_anno$entrez_selected)],
        sig = "NS",
        stringsAsFactors = FALSE
    )
   dea_mat$Feature[is.na(dea_mat$Feature)|dea_mat$Feature==""] <- dea_mat$Gene_ID[is.na(dea_mat$Feature)|dea_mat$Feature==""]
    dea_mat$sig[dea_mat$Gene_ID %in% hits] <- "sig"
    
    # Add loading scores to dea_mat for matching genes
    for(i in 1:nrow(loading_lookup)) {
        idx <- which(dea_mat$Gene_ID == loading_lookup$entrez[i])
        if(length(idx) > 0) {
            dea_mat$Log2FC[idx] <- loading_lookup$loading[i]
            dea_mat$t.stat[idx] <- loading_lookup$loading[i]
        }
    }
 dea_mat = dea_mat[,c("Feature", "Log2FC" ,  "t.stat","Gene_ID", "Description", "sig")]
    write.csv(dea_mat, "dea_results.csv", row.names = FALSE)
    
    rcmd <<- paste0('performORA("', library, '", "', fdr, '", "', collapse, '", "', mode, '")')
    result <- performORA(library, fdr, collapse, mode)
   
    return(result)
}

GetGSEAMultiOmics <- function(library, rank, fdr, collapse, mode, omics_for_pathway = NULL, use_network_only = FALSE, circos_features = NULL) {
    # Precomputed WGCNA / MOFA pathway routing (same gsea Java step as DIABLO).
    # Must run BEFORE the screen_tstat.qs branch, which would otherwise mis-split the token.
    if (!is.null(omics_for_pathway) && grepl("^PC(WGCNA|MOFA),", omics_for_pathway)) {
        return(.pc_run_pathway(omics_for_pathway, library, fdr, collapse, mode, analysis = "gsea", rank = rank))
    }
    # ---- GSEA ranks ALL features by the SIGNED limma t of the contrast (screen_tstat.qs) ----------
    # GSEA needs a full-genome ranking; sparse DIABLO loadings can't provide one (0 for every non-
    # selected feature). So we rank the whole block by the moderated t computed at the screen step
    # (covariate-adjusted when covariates were selected). This is the group-CONTRAST GSEA -- a
    # different, powered question from the ORA of the DIABLO signature. circos_features is intentionally
    # ignored here (GSEA uses all measured genes). Falls through to the legacy path if screen_tstat absent.
    if (file.exists("screen_tstat.qs")) {
        sv <- if (!is.null(omics_for_pathway) && omics_for_pathway != "") strsplit(omics_for_pathway, ",")[[1]] else c("proc_rnaseq", "proc_nanostring_merge", "proc_prot_combat")
        out_anno     <- qs::qread(paste0(other.tables.path,"libraries/gene_anno_full.qs"))
        screen_tstat <- qs::qread("screen_tstat.qs")
        tvec <- numeric(0)
        for (view in intersect(names(screen_tstat), sv)) tvec <- c(tvec, screen_tstat[[view]])
        if (length(tvec) > 0) {
            feat_ids <- names(tvec)
            if (file.exists("feature_anno.qs")) {
                .fa <- qs::qread("feature_anno.qs")
                gene_ids <- unname(setNames(as.character(.fa$entrez), as.character(.fa$feature))[feat_ids])
            } else {
                gene_ids <- feat_ids
                .e <- grepl("ENSG", gene_ids)
                gene_ids[.e] <- out_anno$entrez_selected[match(gene_ids[.e], out_anno$ensembl_gene_id)]
            }
            keep <- !is.na(gene_ids) & gene_ids != "" & !is.na(tvec)
            gene_ids <- gene_ids[keep]; tvec <- as.numeric(tvec[keep]); feat_ids <- feat_ids[keep]
            ord <- order(abs(tvec), decreasing = TRUE)   # one row per gene: strongest |t|
            gene_ids <- gene_ids[ord]; tvec <- tvec[ord]; feat_ids <- feat_ids[ord]
            dup <- duplicated(gene_ids); gene_ids <- gene_ids[!dup]; tvec <- tvec[!dup]; feat_ids <- feat_ids[!dup]
            sym <- out_anno$hgnc_symbol[match(gene_ids, out_anno$entrez_selected)]
            sym[is.na(sym) | sym == ""] <- feat_ids[is.na(sym) | sym == ""]
            dea_mat <- data.frame(
                Feature = sym, Log2FC = tvec, t.stat = tvec, Gene_ID = as.character(gene_ids),
                Description = out_anno$description[match(gene_ids, out_anno$entrez_selected)],
                sig = "NS", stringsAsFactors = FALSE
            )
            write.csv(dea_mat, "dea_results.csv", row.names = FALSE)
            rcmd <<- paste0('performGSEA("', library, '", "', rank, '", "', fdr, '", "', collapse, '", "', mode, '")')
            return(performGSEA(library, rank, fdr, collapse, mode))
        }
    }
    top_features <- qs::qread("top_features.qs")
    kpList <- qs::qread("omics_kpList.qs")
    
    # Use provided omics or default to transcriptomics
    if(!is.null(omics_for_pathway) && omics_for_pathway != "") {
        selected_views <- strsplit(omics_for_pathway, ",")[[1]]
    } else {
        selected_views <- c("proc_rnaseq", "proc_nanostring_merge", "proc_prot_combat")
    }
    
    # Get ALL features from selected views as background
    background_features <- NULL
    for(view in names(kpList)) {
        if(view %in% selected_views) {
            background_features <- c(background_features, rownames(kpList[[view]]))
        }
    }
    background_features <- unique(background_features)
    
    if(length(background_features) == 0) {
        return("NO;No features found in selected omics for pathway analysis")
    }
    
    # Filter by network features if requested
    if(use_network_only && file.exists("edges_A_B.qs")) {
        edges_A_B <- qs::qread("edges_A_B.qs")
        network_features <- unique(c(edges_A_B$from, edges_A_B$to))
        background_features <- intersect(background_features, network_features)
        
        if(length(background_features) == 0) {
            return("NO;No network features found in selected omics")
        }
    }
    
    # Convert to entrez
    out_anno <- qs::qread(paste0(other.tables.path,"libraries/gene_anno_full.qs")) 
    
    background_features[grepl("ENSG", background_features)] <- 
        out_anno$entrez_selected[match(background_features[grepl("ENSG", background_features)], out_anno$ensembl_gene_id)]
    background_features <- unique(na.omit(background_features))
    
    # Create a lookup table for loading scores by entrez ID
    loading_lookup <- data.frame(
        entrez = character(),
        loading = numeric(),
        stringsAsFactors = FALSE
    )
    
    for(i in 1:nrow(top_features)) {
        gene_id <- top_features$feature[i]
        if(grepl("ENSG", gene_id)) {
            gene_id <- out_anno$entrez_selected[match(gene_id, out_anno$ensembl_gene_id)]
        }
        if(!is.na(gene_id)) {
            loading_lookup <- rbind(loading_lookup, data.frame(
                entrez = as.character(gene_id),
                loading = top_features$mean_abs_loading[i],
                stringsAsFactors = FALSE
            ))
        }
    }
    
    # Create dea_results with loading scores
   dea_mat <- data.frame(
        Feature = out_anno$hgnc_symbol[match(background_features,out_anno$entrez_selected)],
        Gene_ID = background_features,
        Log2FC = 0,
        t.stat = 0,
        Description = out_anno$description[match(background_features,out_anno$entrez_selected)],
        sig = "NS",
        stringsAsFactors = FALSE
    )
   dea_mat$Feature[is.na(dea_mat$Feature)|dea_mat$Feature==""] <- dea_mat$Gene_ID[is.na(dea_mat$Feature)|dea_mat$Feature==""]

    # Assign loading scores to DIABLO-selected features
    for(i in 1:nrow(loading_lookup)) {
        idx <- which(dea_mat$Gene_ID == loading_lookup$entrez[i])
        if(length(idx) > 0) {
            dea_mat$Log2FC[idx] <- loading_lookup$loading[i]
            dea_mat$t.stat[idx] <- loading_lookup$loading[i]
            dea_mat$sig[idx] <- "sig"
        }
    }
     dea_mat = dea_mat[,c("Feature", "Log2FC" ,  "t.stat","Gene_ID", "Description", "sig")]

    write.csv(dea_mat, "dea_results.csv", row.names = FALSE)

    rcmd <<- paste0('performGSEA("', library, '", "', rank, '", "', fdr, '", "', collapse, '", "', mode, '")')
    result <- performGSEA(library, rank, fdr, collapse, mode)

    return(result)
}


################################################################################
# On-demand selection STABILITY (mixOmics perf): reuse the saved DIABLO model,
# re-run the sparse selection under repeated CV, and annotate each feature with
# the fraction of folds it is reselected on its component (0-100%). Writes the
# stability back into all_features.json. Does NOT re-run the main integration.
################################################################################
GetSelectionStability <- function(folds = 5, nrepeat = 3) {
  if(!file.exists("diablo_model.qs") || !file.exists("all_features.qs")) {
    return("NO;No integration results found. Please run the integration first.")
  }
  suppressMessages(library(mixOmics))
  diablo_model <- qs::qread("diablo_model.qs")

  # Stability via perf() is defined for the discriminant model (block.splsda).
  if(!inherits(diablo_model, c("block.splsda", "sgccda"))) {
    return("NO;Selection stability is only available for group-comparison (discriminant) runs.")
  }

  # Clamp folds so every CV fold keeps samples of the smallest class.
  min_class <- suppressWarnings(min(table(diablo_model$Y)))
  folds   <- max(2, min(as.integer(folds), as.integer(min_class)))
  nrepeat <- max(1, as.integer(nrepeat))

  perf_res <- tryCatch(
    suppressWarnings(perf(diablo_model, validation = "Mfold", folds = folds, nrepeat = nrepeat)),
    error = function(e) NULL)
  if(is.null(perf_res) || is.null(perf_res$features$stable)) {
    return("NO;Stability computation failed (too few samples for cross-validation).")
  }

  # perf_res$features$stable: nrep -> block -> comp -> table(feature -> freq within that repeat).
  # Average the per-repeat frequency across repeats for each (feature, component); absent = 0.
  recs <- list(); stab <- perf_res$features$stable
  for(rep_i in seq_along(stab)) {
    for(bn in names(stab[[rep_i]])) {
      comps <- stab[[rep_i]][[bn]]
      for(ci in seq_along(comps)) {
        tb <- comps[[ci]]
        if(length(tb) == 0) next
        recs[[length(recs) + 1]] <- data.frame(
          feature = names(tb), component = ci,
          freq = as.numeric(tb), stringsAsFactors = FALSE)
      }
    }
  }
  if(length(recs) == 0) return("NO;No stability data returned by perf().")
  df  <- do.call(rbind, recs)
  agg <- aggregate(freq ~ feature + component, data = df, FUN = sum)
  agg$stability <- agg$freq / length(stab)          # mean across repeats (absent repeats = 0)

  all_features <- qs::qread("all_features.qs")
  all_features$stability <- agg$stability[match(paste(all_features$feature, all_features$component),
                                                paste(agg$feature, agg$component))]
  all_features$stability[is.na(all_features$stability)] <- 0
  all_features$stability <- round(all_features$stability * 100, 1)   # percent

  qs::qsave(all_features, "all_features.qs")
  jsonlite::write_json(all_features[, c("feature", "label", "omics", "mean_abs_loading",
                                        "loading_signed", "component", "color", "stability")],
                       "all_features.json")
  return(paste0("RES-OK;stability computed (", folds, "-fold x ", nrepeat, "-repeat)"))
}

################################################################################
# Lightweight re-filter: re-slice top features and regenerate feature plot only
# Reads saved all_features.qs, does NOT re-run DIABLO
################################################################################
RefreshTopFeatures <- function(topFeat = 100) {
  library(dplyr)
  library(ggplot2)

  if(!file.exists("all_features.qs")) {
    return("NO;No integration results found. Please run integration first.")
  }

  fig_counts <- qs::qread("fig_counts.qs")
  all_features <- qs::qread("all_features.qs")
  out_anno <- qs::qread(paste0(other.tables.path, "libraries/gene_anno_full.qs")) 

  topFeat <- as.integer(topFeat)
  topFeat <- max(10, min(topFeat, nrow(all_features)))

  top_features <- all_features %>% arrange(desc(mean_abs_loading)) %>% slice_head(n = topFeat)

  top_features$label <- as.character(top_features$feature)
  top_features$omics <- view2omics(top_features$view)
  top_features$label[top_features$label %in% out_anno$ensembl_gene_id] <- out_anno$hgnc_symbol[match(top_features$label[top_features$label %in% out_anno$ensembl_gene_id], out_anno$ensembl_gene_id)]
  top_features$label[top_features$label %in% out_anno$entrez_selected] <- out_anno$hgnc_symbol[match(top_features$label[top_features$label %in% out_anno$entrez_selected], out_anno$entrez_selected)]
  top_features <- top_features %>% mutate(label = factor(label, levels = rev(unique(label))))

  feature_colors <- sapply(unique(top_features$omics), get_omics_color)
  names(feature_colors) <- unique(top_features$omics)

  nm3 <- paste0("sel_feature_1_", fig_counts, ".png")
  p3 <- ggplot(top_features, aes(x = mean_abs_loading, y = label, colour = omics)) +
    geom_segment(aes(x = 0, xend = mean_abs_loading, yend = label), linewidth = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = feature_colors, name = "Omics") +
    labs(title = " ", x = "Loading", y = NULL) +
    coord_cartesian(xlim = c(0, max(top_features$mean_abs_loading, na.rm = TRUE) * 1.08)) +
    theme_minimal(base_size = 9, base_family = "Arial") +
    theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
          axis.title.x = element_text(colour = "black"), plot.title = element_text(face = "bold"), legend.position = "right")
  n_features <- nrow(top_features)
  ggsave(nm3, p3, width = 6, height = max(4, n_features * 0.12), dpi = 300)

  # Update saved top_features so network uses the new selection
  qs::qsave(top_features, "top_features.qs")

  # Build feature count message
  feature_counts <- table(top_features$view)
  feature_msg <- paste(sapply(names(feature_counts), function(nm) {
    paste0(feature_counts[nm], " ", nm)
  }), collapse = " and ")
  msg <- sprintf("Showing top %d features (%s) ranked by DIABLO loading score.", n_features, feature_msg)

  return(paste0("RES-OK;", nm3, ";", msg))
}


################################################################################
################################################################################
##########################  WGCNA FUNCTIONS  ###################################
################################################################################
################################################################################

# Intersect `shared` with the multi-omics donor subset on disk (donors_multiomics.rds),
# kept separate from donors.rds so the multi-omics subset cannot leak into other views.
.apply_multiomics_subset <- function(shared, donors) {
  if (identical(donors, "subset") && file.exists("donors_multiomics.rds")) {
    donor_list <- readRDS("donors_multiomics.rds")
    return(intersect(donor_list, shared))
  }
  shared
}

GetWGCNASummary <- function(omicsType, donors = "all") {
  wgcna_dir <- paste0(other.tables.path, "multi_omics/wgcna_results/", omicsType, "/")

  if (!dir.exists(wgcna_dir)) {
    return(paste0("NO;WGCNA results not found for omics type: ", omicsType))
  }

  # Load precomputed files
  summary_raw <- jsonlite::fromJSON(paste0(wgcna_dir, "summary.json"))
  modules_df  <- read.csv(paste0(wgcna_dir, "modules.csv"), stringsAsFactors = FALSE)
  hub_df      <- read.csv(paste0(wgcna_dir, "hub_features.csv"), stringsAsFactors = FALSE)
  sft_df      <- read.csv(paste0(wgcna_dir, "sft_fit.csv"), stringsAsFactors = FALSE)

  # Build module summary (exclude grey/module 0)
  mod_summary <- list()
  for (m in sort(unique(modules_df$module_number))) {
    mod_feats <- modules_df[modules_df$module_number == m, ]
    mod_hubs  <- hub_df[hub_df$module == m, ]
    top_hub <- if (nrow(mod_hubs) > 0) {
      mod_hubs[which.max(mod_hubs$kME), ]
    } else {
      NULL
    }
    mod_summary[[length(mod_summary) + 1]] <- list(
      module_number = m,
      module_color  = mod_feats$module_color[1],
      size          = nrow(mod_feats),
      top_hub_id    = if (!is.null(top_hub)) top_hub$feature_id else NA,
      top_hub_name  = if (!is.null(top_hub)) top_hub$feature_name else NA,
      top_hub_kME   = if (!is.null(top_hub)) top_hub$kME else NA,
      avg_kME       = if (nrow(mod_hubs) > 0) round(mean(mod_hubs$kME), 4) else NA
    )
  }

  combined <- list(
    summary       = summary_raw,
    module_summary = mod_summary,
    hub_features  = hub_df,
    sft           = sft_df
  )

  out_file <- paste0("wgcna_summary_", omicsType, ".json")
  jsonlite::write_json(combined, out_file, auto_unbox = TRUE, pretty = TRUE)

  n_mod <- summary_raw$n_modules
  n_feat <- summary_raw$n_features
  n_don <- summary_raw$n_donors
  pw <- summary_raw$power
  msg <- sprintf(
    "WGCNA detected %d modules from %d features across %d donors (power=%d, %s network).",
    n_mod, n_feat, n_don, pw, summary_raw$network_type
  )

  return(paste0("RES-OK;", out_file, ";", msg))
}


################################################################################
GetWGCNATrait <- function(
    omicsType,
    varGroup,
    analysisVar,
    class1 = NULL,
    class2 = NULL,
    primaryType = "disc",
    donors = "all") {

  wgcna_dir <- paste0(other.tables.path, "multi_omics/wgcna_results/", omicsType, "/")

  # Load eigengenes (rows = donors, cols = ME1, ME2, ...)
  eigen_df <- read.csv(paste0(wgcna_dir, "eigengenes.csv"), row.names = 1, stringsAsFactors = FALSE)
  summary_raw <- jsonlite::fromJSON(paste0(wgcna_dir, "summary.json"))

  # Load metadata
  meta <- read.csv(paste0(other.tables.path, "display_data/metadata_sum_norm.csv"), stringsAsFactors = FALSE)

  # Inject cluster (or any proc_metadata-only) column if selected var is missing from the CSV
  if (!is.null(analysisVar) && !(analysisVar %in% colnames(meta))) {
    library(RSQLite)
    meta_db <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
    pheno <- dbReadTable(meta_db, "proc_metadata")
    dbDisconnect(meta_db)
    if (analysisVar %in% colnames(pheno)) {
      extra <- pheno[, c("record_id", analysisVar)]
      meta <- merge(meta, extra, by = "record_id", all.x = TRUE)
    } else {
      return(paste0("NO;Variable '", analysisVar,
                    "' not found in metadata or donor table."))
    }
  }

  # Find shared donors (multi-omics subset lives in donors_multiomics.rds, not donors.rds)
  shared <- intersect(rownames(eigen_df), meta$record_id)
  shared <- .apply_multiomics_subset(shared, donors)

  if (length(shared) < 10) {
    return("NO;Too few shared donors between WGCNA eigengenes and metadata.")
  }

  eigen_sub <- eigen_df[shared, , drop = FALSE]
  meta_sub  <- meta[match(shared, meta$record_id), ]

  # Get display name for the variable
  proc_variable <- read.csv(paste0(other.tables.path, "display_interface/proc_variable_summary_v2.csv"), stringsAsFactors = FALSE)
  trait_display <- proc_variable$display[proc_variable$column == analysisVar]
  if (length(trait_display) == 0) trait_display <- analysisVar

  # Get module info
  modules_info <- summary_raw$modules

  # Compute correlations
  results <- list()

  if (primaryType == "disc") {
    if (length(class2) == 0) class2 <- "NULL"
    meta_filt <- meta_sub[meta_sub[[analysisVar]] %in% c(class1, class2), ]
    shared_filt <- meta_filt$record_id
    eigen_filt <- eigen_sub[shared_filt, , drop = FALSE]
    y_numeric <- ifelse(meta_filt[[analysisVar]] == class1, 0, 1)

    for (me in colnames(eigen_filt)) {
      ct <- cor.test(eigen_filt[[me]], y_numeric, method = "pearson")
      mod_num <- gsub("ME", "", me)
      mod_color <- if (!is.null(modules_info[[mod_num]])) modules_info[[mod_num]]$color else "grey"
      mod_size  <- if (!is.null(modules_info[[mod_num]])) modules_info[[mod_num]]$size else 0
      results[[length(results) + 1]] <- list(
        module      = me,
        color       = mod_color,
        size        = mod_size,
        correlation = round(ct$estimate, 4),
        pvalue      = ct$p.value
      )
    }
  } else {
    meta_filt <- meta_sub[!is.na(meta_sub[[analysisVar]]), ]
    shared_filt <- meta_filt$record_id
    eigen_filt <- eigen_sub[shared_filt, , drop = FALSE]
    y_val <- as.numeric(meta_filt[[analysisVar]])

    for (me in colnames(eigen_filt)) {
      ct <- cor.test(eigen_filt[[me]], y_val, method = "pearson")
      mod_num <- gsub("ME", "", me)
      mod_color <- if (!is.null(modules_info[[mod_num]])) modules_info[[mod_num]]$color else "grey"
      mod_size  <- if (!is.null(modules_info[[mod_num]])) modules_info[[mod_num]]$size else 0
      results[[length(results) + 1]] <- list(
        module      = me,
        color       = mod_color,
        size        = mod_size,
        correlation = round(ct$estimate, 4),
        pvalue      = ct$p.value
      )
    }
  }

  # Apply FDR correction
  pvals <- sapply(results, function(x) x$pvalue)
  fdrs  <- p.adjust(pvals, method = "BH")
  for (i in seq_along(results)) {
    results[[i]]$fdr <- round(fdrs[i], 6)
    results[[i]]$pvalue <- round(results[[i]]$pvalue, 6)
  }

  out <- list(
    trait         = analysisVar,
    trait_display = trait_display,
    modules       = results,
    n_donors      = length(shared_filt),
    primaryType   = primaryType,
    class1        = class1,
    class2        = class2
  )

  out_file <- paste0("wgcna_trait_", omicsType, ".json")
  jsonlite::write_json(out, out_file, auto_unbox = TRUE, pretty = TRUE)

  n_sig <- sum(fdrs < 0.05)
  msg <- sprintf(
    "Module-trait analysis: %d of %d modules significantly correlated with %s (FDR < 0.05, %d donors).",
    n_sig, length(results), trait_display, length(shared_filt)
  )

  return(paste0("RES-OK;", out_file, ";", msg))
}


################################################################################
GetWGCNANetwork <- function(omicsType, moduleNum, donors = "all") {
  # Network edges and nodes are precomputed from the full WGCNA model; module
  # membership doesn't change with a donor subset, so `donors` is accepted for
  # signature consistency but not applied here.

  wgcna_dir <- paste0(other.tables.path, "multi_omics/wgcna_results/", omicsType, "/")

  edges_df   <- read.csv(paste0(wgcna_dir, "network_edges.csv"), stringsAsFactors = FALSE)
  modules_df <- read.csv(paste0(wgcna_dir, "modules.csv"), stringsAsFactors = FALSE)
  hub_df     <- read.csv(paste0(wgcna_dir, "hub_features.csv"), stringsAsFactors = FALSE)
  summary_raw <- jsonlite::fromJSON(paste0(wgcna_dir, "summary.json"))

  moduleNum <- as.integer(moduleNum)

  # Filter edges for this module
  mod_edges <- edges_df[edges_df$module == moduleNum, , drop = FALSE]

  if (nrow(mod_edges) == 0) {
    return(paste0("NO;No network edges found for module ", moduleNum))
  }

  # Get unique node IDs from edges
  node_ids <- unique(c(mod_edges$source, mod_edges$target))

  # Build node info from modules_df and hub_df
  mod_info <- modules_df[modules_df$feature_id %in% node_ids, ]
  hub_info <- hub_df[hub_df$feature_id %in% node_ids, ]

  nodes <- lapply(node_ids, function(nid) {
    mi <- mod_info[mod_info$feature_id == nid, ]
    hi <- hub_info[hub_info$feature_id == nid, ]
    list(
      id    = nid,
      label = if (nrow(mi) > 0 && !is.na(mi$feature_name[1])) mi$feature_name[1] else nid,
      type  = omicsType,
      kME   = if (nrow(hi) > 0) hi$kME[1] else 0,
      connectivity = if (nrow(hi) > 0) hi$intramodular_connectivity[1] else 0
    )
  })

  edges <- lapply(seq_len(nrow(mod_edges)), function(i) {
    list(
      source = mod_edges$source[i],
      target = mod_edges$target[i],
      weight = mod_edges$weight[i]
    )
  })

  mod_color <- if (!is.null(summary_raw$modules[[as.character(moduleNum)]])) {
    summary_raw$modules[[as.character(moduleNum)]]$color
  } else {
    "grey"
  }

  net <- list(
    nodes = nodes,
    edges = edges,
    label = paste0("Module ", moduleNum, " (", mod_color, ")"),
    module_color = mod_color,
    n_nodes = length(nodes),
    n_edges = length(edges)
  )

  out_file <- paste0("wgcna_network_", omicsType, "_", moduleNum, ".json")
  jsonlite::write_json(net, out_file, auto_unbox = TRUE, pretty = TRUE)

  msg <- sprintf(
    "Module %d (%s): %d nodes, %d edges.",
    moduleNum, mod_color, length(nodes), length(edges)
  )

  return(paste0("RES-OK;", out_file, ";", msg))
}


################################################################################
GetWGCNACross <- function(omicsType1, omicsType2, donors = "all") {

  dir1 <- paste0(other.tables.path, "multi_omics/wgcna_results/", omicsType1, "/")
  dir2 <- paste0(other.tables.path, "multi_omics/wgcna_results/", omicsType2, "/")

  if (!dir.exists(dir1)) return(paste0("NO;WGCNA results not found for: ", omicsType1))
  if (!dir.exists(dir2)) return(paste0("NO;WGCNA results not found for: ", omicsType2))

  eigen1 <- read.csv(paste0(dir1, "eigengenes.csv"), row.names = 1, stringsAsFactors = FALSE)
  eigen2 <- read.csv(paste0(dir2, "eigengenes.csv"), row.names = 1, stringsAsFactors = FALSE)
  sum1   <- jsonlite::fromJSON(paste0(dir1, "summary.json"))
  sum2   <- jsonlite::fromJSON(paste0(dir2, "summary.json"))

  # Shared donors (intersect with multi-omics subset if active)
  shared <- intersect(rownames(eigen1), rownames(eigen2))
  shared <- .apply_multiomics_subset(shared, donors)

  if (length(shared) < 10) {
    return("NO;Too few shared donors between the two omics types.")
  }

  e1 <- eigen1[shared, , drop = FALSE]
  e2 <- eigen2[shared, , drop = FALSE]

  # Correlation matrix
  cor_mat <- matrix(NA, nrow = ncol(e1), ncol = ncol(e2))
  p_mat   <- matrix(NA, nrow = ncol(e1), ncol = ncol(e2))
  rownames(cor_mat) <- colnames(e1)
  colnames(cor_mat) <- colnames(e2)

  for (i in seq_len(ncol(e1))) {
    for (j in seq_len(ncol(e2))) {
      ct <- cor.test(e1[, i], e2[, j], method = "pearson")
      cor_mat[i, j] <- round(ct$estimate, 4)
      p_mat[i, j]   <- round(ct$p.value, 6)
    }
  }

  # Get module colors
  colors1 <- sapply(colnames(e1), function(me) {
    mn <- gsub("ME", "", me)
    if (!is.null(sum1$modules[[mn]])) sum1$modules[[mn]]$color else "grey"
  })
  colors2 <- sapply(colnames(e2), function(me) {
    mn <- gsub("ME", "", me)
    if (!is.null(sum2$modules[[mn]])) sum2$modules[[mn]]$color else "grey"
  })

  out <- list(
    omics1_modules = colnames(e1),
    omics2_modules = colnames(e2),
    omics1_colors  = as.character(colors1),
    omics2_colors  = as.character(colors2),
    correlations   = lapply(1:nrow(cor_mat), function(i) as.numeric(cor_mat[i, ])),
    pvalues        = lapply(1:nrow(p_mat), function(i) as.numeric(p_mat[i, ])),
    n_shared_donors = length(shared),
    omics1_name    = omicsType1,
    omics2_name    = omicsType2
  )

  out_file <- "wgcna_cross.json"
  jsonlite::write_json(out, out_file, auto_unbox = TRUE, pretty = TRUE)

  msg <- sprintf(
    "Cross-omics correlation: %s (%d modules) vs %s (%d modules) on %d shared donors.",
    omicsType1, ncol(e1), omicsType2, ncol(e2), length(shared)
  )

  return(paste0("RES-OK;", out_file, ";", msg))
}


################################################################################
GetWGCNACrossNetwork <- function(omicsType1, omicsType2, topN = 5, sigModules1 = NULL, sigModules2 = NULL, donors = "all") {

  wgcna_dir1 <- paste0(other.tables.path, "multi_omics/wgcna_results/", omicsType1, "/")
  wgcna_dir2 <- paste0(other.tables.path, "multi_omics/wgcna_results/", omicsType2, "/")
  expr_dir   <- paste0(other.tables.path, "processed_comprehensive/omics/")

  if (!dir.exists(wgcna_dir1)) return(paste0("NO;WGCNA results not found for: ", omicsType1))
  if (!dir.exists(wgcna_dir2)) return(paste0("NO;WGCNA results not found for: ", omicsType2))

  # --- Map omics type to expression file and ID column ---
  get_expr_info <- function(omicsType) {
    if (omicsType == "rnaseq") {
      list(file = "rnaseq_filtered.csv", id_col = "accession", meta_cols = 5)
    } else if (omicsType == "proteomics") {
      list(file = "proteomics.csv", id_col = "gene_id", meta_cols = 3)
    } else if (omicsType == "nanostring") {
      list(file = "nanostring.csv", id_col = "gene_id", meta_cols = 3)
    } else if (omicsType %in% c("pseudobulk_alpha", "pseudobulk_beta",
                                  "patchseq_alpha", "patchseq_beta")) {
      # patchseq_* kept as legacy aliases (data was renamed to pseudobulk_*.csv).
      cell <- sub("^(pseudobulk|patchseq)_", "", omicsType)
      list(file = paste0("pseudobulk_", cell, ".csv"),
           id_col = "gene_id", meta_cols = 3)
    } else if (omicsType %in% c("flux_HG", "flux_LG", "flux_ratio")) {
      list(file = paste0(omicsType, ".csv"), id_col = "rxn", meta_cols = 4)
    } else {
      NULL
    }
  }

  info1 <- get_expr_info(omicsType1)
  info2 <- get_expr_info(omicsType2)
  if (is.null(info1)) return(paste0("NO;Unsupported omics type: ", omicsType1))
  if (is.null(info2)) return(paste0("NO;Unsupported omics type: ", omicsType2))

  # --- Load hub features for both omics ---
  hub1 <- read.csv(paste0(wgcna_dir1, "hub_features.csv"), stringsAsFactors = FALSE)
  hub2 <- read.csv(paste0(wgcna_dir2, "hub_features.csv"), stringsAsFactors = FALSE)
  sum1 <- jsonlite::fromJSON(paste0(wgcna_dir1, "summary.json"))
  sum2 <- jsonlite::fromJSON(paste0(wgcna_dir2, "summary.json"))

  # Take top N hub features per module (by kME), excluding module 0
  get_top_hubs <- function(hub_df, topN) {
    hub_df <- hub_df[hub_df$module != 0, ]
    do.call(rbind, lapply(split(hub_df, hub_df$module), function(df) {
      df <- df[order(-df$kME), ]
      head(df, topN)
    }))
  }

  # Filter to trait-significant modules if provided (comma-separated string)
  if (!is.null(sigModules1) && nchar(sigModules1) > 0 && sigModules1 != "NA") {
    sig1 <- as.integer(unlist(strsplit(sigModules1, ",")))
    hub1 <- hub1[hub1$module %in% sig1, ]
  }
  if (!is.null(sigModules2) && nchar(sigModules2) > 0 && sigModules2 != "NA") {
    sig2 <- as.integer(unlist(strsplit(sigModules2, ",")))
    hub2 <- hub2[hub2$module %in% sig2, ]
  }

  top1 <- get_top_hubs(hub1, topN)
  top2 <- get_top_hubs(hub2, topN)

  if (is.null(top1) || is.null(top2) || nrow(top1) == 0 || nrow(top2) == 0) {
    return("NO;No hub features found for trait-significant modules. Try relaxing the significance threshold.")
  }

  # --- Load expression data ---
  expr1_raw <- read.csv(paste0(expr_dir, info1$file), stringsAsFactors = FALSE, check.names = FALSE)
  expr2_raw <- read.csv(paste0(expr_dir, info2$file), stringsAsFactors = FALSE, check.names = FALSE)

  # Set row names to feature IDs
  rownames(expr1_raw) <- as.character(expr1_raw[[info1$id_col]])
  rownames(expr2_raw) <- as.character(expr2_raw[[info2$id_col]])

  # Extract donor columns only
  donors1 <- colnames(expr1_raw)[-(1:info1$meta_cols)]
  donors2 <- colnames(expr2_raw)[-(1:info2$meta_cols)]

  # Filter to hub features that exist in expression data
  feat1_ids <- as.character(top1$feature_id)
  feat2_ids <- as.character(top2$feature_id)
  feat1_ids <- feat1_ids[feat1_ids %in% rownames(expr1_raw)]
  feat2_ids <- feat2_ids[feat2_ids %in% rownames(expr2_raw)]

  if (length(feat1_ids) == 0 || length(feat2_ids) == 0) {
    return("NO;Hub feature IDs not found in expression data.")
  }

  # Find shared donors (intersect with multi-omics subset if active)
  shared_donors <- intersect(donors1, donors2)
  shared_donors <- .apply_multiomics_subset(shared_donors, donors)
  if (length(shared_donors) < 10) {
    return("NO;Too few shared donors between the two omics types.")
  }

  # Extract expression matrices (features x shared donors)
  mat1 <- as.matrix(expr1_raw[feat1_ids, shared_donors, drop = FALSE])
  mat2 <- as.matrix(expr2_raw[feat2_ids, shared_donors, drop = FALSE])

  # Convert to numeric (safety)
  mat1 <- apply(mat1, 2, as.numeric)
  mat2 <- apply(mat2, 2, as.numeric)
  rownames(mat1) <- feat1_ids
  rownames(mat2) <- feat2_ids

  # --- Compute cross-omics correlations ---
  edges <- list()
  for (i in seq_len(nrow(mat1))) {
    for (j in seq_len(nrow(mat2))) {
      x <- mat1[i, ]
      y <- mat2[j, ]
      valid <- !is.na(x) & !is.na(y)
      if (sum(valid) < 10) next
      ct <- cor.test(x[valid], y[valid], method = "pearson")
      if (ct$p.value < 0.05 && abs(ct$estimate) > 0.3) {
        edges[[length(edges) + 1]] <- list(
          source  = rownames(mat1)[i],
          target  = rownames(mat2)[j],
          weight  = round(ct$estimate, 4),
          pvalue  = round(ct$p.value, 6),
          n_donors = sum(valid)
        )
      }
    }
  }

  # Limit edges: keep top 300 by |correlation|
  if (length(edges) > 300) {
    abs_cors <- sapply(edges, function(e) abs(e$weight))
    keep_idx <- order(-abs_cors)[1:300]
    edges <- edges[keep_idx]
  }

  # --- Build node info ---
  # Collect unique node IDs from edges
  edge_sources <- sapply(edges, function(e) e$source)
  edge_targets <- sapply(edges, function(e) e$target)
  node1_ids <- unique(edge_sources)
  node2_ids <- unique(edge_targets)

  build_nodes <- function(ids, hub_df, summary_raw, omicsType) {
    lapply(ids, function(fid) {
      hi <- hub_df[hub_df$feature_id == fid, ]
      mod_num <- if (nrow(hi) > 0) hi$module[1] else 0
      mod_color <- if (!is.null(summary_raw$modules[[as.character(mod_num)]])) {
        summary_raw$modules[[as.character(mod_num)]]$color
      } else { "grey" }
      list(
        id        = fid,
        label     = if (nrow(hi) > 0 && !is.na(hi$feature_name[1])) hi$feature_name[1] else fid,
        omics     = omicsType,
        module    = mod_num,
        color     = mod_color,
        kME       = if (nrow(hi) > 0) round(hi$kME[1], 4) else 0
      )
    })
  }

  nodes1 <- build_nodes(node1_ids, top1, sum1, omicsType1)
  nodes2 <- build_nodes(node2_ids, top2, sum2, omicsType2)

  net <- list(
    nodes        = c(nodes1, nodes2),
    edges        = edges,
    omics1_name  = omicsType1,
    omics2_name  = omicsType2,
    n_nodes      = length(nodes1) + length(nodes2),
    n_edges      = length(edges),
    n_shared_donors = length(shared_donors),
    n_hub_features1 = length(feat1_ids),
    n_hub_features2 = length(feat2_ids)
  )

  out_file <- "wgcna_crossnet.json"
  jsonlite::write_json(net, out_file, auto_unbox = TRUE, pretty = TRUE)

  msg <- sprintf(
    "Cross-omics network: %d nodes (%d %s + %d %s), %d significant edges (|r|>0.3, p<0.05, %d shared donors).",
    net$n_nodes, length(nodes1), omicsType1, length(nodes2), omicsType2,
    length(edges), length(shared_donors)
  )

  return(paste0("RES-OK;", out_file, ";", msg))
}


################################################################################
# ------------------------------------------------------------------------------
# Fast targeted expression read for very large omics matrices (e.g. CpG-level
# methylation, ~729k rows / 678 MB): pull ONLY the requested feature rows via grep
# instead of loading the whole file just to correlate the handful of top-N
# features. Safe only for distinctive column-1 feature ids (CpG cg########).
# ------------------------------------------------------------------------------
.circos_read_expr_subset <- function(fpath, ids) {
  ids <- unique(as.character(ids)); ids <- ids[nzchar(ids)]
  if (length(ids) == 0) return(NULL)
  hdr <- readLines(fpath, n = 1)
  idf <- tempfile(); writeLines(ids, idf); on.exit(unlink(idf), add = TRUE)
  rows <- tryCatch(system2("grep", c("-F", "-w", "-f", shQuote(idf), shQuote(fpath)),
                           stdout = TRUE), error = function(e) character(0))
  if (length(rows) == 0) return(NULL)
  read.csv(text = c(hdr, rows), stringsAsFactors = FALSE, check.names = FALSE)
}

# MOFA cross-omics CIRCOS data (real correlations, not the w1*w2 co-loading).
# For one factor: take the top-N features per omics by |weight| on that factor,
# load their raw expression, and compute the ACTUAL pairwise cross-omics Pearson
# correlation (+ p + BH FDR) across shared donors -- exactly like build_diablo_circos.
# Returns JSON { factor, nodes:[{id,type,label}], edges:[{source,target,weight,pval,padj,sign,n}] }
# consumed by the frontend circos renderer (drawFeatureCircos). Served via
# DoPrecomputed(step="mofa_circos"): traits=factor, category=topN.
################################################################################
GetMofaCircos <- function(base_path, combo, factor = "Factor1", topN = "15", donors = "all") {
  topN <- suppressWarnings(as.integer(topN)); if (is.na(topN) || topN < 3) topN <- 15
  mofa_dir <- file.path(base_path, "mofa_results", combo)
  if (!dir.exists(mofa_dir))
    return(jsonlite::toJSON(list(error = paste0("MOFA results not found: ", combo)), auto_unbox = TRUE))
  summ_path <- file.path(mofa_dir, "summary.json")
  if (!file.exists(summ_path))
    return(jsonlite::toJSON(list(error = "MOFA summary.json not found"), auto_unbox = TRUE))
  omics_types <- jsonlite::fromJSON(summ_path)$omics_types
  if (is.null(omics_types) || length(omics_types) < 2)
    return(jsonlite::toJSON(list(error = "Fewer than two omics in this MOFA model"), auto_unbox = TRUE))

  expr_dir <- paste0(base_path, "omics/")
  get_expr_info <- function(omicsType) {
    if (omicsType == "rnaseq") list(file = "rnaseq_filtered.csv", id_col = "accession", meta_cols = 5)
    else if (omicsType == "proteomics") list(file = "proteomics_v2.csv", id_col = "symbol", meta_cols = 4, fallback_col = "Protein_Group")
    else if (omicsType == "nanostring") list(file = "nanostring.csv", id_col = "symbol", meta_cols = 4)
    else if (omicsType == "methylation") list(file = "methylation_M.csv", id_col = "feature_id", meta_cols = 8)
    else if (omicsType %in% c("pseudobulk_alpha", "pseudobulk_beta", "patchseq_alpha", "patchseq_beta")) {
      cell <- sub("^(pseudobulk|patchseq)_", "", omicsType)
      list(file = paste0("pseudobulk_", cell, ".csv"), id_col = "gene_id", meta_cols = 3)
    }
    else if (omicsType %in% c("flux_HG", "flux_LG", "flux_ratio")) list(file = paste0(omicsType, ".csv"), id_col = "rxn", meta_cols = 4)
    else NULL
  }

  # Per-omics: top-N features on this factor (by |weight|) + their expression rows.
  per_omics <- list()
  for (om in omics_types) {
    wpath <- file.path(mofa_dir, paste0("weights_", om, ".csv"))
    info  <- get_expr_info(om)
    if (!file.exists(wpath) || is.null(info) || !file.exists(paste0(expr_dir, info$file))) next
    w <- read.csv(wpath, stringsAsFactors = FALSE, check.names = FALSE)
    if (!(factor %in% colnames(w))) next
    w$.wt <- suppressWarnings(as.numeric(w[[factor]]))
    w <- w[!is.na(w$.wt), , drop = FALSE]
    w <- w[order(-abs(w$.wt)), , drop = FALSE]
    top <- utils::head(w, topN)
    if (nrow(top) == 0) next
    # Big matrices (CpG-level methylation) read only the top-N rows; others read whole.
    expr <- if (identical(om, "methylation"))
              .circos_read_expr_subset(paste0(expr_dir, info$file), top$feature_id)
            else
              read.csv(paste0(expr_dir, info$file), stringsAsFactors = FALSE, check.names = FALSE)
    if (is.null(expr) || nrow(expr) == 0) next
    # Rebuild the feature id EXACTLY as the pipeline built weights$feature_id, so
    # rownames are unique AND line up with the weights (no crash, nothing dropped):
    #   1. blank id -> fallback_col  (proteomics_v2: 2 rows have no symbol -> Protein_Group)
    #   2. make.unique on the same file/row order (proteomics_v2: 7 symbols shared by
    #      two protein groups -> "SYM", "SYM_dup1")
    expr_ids <- as.character(expr[[info$id_col]])
    if (!is.null(info$fallback_col) && info$fallback_col %in% colnames(expr)) {
      blank <- is.na(expr_ids) | expr_ids == "" | expr_ids == "NA"   # same test as the pipeline
      if (any(blank)) expr_ids[blank] <- as.character(expr[[info$fallback_col]])[blank]
    }
    if (any(duplicated(expr_ids))) expr_ids <- make.unique(expr_ids, sep = "_dup")
    rownames(expr) <- expr_ids
    dons <- colnames(expr)[-(1:info$meta_cols)]
    fids <- as.character(top$feature_id); fids <- fids[fids %in% rownames(expr)]
    if (length(fids) == 0) next
    labs <- setNames(as.character(top$feature_name), as.character(top$feature_id))
    per_omics[[om]] <- list(expr = expr, dons = dons, fids = fids, labs = labs)
  }
  oms <- names(per_omics)
  if (length(oms) < 2)
    return(jsonlite::toJSON(list(error = "Not enough omics with weights + expression for this factor"), auto_unbox = TRUE))

  # All cross-omics feature pairs -> real Pearson r on shared donors.
  edges <- list(); nodes_seen <- list()
  for (ai in seq_len(length(oms) - 1)) for (bi in (ai + 1):length(oms)) {
    A <- per_omics[[oms[ai]]]; B <- per_omics[[oms[bi]]]
    shared <- .apply_multiomics_subset(intersect(A$dons, B$dons), donors)
    if (length(shared) < 10) next
    matA <- apply(as.matrix(A$expr[A$fids, shared, drop = FALSE]), 2, as.numeric); rownames(matA) <- A$fids
    matB <- apply(as.matrix(B$expr[B$fids, shared, drop = FALSE]), 2, as.numeric); rownames(matB) <- B$fids
    for (i in seq_len(nrow(matA))) for (j in seq_len(nrow(matB))) {
      x <- matA[i, ]; y <- matB[j, ]; ok <- !is.na(x) & !is.na(y)
      if (sum(ok) < 10 || sd(x[ok]) == 0 || sd(y[ok]) == 0) next
      ct <- try(cor.test(x[ok], y[ok], method = "pearson"), silent = TRUE)
      if (inherits(ct, "try-error")) next
      sid <- paste0(oms[ai], "::", A$fids[i]); tid <- paste0(oms[bi], "::", B$fids[j])
      edges[[length(edges) + 1]] <- list(source = sid, target = tid,
        weight = round(unname(ct$estimate), 4), pval = ct$p.value, n = sum(ok),
        sign = if (unname(ct$estimate) >= 0) "pos" else "neg")
      nodes_seen[[sid]] <- list(id = sid, type = oms[ai], label = unname(A$labs[A$fids[i]]))
      nodes_seen[[tid]] <- list(id = tid, type = oms[bi], label = unname(B$labs[B$fids[j]]))
    }
  }
  if (length(edges) == 0)
    return(jsonlite::toJSON(list(factor = factor, nodes = list(), edges = list(),
                                 message = "No cross-omics correlations for this factor."), auto_unbox = TRUE))

  # BH FDR across every returned pair (the frontend slider filters |r| on top).
  padj <- p.adjust(vapply(edges, function(e) e$pval, numeric(1)), method = "BH")
  for (k in seq_along(edges)) { edges[[k]]$pval <- signif(edges[[k]]$pval, 4); edges[[k]]$padj <- signif(padj[k], 4) }

  jsonlite::toJSON(list(factor = factor, nodes = unname(nodes_seen), edges = edges),
                   auto_unbox = TRUE, na = "null", digits = 6)
}


################################################################################
# Precomputed WGCNA / MOFA pathway enrichment
# ------------------------------------------------------------------------------
# Routed through the SAME ora/gsea Java step DIABLO uses (see the early dispatch
# added to GetORAMultiOmics / GetGSEAMultiOmics). The frontend encodes the source
# and the feature selection in the `omics_for_pathway` parameter (comma-delimited;
# every character passes the servlet whitelist) and leaves circos_features empty:
#     WGCNA : "PCWGCNA,<combo>,<moduleNumber>"
#     MOFA  : "PCMOFA,<combo>,<factor>,<topN>"
# We build dea_results.csv in the exact schema performORA()/performGSEA() expect
# (Gene_ID = ENTREZ, column 2 = signed score) from the precomputed
# wgcna_combined_results/<combo>/modules.csv and mofa_results/<combo>/weights_*.csv,
# then run the shared enrichment engine -- no new Java endpoint, no engine change.
################################################################################

# Gene omics whose feature ids resolve to Entrez. flux / metabolite are a
# different pathway domain (GEM reaction / KEGG metabolite) and are excluded here.
.pc_gene_omics <- c("rnaseq", "proteomics", "nanostring",
                    "pseudobulk_alpha", "pseudobulk_beta",
                    "patchseq_alpha", "patchseq_beta")

# Vectorised feature-id -> Entrez, following the WORKING DIABLO convention:
#   rnaseq / nanostring / pseudobulk : feature id is an ENSG accession -> Entrez
#   proteomics                       : feature id IS the Entrez gene id
#   fallback (any omics)             : HGNC symbol -> Entrez
.pc_id_to_entrez <- function(raw, omics, sym) {
  ga <- qs::qread(paste0(other.tables.path, "libraries/gene_anno_full.qs"))
  ens2ent <- setNames(as.character(ga$entrez_selected), as.character(ga$ensembl_gene_id))
  sym2ent <- setNames(as.character(ga$entrez_selected), as.character(ga$hgnc_symbol))
  ent_all <- unique(as.character(ga$entrez_selected))
  raw <- as.character(raw); omics <- as.character(omics); sym <- as.character(sym)
  e <- rep(NA_character_, length(raw))
  is_ens <- grepl("^ENSG", raw)
  e[is_ens] <- ens2ent[raw[is_ens]]
  is_prot <- is.na(e) & omics == "proteomics" & raw %in% ent_all
  e[is_prot] <- raw[is_prot]
  need <- (is.na(e) | e == "") & !is.na(sym) & sym != ""
  e[need] <- sym2ent[sym[need]]
  unname(e)
}

# WGCNA selection: one module's members (foreground) + all analysed features of
# that module's omics in the combo (background = the standard WGCNA universe).
.pc_wgcna_selection <- function(base_path, combo, moduleNum) {
  mfile <- file.path(base_path, "wgcna_combined_results", combo, "modules.csv")
  if (!file.exists(mfile)) return(list(err = paste0("WGCNA modules.csv not found for combo: ", combo)))
  mods <- read.csv(mfile, stringsAsFactors = FALSE)
  modn <- suppressWarnings(as.integer(moduleNum))
  if (is.na(modn)) return(list(err = paste0("Invalid module: ", moduleNum)))
  hit <- mods[mods$module_number == modn, , drop = FALSE]
  if (nrow(hit) == 0) return(list(err = paste0("No features in module ", modn)))
  om <- unique(hit$omics_source); om <- om[om %in% .pc_gene_omics]
  if (length(om) == 0)
    return(list(err = paste0("Module M", modn, " is a non-gene omics module; pathway enrichment applies to gene omics (RNA/protein) only.")))
  bg <- mods[mods$omics_source %in% om, , drop = FALSE]
  list(hit_raw = hit$original_id, hit_om = hit$omics_source, hit_sym = hit$feature_name, hit_score = NULL,
       bg_raw  = bg$original_id,  bg_om  = bg$omics_source,  bg_sym  = bg$feature_name,  bg_score = NULL)
}

# MOFA selection: top-N features per gene omics by |weight| on the factor
# (foreground) + every weighted gene feature (background). The signed factor
# weight is carried as the per-gene score so GSEA can rank on it.
.pc_mofa_selection <- function(base_path, combo, factor, topN) {
  mdir <- file.path(base_path, "mofa_results", combo)
  sfile <- file.path(mdir, "summary.json")
  if (!file.exists(sfile)) return(list(err = paste0("MOFA summary.json not found for combo: ", combo)))
  oms <- jsonlite::fromJSON(sfile)$omics_types
  oms <- oms[oms %in% .pc_gene_omics]
  if (length(oms) == 0) return(list(err = "No gene omics in this MOFA model for pathway enrichment."))
  hr <- c(); ho <- c(); hs <- c(); hc <- c(); br <- c(); bo <- c(); bs <- c(); bc <- c()
  for (om in oms) {
    wpath <- file.path(mdir, paste0("weights_", om, ".csv"))
    if (!file.exists(wpath)) next
    w <- read.csv(wpath, stringsAsFactors = FALSE, check.names = FALSE)
    if (!(factor %in% colnames(w)) || !("feature_id" %in% colnames(w))) next
    w$.wt <- suppressWarnings(as.numeric(w[[factor]]))
    w <- w[!is.na(w$.wt), , drop = FALSE]
    if (nrow(w) == 0) next
    sym <- if ("feature_name" %in% colnames(w)) w$feature_name else w$feature_id
    br <- c(br, w$feature_id); bo <- c(bo, rep(om, nrow(w))); bs <- c(bs, sym); bc <- c(bc, w$.wt)
    top <- utils::head(w[order(-abs(w$.wt)), , drop = FALSE], topN)
    tsym <- if ("feature_name" %in% colnames(top)) top$feature_name else top$feature_id
    hr <- c(hr, top$feature_id); ho <- c(ho, rep(om, nrow(top))); hs <- c(hs, tsym); hc <- c(hc, top$.wt)
  }
  if (length(hr) == 0) return(list(err = paste0("Factor ", factor, " has no gene-omics weights.")))
  list(hit_raw = hr, hit_om = ho, hit_sym = hs, hit_score = hc,
       bg_raw = br, bg_om = bo, bg_sym = bs, bg_score = bc)
}

# Write dea_results.csv (Gene_ID = Entrez) from a selection; return hit/universe counts.
.pc_write_pathway_dea <- function(sel) {
  ga <- qs::qread(paste0(other.tables.path, "libraries/gene_anno_full.qs"))
  bg_ent  <- .pc_id_to_entrez(sel$bg_raw,  sel$bg_om,  sel$bg_sym)
  hit_ent <- .pc_id_to_entrez(sel$hit_raw, sel$hit_om, sel$hit_sym)
  score <- if (!is.null(sel$bg_score) && length(sel$bg_score) == length(bg_ent))
    suppressWarnings(as.numeric(sel$bg_score)) else rep(0, length(bg_ent))
  score[is.na(score)] <- 0
  keep <- !is.na(bg_ent) & bg_ent != ""
  bg_ent <- bg_ent[keep]; score <- score[keep]
  ord <- order(-abs(score)); bg_ent <- bg_ent[ord]; score <- score[ord]   # dedup keeps max |score|
  dup <- duplicated(bg_ent); bg_ent <- bg_ent[!dup]; score <- score[!dup]
  hit_set <- unique(hit_ent[!is.na(hit_ent) & hit_ent != ""])
  if (length(bg_ent) < 10 || length(hit_set) < 2)
    return(list(hits = length(hit_set), universe = length(bg_ent), ok = FALSE))
  sym  <- ga$hgnc_symbol[match(bg_ent, as.character(ga$entrez_selected))]
  desc <- ga$description[match(bg_ent, as.character(ga$entrez_selected))]
  dea <- data.frame(
    Feature = ifelse(is.na(sym) | sym == "", bg_ent, sym),
    Log2FC  = score, t.stat = score, Gene_ID = bg_ent,
    Description = ifelse(is.na(desc), "", desc),
    sig = ifelse(bg_ent %in% hit_set, "sig", "NS"),
    stringsAsFactors = FALSE
  )
  write.csv(dea, "dea_results.csv", row.names = FALSE)
  list(hits = sum(dea$sig == "sig"), universe = nrow(dea), ok = TRUE)
}

# Entry point invoked by GetORAMultiOmics / GetGSEAMultiOmics when the routing
# token is present. analysis = "ora" | "gsea".
.pc_run_pathway <- function(routing, library, fdr, collapse, mode, analysis, rank = "coef") {
  base_path <- paste0(other.tables.path, "processed_comprehensive/")
  parts <- strsplit(routing, ",", fixed = TRUE)[[1]]
  src <- parts[1]
  sel <- if (src == "PCWGCNA") {
    .pc_wgcna_selection(base_path, parts[2], parts[3])
  } else if (src == "PCMOFA") {
    tN <- suppressWarnings(as.integer(parts[4])); if (is.na(tN) || tN < 5) tN <- 100
    .pc_mofa_selection(base_path, parts[2], parts[3], tN)
  } else list(err = paste0("Unknown pathway source: ", src))
  if (!is.null(sel$err)) return(paste0("RES-NO; ", sel$err))
  info <- .pc_write_pathway_dea(sel)
  if (!isTRUE(info$ok))
    return(paste0("RES-NO; too few genes mapped for enrichment (foreground = ",
                  info$hits, ", background = ", info$universe, ")."))
  if (analysis == "ora") {
    rcmd <<- paste0('performORA("', library, '", "', fdr, '", "', collapse, '", "', mode, '")')
    performORA(library, fdr, collapse, mode)
  } else {
    rcmd <<- paste0('performGSEA("', library, '", "', rank, '", "', fdr, '", "', collapse, '", "', mode, '")')
    performGSEA(library, rank, fdr, collapse, mode)
  }
}


################################################################################
# WGCNA Pathway Enrichment Analysis
# Runs ORA on features from a specific WGCNA module
# Input: omicsType, moduleNum, library, fdr, collapse
# Writes dea_results.csv in the format expected by performORA()/performGSEA()
# Then calls performORA() or performGSEA()
################################################################################
GetWGCNAPathway <- function(omicsType, moduleNum, library, fdr, collapse, analysisType = "ora", donors = "all") {
  # Module membership is a fixed property of the precomputed WGCNA model; pathway
  # enrichment runs over the module's feature list, so donor subset doesn't apply here.

  wgcna_dir <- paste0(other.tables.path, "multi_omics/wgcna_results/", omicsType, "/")
  if (!dir.exists(wgcna_dir)) return(paste0("NO;WGCNA results not found for: ", omicsType))

  # Read module assignments and hub features
  modules_df <- read.csv(paste0(wgcna_dir, "modules.csv"), stringsAsFactors = FALSE)
  hub_df     <- read.csv(paste0(wgcna_dir, "hub_features.csv"), stringsAsFactors = FALSE)
  sum_json   <- jsonlite::fromJSON(paste0(wgcna_dir, "summary.json"))

  mod_num <- as.integer(moduleNum)

  # Get module color for the message
  mod_info <- sum_json$modules[[as.character(mod_num)]]
  mod_color <- if (!is.null(mod_info)) mod_info$color else "unknown"
  mod_size  <- if (!is.null(mod_info)) mod_info$size else sum(modules_df$module_number == mod_num)

  # Get features in this module (hits)
  hit_features <- modules_df$feature_id[modules_df$module_number == mod_num]
  # All features across all modules (background/universe)
  all_features <- modules_df$feature_id

  if (length(hit_features) == 0) {
    return(paste0("NO;No features found in module ", mod_num))
  }

  # Get kME values for ranking (for GSEA)
  hub_kme <- hub_df[hub_df$module == mod_num, c("feature_id", "kME")]
  kme_lookup <- setNames(hub_kme$kME, hub_kme$feature_id)

  # Load gene annotation for ID conversion
  gene_anno_file <- paste0(other.tables.path, "libraries/gene_anno_full.qs")
  if (!file.exists(gene_anno_file)) {
    return("NO;Gene annotation file not found.")
  }
  gene_anno <- qs::qread(gene_anno_file)

  # Convert feature IDs to Entrez IDs
  # For rnaseq: feature_id = ENSG (accession); for proteomics: feature_id = numeric gene_id
  id_to_entrez <- setNames(as.character(gene_anno$entrezgene_id), as.character(gene_anno$ensembl_gene_id))

  # Try direct mapping first (for ENSG IDs)
  hit_entrez <- na.omit(id_to_entrez[hit_features])
  all_entrez <- na.omit(id_to_entrez[all_features])

  # If direct mapping fails (proteomics uses numeric gene_id), try alternative
  if (length(hit_entrez) < 2) {
    # Try symbol-based mapping
    feat_names <- modules_df$feature_name[modules_df$module_number == mod_num]
    sym_to_entrez <- setNames(as.character(gene_anno$entrezgene_id), as.character(gene_anno$hgnc_symbol))
    hit_entrez <- na.omit(sym_to_entrez[feat_names])
    all_names <- modules_df$feature_name
    all_entrez <- na.omit(sym_to_entrez[all_names])
  }

  if (length(hit_entrez) < 2) {
    return(paste0("NO;Could not map enough features to gene IDs. Only ", length(hit_entrez), " mapped."))
  }

  # Build dea_results.csv in the format performORA/performGSEA expects
  dea_mat <- data.frame(
    Feature = all_entrez,
    Log2FC  = 0,
    t.stat  = 0,
    Gene_ID = names(all_entrez),
    Description = "",
    sig = "NS",
    stringsAsFactors = FALSE
  )

  # Mark hit features
  dea_mat$sig[dea_mat$Feature %in% hit_entrez] <- "sig"

  # Set kME as Log2FC for ranking (GSEA uses this)
  for (i in seq_len(nrow(dea_mat))) {
    fid <- dea_mat$Gene_ID[i]
    if (fid %in% names(kme_lookup)) {
      dea_mat$Log2FC[i] <- kme_lookup[fid]
      dea_mat$t.stat[i] <- kme_lookup[fid]
    }
  }

  # For non-module features, set small random values for GSEA ranking
  ns_idx <- dea_mat$sig == "NS"
  dea_mat$Log2FC[ns_idx] <- runif(sum(ns_idx), -0.01, 0.01)
  dea_mat$t.stat[ns_idx] <- dea_mat$Log2FC[ns_idx]

  write.csv(dea_mat, "dea_results.csv", row.names = FALSE)

  # Run the analysis
  if (analysisType == "ora") {
    res <- performORA(funcLib = library, fdr = as.numeric(fdr), collapse = collapse, mode = "tool")
  } else {
    res <- performGSEA(funcLib = library, rank.stat = "coef", fdr = as.numeric(fdr), collapse = collapse, mode = "tool")
  }

  # Prepend module info to the result message
  if (grepl("^RES-OK", res)) {
    parts <- unlist(strsplit(res, ";"))
    n_sig <- parts[2]
    msg <- sprintf("Module M%d (%s, %d features from %s): %s significant pathways found.",
                   mod_num, mod_color, mod_size, omicsType, n_sig)
    return(paste0("RES-OK;", n_sig, ";", msg))
  }
  return(res)
}


################################################################################
GetWGCNAEigengeneBox <- function(omicsType, moduleNum, analysisVar, class1, class2, primaryType, donors = "all") {

  wgcna_dir <- paste0(other.tables.path, "multi_omics/wgcna_results/", omicsType, "/")
  eigen_df  <- read.csv(paste0(wgcna_dir, "eigengenes.csv"), row.names = 1, stringsAsFactors = FALSE)
  meta      <- read.csv(paste0(other.tables.path, "display_data/metadata_sum_norm.csv"), stringsAsFactors = FALSE)

  # Inject cluster (or any proc_metadata-only) column if selected var is missing from the CSV
  if (!is.null(analysisVar) && !(analysisVar %in% colnames(meta))) {
    library(RSQLite)
    meta_db <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
    pheno <- dbReadTable(meta_db, "proc_metadata")
    dbDisconnect(meta_db)
    if (analysisVar %in% colnames(pheno)) {
      extra <- pheno[, c("record_id", analysisVar)]
      meta <- merge(meta, extra, by = "record_id", all.x = TRUE)
    } else {
      return(paste0("NO;Variable '", analysisVar,
                    "' not found in metadata or donor table."))
    }
  }

  me_col <- paste0("ME", moduleNum)
  if (!me_col %in% colnames(eigen_df)) {
    return(paste0("NO;Module ", moduleNum, " not found in eigengenes."))
  }

  shared <- intersect(rownames(eigen_df), meta$record_id)
  shared <- .apply_multiomics_subset(shared, donors)
  eigen_sub <- eigen_df[shared, me_col, drop = FALSE]
  meta_sub  <- meta[match(shared, meta$record_id), ]

  if (primaryType == "disc") {
    keep <- meta_sub[[analysisVar]] %in% c(class1, class2)
    meta_sub  <- meta_sub[keep, ]
    eigen_sub <- eigen_sub[keep, , drop = FALSE]
  }

  df_out <- data.frame(
    record  = meta_sub$record_id,
    meta    = meta_sub[[analysisVar]],
    feature = eigen_sub[[me_col]]
  )

  write.csv(df_out, "df.csv", row.names = FALSE)
  return("RES-OK")
}


################################################################################
# Precomputed multi-omics dispatcher
# Called from PrecomputedMultiomics.java via RCenter.doPrecomputed()
# Handles all precomputed WGCNA/MOFA data serving and computation
################################################################################
DoPrecomputed <- function(step, method="NULL", combo="NULL", file="NULL",
                           traits="NULL", category="NULL", primaryType="NULL",
                           class1="NULL", class2="NULL", donors="NULL") {

  base_path <- paste0(other.tables.path, "processed_comprehensive/")

  if (step == "combos") {
    return(.pc_list_combos(base_path))
  } else if (step == "wgcna_file") {
    return(.pc_serve_file(base_path, "wgcna_combined_results", combo, file))
  } else if (step == "mofa_file") {
    return(.pc_serve_file(base_path, "mofa_results", combo, file))
  } else if (step == "phenotype") {
    return(.pc_serve_phenotype(base_path, category, donors))
  } else if (step == "traitcorr") {
    return(.pc_trait_corr(base_path, method, combo, traits, category, primaryType, class1, class2, donors))
  } else if (step == "traitcorr_multigroup") {
    return(.pc_multigroup_corr(base_path, method, combo, traits, category, donors))
  } else if (step == "overview_heatmap") {
    return(.pc_overview_heatmap(base_path, method, combo, category, traits, donors))
  } else if (step == "mofa_circos") {
    # traits = factor (e.g. "Factor1"), category = topN features/omics
    return(GetMofaCircos(base_path, combo, traits, category, donors))
  } else {
    return(jsonlite::toJSON(list(error = paste0("Unknown step: ", step)), auto_unbox = TRUE))
  }
}

# ---- Helper: list all combo analyses ----
.pc_list_combos <- function(base_path) {
  wgcna_dir <- file.path(base_path, "wgcna_combined_results")
  mofa_dir <- file.path(base_path, "mofa_results")

  combos <- unique(c(
    if (dir.exists(wgcna_dir)) list.dirs(wgcna_dir, recursive = FALSE, full.names = FALSE) else character(0),
    if (dir.exists(mofa_dir)) list.dirs(mofa_dir, recursive = FALSE, full.names = FALSE) else character(0)
  ))

  result <- lapply(combos, function(combo) {
    entry <- list(combo_name = combo)
    # WGCNA summary
    ws_path <- file.path(wgcna_dir, combo, "summary.json")
    if (file.exists(ws_path)) {
      ws <- jsonlite::fromJSON(ws_path)
      entry$has_wgcna <- TRUE
      entry$wgcna_n_modules <- ws$n_modules
      entry$wgcna_n_donors <- ws$n_donors
      entry$wgcna_n_features <- ws$n_features
      entry$wgcna_sft_quality <- ws$sft_quality
      entry$wgcna_sft_r2 <- ws$sft_r2
      entry$wgcna_power <- ws$power
    } else {
      entry$has_wgcna <- FALSE
    }
    # MOFA summary
    ms_path <- file.path(mofa_dir, combo, "summary.json")
    if (file.exists(ms_path)) {
      ms <- jsonlite::fromJSON(ms_path)
      entry$has_mofa <- TRUE
      entry$mofa_n_donors <- ms$n_donors_total
      entry$mofa_n_factors <- ms$n_factors_active
    } else {
      entry$has_mofa <- FALSE
    }
    entry
  })

  return(jsonlite::toJSON(result, auto_unbox = TRUE))
}

# ---- Helper: serve a file as text ----
.pc_serve_file <- function(base_path, subdir, combo, file) {
  fpath <- file.path(base_path, subdir, combo, file)
  if (!file.exists(fpath)) {
    return(jsonlite::toJSON(list(error = paste0("File not found: ", combo, "/", file)), auto_unbox = TRUE))
  }
  return(paste(readLines(fpath, warn = FALSE), collapse = "\n"))
}

# ---- Helper: serve phenotype CSV ----
.pc_serve_phenotype <- function(base_path, category, donors) {
  if (category == "NULL" || category == "") category <- "core"
  pheno_file <- switch(category,
    "all" = "all_features.csv", "isolation" = "features_isolation.csv",
    "gsis" = "features_gsis.csv", "electrophysiology" = "features_electrophysiology.csv",
    "grs" = "features_grs.csv", "cell" = "features_cell.csv",
    "seahorse" = "features_seahorse.csv", "perifusion" = "features_perifusion.csv",
    "lip_extract" = "features_lip_extract.csv", "features_core.csv")
  fpath <- file.path(base_path, "phenotype", pheno_file)
  if (!file.exists(fpath)) {
    return(jsonlite::toJSON(list(error = paste0("Phenotype file not found: ", pheno_file)), auto_unbox = TRUE))
  }
  lines <- readLines(fpath, warn = FALSE)
  if (donors == "NULL" || donors == "" || donors == "all") {
    return(paste(lines, collapse = "\n"))
  }
  # Filter by donor
  donor_set <- strsplit(donors, ",")[[1]]
  header <- lines[1]
  filtered <- lines[-1][sapply(lines[-1], function(l) {
    sid <- gsub('"', '', strsplit(l, ",")[[1]][1])
    sid %in% donor_set
  })]
  return(paste(c(header, filtered), collapse = "\n"))
}

# ---- Helper: load scores (eigengenes or factors) and phenotype, return merged ----
# `donors` is "all" (default) or "subset" — when "subset", the result is restricted
# to donors listed in donors_multiomics.rds (which is the multi-omics-page-only
# subset file, isolated from omics/phenotype views' donors.rds).
.pc_load_scores_pheno <- function(base_path, method, combo, traits, category, donors = "all") {
  if (method == "wgcna") {
    scores_path <- file.path(base_path, "wgcna_combined_results", combo, "eigengenes.csv")
  } else {
    scores_path <- file.path(base_path, "mofa_results", combo, "factors.csv")
  }
  if (!file.exists(scores_path)) return(NULL)

  scores <- read.csv(scores_path, stringsAsFactors = FALSE, check.names = FALSE)

  if (method == "wgcna") {
    donor_ids <- scores[[1]]
    score_cols <- colnames(scores)[-1]
    scores_mat <- scores[, -1, drop = FALSE]
  } else {
    donor_ids <- scores$donor
    score_cols <- setdiff(colnames(scores), "donor")
    scores_mat <- scores[, score_cols, drop = FALSE]
  }

  if (category == "NULL" || category == "") category <- "core"
  pheno_file <- switch(category,
    "all" = "all_features.csv", "isolation" = "features_isolation.csv",
    "gsis" = "features_gsis.csv", "electrophysiology" = "features_electrophysiology.csv",
    "grs" = "features_grs.csv", "cell" = "features_cell.csv",
    "seahorse" = "features_seahorse.csv", "perifusion" = "features_perifusion.csv",
    "lip_extract" = "features_lip_extract.csv", "features_core.csv")
  pheno_path <- file.path(base_path, "phenotype", pheno_file)
  if (!file.exists(pheno_path)) return(NULL)

  pheno <- read.csv(pheno_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (traits == "NULL" || traits == "") traits <- "donorage"
  trait <- strsplit(traits, ",")[[1]][1]

  common <- intersect(donor_ids, pheno$sample_id)
  common <- .apply_multiomics_subset(common, donors)
  if (length(common) < 5) return(NULL)

  idx_s <- match(common, donor_ids)
  idx_p <- match(common, pheno$sample_id)

  list(scores_mat = scores_mat[idx_s, , drop = FALSE],
       score_cols = score_cols,
       trait_vals = pheno[[trait]][idx_p],
       trait = trait,
       n = length(common))
}

# ---- Trait correlation (continuous or binary categorical) ----
.pc_trait_corr <- function(base_path, method, combo, traits, category, primaryType, class1, class2, donors = "all") {
  d <- .pc_load_scores_pheno(base_path, method, combo, traits, category, donors)
  if (is.null(d)) return(jsonlite::toJSON(list(correlations = list(), n_donors = 0), auto_unbox = TRUE))

  # Usable-donor count for this trait/comparison. Reported so the UI can explain
  # an empty result ("only N donors...") instead of showing a blank plot; every
  # factor is skipped below when its usable n < 5.
  if (primaryType == "disc" && class1 != "NULL" && class2 != "NULL") {
    n_donors <- sum(as.character(d$trait_vals) %in% c(class1, class2), na.rm = TRUE)
  } else {
    n_donors <- sum(!is.na(suppressWarnings(as.numeric(d$trait_vals))))
  }

  correlations <- list()

  if (primaryType == "cont" || primaryType == "NULL") {
    trait_num <- suppressWarnings(as.numeric(d$trait_vals))
    valid <- !is.na(trait_num)
    for (comp in d$score_cols) {
      comp_vals <- as.numeric(d$scores_mat[[comp]])
      use <- valid & !is.na(comp_vals)
      n <- sum(use)
      if (n < 5) next
      ct <- cor.test(comp_vals[use], trait_num[use], method = "pearson")
      correlations[[length(correlations) + 1]] <- list(
        component = comp, trait = d$trait,
        r = round(ct$estimate, 4), p = ct$p.value, n = n)
    }
  } else if (primaryType == "disc" && class1 != "NULL" && class2 != "NULL") {
    cat_vals <- as.character(d$trait_vals)
    for (comp in d$score_cols) {
      comp_vals <- as.numeric(d$scores_mat[[comp]])
      use <- !is.na(comp_vals) & (cat_vals == class1 | cat_vals == class2)
      n <- sum(use)
      if (n < 5) next
      x <- comp_vals[use]
      y <- ifelse(cat_vals[use] == class1, 1.0, 0.0)
      ct <- cor.test(x, y, method = "pearson")
      correlations[[length(correlations) + 1]] <- list(
        component = comp, trait = d$trait,
        r = round(ct$estimate, 4), p = ct$p.value, n = n,
        class1 = class1, class2 = class2,
        n1 = sum(y == 1), n2 = sum(y == 0))
    }
  }

  return(jsonlite::toJSON(list(method = method, combo = combo,
                                primaryType = primaryType,
                                n_donors = n_donors,
                                correlations = correlations), auto_unbox = TRUE))
}

# ---- Multi-group categorical correlation (Kruskal-Wallis + Cohen's d) ----
.pc_multigroup_corr <- function(base_path, method, combo, traits, category, donors = "all") {
  d <- .pc_load_scores_pheno(base_path, method, combo, traits, category, donors)
  if (is.null(d)) return(jsonlite::toJSON(list(correlations = list()), auto_unbox = TRUE))

  cat_vals <- as.character(d$trait_vals)
  valid_cat <- !is.na(cat_vals) & cat_vals != "" & cat_vals != "NA"

  correlations <- list()

  for (comp in d$score_cols) {
    comp_vals <- as.numeric(d$scores_mat[[comp]])
    use <- valid_cat & !is.na(comp_vals)
    groups_factor <- factor(cat_vals[use])
    if (nlevels(groups_factor) < 2) next
    if (sum(use) < 5) next

    vals <- comp_vals[use]
    overall_mean <- mean(vals)
    overall_sd <- sd(vals)
    if (overall_sd == 0) next

    kw <- tryCatch(kruskal.test(vals ~ groups_factor),
                   error = function(e) list(statistic = NA, p.value = 1))

    group_stats <- lapply(levels(groups_factor), function(g) {
      g_vals <- vals[groups_factor == g]
      list(group = g, mean = round(mean(g_vals), 4),
           n = length(g_vals),
           cohen_d = round((mean(g_vals) - overall_mean) / overall_sd, 4))
    })

    correlations[[length(correlations) + 1]] <- list(
      component = comp, trait = d$trait,
      overall_mean = round(overall_mean, 4),
      overall_sd = round(overall_sd, 4),
      kw_p = kw$p.value,
      kw_h = round(as.numeric(kw$statistic), 2),
      groups = group_stats)
  }

  return(jsonlite::toJSON(list(method = method, combo = combo,
                                mode = "multigroup",
                                correlations = correlations), auto_unbox = TRUE))
}

# ---- Overview heatmap: factors/modules x phenotype-category, matrix form ----
# Reads the pre-computed phenotype-correlation file written by the pipeline
# (factor_phenotype_correlations.csv for MOFA, module_phenotype_correlations.csv
# for WGCNA), filters to a phenotype category (default = "core"), picks the
# preferred correlation method per phenotype, and reshapes to a matrix the
# frontend can render directly as a heatmap.
#
# Method selection per phenotype type:
#   numeric            -> bicor (WGCNA) / pearson (MOFA)
#   binary_categorical -> pearson (point-biserial; bicor isn't run on 0/1)
#   multilevel_cat     -> anova_eta (sqrt(eta^2), unsigned 0-1)
#
# Output JSON shape:
#   { method, combo, category,
#     rows: ["Factor1", ...],      cols: ["donorage", ...],
#     col_types: ["numeric", "binary_categorical", "multilevel_categorical", ...],
#     col_methods: ["pearson", "pearson", "anova_eta", ...],
#     statistic, p_value, p_adj, n   (all 2D arrays, rows x cols, NA -> null)
#   }
.pc_overview_heatmap <- function(base_path, method, combo, category,
                                   traits = "NULL", donors = "all") {
  if (method == "NULL" || method == "") method <- "mofa"
  if (category == "NULL" || category == "") category <- "core"

  if (method == "wgcna") {
    pcorr_path <- file.path(base_path, "wgcna_combined_results", combo,
                             "module_phenotype_correlations.csv")
  } else {
    pcorr_path <- file.path(base_path, "mofa_results", combo,
                             "factor_phenotype_correlations.csv")
  }
  if (!file.exists(pcorr_path)) {
    return(jsonlite::toJSON(list(
      error = paste0("phenotype correlations not found for ",
                     method, "/", combo)), auto_unbox = TRUE))
  }
  pcorr <- read.csv(pcorr_path, stringsAsFactors = FALSE, check.names = FALSE)

  # When a multi-omics donor subset is active, recompute every (score, phenotype)
  # stat live from the scores + phenotype CSVs instead of returning the
  # precomputed full-cohort values. The precomputed CSV is still used here as
  # the source of truth for each phenotype's type (numeric / binary_categorical /
  # multilevel_categorical) — only the statistic, p-value, p_adj and n are
  # recomputed on the subset.
  if (identical(donors, "subset") && file.exists("donors_multiomics.rds")) {
    return(.pc_overview_heatmap_subset(base_path, method, combo, category,
                                       traits, pcorr))
  }

  # Resolve which phenotype columns to render.
  # Priority: explicit `traits` (comma-separated column names) overrides
  # `category`. This lets the frontend's phenotype-selection dialog send an
  # arbitrary subset across multiple categories.
  if (traits != "NULL" && traits != "") {
    pheno_cols <- trimws(strsplit(traits, ",", fixed = TRUE)[[1]])
    pheno_cols <- pheno_cols[nzchar(pheno_cols)]
  } else {
    pheno_file <- switch(category,
      "all"             = "all_features.csv",
      "isolation"       = "features_isolation.csv",
      "gsis"            = "features_gsis.csv",
      "electrophysiology"= "features_electrophysiology.csv",
      "grs"             = "features_grs.csv",
      "cell"            = "features_cell.csv",
      "seahorse"        = "features_seahorse.csv",
      "perifusion"      = "features_perifusion.csv",
      "lip_extract"     = "features_lip_extract.csv",
      "features_core.csv")
    pheno_path <- file.path(base_path, "phenotype", pheno_file)
    if (!file.exists(pheno_path)) {
      return(jsonlite::toJSON(list(
        error = paste0("phenotype file not found: ", pheno_file)),
        auto_unbox = TRUE))
    }
    pheno_cols <- setdiff(
      colnames(read.csv(pheno_path, nrow = 1, check.names = FALSE)),
      "sample_id")
  }

  pcorr <- pcorr[pcorr$phenotype %in% pheno_cols, , drop = FALSE]
  if (nrow(pcorr) == 0) {
    return(jsonlite::toJSON(list(
      method = method, combo = combo, category = category,
      rows = list(), cols = list(),
      col_types = list(), col_methods = list(),
      statistic = list(), p_value = list(), p_adj = list(), n = list(),
      message = "no correlations available for this category"),
      auto_unbox = TRUE))
  }

  # One method per matrix (no mixing) — Pearson for everything except
  # multi-level categoricals, which need ANOVA.  We deliberately do NOT mix
  # pearson + bicor on the same heatmap because that confuses readers.
  # bicor is still available row-by-row in factor_phenotype_correlations.csv
  # for users who want it.
  pick_method <- function(rows_for_pheno) {
    type <- rows_for_pheno$type[1]
    if (type == "multilevel_categorical") return("anova_eta")
    "pearson"
  }

  # Sort score names naturally: Factor1 < Factor2 < ... or ME1 < ME2 < ...
  scores <- unique(pcorr$score)
  num_part <- suppressWarnings(as.numeric(gsub("[^0-9]", "", scores)))
  scores <- scores[order(num_part, scores)]

  # Preserve the phenotype CSV's original column order
  cols_present <- pheno_cols[pheno_cols %in% pcorr$phenotype]
  n_rows <- length(scores); n_cols <- length(cols_present)

  stat_mat <- matrix(NA_real_, n_rows, n_cols)
  p_mat    <- matrix(NA_real_, n_rows, n_cols)
  padj_mat <- matrix(NA_real_, n_rows, n_cols)
  n_mat    <- matrix(NA_integer_, n_rows, n_cols)
  col_methods <- character(n_cols)
  col_types   <- character(n_cols)

  for (j in seq_along(cols_present)) {
    pn <- cols_present[j]
    rows_pn <- pcorr[pcorr$phenotype == pn, , drop = FALSE]
    chosen  <- pick_method(rows_pn)
    col_methods[j] <- chosen
    col_types[j]   <- rows_pn$type[1]
    rows_pn <- rows_pn[rows_pn$method == chosen, , drop = FALSE]
    for (i in seq_along(scores)) {
      r <- rows_pn[rows_pn$score == scores[i], , drop = FALSE]
      if (nrow(r) == 1) {
        stat_mat[i, j] <- r$statistic
        p_mat[i, j]    <- r$p_value
        padj_mat[i, j] <- r$p_adj
        n_mat[i, j]    <- r$n
      }
    }
  }

  # Effective donor count (best-covered phenotype), then drop phenotype columns
  # with no computable value (all-NA) so the heatmap only shows real data.
  nd <- suppressWarnings(max(n_mat, na.rm = TRUE)); if (!is.finite(nd)) nd <- NA
  keep <- which(colSums(!is.na(stat_mat)) > 0)
  cols_present <- cols_present[keep]
  col_types    <- col_types[keep]
  col_methods  <- col_methods[keep]
  stat_mat     <- stat_mat[, keep, drop = FALSE]
  p_mat        <- p_mat[, keep, drop = FALSE]
  padj_mat     <- padj_mat[, keep, drop = FALSE]
  n_mat        <- n_mat[, keep, drop = FALSE]

  out <- list(
    method      = jsonlite::unbox(method),
    combo       = jsonlite::unbox(combo),
    category    = jsonlite::unbox(category),
    n_donors    = jsonlite::unbox(nd),
    rows        = scores,
    cols        = cols_present,
    col_types   = col_types,
    col_methods = col_methods,
    statistic   = stat_mat,
    p_value     = p_mat,
    p_adj       = padj_mat,
    n           = n_mat
  )
  # auto_unbox=FALSE keeps rows/cols/col_types/col_methods as arrays even when
  # length 1. Scalars (method/combo/category) are explicitly unboxed above.
  jsonlite::toJSON(out, na = "null", auto_unbox = FALSE,
                   matrix = "rowmajor", digits = 6)
}


# ---- Overview heatmap recomputed on the multi-omics donor subset ----
# Mirrors .pc_overview_heatmap's output JSON exactly, but every (score, phenotype)
# statistic/p_value/p_adj/n is recomputed from the score + phenotype CSVs after
# restricting to donors_multiomics.rds. Phenotype types come from `pcorr` (the
# precomputed correlations table), so test selection stays consistent with the
# full-cohort heatmap:
#   numeric, binary_categorical  -> pearson  (signed)
#   multilevel_categorical        -> anova_eta = sqrt(eta^2)   (unsigned 0..1)
.pc_overview_heatmap_subset <- function(base_path, method, combo, category,
                                        traits, pcorr) {
  # ---- load scores (eigengenes for WGCNA, factors for MOFA) ----
  if (method == "wgcna") {
    scores_path <- file.path(base_path, "wgcna_combined_results", combo, "eigengenes.csv")
  } else {
    scores_path <- file.path(base_path, "mofa_results", combo, "factors.csv")
  }
  if (!file.exists(scores_path)) {
    return(jsonlite::toJSON(list(
      error = paste0("scores file not found: ", scores_path)),
      auto_unbox = TRUE))
  }
  sc <- read.csv(scores_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (method == "wgcna") {
    donor_ids  <- sc[[1]]
    score_cols <- colnames(sc)[-1]
    scores_mat <- sc[, -1, drop = FALSE]
  } else {
    donor_ids  <- sc$donor
    score_cols <- setdiff(colnames(sc), "donor")
    scores_mat <- sc[, score_cols, drop = FALSE]
  }
  rownames(scores_mat) <- donor_ids

  # ---- resolve phenotype columns (same rules as .pc_overview_heatmap) ----
  if (traits != "NULL" && traits != "") {
    pheno_cols <- trimws(strsplit(traits, ",", fixed = TRUE)[[1]])
    pheno_cols <- pheno_cols[nzchar(pheno_cols)]
    pheno_file <- "all_features.csv"
  } else {
    pheno_file <- switch(category,
      "all"              = "all_features.csv",
      "isolation"        = "features_isolation.csv",
      "gsis"             = "features_gsis.csv",
      "electrophysiology"= "features_electrophysiology.csv",
      "grs"              = "features_grs.csv",
      "cell"             = "features_cell.csv",
      "seahorse"         = "features_seahorse.csv",
      "perifusion"       = "features_perifusion.csv",
      "lip_extract"      = "features_lip_extract.csv",
      "features_core.csv")
    pheno_cols <- NULL  # filled in after pheno is read
  }
  pheno_path <- file.path(base_path, "phenotype", pheno_file)
  if (!file.exists(pheno_path)) {
    return(jsonlite::toJSON(list(
      error = paste0("phenotype file not found: ", pheno_file)),
      auto_unbox = TRUE))
  }
  pheno <- read.csv(pheno_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (is.null(pheno_cols)) {
    pheno_cols <- setdiff(colnames(pheno), "sample_id")
  } else {
    pheno_cols <- intersect(pheno_cols, colnames(pheno))
  }

  # restrict to phenotypes that exist in the precomputed table (so type lookup works)
  pheno_cols <- pheno_cols[pheno_cols %in% pcorr$phenotype]
  if (length(pheno_cols) == 0) {
    return(jsonlite::toJSON(list(
      method = method, combo = combo, category = category,
      rows = list(), cols = list(),
      col_types = list(), col_methods = list(),
      statistic = list(), p_value = list(), p_adj = list(), n = list(),
      message = "no correlations available for this category"),
      auto_unbox = TRUE))
  }

  # ---- restrict donors to the subset, then to the score×pheno overlap ----
  common <- intersect(donor_ids, pheno$sample_id)
  common <- .apply_multiomics_subset(common, "subset")
  if (length(common) < 5) {
    return(jsonlite::toJSON(list(
      method = method, combo = combo, category = category,
      rows = list(), cols = list(),
      col_types = list(), col_methods = list(),
      statistic = list(), p_value = list(), p_adj = list(), n = list(),
      message = "too few donors in subset for correlations"),
      auto_unbox = TRUE))
  }
  scores_sub <- scores_mat[match(common, donor_ids), , drop = FALSE]
  pheno_sub  <- pheno[match(common, pheno$sample_id), , drop = FALSE]

  # ---- sort score names naturally (Factor1 < Factor2, ME1 < ME2) ----
  num_part <- suppressWarnings(as.numeric(gsub("[^0-9]", "", score_cols)))
  score_cols <- score_cols[order(num_part, score_cols)]

  n_rows <- length(score_cols)
  n_cols <- length(pheno_cols)
  stat_mat <- matrix(NA_real_,    n_rows, n_cols)
  p_mat    <- matrix(NA_real_,    n_rows, n_cols)
  n_mat    <- matrix(NA_integer_, n_rows, n_cols)
  col_methods <- character(n_cols)
  col_types   <- character(n_cols)

  # ---- per-phenotype: look up type from pcorr, then compute on subset ----
  for (j in seq_along(pheno_cols)) {
    pn      <- pheno_cols[j]
    rows_pn <- pcorr[pcorr$phenotype == pn, , drop = FALSE]
    ptype   <- if (nrow(rows_pn) > 0) rows_pn$type[1] else "numeric"
    chosen  <- if (ptype == "multilevel_categorical") "anova_eta" else "pearson"
    col_methods[j] <- chosen
    col_types[j]   <- ptype

    y_raw <- pheno_sub[[pn]]
    if (chosen == "pearson") {
      # numeric or binary (binary already encoded 0/1 in features_*.csv)
      y_num <- suppressWarnings(as.numeric(y_raw))
      for (i in seq_along(score_cols)) {
        x <- as.numeric(scores_sub[[score_cols[i]]])
        use <- !is.na(x) & !is.na(y_num)
        if (sum(use) < 5 || sd(x[use]) == 0 || sd(y_num[use]) == 0) next
        ct <- tryCatch(cor.test(x[use], y_num[use], method = "pearson"),
                       error = function(e) NULL)
        if (is.null(ct)) next
        stat_mat[i, j] <- ct$estimate
        p_mat[i, j]    <- ct$p.value
        n_mat[i, j]    <- sum(use)
      }
    } else {
      # anova_eta: unsigned 0..1
      g_raw <- as.character(y_raw)
      for (i in seq_along(score_cols)) {
        x <- as.numeric(scores_sub[[score_cols[i]]])
        use <- !is.na(x) & !is.na(g_raw) & g_raw != "" & g_raw != "NA"
        if (sum(use) < 5) next
        grp <- factor(g_raw[use])
        if (nlevels(grp) < 2) next
        fit <- tryCatch(stats::aov(x[use] ~ grp), error = function(e) NULL)
        if (is.null(fit)) next
        s   <- summary(fit)[[1]]
        ss  <- s[["Sum Sq"]]
        if (sum(ss, na.rm = TRUE) == 0) next
        eta <- ss[1] / sum(ss)
        stat_mat[i, j] <- sqrt(max(0, eta))
        p_mat[i, j]    <- s[["Pr(>F)"]][1]
        n_mat[i, j]    <- sum(use)
      }
    }
  }

  # BH-adjust across the whole matrix (same convention as the precomputed CSV)
  pvec <- as.vector(p_mat)
  padj_vec <- p.adjust(pvec, method = "BH")
  padj_mat <- matrix(padj_vec, nrow = n_rows, ncol = n_cols)

  # Drop phenotype columns that are all-NA in this subset (constant phenotype such
  # as Sex within a single-sex subset, <2 groups, etc.) so the heatmap only shows
  # phenotypes that actually vary in what the user subset to.
  keep <- which(colSums(!is.na(stat_mat)) > 0)
  if (length(keep) == 0) {
    return(jsonlite::toJSON(list(
      method = method, combo = combo, category = category,
      rows = list(), cols = list(), col_types = list(), col_methods = list(),
      statistic = list(), p_value = list(), p_adj = list(), n = list(),
      n_donors = length(common),
      message = "no phenotypes vary within this subset"),
      auto_unbox = TRUE))
  }
  pheno_cols  <- pheno_cols[keep]
  col_types   <- col_types[keep]
  col_methods <- col_methods[keep]
  stat_mat    <- stat_mat[, keep, drop = FALSE]
  p_mat       <- p_mat[, keep, drop = FALSE]
  padj_mat    <- padj_mat[, keep, drop = FALSE]
  n_mat       <- n_mat[, keep, drop = FALSE]

  out <- list(
    method      = jsonlite::unbox(method),
    combo       = jsonlite::unbox(combo),
    category    = jsonlite::unbox(category),
    n_donors    = jsonlite::unbox(length(common)),
    rows        = score_cols,
    cols        = pheno_cols,
    col_types   = col_types,
    col_methods = col_methods,
    statistic   = stat_mat,
    p_value     = p_mat,
    p_adj       = padj_mat,
    n           = n_mat
  )
  jsonlite::toJSON(out, na = "null", auto_unbox = FALSE,
                   matrix = "rowmajor", digits = 6)
}


################################################################################
# Legacy functions kept below for backward compatibility
################################################################################

# Precomputed trait correlation (replaces Java-based traitcorr)
# Returns JSON string with correlations between eigengenes/factors and a trait
################################################################################
PrecomputedTraitCorr <- function(method, combo, trait, category, primaryType,
                                  class1 = "NULL", class2 = "NULL") {

  base_path <- paste0(other.tables.path, "processed_comprehensive/")

  # Load eigengene/factor scores
  if (method == "wgcna") {
    scores_path <- file.path(base_path, "wgcna_combined_results", combo, "eigengenes.csv")
  } else {
    scores_path <- file.path(base_path, "mofa_results", combo, "factors.csv")
  }

  if (!file.exists(scores_path)) {
    return(jsonlite::toJSON(list(error = paste0("Scores file not found: ", scores_path)),
                            auto_unbox = TRUE))
  }

  scores <- read.csv(scores_path, stringsAsFactors = FALSE, check.names = FALSE)

  # Standardize donor ID column
  if (method == "wgcna") {
    # WGCNA: first column is unnamed row names
    donor_ids <- scores[[1]]
    score_cols <- colnames(scores)[-1]
    scores_mat <- scores[, -1, drop = FALSE]
  } else {
    # MOFA: last column is "donor"
    donor_ids <- scores$donor
    score_cols <- setdiff(colnames(scores), "donor")
    scores_mat <- scores[, score_cols, drop = FALSE]
  }

  # Load phenotype file
  pheno_file <- switch(category,
    "all" = "all_features.csv",
    "isolation" = "features_isolation.csv",
    "gsis" = "features_gsis.csv",
    "electrophysiology" = "features_electrophysiology.csv",
    "grs" = "features_grs.csv",
    "cell" = "features_cell.csv",
    "seahorse" = "features_seahorse.csv",
    "perifusion" = "features_perifusion.csv",
    "lip_extract" = "features_lip_extract.csv",
    "features_core.csv")
  pheno_path <- file.path(base_path, "phenotype", pheno_file)

  if (!file.exists(pheno_path)) {
    return(jsonlite::toJSON(list(error = paste0("Phenotype file not found: ", pheno_path)),
                            auto_unbox = TRUE))
  }

  pheno <- read.csv(pheno_path, stringsAsFactors = FALSE, check.names = FALSE)

  # Match donors
  common <- intersect(donor_ids, pheno$sample_id)
  if (length(common) < 5) {
    return(jsonlite::toJSON(list(correlations = list()), auto_unbox = TRUE))
  }

  idx_scores <- match(common, donor_ids)
  idx_pheno <- match(common, pheno$sample_id)

  correlations <- list()

  if (primaryType == "cont") {
    # Continuous: Pearson correlation
    trait_vals <- suppressWarnings(as.numeric(pheno[[trait]][idx_pheno]))
    valid <- !is.na(trait_vals)

    for (comp in score_cols) {
      comp_vals <- as.numeric(scores_mat[[comp]][idx_scores])
      use <- valid & !is.na(comp_vals)
      n <- sum(use)
      if (n < 5) next

      ct <- cor.test(comp_vals[use], trait_vals[use], method = "pearson")
      correlations[[length(correlations) + 1]] <- list(
        component = comp,
        trait = trait,
        r = round(ct$estimate, 4),
        p = ct$p.value,
        n = n
      )
    }

  } else if (primaryType == "disc" && class1 != "NULL" && class2 != "NULL") {
    # Categorical binary: point-biserial correlation
    cat_vals <- as.character(pheno[[trait]][idx_pheno])

    for (comp in score_cols) {
      comp_vals <- as.numeric(scores_mat[[comp]][idx_scores])
      use <- !is.na(comp_vals) & (cat_vals == class1 | cat_vals == class2)
      n <- sum(use)
      if (n < 5) next

      x <- comp_vals[use]
      y <- ifelse(cat_vals[use] == class1, 1.0, 0.0)
      n1 <- sum(y == 1)
      n2 <- sum(y == 0)

      ct <- cor.test(x, y, method = "pearson")
      correlations[[length(correlations) + 1]] <- list(
        component = comp,
        trait = trait,
        r = round(ct$estimate, 4),
        p = ct$p.value,
        n = n,
        class1 = class1,
        class2 = class2,
        n1 = n1,
        n2 = n2
      )
    }
  }

  result <- list(
    method = method,
    combo = combo,
    primaryType = primaryType,
    correlations = correlations
  )

  return(jsonlite::toJSON(result, auto_unbox = TRUE))
}


################################################################################
# Multi-group trait correlation for Sankey visualization
# Returns per-group eigengene/factor statistics with Kruskal-Wallis test
################################################################################
PrecomputedMultigroupCorr <- function(method, combo, trait, category) {

  base_path <- paste0(other.tables.path, "processed_comprehensive/")

  # Load eigengene/factor scores
  if (method == "wgcna") {
    scores_path <- file.path(base_path, "wgcna_combined_results", combo, "eigengenes.csv")
  } else {
    scores_path <- file.path(base_path, "mofa_results", combo, "factors.csv")
  }

  if (!file.exists(scores_path)) {
    return(jsonlite::toJSON(list(error = paste0("Scores file not found: ", scores_path)),
                            auto_unbox = TRUE))
  }

  scores <- read.csv(scores_path, stringsAsFactors = FALSE, check.names = FALSE)

  if (method == "wgcna") {
    donor_ids <- scores[[1]]
    score_cols <- colnames(scores)[-1]
    scores_mat <- scores[, -1, drop = FALSE]
  } else {
    donor_ids <- scores$donor
    score_cols <- setdiff(colnames(scores), "donor")
    scores_mat <- scores[, score_cols, drop = FALSE]
  }

  # Load phenotype
  pheno_file <- switch(category,
    "all" = "all_features.csv",
    "isolation" = "features_isolation.csv",
    "gsis" = "features_gsis.csv",
    "electrophysiology" = "features_electrophysiology.csv",
    "grs" = "features_grs.csv",
    "cell" = "features_cell.csv",
    "seahorse" = "features_seahorse.csv",
    "perifusion" = "features_perifusion.csv",
    "lip_extract" = "features_lip_extract.csv",
    "features_core.csv")
  pheno_path <- file.path(base_path, "phenotype", pheno_file)

  if (!file.exists(pheno_path)) {
    return(jsonlite::toJSON(list(error = paste0("Phenotype file not found: ", pheno_path)),
                            auto_unbox = TRUE))
  }

  pheno <- read.csv(pheno_path, stringsAsFactors = FALSE, check.names = FALSE)

  # Match donors
  common <- intersect(donor_ids, pheno$sample_id)
  if (length(common) < 5) {
    return(jsonlite::toJSON(list(correlations = list()), auto_unbox = TRUE))
  }

  idx_scores <- match(common, donor_ids)
  idx_pheno <- match(common, pheno$sample_id)

  cat_vals <- as.character(pheno[[trait]][idx_pheno])
  valid_cat <- !is.na(cat_vals) & cat_vals != "" & cat_vals != "NA"

  correlations <- list()

  for (comp in score_cols) {
    comp_vals <- as.numeric(scores_mat[[comp]][idx_scores])
    use <- valid_cat & !is.na(comp_vals)

    groups_factor <- factor(cat_vals[use])
    if (nlevels(groups_factor) < 2) next
    n_total <- sum(use)
    if (n_total < 5) next

    vals <- comp_vals[use]
    overall_mean <- mean(vals)
    overall_sd <- sd(vals)
    if (overall_sd == 0) next

    # Kruskal-Wallis test
    kw <- tryCatch(
      kruskal.test(vals ~ groups_factor),
      error = function(e) list(statistic = NA, p.value = 1)
    )

    # Per-group statistics
    group_stats <- lapply(levels(groups_factor), function(g) {
      g_vals <- vals[groups_factor == g]
      g_mean <- mean(g_vals)
      g_n <- length(g_vals)
      g_cohen_d <- (g_mean - overall_mean) / overall_sd
      list(
        group = g,
        mean = round(g_mean, 4),
        n = g_n,
        cohen_d = round(g_cohen_d, 4)
      )
    })

    correlations[[length(correlations) + 1]] <- list(
      component = comp,
      trait = trait,
      overall_mean = round(overall_mean, 4),
      overall_sd = round(overall_sd, 4),
      kw_p = kw$p.value,
      kw_h = round(as.numeric(kw$statistic), 2),
      groups = group_stats
    )
  }

  result <- list(
    method = method,
    combo = combo,
    mode = "multigroup",
    correlations = correlations
  )

  return(jsonlite::toJSON(result, auto_unbox = TRUE))
}
