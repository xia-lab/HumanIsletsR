# R Functions for Metabolic Atlas data layers
# Author: Yao Lu
# Computes per-pathway or per-reaction scores for the Metabolic Metro Map overlay

################################################################################

getAtlasLayer <- function(omicsType, layerType, donorId, condition, analysisVar,
                          ref, contrast, fixedEffects, fdr, donors, version="v2") {

  library(RSQLite)
  library(dplyr)

  fdr <- as.numeric(fdr)
  if (fixedEffects == "NA") fixedEffects <- NULL else fixedEffects <- strsplit(fixedEffects, ";")[[1]]
  if (donorId == "" || donorId == "NA") donorId <- NULL
  if (analysisVar == "NA") analysisVar <- NULL
  if (ref == "NA") ref <- NULL
  if (contrast == "NA") contrast <- NULL

  # Connect to omics database (for feature tables)
  mydb <- dbConnect(SQLite(), paste0(sqlite.path, "HI_omics_v2.sqlite"))
  on.exit(dbDisconnect(mydb), add = TRUE)

  # Load metadata from HI_tables.sqlite (separate database)
  meta_db <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
  pheno <- dbReadTable(meta_db, "proc_metadata")
  dbDisconnect(meta_db)
  rownames(pheno) <- pheno$record_id

  # Load KEGG pathway library for gene-pathway mapping
  lib.path <- paste0(other.tables.path, "libraries/")
  gem_lib <- NULL
  if (file.exists(paste0(lib.path, "hsaGem_reaction.qs"))) {
    gem_lib <- qs::qread(paste0(lib.path, "hsaGem_reaction.qs"))
  }

  # ── Load omics data based on type ──────────────────────────────────────
  if (omicsType == "proc_flux") {
    if (condition == "HG_LG_ratio") {
      table.nm <- "proc_flux_ratio"
    } else {
      table.nm <- paste0("proc_flux_", condition)
    }
    feature_table <- dbReadTable(mydb, table.nm)
    rownames(feature_table) <- feature_table$rxn
    feature_ids <- feature_table$rxn
    feature_info <- feature_table[, c(1, 3, 2, 4)]
    feature_table <- feature_table[, -c(1:4)]
    is_flux <- TRUE

  } else if (omicsType == "proc_prot") {
    feature_table <- dbReadTable(mydb, "proc_prot")
    rownames(feature_table) <- feature_table$gene_id
    feature_ids <- feature_table$gene_id
    feature_info <- feature_table[, c(1:3)]
    feature_table <- feature_table[, -c(1:3)]
    is_flux <- FALSE

  } else if (omicsType == "proc_rnaseq") {
    feature_table <- dbReadTable(mydb, "proc_rnaseq")
    rownames(feature_table) <- feature_table$accession
    feature_ids <- feature_table$accession
    feature_info <- feature_table[, c(1:5)]
    feature_table <- feature_table[, -c(1:5)]
    is_flux <- FALSE

  } else if (omicsType == "proc_scrna") {
    library(rhdf5)

    # Parse cell type and glucose from condition (e.g. "Beta_5")
    cond_parts <- strsplit(condition, "_")[[1]]
    cell <- cond_parts[1]       # "Alpha" or "Beta"
    glucose <- cond_parts[2]    # "1", "5", or "10"

    # Load HDF5 file (v2 uses hdf5_V2 directory)
    sc.h5.path <- ifelse(version=="v2", h5.v2.path, h5.path)
    table.path <- paste0(sc.h5.path, "sc_", cell, "_", glucose, ".h5")
    if (!file.exists(table.path)) {
      return(paste0("Error: HDF5 file not found: sc_", cell, "_", glucose, ".h5"))
    }
    feature_table <- h5read(table.path, "data/norm_expression") %>% as.data.frame()
    cells <- h5read(table.path, "meta/cells/cellid")
    genes <- h5read(table.path, "meta/genes") %>% as.data.frame()
    H5close()

    # Set up gene IDs as rownames
    if(version=="v2"){
      # v2: symbol is the primary gene identifier
      feature_info <- data.frame(gene_id = genes[,1], symbol = genes[,1], stringsAsFactors = FALSE)
      if(ncol(genes) >= 2) feature_info$name <- genes[,2]
    } else {
      feature_info <- genes[, c(1, 3, 2)]
      colnames(feature_info) <- c("gene_id", "symbol", "name")
    }
    colnames(feature_table) <- cells
    genes.keep <- !is.na(feature_info$gene_id)
    feature_table <- feature_table[genes.keep, ]
    feature_info <- feature_info[genes.keep, ]
    rownames(feature_table) <- feature_info$gene_id

    # Load ephys_cell metadata for cell-to-donor mapping
    meta_db2 <- dbConnect(SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
    cell_meta <- dbReadTable(meta_db2, "ephys_cell")
    dbDisconnect(meta_db2)

    # Filter metadata to cells present in this HDF5 file
    cell_meta <- cell_meta[cell_meta$cell_id %in% colnames(feature_table), ]
    cell_meta <- cell_meta[!is.na(cell_meta$record_id), ]

    if (nrow(cell_meta) < 10) {
      return("Error: Fewer than 10 cells with metadata for this cell type / glucose")
    }

    # Aggregate cell-level expression to donor-level (mean per donor)
    donor_ids <- unique(cell_meta$record_id)
    donor_means <- sapply(donor_ids, function(d) {
      cell_ids <- cell_meta$cell_id[cell_meta$record_id == d]
      cell_ids <- cell_ids[cell_ids %in% colnames(feature_table)]
      if (length(cell_ids) == 0) return(rep(NA, nrow(feature_table)))
      if (length(cell_ids) == 1) return(as.numeric(feature_table[, cell_ids]))
      rowMeans(feature_table[, cell_ids], na.rm = TRUE)
    })
    colnames(donor_means) <- donor_ids
    feature_table <- as.data.frame(donor_means)
    feature_ids <- rownames(feature_table)
    is_flux <- FALSE

  } else {
    return("Error: Unsupported omics type for atlas layer")
  }

  # ── DONOR-SPECIFIC MODE ────────────────────────────────────────────────
  if (layerType == "donor" && !is.null(donorId)) {

    if (!(donorId %in% colnames(feature_table))) {
      return(paste0("Error: Donor ", donorId, " not found in ", omicsType, " data"))
    }

    donor_vals <- as.numeric(feature_table[, donorId])
    names(donor_vals) <- rownames(feature_table)

    # Compute per-feature z-scores (donor vs population)
    means <- rowMeans(feature_table, na.rm = TRUE)
    sds <- apply(feature_table, 1, sd, na.rm = TRUE)
    sds[sds == 0] <- 1
    zscores <- (donor_vals - means) / sds
    names(zscores) <- rownames(feature_table)
    zscores <- zscores[!is.na(zscores)]

    # For flux: save per-reaction z-scores for segment-click detail
    if (is_flux) {
      rxn_detail <- data.frame(id = names(zscores), value = round(zscores, 4),
                               stringsAsFactors = FALSE)
      write.csv(rxn_detail, "atlas_layer_reactions.csv", row.names = FALSE)
    }

    # Pathway-level test using fgsea (n=1, cannot use limma)
    if (!is.null(gem_lib)) {
      result <- pathway_fgsea_test(zscores, gem_lib, names(zscores), is_flux)
    } else {
      result <- data.frame(id = character(), value = numeric(),
                           pval = numeric(), padj = numeric(),
                           stringsAsFactors = FALSE)
    }
    score_definition <- "zscore_vs_population"
    comparison_label <- donorId

  # ── CLUSTER-SPECIFIC MODE ──────────────────────────────────────────────
  } else if (layerType == "cluster" && !is.null(donorId)) {

    library(limma)

    # donorId holds the cluster name (e.g. "C0")
    cluster_id <- donorId
    cluster_donors <- pheno$record_id[pheno$final_cluster == cluster_id]
    cluster_donors <- intersect(cluster_donors, colnames(feature_table))

    if (length(cluster_donors) < 3) {
      return(paste0("Error: Fewer than 3 donors in cluster ", cluster_id))
    }

    # Build binary grouping: cluster vs rest
    all_donors <- colnames(feature_table)
    group <- ifelse(all_donors %in% cluster_donors, "cluster", "rest")
    group <- factor(group, levels = c("rest", "cluster"))

    # Run limma: cluster vs rest
    design <- model.matrix(~ 0 + group)
    colnames(design) <- c("rest", "cluster")
    contrast.matrix <- makeContrasts(cluster - rest, levels = design)

    feat_mat <- as.matrix(feature_table)

    # Remove features with zero/near-zero variance (cause NA coefficients in limma)
    row_vars <- apply(feat_mat, 1, var, na.rm = TRUE)
    valid_rows <- !is.na(row_vars) & row_vars > 1e-10
    feat_mat <- feat_mat[valid_rows, , drop = FALSE]

    fit <- lmFit(feat_mat, design)
    fit_c <- contrasts.fit(fit, contrast.matrix)
    fit_c <- eBayes(fit_c)

    # For flux: save per-reaction limma results for segment-click detail
    if (is_flux) {
      res <- topTable(fit_c, number = Inf, sort.by = "none")
      rxn_detail <- data.frame(
        id = rownames(feat_mat),
        value = round(res$logFC, 4),
        pval = res$P.Value,
        padj = res$adj.P.Val,
        stringsAsFactors = FALSE
      )
      rxn_detail <- rxn_detail[!is.na(rxn_detail$value), ]
      write.csv(rxn_detail, "atlas_layer_reactions.csv", row.names = FALSE)
    }

    # Pathway-level test using camera
    if (!is.null(gem_lib)) {
      result <- pathway_camera_test(fit, contrast.matrix,
                                    feat_mat, gem_lib, is_flux)
    } else {
      result <- data.frame(id = character(), value = numeric(),
                           pval = numeric(), padj = numeric(),
                           stringsAsFactors = FALSE)
    }
    score_definition <- "log2fc_cluster_vs_rest"
    comparison_label <- cluster_id

  # ── COMPARISON MODE ────────────────────────────────────────────────────
  } else if (layerType == "comparison") {

    if (is.null(analysisVar) || is.null(ref) || is.null(contrast)) {
      if (is.null(analysisVar)) {
        return("Error: Missing comparison parameter (analysisVar)")
      }
    }

    library(limma)

    # Align donors between phenotype and feature tables
    common_donors <- intersect(colnames(feature_table), pheno$record_id)
    if (length(common_donors) < 10) {
      return("Error: Fewer than 10 donors with matching data")
    }

    feat_mat <- feature_table[, common_donors]
    pheno_sub <- pheno[match(common_donors, pheno$record_id), ]

    if (!(analysisVar %in% colnames(pheno_sub))) {
      return(paste0("Error: Variable '", analysisVar, "' not found in phenotype data"))
    }

    is_continuous <- is.null(ref) || is.null(contrast)

    if (!is_continuous) {
      groups <- pheno_sub[[analysisVar]]
      group1_idx <- which(groups == contrast)
      group2_idx <- which(groups == ref)

      if (length(group1_idx) < 5 || length(group2_idx) < 5) {
        return("Error: Fewer than 5 donors in one or both comparison groups")
      }

      # limma differential analysis
      group_labels <- rep(NA, ncol(feat_mat))
      group_labels[group1_idx] <- "grpB"
      group_labels[group2_idx] <- "grpA"
      keep <- !is.na(group_labels)
      feat_sub <- as.matrix(feat_mat[, keep])
      # Remove features with zero/near-zero variance
      row_vars <- apply(feat_sub, 1, var, na.rm = TRUE)
      valid_rows <- !is.na(row_vars) & row_vars > 1e-10
      feat_sub <- feat_sub[valid_rows, , drop = FALSE]

      group_factor <- factor(group_labels[keep], levels = c("grpA", "grpB"))

      design <- model.matrix(~ 0 + group_factor)
      colnames(design) <- c("grpA", "grpB")
      contrast.matrix <- makeContrasts(grpB - grpA, levels = design)

      fit <- lmFit(feat_sub, design)
      fit <- eBayes(fit)

      # Save per-feature limma results for flux segment-click detail
      if (is_flux) {
        fit_c <- contrasts.fit(fit, contrast.matrix)
        fit_c <- eBayes(fit_c)
        res <- topTable(fit_c, number = Inf, sort.by = "none")
        rxn_detail <- data.frame(
          id = rownames(feat_sub),
          value = round(res$logFC, 4),
          pval = res$P.Value,
          padj = res$adj.P.Val,
          stringsAsFactors = FALSE
        )
        rxn_detail <- rxn_detail[!is.na(rxn_detail$value), ]
        write.csv(rxn_detail, "atlas_layer_reactions.csv", row.names = FALSE)
      }

      # Pathway-level test using camera
      if (!is.null(gem_lib)) {
        result <- pathway_camera_test(fit, contrast.matrix, feat_sub,
                                      gem_lib, is_flux)
      } else {
        result <- data.frame(id = character(), value = numeric(),
                             pval = numeric(), padj = numeric(),
                             stringsAsFactors = FALSE)
      }
      score_definition <- "log2fc_grpB_vs_grpA"
      comparison_label <- paste0(contrast, "_vs_", ref)

    } else {
      # Continuous metadata: regression slope per unit change
      predictor <- pheno_sub[[analysisVar]]
      if (!is.numeric(predictor)) predictor <- suppressWarnings(as.numeric(predictor))
      if (all(is.na(predictor))) {
        return(paste0("Error: Variable '", analysisVar, "' is not numeric for continuous comparison"))
      }

      design_df <- data.frame(predictor = predictor)
      if (!is.null(fixedEffects) && length(fixedEffects) > 0) {
        keep_cov <- fixedEffects[fixedEffects %in% colnames(pheno_sub)]
        if (length(keep_cov) > 0) {
          design_df <- cbind(design_df, pheno_sub[, keep_cov, drop = FALSE])
        }
      }

      complete <- complete.cases(design_df)
      if (sum(complete) < 10) {
        return("Error: Fewer than 10 donors with complete metadata for continuous comparison")
      }

      feat_sub <- as.matrix(feat_mat[, complete])
      # Remove features with zero/near-zero variance
      row_vars <- apply(feat_sub, 1, var, na.rm = TRUE)
      valid_rows <- !is.na(row_vars) & row_vars > 1e-10
      feat_sub <- feat_sub[valid_rows, , drop = FALSE]

      design <- model.matrix(~ ., data = design_df[complete, , drop = FALSE])

      fit <- lmFit(feat_sub, design)
      fit <- eBayes(fit)
      coef_name <- colnames(design)[2]  # predictor

      # Save per-feature limma results for flux segment-click detail
      if (is_flux) {
        res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
        rxn_detail <- data.frame(
          id = rownames(feat_sub),
          value = round(res$logFC, 4),
          pval = res$P.Value,
          padj = res$adj.P.Val,
          stringsAsFactors = FALSE
        )
        rxn_detail <- rxn_detail[!is.na(rxn_detail$value), ]
        write.csv(rxn_detail, "atlas_layer_reactions.csv", row.names = FALSE)
      }

      # Pathway-level test using camera (coef_name is column index 2)
      if (!is.null(gem_lib)) {
        result <- pathway_camera_test(fit, 2, feat_sub, gem_lib, is_flux,
                                      is_coef = TRUE)
      } else {
        result <- data.frame(id = character(), value = numeric(),
                             pval = numeric(), padj = numeric(),
                             stringsAsFactors = FALSE)
      }
      score_definition <- paste0("coef_", analysisVar)
      comparison_label <- analysisVar
    }

  } else {
    return("Error: Invalid layerType. Use 'donor', 'cluster', or 'comparison'")
  }

  # Add metadata columns
  if (nrow(result) > 0) {
    if (!("pval" %in% names(result))) result$pval <- NA
    if (!("padj" %in% names(result))) result$padj <- NA
    result$score_definition <- score_definition
    result$comparison <- comparison_label
    result$direction <- ifelse(result$value > 0, "up", "down")
  }

  # Write output CSV
  write.csv(result, "atlas_layer.csv", row.names = FALSE)

  return(paste0("RES-OK;", nrow(result)))
}


################################################################################
# Helper: Build GEM subsystem pathway list for flux (reaction-level) data
################################################################################

build_gem_rxn_pathways <- function(gem_lib, feature_ids) {
  rxn_tab <- gem_lib$reacTab
  pathway_names <- unique(rxn_tab[, c("sysID", "pathwayVis")])

  # Split reaction IDs by subsystem
  pathway_list <- split(rxn_tab$humanGEM, rxn_tab$sysID)

  # Filter to reactions present in feature data
  pathway_list <- lapply(pathway_list, function(rxns) intersect(rxns, feature_ids))
  pathway_list <- pathway_list[sapply(pathway_list, length) >= 3]

  list(sets = pathway_list, names = pathway_names)
}

################################################################################
# Helper: Build GEM subsystem pathway list for gene-level data (via GPR rules)
################################################################################

build_gem_gene_pathways <- function(gem_lib, feature_ids) {
  rxn_tab <- gem_lib$reacTab
  rules <- gem_lib$rules
  pathway_names <- unique(rxn_tab[, c("sysID", "pathwayVis")])

  # Map each reaction to its genes from GPR rules
  # rules is a list; each element is a data.frame with 'id' column (gene IDs)
  rxn_genes <- lapply(rules, function(r) {
    ids <- r$id
    # Handle compound IDs (ENSG001 & ENSG002) by splitting
    all_genes <- unlist(lapply(ids, function(x) {
      if (grepl("&", x)) trimws(strsplit(x, "&")[[1]]) else x
    }))
    intersect(unique(all_genes), feature_ids)
  })

  # Strip compound names (MAR00001.ENSG...) to get reaction IDs
  rxn_ids <- sub("\\..*$", "", names(rxn_genes))

  # Build subsystem → gene sets
  # For each subsystem, collect all genes from its reactions
  pathway_list <- list()
  for (i in seq_along(rxn_genes)) {
    rxn_id <- rxn_ids[i]
    # Find which subsystem this reaction belongs to
    sys_idx <- match(rxn_id, rxn_tab$humanGEM)
    if (is.na(sys_idx)) next
    sys_id <- rxn_tab$sysID[sys_idx]
    if (is.null(pathway_list[[sys_id]])) pathway_list[[sys_id]] <- character()
    pathway_list[[sys_id]] <- c(pathway_list[[sys_id]], rxn_genes[[i]])
  }
  # Deduplicate genes per pathway
  pathway_list <- lapply(pathway_list, unique)
  pathway_list <- pathway_list[sapply(pathway_list, length) >= 3]

  list(sets = pathway_list, names = pathway_names)
}

################################################################################
# Helper: fgsea pathway test (for donor mode, n=1)
################################################################################

pathway_fgsea_test <- function(feature_ranks, gem_lib, feature_ids, is_flux) {
  library(fgsea)

  if (is_flux) {
    pw_info <- build_gem_rxn_pathways(gem_lib, feature_ids)
  } else {
    pw_info <- build_gem_gene_pathways(gem_lib, feature_ids)
  }
  pathway_list <- pw_info$sets
  pathway_names <- pw_info$names

  if (length(pathway_list) == 0) {
    return(data.frame(id = character(), value = numeric(),
                      pval = numeric(), padj = numeric(),
                      stringsAsFactors = FALSE))
  }

  set.seed(42)
  res <- fgsea(pathways = pathway_list, stats = feature_ranks,
               minSize = 3, maxSize = 500)

  if (nrow(res) == 0) {
    return(data.frame(id = character(), value = numeric(),
                      pval = numeric(), padj = numeric(),
                      stringsAsFactors = FALSE))
  }

  # Compute mean score per pathway for map coloring/sizing
  pw_means <- sapply(res$pathway, function(pw) {
    members <- intersect(pathway_list[[pw]], names(feature_ranks))
    if (length(members) == 0) return(0)
    mean(feature_ranks[members], na.rm = TRUE)
  })

  result <- data.frame(
    id = pathway_names$pathwayVis[match(res$pathway, pathway_names$sysID)],
    value = round(pw_means, 4),
    pval = res$pval,
    padj = res$padj,
    n_total = res$size,
    stringsAsFactors = FALSE
  )
  result <- result[!is.na(result$id), ]
  return(result)
}

################################################################################
# Helper: camera pathway test (for cluster/comparison modes with limma fit)
################################################################################

pathway_camera_test <- function(fit, contrast_or_coef, feature_table, gem_lib,
                                is_flux, is_coef = FALSE) {
  library(limma)

  all_features <- rownames(feature_table)

  if (is_flux) {
    pw_info <- build_gem_rxn_pathways(gem_lib, all_features)
  } else {
    pw_info <- build_gem_gene_pathways(gem_lib, all_features)
  }
  pathway_list <- pw_info$sets
  pathway_names <- pw_info$names

  if (length(pathway_list) == 0) {
    return(data.frame(id = character(), value = numeric(),
                      pval = numeric(), padj = numeric(),
                      stringsAsFactors = FALSE))
  }

  # Convert pathway gene/reaction lists to row indices
  pathway_indices <- lapply(pathway_list, function(ids) {
    idx <- match(ids, all_features)
    idx[!is.na(idx)]
  })
  pathway_indices <- pathway_indices[sapply(pathway_indices, length) >= 3]

  # Remove rows with any NA from feature_table before camera
  row_ok <- apply(feature_table, 1, function(x) all(is.finite(x)))
  if (any(!row_ok)) {
    feature_table <- feature_table[row_ok, , drop = FALSE]
    all_features <- rownames(feature_table)
    pathway_indices <- lapply(pathway_list, function(ids) {
      idx <- match(ids, all_features)
      idx[!is.na(idx)]
    })
    pathway_indices <- pathway_indices[sapply(pathway_indices, length) >= 3]
    # Re-fit limma on cleaned data
    fit <- lmFit(feature_table, fit$design)
    fit <- eBayes(fit)
  }

  if (length(pathway_indices) == 0) {
    return(data.frame(id = character(), value = numeric(),
                      pval = numeric(), padj = numeric(),
                      stringsAsFactors = FALSE))
  }

  # Run camera
  cam_res <- tryCatch({
    if (is_coef) {
      camera(feature_table, index = pathway_indices,
             design = fit$design, coef = contrast_or_coef)
    } else {
      camera(feature_table, index = pathway_indices,
             design = fit$design, contrast = contrast_or_coef)
    }
  }, error = function(e) {
    print(paste("camera error:", e$message))
    return(NULL)
  })

  if (is.null(cam_res) || nrow(cam_res) == 0) {
    return(data.frame(id = character(), value = numeric(),
                      pval = numeric(), padj = numeric(),
                      stringsAsFactors = FALSE))
  }

  # Get per-feature statistics for computing pathway mean scores
  if (is_coef) {
    feature_stats <- fit$coefficients[, contrast_or_coef]
  } else {
    fit2 <- contrasts.fit(fit, contrast_or_coef)
    feature_stats <- fit2$coefficients[, 1]
  }
  names(feature_stats) <- all_features

  # Compute mean statistic per pathway for map coloring
  matched_sysIDs <- intersect(rownames(cam_res), names(pathway_list))
  pw_means <- sapply(matched_sysIDs, function(pw) {
    members <- intersect(pathway_list[[pw]], names(feature_stats))
    if (length(members) == 0) return(0)
    mean(feature_stats[members], na.rm = TRUE)
  })

  # Apply camera direction sign to the mean score
  for (pw in matched_sysIDs) {
    if (cam_res[pw, "Direction"] == "Down") {
      pw_means[pw] <- -abs(pw_means[pw])
    } else {
      pw_means[pw] <- abs(pw_means[pw])
    }
  }

  result <- data.frame(
    id = pathway_names$pathwayVis[match(matched_sysIDs, pathway_names$sysID)],
    value = round(pw_means, 4),
    pval = cam_res[matched_sysIDs, "PValue"],
    padj = cam_res[matched_sysIDs, "FDR"],
    n_total = cam_res[matched_sysIDs, "NGenes"],
    stringsAsFactors = FALSE
  )
  result <- result[!is.na(result$id), ]
  return(result)
}
