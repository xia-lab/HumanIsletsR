# R Functions for metadata analysis on HumanIslets web tool V2
# Author: Yao Lu

# Pre-load all packages at source time (avoids 4+ sec delay per function call)
suppressPackageStartupMessages({
  #library(pcaMethods)
  library(ggplot2)
  library(plotly)
  library(dplyr)
  library(jsonlite)
  library(ggrepel)
})

################################################################################

getHeatmapData<- function(reset="no"){
   print(c(reset,"reset==="))
  # If reset is "yes", ignore any existing donor/phenotype filters
  if(reset == "yes" && file.exists("donors.rds")){
    file.remove("donors.rds")
  }
  if(reset == "yes" && file.exists("rmphemo.rds")){
    file.remove("rmphemo.rds")
  }

  dat<-read.csv(paste0(other.tables.path,"display_data/metadata_sum_raw.csv"))

  if(file.exists("donors.rds")){
    donors <- readRDS("donors.rds")
   dat <- dat[dat$record_id %in%donors, ]

  }

  if(file.exists("rmphemo.rds")){
  rmphemo <- readRDS("rmphemo.rds")
     dat <- dat[, !names(dat) %in% rmphemo]
  }
 if(!is.data.frame(dat)){
      res = "NO"

 }else{

 datnm <- paste0("heatmapData",sample(1:10000, 1),".csv")
 write.csv(dat,datnm,row.names = F)
  res <- paste0("RES-OK;",datnm);

  }

  
  print(res)
  return(res)
}

################################################################################

getChordData<- function(corrcutoff=0.5,pvalcutoff=0.05,intragroup,reset="no"){
   print(c(reset,"reset==="))
  corrcutoff <- as.numeric(corrcutoff);
  pvalcutoff <- as.numeric(pvalcutoff);

  # If reset is "yes", ignore any existing donor/phenotype filters
  if(reset == "yes" && file.exists("donors.rds")){
    file.remove("donors.rds")
  }
  if(reset == "yes" && file.exists("rmphemo.rds")){
    file.remove("rmphemo.rds")
  }

  # Always load meta_group early so it's available for all code paths
  meta_group = read.csv(paste0(other.tables.path,"display_interface/meta_groups.csv"))

  if(file.exists("donors.rds")){
    donors <- readRDS("donors.rds")
library(Hmisc)
     metadat<-read.csv(paste0(other.tables.path,"display_data/metadata_sum_norm.csv"))
 df_clean <- metadat[metadat$record_id %in%donors,-1]

df_clean[] <- lapply(df_clean, function(x) {
  if (is.character(x) || is.factor(x)) as.numeric(factor(x)) else x
})
df_clean <- df_clean[, colSums(!is.na(df_clean)) >= 6]

rc <- rcorr(as.matrix(df_clean), type = "spearman")

cor_df = reshape2::melt(rc$r)
pval_df = reshape2::melt(rc$P)
npairs_df = reshape2::melt(rc$n)
cor_df$pval = pval_df$value
cor_df$npairs = npairs_df$value

# Convert factor columns to character to prevent c() integer coercion
cor_df$Var1 <- as.character(cor_df$Var1)
cor_df$Var2 <- as.character(cor_df$Var2)

library(dplyr)
dat <- cor_df %>%
  filter(Var1 != Var2) %>%  # remove self-correlations
  rowwise() %>%
  mutate(pair_id = paste(sort(c(Var1, Var2)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair_id, .keep_all = TRUE) %>%  # keep only one of A-B or B-A
  select(Var1, Var2, value, pval,npairs)
   dat= dat[dat$Var1!=dat$Var2 & !is.na(dat$pval) & dat$npairs>5,]

  }else{
   dat<-read.csv(paste0(other.tables.path,"display_data/metadata_corr.csv"))
  }

  # Ensure Var1/Var2 are character (not factor) for reliable matching
  dat$Var1 <- as.character(dat$Var1)
  dat$Var2 <- as.character(dat$Var2)

  dat<-dat[abs(dat$value)>corrcutoff & dat$pval<pvalcutoff,]

   if(nrow(dat)==0){
   res<- "NO"

  }else{
    if(intragroup=="no"){
 dat$type1 <-meta_group$group[match(dat$Var1,meta_group$column)]
dat$type2 <- meta_group$group[match(dat$Var2,meta_group$column)]
dat <- dat[!is.na(dat$type1) & !is.na(dat$type2) & dat$type1!=dat$type2,]

  }
  }

 if(file.exists("rmphemo.rds")){
      rmphemo <- readRDS("rmphemo.rds")
     dat <- dat[(!dat$Var1 %in%rmphemo) & (!dat$Var2 %in%rmphemo ), ]

  }
  if(nrow(dat)==0){
   res<- "NO"

  }else{
  names(dat)[1:2] <- c("source","target")


  node <- unique(c(dat$source,dat$target))
  node <- node[!is.na(node)]
  node_groups = meta_group[meta_group$column %in% node, ]

  # Fallback: if some nodes not found in meta_group, add them with defaults
  missing <- setdiff(node, node_groups$column)
  if(length(missing) > 0) {
    print(paste("WARNING: nodes not found in meta_groups.csv:", paste(missing, collapse=", ")))
    fallback <- data.frame(
      group = sapply(missing, function(n) {
        # Try to get group from type1/type2 columns if available
        g <- NA
        if("type1" %in% names(dat)) g <- dat$type1[dat$source == n][1]
        if(is.na(g) && "type2" %in% names(dat)) g <- dat$type2[dat$target == n][1]
        if(is.na(g)) g <- "Other"
        return(as.character(g))
      }),
      column = missing,
      display = missing,
      type = rep("cont", length(missing)),
      groupID = rep("other", length(missing)),
      stringsAsFactors = FALSE
    )
    node_groups <- rbind(node_groups, fallback)
  }

  datnm <- paste0("chordData",sample(1:10000, 1),".json")


 jsonlite::write_json(list(nodes = node_groups, links = dat), datnm, pretty = TRUE)

  res <- paste0("RES-OK;",datnm);

 }



  return(res)
}

################################################################################

getPairData<- function(meta1,meta2){
library(ggplot2)
library(plotly)
   dat<-read.csv(paste0(other.tables.path,"display_data/metadata_sum_norm.csv"))
 
dat <- unique(na.omit(dat[,c("record_id",meta1,meta2)]));
 
  meta_group = read.csv(paste0(other.tables.path,"display_interface/meta_groups.csv"))
 
    type1 = meta_group$type[meta_group$column == meta1]
  type2 = meta_group$type[meta_group$column == meta2]

datnm <- paste0("pair",sample(1:10000, 1),".json")

  if(type1 == "cont" & type2 == "cont" ){
  
labels=meta_group$display[match(names(dat)[2:3],meta_group$column)]

p <- ggplot(dat, aes(x = .data[[meta1]], y = .data[[meta2]])) +
  geom_point(
    aes(text = paste0(
      "Donor: ", .data[["record_id"]],
      "<br>", labels[1], ": ", sprintf("%.2f", .data[[meta1]]),
      "<br>", labels[2], ": ", sprintf("%.2f", .data[[meta2]])
    )),
    size = 2,  color = "black" 
  ) +
  geom_smooth(method = "lm", se = TRUE, colour = NA, fill = "blue", alpha = 0.2) +
  geom_smooth(method = "lm", se = FALSE, colour = "blue", linewidth = 1)  +
  labs(title = "", x = labels[1], y = labels[2]) +
  theme_minimal()
pl <- ggplotly(p, tooltip = "text")

json <- plotly_json(pl, jsonedit = FALSE, pretty = TRUE)

cat(json, file = datnm)

  
   type = "line"
   }else if(type1 == "disc" | type2 == "cont" ){
cont_var =  names(dat)[3]
disc_var = names(dat)[2]
   
labels=meta_group$display[match(names(dat)[2:3],meta_group$column)]
 
p <- ggplot(dat, aes(x = .data[[disc_var]], y = .data[[cont_var]])) +
  geom_boxplot(aes(fill = .data[[disc_var]]), outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = .data[[cont_var]], text  = paste0(
    "Donor: ", record_id,
    "<br>", labels[1], ": ", .data[[disc_var]],
    "<br>", labels[2], ": ", round(.data[[cont_var]], 2)
  )),
              width = 0.2, size = 2, shape = 16) +   # shape 16 = solid, no border
  scale_color_gradient(low = "lightblue", high = "darkblue", name = labels[1]) +
  ggsci::scale_fill_npg(name = labels[2]) +
  labs(x = "", y =  labels[2], fill = disc_var) +
  theme_minimal()

pl=ggplotly(p,tooltip = "text")

pl[["x"]][["layout"]][["legend"]][["title"]][["text"]] =  labels[2]


json <- plotly_json(pl, jsonedit = FALSE, pretty = TRUE)

cat(json, file = datnm)

    type = "violin"
   }else if(type1 == "cont" & type2 == "disc" ){
      dat <- dat[,c(1,3,2)]
cont_var =  names(dat)[3]
disc_var = names(dat)[2]
labels=meta_group$display[match(names(dat)[2:3],meta_group$column)]
 
p <- ggplot(dat, aes(x = .data[[disc_var]], y = .data[[cont_var]])) +
  geom_boxplot(aes(fill = .data[[disc_var]]), outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = .data[[cont_var]], text  = paste0(
    "Donor: ", record_id,
    "<br>", labels[1], ": ", .data[[disc_var]],
    "<br>", labels[2], ": ", round(.data[[cont_var]], 2)
  )),
              width = 0.2, size = 2, shape = 16) +   # shape 16 = solid, no border
  scale_color_gradient(low = "lightblue", high = "darkblue", name = labels[1]) +
  ggsci::scale_fill_npg(name = labels[2]) +
  labs(x = "", y =  labels[2], fill = disc_var) +
  theme_minimal()

pl=ggplotly(p,tooltip = "text")

pl[["x"]][["layout"]][["legend"]][["title"]][["text"]] =  labels[1]


json <- plotly_json(pl, jsonedit = FALSE, pretty = TRUE)

cat(json, file = datnm)
    type = "violin"
   }else{
library(dplyr)
labels=meta_group$display[match(names(dat)[2:3],meta_group$column)]

p <- dat %>%
  arrange(meta1, meta2) %>%
  group_by(meta1) %>%
  mutate(ypos = row_number()) %>%
  ggplot(aes(
    x = meta1,
    y = ypos,
    fill = meta2,
    text = paste0(
    "Donor: ", record_id,
    "<br>", labels[1], ": ", .data[[meta1]],
    "<br>", labels[2], ": ", .data[[meta2]]
    )
  )) +
  geom_tile(width = 0.6, height = 1, color = "white") +
  ggsci::scale_fill_npg() +
  labs(x = "", y = "Donors", fill = labels[2]) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

pl <- ggplotly(p, tooltip = "text")

cat(json, file = datnm)


    type="stack"
  }


 # datnm <- paste0("pair",sample(1:10000, 1),".csv")
  
  #write.csv(dat,datnm,row.names = F)



  res <- paste0("RES-OK;",datnm,";",type, paste0(labels[1], " vs ",labels[2]) );
  
  return(res)
}


################################################################################

#' Phenotype PCA Analysis using Probabilistic PCA
#'
#' Performs Probabilistic PCA (PPCA) on phenotype data, handles missing values,
#' and returns interactive plot data for the web interface.
#'
#' @param colorBy Variable to color points by (e.g., "diagnosis", "purity", "donorsex")
#' @param pcX Principal component for X-axis (default: 1)
#' @param pcY Principal component for Y-axis (default: 2)
#' @param nComponents Number of principal components to compute (default: 10)
#' @param minSamples Minimum samples required for a feature (default: 50)
#' @param highlightOutliers Number of outlier donors to highlight (default: 10)
#' @param reset Whether to reset donor/phenotype filters
#' @return JSON string with plot data and PCA results

getPhenotypePCA <- function(colorBy = "diagnosis",
                            pcX = 1,
                            pcY = 2,
                            nComponents = 10,
                            minSamples = 50,
                            highlightOutliers = 10,
                            reset = "no",
                            shapeBy = "") {

  # Use non-interactive graphics backend for headless server
  if(capabilities("cairo")) {
    options(bitmapType = "cairo")
  }
  # Ensure a graphics device is available (needed by ggplotly on headless servers)
  if(length(dev.list()) == 0) {
    pdf(NULL)
  }
 
  pcX <- as.numeric(pcX)
  pcY <- as.numeric(pcY)
  nComponents <- as.numeric(nComponents)
  minSamples <- as.numeric(minSamples)
  highlightOutliers <- as.numeric(highlightOutliers)

  # Handle filter reset
  if(reset == "yes" && file.exists("donors.rds")) {
    file.remove("donors.rds")
  }
  if(reset == "yes" && file.exists("rmphemo.rds")) {
    file.remove("rmphemo.rds")
  }

  # Load phenotype data (normalized for PCA)
  dat <- read.csv(paste0(other.tables.path, "display_data/metadata_sum_norm.csv"))
  meta_group <- read.csv(paste0(other.tables.path, "display_interface/meta_groups.csv"))

  # Apply donor filter if exists
  if(file.exists("donors.rds")) {
    donors <- readRDS("donors.rds")
    dat <- dat[dat$record_id %in% donors, ]
  }

  # Apply phenotype filter if exists
  if(file.exists("rmphemo.rds")) {
    rmphemo <- readRDS("rmphemo.rds")
    dat <- dat[, !names(dat) %in% rmphemo]
  }

  if(nrow(dat) < 10) {
    return("NO;Insufficient samples for PCA analysis")
  }

  # Identify numeric columns only (continuous phenotypes)
  numeric_cols <- meta_group$column[meta_group$type == "cont" & meta_group$column %in% names(dat)] 

  # Keep colorBy column for later
  colorVar <- NULL
  if(colorBy %in% names(dat)) {
    colorVar <- dat[[colorBy]]
  }

  # Keep shapeBy column for later (optional)
  shapeVar <- NULL
  if(nzchar(shapeBy) && shapeBy %in% names(dat)) {
    shapeVar <- dat[[shapeBy]]
  }

  # Prepare data matrix for PCA
  record_ids <- dat$record_id
  pca_data <- dat[, numeric_cols, drop = FALSE]

  # Remove columns with too few samples
  valid_cols <- colSums(!is.na(pca_data)) >= minSamples
  pca_data <- pca_data[, valid_cols, drop = FALSE]
 
  if(ncol(pca_data) < 5) { 
    return("NO;Insufficient features for PCA analysis after filtering")
  }

  # Scale the data (center and scale)
  pca_matrix <- as.matrix(pca_data)

  # Calculate column means and sds ignoring NA
  col_means <- colMeans(pca_matrix, na.rm = TRUE)
  col_sds <- apply(pca_matrix, 2, sd, na.rm = TRUE)
  col_sds[col_sds == 0] <- 1  # Avoid division by zero

  # Center and scale
  pca_matrix <- sweep(pca_matrix, 2, col_means, "-")
  pca_matrix <- sweep(pca_matrix, 2, col_sds, "/")

  # Perform Probabilistic PCA (handles missing values)
  nComponents <- min(nComponents, ncol(pca_matrix) - 1, nrow(pca_matrix) - 1)
 
  ppca_result <- tryCatch({
    pcaMethods::pca(pca_matrix,
        method = "ppca",
        nPcs = nComponents,
        seed = 42,
        maxIterations = 500)
  }, error = function(e) { 
    return(NULL)
  })

  if(is.null(ppca_result)) { 
    return("NO;PPCA computation failed")
  }
 

  # Extract scores and loadings
  pca_scores <- pcaMethods::scores(ppca_result)
  pca_loadings <- pcaMethods::loadings(ppca_result)
  var_explained <- ppca_result@R2cum
  var_each <- c(var_explained[1], diff(var_explained)) * 100

  # Create results data frame
  pca_df <- data.frame(
    record_id = record_ids,
    PC1 = pca_scores[, 1],
    PC2 = pca_scores[, 2],
    stringsAsFactors = FALSE
  )

  # Add more PCs if available
  for(i in 3:min(nComponents, ncol(pca_scores))) {
    pca_df[[paste0("PC", i)]] <- pca_scores[, i]
  }

  # Add color variable
  if(!is.null(colorVar)) {
    pca_df$colorVar <- colorVar
  } else {
    pca_df$colorVar <- "All"
  }

  # Add shape variable (optional)
  if(!is.null(shapeVar)) {
    pca_df$shapeVar <- shapeVar
  }

  # Calculate distance from centroid for outlier detection
  pc_cols <- paste0("PC", 1:min(4, ncol(pca_scores)))
  centroid <- colMeans(pca_df[, pc_cols], na.rm = TRUE)
  pca_df$distance <- sqrt(rowSums(sweep(pca_df[, pc_cols], 2, centroid)^2, na.rm = TRUE))

  # Identify outliers
  pca_df <- pca_df[order(-pca_df$distance), ]
  pca_df$is_outlier <- FALSE
  pca_df$is_outlier[1:min(highlightOutliers, nrow(pca_df))] <- TRUE

  # Get display names for colorBy
  colorByDisplay <- colorBy
  if(colorBy %in% meta_group$column) {
    colorByDisplay <- meta_group$display[meta_group$column == colorBy]
  }

  # Get display name for shapeBy
  shapeByDisplay <- shapeBy
  if(nzchar(shapeBy) && shapeBy %in% meta_group$column) {
    shapeByDisplay <- meta_group$display[meta_group$column == shapeBy]
  }

  # Prepare axis labels
  xlab <- paste0("PC", pcX, " (", round(var_each[pcX], 1), "%)")
  ylab <- paste0("PC", pcY, " (", round(var_each[pcY], 1), "%)")

  # Get PC columns for selected axes
  pcX_col <- paste0("PC", pcX)
  pcY_col <- paste0("PC", pcY)

  # Check if columns exist
  if(!pcX_col %in% names(pca_df) || !pcY_col %in% names(pca_df)) {
    return("NO;Selected PC components not available")
  }

  # Build interactive plot
  pca_df$x <- pca_df[[pcX_col]]
  pca_df$y <- pca_df[[pcY_col]]

  # Create hover text (include shape line only if shapeVar is set)
  if(!is.null(shapeVar)) {
    pca_df$hover_text <- paste0(
      "Donor: ", pca_df$record_id,
      "<br>", colorByDisplay, ": ", pca_df$colorVar,
      "<br>", shapeByDisplay, ": ", pca_df$shapeVar,
      "<br>PC", pcX, ": ", round(pca_df$x, 3),
      "<br>PC", pcY, ": ", round(pca_df$y, 3)
    )
  } else {
    pca_df$hover_text <- paste0(
      "Donor: ", pca_df$record_id,
      "<br>", colorByDisplay, ": ", pca_df$colorVar,
      "<br>PC", pcX, ": ", round(pca_df$x, 3),
      "<br>PC", pcY, ": ", round(pca_df$y, 3)
    )
  }

  # Determine if colorVar is numeric or categorical
  is_numeric_color <- is.numeric(pca_df$colorVar) && length(unique(pca_df$colorVar)) > 10

  # Shape palette (used only if shapeVar present). 16 codes cover up to 16 levels.
  shape_palette <- c(16, 17, 15, 18, 8, 4, 3, 5, 7, 9, 10, 11, 12, 13, 14, 19)

  if(is_numeric_color) {
    # Continuous color scale
    p <- ggplot(pca_df, aes(x = x, y = y, color = colorVar, text = hover_text)) +
      scale_color_viridis_c(name = colorByDisplay)
  } else {
    # Categorical color (Okabe-Ito colour-blind safe palette)
    okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
    p <- ggplot(pca_df, aes(x = x, y = y, color = as.factor(colorVar), text = hover_text)) +
      scale_color_manual(values = okabe_ito, name = colorByDisplay)
  }

  if(!is.null(shapeVar)) {
    n_shape_levels <- length(unique(stats::na.omit(pca_df$shapeVar)))
    shape_vals <- shape_palette[seq_len(min(n_shape_levels, length(shape_palette)))]
    p <- p + geom_point(aes(shape = as.factor(shapeVar)), size = 3, alpha = 0.7) +
      scale_shape_manual(values = shape_vals, name = shapeByDisplay)
  } else {
    p <- p + geom_point(size = 3, alpha = 0.7)
  }

  p <- p +
    labs(x = xlab, y = ylab, title = paste0("Phenotype PCA (n=", nrow(pca_df), ")")) +
    theme_minimal() +
    theme(legend.position = "right")

  # Add outlier labels
  outlier_df <- pca_df[pca_df$is_outlier, ]
  if(nrow(outlier_df) > 0) {
    p <- p + geom_text(data = outlier_df,
                       aes(label = record_id),
                       size = 3,
                       vjust = -0.5,
                       hjust = 0.5,
                       color = "black",
                       fontface = "bold")
  }

  # Convert to plotly
 
  pl <- tryCatch({
    ggplotly(p, tooltip = "text")
  }, error = function(e) {
 
    return(NULL)
  })
  if(is.null(pl)) return("NO;Plotly conversion failed")

  # Save interactive plot JSON
  plot_filename <- paste0("pca_plot_", sample(1:100000, 1), ".json")
 
  json_plot <- plotly_json(pl, jsonedit = FALSE, pretty = TRUE)
  cat(json_plot, file = plot_filename)

  # Prepare loadings data for feature importance
  loadings_df <- data.frame(
    feature = rownames(pca_loadings),
    stringsAsFactors = FALSE
  )
  for(i in 1:ncol(pca_loadings)) {
    loadings_df[[paste0("PC", i)]] <- pca_loadings[, i]
  }

  # Add display names
  loadings_df$display <- meta_group$display[match(loadings_df$feature, meta_group$column)]
  loadings_df$display[is.na(loadings_df$display)] <- loadings_df$feature[is.na(loadings_df$display)]

  # Calculate importance (absolute loading values)
  loadings_df$importance <- abs(loadings_df[[paste0("PC", pcX)]]) + abs(loadings_df[[paste0("PC", pcY)]])
  loadings_df <- loadings_df[order(-loadings_df$importance), ]

  # Save loadings
  loadings_filename <- paste0("pca_loadings_", sample(1:100000, 1), ".csv")
  write.csv(loadings_df, loadings_filename, row.names = FALSE)

  # Save scores for download
  scores_filename <- paste0("pca_scores_", sample(1:100000, 1), ".csv")
  write.csv(pca_df[, c("record_id", paste0("PC", 1:min(nComponents, ncol(pca_scores))), "colorVar", "distance", "is_outlier")],
            scores_filename, row.names = FALSE)

  # Prepare variance explained data
  var_df <- data.frame(
    PC = 1:length(var_each),
    variance = var_each,
    cumulative = var_explained * 100
  )
  var_filename <- paste0("pca_variance_", sample(1:100000, 1), ".csv")
  write.csv(var_df, var_filename, row.names = FALSE)

  # Prepare donor profiles (vectorized: compute all z-scores at once)
  feat_cols <- numeric_cols[valid_cols]
  z_matrix <- sweep(as.matrix(dat[, feat_cols, drop = FALSE]), 2, col_means, "-")
  z_matrix <- sweep(z_matrix, 2, col_sds, "/")
  rownames(z_matrix) <- dat$record_id

  # Pre-compute display names lookup
  display_lookup <- meta_group$display[match(feat_cols, meta_group$column)]
  display_lookup[is.na(display_lookup)] <- feat_cols[is.na(display_lookup)]

  # Build profiles: for each donor, pick top 10 features by |z-score|
  donor_profiles <- list()
  for(i in 1:nrow(pca_df)) {
    donor_id <- pca_df$record_id[i]
    zs <- z_matrix[donor_id, ]
    valid <- !is.na(zs)
    if(sum(valid) > 0) {
      zs <- zs[valid]
      top_idx <- head(order(-abs(zs)), 10)
      donor_profiles[[donor_id]] <- data.frame(
        feature = feat_cols[valid][top_idx],
        z_score = unname(zs[top_idx]),
        display = display_lookup[valid][top_idx],
        stringsAsFactors = FALSE
      )
    }
  }

  profiles_filename <- paste0("pca_profiles_", sample(1:100000, 1), ".json")
 
  write_json(donor_profiles, profiles_filename, auto_unbox = TRUE)

  # Return result
  res <- paste0("RES-OK;",
                plot_filename, ";",
                loadings_filename, ";",
                scores_filename, ";",
                var_filename, ";",
                profiles_filename, ";",
                nrow(pca_df), ";",
                ncol(pca_data))
 
  print(res)
  return(res)
}

################################################################################

#' Generate static PCA plot for download
#'
#' Creates publication-quality static PCA plots in multiple formats
#'
#' @param colorBy Variable to color points by
#' @param pcX Principal component for X-axis
#' @param pcY Principal component for Y-axis
#' @param format Output format: "pdf", "png", or "svg"
#' @param width Plot width in inches
#' @param height Plot height in inches

getPhenotypePCAStatic <- function(colorBy = "diagnosis",
                                   pcX = 1,
                                   pcY = 2,
                                   format = "pdf",
                                   width = 8,
                                   height = 6,
                                   showOutlierLabels = "yes",
                                   highlightOutliers = 10,
                                   shapeBy = "") {

  
  pcX <- as.numeric(pcX)
  pcY <- as.numeric(pcY)
  width <- as.numeric(width)
  height <- as.numeric(height)
  highlightOutliers <- as.numeric(highlightOutliers)

  # Load data (normalized for PCA)
     dat <- read.csv(paste0(other.tables.path, "display_data/metadata_sum_norm.csv"))
  meta_group <- read.csv(paste0(other.tables.path, "display_interface/meta_groups.csv"))
    # Apply filters
  if(file.exists("donors.rds")) {
    donors <- readRDS("donors.rds")
    dat <- dat[dat$record_id %in% donors, ] 
  }
  if(file.exists("rmphemo.rds")) {
    rmphemo <- readRDS("rmphemo.rds")
    dat <- dat[, !names(dat) %in% rmphemo] 
  }

  if(nrow(dat) < 10) {
 
    return("NO;Insufficient samples")
  }

  # Prepare data
  numeric_cols <- meta_group$column[meta_group$type == "cont" & meta_group$column %in% names(dat)]
 
  colorVar <- if(colorBy %in% names(dat)) dat[[colorBy]] else "All"
  shapeVar <- if(nzchar(shapeBy) && shapeBy %in% names(dat)) dat[[shapeBy]] else NULL

  record_ids <- dat$record_id
  pca_data <- dat[, numeric_cols, drop = FALSE]
  valid_cols <- colSums(!is.na(pca_data)) >= 50
  pca_data <- pca_data[, valid_cols, drop = FALSE]
  
  if(ncol(pca_data) < 5) { 
    return("NO;Insufficient features")
  }

  # Scale data
  pca_matrix <- as.matrix(pca_data)
  col_means <- colMeans(pca_matrix, na.rm = TRUE)
  col_sds <- apply(pca_matrix, 2, sd, na.rm = TRUE)
  col_sds[col_sds == 0] <- 1
  pca_matrix <- sweep(pca_matrix, 2, col_means, "-")
  pca_matrix <- sweep(pca_matrix, 2, col_sds, "/")
 

  # PPCA
  nComponents <- min(10, ncol(pca_matrix) - 1, nrow(pca_matrix) - 1)
 
  ppca_result <- tryCatch({
    pcaMethods::pca(pca_matrix, method = "ppca", nPcs = nComponents, seed = 42)
  }, error = function(e) {
 
    return(NULL)
  })

  if(is.null(ppca_result)) {
 
    return("NO;PCA computation failed")
  }
 

  pca_scores <- pcaMethods::scores(ppca_result)
  var_explained <- ppca_result@R2cum
  var_each <- c(var_explained[1], diff(var_explained)) * 100
  
  # Create data frame
    pca_df <- data.frame(
    record_id = record_ids,
    x = pca_scores[, pcX],
    y = pca_scores[, pcY],
    colorVar = colorVar,
    stringsAsFactors = FALSE
  )
  if(!is.null(shapeVar)) {
    pca_df$shapeVar <- shapeVar
  }
 

  # Calculate outliers
  centroid <- c(mean(pca_df$x, na.rm = TRUE), mean(pca_df$y, na.rm = TRUE))
  pca_df$distance <- sqrt((pca_df$x - centroid[1])^2 + (pca_df$y - centroid[2])^2)
  pca_df <- pca_df[order(-pca_df$distance), ]
  pca_df$is_outlier <- FALSE
  pca_df$is_outlier[1:min(highlightOutliers, nrow(pca_df))] <- TRUE
 

  # Labels
  colorByDisplay <- meta_group$display[meta_group$column == colorBy]
  if(length(colorByDisplay) == 0) colorByDisplay <- colorBy
  shapeByDisplay <- shapeBy
  if(nzchar(shapeBy) && shapeBy %in% meta_group$column) {
    shapeByDisplay <- meta_group$display[meta_group$column == shapeBy]
  }
  xlab <- paste0("PC", pcX, " (", round(var_each[pcX], 1), "%)")
  ylab <- paste0("PC", pcY, " (", round(var_each[pcY], 1), "%)")


  # Build plot
  is_numeric_color <- is.numeric(pca_df$colorVar) && length(unique(pca_df$colorVar)) > 10

  shape_palette <- c(16, 17, 15, 18, 8, 4, 3, 5, 7, 9, 10, 11, 12, 13, 14, 19)

  if(is_numeric_color) {
    p <- ggplot(pca_df, aes(x = x, y = y, color = colorVar)) +
      scale_color_viridis_c(name = colorByDisplay)
  } else {
    # Okabe-Ito colour-blind safe palette
    okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
    p <- ggplot(pca_df, aes(x = x, y = y, color = as.factor(colorVar))) +
      scale_color_manual(values = okabe_ito, name = colorByDisplay)
  }

  if(!is.null(shapeVar)) {
    n_shape_levels <- length(unique(stats::na.omit(pca_df$shapeVar)))
    shape_vals <- shape_palette[seq_len(min(n_shape_levels, length(shape_palette)))]
    p <- p + geom_point(aes(shape = as.factor(shapeVar)), size = 3, alpha = 0.7) +
      scale_shape_manual(values = shape_vals, name = shapeByDisplay)
  } else {
    p <- p + geom_point(size = 3, alpha = 0.7)
  }

  p <- p +
    labs(x = xlab, y = ylab,
         title = paste0("Phenotype PCA (n=", nrow(pca_df), ")")) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
 

  # Add outlier labels
  if(showOutlierLabels == "yes") {
    outlier_df <- pca_df[pca_df$is_outlier, ]
    if(nrow(outlier_df) > 0) {
      p <- p + geom_text_repel(
        data = outlier_df,
        aes(label = record_id),
        size = 3,
        fontface = "bold",
        color = "black",
        box.padding = 0.5,
        max.overlaps = 20
      ) 
    }
  }

  # Save plot
  filename <- paste0("pca_static_", sample(1:100000, 1), ".", format) 

  tryCatch({
    if(format == "pdf") {
      ggsave(filename, p, width = width, height = height, device = "pdf")
    } else if(format == "png") {
      ggsave(filename, p, width = width, height = height, dpi = 300, device = "png")
    } else if(format == "svg") {
      ggsave(filename, p, width = width, height = height, device = "svg")
    } 
  }, error = function(e) {
 paste("ERROR saving plot:", e$message)
  })

  res <- paste0("RES-OK;", filename) 
  return(res)
}

################################################################################


getPheno_fun <- function (){
   
require(jsonlite)
filter <- read_json("pheno.json"); 
filter = unlist(filter)
filter = filter[filter=="FALSE"]
rm = names(filter)
saveRDS(rm,"rmphemo.rds")

   
  res <- paste0(length(rm), "RES-OK");
  
  return(res)
}

################################################################################

