setPaths <- function(){

  # set sqlite file path
  if(file.exists("/home/glassfish/sqlite/")){
    sqlite.path <<- "/home/glassfish/sqlite/"; # public server
    h5.path <<- "/home/glassfish/hdf5/";
    h5.v2.path <<- "/home/glassfish/hdf5/hdf5_V2/";
    other.tables.path <<- "/home/glassfish/resources/humanislets/"
  }else if(file.exists("/Users/xia/Dropbox/sqlite/")){
    sqlite.path <<- "/Users/xia/Dropbox/sqlite/"; #xia local
    h5.path <<- "/Users/xia/Dropbox/sqlite/";
    h5.v2.path <<- "/Users/xia/Dropbox/hdf5_V2/";
  } else if(file.exists("/Users/lzy/humanislet/sqlite")){
    sqlite.path <<- "/Users/lzy/humanislet/sqlite/"; #ly local
    h5.path <<- "/Users/lzy/humanislet/hdf5/";
    h5.v2.path <<- "/Users/lzy/humanislet/hdf5_V2/";
    other.tables.path <<- "/Users/lzy/humanislet/humanislets/"
  }else {
    sqlite.path <<- "/Users/jessicaewald/sqlite/"; # ewald local
    h5.path <<- "/Users/jessicaewald/hdf5/"
    h5.v2.path <<- "/Users/jessicaewald/hdf5_V2/"
    other.tables.path <<- "/Users/jessicaewald/Desktop/RestTest/resources/humanislets/"
  }

}

## ---- Multi-omics cluster with predicted donors (Omics View / Multi-omics View "Cluster donors" choice) ----
## Multi-omics cluster (cluster_v2) for ONE calculation. The tables hold every donor's cluster: the omics donors' and
## the predicted ones, with cluster_v2_source (omics | predicted) and cluster_v2_reliability (high | medium | low).
## A calculation keeps the donors of the selected set and blanks the others' cluster. The page sends "cluster_v2.<level>"
## (level = high | highmed | all) to add predicted donors at that reliability or better; plain "cluster_v2" or no choice =
## "Omics donors only", the default everywhere (user 2026-10-09).
.cluster_v2_level <- function(var){
  if(is.null(var) || length(var) != 1 || is.na(var)) return(NULL)
  m <- regmatches(var, regexec("^cluster_v2[.](high|highmed|all)$", var))[[1]]
  if(length(m) == 2) m[2] else NULL
}
.cluster_v2_base <- function(var) if(is.null(.cluster_v2_level(var))) var else "cluster_v2"
.cluster_v2_reliabilities <- function(level) switch(level, high = "high", highmed = c("high", "medium"), all = c("high", "medium", "low"))
.cluster_v2_table <- function(){
  f <- paste0(other.tables.path, "donor_clustering_pipeline_v2/cluster_live_v2/donor_cluster.csv")
  if(!file.exists(f)) return(NULL)
  read.csv(f, stringsAsFactors = FALSE, colClasses = c(donor_id = "character"))
}
## meta$cluster_v2 for the selected set: omics donors + predicted donors at the level (level NULL = omics donors only).
## Uses the table's cluster_v2_source / cluster_v2_reliability; a table without them is filled from donor_cluster.csv.
.fill_cluster_v2 <- function(meta, level = NULL, id.col = "record_id"){
  if(is.null(meta) || !(id.col %in% colnames(meta))) return(meta)
  rel <- if(is.null(level)) character(0) else .cluster_v2_reliabilities(level)
  if(all(c("cluster_v2", "cluster_v2_source", "cluster_v2_reliability") %in% colnames(meta))){
    keep <- meta$cluster_v2_source %in% "omics" | (meta$cluster_v2_source %in% "predicted" & meta$cluster_v2_reliability %in% rel)
    meta$cluster_v2[!keep] <- NA
    return(meta)
  }
  if(is.null(level) && !("cluster_v2" %in% colnames(meta))) return(meta)
  dc <- .cluster_v2_table(); if(is.null(dc)) return(meta)
  dc <- dc[dc$source == "omics" | (dc$source == "predicted" & dc$reliability %in% rel), ]
  meta$cluster_v2 <- dc$cluster[match(as.character(meta[[id.col]]), dc$donor_id)]
  meta
}
## a phenotype table read for ONE calculation (all_features.csv & co., id column sample_id): Multi-omics cluster of the
## selected donors (level NULL = omics donors only); the two bookkeeping columns are never offered as phenotypes
.cluster_v2_pheno <- function(pheno, level = NULL, id.col = "sample_id"){
  pheno <- .fill_cluster_v2(pheno, level, id.col)
  pheno[, setdiff(colnames(pheno), c("cluster_v2_source", "cluster_v2_reliability")), drop = FALSE]
}
