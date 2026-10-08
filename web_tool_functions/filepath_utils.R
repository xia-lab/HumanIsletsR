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
## Clusters exist for the donors with RNA-seq or proteomics; other donors only have a cluster PREDICTED from phenotypes
## (cluster_live_v2/donor_cluster.csv: source, reliability). The page sends the variable as "cluster_v2.<level>"
## (level = high | highmed | all) to add the predicted donors at that reliability or better; plain "cluster_v2" =
## clustered donors only. Predictions use phenotypes, not omics, so this is for omics analyses only.
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
## set meta$cluster_v2 = clustered donors + predicted donors at the level (no level: meta unchanged)
.fill_cluster_v2 <- function(meta, level, id.col = "record_id"){
  if(is.null(level) || is.null(meta) || !(id.col %in% colnames(meta))) return(meta)
  dc <- .cluster_v2_table(); if(is.null(dc)) return(meta)
  dc <- dc[dc$source == "omics" | (dc$source == "predicted" & dc$reliability %in% .cluster_v2_reliabilities(level)), ]
  meta$cluster_v2 <- dc$cluster[match(as.character(meta[[id.col]]), dc$donor_id)]
  meta
}
