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
