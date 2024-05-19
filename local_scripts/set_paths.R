# Set important file paths

setPaths <- function(){
  
  # set sqlite file path
  if(file.exists("/home/glassfish/sqlite/")){
    sqlite.path <<- "/home/glassfish/sqlite/"; # public server
    h5.path <<- "/home/glassfish/hdf5/";
    other.tables.path <<- "/home/glassfish/resources/humanislets/"
  }else if(file.exists("/Users/xia/Dropbox/sqlite/")){
    sqlite.path <<- "/Users/xia/Dropbox/sqlite/"; #xia local
    h5.path <<- "/Users/xia/Dropbox/sqlite/";
  } else if(file.exists("/Users/lzy/humanislet/sqlite")){
    sqlite.path <<- "/Users/lzy/humanislet/sqlite/"; #ly local
    h5.path <<- "/Users/lzy/humanislet/hdf5/";
    other.tables.path <<- "/Users/lzy/humanislet/humanislets/"
  }else if (file.exists("/Users/jessicaewald/sqlite/")) {
    sqlite.path <<- "/Users/jessicaewald/sqlite/"; # ewald local
    h5.path <<- "/Users/jessicaewald/hdf5/"
    other.tables.path <<- "/Users/jessicaewald/Desktop/RestTest/resources/humanislets/"
    restapi.path <<- "/Users/jessicaewald/NetbeansProjects/restxialab/"
  } else {
    print("Please set your local paths for this computer!")
  }
  
}