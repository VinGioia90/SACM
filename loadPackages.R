#########################################################
# Function for installing and load the packages         #
# (The packages are installed only if they are missing) #
#########################################################

instload_packages <- function(){
  neededPackages <- c("parallel", "Rcpp", "TMB", "microbenchmark",
                      "bamlss",  "scoringRules", "boot", "devtools",
                      "mgcViz", "stringr", "BMisc", "electBook",
                      "Gmisc", "glue", "fields", "mvnfast",
                      "htmlTable", "grid", "magrittr",
                      "ggpubr", "gridExtra", "ggplot2", "lattice",
                      "dplyr", "plyr", "tidyr", "scales", "ggh4x",
                      "ggnewscale", "lubridate")#, "geojsonio",
                     # "rgeos", "maptools", "gdata")

  inst_pack <- installed.packages()

  for ( j in 1 : length(neededPackages) ) {
    if ( !(neededPackages[j] %in% inst_pack[, 1]) ) {
      install.packages(neededPackages[j])
    }
    print(paste("Loading", neededPackages[j]))
    library(neededPackages[j], character.only = TRUE)
  }

  if ( !("SCM" %in% inst_pack[,1]) ) {
   devtools::install_github("VinGioia90/SCM@Master")
  }
  print("Loading SCM")
  library(SCM)
  if ( !("electBook" %in% inst_pack[,1]) ) {
    devtools::install_github("mfasiolo/electBook@Master")
  }
  print("Loading electBook")
  library(electBook)
  if ( !("mvnchol" %in% inst_pack[,1]) ) {
    devtools::install_github("meteosimon/mvnchol")
  }
  print("Loading mvnchol")
  library(mvnchol)
}


