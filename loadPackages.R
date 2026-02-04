#########################################################
# Function for installing and load the packages         #
# (The packages are installed only if they are missing) #
#########################################################
instload_packages <- function() {

  cran_packages <- c(
    "remotes", "devtools", "htmlTable", "parallel",  "TMB",
    "microbenchmark", "bamlss", "scoringRules", "boot", "mgcViz",
    "stringr", "BMisc", "Gmisc", "fields", "mvnfast", "grid",
    "magrittr", "ggpubr", "gridExtra", "ggplot2", "lattice",
    "dplyr", "plyr", "tidyr", "scales", "ggh4x", "ggnewscale",
    "lubridate",  "glue"
  )

  github_packages <- c(
    "SCM"       = "VinGioia90/SCM",
    "electBook" = "mfasiolo/electBook",
    "mvnchol"   = "meteosimon/mvnchol"
  )

  installed <- rownames(installed.packages())

  # Install missing CRAN packages
  missing_cran <- setdiff(cran_packages, installed)
  if ( length(missing_cran) ) {
    install.packages(missing_cran)
  }

  # Install missing GitHub packages
  for (pkg in names(github_packages)) {
    if ( !(pkg %in% installed) ) {
      remotes::install_github(github_packages[pkg])
    }
  }

  # Load all packages
  all_packages <- c(cran_packages,
                    "Rcpp", "SparseChol",
                    names(github_packages))
  for (pkg in all_packages) {
    print(paste("Loading ", pkg))
    library(pkg, character.only = TRUE)
  }
}


