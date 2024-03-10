

read_packages <- function(path) {
  packages <- readLines(path)
  
  # Check and install packages
  install_missing_packages <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
      install.packages(package, dependencies = TRUE, repos = "https://cran.rstudio.com/")
    }
  }
  
  # Install missing packages
  invisible(lapply(packages, install_missing_packages))
  
  # Load all packages
  invisible(lapply(packages, library, character.only = TRUE))
}
