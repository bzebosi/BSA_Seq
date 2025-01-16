# Function install and load multiple packages at once
Install_multi_package_bz <- function(packages) {
  # Remove duplicates
  packages <- unique(packages)
  
  # install BiocManager if not installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # Filter out already installed packages
  not_installed_yet <- packages[!(packages %in% rownames(installed.packages()))]
  
  if (length(not_installed_yet) > 0) {
    cran_packages <- not_installed_yet[!(not_installed_yet %in% BiocManager::available())]
    bioc_packages <- not_installed_yet[not_installed_yet %in% BiocManager::available()]
    
    # Install CRAN packages
    if (length(cran_packages) > 0) {
      message("Installing CRAN packages: ", paste(cran_packages, collapse = ", "))
      install.packages(cran_packages, dependencies = TRUE)
    }
    
    # Install Bioconductor packages
    if (length(bioc_packages) > 0) {
      message("Installing Bioconductor packages: ", paste(bioc_packages, collapse = ", "))
      BiocManager::install(bioc_packages)
    }
  } else {
    message("All packages are already installed.")
  }
  
  # Load all packages
  invisible(lapply(packages, function(pkg) {
    suppressMessages(require(pkg, character.only = TRUE))
  }))
  
  message("All packages are loaded successfully.")
}


