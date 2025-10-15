# Load all required functions from GitHub
# source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/multi_package_installer.R")
# source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/bsa_theme.R")
# source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/import_vcfdata.R")
# source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/analyze_vcfdata.R")
# source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/window_peak_interval.R")
# source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/window_homozygosity_compute.R")
# source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/window_homozygosity_plot.R")
# source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/window_homozygosity_single.R")
# source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/window_homozygosity_multi.R")

src_base <- "https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/"
for (f in c(
  "multi_package_installer.R","bsa_theme.R","import_vcfdata.R","analyze_vcfdata.R","peak_interval.R",
  "window_homozygosity_compute.R", "window_homozygosity_plot.R", 
  "window_homozygosity_single.R", "window_homozygosity_multi.R"
  )) source(file.path(src_base, f), local = TRUE)

message("BSA-Seq environment ready. Call window_homozygosity_multi() when you’re set.")

# required packages
packages <- c(
  "reshape2", "readxl", "BiocManager", "zoo", "plyr", "GlobalOptions", "shape", "scales",
  "tidyverse", "openxlsx", "stringr", "IRanges", "magrittr", "data.table", "naturalsort", "locfit", "rlang"
)

# Install required packages
Install_multi_package_bz(packages)