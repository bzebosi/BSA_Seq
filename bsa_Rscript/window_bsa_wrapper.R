src_base <- "https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/"
for (f in c(
  "multi_package_installer.R","bsa_theme.R","import_vcfdata.R","analyze_vcfdata.R",
  "window_peak_interval.R", "window_bsa_plot.R", "window_bsa_pipeline.R", 
  "window_bsa_compute_one.R", "window_bsa_compute_all.R"
)) source(file.path(src_base, f), local = TRUE)

message("BSA-Seq environment ready. Call window_bsa_pipeline() when you’re set.")

# required packages
packages <- c(
  "reshape2", "readxl", "BiocManager", "zoo", "plyr", "GlobalOptions", "shape", "scales",
  "tidyverse", "openxlsx", "stringr", "IRanges", "magrittr", "data.table", "naturalsort", "locfit", "rlang"
)

# Install required packages
Install_multi_package_bz(packages)