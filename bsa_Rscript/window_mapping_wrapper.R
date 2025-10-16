src_base <- "https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/"
for (f in c(
  "multi_package_installer.R","bsa_theme.R","import_vcfdata.R","analyze_vcfdata.R",
  "peak_interval.R", "window_mapping_compute.R", "window_mapping_compute_auto.R", 
  "window_mapping_plot.R", "window_mapping_pipeline.R", "window_mapping_pipeline_auto.R"
)) source(file.path(src_base, f), local = TRUE)

message("BSA-Seq environment ready. Call window_bsa_auto() when you’re set.")

# required packages
packages <- c(
  "reshape2", "readxl", "BiocManager", "zoo", "plyr", "GlobalOptions", "shape", "scales",
  "tidyverse", "openxlsx", "stringr", "IRanges", "magrittr", "data.table", "naturalsort", "locfit", "rlang"
)

# Install required packages
Install_multi_package_bz(packages)