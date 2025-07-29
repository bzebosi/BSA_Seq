# Load all required functions from GitHub
source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/multi_package_installer.R")
source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/bsa_theme.R")
source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/import_vcfdata.R")
source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/analyze_vcfdata.R")
source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/plot_vcfdata.R")
source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/run_af_only.R")
source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/run_ed_only.R")
source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/run_g_only.R")
source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/run_histogram_only.R")
source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/run_simpval_only.R")
source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/run_bsa_pipeline.R")
source("https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/run_bsa_auto.R")

# Install required packages
packages <- c(
  "reshape2", "readxl", "BiocManager", "zoo", "plyr", "GlobalOptions", "shape", "scales",
  "tidyverse", "openxlsx", "stringr", "IRanges", "magrittr", "data.table", "naturalsort", "locfit", "rlang"
)
Install_multi_package_bz(packages)