# Run required dependencies and available on github
source(file = "https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/multi_package_installer.R")
source(file = "https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/bsa_theme.R")
source(file = "https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/find_locus.R")

data <- "/Users/zebosi/Documents/BSA/rlz1/data/"
plots_dir <- "/Users/zebosi/Documents/BSA/rlz11/rlz1_plots"
wt_list <- c("rzl4", "MB3", "MB4")
mt_list <- c("MB3", "MB4", "MB5")
inbreds <- c("B73")
prefixes <- c("b73")

for (i in seq_along(inbreds)) {
  inbx <- inbreds[i]
  pfix <- prefixes[i]
  print(paste("Processing inbred:", inbx))
  
  # Combine the lists into one and generate unique pairs
  all_list <- unique(c(wt_list, mt_list)) # Combine and ensure no duplicates
  unique_pairs <- combn(all_list, 2, simplify = FALSE) # Generate unique pairs
  
  # Iterate over the unique pairs
  for (pair in unique_pairs) {
    wt <- pair[1]
    mt <- pair[2]
    print(c(wt, mt))
    # Create subdirectory for plots
    sub_plots <- file.path(plots_dir, paste0(inbx, "_", wt, "_", mt))
    print(sub_plots)
    if (!dir.exists(sub_plots)) dir.create(sub_plots, recursive = TRUE)
    # Print status
    message(paste("Processing WT:", wt, "and MT:", mt, "in", sub_plots))
    # Call the visualize_vcfdata function and assign the result
    assign(inbx, visualize_vcfdata(vcf_dir=data, plots_dir = sub_plots, 
                                   inbred = inbx, prefix=pfix, Genotype=c(wt=wt, mt=mt), nn_prop = 0.12, 
                                   width = 26, height = 8, dpi = 300))
    
    message(paste("Completed processing for WT:", wt, "and MT:", mt, "for", inbx, "in", sub_plots))
  }
}

