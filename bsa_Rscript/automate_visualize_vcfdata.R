automate_visualize_vcfdata <- function(data, plots_dir, output_dir, wt_list, mt_list, prefixes, pattern = pattern, save_results = TRUE, 
                            min_DP = 1, min_QUAL = 1, threshold = -log10(0.05) * 10, n_prop = 0.1, 
                            hwidth = 30, hheight = 16, width = 26, height = 8, dpi = 300, device = "tiff", plot_type = "pval") {
  # Initialize a list to store data for each WT-MT combination
  all_results <- list()
  for (i in seq_along(prefixes)) {
    pfix <- prefixes[i]
    # print message for the inbred
    print(paste("Processing inbred:", pfix))
    
    # Check if single WT-MT combination or multiple combinations
    if (length(wt_list) == 1 && length(mt_list) == 1) {
      # Directly process the single WT-MT combination
      wt <- wt_list[1]
      mt <- mt_list[1]
      # Create a sub-directory for the WT-MT combination for the plots
      sub_plots <- file.path(plots_dir, paste0(pfix, "_", wt, "_", mt))
      if (!dir.exists(sub_plots)) dir.create(sub_plots, recursive = TRUE)
      # Print message
      message(paste("Processing single WT:", wt, "and single MT:", mt, "in", sub_plots))
      
      # Call the visualize_vcfdata function and assign the result
      result <- paste0(pfix,"_", wt, "_", mt)
      print(result)
      # Call the visualize_vcfdata function and assign the result
      result_data <- load_analze_visualize_vcfdata(vcf_dir=vcf_dir, prefix=prefix, pattern = pattern, output_dir = output_dir, plots_dir = plots_dir , save_results = save_results,
                                                   Genotype=c(wt=wt, mt=mt), min_DP = min_DP, min_QUAL = min_QUAL, threshold = threshold, n_prop = n_prop, 
                                                   hwidth = hwidth, hheight = hheight, width = width, height = height, dpi = dpi, device = device, plot_type = plot_type)
      # Store the result in the list
      if (!is.null(result_data)) {
        all_results[[result]] <- result_data
        print(paste("Stored result for:", result))
      } else {
        warning(paste("Failed to process WT:", wt, "and MT:", mt))
      }
      # Print message
      message(paste("Completed processing for WT:", wt, "and MT:", mt, "for", pfix, "in", sub_plots))
    } else{
      # Combine and ensure no duplicates
      all_list <- unique(c(wt_list, mt_list))
      if (length(all_list) < 2) {
        warning("Unique elements less than in wt_list and mt_list to create pairs.")
        next}
      # Combine the lists into one and generate unique pairs
      unique_pairs <- combn(all_list, 2, simplify = FALSE)
      # loop through the unique pairs
      for (pair in unique_pairs) {
        wt <- pair[1]
        mt <- pair[2]
        sub_plots <- file.path(plots_dir, paste0(pfix, "_", wt, "_", mt))
        if (!dir.exists(sub_plots)) dir.create(sub_plots, recursive = TRUE)
        # Print message
        message(paste("Processing WT:", wt, "and MT:", mt, "in", sub_plots))
        
        # Call the visualize_vcfdata function and assign the result
        result <- paste0(pfix,"_", wt, "_", mt)
        print(result)
        result_data <- load_analze_visualize_vcfdata(vcf_dir=vcf_dir, prefix=prefix, pattern = pattern, output_dir = output_dir, plots_dir = plots_dir , save_results = save_results,
                                                     Genotype=c(wt=wt, mt=mt), min_DP = min_DP, min_QUAL = min_QUAL, threshold = threshold, n_prop = n_prop, 
                                                     hwidth = hwidth, hheight = hheight, width = width, height = height, dpi = dpi, device = device, plot_type = plot_type)
        # Store the result in the list
        if (!is.null(result_data)) {
          all_results[[result]] <- result_data
          print(paste("Stored result for:", result))
        } else {
          warning(paste("Failed to process WT:", wt, "and MT:", mt))
        }
      }
    }
  }
  return(all_results)
}


# vcf_dir <- "path/data"
# plots_dir <- "path/plots"
# output_dir <-  "path/post_out"
# wt_list <- c("MB8")
# mt_list <- c("MB9")
# prefix <- c("p39", "b73")
# pattern <- "snps\\.txt$"


# results <- automate_visualize_vcfdata(vcf_dir=vcf_dir, plots_dir=plots_dir , output_dir=output_dir, wt_list=wt_list, mt_list=mt_list, prefix=prefix, pattern = pattern, save_results = TRUE, 
#                                       min_DP = 1, min_QUAL = 1, threshold = -log10(0.05) * 10, n_prop = 0.1, 
#                                       hwidth = 30, hheight = 16, width = 26, height = 8, dpi = 300, device = "tiff", plot_type = "pval")
