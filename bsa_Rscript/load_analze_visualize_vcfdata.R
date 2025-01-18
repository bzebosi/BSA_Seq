# Function to process and visualize VCF data
load_analze_visualize_vcfdata <- function(vcf_dir, prefix, pattern = "*snps.txt", output_dir = "analysis_output", plots_dir = "plots", save_results = TRUE,
                           Genotype=c(wt="wt", mt="mt"), min_DP = 1, min_QUAL = 1, threshold = -log10(0.05) * 10, n_prop = 0.1,
                           hwidth = 30, hheight = 16, width = 26, height = 8, dpi = 300, device = "tiff", plot_type = "all" ) {
  
  
  # Helper Function : Load and Process VCF data
  message("Launching the load vcfdata function pipeline...")
  geno_data <- load_vcfdata(vcf_dir, prefix, pattern, Genotype, min_DP, min_QUAL)
  
  # Merge and analyze SNPs
  message("Starting merge_analyze snps function pipeline...")
  results <- merge_analyze_snps(geno_data, prefix, threshold, save_results, output_dir)
  
  wt_mt <- results$wt_mt
  ant_wt <- results$ant_wt
  ant_mt <- results$ant_mt
  significant_snps <- results$significant_snps
  

  
  message("Generating plots...")
  # Plotting helper function
  
  
  # Generate and save plots for WT, MT, and AF_diff
  pt1 <- plot_config(wt_mt, ant_wt, ant_mt, wt, mt, plot_type)
  # Loop through the pt1 and call plot_results
  for (plots in pt1) {
    genarate_vcfplots(
      data = plots$data, column = plots$column, plot_title = plots$plot_title,
      y_title = plots$y_title, file_suffix = plots$file_suffix, is_histogram = plots$is_histogram,
      ylim = plots$ylim, is_smooth = plots$is_smooth, threshold = plots$threshold,
      lcolor = plots$lcolor, lsize = plots$lsize, psize = plots$psize, rollmedian = plots$rollmedian, nn_prop=nn_prop, 
      prefix = prefix, hwidth = hwidth, hheight = hheight, width = width, height = height, dpi = dpi, device = device, plots_dir = plots_dir, plot_type = plot_type
    )
  }
}


# data <- "path/data"
# plots_dir <- "path/lots"
# output_dir <-  "path/post_analyzisout"
# wt <- c("MB1")
# mt <- c("MB2")
# prefixes <- c("b73")

# gen <- load_analze_visualize_vcfdata(vcf_dir=data, prefix=prefixes, pattern = "*snps.txt", output_dir = output_dir, plots_dir = plots_dir , save_results = TRUE,
#                      Genotype=c(wt=wt, mt=mt), min_DP = 1, min_QUAL = 1, threshold = -log10(0.05) * 10, n_prop = 0.1, 
#                      hwidth = 30, hheight = 16, width = 26, height = 8, dpi = 300, device = "tiff", plot_type = "pval") 
