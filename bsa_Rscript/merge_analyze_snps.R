# Function merge_analyze_snps merges wild-type and mutant SNP datasets, and calculates allele frequencies, Euclidean Distance (ED), 
# identifies significant snps associated mutation 

merge_analyze_snps <- function(geno_data, prefix, threshold, save_results, output_dir){
  # Parameter Descriptions:
  # geno_data: A list containing two data.tables: wild-type SNP data (geno_data$wt) and : mutant SNP data (geno_data$mt).
  # prefix: A string used as a prefix (e.g. b73) for output file names. prefix: 
  # threshold: A numeric value representing the -10 * log10(adjusted p-value) significance threshold for identifying significant SNPs.
  # save_results: save results or not logical input (TRUE OR FALSE)
  # output_dir: Specify path where to save results if  save_results is TRUE.
  
  if (is.null(geno_data$wt) || is.null(geno_data$mt)) {stop("Both wild-type and mutant datasets are required.")}
  wt_data <- geno_data$wt
  mt_data <- geno_data$mt
  
  # Merge the two filtered datasets by CHROM and POS
  message("Merging wild-type and mutant datasets...")
  wt_mt <- merge(wt_data, mt_data, by = c("CHROM", "POS"))
  wt_mt[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
  wt_mt <- wt_mt[order(CHROM, POS)]
  
  # Anti-join to find unique positions in mutant and wild-type datasets
  message("Identifying unique SNPs in mutant and wild-type datasets...")
  ant_mt <- anti_join(mt_data, wt_data, by = c("CHROM", "POS"))
  ant_mt[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
  ant_mt <- ant_mt[order(CHROM, POS)]
  
  ant_wt <- anti_join(wt_data, mt_data, by = c("CHROM", "POS"))
  ant_wt[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
  ant_wt <- ant_wt[order(CHROM, POS)]
  
  # Calculate allele frequency differences
  message("Calculating allele frequency differences and statistical metrics...")
  wt_mt[, AF_diff := wt_AF - mt_AF]
  
  # Calculate Euclidean Distance (ED) and its fourth power (ED4)
  wt_mt[, ED := sqrt((mt_Fref + mt_Rref - wt_Fref - wt_Rref)^2 + (mt_Falt + mt_Ralt - wt_Falt - wt_Ralt)^2)]
  wt_mt[, ED4 := ED^4]
  
  # Calculate p-value using Fisher's exact test for each position
  wt_mt[, p.value := mapply(function(a_ref, a_alt, b_ref, b_alt) {
    tryCatch({
      fisher.test(matrix(c(a_alt, a_ref, b_alt, b_ref), nrow = 2))$p.value
    }, error = function(e) NA)
  }, mt_Fref + mt_Rref, mt_Falt + mt_Ralt, wt_Fref + wt_Rref, wt_Falt + wt_Ralt)]
  
  # Adjust p-values for multiple testing
  wt_mt[, adj.p.value := p.adjust(p.value, method = "bonferroni")]
  
  # Transform adjusted p-value to -10 * log10(adj.p.value)
  wt_mt[, log_adj.pval := -10 * log10(adj.p.value)]
  
  # Identify significant SNPs based on the threshold
  significant_snps <- wt_mt[log_adj.pval > threshold]
  # Log the number of significant SNPs
  if (nrow(significant_snps) > 0) {
    message(paste("Significant SNPs:", nrow(significant_snps)))
  } else {message("No significant SNPs found above the threshold.")}
  
  # Bundle results into a list
  results <- list(wt_mt = wt_mt, ant_mt = ant_mt, ant_wt = ant_wt, significant_snps = significant_snps)
  
  # Save results as an RDS file
  if(save_results){
    message("Saving results to output directory...")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    rds_file <- file.path(output_dir, paste0(prefix, "_results.rds"))
    saveRDS(results, rds_file)
    
    # Save results
    for (name in names(results)) {
      csv_file <- file.path(output_dir, paste0(prefix, "_", name, ".csv"))
      fwrite(results[[name]], file = csv_file)
    }
    message("Results saved successfully in ", output_dir)
  }
  message("SNP analysis pipeline completed successfully and results lists created.")
  
  return(results)
}

# Example Usage
# prefix <- c("b73")
# output_dir <- "path/analysis_output"
# geno_data <- geno_data # generate used load_vcfdata function
# save_results <- TRUE
# threshold = -log10(0.05)
# results <- merge_analyze_snps(geno_data, prefix, threshold, save_results, output_dir)
