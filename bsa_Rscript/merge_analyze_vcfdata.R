#' Merge and Analyze VCF Data
#' @description Merges processed wild-type and mutant VCF data and identifies unique SNPs.
#' @param geno_data A list of processed VCF data for each genotype (output of `import_vcfdata`).
#' @param prefix Prefix for output file names.
#' @param save_results Logical. If TRUE, saves results to the output directory.
#' @param output_dir Directory to save results.
#' @return A list containing merged SNP data and unique SNPs for each genotype.
#' @export
merge_analyze_vcfdata <- function(geno_data, prefix, save_results = FALSE, output_dir = "post_analysis", threshold = -log10(0.05) * 10){
  # Ensure required packages are installed
  # Load required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required but not installed.")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required but not installed.")
  
  # Ensure required genotypes are available
  if (is.null(geno_data$wt) || is.null(geno_data$mt)) {stop("Both wild-type and mutant datasets are required.")}
  wt_data <- geno_data$wt
  mt_data <- geno_data$mt
  
  # Merge the two filtered datasets by CHROM and POS
  message("Merging wild-type and mutant datasets by shared snp positions...")
  wt_mt <- merge(wt_data, mt_data, by = c("CHROM", "POS"))
  wt_mt[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
  wt_mt <- wt_mt[order(CHROM, POS)]
  message("Number of shared SNPs (wt_mt): ", nrow(wt_mt))
  
  # Anti-join to find unique positions in mutant and wild-type datasets
  id_unique_snps <- function(data1, data2, label) {
    message(paste("Keeping snps unique to", label, "dataset..."))
    unique_snps <- anti_join(data1, data2, by = c("CHROM", "POS"))
    unique_snps[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
    unique_snps <- unique_snps[order(CHROM, POS)]
    message("Number of SNPs unique to ", label, " dataset: ", nrow(unique_snps))
    return(unique_snps)
  }
  
  # Identify unique SNPs
  ant_wt <- id_unique_snps(wt_data, mt_data, "wildtype")
  ant_mt <- id_unique_snps(mt_data, wt_data, "mutant")
  
  get_ems <-function(data, ref_col, alt_col, label){
    message(paste("Filtering potential canonical ems (G/A to C/T)", ref_col,"/",alt_col, "in the ", label, " dataset..."))
    ems_variants <- data[(get(ref_col) == "G" & get(alt_col) == "A") | (get(ref_col) == "C" & get(alt_col) == "T")]
    ems_variants[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
    ems_variants <- ems_variants[order(CHROM, POS)]
    message("Number of potential canonical EMS SNPs (G/A to C/T) unique to " , label, " dataset: ", nrow(ems_variants))
    return(ems_variants)
  }
  
  ant_wt_ems <- get_ems(wt_data, "wt_REF", "wt_ALT", "wildtype")
  ant_mt_ems <- get_ems(mt_data, "mt_REF", "mt_ALT", "mutant" )
  
  
  # Calculate allele frequency differences
  message("Calculating allele frequency differences and statistical metrics...")
  wt_mt[, AF_diff := wt_AF - mt_AF]
  
  # Calculate Euclidean Distance (ED) and its fourth power (ED4)
  message("Calculating allele Euclidean Distance (ED) and its fourth power (ED4)...")
  # Calculate Euclidean Distance (ED)
  message("Calculating Euclidean Distance (ED)...")
  wt_mt[, ED := sqrt(2) * abs(wt_AF - mt_AF)]
  
  # Raise ED to the fourth power (ED4) to amplify differentiation signal
  message("Calculating ED to the fourth power (ED4)...")
  wt_mt[, ED4 := ED^4]
  
  # Calculate G-statistics for each SNP
  message("Calculating G-statistics for each SNP...")
  compute_Gstat <- function(a_ref, a_alt, b_ref, b_alt) {
    total <- (a_ref + b_ref + a_alt + b_alt)
    
    if (total == 0) {return(NA)}
    e1 <- (a_ref + b_ref) * (a_ref + a_alt) / total
    e2 <- (a_ref + b_ref) * (b_ref + b_alt) / total
    e3 <- (a_alt + b_alt) * (a_ref + a_alt) / total
    e4 <- (a_alt + b_alt) * (b_ref + b_alt) / total
    
    obs <- c(a_ref, b_ref, a_alt, b_alt)
    exp <- c(e1, e2, e3, e4)
    
    # Replace zeros in obs/exp to prevent log(0) issues
    obs[obs == 0] <- 1e-10
    exp[exp == 0] <- 1e-10
    
    return(2 * sum(obs * log(obs / exp), na.rm = TRUE))
  }
  
  wt_mt[, Gstat := mapply(compute_Gstat, a_ref=(mt_Fref + mt_Rref), a_alt=(mt_Falt + mt_Ralt), 
                          b_ref=(wt_Fref + wt_Rref), b_alt=(wt_Falt + wt_Ralt))]
  
  # Compute Chi-square p-value based on G-statistic
  wt_mt[, chi_pval := p.adjust(1 - pchisq(Gstat, df = 1), , method = "bonferroni")]
  wt_mt[, chi_log_adj.pval := -10 * log10(chi_pval)]
  
  
  # Calculate p-value using Fisher's exact test for each position
  message("Calculating p-value using Fisher's exact test for each position...")
  fisher_test_pval <- function(a_ref, a_alt, b_ref, b_alt) {
    if (anyNA(c(a_ref, a_alt, b_ref, b_alt))) return(NA)
    tryCatch(fisher.test(matrix(c(a_alt, a_ref, b_alt, b_ref), nrow = 2))$p.value, error = function(e) NA)
  }
  
  
  wt_mt[, p.value := mapply(fisher_test_pval, a_ref=(mt_Fref + mt_Rref), a_alt=(mt_Falt + mt_Ralt), 
                            b_ref=(wt_Fref + wt_Rref), b_alt=(wt_Falt + wt_Ralt))]
  
  # Adjust p-values for multiple testing
  wt_mt[, adj.p.value := p.adjust(p.value, method = "bonferroni")]
  
  # Transform adjusted p-value to -10 * log10(adj.p.value)
  wt_mt[, log_adj.pval := -10 * log10(adj.p.value)]
  
  # Identify significant SNPs based on the threshold
  fisher_significant_snps <- wt_mt[log_adj.pval > threshold]
  # Log the number of significant SNPs
  if (nrow(fisher_significant_snps) > 0) {
    message(paste("Significant SNPs based on fisher's exact test:", nrow(fisher_significant_snps)))
  } else {message("No significant SNPs based on fisher's exact test found above the threshold.")}
  
  # Identify significant SNPs based on the threshold
  chi_significant_snps <- wt_mt[chi_log_adj.pval > threshold]
  # Log the number of significant SNPs
  if (nrow(chi_significant_snps) > 0) {
    message(paste("Significant SNPs based on chisquare's test :", nrow(chi_significant_snps)))
  } else {message("No significant SNPs based on chisquare's test found above the threshold.")}
  

  # Bundle results into a list
  analyzed_merged_genodata <- list(wt_mt = wt_mt, ant_mt = ant_mt, ant_wt = ant_wt, 
                                   ant_mt_ems = ant_mt_ems, ant_wt_ems = ant_wt_ems, 
                                   fisher_significant_snps = fisher_significant_snps, 
                                   chi_significant_snps = chi_significant_snps)
  
  # Save results as an RDS file
  if(save_results){
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (is.null(output_dir)) stop("Output directory must be specified when save_results = TRUE.")
    
    message("Saving results to output directory...")
    rds_file <- file.path(output_dir, paste0(prefix, "_results.rds"))
    saveRDS(analyzed_merged_genodata, rds_file)
    
    # Save results
    for (name in names(analyzed_merged_genodata)) {
      csv_file <- file.path(output_dir, paste0(prefix, "_", name, ".csv"))
      fwrite(merged_results[[name]], file = csv_file)
    }
    message("Results saved successfully in ", output_dir)
  }
  
  return(analyzed_merged_genodata)
}

# vcf_dir <- "/Users/zebosi/Documents/BSA/S7_6508K/data/"
# plots_dir <- "/Users/zebosi/Documents/BSA/S7_6508K/plots"
# output_dir <-  "/Users/zebosi/Documents/BSA/S7_6508K/postanalysis"
# pattern <- "snps\\.txt$"
# wt <- c("S7A6508K")
# mt <- c("S7B6508K")
# prefix<-c("b73")
# min_DP <- 10
# min_QUAL <- 30
# geno_data <- import_vcfdata(vcf_dir, prefix, pattern, Genotypes=list(wt=wt, mt=mt), min_DP=5, min_QUAL=5)
# merg_data <- merge_analyze_vcfdata(geno_data, prefix, save_results = FALSE, output_dir = "post_analysis", threshold = -log10(0.05) * 10)


