# Run required dependencies and available on github
source(file = "https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/multi_package_installer.R")
source(file = "https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/bsa_theme.R")

#' Import VCF Data
#' @description Imports and processes VCF files for specific genotypes.
#' @param vcf_dir Directory containing the VCF files.
#' @param prefix Common prefix for the VCF files.
#' @param pattern File name pattern for VCF files.
#' @param Genotypes A named list of genotype suffixes.
#' @param min_DP Minimum depth for filtering.
#' @param min_QUAL Minimum quality score for filtering.
#' @return A list of filtered and processed data tables.
#' @export
import_vcfdata <- function (vcf_dir, prefix, pattern, Genotypes=list(wt="wildtype", mt="mutant"), min_DP=5, min_QUAL=5){
  
  vcf_list <- list.files(path = vcf_dir, pattern = pattern, full.names = TRUE)
  if (length(vcf_list) == 0) {stop("No vcf files found.")}
  
  geno_data <- list()
  
  for (genotype in names(Genotypes)) {
    file_pattern <- paste0("^", prefix, "_", Genotypes[[genotype]])
    geno_file <- vcf_list[grepl(file_pattern, basename(vcf_list))]
    
    if (length(geno_file) == 0){
      stop(paste("File for", genotype, "not found."))
    } else if (length(geno_file) > 1) {
      message(paste("Mutiple files matched", genotype, "please rename:"))
      
      # Select a file of interest from the menu
      selection <- menu(geno_file, title = paste("Select a file for", genotype))
      
      if (selection == 0) {
        stop("No selection made. Please restart and choose a valid file.")
      }
    }
    
    data <- tryCatch({fread(geno_file)}, error = function(e) {
      stop(paste("Error reading file for", genotype, ":", geno_file, ":", e$message))})
    
    if (is.null(data) || nrow(data) == 0) {
      stop(paste("File for", genotype, "is either empty or missing"))
    }
    
    message(paste("Loaded", genotype, "data with", nrow(data), "rows."))
    
    # Rename columns and convert to data.table if it is not already
    col_names <- c("CHROM", "POS", "REF", "ALT", "QUAL", "DP", "Fref", "Rref", "Falt", "Ralt")
    colnames(data) <- col_names
    if (!is.data.table(data)) data <- as.data.table(data)
    
    # Keep only SNPs
    data <- data[nchar(REF) == 1 & nchar(ALT) == 1]
    
    # Filter rows based on depth (DP) and quality (QUAL)
    data <- data[!is.na(DP) & DP > min_DP & QUAL >= min_QUAL]
    
    # Calculate allele frequency (AF)
    data[, AF := (Falt + Ralt) / (Fref + Rref + Falt + Ralt)]
    
    # Rename columns with genotype prefix
    setnames(data, old = setdiff(names(data), c("CHROM", "POS")), new = paste0(genotype, "_", setdiff(names(data), c("CHROM", "POS"))))
    
    # Store processed data in the list
    geno_data[[genotype]] <- data
  }
  
  # Ensure both wild-type and mutant data are present
  if (is.null(geno_data$wt) || is.null(geno_data$mt)) {stop("Both wild-type and mutant data are required.")}
  
  message("successfully created list for results list")
  return(geno_data)
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

# geno_data<-import_vcfdata(vcf_dir, prefix, pattern, Genotypes=list(wt=wt, mt=mt), min_DP=5, min_QUAL=5)
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
#geno_data <- import_vcfdata(vcf_dir, prefix, pattern, Genotypes=list(wt=wt, mt=mt), min_DP=5, min_QUAL=5)
# merg_data <- merge_analyze_vcfdata(geno_data, prefix, save_results = FALSE, output_dir = "post_analysis", threshold = -log10(0.05) * 10)

#' Generate VCF Data Plots
#' @description Creates SNP-based plots (scatter, smooth, histogram)
#' @param data Data frame containing VCF SNP data.
#' @param prefix Prefix for file naming.
#' @param column Name of column to plot on the y-axis.
#' @param y_title Y-axis label.
#' @param plot_title Title of the plot.
#' @param file_suffix Suffix for saved file.
#' @param ylim Y-axis limits.
#' @param is_smooth Apply smoothing (Loess).
#' @param is_rollmedian Apply rolling median.
#' @param rollmed Window size for rolling median.
#' @param is_histogram If TRUE, create a histogram.
#' @param threshold Optional threshold line.
#' @param hwidth Histogram width.
#' @param hheight Histogram height.
#' @param width Plot width.
#' @param height Plot height.
#' @param dpi Resolution.
#' @param device Output format (e.g., "tiff", "png").
#' @param plots_dir Directory to save plots.
#' @param nn_prop Smoothing parameter.
#' @export
generate_vcfplots <- function(data, prefix, column, y_title, plot_title, file_suffix, ylim = NULL, is_smooth = FALSE, is_rollmedian = TRUE, bwidth = 1000000,
                              is_histogram=FALSE, threshold = NULL, rollmedian = 501, hwidth, hheight, width, height, dpi, device, plots_dir, nn_prop){
  # Ensure plots directory exists
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  # Helper mini-function : Captalize the first letter 
  capitalize_first <- function(text) {
    substr(text, 1, 1) <- toupper(substr(text, 1, 1)) 
    return(text)
  }
  
  # Capitalize inbred name from prefix
  inbred <- capitalize_first(prefix)
  
  # Define color panel for chromosomes
  num_chromosomes <- length(unique(data$CHROM))
  color_panel <- rep(c("blue", "red"), length.out = num_chromosomes)
  
  # Generate scatter or smooth plot or histogram
  if(!is_histogram){
    plot <- ggplot(data, aes(x = POS, y = !!sym(column), color = CHROM)) +
      facet_grid(. ~ CHROM, scales = "free_x", space = "free_y") +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
      scale_x_continuous(breaks = seq(min(data$POS), max(data$POS), by = 20000000), labels = scales::comma) +
      scale_color_manual(values = color_panel) + guides(color = FALSE) +
      labs(title = paste0("Aligned to ", inbred, " : ", " ", plot_title), x = "\n Chromosome Position (bp) \n", y = paste0("\n", y_title, "\n")) + bsa_theme()
    
    # Add geom_point only if is_smooth is FALSE
    if (is_smooth) {
      plot <- plot + stat_smooth(method = "locfit", formula = y ~ lp(x, nn = nn_prop))
    } else if (is_rollmedian) {
      plot <- plot + geom_point(size = 0.5) + 
        geom_line(aes(y = rollmedian(!!sym(column), rollmedian, na.pad = TRUE)), color = "black", linewidth = 2.5)
    } else {
      plot <- plot + geom_point(size = 2)
    }
    
    # Add horizontal threshold line if applicable
    if (!is.null(threshold)) {plot <- plot + geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 2.5)}
    
  } else {
    plot <- ggplot(data, aes(x=POS, fill=CHROM)) + geom_histogram(binwidth=bwidth,  alpha = 0.7) +
      facet_wrap(~ CHROM, ncol = 5, scales = "free_x") + scale_fill_manual(values = color_panel) +
      theme(axis.ticks.x  = element_line(color = "black", linewidth = 1),
            axis.minor.ticks.length = unit(0.5, "cm"), axis.text.x = element_blank()) + guides(fill = FALSE) +
      labs(title = paste0("Aligned to ", inbred, " : ", " ", plot_title), x = "\n Chromosome Position (bp) \n", y = paste0("\n", y_title, "\n")) + bsa_theme()
  }
  
  # Apply fixed y-axis limits only if ylim is not NULL
  if (!is.null(ylim)) {plot <- plot + coord_cartesian(ylim = ylim)}
  
  # Construct file path
  file_path <- file.path(plots_dir, paste0(inbred, "_", file_suffix, ".", device))
  
  # Save plot
  ggsave(filename = file_path, plot = plot, device = device, 
         width = ifelse(is_histogram, hwidth, width), 
         height = ifelse(is_histogram, hheight, height), dpi = dpi)
  
  message(paste0("Plot saved to: ", file_path))
}
#' Generate Allele Frequency (AF) Plot Configurations
#' @description Generates configurations for allele frequency (AF) plots used in Bulk Segregant Analysis (BSA).
#' @param wt_mt A data frame containing shared SNPs between wild-type and mutant.
#' @param ant_wt A data frame containing SNPs unique to wild-type.
#' @param ant_mt A data frame containing SNPs unique to mutant.
#' @param ant_wt_ems A data frame containing EMS-induced SNPs unique to wild-type.
#' @param ant_mt_ems A data frame containing EMS-induced SNPs unique to mutant.
#' @param wt Wild-type identifier.
#' @param mt Mutant identifier.
#' @param rollmedian Integer. Window size for rolling median smoothing.
#' @param is_histogram Logical. If TRUE, enables histogram-based plots.
#' @param threshold Optional. Threshold value for significance filtering.
#' @param ylim Optional. Limits for the y-axis.
#' @return A list of plot configurations.
#' @export
af_configs <- function(wt_mt, ant_wt, ant_mt, ant_wt_ems, ant_mt_ems, wt, mt, rollmedian = 501, is_histogram = FALSE, ylim = NULL) {
  
  # Define common plot parameters
  parameters <- list(
    list(name = "wt_AF", data = wt_mt[wt_AF < 1], title = paste0(wt, " SNP-Index (Shared SNPs with ", mt, ")"), y_title = paste0(wt, "  AF"), plotid = paste0(wt, "_AF_shared_", mt)),
    list(name = "mt_AF", data = wt_mt[mt_AF < 1], title = paste0(mt, " SNP-Index (Shared SNPs with ", wt, ")"), y_title = paste0(mt, "  AF"), plotid = paste0(mt, "_AF_shared_", wt)),
    list(name = "AF_diff", data = wt_mt[AF_diff < 1], title = paste0("AF diff (Δ SNP-Index): ", wt, " - ", mt), y_title = "AF-diff", plotid = paste0(wt, "_AF_diff_", mt)),
    list(name = "wt_AF", data = ant_wt[wt_AF < 1], title = paste0(wt, " SNP-Index (Unique SNPs)"), y_title = paste0(wt, "  AF"), plotid = paste0(wt, "_AF_unique")),
    list(name = "mt_AF", data = ant_mt[mt_AF < 1], title = paste0(mt, " SNP-Index (Unique SNPs)"), y_title = paste0(mt, "  AF"), plotid = paste0(mt, "_AF_unique")),
    list(name = "wt_AF", data = ant_wt_ems[wt_AF < 1], title = paste0(wt, " SNP-Index (Unique ems SNPs)"), y_title = paste0(wt, "  AF"), plotid = paste0(wt, "_AF_unique_ems")),
    list(name = "mt_AF", data = ant_mt_ems[mt_AF < 1], title = paste0(mt, " SNP-Index (Unique ems SNPs)"), y_title = paste0(mt, "  AF"), plotid = paste0(mt, "_AF_unique_ems"))
  )
  
  # Generate configurations for plots
  allparameters <- list()
  for (pmt in parameters ) {
    if (is.null(pmt$data) || nrow(pmt$data) == 0) next
    
    allparameters <- append(allparameters, list(
      list(data = pmt$data, column = pmt$name, plot_title = pmt$title, file_suffix = paste0("af_locfit_", pmt$plotid), y_title = pmt$y_title,
           is_smooth = TRUE, rollmedian = 0, ylim = ylim, is_rollmedian = FALSE, is_histogram = is_histogram, threshold = NULL)))

    if (rollmedian > 0){
      allparameters <- append(allparameters, list(
        list(data = pmt$data, column = pmt$name, plot_title = pmt$title, file_suffix = paste0("af_rollmedian_", pmt$plotid),y_title = pmt$y_title,
             is_smooth = FALSE, ylim = ylim, rollmedian = rollmedian, is_rollmedian = TRUE, is_histogram = is_histogram, threshold = NULL)))
    }
  }
  return(allparameters)
}
#' Generate G-Statistic Plot Configurations
#' @description This function generates configuration settings for plotting 
#' @G-statistic values from SNP data in Bulk Segregant Analysis (BSA).
#' @param wt_mt A data frame containing SNP data for wild-type and mutant comparisons.
#' @param wt A character string representing the wild-type genotype.
#' @param mt A character string representing the mutant genotype.
#' @param rollmedian Integer. Window size for rolling median smoothing.
#' @param is_histogram Logical. If TRUE, generates histogram plots.
#' @param ylim Optional. A numeric vector specifying the y-axis limits.
#' @return A list of plot configurations for G-statistics.
#' @export
gstat_configs <- function(wt_mt, wt, mt, rollmedian = 501, is_histogram = FALSE, ylim = NULL) {
  
  # Define common plot parameters
  parameters <- list(
    list(name = "Gstat", data = wt_mt[Gstat > 0.01], title = paste0(wt, " & ", mt, " SmoothG"), y_title = "SmoothG", plotid = paste0("gstat_", wt, "_", mt))
  )
  
  # Generate configurations for plots
  allparameters <- list()
  for (pmt in parameters) {
    if (is.null(pmt$data) || nrow(pmt$data) == 0) next
    
    allparameters <- append(allparameters, list(
      list(data = pmt$data, column = pmt$name, plot_title = pmt$title, file_suffix = paste0("locfit_", pmt$plotid), y_title = pmt$y_title,
           is_smooth = TRUE, rollmedian = 0, ylim = ylim, is_rollmedian = FALSE, is_histogram = is_histogram, threshold = NULL)))
    
    
    if (rollmedian > 0) {
      allparameters <- append(allparameters, list(
        list(data = pmt$data, column = pmt$name, plot_title = pmt$title, file_suffix = paste0("rollmedian_", pmt$plotid),y_title = pmt$y_title,
             is_smooth = FALSE, ylim = ylim, rollmedian = rollmedian, is_rollmedian = TRUE, is_histogram = is_histogram, threshold = NULL)))
    }
  }
  return(allparameters)
}
#' Generate Euclidean Distance (ED) Plot Configurations
#' @description Configures plots for Euclidean Distance (ED) and ED^4 from SNP-index data.
#' @param wt_mt Data frame containing SNP data for wild-type and mutant.
#' @param wt Wild-type identifier.
#' @param mt Mutant identifier.
#' @param rollmedian Integer. Window size for rolling median smoothing.
#' @param is_histogram Logical. If TRUE, generate a histogram.
#' @param threshold Numeric. Optional threshold line.
#' @param ylim Numeric vector. Y-axis limits.
#' @return A list of plot configurations.
#' @export
ed_configs <- function(wt_mt, wt, mt, rollmedian = 501, is_histogram = FALSE, ylim = NULL) {
  
  # Define common plot parameters
  parameters <- list(
    list(name = "ED", data = wt_mt[ED > 0.01], title = paste0(wt, " & ", mt, " ED"), y_title = "ED", plotid = paste0("ED_", wt, "_", mt)),
    list(name = "ED4", data = wt_mt[ED4 > 0.01], title = paste0(wt, " & ", mt, " ED4"), y_title = "ED4", plotid = paste0("ED4_", wt, "_", mt))
  )
  
  # Generate configurations for plots
  allparameters <- list()
  for (pmt in parameters) {
    if (is.null(pmt$data) || nrow(pmt$data) == 0) next
    
    allparameters <- append(allparameters, list(
      list(data = pmt$data, column = pmt$name, plot_title = pmt$title, file_suffix = paste0("locfit_", pmt$plotid), y_title = pmt$y_title,
           is_smooth = TRUE, rollmedian = 0, ylim = ylim, is_rollmedian = FALSE, is_histogram = is_histogram, threshold = NULL)))
    
    if (rollmedian > 0) {
      allparameters <- append(allparameters, list(
        list(data = pmt$data, column = pmt$name, plot_title = pmt$title, file_suffix = paste0("rollmedian_", pmt$plotid),y_title = pmt$y_title,
             is_smooth = FALSE, ylim = ylim, rollmedian = rollmedian, is_rollmedian = TRUE, is_histogram = is_histogram, threshold = NULL)))
    }
  }
  return(allparameters)
}
#' Generate P-Value Plot Configurations
#' @description Generates configurations for p-value plots in Bulk Segregant Analysis (BSA).
#' @param wt_mt A data frame containing SNP-index and statistical data for wild-type vs. mutant.
#' @param wt Wild-type identifier.
#' @param mt Mutant identifier.
#' @param rollmedian Integer. Window size for rolling median smoothing.
#' @param is_histogram Logical. If TRUE, enables histogram-based plots.
#' @param threshold Numeric. Significance threshold (-log10 transformed p-value).
#' @param ylim Optional. Limits for the y-axis.
#' @return A list of plot configurations.
#' @export
pval_configs <- function(wt_mt, wt, mt, rollmedian = 0, threshold = -log10(0.05) * 10, ylim = NULL) {
  
  # Define common plot parameters
  parameters <- list(
    list(name = "log_adj.pval", data = wt_mt[AF_diff < 1], title = paste0("-log10(p-value) | Fisher's Exact Test : " , wt, " & ", mt),
         y_title = "-log10(p-value)", plotid = paste0("fisher_log_adj.pval_", wt, "_", mt)),
    list(name = "chi_log_adj.pval", data = wt_mt[AF_diff < 1], title = paste0("-log10(p-value) | Chisquare Test: " , wt, " & ", mt),
         y_title = "-log10(p-value)", plotid = paste0("chi_log_adj.pval_", wt, "_", mt))
  )
  
  # Generate configurations for plots
  allparameters <- list()
  for (pmt in parameters) {
    if (is.null(pmt$data) || nrow(pmt$data) == 0) next
    
    allparameters <- append(allparameters, list(
      list(data = pmt$data, column = pmt$name, plot_title = pmt$title, file_suffix = paste0("standard_", pmt$plotid), y_title = pmt$y_title,
           is_smooth = FALSE, rollmedian = 0, ylim = ylim, is_rollmedian = FALSE, is_histogram = FALSE, threshold = threshold)))
  }
  return(allparameters)
}
#' Generate Histogram Plot Configurations for SNP Analysis
#' @description Generates histogram plot configurations for visualizing homozygous SNPs per megabase (Mb).
#' @param ant_wt Data frame containing SNPs unique to the wild-type.
#' @param ant_mt Data frame containing SNPs unique to the mutant.
#' @param ant_wt_ems Data frame containing EMS-induced SNPs unique to the wild-type.
#' @param ant_mt_ems Data frame containing EMS-induced SNPs unique to the mutant.
#' @param wt Wild-type identifier.
#' @param mt Mutant identifier.
#' @param is_histogram Logical. If TRUE, enables histogram-based plots.
#' @param ylim Optional. Limits for the y-axis.
#' @return A list of histogram plot configurations.
#' @export
hist_configs <- function(ant_wt, ant_mt, ant_wt_ems, ant_mt_ems, wt, mt, is_histogram = TRUE, ylim = NULL) {
  
  # Define common plot parameters
  parameters <- list(
    list(name = "wt_AF", data = ant_wt[wt_AF>=0.99], title = paste0(wt, " unique snps only"), y_title = paste0("Homozygous SNPs per Mb"), plotid = paste0(wt, "_AF_unique")),
    list(name = "mt_AF", data = ant_mt[mt_AF>=0.99], title = paste0(mt, " unique snps only"), y_title = paste0("Homozygous SNPs per Mb"), plotid = paste0(mt, "_AF_unique")),
    list(name = "wt_AF", data = ant_wt_ems[wt_AF>=0.99], title = paste0(wt, " unique ems snps only"), y_title = paste0("Homozygous SNPs per Mb"), plotid = paste0(wt, "_AF_unique_ems")),
    list(name = "mt_AF", data = ant_mt_ems[mt_AF>=0.99], title = paste0(mt, " unique ems snps only"), y_title = paste0("Homozygous SNPs per Mb"), plotid = paste0(mt, "_AF_unique_ems"))
  )
  # Generate configurations for plots
  allparameters <- list()
  for (pmt in parameters) {
    if (is.null(pmt$data) || nrow(pmt$data) == 0) next
    
    allparameters <- append(allparameters, list(
      list(data = pmt$data, column = pmt$name, plot_title = pmt$title, file_suffix = paste0("histogram_", pmt$plotid), y_title = pmt$y_title,
           is_smooth = FALSE, rollmedian = 0, ylim = ylim, is_rollmedian = FALSE, is_histogram = TRUE, threshold = NULL)))
  }
  return(allparameters)
}
#' Automate VCF Plot Generation
#' @title Automate VCF Data Plot Generation
#' @description Generates various SNP-based plots (AF, P-value, G-statistic, ED, Histogram) for BSA analysis.
#' @param wt_mt Data frame containing SNP-index and statistical data for wild-type vs. mutant.
#' @param ant_wt Data frame containing SNPs unique to the wild-type.
#' @param ant_mt Data frame containing SNPs unique to the mutant.
#' @param ant_wt_ems Data frame containing EMS-induced SNPs unique to the wild-type.
#' @param ant_mt_ems Data frame containing EMS-induced SNPs unique to the mutant.
#' @param wt Wild-type identifier.
#' @param mt Mutant identifier.
#' @param rollmedian Integer. Window size for rolling median smoothing.
#' @param threshold Numeric. Significance threshold (-log10 transformed p-value).
#' @param ylim Optional. Limits for the y-axis.
#' @param plots_dir Character. Directory to save plots.
#' @param prefix Character. Prefix for output filenames.
#' @param device Character. Output file format (e.g., "tiff", "png").
#' @param width, height, hwidth, hheight, dpi Numeric. Plot dimensions and resolution.
#' @param nn_prop Numeric. Smoothing parameter.
#' @param plot_types Character vector. List of plots to generate (options: "af", "pval", "gstat", "ed", "histogram").
#' @export
automate_vcf_plots <- function(wt_mt, ant_wt, ant_mt, ant_wt_ems, ant_mt_ems, wt, mt, rollmedian = 501, threshold = -log10(0.05) * 10, ylim = NULL, 
                               plots_dir = "plots", prefix = "b73", device = "tiff", width = 26, height = 8, hwidth = 30, hheight = 16, dpi = 300, nn_prop=0.15,
                               plot_types = c("af","pval","gstat","ed", "histogram")) {

  all_plot_data <- list()
  
  # Normalize input for case insensitivity
  plot_types <- tolower(plot_types)
  
  # Define valid plot types
  valid_plot_types <- c("af", "pval", "gstat", "ed", "histogram")
  
  # Check for invalid plot types
  if (!all(plot_types %in% valid_plot_types)) {
    invalid_types <- plot_types[!plot_types %in% valid_plot_types]
    warning(paste("Invalid plot types detected:", paste(invalid_types, collapse = ", "), "- Skipping these."))
    plot_types <- plot_types[plot_types %in% valid_plot_types]
  }
  
  # Generate AF plots if enabled
  if ("af" %in% plot_types) {
    message("Generating AF configurations...")
    af_plot_data <- af_configs(wt_mt, ant_wt, ant_mt, ant_wt_ems, ant_mt_ems, wt, mt, rollmedian)
    all_plot_data <- c(all_plot_data, af_plot_data)
  }

  # Generate Gstat plots if enabled
  if ("gstat" %in% plot_types) {
    message("Generating Gstat configurations...")
    gstat_plot_data <- gstat_configs(wt_mt, wt, mt, rollmedian)
    all_plot_data <- c(all_plot_data, gstat_plot_data)
  }

  # Generate ED plots if enabled
  if ("ed" %in% plot_types) {
    message("Generating ED configurations...")
    ed_plot_data <- ed_configs(wt_mt, wt, mt, rollmedian)
    all_plot_data <- c(all_plot_data, ed_plot_data)
  }

  # Generate pvalue plots if enabled
  if ("pval" %in% plot_types) {
    message("Generating pvalue configurations...")
    pval_plot_data <- pval_configs(wt_mt, wt, mt, threshold)
    all_plot_data <- c(all_plot_data, pval_plot_data)
  }

  # Generate pvalue plots if enabled
  if ("histogram" %in% plot_types) {
    message("Generating histogram configurations...")
    hist_plot_data <- hist_configs(ant_wt, ant_mt, ant_wt_ems, ant_mt_ems, wt, mt, is_histogram = TRUE, ylim = NULL)
    all_plot_data <- c(all_plot_data, hist_plot_data)
  }

  # check if there's data to plot
  if (length(all_plot_data) == 0) {
    message("No plots were selected. Please enable at least one plot type.")
    return()
  }

  for (plots in all_plot_data) {
  
    message(paste0("Generating plot: ", plots$plot_title))
    generate_vcfplots(data = plots$data, column = plots$column, y_title = plots$y_title, plot_title = plots$plot_title, file_suffix = plots$file_suffix, ylim = plots$ylim, 
                      is_smooth = plots$is_smooth, is_rollmedian = plots$is_rollmedian,is_histogram = plots$is_histogram,
                      threshold = plots$threshold, nn_prop=nn_prop, width=width, height=height, hwidth=hwidth, hheight=hheight, dpi=dpi,  device=device, 
                      prefix = prefix, plots_dir = plots_dir)
  }
  message("All selected plots successfully generated and saved!")
}     
#' @title Process and Visualize VCF Data
#' @description Loads, processes, and visualizes SNP data from VCF files.
#' @param vcf_dir Character. Directory containing VCF files.
#' @param prefix Character. Prefix used for file naming.
#' @param pattern Character. File pattern to match (default: "*snps.txt").
#' @param output_dir Character. Directory for analysis results.
#' @param plots_dir Character. Directory for saving plots.
#' @param save_results Logical. Whether to save SNP analysis results.
#' @param GenotypeS Named vector. Defines wild-type (`wt`) and mutant (`mt`) identifiers.
#' @param min_DP Integer. Minimum depth filter for SNPs.
#' @param min_QUAL Integer. Minimum quality score filter for SNPs.
#' @param threshold Numeric. Threshold for statistical significance (-log10 transformed p-value).
#' @param rollmedian Integer. Window size for rolling median smoothing.
#' @param n_prop Numeric. Smoothing parameter.
#' @param hwidth, hheight, width, height, dpi Numeric. Plot dimensions and resolution.
#' @param device Character. Output format (e.g., "tiff", "png").
#' @param plot_types Character vector. List of plot types to generate (options: "af", "pval", "gstat", "histogram").
#' @return A list containing SNP analysis results and generated plots.
#' @export
visualize_vcfdata <- function(vcf_dir, prefix, pattern, output_dir, plots_dir, 
                              save_results = FALSE, Genotypes = c(wt = "wt", mt = "mt"), min_DP, min_QUAL, plot_data = TRUE,
                              threshold = -log10(0.05) * 10, n_prop = 0.1, rollmedian = 501, hwidth = 30, hheight = 16, ylim = NULL,
                              width = 26, height = 8, dpi = 300, device = "tiff", plot_types = c("af", "pval","ed", "gstat", "histogram")) {
  
  message("Starting VCF Data Processing Pipeline...")
  
  # Validate input directories
  if (!dir.exists(vcf_dir)) stop("Error: Specified VCF directory does not exist.")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  # Validate Genotype identifiers
  if (!all(c("wt", "mt") %in% names(Genotypes))) {
    stop("Error: GenotypeS must be a named vector with 'wt' (wild-type) and 'mt' (mutant).")
  }
  
  
  # Step 1: Load and Process VCF Data
  message("Loading VCF data...")
  geno_data <- import_vcfdata(vcf_dir, prefix, pattern, Genotypes, min_DP, min_QUAL)
  
  # Step 2: Merge and Analyze SNPs
  message("Merging and Analyzing SNPs...")
  results <- merge_analyze_vcfdata(geno_data, prefix, save_results, output_dir)

  if (plot_data){
  # Extract SNP datasets from results
  wt_mt <- results$wt_mt
  ant_wt <- results$ant_wt
  ant_mt <- results$ant_mt
  ant_wt_ems <- results$ant_wt_ems
  ant_mt_ems <- results$ant_mt_ems
    
  automate_vcf_plots(
    wt_mt, ant_wt, ant_mt, ant_wt_ems, ant_mt_ems, wt = Genotypes[["wt"]], mt = Genotypes[["mt"]], rollmedian = rollmedian, threshold = threshold, ylim = ylim, 
    plots_dir = plots_dir, prefix = prefix, device = device, width = width, height = height, hwidth = hwidth, hheight = hheight, dpi = dpi, 
    nn_prop = n_prop, plot_types = plot_types
    )
  }
  message("CF Data Processing and Visualization Completed!")
  return(results)
}
#' Process VCF Data and Generate Visualizations
#' @description This function processes VCF files for a list of wild-type and mutant genotypes across multiple prefixes (inbreds).
#' @description It generates allele frequency plots, p-value distributions, and other visualizations.
#' @param vcf_dir Directory where the VCF files are stored.
#' @param output_dir Directory to save processed data.
#' @param plots_dir Directory to store plots.
#' @param wt_list Vector of wild-type genotypes.
#' @param mt_list Vector of mutant genotypes.
#' @param prefix_list Vector of prefixes.
#' @param pattern Regex pattern to match VCF files.
#' @param save_results Logical. If TRUE, saves processed results.
#' @param min_DP Minimum depth for filtering.
#' @param min_QUAL Minimum quality score.
#' @param threshold Threshold for highlighting significant values.
#' @param rollmedian Rolling median window size.
#' @param hheight height dimensions for histogram
#' @param hwidth width dimensions for histogram
#' @param height height dimensions for other graphs 
#' @param width width dimensions for other graphs
#' @param device Output format (e.g., "tiff", "png").
#' @param plot_types A vector specifying which plots to generate.
#' @return A list of processed results.
#' @export

auto_visualize_vcfdata <- function(vcf_dir, output_dir, plots_dir, wt_list, mt_list, prefix_list, pattern, save_results = FALSE, plot_data = TRUE,
                                   min_DP = 10, min_QUAL = 30, threshold = -log10(0.05) * 10, n_prop = 0.1, rollmedian = 501, ylim = NULL,
                                   hwidth = 30, hheight = 18, width = 40, height = 13, dpi = 300, device = "tiff", plot_types = c("af", "pval", "gstat", "ed","histogram")){
  all_results <- list()
  for (pfix in prefix_list) {
    print(paste("Processing inbred:", pfix))
    
    # if single wiltype and mutant case, process as pair
    if (length(wt_list) == 1 && length(mt_list) == 1) {
      wt <- wt_list[1]
      mt <- mt_list[1]
      Genotypes <- c(wt = wt, mt = mt)
      
      # Create a sub-directory for the WT-MT combination for the plots
      sub_plots <- file.path(plots_dir, paste0(pfix, "_", wt, "_", mt))
      if (!dir.exists(sub_plots)) dir.create(sub_plots, recursive = TRUE)
      
      message(paste("Processing single wildtype:", wt, "and single mutant:", mt, "in", sub_plots))
      result <- paste0(pfix, "_", wt, "_", mt)
      print(result)
      
      # Create sub-directory for each WT-MT pair
      result_data <- visualize_vcfdata(vcf_dir = vcf_dir, prefix = pfix, pattern = pattern, output_dir = output_dir, plots_dir = sub_plots, save_results = save_results, 
                                       Genotypes = Genotypes, min_DP = min_DP, min_QUAL = min_QUAL, threshold = threshold, n_prop = n_prop, rollmedian = rollmedian, plot_data = plot_data,
                                       ylim = ylim, hwidth = hwidth, hheight = hheight, width = width, height = height, dpi = dpi, device = device, plot_types = plot_types)
      
      # Store the result in the list
      if (!is.null(result_data)) {
        all_results[[result]] <- result_data
        print(paste("Stored result for:", result))
      } else {warning(paste("Failed to process wildtype:", wt, "and mutant:", mt))}
      message(paste("Completed processing for wildtype:", wt, "and mutant:", mt, "for", pfix, "in", sub_plots))
      
    } else{
      # Combine and ensure no duplicates
      all_list <- unique(c(wt_list, mt_list))
      if (length(all_list) < 2) {
        warning("No enough unique elements in in wt_list and mt_list to create pairs.")
        next}
      
      unique_pairs <- combn(all_list, 2, simplify = FALSE)
      
      for (pair in unique_pairs) {
        wt <- pair[1]
        mt <- pair[2]
        Genotypes <- c(wt = wt, mt = mt)
        
        # Create sub-directory for each WT-MT pair
        sub_plots <- file.path(plots_dir, paste0(pfix, "_", wt, "_", mt))
        if (!dir.exists(sub_plots)) dir.create(sub_plots, recursive = TRUE)
        
        message(paste("Processing wildtype:", wt, "and mutant:", mt, "in", sub_plots))
        
        result <- paste0(pfix,"_", wt, "_", mt)
        print(result)
        
        # Call visualize_vcfdata function and store the result
        result_data <- visualize_vcfdata(vcf_dir = vcf_dir, prefix = pfix, pattern = pattern, output_dir = output_dir, plots_dir = sub_plots, save_results = save_results, 
                                         Genotypes = Genotypes, min_DP = min_DP, min_QUAL = min_QUAL, threshold = threshold, n_prop = n_prop, rollmedian = rollmedian, plot_data = plot_data,
                                         ylim = ylim, hwidth = hwidth, hheight = hheight, width = width, height = height, dpi = dpi, device = device, plot_types = plot_types)
        
        # Store the result in the list
        if (!is.null(result_data)) {
          all_results[[result]] <- result_data
          print(paste("Stored result for:", result))
        } else {
          warning(paste("Failed to process wildtype:", wt, "and mutant:", mt))
        }
      }
    }
  }
  print(paste("Total stored results:", length(all_results)))
  return(all_results)
}


#vcf_dir = "/Users/zebosi/Documents/BSA/S7_6508K/data/"
#output_dir = "/Users/zebosi/Documents/BSA/S7_6508K/postanalysis"
#plots_dir = "/Users/zebosi/Documents/BSA/S7_6508K/plots"
#wt_list <- c("S7A6508K","S7B6508K")
#mt_list <- c("S7B6508K", "S7C6508K")
#prefix_list = c("b73")
#pattern = "snps\\.txt$"
#threshold = -log10(0.05) * 10
#rollmedian = 501




#results_data <- auto_visualize_vcfdata(plot_data = FALSE,
#                                       vcf_dir = vcf_dir, output_dir = output_dir, plots_dir = plots_dir,
#                                       wt_list = wt_list, mt_list = mt_list, prefix_list = prefix_list,
#                                       pattern = pattern, threshold = threshold, rollmedian = rollmedian
#)

