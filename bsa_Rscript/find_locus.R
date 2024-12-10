# Run required dependencies and available on github
source(file = "https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/multi_package_installer.R")
source(file = "https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/bsa_theme.R")

# functions load wildtyp and mutant vcfdata
visualize_vcfdata <- function(vcf_dir, inbred, prefix, Genotype=c(wt="wt", mt="mt"),
                          nn_prop = 0.1, plots_dir = "plots", width = 26, height = 8, dpi = 300) {
  # Generate a list of VCF files in the directory
  vcf_list <- list.files(path = vcf_dir, pattern = "*.txt", full.names = TRUE)
  
  # Initialize a list to hold data for each genotype
  geno_data <- list()
  
  # Loop through Genotype to find and process matching files
  for (genotype in names(Genotype)) {
    file_pattern <- paste0("^", prefix, "_", Genotype[[genotype]])
    geno_file <- vcf_list[grepl(file_pattern, basename(vcf_list))]
    
    if (length(geno_file) > 0) {
      # Load the genotype data
      data <- read.delim(geno_file, header = TRUE)
      
      # Print a summary of loaded data
      print(paste("Loaded", genotype, "data:", nrow(data), "rows"))
      
      # Define expected column names
      col_names <- c("CHROM", "POS", "REF", "ALT", "QUAL", "DP", "Fref", "Rref", "Falt", "Ralt")
      colnames(data) <- col_names
      
      # Convert to data.table if necessary
      if (!is.data.table(data)) data <- as.data.table(data)
      
      # Filter based on DP and QUAL
      data <- data[!is.na(DP) & DP > 20 & QUAL >= 50]
      
      # Calculate allele frequencies for wild-type and mutant
      data[, AF := (Falt + Ralt) / (Fref + Rref + Falt + Ralt)]
      
      setnames(data, old = setdiff(names(data), c("CHROM", "POS")), 
               new = paste0(genotype,"_", setdiff(names(data), c("CHROM", "POS"))))
      
      # Add the processed data to the list
      geno_data[[genotype]] <- data
    }
  }
  
  if (!is.null(geno_data$wt) && !is.null(geno_data$mt)){
    wt_data <- geno_data$wt
    mt_data <- geno_data$mt
    
    # Merge the two filtered datasets by CHROM and POS
    wt_mt <- merge(wt_data, mt_data, by = c("CHROM", "POS"), suffixes = c("wt_", "mt_"))
    wt_mt[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
    wt_mt <- wt_mt[order(CHROM, POS)]
    
    # anti join and keep position only unique to mt dataset
    ant_mt <- anti_join(mt_data, wt_data, by = c("CHROM", "POS"))
    ant_mt[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
    ant_mt <- ant_mt[order(CHROM, POS)]
    
    # anti join and keep position only unique to wt dataset 
    ant_wt <- anti_join(wt_data, mt_data, by = c("CHROM", "POS"))
    ant_wt[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
    ant_wt <- ant_wt[order(CHROM, POS)]
    
    # Calculate allele frequency differences
    wt_mt[, AF_diff := wt_AF - mt_AF]
    
    # Calculate Euclidean Distance (ED) and its fourth power (ED4)
    wt_mt[, ED := sqrt((mt_Fref + mt_Rref - wt_Fref - wt_Rref)^2 + (mt_Falt + mt_Ralt - wt_Falt - wt_Ralt)^2)]
    wt_mt[, ED4 := ED^4]
    
    #Calculate p-value using Fisher's exact test for each position
    wt_mt[, p.value := mapply(function(a_ref, a_alt, b_ref, b_alt) {
      fisher.test(matrix(c(a_alt, a_ref, b_alt, b_ref), nrow = 2))$p.value
    }, mt_Fref + mt_Rref, mt_Falt + mt_Ralt, wt_Fref + wt_Rref, wt_Falt + wt_Ralt)]
    
    # Adjust p-values for multiple testing
    wt_mt[, adj.p.value := p.adjust(p.value, method = "bonferroni")]
    
    # Transform adjusted p-value to -10 * log10(adj.p.value)
    wt_mt[, log_adj.pval := -10 * log10(adj.p.value)]
    
  }
  
  # Create directory if it doesn't exist
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  base_plot <- function(data, aes_y, y_label) {
    ggplot(data, aes(x = POS, y = !!sym(aes_y), color = CHROM)) +
      facet_grid(. ~ CHROM, space = "free_x", scales = "free_x") +
      
      scale_color_manual(values = rep(c("blue", "red"), 5)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.0001))+
      guides(color = FALSE) + bsa_theme() +
      labs(title = "", x = "\n Chromosome", y = y_label)
  }
  
  plot1 <- function(data, column, plot_title, file_suffix, ylim = c(0, 1), lcolor = "black", lsize = 1, psize = 0.5) {
    # Construct the file path for saving the plot
    file_path <- file.path(plots_dir, paste0(inbred, "_", file_suffix, ".tiff"))
    
    actual_line_color <- ifelse(lcolor == "NA", NA, lcolor)
    
    # Generate and save the plot
    ggsave(
      filename = file_path,
      plot = base_plot(data, column, paste0(inbred, "\n", plot_title, "\n")) +
        geom_point(size=psize) + coord_cartesian(ylim = ylim) + 
        geom_line(aes(y = rollmedian(!!sym(column), 501, na.pad = TRUE)), color = actual_line_color, size = lsize),
      device = "tiff", width = width, height = height, dpi = dpi
    )
    message(paste0("Plot for '", inbred, "' saved to: ", file_path))
  }
  
  # Generate and save plots for WT, MT, and AF_diff
  pt1 <- list(
    list(data = wt_mt[wt_AF < 1], column = "wt_AF", plot_title = "WT SNP Index", file_suffix = "wt_AF_untransformed", ylim = c(0.1, 1),lcolor = "black", lsize = 1, psize=0.5),
    list(data = wt_mt[mt_AF < 1], column = "mt_AF", plot_title = "MT SNP Index", file_suffix = "mt_AF_untransformed", ylim = c(0.1, 1),lcolor = "black", lsize = 1, psize=0.5),
    list(data = wt_mt[AF_diff < 1], column = "AF_diff", plot_title = "Δ SNP Index", file_suffix = "AF_diff_untranformed", ylim = c(-1, 1),lcolor = "black", lsize = 1, psize=0.5),
    list(data = wt_mt[log_adj.pval >= 0], column = "log_adj.pval", plot_title = "-log10(p-value)", file_suffix = "log_adj.pval", ylim = c(0, 60),lcolor = "NA", lsize = 1, psize=0.5),
    list(data = ant_mt[mt_AF < 1], column = "mt_AF", plot_title = "Mutant SNP Index", file_suffix = "ant_mt_AF_untranformed", ylim = c(0.1, 1),lcolor = "black", lsize = 1, psize=0.5),
    list(data = wt_mt[ED > 1], column = "ED", plot_title = "ED", file_suffix = "ED_untranformed", ylim = c(1, 50),lcolor = "black", lsize = 1, psize=0.5),
    list(data = wt_mt[ED4 > 1], column = "ED4", plot_title = "ED4", file_suffix = "ED4_untranformed", ylim = c(1,250000),lcolor = "black", lsize = 1, psize=0.5)
  )
  # Loop through the pt1 and call plot1
  for (plots in pt1) {
    plot1(
      data = plots$data, column = plots$column,
      plot_title = plots$plot_title, file_suffix = plots$file_suffix,
      ylim = plots$ylim, lcolor = plots$lcolor, lsize = plots$lsize, psiz = plots$psize
    )
  }
  
  
  plot2 <- function(data, column, plot_title, file_suffix, ylim = c(0, 1)) {
    # Construct the file path for saving the plot
    file_path <- file.path(plots_dir, paste0(inbred, "_", file_suffix, ".tiff"))
    
    # Generate and save the plot
    ggsave(
      filename = file_path,
      plot = base_plot(data, column, paste0(inbred, "\n", plot_title, "\n")) +
        stat_smooth(method = "locfit", formula = y ~ lp(x, nn = nn_prop))+ coord_cartesian(ylim = ylim),
      device = "tiff", width = width, height = height, dpi = dpi
    )
    message(paste0("Plot for '", inbred, "' saved to: ", file_path))
  }
  
  pt2 <- list(
    list(data=wt_mt[mt_AF < 1], column="mt_AF", plot_title ="Mutant SNP Index", file_suffix ="mt_AF_tranformed", ylim = c(0, 1)),
    list(data=wt_mt[wt_AF < 1], column="wt_AF", plot_title ="WT SNP Index", file_suffix ="wt_AF_tranformed", ylim = c(0, 1)),
    list(data=ant_mt[mt_AF < 1], column="mt_AF", plot_title ="Mutant SNP Index", file_suffix ="ant_mt_AF_tranformed", ylim = c(0, 1)),
    list(data=ant_wt[wt_AF < 1], column="wt_AF", plot_title ="WT SNP Index", file_suffix ="ant_wt_AF_tranformed", ylim = c(0, 1)),
    list(data=wt_mt[AF_diff < 1], column="AF_diff", plot_title ="Δ SNP Index", file_suffix ="AF_diff_tranformed", ylim = c(-1, 1)),
    list(data=wt_mt[ED > 1], column="ED", plot_title ="ED", file_suffix ="ED_tranformed", ylim = c(0.0, 50)),
    list(data=wt_mt[ED4 > 1], column="ED4", plot_title ="ED4", file_suffix ="ED4_tranformed", ylim = c(0.0, 250000))
  )
  for (graphs in pt2) {
    plot2(
      data = graphs$data, column = graphs$column,
      plot_title = graphs$plot_title, file_suffix = graphs$file_suffix, ylim = graphs$ylim
    )
  }
  
  return(list(wt_mt = wt_mt, ant_mt = ant_mt, ant_wt = ant_wt))
}




