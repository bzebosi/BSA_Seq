# Find location of sourcecodes 
source(file = "~/Documents/BSA_Ts3/customized_Rscripts/multi_package_installer.R")
source(file = "~/Documents/BSA_Ts3/customized_Rscripts/bsa_theme.R")
# function to load and visualize data
visualize_vcfdata <- function(vcf_dir, inbred, prefix, nn_prop = 0.1, wt, mt,
                        plots_dir = "plots", width = 26, height = 8, dpi = 300){
  # Generate list of VCF files in the directory
  vcf_list <- list.files(path = vcf_dir, pattern = "*.txt", full.names = TRUE)
  
  # Find files that match the prefix and contain "mt" (mutant) or "wt" (wild-type)
  mt_file <- vcf_list[grepl(paste0("^", prefix, "_", mt), basename(vcf_list))]
  wt_file <- vcf_list[grepl(paste0("^", prefix, "_", wt), basename(vcf_list))]
  
  # Check if both files exist
  if (length(mt_file) > 0 && length(wt_file) > 0) {
    
    # Load mutant and wild-type data
    mt_data <- read.delim(mt_file, header = TRUE)
    wt_data <- read.delim(wt_file, header = TRUE)
    
    # Define expected column names
    col_names <- c("CHROM","POS","REF","ALT","QUAL","DP","Fref","Rref","Falt","Ralt")
    colnames(mt_data) <- col_names
    colnames(wt_data) <- col_names
    
    # Convert to data.table if necessary
    if (!is.data.table(wt_data)) wt_data <- as.data.table(wt_data)
    if (!is.data.table(mt_data)) mt_data <- as.data.table(mt_data)
    
    # Filter based on DP and QUAL for both wild-type and mutant data
    wt_data <- wt_data[!is.na(DP) & DP > 20 & QUAL >= 200]
    mt_data <- mt_data[!is.na(DP) & DP > 20 & QUAL >= 200]
  }
  if (length(wt_data) > 0 && length(mt_data) > 0) {
    setnames(wt_data, old = setdiff(names(wt_data), c("CHROM", "POS")), 
             new = paste0("wt_", setdiff(names(wt_data), c("CHROM", "POS"))))
    setnames(mt_data, old = setdiff(names(mt_data), c("CHROM", "POS")), 
             new = paste0("mt_", setdiff(names(mt_data), c("CHROM", "POS"))))
    
    # Merge the two filtered datasets by CHROM and POS
    wt_mt <- merge(wt_data, mt_data, by = c("CHROM", "POS"), suffixes = c("wt_", "mt_"))
    wt_mt[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
    wt_mt <- wt_mt[order(CHROM, POS)]
    
    # anti join keep 
    ant_mt <- anti_join(mt_data, wt_data, by = c("CHROM", "POS"))
    ant_mt[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
    ant_mt <- ant_mt[order(CHROM, POS)]
    
    # anti join keep 
    ant_wt <- anti_join(wt_data, mt_data, by = c("CHROM", "POS"))
    ant_wt[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
    ant_wt <- ant_wt[order(CHROM, POS)]
    
    # Calculate allele frequencies for wild-type and mutant
    wt_mt[, wt_AFQ := (wt_Falt + wt_Ralt) / (wt_Fref + wt_Rref + wt_Falt + wt_Ralt)]
    wt_mt[, mt_AFQ := (mt_Falt + mt_Ralt) / (mt_Fref + mt_Rref + mt_Falt + mt_Ralt)]
    
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
    
    # Dynamically create the final columns to keep
    keepcols <- c("CHROM", "POS", "wt_REF", "wt_ALT", 
                  "mt_REF", "mt_ALT", "wt_QUAL", "mt_QUAL", "wt_DP", "mt_DP", 
                  "wt_AFQ", "mt_AFQ", "AF_diff", "ED", "ED4", "adj.p.value", "log_adj.pval")
    
    # Keep only the specified columns
    wt_mt <- wt_mt[, keepcols, with = FALSE]
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
  
  # Plot and save AD8_unAF_transformed
  ggsave(
    filename = file.path(plots_dir, paste0(inbred, "_wt_AF_untranformed.tiff")),
    plot = base_plot(wt_mt[wt_AFQ < 1], "wt_AFQ", paste0(inbred, "\n WT SNP Index \n"))+
      geom_point() +  coord_cartesian(ylim = c(0, 1)) +
      geom_line(aes(y = rollmedian(wt_AFQ, 401, na.pad = TRUE)), color = "black"),
    device = "tiff", width = width, height = height, dpi = dpi)
  message(paste0("Plot for '", inbred, "' saved to: ", file.path(plots_dir, paste0(inbred, "_wt_AF_untranformed.tiff"))))
  
  # Plot and save mt_unAF_transformed
  ggsave(
    filename = file.path(plots_dir, paste0(inbred, "_mt_AF_untranformed.tiff")),
    plot = base_plot(wt_mt[mt_AFQ < 1], "mt_AFQ", paste0(inbred, "\n Mutant SNP Index \n"))+
      geom_point() +  coord_cartesian(ylim = c(0, 1)) +
      geom_line(aes(y = rollmedian(mt_AFQ, 401, na.pad = TRUE)), color = "black"),
    device = "tiff", width = width, height = height, dpi = dpi)
  message(paste0("Plot for '", inbred, "' saved to: ", file.path(plots_dir, paste0(inbred, "_mt_AF_untranformed.tiff"))))
  
  # Plot and save AF_diff_untransformed
  ggsave(
    filename = file.path(plots_dir, paste0(inbred, "_AF_diff_untranformed.tiff")),
    plot = base_plot(wt_mt[AF_diff < 1], "AF_diff", paste0(inbred, "\n Δ SNP Index \n"))+
      geom_point() +  coord_cartesian(ylim = c(0, 1)) +
      geom_line(aes(y = rollmedian(AF_diff, 401, na.pad = TRUE)), color = "black"),
    device = "tiff", width = width, height = height, dpi = dpi)
  message(paste0("Plot for '", inbred, "' saved to: ", file.path(plots_dir, paste0(inbred, "_AF_diff_untranformed.tiff"))))
  
  # Plot and save wt_AF_transformed
  ggsave(
    filename = file.path(plots_dir, paste0(inbred, "_wt_tranformed.tiff")),
    plot = base_plot(wt_mt[wt_AFQ  < 1], "wt_AFQ", paste0(inbred, "\n WT SNP Index \n"))+
      stat_smooth(method = "locfit", formula = y ~ lp(x, nn = nn_prop))+
      coord_cartesian(ylim = c(0.5, 1)),
    device = "tiff", width = width, height = height, dpi = dpi)
  message(paste0("Plot for '", inbred, "' saved to: ", file.path(plots_dir, paste0(inbred, "_wt_tranformed.tiff"))))
  
  # Plot and save mutant_AF_transformed
  ggsave(
    filename = file.path(plots_dir, paste0(inbred, "_mt_tranformed.tiff")),
    plot = base_plot(wt_mt[mt_AFQ  < 1], "mt_AFQ", paste0(inbred, "\n Mutant SNP Index \n"))+
      stat_smooth(method = "locfit", formula = y ~ lp(x, nn = nn_prop))+
      coord_cartesian(ylim = c(0.5, 1)),
    device = "tiff", width = width, height = height, dpi = dpi)
  message(paste0("Plot for '", inbred, "' saved to: ", file.path(plots_dir, paste0(inbred, "_mt_tranformed.tiff"))))
  
  
  # Plot and save _AF_diff_transformed
  ggsave(
    filename = file.path(plots_dir, paste0(inbred, "_AF_diff_tranformed.tiff")),
    plot = base_plot(wt_mt[AF_diff  < 1], "AF_diff", paste0(inbred, "\n Δ SNP Index \n"))+
      stat_smooth(method = "locfit", formula = y ~ lp(x, nn = nn_prop))+
      coord_cartesian(ylim = c(0.000, 0.1)),
    device = "tiff", width = width, height = height, dpi = dpi)
  message(paste0("Plot for '", inbred, "' saved to: ", file.path(plots_dir, paste0(inbred, "_AF_diff_tranformed.tiff"))))
  
  # Plot and save ED transformed
  ggsave(
    filename = file.path(plots_dir, paste0(inbred, "_ED_smooth.tiff")),
    plot = base_plot(wt_mt[ED > 1], "ED", paste0(inbred, "\n ED \n"))+
      stat_smooth(method = "locfit", formula = y ~ lp(x, nn = nn_prop)),
    device = "tiff", width = width, height = height, dpi = dpi)
  message(paste0("Plot for '", inbred, "' saved to: ", file.path(plots_dir, paste0(inbred, "_ED_smooth.tiff"))))
  
  # Plot and save ED4 transformed
  ggsave(
    filename = file.path(plots_dir, paste0(inbred, "_ED4_smooth.tiff")),
    plot = base_plot(wt_mt[ED4 > 1], "ED4", paste0(inbred, "\n ED4 \n"))+
      stat_smooth(method = "locfit", formula = y ~ lp(x, nn = nn_prop)),
    device = "tiff", width = width, height = height, dpi = dpi)
  message(paste0("Plot for '", inbred, "' saved to: ", file.path(plots_dir, paste0(inbred, "_ED4_smooth.tiff"))))

  # Plot and save log_adj.pval transformed
  ggsave(
    filename = file.path(plots_dir, paste0(inbred, "_log_adj.pval.tiff")),
    plot = base_plot(wt_mt[log_adj.pval >= 0], "log_adj.pval", paste0(inbred, "\n -log10(p-value) \n"))+
      geom_point(),
    device = "tiff", width = width, height = height, dpi = dpi)
  message(paste0("Plot for '", inbred, "' saved to: ", file.path(plots_dir, paste0(inbred, "_log_adj.pval.tiff"))))
  
  
  # Plot and save AD7 only
  ggsave(
    filename = file.path(plots_dir, paste0(inbred, "_ant_mt_AF_untranformed.tiff")),
    plot = base_plot(ant_mt[mt_AF < 1], "mt_AF", paste0(inbred, "\n Mutant SNP Index \n"))+
      geom_point() +  coord_cartesian(ylim = c(0, 1)) +
      geom_line(aes(y = rollmedian(mt_AF, 501, na.pad = TRUE)), color = "black"),
    device = "tiff", width = width, height = height, dpi = dpi)
  message(paste0("Plot for '", inbred, "' saved to: ", file.path(plots_dir, paste0(inbred, "_ant_mt_AF_untranformed.tiff"))))
  
  # Plot and save AD7 only
  ggsave(
    filename = file.path(plots_dir, paste0(inbred, "_ant_mt_AF_tranformed.tiff")),
    plot = base_plot(ant_mt[mt_AF < 1], "mt_AF", paste0(inbred, "\n Mutant SNP Index \n"))+
      stat_smooth(method = "locfit", formula = y ~ lp(x, nn = nn_prop)),
    device = "tiff", width = width, height = height, dpi = dpi)
  message(paste0("Plot for '", inbred, "' saved to: ", file.path(plots_dir, paste0(inbred, "_ant_mt_AF_tranformed.tiff"))))
  
  # Plot and save AD7 only
  ggsave(
    filename = file.path(plots_dir, paste0(inbred, "_ant_wt_AF_untranformed.tiff")),
    plot = base_plot(ant_wt[wt_AF < 1], "wt_AF", paste0(inbred, "\n WT SNP Index \n"))+
      geom_point() +  coord_cartesian(ylim = c(0, 1)) +
      geom_line(aes(y = rollmedian(wt_AF, 501, na.pad = TRUE)), color = "black"),
    device = "tiff", width = width, height = height, dpi = dpi)
  message(paste0("Plot for '", inbred, "' saved to: ", file.path(plots_dir, paste0(inbred, "_ant_wt_AF_untranformed.tiff"))))
  
  # Plot and save AD7 only
  ggsave(
    filename = file.path(plots_dir, paste0(inbred, "_ant_wt_AF_tranformed.tiff")),
    plot = base_plot(ant_wt[wt_AF < 1], "wt_AF", paste0(inbred, "\n WT SNP Index \n"))+
      stat_smooth(method = "locfit", formula = y ~ lp(x, nn = nn_prop)),
    device = "tiff", width = width, height = height, dpi = dpi)
  message(paste0("Plot for '", inbred, "' saved to: ", file.path(plots_dir, paste0(inbred, "_ant_wt_AF_tranformed.tiff"))))
  
  # Return both datasets as a list
  return(list(wt_mt = wt_mt, ant_mt = ant_mt, ant_wt = ant_wt))
}

data="~/Documents/BSA_Ts3/Ts3/data"
plots="~/Documents/BSA_Ts3/Ts3/plots"
inbreds <- c("A188")
prefixes <- c("a188")
wt <- c("AD8")
mt <- c("AD7")

# Loop through each genotype and prefix, and call load_visualize_vcfdata
for (i in seq_along(inbreds)) {
  inbx <- inbreds[i]
  pfix <- prefixes[i]
  print(paste("Processing inbred:", inbx))
  assign(inbx, visualize_vcfdata(vcf_dir=data, plots_dir = plots, 
                                 inbred = inbx, prefix=pfix, wt=wt, mt=mt,
                                 nn_prop = 0.12, 
                                 width = 26, height = 8, dpi = 300)
         )
  message(paste("Completed processing for:", inbred))
}
