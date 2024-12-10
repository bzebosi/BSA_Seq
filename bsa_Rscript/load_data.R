# Run required dependencies and available on github
source(file = "https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/multi_package_installer.R")
source(file = "https://raw.githubusercontent.com/bzebosi/BSA_Seq/main/bsa_Rscript/bsa_theme.R")

# functions load wildtyp and mutant vcfdata
loadvcfdata <- function(vcf_dir, inbred, prefix, Genotype=c(wt="wt", mt="mt")) {
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
  
  return(list(wt_mt = wt_mt, ant_mt = ant_mt, ant_wt = ant_wt))
}
