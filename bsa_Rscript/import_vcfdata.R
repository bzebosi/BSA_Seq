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
