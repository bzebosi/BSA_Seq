load_vcfdata <- function (vcf_dir, prefix, pattern, Genotype, min_DP, min_QUAL){
  # This function processes VCF files for specific genotypes, filters the data based on depth (DP) and quality score (QUAL),
  # and computes allele frequencies (AF). It returns a list containing processed data for each genotype.
  
  # Description of Parameters:
  # vcf_dir: Directory with the VCF files. Must be a valid path (e.g., "path/to/data").
  # prefix: Common prefix for the VCF files.(e.g., "b73"). Used to identify files for specific genotypes.
  # pattern: File name pattern to identify VCF files (e.g., "snps\\.txt$"). 
  # Genotype: list where the names are the genotypes of interest (e.g., "wt", "mt") and the values are file suffixes.
  #           corresponding to those genotypes (e.g., list(wt = "MB1", mt = "MB2")).
  # min_DP: Minimum DP required to include a variant in the analyss and Numeric value (e.g., 10).
  # min_QUAL: Minimum QUAL required to include a variant in the analysis. Numeric value (e.g., 30).
  
  # List VCF files
  vcf_list <- list.files(path = vcf_dir, pattern = pattern, full.names = TRUE)
  if (length(vcf_list) == 0) {stop("No vcf files found.")}
  
  # initialize a list to store genotype data
  geno_data <- list()
  
  # Process each genotype
  for (genotype in names(Genotype)) {
    # Generate file pattern for the genotype
    file_pattern <- paste0("^", prefix, "_", Genotype[[genotype]])
    geno_file <- vcf_list[grepl(file_pattern, basename(vcf_list))]
    
    if (length(geno_file) > 0) {
      data <- tryCatch({fread(geno_file)}, error = function(e) {
        stop(paste("Error reading file for", genotype, ":", geno_file, ":", e$message))})
      
      # check data and ensure the file is not emty or missing
      if (is.null(data) || nrow(data) == 0) {stop(paste("File for", genotype, "is either empty or missing"))}
      message(paste("Loaded", genotype, "data with", nrow(data), "rows."))
      
      # Rename columns
      col_names <- c("CHROM", "POS", "REF", "ALT", "QUAL", "DP", "Fref", "Rref", "Falt", "Ralt")
      colnames(data) <- col_names
      
      # Convert to data.table if it is not already
      if (!is.data.table(data)) data <- as.data.table(data)
      
      # Filter rows based on depth (DP) and quality (QUAL)
      data <- data[!is.na(DP) & DP > min_DP & QUAL >= min_QUAL]
      
      # Calculate allele frequency (AF)
      data[, AF := (Falt + Ralt) / (Fref + Rref + Falt + Ralt)]
      
      # Rename columns with genotype prefix
      setnames(data, old = setdiff(names(data), c("CHROM", "POS")), new = paste0(genotype, "_", setdiff(names(data), c("CHROM", "POS"))))
      
      # Store processed data in the list
      geno_data[[genotype]] <- data
    } else {stop(paste("File for", genotype, "not found"))}
  }
  
  # Ensure both wild-type and mutant data are present
  if (is.null(geno_data$wt) || is.null(geno_data$mt)) {stop("Both wild-type and mutant data are required.")}
  
  message("successfully created list for results list")
  return(geno_data)
}

# Example Usage
# vcfiles : "b73_MB1_snps.txt""b73_MB2_snps.txt"
# vcf_dir <- "path/data"
# dir.exists(vcf_dir)
# prefix <- c("b73")
# pattern <- "snps\\.txt$"
# vcf_list <- list.files(path = vcf_dir, pattern, full.names = TRUE)
# print(vcf_list)
# min_DP <- 10
# min_QUAL <- 30
# wt <- "MB1"
# mt <- "MB2"
# geno_data <- load_vcfdata (vcf_dir, prefix, pattern, Genotype=c(wt=wt, mt=mt), min_DP, min_QUAL)
# str(geno_data)
