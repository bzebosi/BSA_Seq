load_vcfdata <- function (vcf_dir, prefix, pattern, Genotype, min_DP, min_QUAL){
  
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
      
      # check data
      if (is.null(data) || nrow(data) == 0) {stop(paste("File for", genotype, "is either empty or missing"))}
      message(paste("Loaded", genotype, "data with", nrow(data), "rows."))
      
      # Rename columns
      col_names <- c("CHROM", "POS", "REF", "ALT", "QUAL", "DP", "Fref", "Rref", "Falt", "Ralt")
      colnames(data) <- col_names
      
      # Convert to data.table
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
