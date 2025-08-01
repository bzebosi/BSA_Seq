#' Import VCF Data
#' @description Imports and processes VCF files for specific genotypes.
#' @param vcf_dir Directory that contains the TSV / VCF-like files.
#' @param prefix Common prefix for the VCF files.
#' @param pattern File name pattern for VCF files.
#' @param Genotypes A named list of genotype suffixes.
#' @param min_DP Minimum read depth to keep a record.
#' @param min_QUAL Minimum QUAL score to keep a record.
#' @param only_mutant Logical; if TRUE read only the mutant file.
#' @return A list of filtered and processed data tables.
#' @examples
#' \dontrun{
#' geno <- import_vcfdata(
#'   vcf_dir = "/path/to/files", prefix = "b73", pattern = "snps\\.tsv$",
#'   Genotypes = list(wt = "S7A", mt = "S7B"), min_DP = 10, min_QUAL = 30)
#' }
#' @export

import_vcfdata <- function(vcf_dir, prefix, pattern, Genotypes = list(wt = "wildtype", mt = "mutant"),
                           min_DP = 5, min_QUAL = 5, only_mutant = FALSE) {
  
  # Get matching VCF files
  vcf_list <- list.files(path = vcf_dir, pattern = pattern, full.names = TRUE)
  if (length(vcf_list) == 0) stop("No vcf files found.")
  
  # Decide whether to load only mutant or both genotypes
  selected_genotypes <- if (only_mutant) list(mt = Genotypes[["mt"]]) else Genotypes
  geno_data <- list()
  
  for (genotype in names(selected_genotypes)) {
    file_pattern <- paste0("^", prefix, "_", selected_genotypes[[genotype]])
    matched_file <- vcf_list[grepl(file_pattern, basename(vcf_list))]
    
    if (length(matched_file) == 0) {
      stop(paste("File for", genotype, "not found."))
    } else if (length(matched_file) > 1) {
      warning(paste("Multiple files matched for", genotype, "- using the first one:", basename(matched_file[1])))
    }
    
    geno_file <- matched_file[1]
    message(paste("Reading file for", genotype, ":", basename(geno_file)))
    
    # Read file
    data <- tryCatch({
      fread(geno_file)
    }, error = function(e) {
      stop(paste("Error reading file for", genotype, ":", geno_file, ":", e$message))
    })
    
    if (is.null(data) || nrow(data) == 0) {stop(paste("File for", genotype, "is empty or invalid"))}
    
    message(paste("Loaded", genotype, "data with", nrow(data), "rows."))
    
    # Rename columns
    if (ncol(data) >= 10) {
      expectedcol_names <- c("CHROM", "POS", "REF", "ALT", "QUAL", "DP", "Fref", "Rref", "Falt", "Ralt")
      setnames(data, old = colnames(data)[1:10], new = expectedcol_names)
      if (!is.data.table(data)) data <- as.data.table(data)
    } else {
      stop("The input file does not contain the expected number of columns (10).")
    }
    
    # Keep only SNPs (exclude indels)
    data <- data[nchar(REF) == 1 & nchar(ALT) == 1]
    
    if (!is.data.table(data)) data <- as.data.table(data)
    
    # Filter by depth and quality
    data <- data[!is.na(DP) & DP > min_DP & QUAL >= min_QUAL]
    message(paste("Filtered", genotype, "data with", nrow(data), "rows based on DP and QUAL."))
    
    # Compute allele frequency (AF)
    data[, AF := fifelse(
      (Fref + Rref + Falt + Ralt) > 0, 
      (Falt + Ralt) / (Fref + Rref + Falt + Ralt), NA_real_) ]
    
    # Sort and prefix columns
    setorder(data, CHROM, POS)
    cols_to_prefix <- setdiff(names(data), c("CHROM", "POS"))
    setnames(data, old = cols_to_prefix, new = paste0(genotype, "_", cols_to_prefix))
    
    geno_data[[genotype]] <- data
  }
  
  # Ensure both genotypes present if not mutant-only mode
  if (!only_mutant && (is.null(geno_data$wt) || is.null(geno_data$mt))) {
    stop("Both wild-type and mutant data are required.")
  }
  
  message("Successfully created list for results.")
  return(geno_data)
}
