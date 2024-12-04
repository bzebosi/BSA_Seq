source(file = "~/Documents/BSA_Ts3/customized_Rscripts/multi_package_installer.R")
source(file = "~/Documents/BSA_Ts3/customized_Rscripts/bsa_theme.R")
# function to load and visualize data
loadvcfdata <- function(vcf_dir, inbred, prefix, wt, mt){
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
    col_names <- c("CHROM","POS","REF","ALT","QUAL","DP","Fref","Rref","Falt","Ralt","GT","Tref","Talt","AF")
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
  # Return both datasets as a list
  return(result = list(wt_mt = wt_mt, ant_mt = ant_mt, ant_wt = ant_wt))
}

data="~/Documents/BSA_Ts3/Ts3/data/"
inbreds <- c("B73")
prefixes <- c("b73")
wt <- c("AD8")
mt <- c("AD7")

for (i in seq_along(inbreds)) {
  inbx <- inbreds[i]
  pfix <- prefixes[i]
  print(paste("Processing inbred:", inbx))
  assign(pfix, loadvcfdata(vcf_dir=data, inbred=inbx, prefix=pfix, wt=wt, mt=mt))
}