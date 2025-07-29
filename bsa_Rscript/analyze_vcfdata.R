analyze_vcfdata <- function(geno_data, prefix, save_results = FALSE, output_dir = "post_analysis", only_mutant = FALSE) {
  # Helper function : Filter and Extract EMS SNPs from mutant
  get_ems <- function(data, ref_col, alt_col, label) {
    message(paste("Filtering EMS SNPs in", label))
    ems_variants <- data[(get(ref_col) == "G" & get(alt_col) == "A") | (get(ref_col) == "C" & get(alt_col) == "T")]
    ems_variants[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
    ems_variants <- ems_variants[order(CHROM, POS)]
    message("Number of Total SNPs in ", label, ":", nrow(data))
    message("Number of EMS SNPs in ", label, ":", nrow(ems_variants))
    return(ems_variants)
  }
  
  # Helper Function :  anti-join and identity unique SNPs
  id_unique_snps <- function(data1, data2, label) {
    message(paste("Identifying SNPs unique to", label))
    unique_snps <- anti_join(data1, data2, by = c("CHROM", "POS"))
    unique_snps[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
    unique_snps <- unique_snps[order(CHROM, POS)]
    return(unique_snps)
  }

  # Process mutant-only mode
  if (only_mutant) {
    if (is.null(geno_data$mt)) stop("Mutant data is required.")
    mt_data <- geno_data$mt
    mt_data[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
    mt_data <- mt_data[order(CHROM, POS)]
    
    ant_mt_ems <- get_ems(mt_data, "mt_REF", "mt_ALT", "mutant")
    result <- list(ant_mt = mt_data, ant_mt_ems = ant_mt_ems)
    
    if (save_results) {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      saveRDS(result, file.path(output_dir, paste0(prefix, "_mutant_only_results.rds")))
      fwrite(mt_data, file.path(output_dir, paste0(prefix, "_ant_mt.csv")))
      fwrite(ant_mt_ems, file.path(output_dir, paste0(prefix, "_ant_mt_ems.csv")))
      message("Mutant-only results saved to ", output_dir)
    }
    return(result)
  } else {
      # === WT vs MT analysis ===
      if (is.null(geno_data$wt) || is.null(geno_data$mt)) {stop("Both wild-type and mutant datasets are required.")}
      wt_data <- geno_data$wt
      mt_data <- geno_data$mt
      
      # Merge wild-type and mutant by CHROM and POS
      message("Merging wild-type and mutant datasets...")
      wt_mt <- merge(wt_data, mt_data, by = c("CHROM", "POS"))
      wt_mt[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
      wt_mt <- wt_mt[order(CHROM, POS)]
      message("Number of shared SNPs (wt_mt): ", nrow(wt_mt))
      
      # Unique SNPs and EMS filtering
      ant_wt <- id_unique_snps(wt_data, mt_data, "wildtype")
      ant_mt <- id_unique_snps(mt_data, wt_data, "mutant")
      ant_wt_ems <- get_ems(ant_wt, "wt_REF", "wt_ALT", "wildtype")
      ant_mt_ems <- get_ems(ant_mt, "mt_REF", "mt_ALT", "mutant")
      
      #  # Statistical metrics : Allele frequency diff, ED, G-statistics, p-values
      message("Calculating allele frequency differences (AFD)")
      wt_mt[, AFD := abs(wt_AF - mt_AF)]
      
      # Calculate Euclidean Distance (ED) and its fourth power (ED4)
      message("Calculating allele Euclidean Distance (ED) and its fourth power (ED4)")
      wt_mt[, ED := sqrt(2) * abs(wt_AF - mt_AF)]
      wt_mt[, ED4 := ED^4]
      
      # Calculate G-statistics for each SNP
      message("Calculating G-statistics for each SNP...")
      compute_G <- function(mt_Fref, mt_Rref, mt_Falt, mt_Ralt,
                                 wt_Fref, wt_Rref, wt_Falt, wt_Ralt) {
        # Compute allele counts
        n_Am <- mt_Fref + mt_Rref   # reference allele in mutant bulk
        n_am <- mt_Falt + mt_Ralt   # alternate allele in mutant bulk
        n_Aw <- wt_Fref + wt_Rref   # reference allele in wildtype bulk
        n_aw <- wt_Falt + wt_Ralt   # alternate allele in wildtype bulk
        
        total <- n_Am + n_am + n_Aw + n_aw
        if (total == 0) return(NA_real_)
        
        # Expected counts under independence
        e1 <- (n_Am + n_Aw) * (n_Am + n_am) / total
        e2 <- (n_Am + n_Aw) * (n_Aw + n_aw) / total
        e3 <- (n_am + n_aw) * (n_Am + n_am) / total
        e4 <- (n_am + n_aw) * (n_Aw + n_aw) / total
        
        obs <- c(n_Am, n_Aw, n_am, n_aw)
        exp <- c(e1, e2, e3, e4)
        
        # Avoid division by zero or log(0)
        obs[obs == 0] <- 1e-10
        exp[exp == 0] <- 1e-10
        
        return(2 * sum(obs * log(obs / exp), na.rm = TRUE))
      }
      
      wt_mt[, G := mapply(compute_G,
                          mt_Fref, mt_Rref, mt_Falt, mt_Ralt,
                          wt_Fref, wt_Rref, wt_Falt, wt_Ralt)]
      
        
      result <- list(wt_mt = wt_mt, ant_mt = ant_mt, ant_wt = ant_wt, ant_mt_ems = ant_mt_ems, ant_wt_ems = ant_wt_ems)
      
      # === Save output ===
      if (save_results) {
        if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
        saveRDS(result, file.path(output_dir, paste0(prefix, "_results.rds")))
        for (name in names(result)) {
          fwrite(result[[name]], file.path(output_dir, paste0(prefix, "_", name, ".csv")))
        }
        message("All results to", output_dir)
      }
      return(result)
  }
}
