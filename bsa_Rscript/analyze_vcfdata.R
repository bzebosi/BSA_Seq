#' Analyze BSA-Seq SNP Data
#' Performs variant comparison and statistical analysis between mutant and wild-type SNP tables 
#' from Bulk Segregant Analysis (BSA-Seq) experiments. Computes allele frequency difference (AFD), 
#' Euclidean distance (ED), ED⁴, and G-statistics. Also identifies EMS-type SNPs and unique variants 
#' in each bulk. Supports mutant-only or WT-vs-MT modes.
#' @param geno_data List from import_vcfdata() with $wt and/or $mt.
#' @param prefix Character prefix to label output files.
#' @param save_results Logical. If TRUE, results are saved to disk in CSV and RDS format.
#' @param output_dir Directory where results will be saved.
#' @param only_mutant Logical. If TRUE, assumes mutant-only analysis (no wild-type comparison).
#' @param bsa_metrics Character vector among c("afd","ed","g","all").
#'        - NULL/empty => compute nothing; "all" => afd+ed+g.
#' @return List with wt_mt, ant_mt, ant_wt, ant_mt_ems, ant_wt_ems (subset depends on mode).
#' @export
analyze_vcfdata <- function(geno_data, prefix, save_results = FALSE, output_dir = "post_analysis", 
                            bsa_metrics = c("afd", "g","ed", "all"), only_mutant = FALSE) {
  
  naturalsort <- function(x) {
    if (requireNamespace("stringr", quietly = TRUE)) stringr::str_sort(x, numeric = TRUE) else sort(x)
  }
  
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
    unique_snps <- dplyr::anti_join(data1, data2, by = c("CHROM", "POS"))
    unique_snps <- data.table::as.data.table(unique_snps)
    unique_snps[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
    unique_snps <- unique_snps[order(CHROM, POS)]
    return(unique_snps)
  }
  
  # normalize metric selection
  if (is.null(bsa_metrics) || length(bsa_metrics) == 0L) {
    bsa_metrics <- character(0)  # compute nothing
  } else {
    bsa_metrics <- unique(tolower(bsa_metrics))
    if ("all" %in% bsa_metrics) bsa_metrics <- c("afd","ed","g")
  }
  
  # Process mutant-only mode
  if (only_mutant) {
    if (is.null(geno_data$mt)) stop("Mutant data is required.")
    mt_data <- data.table::as.data.table(geno_data$mt)
    mt_data[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
    mt_data <- mt_data[order(CHROM, POS)]
    
    ant_mt_ems <- get_ems(mt_data, "mt_REF", "mt_ALT", "mutant")
    result <- list(ant_mt = mt_data, ant_mt_ems = ant_mt_ems)
    
    if (save_results) {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      saveRDS(result, file.path(output_dir, paste0(prefix, "_mutant_only_results.rds")))
      data.table::fwrite(mt_data, file.path(output_dir, paste0(prefix, "_ant_mt.csv")))
      data.table::fwrite(ant_mt_ems, file.path(output_dir, paste0(prefix, "_ant_mt_ems.csv")))
      message("Mutant-only results saved to ", output_dir)
    }
    return(result)
  } else {
      # === WT vs MT analysis ===
      if (is.null(geno_data$wt) || is.null(geno_data$mt)) {stop("Both wild-type and mutant datasets are required.")}
      wt_data <- data.table::as.data.table(geno_data$wt)
      mt_data <- data.table::as.data.table(geno_data$mt)
      
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
      if ("afd" %in% bsa_metrics) {
        message("Calculating AFD")
        if (!all(c("wt_AF","mt_AF") %in% names(wt_mt)))
          stop("AFD requires columns wt_AF and mt_AF in wt_mt.")
        wt_mt[, AFD := abs(wt_AF - mt_AF)]
      }
      
      # Calculate Euclidean Distance (ED) and its fourth power (ED4)
      if ("ed" %in% bsa_metrics) {
        message("Calculating ED and ED4")
        if (!all(c("wt_AF","mt_AF") %in% names(wt_mt)))
          stop("ED requires columns wt_AF and mt_AF in wt_mt.")
        
        wt_mt[, ED := sqrt(2) * abs(wt_AF - mt_AF)]
        wt_mt[, ED4 := ED^4]
      }
      
      # Calculate G-statistics for each SNP
      if ("g" %in% bsa_metrics) {
        message("Calculating G-stats")
        need_g <- c("mt_Fref","mt_Rref","mt_Falt","mt_Ralt","wt_Fref","wt_Rref","wt_Falt","wt_Ralt")
        if (!all(need_g %in% names(wt_mt)))
          stop("G requires columns: ", paste(need_g, collapse = ", "), " in wt_mt.")
        
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
        
        wt_mt[, G := mapply(compute_G,mt_Fref, mt_Rref, mt_Falt, mt_Ralt, wt_Fref, wt_Rref, wt_Falt, wt_Ralt)]
      }
        
      result <- list(wt_mt = wt_mt, ant_mt = ant_mt, ant_wt = ant_wt, ant_mt_ems = ant_mt_ems, ant_wt_ems = ant_wt_ems)
      
      # === Save output ===
      if (save_results) {
        if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
        saveRDS(result, file.path(output_dir, paste0(prefix, "_results.rds")))
        for (name in names(result)) {
          data.table::fwrite(result[[name]], file.path(output_dir, paste0(prefix, "_", name, ".csv")))
        }
        message("All results to", output_dir)
      }
      return(result)
  }
}


# vcf_dir = "/Users/zebosi/Documents/osu_postdoc/BSA/S7_6508K/data/snps"
# output_dir = "/Users/zebosi/Documents/osu_postdoc/BSA/S7_6508K/post_analysis"
# wt <- c("S7A6508K") or wt <- NULL (if only_mutant = TRUE)
# mt <- c("S7B6508K")
# pattern = "snps\\.tsv$"
# min_DP=10
# min_QUAL=10
# prefix = c("b73")
# a <- import_vcfdata(vcf_dir, prefix, pattern, Genotypes = list(wt = wt, mt = mt),
#                    min_DP, min_QUAL, only_mutant = FALSE)
#d <- analyze_vcfdata(a, prefix, save_results = FALSE, bsa_metrics = "all", output_dir, only_mutant = FALSE)
