#' Sliding-window homozygosity (with optional peak intervals)
#' Computes windowed homozygosity along chromosomes from SNP AFs, plus
#' optional smoothed/rolled summaries and peak intervals.
#' @param data data.frame/data.table with at least CHROM, POS, and an AF column.
#' @param af_col Character. Name of the allele-frequency column to use (e.g. "mt_AF" or "wt_AF").
#' @param window_size,step_size Numeric window size and step (bp).
#' @param rollmedian Integer width for rolling median over window series.
#' @param nn_prop Numeric; nearest-neighbor fraction for locfit smoothing.
#' @param af_min Numeric in (0,1). SNPs with AF >= af_min (or <= 1-af_min) count as homozygous.
#' @param find_intervals Logical; if TRUEc, call window_peak_interval()} on the window table.
#' @param offhold Numeric in (0,1). Interval cutoff = peak * offhold.
#' @param min_vsize Integer. Minimum genomic width (bp) to keep an interval.
#' @param use_cols Character vector of score columns (in the window table) to scan for peaks.
#' use_cols= c("Homozygosity_lft","Homozygosity_ma","Homozygosity")
#' @return If find_intervals=FALSE: list(windows = data.table).
#' @seealso window_peak_interval
#' @export
#' 
window_homozygosity_compute <- function(data, af_col = NULL, window_size = 5e5, step_size = 1e5, rollmedian = 21L, nn_prop = 0.1, af_min = 0.99, 
                                find_intervals = TRUE, offhold = 0.80, min_vsize = 1e6, use_cols= c("Homozygosity_lft","Homozygosity_ma","Homozygosity")) {
 
  naturalsort <- function(x) {
    if (requireNamespace("stringr", quietly = TRUE)) stringr::str_sort(x, numeric = TRUE) else sort(x)
  }
  
  if (!all(c("CHROM","POS") %in% names(data))) stop("Missing required columns: CHROM and POS")
  if (!(af_col %in% names(data))) stop("af_col not found in data (use 'mt_AF' or 'wt_AF')")
  has_locfit <- requireNamespace("locfit", quietly = TRUE)
  
  datax <- data.table::as.data.table(data)
  datax[, POS := as.integer(POS)]
  datax[, Hom := as.integer(get(af_col) >= af_min | get(af_col) <= (1 - af_min))]
  
  ma_col     <- paste0("Homozygosity", "_rmd")  # rolling median of the window medians
  smooth_col <- paste0("Homozygosity", "_lft")
  
  results <- list()
  for (chr in unique(datax$CHROM)) {
    chr_snps <- datax[CHROM == chr]
    if (nrow(chr_snps) == 0L) next
    chr_len  <- max(chr_snps$POS, na.rm = TRUE)
    half    <- floor(window_size / 2)
    centers <- seq(half, chr_len - half, by = step_size)
    
    if (length(centers) == 0L) next
    out_rows <- lapply(centers, function(center) {
      start <- max(1L, center - half)
      end   <- min(chr_len, center + half)
      win   <- chr_snps[POS >= start & POS <= end]
      total <- nrow(win)
      hom   <- if (total > 0L) sum(win$Hom) else 0L
      hz    <- if (total > 0L) hom / total else NA_real_
      data.table::data.table(
        CHROM = chr, start = start, end = end, POS = center,
        N = as.integer(total), Hom = as.integer(hom), Homozygosity = hz
      )
    })
    
    k <- max(1L, as.integer(rollmedian))
    chr_res <- data.table::rbindlist(out_rows, use.names = TRUE, fill = TRUE)
    data.table::setorder(chr_res, POS)
    
    chr_res[, Homozygosity_ma := zoo::rollapply(
      Homozygosity, width = k, FUN = stats::median, na.rm = TRUE, fill = NA_real_
    )]
    
    x <- chr_res$POS
    y <- chr_res$Homozygosity
    keep <- is.finite(x) & is.finite(y)
    
    if (has_locfit && sum(keep) >= 5) {
      fit <- locfit::locfit.raw(x[keep], y[keep], alpha = c(nn = nn_prop))
      chr_res[, (smooth_col) := as.numeric(predict(fit, newdata = x))]
      chr_res[, (smooth_col) := pmax(0, pmin(1, get(smooth_col)))]
      chr_res[, (smooth_col) := pmax(0, get(smooth_col))] 
    } else {
      if (!has_locfit) message("locfit not installed; Homozygosity_smooth set to NA for CHROM = ", chr)
      chr_res[, (smooth_col) := NA_real_]
    }
    results[[length(results) + 1L]] <- chr_res
  }
  
  out <- data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
  out[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
  data.table::setorder(out, CHROM, POS)
  
  # === optional: compute genomic peak intervals from the window table ===
  if (isTRUE(find_intervals)) {
    if (!exists("window_peak_interval", mode = "function")) {
      stop("window_peak_interval=TRUE but window_peak_interval() is not available in the environment.")
    }
    
    intervals <- window_peak_interval(data=out, offhold = offhold, min_vsize= min_vsize, use_cols = use_cols)

    return(list(windows = out, intervals = intervals))
  } else {
    return(list(windows = out))
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
# d <- analyze_vcfdata(a, prefix, save_results = FALSE, bsa_metrics = "all", output_dir, only_mutant = FALSE)

#, win <- compute_window_homozygosity(
#'  data = d$ant_mt, af_col = "mt_AF",
#'  window_size = 5e5, step_size = 1e5, find_intervals = TRUE,
#'  rollmedian = 101L, af_min = 0.99, nn_prop = 0.1, offhold = 0.80, min_vsize = 1e6,
#'   use_cols = c("Homozygosity_lft","Homozygosity_ma","Homozygosity")
# )




