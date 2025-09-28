#' Sliding window computation for BSA metrics
#' @description
#' Compute median values of a metric in sliding windows across chromosomes.
#' Also adds a rolling median and optional locfit smoothing, and can call
#' `window_peak_interval()` to detect peak regions.
#' @param data Data.frame/data.table with `CHROM`, `POS`, and the metric column.
#' @param metric_col Name of the numeric metric column to use.
#' @param window_size Window size in bp (default 2e6).
#' @param step_size Step size in bp (default 1e5).
#' @param rollmedian Rolling median window (default 100).
#' @param nn_prop Proportion for locfit smoothing (default 0.1).
#' @param find_intervals Logical, run `window_peak_interval()` (default TRUE).
#' @param offhold Peak cutoff fraction (default 0.8).
#' @param min_vsize Minimum bp size for intervals (default 1e6).
#' @param use_cols Columns to use for interval calling. Defaults to
#'   `c(metric_wmd, metric_rmd, metric_lft)`.
#' @return A list with:
#' \itemize{
#'   \item `windows`: window-level summary table
#'   \item `intervals`: peak interval table (if `find_intervals=TRUE`)
#' }

window_bsa_compute_one <- function(data, metric_col = NULL, window_size = 2e6, step_size = 1e5, rollmedian = 100L, nn_prop = 0.1,
                               find_intervals = TRUE, offhold = 0.80, min_vsize = 1e6, use_cols = c()){
  
  naturalsort <- function(x) {
    if (requireNamespace("stringr", quietly = TRUE)) stringr::str_sort(x, numeric = TRUE) else sort(x)
  }
  
  if (!all(c("CHROM","POS") %in% names(data))) stop("Missing required columns: CHROM and POS")
  if (is.null(metric_col) || !(metric_col %in% names(data))) stop("metric_col not found in data")
  if (!is.numeric(data[[metric_col]])) stop("bsa_col must be numeric")

  has_locfit <- requireNamespace("locfit", quietly = TRUE)
  
  datax <- data.table::as.data.table(data)
  datax[, POS := as.integer(POS)]
  
  stat_col   <- paste0(metric_col, "_wmd")  # window median
  ma_col     <- paste0(metric_col, "_rmd")  # rolling median of the window medians
  smooth_col <- paste0(metric_col, "_lft")  # locfit
  
  k <- max(1L, as.integer(rollmedian))
  # pre-compute chromosome ordering from input (safer)
  all_chr <- naturalsort(unique(as.character(datax$CHROM)))
  
  results <- list()
  
  # per-chromosome window computation 
  for (chr in unique(datax$CHROM)) {
    chr_snps <- datax[CHROM == chr]
    if (!nrow(chr_snps)) next                      # safety
    chr_len  <- max(chr_snps$POS, na.rm = TRUE)
    
    half    <- floor(window_size / 2)
    centers <- seq(half, chr_len - half, by = step_size)
    if (length(centers) == 0L) next
    
    out_rows <- lapply(centers, function(center) {
      start <- max(1L, center - half)
      end   <- min(chr_len, center + half)
      win   <- chr_snps[POS >= start & POS <= end]
      total <- nrow(win)
      val <- if (total > 0L) stats::median(as.numeric(win[[metric_col]]), na.rm = TRUE) else NA_real_
      data.table::data.table(CHROM = chr, start = start, end = end, POS = center,
                             N = as.integer(total), Stat = val)
    })
    
    chr_res <- data.table::rbindlist(out_rows, use.names = TRUE, fill = TRUE)
    data.table::setorder(chr_res, POS)
    data.table::setnames(chr_res, "Stat", stat_col)
    
    # rolling median across window medians
    chr_res[, (ma_col) := zoo::rollapply(get(stat_col), width = k,
                                         FUN = stats::median, na.rm = TRUE, fill = NA_real_)]
    
    # locfit smoothing
    x <- chr_res$POS
    y <- chr_res[[stat_col]]
    keep <- is.finite(x) & is.finite(y)
    
    if (has_locfit && sum(keep) >= 5) {
      fit <- locfit::locfit.raw(x[keep], y[keep], alpha = c(nn = nn_prop))
      chr_res[, (smooth_col) := as.numeric(predict(fit, newdata = x))]
      chr_res[, (smooth_col) := pmax(0, get(smooth_col))]  # always clamp non-negative
    } else {
      if (!has_locfit) message("locfit not installed; Smooth set to NA for CHROM = ", chr)
      chr_res[, (smooth_col) := NA_real_]
    }
    results[[length(results) + 1L]] <- chr_res
  }
  
  
  # ---- combine chromosomes; handle empty edge case ----
  if (!length(results)) {
    out <- data.table::data.table(
      CHROM = character(), start = integer(), end = integer(),POS = integer(), N = integer()
    )
    out[[stat_col]] <- numeric(); out[[ma_col]] <- numeric(); out[[smooth_col]] <- numeric()
    if (find_intervals) {
      return(list(windows = out, intervals = out[0]))
    } else {
      return(list(windows = out))
    }
  }
  
  out <- data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
  out[, CHROM := factor(CHROM, levels = all_chr)]
  data.table::setorder(out, CHROM, POS)
  
  # === optional: compute genomic peak intervals from the window table ===
  if (find_intervals) {
    if (!exists("window_peak_interval", mode = "function")) {
      stop("window_peak_interval=TRUE but window_peak_interval() is not available in the environment.")
    }
    
    use_cols_final <- if (length(use_cols)) use_cols else c(stat_col, ma_col, smooth_col)
    intervals <- window_peak_interval(data=out, offhold = offhold, min_vsize= min_vsize, use_cols = use_cols_final)
    
    return(list(windows = out, intervals = intervals))
  } else {
    return(list(windows = out))
  }
}

#' vcf_dir = "/Users/zebosi/Documents/osu_postdoc/BSA/S7_6508K/data/snps"
#' output_dir = "/Users/zebosi/Documents/osu_postdoc/BSA/S7_6508K/interval_analysis"
#' plots_dir = "/Users/zebosi/Documents/osu_postdoc/BSA/S7_6508K/ppl_analysis"
#' wt <- c("S7A6508K")
#' mt <- c("S7B6508K")
#' pattern = "snps\\.tsv$"
#' threshold = -log10(0.05) * 10
#' min_DP=10
#' min_QUAL=10
#' prefix = c("b73")
#' 
#' a <- import_vcfdata(vcf_dir, prefix, pattern, Genotypes = list(wt = wt, mt = mt),
#'                    min_DP = 5, min_QUAL = 5, only_mutant = FALSE)
#' d <- analyze_vcfdata(a, prefix, save_results = FALSE, bsa_metrics = "all", output_dir = "post_analysis", only_mutant = FALSE)
#' win <- window_bsa_compute_one(d$wt_mt, metric_col = "ED4", window_size = 2e6, step_size = 1e5, 
#'                          rollmedian = 100L, nn_prop = 0.1, find_intervals = TRUE, offhold = 0.90, min_vsize = 1e6)
