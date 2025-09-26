window_bsa_compute <- function(data, bsa_col = NULL, window_size = 2e6, step_size = 1e5, rollmedian = 100L,
                               find_interval = TRUE, offhold = 0.80, min_vsize = 1e6, use_cols = c("lft","ma")){
  
  naturalsort <- function(x) {
    if (requireNamespace("stringr", quietly = TRUE)) stringr::str_sort(x, numeric = TRUE) else sort(x)
  }
  
  if (!all(c("CHROM","POS") %in% names(data))) stop("Missing required columns: CHROM and POS")
  if (!(af_col %in% names(data))) stop("af_col not found in data (use 'mt_AF' or 'wt_AF')")
  has_locfit <- requireNamespace("locfit", quietly = TRUE)
  
  datax <- data.table::as.data.table(data)
  datax[, POS := as.integer(POS)]
  
  stat_col   <- paste0(metric_col, "_wmd")  # window median
  ma_col     <- paste0(metric_col, "_rmd")  # rolling median of the window medians
  smooth_col <- paste0(metric_col, "_lft")  # locfit
  
  k <- max(1L, as.integer(ma_k))
  
  results <- list()
  
  for (chr in unique(datax$CHROM)) {
    chr_snps <- datax[CHROM == chr]
    if (!nrow(chr_snps)) next                      # safety
    chr_len  <- max(chr_snps$POS, na.rm = TRUE)
    
    half    <- floor(window_bp / 2)
    centers <- seq(half, chr_len - half, by = step_bp)
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
    chr_res[, (ma_col) := zoo::rollapply(get(stat_col), width = k,
                                         FUN = stats::median, na.rm = TRUE, fill = NA_real_)]
    
    x <- chr_res$POS
    y <- chr_res[[stat_col]]
    keep <- is.finite(x) & is.finite(y)
    
    if (has_locfit && sum(keep) >= 5) {
      fit <- locfit::locfit.raw(x[keep], y[keep], alpha = c(nn = n_p))
      chr_res[, (smooth_col) := as.numeric(predict(fit, newdata = x))]
      chr_res[, (smooth_col) := pmax(0, get(smooth_col))]  # always clamp non-negative
    } else {
      if (!has_locfit) message("locfit not installed; Smooth set to NA for CHROM = ", chr)
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
