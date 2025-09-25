#' Peak intervals from window scores
#'
#' @description
#' For each chromosome and each column in use_cols, find the global peak,
#' expand left/right while score ≥ offhold * peak, and return that interval.
#' Intervals shorter than min_vsize (bp) are dropped.
#'
#' @param data data.frame/data.table with CHROM,start,end,POS` and score columns.
#' @param offhold numeric in (0,1]; cutoff as a fraction of the peak (default 0.80).
#' @param min_vsize integer; minimum interval length in bp.
#' @param use_cols character; score columns to scan (e.g. "Homozygosity_lft").
#'
#' @return data.table with:
#' `CHROM, score_col, cutoff, start, end, width_bp, peak_pos, peak_score,
#'  mean_score, max_score, min_score, area (where area = max_score * width_bp).
#' @export
window_peak_interval <- function(data, offhold=0.80, min_vsize = 0L, use_cols = c()) {
  # input checks
  datax <- as.data.table(data)
  need <- c("CHROM","start","end","POS")
  if (!all(need %in% names(datax))) stop("Input needs CHROM, start, end, POS")
  
  cols_use <- intersect(use_cols, names(datax))
  if (!length(cols_use)) stop("None of use_cols found in data")
  
  # empty result (data.frame)
  out <- data.table(CHROM=character(), score_col=character(), cutoff=numeric(),
    start=integer(), end=integer(), width_bp=integer(),peak_pos=integer(), 
    peak_score=numeric(), mean_score=numeric(), max_score=numeric(), area=numeric())
  
  results <- list()
  
  for (chr in unique(datax$CHROM)) {
    dt <- datax[CHROM == chr][order(POS)]
    if (!nrow(dt)) next
    
    for (col in cols_use) {
      s <- dt[[col]]
      if (all(is.na(s))) next
      
      # peak and cutoff
      s2 <- ifelse(is.na(s), -Inf, s)
      peak_idx <- which.max(s2)
      peak_val <- s2[peak_idx]
      if (!is.finite(peak_val)) next
      cutoff <- peak_val * offhold
      
      # expand from peak while score >= cutoff
      L <- R <- peak_idx; n <- nrow(dt)
      while (L > 1L && s2[L - 1L] >= cutoff) L <- L - 1L
      while (R < n  && s2[R + 1L] >= cutoff) R <- R + 1L
      
      # build segment + size filter
      seg <- dt[L:R]
      seg_start <- as.integer(min(seg$start, na.rm = TRUE))
      seg_end   <- as.integer(max(seg$end,   na.rm = TRUE))
      width_bp  <- as.integer(seg_end - seg_start + 1L)
      if (width_bp < as.integer(min_vsize)) next
      
      # collect
      mean_sc <- as.numeric(mean(seg[[col]], na.rm = TRUE))
      max_sc  <- as.numeric(max(seg[[col]],  na.rm = TRUE))
      min_sc  <- as.numeric(min(seg[[col]],  na.rm = TRUE))
                            
      results[[length(results) + 1L]] <- data.frame(CHROM = seg$CHROM[1], 
        score_col= col, cutoff = cutoff, start = seg_start,end = seg_end, 
        width_bp = width_bp, peak_pos = as.integer(dt$POS[peak_idx]),
        peak_score = as.numeric(s[peak_idx]),
        mean_score = mean_sc,
        max_score  = max_sc,
        min_score  = min_sc,
        area = max_sc * width_bp,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (!length(results)) return(out)
  rbindlist(results, use.names = TRUE)
}


#' window_peak_interval(data, offhold=0.8, min_vsize=1e6,
#'                    use_cols=c("Homozygosity_lft","Homozygosity_ma"))