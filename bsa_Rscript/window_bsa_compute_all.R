#' Window-based BSA analysis for multiple metrics
#'
#' @description
#' Wrapper around window_bsa_compute() to run sliding-window summaries for:
#' - AF (allele frequency) from mutant (mt_AF), wildtype (wt_AF), or both (via ant_*  and ant_*_ems)
#' - Joint metrics ED, ED4, AFD, G from wt_mt.
#' Respects use_ems (use EMS tables) and only_mutant (skip WT + joint metrics).
#' @param data List with tables: ant_mt, ant_wt, optional ant_mt_ems, ant_wt_ems, and wt_mt for joint metrics.
#' @param metrics Character vector: one or more of "af","ed","ed4","afd","g","all".
#' @param af_col Which AF to use: "mt_AF","wt_AF","both". Default "mt_AF".
#' @param use_ems Logical; if TRUE, use *_ems tables for AF. Default FALSE.
#' @param only_mutant Logical; if TRUE, skip WT AF and joint metrics. Default FALSE.
#' @param window_size,step_size,rollmedian,nn_prop,find_intervals,offhold,min_vsize
#' @inheritParams window_bsa_compute
#'   Passed through to window_bsa_compute().
#' @return A list with:

window_bsa_compute_all <- function(data, metrics = c("af","ed","ed4","afd","g","all"),
                               af_col = "mt_AF", use_ems = FALSE, only_mutant = FALSE,
                               window_size = 2e6, step_size = 1e5,rollmedian = 100L, nn_prop = 0.1,
                               find_intervals = TRUE, offhold = 0.80, min_vsize = 1e6){
  
  metrics <- tolower(metrics)
  if ("all" %in% metrics) metrics <- c("af","ed","ed4","afd","g")
  
  out <- list(windows = list(), intervals = list())
  
  af_cols <- if ("both" %in% af_col) c("mt_AF","wt_AF") else af_col
  
  if ("af" %in% metrics) {
    for (fq in af_cols) {
      if (!fq %in% c("mt_AF","wt_AF")) stop("af_col must be 'mt_AF', 'wt_AF', or 'both'.")
      if (only_mutant && fq == "wt_AF") { message("only_mutant=TRUE: skipping wt_AF."); next }
      
      # choose table name based on AF side + EMS flag
      doi <- if (fq == "mt_AF") {
        paste0("ant_mt", if (use_ems) "_ems")
      } else {
        paste0("ant_wt", if (use_ems) "_ems")
      }
      
      dt <- data[[doi]]
      if (is.null(dt)) stop("Missing table: ", doi)
      if (!fq %in% names(dt)) stop("Column ", fq, " not found in ", doi)
      
      res <- window_bsa_compute(data = dt, metric_col = fq, window_size = window_size, 
                                step_size = step_size, rollmedian = rollmedian, 
                                nn_prop = nn_prop, find_intervals = find_intervals, 
                                offhold = offhold, min_vsize = min_vsize, use_cols = c()
      )
      
      out$windows[[fq]]   <- res$windows
      out$intervals[[fq]] <- if ("intervals" %in% names(res)) res$intervals else NULL
    }
  }
  
  # ---- Joint metrics (from wt_mt) ----
  joint <- unique(metrics[metrics %in% c("ed","ed4","afd","g")])
  if (!only_mutant && length(joint)) {
    wt_mt <- data$wt_mt
    if (is.null(wt_mt)) stop("Missing table: wt_mt")
    
    col_map <- c(ed = "ED", ed4 = "ED4", afd = "AFD", g = "G")
    for (m in joint) {
      metric_col <- unname(col_map[[m]])
      if (!metric_col %in% names(wt_mt)) 
        stop("Column ", metric_col, " not found in wt_mt")
      
      res <- window_bsa_compute(
        data = wt_mt, metric_col = metric_col,
        window_size = window_size, step_size = step_size,
        rollmedian = rollmedian, nn_prop = nn_prop,
        find_intervals = find_intervals, offhold = offhold, min_vsize = min_vsize,
        use_cols = c()
      )
      out$windows[[m]]   <- res$windows
      out$intervals[[m]] <- if ("intervals" %in% names(res)) res$intervals else NULL
    } 
    
  } else if (only_mutant && length(joint)) {
    message("only_mutant=TRUE: skipping joint metrics (", paste(joint, collapse = ", "), ").")
  }
  
  return(out)
}
