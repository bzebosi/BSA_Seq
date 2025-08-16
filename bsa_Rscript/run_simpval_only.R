#' Sliding-Window Fisher’s Exact Test for BSA-Seq
#'
#' For each window along each chromosome, performs Fisher’s exact test on allele counts,
#' adjusts p-values (Bonferroni), and identifies significant intervals.  
#' Optionally saves per-window results, significant SNPs, and merged intervals as CSV or Excel,
#' and plots the –log10(p-value) along the genome.
#'
#' @param wt_mt Data.table of merged WT vs MT SNPs (must contain \code{AFD}, \code{ED}, or \code{ED4}).
#' @param ant_wt Data.table of SNPs unique to WT bulk.
#' @param ant_mt Data.table of SNPs unique to MT bulk.
#' @param ant_wt_ems Data.table of EMS-type SNPs unique to WT bulk.
#' @param ant_mt_ems Data.table of EMS-type SNPs unique to MT bulk.
#' @param wt Character label for wild-type bulk (used in filenames).
#' @param mt Character label for mutant bulk.
#' @param prefix Sample/prefix string (used in output filenames).
#' @param plots_dir Directory in which to save plot images.
#' @param plot_style One of \code{"grid"} or \code{"wrap"} layouts for faceting.
#' @param ylim NULL or numeric \code{c(min, max)} for y-axis limits.
#' @param rollmedian Window size for rolling-median smoothing (ignored in this function).
#' @param threshold Numeric cutoff on –log10(adj.p-value) to flag significance.
#' @param output_dir Directory in which to save CSV/Excel results.
#' @param save_intervals Logical; if \code{TRUE} saves result tables to disk.
#' @param width,height,hwidth,hheight Numeric dimensions (in inches) for saved plots.
#' @param dpi Resolution of saved plot images.
#' @param device Graphics device (\code{"png"}, \code{"pdf"}, etc.).
#' @param plot_data Logical; if \code{TRUE}, generates and saves the sliding-window p-value plot.
#' @param only_mutant Logical; if \code{TRUE}, only the MT-unique tables will be scanned (no WT data).
#' @param window_size Window width in bp (default 1e6).
#' @param step_size Step size in bp between windows (default 1e5).
#' @param is_histogram Ignored here.
#' @param save_excel Logical; if \code{TRUE} writes results into a single Excel workbook.
#' @param stat_method Character vector: which statistics to scan \code{"af"} (AFD) and \code{"ed"} (ED/ED4).
#' @param af_doorstep Numeric AF threshold for defining “high” vs “low”.
#' @param ed_doorstep Numeric ED threshold.
#' @param color_panel Character vector of colors for plotting.
#' @return A named list of length one per statistic run, each with elements:
#'   \itemize{
#'     \item \code{$results}: full sliding-window table,
#'     \item \code{$abovethreshold}: windows exceeding \code{threshold},
#'     \item \code{$intervals}: merged genomic intervals of significance.
#'   }
#' @examples
#' \dontrun{
#' res <- run_simpval_only(
#'   wt_mt=merged_snps, ant_wt=unique_wt, ant_mt=unique_mt, wt="WT", mt="Ts5", prefix="Ts3",
#'   plots_dir="plots/", output_dir="post_analysis/", stat_method=c("af","ed"),
#'   window_size=2e6, step_size=5e5, save_excel=TRUE)
#' }
#' @export
run_simpval_only <- function(wt_mt = NULL, ant_wt = NULL, ant_mt = NULL, ant_wt_ems = NULL, ant_mt_ems = NULL, 
                              wt = "wildtype", mt = "mutant", prefix = "sample", plots_dir = "plots", plot_style = "grid",
                              ylim = NULL, rollmedian = 0, threshold = -log10(0.05) * 10, output_dir = "post_analysis", save_intervals = TRUE, 
                              width = 45, height = 13, hwidth = 36, hheight = 18, dpi = 300, device = "png", plot_data = TRUE,
                              only_mutant = FALSE, window_size = 1000000, step_size = 100000, is_histogram = FALSE, save_excel = TRUE,
                              stat_method = c("af", "ed"), af_doorstep = 0.67, ed_doorstep = 1, color_panel = c("blue", "red")) {
  
  # internal check if the data is not empty.
  
  check_cols <- function(data, colname) {
    if (is.null(data)) { message("data not present"); return(NULL) }
    if (!(colname %in% colnames(data))) { message(paste("Column", colname, "not found in data.")); return(NULL) }
    return(data)
  }
  
  # Inner analysis function
  scan_it <- function(data, af_col, window_size, step_size, method, af_doorstep, ed_doorstep) {
    doorstep <- if (method == "af") af_doorstep else ed_doorstep
    b_high <- sum(data[[af_col]] > doorstep, na.rm = TRUE)
    b_low  <- sum(data[[af_col]] <= doorstep, na.rm = TRUE)
    
    results_out <- list()
    unique_chromosomes <- unique(data$CHROM)
    for (chr in unique_chromosomes) {
      chr_data <- data[CHROM == chr]
      chr_len  <- max(chr_data$POS, na.rm = TRUE)
      starts   <- seq(0, chr_len, by = step_size)
      
      for (start in starts) {
        end  <- start + window_size
        mid  <- (start + end) / 2
        win  <- chr_data[POS >= start & POS <= end]
        if (nrow(win) == 0) next

        hi <- sum(win[[af_col]] > doorstep, na.rm = TRUE)
        lo <- sum(win[[af_col]] <= doorstep, na.rm = TRUE)
        
        mat <- matrix(c(hi, lo, b_high, b_low), nrow = 2)
        p_val <- tryCatch(fisher.test(mat)$p.value, error = function(e) NA_real_)
        
        results_out[[length(results_out) + 1]] <- data.table(
          CHROM = chr, start = start, end = end, POS = mid,
          n_win = nrow(win), high = hi, pval = p_val
        )
      }
    }
    
    # Check if 'pval' column exists
    if (length(results_out) == 0) { message("No valid sliding windows found for: ", af_col); return(data.table())}
    
    
    # Combine all results into a single data table
    results_dt <- rbindlist(results_out, fill = TRUE)
    results_dt[, adj.pval := p.adjust(pval, method = "bonferroni")]
    results_dt[, log.pval := -10 * log10(ifelse(adj.pval == 0, 1e-10, adj.pval))]
    results_dt[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
    results_dt <-  results_dt[order(CHROM, POS)]
    
    # Identify significant windows
    sig_data <- results_dt[log.pval > threshold]
    # Reduce intervals
    if (nrow(sig_data) > 0) {
      reduced_list <- lapply(split(sig_data, sig_data$CHROM), function(subdt) {
        rng <- reduce(IRanges(start = subdt$start, end = subdt$end))
        df <- as.data.frame(rng)
        df$CHROM <- unique(subdt$CHROM)
        return(df[, c("CHROM", "start", "end", "width")])
      })
      
      reduced_intervals <- rbindlist(reduced_list)
    } else {
      reduced_intervals <- data.table(CHROM = character(), start = integer(), end = integer())
    }
    
    return(list(results = results_dt, abovethreshold = sig_data, intervals = reduced_intervals))
  }
  
  # Ensure stat_method is valid
  stat_method <- match.arg(stat_method, choices = c("af", "ed"), several.ok = TRUE)
  if (only_mutant && "ed" %in% stat_method) {
    warning("only_mutant=TRUE: dropping 'ed' from stat_method; ED requires WT+MT.")
    stat_method <- setdiff(stat_method, "ed")
  }
  if (only_mutant && length(stat_method) == 0) {
    message("only_mutant=TRUE: defaulting stat_method to 'af'.")
    stat_method <- "af"
  }
  
  make_parameters <- function() {
    params <- list()
    for (stat in stat_method){
      if (stat == "af"){
        if (only_mutant) {
          params <- append(params,list(
            list(method = "af", label = sprintf("%s_af_unique", mt), column = "mt_AF",  data = check_cols(ant_mt, "mt_AF"), 
                 plot_title = sprintf("log10(p-value) | AF : %s", mt), y_title = "-log10(p-value)", plotid = sprintf("sim_af_fishertest_pvalue_%s", mt)),
            list(method = "af", label = sprintf("%s_af_unique_ems", mt), column = "mt_AF", data = check_cols(ant_mt_ems, "mt_AF"), 
                 plot_title = sprintf("log10(p-value) | AF : %s", mt), y_title = "-log10(p-value)", plotid = sprintf("sim_af_fishertest_pvalue_ems_%s", mt))
          ))
        } else {
          params <- append(params,list(
            list(method = "af", label = sprintf("afd_%s_%s", wt, mt), column = "AFD", data = check_cols(wt_mt, "AFD"), 
                 plot_title = sprintf("log10(p-value) | AFD : %s & %s", wt, mt), y_title = "-log10(p-value)", 
                 plotid = sprintf("sim_afd_fishertest_pvalue_%s_vs_%s", wt, mt))
          ))
        }
      }
      
      if (stat == "ed" && !only_mutant) {
        params <- append(params, list(
          list(method = "ed", label = sprintf("ed4_%s_%s", wt, mt), column = "ED4", data = check_cols(wt_mt, "ED4"), 
               plot_title = sprintf("log10(p-value) | ED4 : %s vs %s", wt, mt), y_title = "-log10(p-value)", plotid = sprintf("sim_ed4_fisher_log_adj.pval_%s_%s", wt, mt))
        ))
      }
    }
    return(params)
  }
  parameters <- make_parameters()
  results_list <- list()
  
  for (pmt in parameters) {
    
    # 0. Skip empty data frames early
    if (!is.null(pmt$data) && nrow(pmt$data) > 0) {
      message("Running: ", pmt$label)
      res <- scan_it(pmt$data, pmt$column, window_size, step_size, pmt$method, af_doorstep, ed_doorstep)
      
      # add significance column
      if (!is.null(res$results) && nrow(res$results) > 0) {
        res$results[, sig := log.pval > threshold]
      } else {
        warning("No sliding-window results for ", pmt$label)
        next
      }
  
      # Plot
      if (plot_data) {
        plot_vcfdata(
          data = res$results, column = "log.pval", y_title = "-log10(p)", plot_title = pmt$y_title,
          file_suffix = pmt$plotid, threshold = threshold, ylim = ylim, is_smooth = FALSE,
          is_rollmedian = FALSE, is_histogram = FALSE, width = width, height = height,
          hwidth = hwidth, hheight = hheight, dpi = dpi, device = device, prefix = prefix,
          plots_dir = plots_dir, plot_style = plot_style, color_panel = color_panel)
      }
      results_list[[pmt$label]] <- res
    }
  }
  
  if (save_intervals) {
    if (save_excel) {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      excel_stub <- if (only_mutant || is.null(wt)) {
        sprintf("%s_%s", prefix, mt)
      } else {
        sprintf("%s_%s_%s", prefix, wt, mt)
      }

      excel_path <- file.path(output_dir, sprintf("%s_%s_%s_sim_pval_results.xlsx", , excel_stub))
      wb <- createWorkbook()
      
      for (label in names(results_list)) {
        res <- results_list[[label]]
        sheet1 <- substr(paste0(label, "_results"), 1, 31)
        sheet2 <- substr(paste0(label, "_thresh"), 1, 31)
        sheet3 <- substr(paste0(label, "_intervals"), 1, 31)
        addWorksheet(wb, sheet1); writeData(wb, sheet1, res$results)
        addWorksheet(wb, sheet2); writeData(wb, sheet2, res$abovethreshold)
        addWorksheet(wb, sheet3); writeData(wb, sheet3, res$intervals)
      }
      saveWorkbook(wb, excel_path, overwrite = TRUE)
      message("Excel saved at: ", excel_path)
    
  } else {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    for (label in names(results_list)) {
      res <- results_list[[label]]
      out_prefix <- file.path(output_dir, sprintf("%s.%s", prefix, label))
      fwrite(res$results, paste0(out_prefix, "_results.csv"))
      fwrite(res$abovethreshold, paste0(out_prefix, "_abovethreshold.csv"))
      fwrite(res$intervals, paste0(out_prefix, "_intervals.csv"))
    }
    message("CSVs saved in: ", output_dir)
    }
  }
  return(results_list)
}