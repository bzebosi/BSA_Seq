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
  if ("ed" %in% stat_method && only_mutant) {
    stop("ed sliding window analysis requires both wildtype and mutant data (only_mutant = FALSE).")
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
          list(method = "ed", label = sprintf("ed_%s_%s", wt, mt), column = "ED", data = check_cols(wt_mt, "ED"), 
               plot_title = sprintf("log10(p-value) | ED : %s vs %s", wt, mt), y_title = "-log10(p-value)", plotid = sprintf("sim_ed_fisher_log_adj.pval_%s_%s", wt, mt)),
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
      excel_path <- file.path(output_dir, sprintf("%s_%s_%s_sim_pval_results.xlsx", prefix, wt, mt))
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