window_bsa_pipeline <- function(vcf_dir, prefix, pattern, Genotypes = list(wt = "wildtype", mt = "mutant"),
                                         output_dir = "post_analysis", min_DP = 5, min_QUAL = 5, only_mutant = FALSE, interval_mutant = FALSE,
                                         use_ems = FALSE, save_results = FALSE, save_interval = TRUE, rollmedian = 501L,
                                         window_size = 5e5, step_size = 1e5, nn_prop = 0.1,
                                         find_intervals = TRUE, offhold = 0.80, min_vsize = 0L, 
                                         use_cols= c("wmd","lft","rmd","all"),
                                         do_plot = TRUE, plots_dir = "window_bsa_plots", device = "png", dpi = 300, 
                                         hwidth = 30, hheight = 18, width = 45, height =15, facet_column = 5, line_size = 5,
                                         bsa_metrics = c("waf","maf","ed","ed4","g","afd", "all"), 
                                         plot_style = c("wrap","grid"), color_panel = c("blue","red")){
  
  
  
  
  
  message("=== Step 1: Importing VCF Data ===")
  geno_data <- import_vcfdata(vcf_dir = vcf_dir, prefix = prefix, pattern = pattern,
                              Genotypes = Genotypes, min_DP = min_DP, min_QUAL = min_QUAL, only_mutant = only_mutant)
  message("=== Step 2: Analyzing VCF Data ===")
  result <- analyze_vcfdata(geno_data = geno_data, prefix = prefix,
                            save_results = save_results, output_dir = output_dir, only_mutant = only_mutant)
  
  # message("=== Step 3: Sliding-window homozygosity ====")

  message("=== Step 3: Sliding-window BSA metrics ===")
  chosen_metrics <- if (is.null(bsa_metrics)) c("waf","maf","ed","ed4","afd","g") else bsa_metrics
  win <- window_bsa_compute_all(
    data = result,
    bsa_metrics = chosen_metrics,
    use_ems = use_ems,only_mutant = only_mutant, window_size = window_size,
    step_size = step_size, rollmedian = rollmedian, nn_prop = nn_prop,
    find_intervals = find_intervals, offhold = offhold, min_vsize = min_vsize
  )
  
  # ---- Simple save: one Excel with interval sheets ----
  # ---- Save intervals to Excel (one sheet per metric) ----
  if (isTRUE(save_interval) && isTRUE(find_intervals)) {
    message("=== Step 4: Saving interval excel ====")
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      warning("save_interval=TRUE but package 'openxlsx' is not installed. Skipping Excel write.")
    } else {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      wt_name <- if (!is.null(Genotypes$wt)) Genotypes$wt else "wt"
      mt_name <- if (!is.null(Genotypes$mt)) Genotypes$mt else "mt"
      
      file_name <- if (isTRUE(interval_mutant)) {
        sprintf("%s_%s_intervals.xlsx", prefix, mt_name)
      } else {
        sprintf("%s_%s_vs_%s_bsa_intervals.xlsx", prefix, wt_name, mt_name)
      }
      
      # If intervals list exists and has any non-empty data frames, write them
      has_any <- FALSE
      wb <- openxlsx::createWorkbook()
      if (!is.null(win$intervals) && length(win$intervals)) {
        for (label in names(win$intervals)) {
          iv <- win$intervals[[label]]
          if (is.null(iv)) next
          if (is.data.frame(iv) && nrow(iv) > 0L) {
            has_any <- TRUE
            sheet <- substr(paste0(label, "_intervals"), 1, 31)  # Excel sheet name limit
            openxlsx::addWorksheet(wb, sheet)
            openxlsx::writeData(wb, sheet, iv)
          }
        }
      }
      
      xlsx_path <- file.path(output_dir, file_name)
      if (isTRUE(has_any)) {
        openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
        message("Saved intervals to: ", xlsx_path)
      } else {
        message("No interval rows to save (all empty).")
      }
    }
  }

  if (isTRUE(do_plot)) {
    message("=== Step 5: Plotting ===")
    if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
    
    # Labels for plot titles
    mt_label <- if (!is.null(Genotypes$mt)) Genotypes$mt else "mutant"
    wt_label <- if (!is.null(Genotypes$wt)) Genotypes$wt else "wildtype"
    
    # Match a single plot_style choice
    plot_style_choice <- match.arg(plot_style)
    
    # window_bsa_plot expects 'use_col' (singular); map your 'use_cols'
    use_col <- if (identical(use_cols, "all")) "all" else use_cols
    
    window_bsa_plot(data=win, prefix = prefix, mt = mt_label, wt = wt_label, only_mutant = only_mutant, 
                    use_ems = use_ems, bsa_metrics = chosen_metrics, use_col = use_col, 
                    hwidth = hwidth, hheight = hheight, width = width, height = height, plots_dir = plots_dir, 
                    facet_column = facet_column, line_size = line_size, plot_style = plot_style_choice, color_panel = color_panel)
  }
  
  return(list(geno_data = geno_data, results = result, out = win))
}
