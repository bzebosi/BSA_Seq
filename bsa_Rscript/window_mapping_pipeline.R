window_mapping_pipeline <- function(vcf_dir, prefix, pattern, Genotypes = list(wt = "wildtype", mt = "mutant"),
                                         output_dir = "post_analysis", min_DP = 5, min_QUAL = 5, only_mutant = FALSE,
                                         use_ems = FALSE, save_results = FALSE, save_interval = TRUE, rollmedian = 501L,
                                         window_size = 5e5, step_size = 1e5, nn_prop = 0.1,
                                         find_intervals = TRUE, offhold = 0.80, min_vsize = 0L, 
                                         use_col= c("wmd","lft","rmd","all"),
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
  all_metrics <- c("waf","maf","ed","ed4","afd","g")
  chosen_metrics <- tolower(unique(bsa_metrics))
  if ("all" %in% chosen_metrics) chosen_metrics <- all_metrics
  chosen_metrics <- intersect(chosen_metrics, all_metrics)
  

  win <- window_mapping_compute_auto(
    data = result, bsa_metrics = chosen_metrics, use_ems = use_ems,only_mutant = only_mutant, 
    window_size = window_size, step_size = step_size, rollmedian = rollmedian, nn_prop = nn_prop,
    find_intervals = find_intervals, offhold = offhold, min_vsize = min_vsize
  )
  
  # ---- Simple save: one Excel with interval sheets ----
  if (isTRUE(save_interval) && isTRUE(find_intervals)) {
    message("=== Step 4: Saving interval excel ===")
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      wt_name <- if (!is.null(Genotypes$wt)) Genotypes$wt else "wt"
      mt_name <- if (!is.null(Genotypes$mt)) Genotypes$mt else "mt"
      
      # Add suffix depending on EMS usage
      ems_suffix <- if (use_ems) "ems" else "all"
      file_name <- if (only_mutant) {
        sprintf("%s_%s_window_bsa_intervals_%s.xlsx", prefix, mt_name, ems_suffix)
      } else {
        sprintf("%s_%s_vs_%s_window_bsa_intervals_%s.xlsx", prefix, wt_name, mt_name, ems_suffix)
      }
      
      wb <- openxlsx::createWorkbook()
      labs <- names(win$intervals)
      if (isTRUE(only_mutant)) labs <- labs[grepl("^mt_", labs)]
      for (label in labs) {
        iv <- win$intervals[[label]]
        if (!is.null(iv) && is.data.frame(iv) && nrow(iv) > 0L) {
          sheet <- substr(paste0(label, "_intervals"), 1, 31)
          openxlsx::addWorksheet(wb, sheet)
          openxlsx::writeData(wb, sheet, iv)
        }
      }
      xlsx_path <- file.path(output_dir, file_name)
      if (length(wb$worksheets) > 0L) {
        openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
        message("Saved intervals to: ", xlsx_path)
      } else {
        message("No interval rows to save (all empty).")
      }
    } else {
      warning("save_interval=TRUE but 'openxlsx' is not installed. Skipping Excel write.")
    }
  }
  

  if (isTRUE(do_plot)) {
    message("=== Step 5: Plotting ===")
    if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
    mt_label <- if (!is.null(Genotypes$mt)) Genotypes$mt else "mutant"
    wt_label <- if (!is.null(Genotypes$wt)) Genotypes$wt else "wildtype"
    
    # Match a single plot_style choice
    plot_style_choice <- match.arg(plot_style)
    
    # window_bsa_plot expects 'use_col' (singular); map your 'use_col'
    use_col <- if (identical(use_col, "all")) "all" else use_col
    
    window_mapping_plot(data=win, prefix = prefix, mt = mt_label, wt = wt_label, only_mutant = only_mutant, 
                    use_ems = use_ems, bsa_metrics = chosen_metrics, use_col = use_col, 
                    hwidth = hwidth, hheight = hheight, width = width, height = height, plots_dir = plots_dir, 
                    facet_column = facet_column, line_size = line_size, plot_style = plot_style_choice, color_panel = color_panel)
  }
  
  return(list(geno_data = geno_data, results = result, out = win))
}
