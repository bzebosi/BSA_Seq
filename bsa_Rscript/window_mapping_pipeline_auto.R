window_mapping_pipeline_auto <- function(
    vcf_dir, pattern, wt_list, mt_list, prefix_list, only_mutant = FALSE, use_ems = FALSE,
    bsa_metrics = c("waf","maf","ed","ed4","g","afd","all"),
    plots_dir = "window_bsa_plots", output_dir = "post_analysis",
    save_results = FALSE, save_interval = TRUE, do_plot = TRUE,
    device = "png", dpi = 300, plot_style = c("wrap","grid"),
    color_panel = c("blue","red"), hwidth = 30, hheight = 18, width = 45, height = 15,
    facet_column = 5, line_size = 5, window_size = 5e5, step_size = 1e5, rollmedian = 501L, nn_prop = 0.1,
    find_intervals = TRUE, offhold = 0.80, min_vsize = 0L,
    use_col = c("wmd","lft","rmd","all")
){
  results <- list()
  seen <- character()
  plot_style_choice <- match.arg(plot_style)
  
  dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (pfx in prefix_list) {
    message("\n==== Running BSA for prefix: ", pfx, " ====")
    
    if (only_mutant) {
      for (mt in mt_list) {
        key <- paste(pfx, mt, sep = "_")
        if (key %in% seen) next
        seen <- c(seen, key)
        
        sub_plots  <- file.path(plots_dir,  key)
        sub_output <- file.path(output_dir, key)
        dir.create(sub_plots,  showWarnings = FALSE, recursive = TRUE)
        dir.create(sub_output, showWarnings = FALSE, recursive = TRUE)
        
        message("Mutant-only: ", key)
        results[[key]] <- window_mapping_pipeline(
          vcf_dir = vcf_dir, prefix = pfx, pattern = pattern, Genotypes = list(mt = mt),
          output_dir = sub_output, only_mutant = TRUE,
          use_ems = use_ems, save_results = save_results, save_interval = save_interval,
          rollmedian = rollmedian, window_size = window_size, step_size = step_size,
          nn_prop = nn_prop, find_intervals = find_intervals, offhold = offhold, min_vsize = min_vsize,
          use_col = use_col, do_plot = do_plot, plots_dir = sub_plots,
          device = device, dpi = dpi, hwidth = hwidth, hheight = hheight, width = width, height = height,
          facet_column = facet_column, line_size = line_size,
          bsa_metrics = bsa_metrics, plot_style = plot_style_choice, color_panel = color_panel
        )
      }
      
    } else {
      for (wt in wt_list) for (mt in mt_list) {
        if (wt == mt) next
        pair <- paste(sort(c(wt, mt)), collapse = "_")
        key  <- paste(pfx, pair, sep = "_")
        if (key %in% seen) next
        seen <- c(seen, key)
        
        sub_plots  <- file.path(plots_dir,  key)
        sub_output <- file.path(output_dir, key)
        dir.create(sub_plots,  showWarnings = FALSE, recursive = TRUE)
        dir.create(sub_output, showWarnings = FALSE, recursive = TRUE)
        
        message("run: ", key, " (wt=", wt, ", mt=", mt, ")")
        results[[key]] <- window_mapping_pipeline(
          vcf_dir = vcf_dir, prefix = pfx, pattern = pattern,
          Genotypes = list(wt = wt, mt = mt),
          output_dir = sub_output, only_mutant = FALSE,
          use_ems = use_ems, save_results = save_results, save_interval = save_interval,
          rollmedian = rollmedian, window_size = window_size, step_size = step_size,
          nn_prop = nn_prop, find_intervals = find_intervals, offhold = offhold, min_vsize = min_vsize,
          use_col = use_col, do_plot = do_plot, plots_dir = sub_plots,
          device = device, dpi = dpi, hwidth = hwidth, hheight = hheight, width = width, height = height,
          facet_column = facet_column, line_size = line_size,
          bsa_metrics = bsa_metrics, plot_style = plot_style_choice, color_panel = color_panel
        )
      }
    }
  }
  
  message("\n==== window_bsa_auto finished ====")
  invisible(results)
}
