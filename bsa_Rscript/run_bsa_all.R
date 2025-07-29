run_bsa_all <- function(
    vcf_dir, pattern,
    wt_list, mt_list, prefix_list,
    min_DP = 5, min_QUAL = 5, only_mutant = FALSE,
    plots_dir = "plots", output_dir = "post_analysis",
    save_results = FALSE, save_intervals = TRUE, plot_data = TRUE, save_excel = TRUE,
    rollmedian = 501, ylim = NULL, device = "png", 
    width = 45, height = 13, hwidth = 30, hheight = 18, dpi = 300,
    nn_prop = 0.1, plot_mode = "both", plot_style = "grid",
    af_min = 0.99, ed_min = 1e-5, threshold = -log10(0.05) * 10, bwidth = 1e6,
    plot_types = c("af", "gstat", "ed", "histogram", "pval"), window_size = 2000000, step_size = 1000000,
    stat_method = c("af", "ed"), af_doorstep = 0.67, ed_doorstep = 1, color_panel = c("blue", "red")
) {
  
  results <- list()
  seen_combos <- character()
  
  for (pfx in prefix_list) {
    message("\n==== Running BSA pipeline for prefix: ", pfx, " ====")
    # ensure output folders exist
    dir.create(plots_dir,   showWarnings=FALSE, recursive=TRUE)
    dir.create(output_dir,  showWarnings=FALSE, recursive=TRUE)
    
    if (only_mutant) {
      # mutant‐only mode: loop only over mt_list
      for (mt in mt_list) {
        combo_key <- paste(pfx, mt, sep = "_")
        if (combo_key %in% seen_combos) next
        seen_combos <- c(seen_combos, combo_key)
        sub_plots <- file.path(plots_dir, combo_key)
        dir.create(sub_plots, showWarnings=FALSE, recursive=TRUE)
        message("Running mutant‐only: ", combo_key)
        
        res <- run_bsa_pipeline(
          vcf_dir=vcf_dir, prefix=pfx, pattern=pattern,Genotypes=list(mt = mt),
          min_DP=min_DP, min_QUAL=min_QUAL, only_mutant=TRUE,
          plots_dir=sub_plots, output_dir=output_dir, save_results=save_results,
          save_intervals= save_intervals, plot_data=plot_data, save_excel=save_excel, rollmedian=rollmedian,
          ylim=ylim, device=device, width=width, height=height, hwidth=hwidth, hheight= hheight, dpi=dpi,
          nn_prop=nn_prop, plot_mode=plot_mode, plot_style=plot_style,
          af_min=af_min, ed_min=ed_min, threshold=threshold,bwidth=bwidth,plot_types=plot_types,
          window_size=window_size,step_size=step_size,stat_method=stat_method,
          af_doorstep=af_doorstep,ed_doorstep=ed_doorstep,color_panel=color_panel
        )
        results[[combo_key ]] <- res
      }
    } else {
      # WT vs MT mode: nested loops
      for (wt in wt_list) {
        for (mt in mt_list) {
        if (wt == mt) next
          pair_name <- paste(sort(c(wt, mt)), collapse = "_")
          combo_key <- paste(pfx, pair_name, sep = "_")
          
          if (combo_key %in% seen_combos) next
          seen_combos <- c(seen_combos, combo_key)
        
          sub_plots <- file.path(plots_dir, combo_key)
          dir.create(sub_plots, showWarnings=FALSE, recursive=TRUE)
        
          message("run:", combo_key, "  (WT = ", wt, ",  MT = ", mt, ")\n")
          
          res <- run_bsa_pipeline(
            vcf_dir=vcf_dir, prefix=pfx, pattern=pattern,Genotypes = list(wt = wt, mt = mt),
            min_DP=min_DP, min_QUAL=min_QUAL, only_mutant=FALSE,
            plots_dir=sub_plots, output_dir=output_dir, save_results=save_results,
            save_intervals= save_intervals, plot_data=plot_data, save_excel=save_excel, rollmedian=rollmedian,
            ylim=ylim, device=device, width=width, height=height, hwidth=hwidth, hheight= hheight, dpi=dpi,
            nn_prop=nn_prop, plot_mode=plot_mode, plot_style=plot_style,
            af_min=af_min, ed_min=ed_min, threshold=threshold,bwidth=bwidth,plot_types=plot_types,
            window_size=window_size,step_size=step_size,stat_method=stat_method,
            af_doorstep=af_doorstep,ed_doorstep=ed_doorstep,color_panel=color_panel
          )
          results[[combo_key]] <- res
          }
        }
      }
  }
  message("\n==== run_bsa_all finished ====")
  return(invisible(results))
}
