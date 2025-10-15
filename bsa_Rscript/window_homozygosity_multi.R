#' Run window-based homozygosity analysis for multiple prefixes/genotype sets
#' @description
#' A wrapper around `window_homozygosity_pipeline()` to automatically run 
#' sliding-window homozygosity analysis across multiple prefixes and 
#' genotype combinations. Creates separate subfolders for each run and 
#' supports both mutant-only and WT vs MT comparisons.
#' @param vcf_dir Directory containing the VCF/TSV files.
#' @param pattern Regex pattern to match input files.
#' @param wt_list Character vector of wildtype genotype identifiers.
#' @param mt_list Character vector of mutant genotype identifiers.
#' @param prefix_list Character vector of prefixes to loop over.
#' @param output_dir Directory to save outputs (default = "post_analysis").
#' @param min_DP Minimum read depth filter.
#' @param min_QUAL Minimum quality filter.
#' @param only_mutant Logical; run in mutant-only mode.
#' @param use_ems Logical; include EMS-type SNPs.
#' @param save_results Logical; save intermediate analysis results.
#' @param save_interval Logical; save interval results to Excel.
#' @param rollmedian Window size for rolling median smoothing.
#' @param af_col Allele frequency column to use (default chooses automatically).
#' @param window_size Sliding window size (bp).
#' @param step_size Step size for sliding windows (bp).
#' @param nn_prop Proportion of neighbors for locfit smoothing.
#' @param af_min Allele frequency cutoff for homozygosity.
#' @param offhold Proportion threshold for peak cutoff.
#' @param min_vsize Minimum interval size (bp).
#' @param find_intervals Logical; whether to detect genomic peak intervals.
#' @param do_plot Logical; whether to generate plots.
#' @param bsa_metrics Character vector of BSA metrics to calculate.
#' @param plot_mode Which plot types ("hist","line","all").
#' @param use_col Columns to use for interval detection ("wmd","lft","rmd","all").
#' @param plots_dir Directory for plot outputs.
#' @param device Plot device (e.g., "png", "pdf").
#' @param dpi Resolution for plots.
#' @param hwidth,hheight Dimensions for wrapped plots.
#' @param width,height Dimensions for grid plots.
#' @param bwidth Histogram bin width (bp).
#' @param facet_column Number of columns for faceting.
#' @param line_size Line size for plots.
#' @param plot_style Faceting style ("wrap" or "grid").
#' @param color_panel Vector of colors for chromosomes.
#' @export
window_homozygosity_multi <- function(
    vcf_dir, pattern, wt_list, mt_list, prefix_list, output_dir = "post_analysis", plots_dir = "plots_dir",
    min_DP = 5, min_QUAL = 5, only_mutant = FALSE, use_ems = TRUE, save_results = FALSE, save_interval = TRUE,
    rollmedian = 501L, af_col = NULL, window_size = 5e5, step_size = 1e5, nn_prop = 0.1, af_min = 0.99, 
    offhold = 0.80, min_vsize = 0L, find_intervals = TRUE, do_plot = TRUE,
    bsa_metrics = c("waf","maf","homozygosity", "all"), plot_mode = c("hist","line","all"), 
    use_col = c("wmd","lft","rmd","all"), device = "png", dpi = 300, hwidth = 30, hheight = 18, width = 45, height = 15,
    bwidth = 1000000, facet_column = 5, line_size = 5, plot_style = c("wrap","grid"), color_panel = c("blue","red")
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
      sub_plots <- file.path(sub_plots, "window_homo_plots")
      sub_output <- file.path(output_dir, key)
      sub_output <- file.path(sub_output, "window_homo_analysis")
      dir.create(sub_plots,  showWarnings = FALSE, recursive = TRUE)
      dir.create(sub_output, showWarnings = FALSE, recursive = TRUE)
      
      message("Mutant-only: ", key)
      results[[key]] <-  window_homozygosity_single(
        vcf_dir = vcf_dir, prefix = pfx, pattern = pattern, Genotypes = list(mt = mt), output_dir = sub_output, 
        min_DP = min_DP, min_QUAL = min_QUAL, only_mutant = TRUE, use_ems = use_ems, 
        save_results = save_results, save_interval = save_interval,rollmedian = rollmedian, af_col = af_col, 
        window_size = window_size, step_size = step_size, nn_prop = nn_prop, af_min = af_min,
        find_intervals = find_intervals, offhold = offhold, min_vsize = min_vsize, bsa_metrics = bsa_metrics, 
        plot_mode = plot_mode, use_col = use_col, do_plot = do_plot, plots_dir = sub_plots,
        device = device, dpi = dpi, hwidth = hwidth, hheight = hheight, width = width, height = height,
        bwidth = bwidth, facet_column = facet_column, line_size = line_size, plot_style = plot_style_choice, color_panel = color_panel
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
      results[[key]] <- window_homozygosity_single(
        vcf_dir = vcf_dir, prefix = pfx, pattern = pattern, Genotypes = list(mt = mt, wt = wt), output_dir = sub_output, 
        min_DP = min_DP, min_QUAL = min_QUAL, only_mutant = FALSE, use_ems = use_ems, 
        save_results = save_results, save_interval = save_interval, rollmedian = rollmedian, af_col = af_col, 
        window_size = window_size, step_size = step_size, nn_prop = nn_prop, af_min = af_min,
        find_intervals = find_intervals, offhold = offhold, min_vsize = min_vsize, bsa_metrics = bsa_metrics, 
        plot_mode = plot_mode, use_col = use_col, do_plot = do_plot, plots_dir = sub_plots,
        device = device, dpi = dpi, hwidth = hwidth, hheight = hheight, width = width, height = height,
        bwidth = bwidth, facet_column = facet_column, line_size = line_size, plot_style = plot_style_choice, color_panel = color_panel
      )
    }
  }
  }
  return(results)
}
