#' Run Bulk Segregant Analysis (BSA) Full Pipeline
#'
#' Wrapper function to import VCF SNP data, compute SNP statistics, and generate multiple
#' plot types for bulk segregant analysis (AFD, G-statistics, ED, histograms, Fisher test p-values).
#' Supports mutant-only and wildtype-mutant comparison designs.
#'
#' @param vcf_dir     Folder with input VCF files.
#' @param prefix      Output prefix used for saving results.
#' @param pattern     Regex pattern to match VCF filenames.
#' @param Genotypes   List with named entries \code{wt} and \code{mt}.
#' @param min_DP,min_QUAL Minimum depth and quality filter.
#' @param only_mutant Whether to run mutant-only analysis.
#' @param plots_dir,output_dir Paths to save plots and results.
#' @param save_results,save_intervals,save_excel Save control flags.
#' @param plot_data Whether to generate plots.
#' @param rollmedian Window size for rolling median smoothing.
#' @param ylim Y-axis limits for all plots.
#' @param device Output plot format (e.g., \code{"png"}).
#' @param width,height,hwidth,hheight,dpi Plot dimensions and resolution.
#' @param nn_prop Locfit smoothing parameter.
#' @param plot_mode Either \code{"rollmedian"}, \code{"locfit"}, or \code{"both"}.
#' @param plot_style Either \code{"grid"} or \code{"wrap"}.
#' @param af_min,ed_min Thresholds for homozygosity plots.
#' @param threshold Y cutoff for p-value plots.
#' @param bwidth Bin width for histograms.
#' @param plot_types Which plot types to run (AF, ED, G-stat, histogram, pval).
#' @param window_size,step_size Sliding window size and step (bp).
#' @param stat_method Statistic(s) to use for Fisher p-values: \code{"af"}, \code{"ed"}, or both.
#' @param af_doorstep,ed_doorstep Cutoffs used for Fisher window counts.
#' @param color_panel  Vector of colors for chromosomes.
#' @return A named list with:
#' \item{geno_data}{VCF SNP input}
#' \item{analysis_result}{Output from \code{analyze_vcfdata()}}
#' \item{sim_pval_results}{Sliding-window p-value results}
#'
#' @examples
#' \dontrun{
#' run_bsa_pipeline(
#'   vcf_dir = "vcf_files", prefix = "b73_ts1", pattern = "ts1",
#'   Genotypes = list(wt = "WT", mt = "ts1"),
#'   only_mutant = FALSE, plot_types = c("af", "pval", "histogram"),
#'   threshold = -log10(0.05) * 10, stat_method = "af"
#' )
#' }
#'
#' @export
run_bsa_pipeline <- function(
    vcf_dir, prefix, pattern,
    Genotypes = list(wt = "wildtype", mt = "mutant"),
    min_DP = 5, min_QUAL = 5, only_mutant = FALSE,
    plots_dir = "plots", output_dir = "post_analysis",
    save_results = FALSE, save_intervals = TRUE, plot_data = TRUE, save_excel = TRUE,
    rollmedian = 501, ylim = NULL, device = "png", 
    width = 45, height = 13, hwidth = 30, hheight = 18, dpi = 300,
    nn_prop = 0.1, plot_mode = "both", plot_style = "grid",
    af_min = 0.99, ed_min = 1e-5, threshold = -log10(0.05) * 10, bwidth = 1e6,
    plot_types = c("af", "gstat", "ed", "histogram", "pval"),
    window_size = 2000000, step_size = 1000000,
    stat_method = c("af", "ed"), af_doorstep = 0.67, ed_doorstep = 1,
    color_panel = c("blue", "red")
) {
  message("=== Step 1: Importing VCF Data ===")
  geno_data <- import_vcfdata(vcf_dir = vcf_dir, prefix = prefix, pattern = pattern,
    Genotypes = Genotypes, min_DP = min_DP, min_QUAL = min_QUAL, only_mutant = only_mutant)
  
  message("=== Step 2: Analyzing VCF Data ===")
  result <- analyze_vcfdata(geno_data = geno_data, prefix = prefix,
    save_results = save_results,output_dir = output_dir, only_mutant = only_mutant)
  
  # Extract SNP datasets from results
  wt_mt <- result$wt_mt
  ant_wt <- result$ant_wt
  ant_mt <- result$ant_mt
  ant_wt_ems <- result$ant_wt_ems
  ant_mt_ems <- result$ant_mt_ems
  
  message("=== Step 3: Running All BSA Plots ===")
  # Normalize input for case insensitivity
  plot_types <- tolower(plot_types)
  
  # Define valid plot types
  valid_plot_types <- c("af", "gstat", "ed", "histogram", "pval")
  
  # Check for invalid plot types
  if (!all(plot_types %in% valid_plot_types)) {
    invalid_types <- plot_types[!plot_types %in% valid_plot_types]
    warning(paste("Invalid plot types detected:", paste(invalid_types, collapse = ", "), "- Skipping these."))
    plot_types <- plot_types[plot_types %in% valid_plot_types]
  }
  
  
  sim_pval_results <- NULL
  if ("pval" %in% plot_types) {
    message("==== Running Sliding Window Fisher p-value Plots ====")
    sim_pval_results <- run_simpval_only(wt_mt, ant_wt, ant_mt, ant_wt_ems, ant_mt_ems,
      wt=if(only_mutant) NULL else Genotypes$wt, mt = Genotypes$mt,
      prefix = prefix, plots_dir = plots_dir, output_dir = output_dir, ylim = ylim, only_mutant = only_mutant, 
      device = device, width = width, height = height, hwidth = hwidth, hheight = hheight, dpi = dpi, 
      plot_style = plot_style, threshold = threshold, save_intervals = save_intervals, save_excel = save_excel, plot_data = plot_data,
      window_size = window_size, step_size = step_size,stat_method = stat_method, af_doorstep = af_doorstep, ed_doorstep = ed_doorstep,
      color_panel = color_panel
    )
  }
  
  if ("af" %in% plot_types) {
    message("==== Generating AF plots ====")
    run_af_only(
      wt_mt, ant_wt, ant_mt, ant_wt_ems, ant_mt_ems,
      wt=if(only_mutant) NULL else Genotypes$wt, mt=Genotypes$mt, prefix=prefix,
      plots_dir=plots_dir, rollmedian=rollmedian, nn_prop=nn_prop, ylim=ylim, only_mutant=only_mutant,
      plot_data=plot_data, plot_mode=plot_mode, plot_style=plot_style,
      width=width, height=height, hwidth=hwidth, hheight=hheight, device=device, dpi=dpi,color_panel=color_panel
    )
  }
  
  if ("ed" %in% plot_types) {
    message("==== Generating ED plots ====")
    run_ed_only(
      wt_mt, wt=if(only_mutant) NULL else Genotypes$wt, mt=Genotypes$mt, prefix=prefix,
      plots_dir=plots_dir, rollmedian=rollmedian, nn_prop=nn_prop, ylim=ylim, only_mutant=only_mutant,
      plot_data=plot_data, plot_mode=plot_mode, plot_style=plot_style,color_panel=color_panel, 
      width=width, height=height, hwidth=hwidth, hheight=hheight, device=device, dpi=dpi)
  }
  
  if ("gstat" %in% plot_types) {
    message("==== Generating gstatistics plots ====")
    run_g_only(
      wt_mt, wt=if(only_mutant) NULL else Genotypes$wt, mt=Genotypes$mt, prefix=prefix,
      plots_dir = plots_dir, rollmedian = rollmedian, ylim = ylim, only_mutant = only_mutant, 
      plot_data = plot_data, plot_style = plot_style, plot_mode = plot_mode, color_panel=color_panel, 
      width=width, height=height, hwidth=hwidth, hheight=hheight, device=device, dpi=dpi)
  }
  
  
  if ("histogram" %in% plot_types) {
    message("==== Generating homozygosity plots ====")
    run_histogram_only(
    wt_mt, ant_wt, ant_mt, ant_wt_ems, ant_mt_ems,
    wt=if(only_mutant) NULL else Genotypes$wt, mt=Genotypes$mt, prefix=prefix,
    plots_dir = plots_dir, ylim = ylim, only_mutant = only_mutant,
    plot_data = plot_data, plot_style = plot_style, color_panel=color_panel,
    is_histogram = TRUE, af_min = af_min, bwidth = bwidth,
    width=width, height=height, hwidth=hwidth, hheight=hheight, device=device, dpi=dpi)
  }
  

  message("=== BSA Pipeline Complete ===")
  
  return(list(
    geno_data = geno_data,
    analysis_result = result,
    sim_pval_results = sim_pval_results
  ))
}
