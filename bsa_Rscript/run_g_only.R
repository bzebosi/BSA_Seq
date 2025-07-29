#' Generate G-statistic Plots for BSA-Seq Data
#'
#' Plots G-statistics from merged wild-type and mutant SNP tables to visualize genomic intervals 
#' with significant allele frequency divergence. Supports locfit and rollmedian smoothing. 
#' Skips plotting in mutant-only mode.
#'
#' @param wt_mt A data.table with merged wild-type and mutant SNPs containing the \code{G} column.
#' @param wt Wild-type label for plot titles.
#' @param mt Mutant label for plot titles.
#' @param prefix Prefix for output plot file names.
#' @param plots_dir Directory to save output plots.
#' @param rollmedian Window size for rolling median smoothing (used if \code{plot_mode} is \code{"rollmedian"} or \code{"both"}).
#' @param ylim Optional y-axis limits.
#' @param only_mutant Logical. If \code{TRUE}, G-statistic plots are skipped.
#' @param device Graphics device to use (e.g., \code{"png"}, \code{"pdf"}).
#' @param plot_data Logical. If \code{TRUE}, plots will be generated.
#' @param width Width of the plot.
#' @param height Height of the plot.
#' @param hwidth Width of the histogram panel.
#' @param hheight Height of the histogram panel.
#' @param dpi Resolution of the output image.
#' @param nn_prop Proportion of nearest neighbors for locfit smoothing.
#' @param g_min Minimum G-statistic value to include SNPs for plotting.
#' @param plot_mode Plotting method. One of \code{"locfit"}, \code{"rollmedian"}, or \code{"both"}.
#' @param plot_style Layout style for plots. Typically \code{"wrap"} or \code{"grid"}.
#' @param color_panel Color vector for plot aesthetics.
#' @return NULL. Plots are saved to the specified directory.
#' @examples
#' ' \dontrun{
#' run_g_only(wt_mt = merged_snps, wt = "WT", mt = "Ts5", prefix = "b73_Ts5",
#' plots_dir = "plots/", plot_mode = "both", g_min = 1e-4)
#' }
#' @export
run_g_only <- function(
    wt_mt = NULL, wt = "wildtype", mt = "mutant", prefix = "sample", plots_dir = "plots", 
    rollmedian = 501, ylim = NULL, only_mutant = FALSE, device = "png",plot_data = TRUE,
    width = 45, height = 12, hwidth = 30, hheight = 18, dpi = 300, nn_prop = 0.1, g_min = 1e-5,
    plot_mode = "both", plot_style = "wrap", color_panel = c("blue", "red")){
  
  
  if (only_mutant) {
    message("Skipping ed plots in mutant-only mode.")
    return(NULL)
  }
  
  # fix: include g_min inside check_cols
  check_cols <- function(data, colname, g_min) {
    if (is.null(data)) return(NULL)
    if (!(colname %in% colnames(data))) return(NULL)
    data <- data[data[[colname]] >= g_min, ]
    if (nrow(data) == 0) return(NULL)
    return(data)
  }
  
  make_parameters <- function() {
    list(
      list(column = "G", data = check_cols(wt_mt, "G", g_min), plot_title = sprintf("G : %s & %s", wt, mt), y_title = "G", plotid = sprintf("g_%s_%s", wt, mt))
    )
  }
  
  plot_with_mode <- function(mode, data, column, plot_title, y_title, plotid) {
    if (mode == "locfit" || mode == "both") {
      plot_vcfdata(
        data = data, column = column, 
        y_title = y_title, plot_title = plot_title, file_suffix = sprintf("locfit_%s", plotid),
        ylim = ylim, is_smooth = TRUE, is_rollmedian = FALSE, rollmedian = 0, is_histogram = FALSE, 
        nn_prop = nn_prop, width = width, height = height, hwidth = hwidth, hheight = hheight, dpi = dpi, device = device, 
        prefix = prefix, plots_dir = plots_dir, plot_style = plot_style, color_panel = color_panel
      )
    }
    
    if ((mode == "rollmedian" || mode == "both") && rollmedian > 0) {
      plot_vcfdata(
        data = data, column = column, 
        y_title = y_title, plot_title = plot_title, file_suffix = sprintf("g_rollmedian_%s", plotid),
        ylim = ylim, is_smooth = FALSE, is_rollmedian = TRUE, rollmedian = rollmedian, is_histogram = FALSE, 
        nn_prop = nn_prop, width = width, height = height, hwidth = hwidth, hheight = hheight, dpi = dpi, device = device, 
        prefix = prefix, plots_dir = plots_dir, plot_style = plot_style, color_panel = color_panel
      )
    }
  }
  
  # Validate plot mode
  valid_modes <- c("locfit", "rollmedian", "both")
  if (!(plot_mode %in% valid_modes)) {
    stop(sprintf("Invalid plot_mode: '%s'. Choose one of: %s", plot_mode, paste(valid_modes, collapse = ", ")))
  }
  
  # Generate plots
  parameters <- make_parameters()
  if (plot_data) {
    for (pmt in parameters) {
      if (!is.null(pmt$data)) {
        message(sprintf("Generating G-statistic plot: %s", pmt$plot_title))
        plot_with_mode(plot_mode, pmt$data, pmt$column, pmt$plot_title, pmt$y_title, pmt$plotid)
      }
    }
  }
  
  message("G-statistic plots completed.")
}


