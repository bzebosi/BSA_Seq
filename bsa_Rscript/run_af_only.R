#' Generate Allele Frequency (AF) and AFD Plots for BSA-Seq Data
#'
#' Creates SNP-index and delta SNP-index (AFD) plots from BSA-Seq data for wild-type and mutant samples.
#' Supports plotting of unique and EMS-specific SNPs, with options for locfit or rollmedian smoothing.
#'
#' @param wt_mt Merged wild-type and mutant SNP table containing the \code{AFD} column.
#' @param ant_wt Unique wild-type SNPs.
#' @param ant_mt Unique mutant SNPs.
#' @param ant_wt_ems EMS-type wild-type SNPs.
#' @param ant_mt_ems EMS-type mutant SNPs.
#' @param wt Label for wild-type sample (used in titles and filenames).
#' @param mt Label for mutant sample (used in titles and filenames).
#' @param prefix Prefix for output filenames.
#' @param plots_dir Directory where plots will be saved.
#' @param rollmedian Rolling median window size (used when \code{plot_mode} includes \code{"rollmedian"}).
#' @param ylim Optional y-axis limits.
#' @param only_mutant Logical. If \code{TRUE}, plots only mutant SNP-index values.
#' @param device Graphics device for saving plots (e.g., \code{"png"}, \code{"pdf"}).
#' @param plot_data Logical. If \code{TRUE}, plots are generated.
#' @param width Width of output plot (in inches).
#' @param height Height of main panel (in inches).
#' @param hwidth Width of histogram or secondary panel (in inches).
#' @param hheight Height of histogram or secondary panel (in inches).
#' @param dpi Plot resolution in dots per inch.
#' @param nn_prop Proportion of data used for locfit smoothing (between 0 and 1).
#' @param plot_mode Smoothing method: one of \code{"locfit"}, \code{"rollmedian"}, or \code{"both"}.
#' @param plot_style Layout of facet panels, e.g., \code{"wrap"} or \code{"grid"}.
#' @param color_panel Color vector for plot elements (e.g., c("blue", "red")).
#' @return NULL. Plots are saved to the specified directory.
#' @examples
#' \dontrun{
#' run_af_only(
#'   wt_mt = merged_data, ant_wt = unique_wt, ant_mt = unique_mt, ant_wt_ems = ems_wt, ant_mt_ems = ems_mt,
#'   wt = "WT", mt = "Ts5", prefix = "b73", plots_dir = "plots", plot_mode = "both")
#' }
#' @export


run_af_only <- function(
    wt_mt = NULL, ant_wt = NULL, ant_mt = NULL, ant_wt_ems = NULL, ant_mt_ems = NULL,
    wt = "wildtype", mt = "mutant", prefix = "sample", plots_dir = "plots", 
    rollmedian = 501, ylim = NULL, only_mutant = FALSE, device = "png", plot_data = TRUE,
    width = 45, height = 13, hwidth = 30, hheight = 18, dpi = 300, nn_prop = 0.1,
    plot_mode = "both", plot_style = "wrap", color_panel = c("blue", "red")) {
  
  
  check_cols <- function(data, colname) {
    if (is.null(data)) return(NULL)
    if (!(colname %in% colnames(data))) return(NULL)
    return(data)
  }
  
  make_parameters <- function() {
    if (only_mutant) {
      list(
        list(column = "mt_AF",  data = check_cols(ant_mt, "mt_AF"), plot_title = sprintf("%s SNP-Index (Unique SNPs)", mt), y_title = "AF", plotid = sprintf("%s_AF_unique", mt)),
        list(column = "mt_AF", data = check_cols(ant_mt_ems, "mt_AF"), plot_title = sprintf("%s SNP-Index (Unique ems SNPs)", mt), y_title = "AF", plotid = sprintf("%s_AF_unique_ems", mt))
      )
    } else {
      list(
        list(column = "AFD", data = check_cols(wt_mt, "AFD"), plot_title = sprintf("AFD (Δ SNP-Index): %s - %s", wt, mt), y_title = "AFD", plotid = sprintf("%s_%s_AFD", wt, mt)),
        list(column = "wt_AF", data = check_cols(ant_wt, "wt_AF"), plot_title = sprintf("%s SNP-Index (Unique SNPs)", wt), y_title =  "AF", plotid = sprintf("%s_AF_unique", wt)),
        list(column = "mt_AF",  data = check_cols(ant_mt, "mt_AF"), plot_title = sprintf("%s SNP-Index (Unique SNPs)", mt), y_title = "AF", plotid = sprintf("%s_AF_unique", mt)),
        list(column = "wt_AF", data = check_cols(ant_wt_ems, "wt_AF"), plot_title = sprintf("%s SNP-Index (Unique ems SNPs)", wt), y_title = "AF", plotid = sprintf("%s_AF_unique_ems", wt)),
        list(column = "mt_AF", data = check_cols(ant_mt_ems, "mt_AF"), plot_title = sprintf("%s SNP-Index (Unique ems SNPs)", mt), y_title = "AF", plotid = sprintf("%s_AF_unique_ems", mt))
      )
    }
  }
  
  plot_with_mode <- function(mode, data, column, plot_title, y_title, plotid) {
    if (mode == "locfit" || mode == "both") {
      plot_vcfdata(
        data = data, column = column, 
        y_title = y_title, plot_title = plot_title, file_suffix = sprintf("af_locfit_%s", plotid),
        ylim = ylim, is_smooth = TRUE, is_rollmedian = FALSE, rollmedian = 0, is_histogram = FALSE, 
        nn_prop = nn_prop, width = width, height = height, hwidth = hwidth, hheight = hheight, dpi = dpi, device = device, 
        prefix = prefix, plots_dir = plots_dir, plot_style = plot_style, color_panel = color_panel
      )
    }
    
    if ((mode == "rollmedian" || mode == "both") && rollmedian > 0) {
      plot_vcfdata(
        data = data, column = column, 
        y_title = y_title, plot_title = plot_title, file_suffix = sprintf("af_rollmedian_%s", plotid),
        ylim = ylim, is_smooth = FALSE, is_rollmedian = TRUE, rollmedian = rollmedian, is_histogram = FALSE, 
        nn_prop = nn_prop, width = width, height = height, hwidth = hwidth, hheight = hheight, dpi = dpi, device = device, 
        prefix = prefix, plots_dir = plots_dir, plot_style = plot_style, color_panel = color_panel
      )
    }
  }
  
  parameters <- make_parameters()
  valid_modes <- c("locfit", "rollmedian", "both")
  if (!(plot_mode %in% valid_modes)) {
    stop(sprintf("Invalid plot_mode: '%s'. Choose one of: %s", plot_mode, paste(valid_modes, collapse = ", ")))
  }
  
  # Generate configurations for plots
  if (plot_data) {
    for (pmt in parameters) {
      if (!is.null(pmt$data) && nrow(pmt$data) > 0) {
        message(sprintf("Generating AF plot: %s", pmt$plot_title))
        plot_with_mode(plot_mode, pmt$data, pmt$column, pmt$plot_title, pmt$y_title, pmt$plotid)
      }
    }
  }
  
  message("AF-only plots completed.")
}



