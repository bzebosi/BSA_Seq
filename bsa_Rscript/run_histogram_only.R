#' Generate Histogram of Unique SNPs for BSA-Seq Data
#' Plots histograms of unique SNP allele frequencies from wild-type and mutant SNP tables.
#' Highlights EMS-type variants and allows filtering by allele frequency threshold.
#' Supports both standard and mutant-only BSA-Seq designs.
#' @param wt_mt Merged SNP table (optional, not used in histogram but included for consistency).
#' @param ant_wt Unique wild-type SNPs.
#' @param ant_mt Unique mutant SNPs.
#' @param ant_wt_ems EMS-type wild-type SNPs.
#' @param ant_mt_ems EMS-type mutant SNPs.
#' @param wt Wild-type label for plot titles.
#' @param mt Mutant label for plot titles.
#' @param prefix Prefix for output filenames.
#' @param plots_dir Directory to save output plots.
#' @param ylim Optional y-axis limits.
#' @param only_mutant Logical. If \code{TRUE}, only mutant-specific plots will be drawn.
#' @param device Graphics device to use (e.g., \code{"png"}, \code{"pdf"}).
#' @param is_histogram Logical. Must be \code{TRUE} to enable histogram plotting.
#' @param plot_data Logical. If \code{TRUE}, plots will be generated.
#' @param width Width of the plot (in inches).
#' @param height Height of the plot (in inches).
#' @param hwidth Width of the histogram panel (in inches).
#' @param hheight Height of the histogram panel (in inches).
#' @param dpi Resolution of the output image.
#' @param plot_style Faceting layout. One of \code{"wrap"} or \code{"grid"}.
#' @param af_min Minimum allele frequency required to include SNPs.
#' @param bwidth Bin width in base pairs (e.g., 1e6 for 1 Mb bins).
#' @param color_panel Color vector used for plotting.
#' @return NULL. Plots are saved to the specified directory.
#' @examples
#' \dontrun{
#' run_histogram_only(
#'   ant_wt = unique_wt_snps, ant_mt = unique_mt_snps, ant_wt_ems = unique_wt_ems,
#'   ant_mt_ems = unique_mt_ems, wt = "WT", mt = "Ts5", prefix = "b73",
#'   plots_dir = "plots/", af_min = 0.99, bwidth = 1e6)
#' }
#' @export
run_histogram_only <- function(
    wt_mt = NULL, ant_wt = NULL, ant_mt = NULL, ant_wt_ems = NULL, ant_mt_ems = NULL,
    wt = "wildtype", mt = "mutant", prefix = "sample", plots_dir = "plots", 
    ylim = NULL, only_mutant = FALSE, device = "png", is_histogram = TRUE, plot_data = TRUE,
    width = 45, height = 15, hwidth = 30, hheight = 18, dpi = 300, 
    plot_style = "wrap", af_min = 0.99, bwidth = 1e6, color_panel = c("blue", "red")) {
  

  # fix: include g_min inside check_cols
  check_cols <- function(data, colname, af_min) {
    if (is.null(data)) return(NULL)
    if (!(colname %in% colnames(data))) return(NULL)
    data <- data[data[[colname]] >= af_min, ]
    if (nrow(data) == 0) return(NULL)
    return(data)
  }
  
  make_parameters <- function() {
    if (only_mutant) {
      list(
        list(column = "mt_AF", data = check_cols(ant_mt, "mt_AF", af_min), plot_title = sprintf("%s unique snps only", mt), y_title = "SNPs / Mb (×10³)", plotid = sprintf("%s_AF_unique", mt)),
        list(column = "mt_AF", data = check_cols(ant_mt_ems, "mt_AF", af_min), plot_title = sprintf("%s unique ems snps only", mt), y_title = "SNPs / Mb (×10³)", plotid = sprintf("%s_AF_unique", mt))
      )
    } else {
      list(
        list(column = "wt_AF", data = check_cols(ant_wt, "wt_AF", af_min), plot_title = sprintf("%s unique snps only", wt), y_title = "SNPs / Mb (×10³)", plotid = sprintf("%s_AF_unique", wt)),
        list(column = "mt_AF", data = check_cols(ant_mt, "mt_AF", af_min), plot_title = sprintf("%s unique snps only", mt), y_title = "SNPs / Mb (×10³)", plotid = sprintf("%s_AF_unique", mt)),
        list(column = "wt_AF", data = check_cols(ant_wt_ems, "wt_AF", af_min), plot_title = sprintf("%s unique ems snps only", wt), y_title = "SNPs / Mb (×10³)", plotid = sprintf("%s_AF_unique_ems", wt)),
        list(column = "mt_AF", data = check_cols(ant_mt_ems, "mt_AF", af_min), plot_title = sprintf("%s unique ems snps only", mt), y_title = "SNPs / Mb (×10³)", plotid = sprintf("%s_AF_unique_ems", mt))
      )
    }
  }
  
  # Generate plots
  parameters <- make_parameters()
  if (plot_data) {
    for (pmt in parameters) {
      if (!is.null(pmt$data)) {
        message(sprintf("Generating histogram: %s", pmt$plot_title))
        plot_vcfdata(
          data = pmt$data, column = pmt$column, 
          y_title = pmt$y_title, plot_title = pmt$plot_title,file_suffix = sprintf("histogram_%s", pmt$plotid),
          ylim = ylim, is_smooth = FALSE, is_rollmedian = FALSE, is_histogram = TRUE, threshold = NULL, bwidth = bwidth,
          width = width, height = height, hwidth = hwidth, hheight = hheight, dpi = dpi, device = device, 
          prefix = prefix, plots_dir = plots_dir, plot_style = plot_style, color_panel = color_panel
        )
      }
    }
  }
  message("Histogram plotting completed.")
} 