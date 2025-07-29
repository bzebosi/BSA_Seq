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