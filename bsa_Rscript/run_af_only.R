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



