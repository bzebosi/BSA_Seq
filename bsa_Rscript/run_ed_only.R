run_ed_only <- function(
    wt_mt = NULL, wt = "wildtype", mt = "mutant", prefix = "sample", plots_dir = "plots", 
    rollmedian = 501, ylim = NULL, only_mutant = FALSE, device = "png", plot_data = TRUE,
    width = 45, height = 12, hwidth = 30, hheight = 18, dpi = 300, nn_prop = 0.1, ed_min = 1e-5,
    plot_mode = "both", plot_style = "wrap", color_panel = c("blue", "red")){
  
  
  if (only_mutant) {
    message("Skipping ed plots in mutant-only mode.")
    return(NULL)
  }
  
  # fix: include ed_min inside check_cols
  check_cols <- function(data, colname, ed_min) {
    if (is.null(data)) return(NULL)
    if (!(colname %in% colnames(data))) return(NULL)
    data <- data[data[[colname]] >= ed_min, ]
    if (nrow(data) == 0) return(NULL)
    return(data)
  }
  
  make_parameters <- function() {
      list(
        list(column = "ED", data = check_cols(wt_mt, "ED" , ed_min), plot_title = sprintf("ED %s & %s", wt, mt), y_title = "ED", plotid = sprintf("ED_%s_%s", wt, mt)), 
        list(column = "ED4", data = check_cols(wt_mt, "ED4", ed_min), plot_title = sprintf("ED4 %s & %s", wt, mt), y_title = "ED4", plotid = sprintf("ED4_%s_%s", wt, mt))
      )
  }
  
  plot_with_mode <- function(mode, data, column, plot_title, y_title, plotid) {
    if (mode == "locfit" || mode == "both") {
      plot_vcfdata(
        data = data, column = column, 
        y_title = y_title, plot_title = plot_title, file_suffix = sprintf("ed_locfit_%s", plotid),
        ylim = ylim, is_smooth = TRUE, is_rollmedian = FALSE, rollmedian = 0, is_histogram = FALSE, 
        nn_prop = nn_prop, width = width, height = height, hwidth = hwidth, hheight = hheight, dpi = dpi, device = device, 
        prefix = prefix, plots_dir = plots_dir, plot_style = plot_style, color_panel = color_panel
      )
    }
    
    if ((mode == "rollmedian" || mode == "both") && rollmedian > 0) {
      plot_vcfdata(
        data = data, column = column, 
        y_title = y_title, plot_title = plot_title, file_suffix = sprintf("ed_rollmedian_%s", plotid),
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
  for (pmt in parameters) {
    if (!is.null(pmt$data)) {
      message(sprintf("Generating ED plot: %s", pmt$plot_title))
      plot_with_mode(plot_mode, pmt$data, pmt$column, pmt$plot_title, pmt$y_title, pmt$plotid)
    }
  }
  
  message("ED-only plots completed.")
}


