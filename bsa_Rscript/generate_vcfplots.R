#' Generate VCF Data Plots
#' @description Creates SNP-based plots (scatter, smooth, histogram)
#' @param data Data frame containing VCF SNP data.
#' @param prefix Prefix for file naming.
#' @param column Name of column to plot on the y-axis.
#' @param y_title Y-axis label.
#' @param plot_title Title of the plot.
#' @param file_suffix Suffix for saved file.
#' @param ylim Y-axis limits.
#' @param is_smooth Apply smoothing (Loess).
#' @param is_rollmedian Apply rolling median.
#' @param rollmed Window size for rolling median.
#' @param is_histogram If TRUE, create a histogram.
#' @param threshold Optional threshold line.
#' @param hwidth Histogram width.
#' @param hheight Histogram height.
#' @param width Plot width.
#' @param height Plot height.
#' @param dpi Resolution.
#' @param device Output format (e.g., "tiff", "png").
#' @param plots_dir Directory to save plots.
#' @param nn_prop Smoothing parameter.
#' @export
generate_vcfplots <- function(data, prefix, column, y_title, plot_title, file_suffix, ylim = NULL, is_smooth = FALSE, is_rollmedian = TRUE, bwidth = 1000000,
                              is_histogram=FALSE, threshold = NULL, rollmedian = 501, hwidth, hheight, width, height, dpi, device, plots_dir, nn_prop){
  # Ensure plots directory exists
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  # Helper mini-function : Captalize the first letter 
  capitalize_first <- function(text) {
    substr(text, 1, 1) <- toupper(substr(text, 1, 1)) 
    return(text)
  }
  
  # Capitalize inbred name from prefix
  inbred <- capitalize_first(prefix)
  
  # Define color panel for chromosomes
  num_chromosomes <- length(unique(data$CHROM))
  color_panel <- rep(c("blue", "red"), length.out = num_chromosomes)
  
  # Generate scatter or smooth plot or histogram
  if(!is_histogram){
    plot <- ggplot(data, aes(x = POS, y = !!sym(column), color = CHROM)) +
      facet_grid(. ~ CHROM, scales = "free_x", space = "free_y") +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
      scale_x_continuous(breaks = seq(min(data$POS), max(data$POS), by = 20000000), labels = scales::comma) +
      scale_color_manual(values = color_panel) + guides(color = FALSE) +
      labs(title = paste0("Aligned to ", inbred, " : ", " ", plot_title), x = "\n Chromosome Position (bp) \n", y = paste0("\n", y_title, "\n")) + bsa_theme()
    
    # Add geom_point only if is_smooth is FALSE
    if (is_smooth) {
      plot <- plot + stat_smooth(method = "locfit", formula = y ~ lp(x, nn = nn_prop))
    } else if (is_rollmedian) {
      plot <- plot + geom_point(size = 0.5) + 
        geom_line(aes(y = rollmedian(!!sym(column), rollmedian, na.pad = TRUE)), color = "black", linewidth = 4)
    } else {
      plot <- plot + geom_point(size = 1)
    }
    
    # Add horizontal threshold line if applicable
    if (!is.null(threshold)) {plot <- plot + geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1.2)}
    
  } else {
    plot <- ggplot(data, aes(x=POS, fill=CHROM)) + geom_histogram(binwidth=bwidth,  alpha = 0.7) +
      facet_wrap(~ CHROM, ncol = 5, scales = "free_x") + scale_fill_manual(values = color_panel) +
      theme(axis.ticks.x  = element_line(color = "black", linewidth = 1),
            axis.minor.ticks.length = unit(0.5, "cm"), axis.text.x = element_blank()) + guides(fill = FALSE) +
      labs(title = paste0("Aligned to ", inbred, " : ", " ", plot_title), x = "\n Chromosome Position (bp) \n", y = paste0("\n", y_title, "\n")) + bsa_theme()
  }
  
  # Apply fixed y-axis limits only if ylim is not NULL
  if (!is.null(ylim)) {plot <- plot + coord_cartesian(ylim = ylim)}
  
  # Construct file path
  file_path <- file.path(plots_dir, paste0(inbred, "_", file_suffix, ".", device))
  
  # Save plot
  ggsave(filename = file_path, plot = plot, device = device, 
         width = ifelse(is_histogram, hwidth, width), 
         height = ifelse(is_histogram, hheight, height), dpi = dpi)
  
  message(paste0("Plot saved to: ", file_path))
}
