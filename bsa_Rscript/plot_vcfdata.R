#' Plot BSA-Seq Data by Chromosome
#'
#' Generates SNP plots or histograms for allele frequency metrics (AF, AFD, ED, G-statistics, etc.)
#' across chromosomes from BSA-Seq data. Supports smoothing via locfit or rolling median. Facets
#' by chromosome using either facet_wrap or facet_grid`.
#' @param data A data.frame or data.table containing SNP-level information including CHROM, POS, and the metric column.
#' @param prefix A prefix (usually sample name or inbred line) used for plot titles and filenames.
#' @param column The name of the numeric column in data to plot on the y-axis (e.g., "AFD", "G", "ED").
#' @param y_title Y-axis label.
#' @param plot_title Plot title.
#' @param file_suffix Suffix to use when naming the output file (before file extension).
#' @param ylim Optional y-axis limits as a numeric vector of length 2.
#' @param is_smooth Logical. If \code{TRUE}, applies locfit smoothing.
#' @param is_rollmedian Logical. If \code{TRUE}, applies rolling median smoothing.
#' @param bwidth Bin width (in bp) for histogram mode.
#' @param is_histogram Logical. If \code{TRUE}, plots histogram instead of scatterplot.
#' @param threshold Optional horizontal line (e.g., significance threshold).
#' @param rollmedian Window size for rolling median smoothing.
#' @param hwidth Width of the plot for histogram or wrap mode.
#' @param hheight Height of the plot for histogram or wrap mode.
#' @param width Width of the plot for grid mode.
#' @param height Height of the plot for grid mode.
#' @param dpi Resolution in dots per inch.
#' @param device Output format for saved plot. One of \code{"png"}, \code{"pdf"}, etc.
#' @param plots_dir Directory where plots should be saved.
#' @param nn_prop Proportion of nearest neighbors for locfit smoothing (between 0 and 1).
#' @param point_size Size of points in scatter plots.
#' @param line_size Thickness of smoothed line or histogram outline.
#' @param alpha_size Alpha transparency for histogram bars.
#' @param facet_column Number of columns for \code{facet_wrap} layout (used if \code{plot_style = "wrap"}).
#' @param plot_style Faceting style. Either \code{"wrap"} (default) or \code{"grid"}.
#' @param remove_x_text Logical. If \code{TRUE}, removes x-axis text and ticks.
#' @param color_panel Vector of colors to use for different chromosomes. Will be repeated to match number of chromosomes.
#' @return (Invisibly) returns the ggplot object. Also saves the plot to disk.
#' @examples
#' \dontrun{
#' plot_vcfdata(data = my_data, prefix = "Ts5", column = "AFD",
#'   y_title = "Δ SNP-index", plot_title = "AFD Plot",
#'   file_suffix = "afd_example", device = "png", plots_dir = "plots/",
#'   is_smooth = TRUE, plot_style = "wrap", nn_prop = 0.1)
#' }
#' @export
plot_vcfdata <- function(data, prefix, column, y_title, plot_title, file_suffix, 
                               ylim = NULL, is_smooth = FALSE, is_rollmedian = TRUE, bwidth = 1000000,
                               is_histogram=FALSE, threshold = NULL, rollmedian = 501, 
                               hwidth, hheight, width, height, dpi, device, plots_dir, 
                               nn_prop, point_size = 4, line_size = 4, alpha_size = 1, facet_column = 5, 
                               plot_style = c("wrap", "grid"), remove_x_text = TRUE, color_panel = c("blue", "red")){
  # Ensure plots directory exists
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  # Helper mini-function : Captalize the first letter 
  capitalize_first <- function(text) {
    substr(text, 1, 1) <- toupper(substr(text, 1, 1)) 
    return(text)
  }
  
  # Capitalize inbred name from prefix
  inbred <- capitalize_first(prefix)
  
  # Add PositionMb column for better x-axis readability  # Add Position in Mb
  data <- data %>% mutate(PositionMb = POS / 1e6)
  bwidth <- bwidth / 1e6 
  
  # Define color panel for chromosomes
  num_chrom <- length(unique(data$CHROM))
  color_panel <- rep(color_panel, length.out = num_chrom)
  
  # plot_style
  plot_style <- match.arg(plot_style)
  
  # Select faceting style
  facet_layer <- if (plot_style == "wrap") {
    facet_wrap(~ CHROM, ncol = facet_column, scales = "free_x")
  } else {
    facet_grid(. ~ CHROM, scales = "free_x", space = "free_x")
  }
  
  # Generate scatter or smooth plot or histogram
  if(!is_histogram){
    plot <- ggplot(data, aes(x = PositionMb, y = !!sym(column), color = CHROM)) +
      facet_layer +
      scale_color_manual(values = color_panel) + guides(color = "none") +
      labs(title = paste0("Aligned to ", inbred, " : ", " ", plot_title), x = "Chromosome Position (Mb)", y = paste0("\n", y_title, "\n")) + 
      bsa_theme() +   theme(
        axis.text.x = if (remove_x_text) element_blank() else element_text(angle = 45, hjust = 1, vjust = 1, margin = margin(t = -1)),
        axis.ticks.x  = if (remove_x_text) element_blank() else element_line(),
        axis.title.x = element_text(color = "black")
      )
    
    # Add geom_point only if is_smooth is FALSE
    if (is_smooth) {
      plot <- plot + stat_smooth(method = "locfit", formula = y ~ lp(x, nn = nn_prop), linewidth = line_size)
    } else if (is_rollmedian) {
      plot <- plot + geom_point(size = point_size) + 
        geom_line(aes(y = rollmedian(!!sym(column), rollmedian, na.pad = TRUE)), color = "black", linewidth = line_size)
    } else {
      plot <- plot + geom_point(size = point_size)
    }
    
    # Add horizontal threshold line if applicable
    if (!is.null(threshold)) {plot <- plot + geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = line_size)}
    
  } else {
    plot <- ggplot(data, aes(x=PositionMb, fill=CHROM)) + geom_histogram(binwidth=bwidth,  alpha = alpha_size,  linewidth = 6) +
      facet_layer + scale_fill_manual(values = color_panel) +
      labs(title = paste0("Aligned to ", inbred, " : ", " ", plot_title), x = "\n Chromosome Position (bp) \n", y = paste0("\n", y_title, "\n")) + 
      bsa_theme() +   theme(
        axis.text.x = if (remove_x_text) element_blank() else element_text(angle = 45, hjust = 1, vjust = 1, margin = margin(t = -1)),
        axis.ticks.x  = if (remove_x_text) element_blank() else element_line(),
        axis.title.x = element_text(color = "black")) + 
      scale_y_continuous(
        labels = scales::label_number(scale = 1e-3, accuracy = 0.01, trim = TRUE)
      )
  }
  
  # Apply fixed y-axis limits only if ylim is not NULL
  if (!is.null(ylim)) {plot <- plot + coord_cartesian(ylim = ylim)}
  
  # Construct file path and save plot
  file_path <- file.path(plots_dir, paste0(inbred, "_", file_suffix, "_", plot_style, ".", device))
  plot_width <- if (plot_style == "wrap") hwidth else width
  plot_height <- if (plot_style == "wrap") hheight else height
  ggsave(filename = file_path, plot = plot, device = device, width = plot_width, height = plot_height, dpi = dpi)
  
  message(paste0("Plot saved to: ", file_path))
  return(invisible(plot))
}




