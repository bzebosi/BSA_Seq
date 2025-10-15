#' Plot Homozygosity and Allele Frequency Distributions
#' Creates histograms of allele frequency and line plots of homozygosity across chromosomes.
#' @param data A list with `results` and `out` tables.
#' @param prefix Reference genome or sample name (e.g., "B73").
#' @param bwidth Binwidth for histograms in base pairs.
#' @param mt, wt Labels for mutant and wildtype samples.
#' @param bsa_metrics Which metrics to plot: "maf", "waf", "homozygosity", or "all".
#' @param use_col Which homozygosity tracks to use: "wmd", "lft", "rmd", or "all".
#' @param plot_mode "hist", "line", or "all".
#' @param use_ems Logical, whether to use EMS-only filtered data.
#' @param af_min Allele frequency filter threshold.
#' @param plot_style Plot layout: "wrap" or "grid".
#' @param plots_dir Output directory for saved plots.
#' @return (Invisibly) a list of plot file names or plot objects.
#' @export
window_homozygosity_plot <- function(data, prefix, bwidth = 1000000, mt = "mutant", wt = "wildtype", 
                                     bsa_metrics = c("waf","maf","homozygosity", "all"),
                                     use_col = c("wmd","lft","rmd","all"), plot_mode = c("hist","line", "all"),
                                     use_ems = FALSE, af_min = 0.99, hwidth = 30, hheight = 18, 
                                     width = 45, height = 15, dpi=300, device = "png", plots_dir = "plots_dir", 
                                     facet_column = 5, line_size = 5, only_mutant = FALSE,
                                     plot_style = c("wrap","grid"), color_panel = c("blue","red")) {
  
  # If only_mutant is TRUE, drop WT label entirely
  if (only_mutant) wt <- NULL
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  # Helper mini-function : Captalize the first letter 
  capitalize_first <- function(text) {
    paste0(toupper(substr(text, 1, 1)), substr(text, 2, nchar(text)))
  }
  
  # Capitalize inbred name from prefix
  inbred <- capitalize_first(prefix)
  
  # plot_style
  plot_style <- match.arg(plot_style)
  
  # Select faceting style
  facet_layer <- if (plot_style == "wrap") {
    facet_wrap(~ CHROM, ncol = facet_column, scales = "free_x")
  } else {
    facet_grid(. ~ CHROM, scales = "free_x", space = "free_x")
  }
  
  # Histogram binwidth: convert once from bp to Mb
  make_hist_plot <- function(df_hist, column, plot_title, y_title, inbred, facet_layer, color_panel, bwidth,
                             plots_dir, device, plot_style, hwidth, hheight, width, height, dpi, file_suffix, af_min) {
    
    if (is.null(df_hist) || !all(c("CHROM","POS") %in% names(df_hist))) return(NULL)
    dfh <- df_hist %>% dplyr::mutate(PosMb = POS / 1e6)
    bwidth <- bwidth / 1e6
    
    # keep only SNPs with AF >= af_min when an AF column is provided
    if (!is.null(column) && column %in% names(dfh)) {
      dfh <- dfh %>% dplyr::filter(!is.na(.data[[column]]), .data[[column]] >= af_min)
    }
    
    # palette sized & named to chroms to avoid "Insufficient values in manual scale"
    chroms <- unique(as.character(dfh$CHROM))
    color_panel <- rep(color_panel, length.out = length(chroms))
    
    hplot <- ggplot(dfh, aes(x = PosMb, fill = CHROM)) +
      geom_histogram(binwidth = bwidth, linewidth = 0.5, aes(weight = 1 / (bwidth * 1000))) +
      facet_layer + scale_fill_manual(values = color_panel) + bsa_theme() +
      labs(title = paste0("Aligned to ", inbred, " : ", plot_title, "\n"),
           x = "\n Chromosome Position (Mbp)", y = paste0(y_title ,"\n")) 
    
    # Construct file path and save plot
    file_path <- file.path(plots_dir, paste0(inbred, "_", file_suffix, "_", plot_style, ".", device))
    plot_width <- if (plot_style == "wrap") hwidth else width
    plot_height <- if (plot_style == "wrap") hheight else height
    ggsave(filename = file_path, plot = hplot, device = device, width = plot_width, height = plot_height, dpi = dpi)
    
    message(paste0("Plot saved to: ", file_path))
    return(hplot)
  }
  
  make_line_plot <- function(df_line, column, plot_title, y_title, inbred, facet_layer, color_panel,
                             plots_dir, device, plot_style, hwidth, hheight, width, height, dpi, file_suffix) {
    
    if (is.null(df_line) || !all(c("CHROM","POS", column) %in% names(df_line))) return(NULL)
    dfl <- df_line %>% dplyr::mutate(PosMb = POS / 1e6)
    
    # palette sized & named to chroms to avoid "Insufficient values in manual scale"
    chroms <- unique(as.character(dfl$CHROM))
    color_panel <- rep(color_panel, length.out = length(chroms))
    
    lplot <- ggplot(dfl, aes(x = PosMb, y = .data[[column]], color = CHROM)) +
      geom_line(linewidth = line_size, show.legend = FALSE) + facet_layer + 
      scale_color_manual(values = color_panel) + guides(color = "none") +
      bsa_theme() +
      labs(title = paste0("Aligned to ", inbred, " : ", plot_title, "\n"),
           x = "\n Chromosome Position (Mbp)", y = paste0(y_title ,"\n"))
    
    # Construct file path and save plot
    file_path  <- file.path(plots_dir, paste0(inbred, "_", file_suffix, "_", plot_style, ".", device))
    plot_width  <- if (plot_style == "wrap") hwidth  else width
    plot_height <- if (plot_style == "wrap") hheight else height
    
    ggsave(filename = file_path, plot = lplot, device = device,
           width = plot_width, height = plot_height, dpi = dpi)
    
    message(paste0("Plot saved to: ", file_path))
    return(lplot)
  }
  
  # Which plot mode(s)?
  mode_choice <- match.arg(plot_mode, several.ok = TRUE)
  do_hist <- mode_choice %in% c("hist","all")
  do_line <- mode_choice %in% c("line","all")
  
  # Normalize requested metrics
  pm <- tolower(unique(bsa_metrics))
  if ("all" %in% pm) {
    if (only_mutant) {
      pm <- c("maf","homozygosity")
    } else {
      pm <- c("waf","maf","homozygosity")
    }
  }
  
  # data locations (per-site for hist; windowed for lines)
  res <- data$results
  out <- data$out
  
  # select EMS or ALL tables
  ant_mt <- if (use_ems) res$ant_mt_ems else res$ant_mt
  ant_wt <- if (use_ems) res$ant_wt_ems else res$ant_wt
  
  # ---- parameter builders ----
  hist_parameters <- function() {
    if (!do_hist) return(list())
    gens <- if (only_mutant) "mt" else c("mt","wt")
    outp <- list()
    for (g in gens) {
      want_metric <- if (g == "mt") "maf" else "waf"
      if (!(want_metric %in% pm)) next
      
      tbl <- if (g == "mt") ant_mt else ant_wt
      if (is.null(tbl)) next
      
      col <- if (g == "mt") "mt_AF" else "wt_AF"
      who <- if (g == "mt") mt else wt
      var <- if (use_ems) "ems" else "all"
      
      outp[[length(outp) + 1L]] <- list(
        ptype = "hist", column = col, data = tbl,
        plot_title = sprintf("%s %s snps only | homo", who, var),
        y_title = "SNPs / Mb (×10³)",
        plotid = sprintf("%s_histogram_%s", who, var)
      )
    }
    outp
  }
  
  
  # track label pretty-name
  tracks <- if (identical(use_col, "all")) c("wmd","lft","rmd") else unique(tolower(use_col))
  tracks <- intersect(tracks, c("wmd","lft","rmd"))
  if (!length(tracks)) tracks <- "wmd"
  track_name <- function(tr) switch(tr, wmd = "winmedian", lft = "locfit", rmd = "rollmedian", tr)
  
  # after (safe & identical behavior for non-mutant-only)
  hom_mt <- {
    x <- if (use_ems) out$homoz_mt_ems else out$homoz_mt
    if (!is.null(x)) x$windows else NULL
  } 
  
  hom_wt <- if (only_mutant) {
    NULL
  } else {
    x <- if (use_ems) out$homoz_wt_ems else out$homoz_wt
    if (!is.null(x)) x$windows else NULL
  }
  
  homo_parameters <- function() {
    if (!do_line || !("homozygosity" %in% pm)) return(list())
    
    gens <- if (only_mutant) "mt" else c("mt","wt")
    outp <- list()
    
    for (g in gens) {
      tbl <- if (g == "mt") hom_mt else hom_wt
      if (is.null(tbl)) next
      
      who <- if (g == "mt") mt else wt
      var <- if (use_ems) "ems" else "all"
      
      for (tr in tracks) {
        # accept either Homozygosity_rmd or Homozygosity_ma
        col <- switch(tr, 
                      wmd = "Homozygosity_wmd", lft = "Homozygosity_lft",
                      rmd = "Homozygosity_rmd")
        
        if (!col %in% names(tbl)) next
        
        outp[[length(outp) + 1L]] <- list(
          ptype = "line", column = col, data = tbl,
          plot_title = sprintf("%s | Homozygosity | %s", who, track_name(tr)),
          y_title = "Homozygosity", plotid   = sprintf("%s_homo_%s_%s", who, tr, var)
        )
      }
    }
    outp
  }
  
  params <- c(hist_parameters(), homo_parameters())
  if (!length(params)) stop("No matching data found for requested metrics/mode.")
  
  # ---- loop: single dispatch by ptype (no AF logic) ----
  out_files <- list()
  for (p in params) {
    if (is.null(p$data)) next
    
    if (identical(p$ptype, "hist")) {
      g <- make_hist_plot(df_hist = p$data, column = p$column, plot_title=p$plot_title, 
                          y_title=p$y_title, inbred=inbred,facet_layer = facet_layer, 
                          color_panel = color_panel,bwidth = bwidth, plots_dir = plots_dir,
                          device = device, plot_style = plot_style, hwidth = hwidth, hheight = hheight,
                          width = width, height = height, dpi = dpi, file_suffix = p$plotid, af_min = af_min
      )
    } else {
      g <- make_line_plot(df_line = p$data, column = p$column, plot_title=p$plot_title, 
                          y_title=p$y_title, inbred=inbred,facet_layer = facet_layer, 
                          color_panel = color_panel,plots_dir = plots_dir, device = device, 
                          plot_style = plot_style, hwidth = hwidth, hheight = hheight,width = width, 
                          height  = height, dpi = dpi, file_suffix = p$plotid
      )
    }
    if (!is.null(g)) out_files[[p$plotid]] <- TRUE
  }
  invisible(out_files)
}

