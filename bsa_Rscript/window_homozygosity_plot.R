window_homozygosity_plot <- function(data, prefix, ylim = NULL, bwidth = 1000000, mt = "mutant", wt = "wildtype", plot_mutant = FALSE,
                         plot_metrics = c("hist","homo","median","locfit","all"), af_min = 0.99, hwidth = 30, hheight = 18, 
                         width = 45, height =15, dpi, device, plots_dir, facet_column = 5, line_size = 5,
                         plot_style = c("wrap","grid"), color_panel = c("blue","red")) {
  
  # Ensure plots directory exists
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  # Helper mini-function : Captalize the first letter 
  capitalize_first <- function(text) {
    substr(text, 1, 1) <- toupper(substr(text, 1, 1)) 
    return(text)
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
  
  homoz <- list(
    mt_all = data$out$homoz_mt$windows,
    mt_ems = data$out$homoz_mt_ems$windows,
    wt_all = if (!plot_mutant) data$out$homoz_wt$windows else NULL,
    wt_ems = if (!plot_mutant) data$out$homoz_wt_ems$windows else NULL
  )
  
  ant_mt_all = data$results$ant_mt
  ant_mt_ems = data$results$ant_mt_ems
  ant_wt_all = if (!plot_mutant) data$results$ant_wt else NULL
  ant_wt_ems = if (!plot_mutant) data$results$ant_wt_ems else NULL
  
  hist_parameters <- function() {
    if (plot_mutant) {
      list(
        list(ptype = "hist", column = "mt_AF", data = ant_mt_all, plot_title = sprintf("%s all snps only | homo", mt), y_title = "SNPs / Mb (×10³)", plotid = sprintf("%s_histogram_all", mt)),
        list(ptype = "hist", column = "mt_AF", data = ant_mt_ems, plot_title = sprintf("%s ems snps only | homo", mt), y_title = "SNPs / Mb (×10³)", plotid = sprintf("%s_histogram_ems", mt))
      )
    } else {
      list(
        list(ptype = "hist", column = "mt_AF", data = ant_mt_all, plot_title = sprintf("%s all snps only | homo", mt), y_title = "SNPs / Mb (×10³)", plotid = sprintf("%s_histogram_all", mt)),
        list(ptype = "hist", column = "mt_AF", data = ant_mt_ems, plot_title = sprintf("%s ems snps only | homo", mt), y_title = "SNPs / Mb (×10³)", plotid = sprintf("%s_histogram_ems", mt)),
        list(ptype = "hist", column = "wt_AF", data = ant_wt_all, plot_title = sprintf("%s all snps only | homo", wt), y_title = "SNPs / Mb (×10³)", plotid = sprintf("%s_histogram_all", wt)),
        list(ptype = "hist", column = "wt_AF", data = ant_wt_ems, plot_title = sprintf("%s ems snps only | homo", wt), y_title = "SNPs / Mb (×10³)", plotid = sprintf("%s_histogram_ems", wt))
      )
    }
  }
  
  hom_parameters <- function() {
    if (plot_mutant) {
      list(
        list(ptype = "line", column = "Homozygosity", data = homoz$mt_all, plot_title = sprintf("%s all snps only", mt), y_title = "Homozygosity", plotid = sprintf("%s_homozygosity_all", mt)),
        list(ptype = "line", column = "Homozygosity", data = homoz$mt_ems, plot_title = sprintf("%s ems snps only", mt), y_title = "Homozygosity", plotid = sprintf("%s_homozygosity_ems", mt))
      )
    } else {
      list(
        list(ptype = "line", column = "Homozygosity", data = homoz$mt_all, plot_title = sprintf("%s all snps only", mt), y_title = "Homozygosity", plotid = sprintf("%s_homozygosity_all", mt)),
        list(ptype = "line", column = "Homozygosity", data = homoz$mt_ems, plot_title = sprintf("%s ems snps only", mt), y_title = "Homozygosity", plotid = sprintf("%s_homozygosity_ems", mt)),
        list(ptype = "line", column = "Homozygosity", data = homoz$wt_all, plot_title = sprintf("%s all snps only", wt), y_title = "Homozygosity", plotid = sprintf("%s_homozygosity_all", wt)),
        list(ptype = "line", column = "Homozygosity", data = homoz$wt_ems, plot_title = sprintf("%s ems snps only", wt), y_title = "Homozygosity", plotid = sprintf("%s_homozygosity_ems", wt))
      )
    }
  }
  
  locfit_parameters <- function() {
    if (plot_mutant) {
      list(
        list(ptype = "line", column = "Homozygosity_lft", data = homoz$mt_all, plot_title = sprintf("%s all snps only | locfit", mt), y_title = "Homozygosity", plotid = sprintf("%s_homozygosity_fit_all", mt)),
        list(ptype = "line", column = "Homozygosity_lft", data = homoz$mt_ems, plot_title = sprintf("%s ems snps only | locfit", mt), y_title = "Homozygosity", plotid = sprintf("%s_homozygosity_fit_ems", mt))
      )
    } else {
      list(
        list(ptype = "line", column = "Homozygosity_lft", data = homoz$mt_all, plot_title = sprintf("%s all snps only | locfit" , mt), y_title = "Homozygosity", plotid = sprintf("%s_homozygosity_fit_all", mt)),
        list(ptype = "line", column = "Homozygosity_lft", data = homoz$mt_ems, plot_title = sprintf("%s ems snps only | locfit", mt), y_title = "Homozygosity", plotid = sprintf("%s_homozygosity_fit_ems", mt)),
        list(ptype = "line", column = "Homozygosity_lft", data = homoz$wt_all, plot_title = sprintf("%s all snps only | locfit", wt), y_title = "Homozygosity", plotid = sprintf("%s_homozygosity_fit_all", wt)),
        list(ptype = "line", column = "Homozygosity_lft", data = homoz$wt_ems, plot_title = sprintf("%s ems snps only | locfit", wt), y_title = "Homozygosity", plotid = sprintf("%s_homozygosity_fit_ems", wt))
      )
    }
  }
  
  rollmedian_parameters <- function() {
    if (plot_mutant) {
      list(
        list(ptype = "line", column = "Homozygosity_ma", data = homoz$mt_all, plot_title = sprintf("%s all snps only | rollmedian", mt), y_title = "Homozygosity", plotid = sprintf("%s_homozygosity_rmedian_all", mt)),
        list(ptype = "line", column = "Homozygosity_ma", data = homoz$mt_ems, plot_title = sprintf("%s ems snps only | rollmedian", mt), y_title = "Homozygosity", plotid = sprintf("%s_homozygosity_rmedian_ems", mt))
      )
    } else {
      list(
        list(ptype = "line", column = "Homozygosity_ma", data = homoz$mt_all, plot_title = sprintf("%s all snps only | rollmedian", mt), y_title = "Homozygosity", plotid = sprintf("%s_homozygosity_rmedian_all", mt)),
        list(ptype = "line", column = "Homozygosity_ma", data = homoz$mt_ems, plot_title = sprintf("%s ems snps only | rollmedian" , mt), y_title = "Homozygosity", plotid = sprintf("%s_homozygosity_rmedian_ems", mt)),
        list(ptype = "line", column = "Homozygosity_ma", data = homoz$wt_all, plot_title = sprintf("%s all snps only | rollmedian", wt), y_title = "Homozygosity", plotid = sprintf("%s_homozygosity_rmedian_all", wt)),
        list(ptype = "line", column = "Homozygosity_ma", data = homoz$wt_ems, plot_title = sprintf("%s ems snps only | rollmedian", wt), y_title = "Homozygosity", plotid = sprintf("%s_homozygosity_rmedian_ems", wt))
      )
    }
  }
  
  # choose what to plot based on plot_metrics
  plot_metrics <- match.arg(plot_metrics,
    choices = c("hist","homo","median","locfit","all"),
    several.ok = TRUE
  )
  params <- list()
  
  if ("all" %in% plot_metrics) {
    params <- c(hist_parameters(), hom_parameters(), locfit_parameters(), rollmedian_parameters())
  } else {
    params <- list()
    if ("homo"   %in% plot_metrics) params <- c(params, hom_parameters())
    if ("hist"   %in% plot_metrics) params <- c(params, hist_parameters())
    if ("median" %in% plot_metrics) params <- c(params, rollmedian_parameters())
    if ("locfit" %in% plot_metrics) params <- c(params, locfit_parameters())
  }
  
  # ---- loop: single dispatch by ptype (no AF logic) ----
  out_files <- list()
  for (p in params) {
    if (is.null(p$data)) next
    
    if (identical(p$ptype, "hist")) {
      g <- make_hist_plot(
        df_hist = p$data, column = p$column, plot_title=p$plot_title, 
        y_title=p$y_title, inbred=inbred,
        facet_layer = facet_layer, color_panel = color_panel,
        bwidth = bwidth, plots_dir = plots_dir,
        device = device, plot_style = plot_style, hwidth = hwidth, hheight = hheight,
        width = width, height  = height, dpi = dpi, file_suffix = p$plotid,af_min = af_min
      )
    } else {
      g <- make_line_plot(
        df_line = p$data, column = p$column, plot_title=p$plot_title, 
        y_title=p$y_title, inbred=inbred,
        facet_layer = facet_layer, color_panel = color_panel,plots_dir = plots_dir,
        device = device, plot_style = plot_style, hwidth = hwidth, hheight = hheight,
        width = width, height  = height, dpi = dpi, file_suffix = p$plotid
      )
    }
    
    if (!is.null(g)) out_files[[p$plotid]] <- TRUE
  }
  invisible(out_files)
}



# --- test run ---
#plot_vcfdata(
#  data = res_both,
#  mt = "S7C6508K",
#  wt = "S7A6508K",
#  prefix = "B73",
#  plots_dir = "/Users/zebosi/Downloads/plots_outiio",
#  plot_metrics = c("hist","median","locfit"),
#  device = "png",
#  plot_style = "grid",
#  plot_mutant = TRUE,
#  width = 45, height = 13, dpi = 300
#)

