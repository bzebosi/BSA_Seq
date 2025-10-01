window_bsa_plot <- function(data, prefix, mt = "mutant", wt = "wildtype", only_mutant = FALSE, use_ems = FALSE,
                           bsa_metrics = c("waf","maf","ed","ed4","g","afd", "all"), use_col = c("wmd","lft","rmd","all"), hwidth = 30, hheight = 18, 
                           width = 45, height =15, dpi=300, device = "png", plots_dir = "plots_dir", facet_column = 5, line_size = 5,
                           plot_style = c("wrap","grid"), color_panel = c("blue","red")) {
  # Ensure plots directory exists
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  # Helper: capitalize first letter
  capitalize_first <- function(text) {
    substr(text, 1, 1) <- toupper(substr(text, 1, 1))
    text
  }
  inbred <- capitalize_first(prefix)
  
  # plot style & facet layer
  plot_style <- match.arg(plot_style)
  facet_layer <- if (plot_style == "wrap") {
    facet_wrap(~ CHROM, ncol = facet_column, scales = "free_x")
  } else {
    facet_grid(. ~ CHROM, scales = "free_x", space = "free_x")
  }
  
  # Simple ggplot line helper
  make_line_plot <- function(df_line, column, plot_title, y_title, inbred, facet_layer, color_panel,
                             plots_dir, device, plot_style, hwidth, hheight, width, height, dpi, file_suffix) {
    if (is.null(df_line) || !all(c("CHROM","POS", column) %in% names(df_line))) return(NULL)
    dfl <- df_line %>% dplyr::mutate(PosMb = POS / 1e6)
    dfl <- dfl %>% dplyr::filter(is.finite(.data[[column]]))
    
    chroms <- unique(as.character(dfl$CHROM))
    color_panel <- rep(color_panel, length.out = length(chroms))
    
    lplot <- ggplot(dfl, aes(x = PosMb, y = .data[[column]], color = CHROM)) +
      geom_line(linewidth = line_size, show.legend = FALSE) +
      facet_layer + scale_color_manual(values = color_panel) + guides(color = "none") +
      bsa_theme() +
      labs(title = paste0("Aligned to ", inbred, " : ", plot_title, "\n"),
           x = "\n Chromosome Position (Mbp)", y = paste0(y_title ,"\n"))
    
    file_path  <- file.path(plots_dir, paste0(inbred, "_", file_suffix, "_", plot_style, ".", device))
    plot_width  <- if (plot_style == "wrap") hwidth  else width
    plot_height <- if (plot_style == "wrap") hheight else height
    ggsave(filename = file_path, plot = lplot, device = device, width = plot_width, height = plot_height, dpi = dpi)
    message(paste0("Plot saved to: ", file_path))
    lplot
  }
  
  # grab windows list
  win <- data$windows
  if (is.null(win) || !length(win)) stop("No window tables found in data$windows.")
  
  # Pick keys (EMS variant if requested and available)
  pick_keys <- function(win, use_ems = FALSE) {
    base <- c("mt_AF","wt_AF")
    ems  <- c("mt_AF_ems","wt_AF_ems")
    if (use_ems && any(ems  %in% names(win))) return(ems [ems  %in% names(win)])
    if (!use_ems && any(base %in% names(win))) return(base[base %in% names(win)])
    # fallback: use whatever exists
    keys <- c(base, ems)
    keys[keys %in% names(win)]
  }
  
  # normalize which metrics to show
  pm <- tolower(unique(bsa_metrics))
  af_keys <- pick_keys(win, use_ems = use_ems)
  if (isTRUE(only_mutant)) af_keys <- af_keys[grepl("^mt_", af_keys)]
  if (!("all" %in% pm)) {
    if (!("maf" %in% pm)) af_keys <- af_keys[!grepl("^mt_", af_keys)]
    if (!("waf" %in% pm)) af_keys <- af_keys[!grepl("^wt_", af_keys)]
  }
  
  # which track suffixes to use?
  tracks <- if (identical(use_col, "all")) c("wmd","lft","rmd") else unique(tolower(use_col))
  
  # ---- parameter builder ----
  af_parameters <- function() {
    out <- list()
    
    for (key in af_keys) {
      df <- win[[key]]
      if (is.null(df)) next
      
      # columns remain mt_AF_* / wt_AF_* even when key ends with _ems
      base <- sub("_ems$","", key)                  # "mt_AF" or "wt_AF"
      who  <- if (grepl("^mt_", base)) mt else wt   # label (uses your mt/wt args)
      ylab <- "AF" # dynamic y-axis label
      
      for (tr in tracks) {
        col <- paste0(base, "_", tr)                # mt_AF_wmd / wt_AF_lft / mt_AF_rmd
        if (!col %in% names(df)) next
        
        out[[length(out)+1L]] <- list(column = col, data = df, y_title = ylab,
          plot_title = sprintf("%s | AF | %s", who, 
                               switch(tr, wmd = " windowmedian", lft = "locfit",
                                      rmd = "rollingmedian", tr)),
          plotid = sprintf("%s_af_%s", who,
                               switch(tr, wmd = "wmd", lft = "locfit",
                                      rmd = "rollmedian", tr))
        )
      }
    }
    out
  }
  
  # joint metric keys requested (from wt_mt-derived windows)
  all_joint_names <- c("ed","ed4","afd","g")
  joint_wanted    <- if ("all" %in% pm) all_joint_names else pm[pm %in% all_joint_names]
  joint_keys      <- intersect(joint_wanted, names(win))
  
  # Joint metrics (ED, ED4, AFD, G)
  joint_parameters <- function() {
    out <- list()
    if (!length(joint_keys)) return(out)
    
    for (m in joint_keys) {
      df <- win[[m]]
      if (is.null(df)) next
      base <- toupper(m)
      ylab <- base
      
      for (tr in tracks) {
        col <- paste0(base, "_", tr)  # e.g. ED_wmd, G_lft
        if (!col %in% names(df)) next
        
        out[[length(out) + 1L]] <- list(column = col, data = df, y_title = ylab,
          plot_title = sprintf("%s | %s (%s vs %s)", base,
                               switch(tr, wmd = "winmedian", lft = "locfit", rmd = "rollmedian", tr), wt, mt),
          plotid = sprintf("%s_%s_%s_vs_%s",
                           m, switch(tr, wmd = "wmd", lft = "lft", rmd = "rmd", tr), wt, mt)
        )
      }
    }
    out
  }
  

  # Build params and render
  params <- c(af_parameters(), joint_parameters())
  if (!length(params)) stop("No matching AF/joint columns found for the requested metrics/tracks.")
  
  out_files <- list()
  for (p in params) {
    g <- make_line_plot(
      df_line = p$data, column = p$column,
      plot_title = p$plot_title, y_title = p$y_title, inbred = inbred,
      facet_layer = facet_layer, color_panel = color_panel,
      plots_dir = plots_dir, device = device, plot_style = plot_style,
      hwidth = hwidth, hheight = hheight, width = width, height = height,
      dpi = dpi, file_suffix = p$plotid
    )
    if (!is.null(g)) out_files[[p$plotid]] <- TRUE
  }
  
  invisible(out_files)
}


#' vcf_dir = "/Users/zebosi/Documents/osu_postdoc/BSA/S7_6508K/data/snps"
#' output_dir = "/Users/zebosi/Documents/osu_postdoc/BSA/S7_6508K/interval_analysis"
#' plots_dir = "/Users/zebosi/Documents/osu_postdoc/BSA/S7_6508K/ppl_analysis"
#' wt <- c("S7A6508K")
#' mt <- c("S7B6508K")
#' pattern = "snps\\.tsv$"
#' min_DP=10
#' min_QUAL=10
#' prefix = c("b73")

#' a <- import_vcfdata(vcf_dir, prefix, pattern, Genotypes = list(wt = wt, mt = mt),
#'                    min_DP = 5, min_QUAL = 5, only_mutant = FALSE)

#' d <- analyze_vcfdata(a, prefix, save_results = FALSE, bsa_metrics = "all", output_dir = "post_analysis", only_mutant = FALSE)
#'


#' win1 <- window_bsa_compute_all(d, metrics = c("all"),
#'                                af_col = "both", use_ems = TRUE, only_mutant = FALSE,
#'                                window_size = 2e6, step_size = 1e5,rollmedian = 100L, nn_prop = 0.1,
#'                                find_intervals = TRUE, offhold = 0.90, min_vsize = 1e6)

ak<- window_bsa_plot(data=win1, prefix, mt = mt, wt = wt, only_mutant = FALSE, use_ems = FALSE,
                     bsa_metrics = c("all"),
                     use_col = c("lft","rmd"), hwidth = 30, hheight = 18, 
                     width = 45, height =15, plots_dir = plots_dir, facet_column = 5, line_size = 5,
                     plot_style = "grid", color_panel = c("blue","red")) 

