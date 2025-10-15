#' @param vcf_dir Directory containing the VCF/TSV files.
#' @param prefix Prefix used in file names.
#' @param pattern Regex pattern to match input files.
#' @param Genotypes Named list with wt and mt identifiers.
#' @param output_dir Directory to save outputs.
#' @param min_DP Minimum read depth filter.
#' @param min_QUAL Minimum quality filter.
#' @param only_mutant Logical; run mutant-only mode.
#' @param use_ems Logical; include EMS-type SNPs.
#' @param save_results Logical; save analysis results.
#' @param save_interval Logical; save interval results to Excel.
#' @param rollmedian Window size for rolling median smoothing.
#' @param af_col Allele frequency column to use.
#' @param bsa_metrics Metrics to calculate (afd, ed, g, all).
#' @param window_size Sliding window size (bp).
#' @param step_size Step size for sliding windows (bp).
#' @param nn_prop Proportion of neighbors for locfit smoothing.
#' @param af_min Allele frequency cutoff for homozygosity.
#' @param find_intervals Logical; find genomic peak intervals.
#' @param offhold Proportion threshold for peak cutoff.
#' @param min_vsize Minimum interval size (bp).
#' @param use_col Columns to use for interval detection.
#' @param do_plot Logical; generate plots.
#' @param plots_dir Directory for plot outputs.
#' @param device Plot device (e.g., "png", "pdf").
#' @param dpi Resolution for plots.
#' @param hwidth,hheight Dimensions for wrapped plots.
#' @param width,height Dimensions for grid plots.
#' @param bwidth Histogram bin width (bp).
#' @param facet_column Number of columns for faceting.
#' @param line_size Line size for plots.
#' @param plot_mode Which plot types ("hist","line","all").
#' @param plot_style Faceting style ("wrap" or "grid").
#' @param color_panel Vector of colors for chromosomes.
#' @return A list with imported data, analysis results, and homozygosity outputs.
#' @export
window_homozygosity_single <- function(
    vcf_dir, prefix, pattern, Genotypes = list(wt = "wildtype", mt = "mutant"),
    output_dir = "post_analysis", min_DP = 5, min_QUAL = 5, only_mutant = FALSE, 
    use_ems = TRUE, save_results = FALSE, save_interval = TRUE, find_intervals = TRUE,
    rollmedian = 501L, af_col = NULL, window_size = 5e5, step_size = 1e5, nn_prop = 0.1, 
    af_min = 0.99, offhold = 0.80, min_vsize = 0L,
    bsa_metrics = c("waf","maf","homozygosity", "all"), plot_mode = c("hist","line","all"), 
    use_col = c("wmd","lft","rmd","all"), do_plot = TRUE, plots_dir = "plots_dir",
    device = "png", dpi = 300, hwidth = 30, hheight = 18, width = 45, height = 15,
    bwidth = 1000000, facet_column = 5, line_size = 5,
    plot_style = c("wrap","grid"), color_panel = c("blue","red")
){
  
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  
  message("=== Step 1: Importing VCF Data ===")
  geno_data <- import_vcfdata(
    vcf_dir = vcf_dir, prefix = prefix, pattern = pattern,
    Genotypes = Genotypes, min_DP = min_DP, min_QUAL = min_QUAL,
    only_mutant = only_mutant
  )
  
  message("=== Step 2: Analyzing VCF Data ===")
  result <- analyze_vcfdata(
    geno_data = geno_data, prefix = prefix,save_results = save_results, 
    output_dir = output_dir, only_mutant = only_mutant
  )
  
  # Extract SNP datasets from results
  ant_wt    <- result$ant_wt
  ant_mt    <- result$ant_mt
  ant_wt_ems <- result$ant_wt_ems
  ant_mt_ems <- result$ant_mt_ems
  
  homo_run <- function(dt, afcol) {
    if (is.null(dt) || !nrow(dt)) return(NULL)
    if (is.null(afcol)) stop("homo_run: 'afcol' must be 'wt_AF' or 'mt_AF'.")
    window_homozygosity_compute(
      data = as.data.table(dt), af_col = afcol,
      window_size = window_size, step_size = step_size,
      rollmedian = rollmedian, nn_prop = nn_prop, af_min = af_min,
      find_intervals = find_intervals, offhold = offhold, min_vsize = min_vsize,
      use_col = use_col
    )
  }
  
  message("=== Step 3: Sliding-window homozygosity ===")
  out <- list()
  
  if (only_mutant) {
    # mutant-only; default AF column if not supplied
    if (is.null(af_col)) af_col <- "mt_AF"
    if (isTRUE(use_ems)) {
      out$homoz_mt_ems <- if (!is.null(ant_mt_ems)) homo_run(ant_mt_ems, af_col) else NULL
    } else {
      out$homoz_mt <- if (!is.null(ant_mt)) homo_run(ant_mt, af_col) else NULL
      out$homoz_mt_ems <- NULL
    }
  } else {
    # WT + MT mode
    if (isTRUE(use_ems)) {
      out$homoz_wt_ems <- if (!is.null(ant_wt_ems)) homo_run(ant_wt_ems, "wt_AF") else NULL
      out$homoz_mt_ems <- if (!is.null(ant_mt_ems)) homo_run(ant_mt_ems, "mt_AF") else NULL
    } else {
      out$homoz_wt <- if (!is.null(ant_wt)) homo_run(ant_wt, "wt_AF") else NULL
      out$homoz_mt <- if (!is.null(ant_mt)) homo_run(ant_mt, "mt_AF") else NULL
    }
  }
  
  # ---- Simple save: one Excel with interval sheets ----
  if (isTRUE(save_interval) && isTRUE(find_intervals)) {
    message("=== Step 4: Saving interval excel ===")
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      warning("save_interval=TRUE but package 'openxlsx' is not installed. Skipping Excel write.")
    } else {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      wt_name <- if (!is.null(Genotypes$wt)) Genotypes$wt else "wt"
      mt_name <- if (!is.null(Genotypes$mt)) Genotypes$mt else "mt"
      
      file_name <- if (isTRUE(only_mutant)) {
        sprintf("%s_%s_homozygosity_intervals.xlsx", prefix, mt_name)
      } else {
        sprintf("%s_%s_vs_%s_homozygosity_intervals.xlsx", prefix, wt_name, mt_name)
      }
      
      wb <- openxlsx::createWorkbook()
      
      for (label in names(out)) {
        x <- out[[label]]
        if (is.null(x)) next
        if (!is.null(x$intervals) && nrow(x$intervals) > 0L) {
          sheet <- substr(paste0(label, "_intervals"), 1, 31)  # Excel sheet name limit
          openxlsx::addWorksheet(wb, sheet)
          openxlsx::writeData(wb, sheet, x$intervals)
        }
      }
      
      xlsx_path <- file.path(output_dir, file_name)
      if (length(wb$worksheets) > 0L) {
        openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
        message("Saved intervals to: ", xlsx_path)
      } else {
        message("No interval rows to save (all empty).")
      }
    }
  }
  
  # ---- Step 5: Plotting ----
  if (isTRUE(do_plot)) {
    message("=== Step 5: Plotting ===")
    if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
    # Labels for plot titles
    mt_label <- if (!is.null(Genotypes$mt)) Genotypes$mt else "mutant"
    wt_label <- if (!is.null(Genotypes$wt)) Genotypes$wt else "wildtype"
    if (isTRUE(only_mutant)) wt_label <- NULL  # drop WT label if only_mutant
    
    # The plotting function expects a list with $results and $out
    plot_data <- list(results = result, out = out)
    
    window_homozygosity_plot(
      data = plot_data, prefix = prefix, mt = mt_label, wt = wt_label,
      plots_dir = plots_dir, device = device, dpi = dpi,
      use_col = use_col, bsa_metrics = bsa_metrics,
      plot_mode = plot_mode, plot_style = plot_style, facet_column = facet_column,
      width = width, height = height, hwidth = hwidth, hheight = hheight,
      line_size = line_size, color_panel = color_panel, bwidth = bwidth,
      only_mutant = only_mutant, af_min = af_min
    )
  }
  
  return(list(geno_data = geno_data, results = result, out = out))
}
