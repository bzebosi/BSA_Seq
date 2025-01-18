# Function generates plot configuration
plot_config <- function(wt_mt, ant_wt, ant_mt, wt, mt, plot_type = "all"){
  # Define plot types and their configurations
  # 'af' for Allele Frequency, 'ed' for Euclidean Distance, 'hist' for Histograms for homozygous snps, and 'pval' for P-values
  # input argument:
  # wt_mt contains wildtype and mutant data
  # ant_wt and ant_mt contains snps unique wildtype and mutant.
  # wt and mt are wildtype and mutant labels
  # plot_type options are "af" , "ed" , "hist", "hist" and "pval"
  # 1. plot_type ="all" or plot_type = c("af", "ed", "hist", "pval"), configurations for all plot types such af, hist, ed and pval
  # 2. plot_type ="af" configurations only include the allele frequency (af) plot type
  # 3. plot_type = c("af", "ed") configurations include the allele frequency (af), Euclidean Distance (ed)  plot type
  plot_types <- list(
    af = list(
      list(data = wt_mt[wt_AF < 1], column = "wt_AF", plot_title = paste0(wt, " AF  (SNP-Index)", " using shared snps with ", " ", mt),
           y_title = paste0(wt," SNP Index"), file_suffix = "P1a_wt_AF_untransformed", ylim = NULL,lcolor = "black", lsize = 1, psize=0.5,
           is_smooth = FALSE, is_histogram=FALSE, threshold = NULL, rollmedian = 401), 
      list(data = wt_mt[wt_AF < 1], column = "wt_AF", plot_title = paste0(wt, " AF  (SNP-Index)", " using shared snps with ", " ", mt),
           y_title = paste0(wt," SNP Index"), file_suffix = "P1b_wt_AF_transformed", ylim = NULL,lcolor = "black", lsize = 1, psize=0.5,
           is_smooth = TRUE, is_histogram=FALSE, threshold = NULL, rollmedian = 0),
      list(data = wt_mt[mt_AF < 1], column = "mt_AF", plot_title = paste0(mt, " AF  (SNP-Index)", " using shared snps with ", " ", wt),
           y_title = paste0(mt," SNP Index"), file_suffix = "P2a_mt_AF_untransformed", ylim = NULL,lcolor = "black", lsize = 1, psize=0.5,
           is_smooth = FALSE, is_histogram=FALSE, threshold = NULL, rollmedian = 401), 
      list(data = wt_mt[mt_AF < 1], column = "mt_AF", plot_title = paste0(mt, " (AF / SNP-Index)", " using shared snps with ", " ", wt),
           y_title = paste0(mt," SNP Index"), file_suffix = "P2b_mt_AF_transformed", ylim = NULL,lcolor = "black", lsize = 1, psize=0.5,
           is_smooth = TRUE, is_histogram=FALSE, threshold = NULL, rollmedian = 0),
      list(data = wt_mt[AF_diff < 1], column = "AF_diff", plot_title = paste0("AF diff (Δ SNP-Index): ",  wt,  " - ",  mt),
           y_title = paste0("AF diff : ", "  ",  "(Δ SNP Index)"), file_suffix = "P3a_AF_diff_untranformed", ylim = NULL,lcolor = "black", lsize = 1, psize=0.5,
           is_smooth = FALSE, is_histogram=FALSE, threshold = NULL, rollmedian = 401), 
      list(data = wt_mt[AF_diff < 1], column = "AF_diff", plot_title = paste0("AF diff (Δ SNP-Index): ",  wt,  " - ",  mt),
           y_title = paste0("AF diff : ", "  ",  "(Δ SNP Index)"), file_suffix = "P3b_AF_diff_tranformed", ylim = NULL,lcolor = "black", lsize = 1, psize=0.5,
           is_smooth = TRUE, is_histogram=FALSE, threshold = NULL, rollmedian = 0), 
      list(data = ant_wt[wt_AF < 1], column = "wt_AF", plot_title = paste0(wt, "AF  (SNP-INDEX)", " using unique snps only "),
           y_title = paste0(wt, " SNP Index"), file_suffix = "P4a_onlywt_AF_untransformed", ylim = NULL, lcolor = "black", lsize = 1, psize=0.5,
           is_smooth = FALSE, is_histogram=FALSE, threshold = NULL, rollmedian = 401), 
      list(data = ant_wt[wt_AF < 1], column = "wt_AF", plot_title = paste0(wt, "AF  (SNP-INDEX)", " using unique snps only "),
           y_title = paste0(wt, " SNP Index"), file_suffix = "P4b_onlywt_AF_transformed", ylim = NULL, lcolor = "black", lsize = 1, psize=0.5,
           is_smooth = TRUE, is_histogram=FALSE, threshold = NULL, rollmedian = 0),
      list(data = ant_mt[mt_AF < 1], column = "mt_AF", plot_title = paste0(mt, "AF  (SNP-INDEX)", " using unique snps only "),
           y_title = paste0(mt, " SNP Index"), file_suffix = "P5a_onlymt_AF_untransformed", ylim = NULL, lcolor = "black", lsize = 1, psize=0.5,
           is_smooth = FALSE, is_histogram=FALSE, threshold = NULL, rollmedian = 401), 
      list(data = ant_mt[mt_AF < 1], column = "mt_AF", plot_title = paste0(mt, " AF ", " using unique snps only "),
           y_title = paste0(mt, " SNP Index"), file_suffix = "P5b_onlymt_AF_transformed", ylim = NULL, lcolor = "black", lsize = 1, psize=0.5,
           is_smooth = TRUE, is_histogram=FALSE, threshold = NULL, rollmedian = 0)
    ), 
    ed = list(
      list(data = wt_mt[ED > 1], column = "ED", plot_title = paste0("ED: ", wt, " vs ", mt),
           y_title = paste0("ED"), file_suffix = "P6a_ED_untranformed", ylim = NULL,lcolor = "black", lsize = 1, psize=0.5,
           is_smooth = FALSE, is_histogram=FALSE, threshold = NULL, rollmedian = 401), 
      list(data = wt_mt[ED > 1], column = "ED", plot_title = paste0("ED: ", wt, " vs ", mt),
           y_title = paste0("ED"), file_suffix = "P6b_ED_tranformed", ylim = NULL, lcolor = "black", lsize = 1, psize=0.5,
           is_smooth = TRUE, is_histogram=FALSE, threshold = NULL, rollmedian = 0),
      list(data = wt_mt[ED4 > 1], column = "ED4", plot_title = paste0("ED4: ", wt, " vs ", mt),
           y_title = paste0("ED4"), file_suffix = "P7a_ED4_untranformed", ylim = NULL,lcolor = "black", lsize = 1, psize=0.5,
           is_smooth = FALSE, is_histogram=FALSE, threshold = NULL, rollmedian = 401), 
      list(data = wt_mt[ED4 > 1], column = "ED4", plot_title = paste0("ED4: ", wt, " vs ", mt),
           y_title = paste0("ED4"), file_suffix = "P7b_ED4_tranformed", ylim = NULL, lcolor = "black", lsize = 1, psize=0.5,
           is_smooth = TRUE, is_histogram=FALSE, threshold = NULL, rollmedian = 0)
    ), 
    hist = list(
      list(data = ant_mt[mt_AF>=1], column = "mt_AF", plot_title = paste0(mt, " unique snps only "),
           y_title = paste0("Homozygous SNPs per Mb"), file_suffix = "P8a_onlymt_AF_histogram", ylim = NULL, lcolor = "black", lsize = 1, psize=0.5,
           is_smooth = FALSE, is_histogram=TRUE, threshold = NULL, rollmedian = 0),
      list(data = ant_wt[wt_AF>=1], column = "wt_AF", plot_title = paste0(wt, " unique snps only with AF >=1"),
           y_title = paste0("Homozygous SNPs per Mb"), file_suffix = "P8b_onlywt_AF_histogram", ylim = NULL, lcolor = "black", lsize = 1, psize=0.5,
           is_smooth = FALSE, is_histogram=TRUE, threshold = NULL, rollmedian = 0)
    ), 
    pval = list(
      list(data = wt_mt[AF_diff < 1], column = "log_adj.pval", plot_title = paste0("-log10(p-value): ", wt, " vs ", mt),
           y_title = "-log10(p-value)", file_suffix = "P6_log_adj.pval", ylim = NULL, lcolor = "black", lsize = 1, psize=0.5,
           is_smooth = FALSE, is_histogram=FALSE, threshold = -log10(0.05)*10, rollmedian = 0)
    )
  )
  selected_plot_types <- if (plot_type == "all") names(plot_types) else plot_type
  

  plot_list <- list()
  for (type in selected_plot_types) {
    plot_list <- c(plot_list, plot_types[[type]])
  }

  return(plot_list)
}

# Example Usage
configs <- plot_config(wt_mt, ant_wt, ant_mt, wt, mt, plot_type = "all")
print(configs)