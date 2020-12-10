
##                                             IMPORTANT !!!                                    ###
###    IF THIS DESCRIPTION IS ALTERED YOU HAVE TO RUN DOCUMENT() AT THE END OF THIS FILE        ###


#'
#' This function makes boxplots of the frequency or cell count of each cluster that's been determined to be differentially abundant using diffcyt.
#'
#' @param diffcyt_result; result of \code{diffcyt()} in differential abundance mode
#' @param sce; what SingleCellExperiment was used for diffcyt?
#' @param counts; logical; if \code{TRUE}, count data is displayed, if \code{FALSE} or \code{NULL} cluster frequencies are shown
#' @param FDR; numeric; false discovery rate used to correct for multiple testing
#' @return a ggplot2 object
#'
#' @export
#'





diffcyt_boxplot <- function(diffcyt_result, sce, counts=NULL, FDR){


results_table <- data.frame(diffcyt::topTable(diffcyt_result, all=T, show_counts = T))
sig_results <-  dplyr::filter(results_table, results_table$p_adj < FDR)

long_sig_results <- tidyr::gather(sig_results, sample_id, count, colnames(sig_results)[4:ncol(sig_results)])
long_sig_results$volunteer <- stringr::str_match(long_sig_results$sample_id, "V[0-9]*")[, 1]
long_sig_results$timepoint <- substr(long_sig_results$sample_id, 12,nchar(long_sig_results$sample_id))

long_sig_results$sample_id <- gsub("counts_", "", long_sig_results$sample_id)
cell_counts <- CATALYST::n_cells(sce)
long_sig_results$frequency <- long_sig_results$count / cell_counts[long_sig_results$sample_id] *100

if (isTRUE(counts))
  long_sig_results$ydim <- long_sig_results$count
else
  long_sig_results$ydim <-long_sig_results$frequency


# long_sig_results$ydim <- ifelse(isTRUE(counts), long_sig_results$count, long_sig_results$frequency); head(long_sig_results$ydim)

sig_plot <- ggplot2::ggplot(long_sig_results, ggplot2::aes(x=factor(timepoint), y=ydim))+
               ggplot2::geom_boxplot(ggplot2::aes(fill=timepoint))+
               ggplot2::geom_point(ggplot2::aes(shape=volunteer))+
               ggplot2::facet_wrap(~cluster_id, scales = "free", ncol=5, labeller=ggplot2::label_wrap_gen())+
               ggplot2::theme_minimal()+
               ggplot2::theme(axis.title = ggplot2::element_blank(),
                              legend.title = ggplot2::element_blank())

return(sig_plot)

}
