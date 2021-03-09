
##                                             IMPORTANT !!!                                    ###
###    IF THIS DESCRIPTION IS ALTERED YOU HAVE TO RUN DOCUMENT() AT THE END OF THIS FILE        ###


#'
#' This function makes boxplots of the frequency or cell count of each cluster that's been determined to be differentially abundant using diffcyt.
#'
#' @param diffcyt_result; result of \code{diffcyt()} in differential abundance mode
#' @param sce; what SingleCellExperiment was used for diffcyt?
#' @param counts; logical; if \code{TRUE}, count data is displayed, if \code{FALSE} or \code{NULL} cluster frequencies are shown
#' @param FDR; numeric; false discovery rate used to correct for multiple testing
#' @param logFC; numeric; symmetrical cutoff for change to be considered relevant both up and down
#' @return a ggplot2 object
#'
#' @export
#'


  # diffcyt_result=da_dod
  # sce=merged_daf
  # FDR=0.05
  # logFC=1


diffcyt_boxplot <- function(diffcyt_result, sce, counts=NULL, FDR, logFC=NULL){

# make a table that contains the output of toptable: cluster names, cell counts and p values of every sample
results_table <- data.frame(diffcyt::topTable(diffcyt_result, all=T, show_counts = T))

# make a different table that has the edgeR output including fold changes; we'll use this to make a list
# of clusters that fulfil our logFC requirement, before removing those with high p values
all_clusters <- SummarizedExperiment::rowData(diffcyt_result$res)

changing_clusters <- subset(all_clusters$cluster_id, all_clusters$logFC > logFC | all_clusters$logFC < -logFC)

results_cut <- subset(results_table, results_table$cluster_id %in% changing_clusters)


sig_results <-  dplyr::filter(results_cut, results_cut$p_adj < FDR)

long_sig_results <- tidyr::gather(sig_results, sample_id, count, colnames(sig_results)[4:ncol(sig_results)])
long_sig_results$volunteer <- stringr::str_match(long_sig_results$sample_id, "v[0-9]*")[, 1]
long_sig_results$timepoint <- substr(long_sig_results$sample_id, 12,nchar(long_sig_results$sample_id))

long_sig_results$sample_id <- gsub("counts_", "", long_sig_results$sample_id)
cell_counts <- CATALYST::n_cells(sce)
long_sig_results$frequency <- long_sig_results$count / cell_counts[long_sig_results$sample_id] *100

if (isTRUE(counts))
  long_sig_results$ydim <- long_sig_results$count

else
  long_sig_results$ydim <-long_sig_results$frequency

fraction <- ifelse(isTRUE(counts), 1, 0.01)

if (isTRUE(counts))
  y_scale <- ggplot2::scale_y_continuous(name = "Number of Cells in Sample")

else
  y_scale <- ggplot2::scale_y_continuous(name = "Fraction of CD3+ T cells", label=scales::label_percent())



# long_sig_results$ydim <- ifelse(isTRUE(counts), long_sig_results$count, long_sig_results$frequency); head(long_sig_results$ydim)
time_col=colorspace::sequential_hcl(5, palette = "Purple Yellow")
sig_plot <- ggplot2::ggplot(long_sig_results, ggplot2::aes(x=factor(timepoint), y=ydim*fraction))+
               ggplot2::geom_boxplot(ggplot2::aes(fill=timepoint), outlier.shape = NA)+
               ggplot2::geom_point(ggplot2::aes(shape=volunteer))+
               ggplot2::facet_wrap(~cluster_id, scales = "free", ncol=7, labeller=ggplot2::label_wrap_gen(width = 20))+
               ggplot2::theme_minimal()+
               ggplot2::scale_x_discrete(name = "Timepoint")+
               y_scale+
               ggplot2::scale_fill_manual(values = c("Baseline"=time_col[4], "C10"=time_col[3], "Diagnosis"=time_col[2], "DoD"=time_col[2], "T6"=time_col[1]))+
               ggplot2::theme(axis.title = ggplot2::element_text(size=10),
                              legend.title = ggplot2::element_blank(),
                              axis.text.x = ggplot2::element_text(size=8, angle = 45),
                              strip.text = ggplot2::element_text(size=9))

return(sig_plot)

}
