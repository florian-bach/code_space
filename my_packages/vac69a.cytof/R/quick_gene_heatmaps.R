##                                             IMPORTANT !!!                                    ###
###    IF THIS DESCRIPTION IS ALTERED YOU HAVE TO RUN DOCUMENT() AT THE END OF THIS FILE        ###


#'
#' This function searches the  DoD_Baseline and T6_Baseline pairwise comparison results for all the search terms supplied
#' in the character vector; the objects DoD_Baseline_sig AND T6_Baseline_sig need to be defined in the environment for this
#' to work

#' @param  search_results character string; a list of lists containing gene names for making heatmaps
#' @param sort_by character string; the name of the comparison used to order genes in the heatmap must be in unique(all_unique_genes$file_name)
#' ("C14_Baseline" "DoD_Baseline" "C56_Baseline" "T6_Baseline"  "T6_DoD")
#' @return a named list of lists containing the data frames used for plotting. The heatmaps will be save in ./figures relative to the
#' working directory
#' @importFrom magrittr %>%

#' @export
#'
#'


quick_gene_heatmaps <- function(search_results, sort_by){


      list_of_dfs <- lapply(search_results, function(x)
        dplyr::filter(all_unique_genes, Symbol %in% x)
      )

      #return(list_of_dfs)

      for (i in 1:length(list_of_dfs)){

        plot_data <- list_of_dfs[[i]]

        plot_data <- dplyr::filter(plot_data, file_name==sort_by & baseMean!=0)

        plot_levels <- plot_data %>%
          dplyr::filter(file_name==sort_by) %>%
          dplyr::arrange(log2FoldChange) %>%
          dplyr::select(Symbol)

      plot_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x=factor(file_name, levels=c("C14_Baseline", "DoD_Baseline", "T6_DoD", "T6_Baseline")),
                                         y=factor(Symbol, levels = c(plot_levels$Symbol))))+
        #uncomment this line to add green box around DEGs
        ggplot2::geom_tile(aes(fill=log2FoldChange, width=0.92, height=0.92), size=0.4, color=ifelse(plot_data$padj<0.05, "darkgreen", "black"))+
        ggplot2::geom_tile(aes(fill=log2FoldChange))+
        ggplot2::scale_fill_gradientn(name="log2FC",
                                        values = scales::rescale(c(-1*max(plot_data$log2FoldChange, na.rm = TRUE),0, max(plot_data$log2FoldChange, na.rm = TRUE))),
                                        colors = c("#0859C6","black","#FFA500"),
                                        limits = c(-1*max(plot_data$log2FoldChange, na.rm = TRUE),
                                                   max(plot_data$log2FoldChange, na.rm = TRUE)
                                                   )
                                        # breaks=c(-2,0,2,4,6)
                                      )+
        #ggplot2::scale_fill_gradient2(low = "#0859C6", mid = "black", high = "#FFA500", midpoint = 0, limits=c(-2, 7), breaks=c(-2,0,2,4,6))+
        ggplot2::theme_void()+
        ggplot2::ggtitle(paste(names(list_of_dfs)[i], "\n", sep=''), )+
        ggplot2::guides(fill=ggplot2::guide_colorbar(nbin=30))+
        #ggplot2::coord_fixed(ratio = 1)+
        ggplot2::theme(
          axis.text = ggplot2::element_text(hjust=1, angle=45, size=7, vjust=1),
          #axis.text.y = ggplot2::element_text(angle=45, hjust = 1, vjust = 1, size = 9),
          plot.title = ggplot2::element_text(hjust=0.5),
          legend.title = ggplot2::element_text(),
          plot.margin=margin(0,0,1,1),
          legend.position = "right"
          #legend.box.margin=margin(0,0,0,0)
          )

      ggplot2::ggsave(paste("./figures/quick_heatmaps/", names(list_of_dfs)[i], ".png", sep=''), width=8, height=2, plot_plot)
      return(plot_plot)
      }


}
