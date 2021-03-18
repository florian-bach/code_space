##                                             IMPORTANT !!!                                    ###
###    IF THIS DESCRIPTION IS ALTERED YOU HAVE TO RUN DOCUMENT() AT THE END OF THIS FILE        ###


#'
#' This function makes pretty UMAP plots. It takes a dataframe as input and uses the UMAP2 and UMAP1 columns as X and Y coordinates and any
#' other named column as the coloring variable. It can facet by one variable.
#'
#' It relies on the ggplot2 package
#'
#' @param df; \code{data.frame} containing cells as rows and features as columns
#' @param color_by; character string; the coloring channel; has to be one of \code{colnames(df)}
#' @param facet_by; character string; the facetting variable; has to be NULL or one of \code{colnames(df)}
#' @return a ggplot2 object
#'
#' @examples
#' sce <- read_small("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/", 500)
#' df <- prep_sce_for_ggplot(sce)
#' cd25_plot <- flo_umap(df, "CD25", "timepoint")
#' @export
#'

vac63c_flo_umap <- function (df, color_by, facet_by = NULL, only_show = NULL){

  if (is.null(facet_by))
    assign("facet_title", ggplot2::element_blank())

  if (!is.null(facet_by))
    assign("facet_title", ggplot2::element_text())

  if (!is.null(only_show))
    assign("df", subset(df, df$timepoint == only_show))


  data <- df[, c("UMAP1", "UMAP2", color_by, facet_by)]
  colnames(data)[3] <- "marker"


  # inferno_mega_lite <- c("#000004", "#8A2267", "#EF802B", "#FFEC89",
  #                        "#FCFFA4")
  UMAP_theme <- ggplot2::theme_minimal()+
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   #legend.position = "none",
                   axis.text = ggplot2::element_blank(),
                   strip.text = facet_title,
                   plot.title = element_text(hjust=0.5, size=8)
    )

  plt <- ggplot2::ggplot(data, ggplot2::aes(x = UMAP2, y = UMAP1, color = marker))+
    ggplot2::guides(colour = guide_colorbar(title="Intensity", barwidth=0.5, barheight=3))+
    ggplot2::geom_point(shape = ".", alpha=0.5) +
    viridis::scale_color_viridis(option="B")+
    UMAP_theme+
    ggplot2::facet_wrap(facet_by, ncol=4)+
    ggplot2::ggtitle(color_by)+
    ggplot2::coord_cartesian(xlim=c(-12.5, 10), ylim=c(-11.2, 9.3))
}


