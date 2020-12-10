##                                             IMPORTANT !!!                                    ###
###    IF THIS DESCRIPTION IS ALTERED YOU HAVE TO RUN DOCUMENT() AT THE END OF THIS FILE        ###


#'
#' This function makes data.frames when passed SingleCellExperiments. This allows making ggplot objects from scratch using sce data.
#'
#' It relies on the ggplot2 package
#'
#' @param sce; a SingleCellExperiment containing UMAP dimensions
#' @return a data.frame containing events as rows and features as columns
#'
#' @examples big_table <- prep_sce_for_ggplot(daf)
#' @export


prep_sce_for_ggplot <- function(sce){

#extract event data, metadata and transpose
df <- data.frame(t(data.frame(SummarizedExperiment::assays(sce)$exprs)))
df <- data.frame(cbind(df, SingleCellExperiment::colData(sce)))

#extract UMAP coordinates
umaps <- data.frame(SingleCellExperiment::reducedDim(sce, "UMAP"))
colnames(umaps) <- c("UMAP1", "UMAP2")

# put it together
df <- data.frame(cbind(df, umaps), stringsAsFactors = F)
return(df)
}
