##                                             IMPORTANT !!!                                    ###
###    IF THIS DESCRIPTION IS ALTERED YOU HAVE TO RUN DOCUMENT() AT THE END OF THIS FILE        ###
#'
#' This function splits a SingleCellExperiment along come conditions and computes summary statistics. This is useful when we want to look at the median expression of some markers on every cluster at a specific timepoint, for example. It almost entirely reuses code working in the background of `CATALYST::plotExprHeatmap()`
#'
#'
#' @param x; a list of matrices made with aggregate.statistics
#' @param fun; which summary statistc to calculate across clusters? has to be one of `"median"`, `"mean"` or `"sum"`
#'
#' @return a single matrix, where each column is a condition (e.g. timepoint) with an additional column for the summary statistic, each row a cluster, and each value a summary statistic calculated with `aggregate.statistics` or `simplify.aggregate`.
#'
#' @examples
#'
#' median_ms <- simplify.aggregate(ms, "median")
#' @seealso
#'
#' \code{\link{aggregate.statistics}}
#' @export


simplify.aggregate <- function(ms, fun){

ms2 <- lapply(ms, function(x)data.frame(x))
ms2 <- Map(cbind, ms2, cluster_id = names(ms))

ms3 <- do.call(rbind, ms2)
ms3$Marker <- unlist(lapply(ms2, rownames))

summary <- apply(ms3[,1:unique(sapply(ms, ncol))], MARGIN = 1, fun)
ms3 <- cbind(ms3, summary)
colnames(ms3)[ncol(ms3)] <- paste(fun, "_expression", sep='')
return(ms3)
}
