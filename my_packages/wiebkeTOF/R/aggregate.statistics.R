##                                             IMPORTANT !!!                                    ###
###    IF THIS DESCRIPTION IS ALTERED YOU HAVE TO RUN DOCUMENT() AT THE END OF THIS FILE        ###
#'
#' This function splits a SingleCellExperiment along come conditions and computes summary statistics. This is useful when we want to look at the median expression of some markers on every cluster at a specific timepoint, for example. It almost entirely reuses code working in the background of `CATALYST::plotExprHeatmap()`
#'
#'
#' @param x; a SingleCellDataset
#' @param markers; character vector: which markers do you wish to compute summary statistics for?
#' @param by; what condition to split the SingleCellExperiment by? has to be one or more of `colnames(colData(sce))`
#' @param clustering; which clustering should be used to split cells? has to be in of `names(cluster_codes(sce))`
#' @param fun; which summary statistc to calculate? has to be one of `"median"`, `"mean"` or `"sum"`
#'
#' @return a list of matrices (one for each cluster), where each column is a condition (e.g. timepoint), each row a cluster, and each value a summary statistic
#'
#' @examples
#'
#' ms <- aggregate.statistics(x=sce, markers=refined_markers, by = c("cluster_id", "timepoint"), fun="median", clustering="meta40")
#'
#' @export


# expression per cluster ####


aggregate.statistics <- function(x, markers, by, clustering, fun = c("median", "mean", "sum")) {

  split.cells <- function(x, by, clustering) {
    stopifnot(is.character(by), by %in% colnames(colData(x)))
    cd <- data.frame(colData(x))
    cd$cluster_id <- cluster_ids(x, k=clustering)
    dt <- data.table::data.table(cd, i = seq_len(ncol(x)))
    dt_split <- split(dt, by = by, sorted = TRUE, flatten = FALSE)
    map_depth(dt_split, length(by), "i")
  }

  fun <- switch(match.arg(fun),
                median = rowMedians, mean = rowMeans, sum = rowSums)
  x_slim <- x[markers,]
  cs <- split.cells(x_slim, by, clustering)
  pb <- map_depth(cs, -1, function(i) {
    if (length(i) == 0) return(numeric(nrow(x)))
    fun(assay(x_slim, "exprs")[, i, drop = FALSE])
  })
  map_depth(pb, -2, function(u) as.matrix(data.frame(
    u, row.names = rownames(x_slim), check.names = FALSE)))
}


