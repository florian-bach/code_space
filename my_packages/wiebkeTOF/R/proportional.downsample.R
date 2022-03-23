
##                                             IMPORTANT !!!                                    ###
###    IF THIS DESCRIPTION IS ALTERED YOU HAVE TO RUN DOCUMENT() AT THE END OF THIS FILE        ###


#'
#' This function performs proportional downsampling of a flowset
#'
#' It randomly selects a wiebke-defined percentage from each fcs file in a flowset.
#'
#' @param fs; a flowCore flowset to be downsampled
#' @param fraction; number between 0 and 1 indicating what fraction of cells to retain in the downsampled setwhat SingleCellExperiment was used for diffcyt?
#' @return a flowCore flowset
#'
#' @note
#'
#' As the downsampling is performed randomly, running it multiple times will give slightly different results. To reproducibly have the same output include manually set the seed, as shown in the example below.
#' @examples
#'
#' # manually set the seed to have reproducible results; use whatever integer you like.
#' set.seed(1234)
#'
#' # keep half of all events
#' smaller_vac69a <- proportional.downsample(vac69a, 0.5)
#'
#' @seealso
#'
#' \code{\link{hard.downsample}}
#'
#' @export
#'

proportional.downsample <- function(fs, fraction){
  flowCore::fsApply(fs, function(ff){
    idx <- sample.int(n=nrow(ff), size=nrow(ff)*fraction)
    ff[idx,]})
}
