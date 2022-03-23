
##                                             IMPORTANT !!!                                    ###
###    IF THIS DESCRIPTION IS ALTERED YOU HAVE TO RUN DOCUMENT() AT THE END OF THIS FILE        ###


#'
#' This function performs hard downsampling of a flowset
#'
#' It randomly selects a wiebke-defined number of cells from each fcs file in a flowset.
#'
#' @param fs; a flowCore flowset to be downsampled
#' @param event_number; a positive integer smaller than or equal to \code{min(fsApply(fs, nrow))}` indicating what number of cells to retain in the downsampled flowset
#' @return a (smaller) flowCore flowset
#'
#' @note
#'
#' As the downsampling is performed randomly, running it multiple times will give slightly different results. To reproducibly have the same output include manually set the seed, as shown in the example below.
#' @examples
#'
#' # manually set the seed to have reproducible results; use whatever integer you like.
#' set.seed(1234)
#'
#' # keep 5000 events from each file
#' smaller_vac69a <- hard.downsample(vac69a, 5000)
#'
#' @seealso
#'
#' \code{\link{proportional.downsample}}
#'
#' @export
#'


hard.downsample <- function(fs, event_number){
  flowCore::fsApply(fs, function(ff){
    idx <- sample.int(nrow(ff), min(event_number, nrow(ff)))
    ff[idx,]
  })
}
