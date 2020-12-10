##                                             IMPORTANT !!!                                    ###
###    IF THIS DESCRIPTION IS ALTERED YOU HAVE TO RUN DOCUMENT() AT THE END OF THIS FILE        ###


#'
#' This function simplifies and streamlines reading a subsampled version of the vac69a cytof dataset into R
#' It relies on the packages flowCore, CATALYST and deplyr and their dependcies (e.g. SingleCellExperiment)
#'
#' The working directory needs to contain all of the following:
#' \itemize{
#' \item The fcs files to be read into the flowSet
#' \item A csv file containing metadata ("file_name", "sample_id", "timepoint", "batch", "volunteer")
#' \item A cvs file containing the panel ("fcs_colname", "antigen", "marker_class")
#' \item A csv file containing the names of markers to be considered for clustering
#' }
#' The seed for clustering is set to 123, 100 clusters, 50 metaclusters
#' The list of markers to consider for clustering is located in refined_markers.csv in the same directory as the fcs files
#'
#' @param path_to_direcory character string; path to the directory containing metadata, panel, fcs files and clustering marker names
#' @param proportional logical; should relative file sizes be maintained
#' @param event_number integer; sampling floor (if proportional = T), or sampling ceiling (if proportional = F)
#' @return a SingleCellExperiment containing all downsampled fcs files in the directory & FlowSOM clustering results
#'
#' @examples
#'
#'# you can call the resulting variable anything other than sce (sce is the internal name) #
#'daf <- read_small("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/", proportional=TRUE, event_number=3000)
#'
#'
#'
#' @export
#'
#'


read_small <- function(path_to_directory, proportional = NULL, event_number){


  ### read in metadata etc. ####
  #path_to_directory <- "~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/"
  #
  # character.check <- is.character(path_to_direcory)
  # space.check <- length(path_to_direcory) == 1
  #
  # if(character.check==FALSE) stop("path_to_direcory isn't a string")
  # if(space.check==FALSE) stop("path_to_direcory contains spaces")

  working_directory <- path_to_directory

  md <- read.csv(paste(working_directory, "meta_data.csv", sep=''), header=T, stringsAsFactors = F)
  md$timepoint <- factor(md$timepoint, levels = c("Baseline", "C8", "C10", "C12", "Diagnosis", "T6"))
  md <- md[with(md, order(md$volunteer, md$timepoint)),]
  md$file_name <- paste(working_directory, md$file_name, sep='')

  #read in panel
  panel <- read.csv(paste(working_directory, "VAC69_PANEL.CSV", sep=''), header = T, stringsAsFactors = F)
  colnames(panel)[2] <- "marker_name"

  ### read in flowfiles using flowCore
  vac69a <- flowCore::read.flowSet(md$file_name)


  #sample.int takes a sample of the specified size from the elements of x using either with or without replacement.
  # set.seed(1234); smaller_vac69a <- fsApply(vac69a, function(ff) {
  #   idx <- sample.int(nrow(ff), min(sampling_ceiling, nrow(ff)))
  #   ff[idx,]  # alt. ff[order(idx),]
  # })

  if(isTRUE(proportional))
    assign("downsample", function(fs, event_number){
      flowCore::fsApply(fs, function(ff){
      idx <- sample.int(nrow(ff), nrow(ff)/min(flowCore::fsApply(fs, nrow))*event_number)
      ff[idx,]
    })
      })

  if(xor(isFALSE(proportional), is.null(proportional)))
    assign("downsample", function(fs, event_number){
      flowCore::fsApply(fs, function(ff){
      idx <- sample.int(nrow(ff), min(event_number, nrow(ff)))
      ff[idx,]
    })
      })


set.seed(1234); smaller_vac69a <- downsample(vac69a, event_number)



  ### CATALYST ####
  #construct daFrame #
  ## md has to have particular properties: file_name=NAMED LIST (chr), ID and everything else in factors
  sce <- CATALYST::prepData(smaller_vac69a, panel, md,

                            md_cols = list(file = "file_name",
                                           id = "sample_id",
                                           factors = c("timepoint", "batch", "volunteer")
                            ),

                            panel_cols = list(channel = "fcs_colname",
                                              antigen = "marker_name",
                                              class = "marker_class"
                            )
  )



  # refined_markers <- c("CD4",
  #                      "CD8",
  #                      "Vd2",
  #                      "Va72",
  #                      "CD38",
  #                      "HLADR",
  #                      "ICOS",
  #                      "CD28",
  #                      "PD1",
  #                      #"TIM3",
  #                      "CD95",
  #                      "BCL2",
  #                      "CD27",
  #                      "Perforin",
  #                      "GZB",
  #                      "CX3CR1",
  #                      "Tbet",
  #                      "CTLA4",
  #                      "Ki67",
  #                      "CD127",
  #                      #"IntegrinB7",
  #                      #"CD56",
  #                      #"CD16",
  #                      "CD161",
  #                      #"CD49d",
  #                      #"CD103",
  #                      "CD25",
  #                      "FoxP3",
  #                      "CD39",
  #                      "CLA",
  #                      #"CXCR5",
  #                      "CD57",
  #                      "CD45RA",
  #                      "CD45RO",
  #                      "CCR7")

  #write.csv(data.frame(refined_markers), paste(working_directory, "refined_markers.csv", sep = ''),  row.names = F)
  refined_markers <- read.csv(paste(working_directory, "refined_markers.csv", sep = ''), stringsAsFactors = F)

  # clustering ####
  set.seed(123);sce <- CATALYST::cluster(sce, features = refined_markers[,1], xdim = 10, ydim = 10, maxK = 50)

  return(sce)
}

#working_directory <- "~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/"
#daf <- import_vac69a_and_cluster("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")
