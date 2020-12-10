##                                             IMPORTANT !!!                                    ###
###    IF THIS DESCRIPTION IS ALTERED YOU HAVE TO RUN DOCUMENT() AT THE END OF THIS FILE        ###


#'
#' This function simplifies and streamlines reading the full vac69a cytof dataset into R
#' It relies on the packages flowCore, CATALYST and deplyr and their dependcies (e.g. SingleCellExperiment)
#'
#' The working directory needs to contain all of the following:
#' 1) The fcs files to be read into the flowSet
#' 2) A csv file containing metadata ("file_name", "sample_id", "timepoint", "batch", "volunteer")
#' 3) A cvs file containing the panel ("fcs_colname", "antigen", "marker_class")
#' 4) A csv file containing the names of markers to be considered for clustering
#'
#' The seed for clustering is set to 123, 100 clusters, 50 metaclusters
#' The list of markers to consider for clustering is located in refined_markers.csv in the same directory as the fcs files
#'
#' @param path_to_directory character string; path to the directory containing metadata, panel, fcs files and clustering marker names
#'
#' @return a SingleCellExperiment containing all fcs files in the directory & FlowSOM clustering results
#'
#' @examples
#'
#'# you can call the resulting variable anything other than sce (sce is the internal name) #
#'daf <- read_full("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")
#'
#'@export
#'
#'
#'
#'

#path_to_directory <- "~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/"

read_full <- function(path_to_directory){


  ### read in metadata etc. ####
  # "~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/"
  #
  # character.check <- is.character(path_to_directory)
  # space.check <- length(path_to_directory) == 1
  #
  # if(character.check==FALSE) stop("path_to_directory isn't a string")
  # if(space.check==FALSE) stop("path_to_directory contains spaces")

  working_directory <- path_to_directory

  md <- read.csv(paste(working_directory, "meta_data.csv", sep=''), header=T, stringsAsFactors = F)
  md$timepoint <- factor(md$timepoint, levels = c("Baseline", "C8", "C10", "C12", "DoD", "T6"))
  md <- md[with(md, order(md$volunteer, md$timepoint)),]
  md$file_name <- paste(working_directory, md$file_name, sep='')

  #read in panel
  panel <- read.csv(paste(working_directory, "VAC69_PANEL.CSV", sep=''), header = T, stringsAsFactors = F)
  colnames(panel)[2] <- "marker_name"

  ### read in flowfiles using flowCore
  vac69a <- flowCore::read.flowSet(md$file_name)


  ### CATALYST ####
  #construct daFrame #
  ## md has to have particular properties: file_name=NAMED LIST (chr), ID and everything else in factors
  sce <- CATALYST::prepData(vac69a, panel, md,

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
