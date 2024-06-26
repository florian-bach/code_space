% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_small.R
\name{read_small}
\alias{read_small}
\title{This function simplifies and streamlines reading a subsampled version of the vac69a cytof dataset into R
It relies on the packages flowCore, CATALYST and deplyr and their dependcies (e.g. SingleCellExperiment)}
\usage{
read_small(path_to_directory, proportional = NULL, event_number)
}
\arguments{
\item{proportional}{logical; should relative file sizes be maintained}

\item{event_number}{integer; sampling floor (if proportional = T), or sampling ceiling (if proportional = F)}

\item{path_to_direcory}{character string; path to the directory containing metadata, panel, fcs files and clustering marker names}
}
\value{
a SingleCellExperiment containing all downsampled fcs files in the directory & FlowSOM clustering results
}
\description{
The working directory needs to contain all of the following:
\itemize{
\item The fcs files to be read into the flowSet
\item A csv file containing metadata ("file_name", "sample_id", "timepoint", "batch", "volunteer")
\item A cvs file containing the panel ("fcs_colname", "antigen", "marker_class")
\item A csv file containing the names of markers to be considered for clustering
}
The seed for clustering is set to 123, 100 clusters, 50 metaclusters
The list of markers to consider for clustering is located in refined_markers.csv in the same directory as the fcs files
}
\examples{

# you can call the resulting variable anything other than sce (sce is the internal name) #
daf <- read_small("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/", proportional=TRUE, event_number=3000)



}
