# source("https://bioconductor.org/biocLite.R")
# biocLite("CATALYST")
library(flowCore)
library(premessa)

setwd("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/post_umap/")

list_of_files <- list.files(path=".", pattern="*.fcs")

t6_files <- list_of_files[grep(pattern = "T+6", list_of_files)]

premessa::concatenate_fcs_files(files.list=list_of_files, output.file = "./concat/all.fcs")

