# source("https://bioconductor.org/biocLite.R")
# biocLite("CATALYST")

setwd("C:/Users/Florian/PhD/cytof/June_2019_titration_fix_opti")

list_of_files <- list.files(path=".", pattern="*.fcs")

flo_set <- read.flowSet(list_of_files)

library(CATALYST)

concatFCS(flo_set, out_path=".")
