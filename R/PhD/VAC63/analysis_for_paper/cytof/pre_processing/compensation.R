# library(flowCore)
# library(CATALYST)


setwd("~/flo_r_biz/")
spillMat <- read.csv("~/flo_r_biz/aug2019_spillmat_final.csv", header = T, row.names = 1)

# batch1 <- read.FCS("not_comped/batch1_normalised_renamed_concat.fcs")
# batch2 <- flowCore::read.FCS("not_comped/batch2_normalised_renamed_concat.fcs")
batch3 <- flowCore::read.FCS("not_comped/batch3_normalised_renamed_concat.fcs")

# compCytof(x=batch1, y=spillMat[,-55], out_path="./comped/", method="nnls")
# system.time(CATALYST::compCytof(x=batch2, y=spillMat[,-55], out_path="./comped/", method="nnls"))
system.time(CATALYST::compCytof(x=batch3, y=spillMat[,-55], out_path="./comped/", method="nnls"))




