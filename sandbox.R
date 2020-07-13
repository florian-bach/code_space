# DIDNT WORK ##

isotope_list2 <- c(CATALYST::isotope_list, list(BCKG=190))


batch1 <- flowCore::read.FCS("~/PhD/cytof/vac63c/normalised_renamed/batch1_normalised_renamed_concat.fcs")

compCytof(x=batch1,
          y=spillMat,
          out_path="~/PhD/cytof/vac63c/normalised_renamed_comped/",
          method="nnls",
          isotope_list = isotope_list2)

####

batch_fs <- flowCore::read.flowSet("~/PhD/cytof/vac63c/normalised_renamed/batch1_normalised_renamed_concat.fcs")

compCytof(x=batch_fs,
          y=spillMat,
          out_path="~/PhD/cytof/vac63c/normalised_renamed_comped/",
          method="nnls",
          isotope_list = isotope_list2)

compy <- function(fs){CATALYST::compCytof(x=fs,
                                          y=spillMat,
                                          out_path="~/PhD/cytof/vac63c/normalised_renamed_comped/",
                                          method="nnls",
                                          isotope_list = isotope_list2)}

flowCore::fsApply(x = batch_fs, FUN = compy)



test <- read.flowSet(files=c("~/PhD/cytof/vac63c/normalised_renamed/debarcoded_not_comped/301_Baseline.fcs", 
                             "~/PhD/cytof/vac63c/normalised_renamed/debarcoded_not_comped/301_C45.fcs"))

comp <- function(ff){compCytof(x=ff, y=spillMat, method, isotope_list=custom_isotope_list, out_path="python")}
fsApply(test, comp)



###

compCytof(x="~/PhD/cytof/vac63c/normalised_renamed/batch1_normalised_renamed_concat.fcs",
          y=spillMat,
          out_path="~/PhD/cytof/vac63c/normalised_renamed_comped/",
          method="nnls",
          isotope_list = isotope_list2)



custom_isotope_list <- c(CATALYST::isotope_list, list(BCKG=190))

fs <- read.FCS("~/PhD/cytof/vac63c/normalised_renamed/batch1_normalised_renamed_concat.fcs")

compCytof(x="~/PhD/cytof/vac63c/normalised_renamed/batch1_normalised_renamed_concat.fcs",
          y=spillMat,
          out_path="~/PhD/cytof/vac63c/normalised_renamed_comped/",
          method="nnls",
          isotope_list = isotope_list2)


all(as.numeric(gsub("[a-zA-Z ]", "", rownames(spillMat))) %in% unlist(isotope_list2))



batch1 <- flowCore::read.FCS("~/PhD/cytof/vac63c/normalised_renamed/batch1_normalised_renamed_concat.fcs")
spillMat <- read.csv("~/PhD/cytof/beads_and_tuning_reports/aug2019/aug2019_spillmat_final.csv", header = T, row.names = 1)

compCytof(x=batch1,
          y=spillMat,
          out_path="~/PhD/cytof/vac63c/normalised_renamed_comped/",
          method="nnls",
          isotope_list = c(CATALYST::isotope_list, list(BCKG=190)))


# isotope_list2 <- c(CATALYST::isotope_list, list(BCKG=190))
# > all(as.numeric(gsub("[a-zA-Z ]", "", colnames(spillMat))) %in% unlist(isotope_list2))
# [1] TRUE
# > all(as.numeric(gsub("[a-zA-Z ]", "", rownames(spillMat))) %in% unlist(isotope_list2))
# [1] TRUE


###############

library(CATALYST)
library(flowCore)



custom_isotope_list <- c(CATALYST::isotope_list,list(BCKG=190))

spillMat <- read.csv("~/PhD/cytof/beads_and_tuning_reports/feb2019/florian_feb2019_spillmat.csv", header = T, row.names = 1)


setwd("~/PhD/cytof/vac63c/normalised_renamed/debarcoded_not_comped/")
vac63c_premessa_table <- premessa::read_parameters(list.files(".", pattern="*.fcs")[1])
vac63c_premessa_table$Remove <- ifelse(vac63c_premessa_table[,1]=="190BCKG", TRUE, FALSE)
vac63c_premessa_table$Parameter <- rownames(vac63c_premessa_table)

premessa::rename_parameters_in_files(".", "renamed", vac63c_premessa_table)

#works??????? piece of shit...
test <- read.FCS("~/PhD/cytof/vac63c/normalised_renamed/batch1_normalised_renamed_concat.fcs")
comped_nnls <- compCytof(x=test, y=spillMat[,-55], method="nnls", out_path="~/PhD/cytof/vac63c/normalised_renamed_comped/")
write.FCS(comped_nnls, "~/PhD/cytof/vac63c/normalised_renamed_comped/batch1_normalised_renamed_comped_concat.fcs")



