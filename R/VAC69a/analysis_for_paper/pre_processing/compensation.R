library(flowCore)
library(CATALYST)

#fix fcs file headers

# library("cytofCore")
# cytofCore.updatePanel(templateFile="~/PhD/cytof/vac69a/reprocessed/panel_update.csv", fcsFolder="~/PhD/cytof/beads_and_tuning_reports/feb2019/")



# get single-stained control samples
ss_exp <- read.FCS("~/PhD/cytof/beads_and_tuning_reports/feb2019/relabeled/beads_feb_2019.fcs")
# specify mass channels stained for

custom_isotope_list <- c(CATALYST::isotope_list,list(BCKG=190))
# custom_isotope_list$Cd <- custom_isotope_list$Cd[-3]


#bc_ms <-  c(89, 102, 104:106, 108, 110:116, 120, 127, 131, 133, 138, 140:156, 158:176, 190:195, 198, 209, 209)
bc_ms <-  c(89, 114, 115, 141:156, 158:176, 198, 209)
# debarcode
re <- assignPrelim(x=ss_exp, y=bc_ms, verbose=TRUE)
re <- estCutoffs(x=re)
re <- applyCutoffs(x=re)
# compute spillover matrix
spillMat <- computeSpillmat(x=re)

write.csv(spillMat, "~/PhD/cytof/beads_and_tuning_reports/feb2019/feb2019_spillmat.csv")


vac69a <- read.flowSet(path = "~/PhD/cytof/vac69a/reprocessed/relabeled/", pattern = "*fcs")
spillMat <- read.csv("~/PhD/cytof/beads_and_tuning_reports/feb2019/feb2019_spillmat.csv", header = T, row.names = 1)

test <- read.FCS("~/PhD/cytof/vac69a/reprocessed/relabeled/V02_DoD.fcs")
comped_nnls <- compCytof(x=test, y=spillMat, method="nnls", isotope_list=custom_isotope_list, out_path="~/PhD/cytof/vac69a/reprocessed/reprocessed_comped/")






vac69a <- read.flowSet(path = "~/PhD/cytof/vac69a/reprocessed/relabeled/", pattern = "*fcs")

comp <- function(ff){compCytof(x=ff, y=spillMat, method, isotope_list=custom_isotope_list, out_path="~/PhD/cytof/vac69a/reprocessed/reprocessed_flow_comped/")}
fsApply(vac69a, comp)

spillMat <- read.csv("~/PhD/cytof/beads_and_tuning_reports/feb2019/feb2019_spillmat.csv", header = T, row.names = 1)
batch1 <- read.FCS("/media/flobuntu/data/cytof_imd/12022019/IMD/normalised/normalised_concat.fcs")

compCytof(x=batch1, y=spillMat, out_path="/media/flobuntu/data/cytof_imd/12022019/IMD/normalised/", method="nnls", isotope_list=custom_isotope_list)




