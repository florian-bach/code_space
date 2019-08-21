library(premessa)

paneleditor_GUI()
setwd("/Users/s1249052/PhD/flow_data/vac69a/20190212")

files_list <- list.files(path=".", pattern="*.fcs")

concatenate_fcs_files(files_list, "concat679.fcs")


debarcoder_GUI()
