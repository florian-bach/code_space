vivax_panel <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/VAC69_PANEL.CSV", header=T)

falciparum_panel <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_panel.csv", header=T)


harmony_panel <- subset(falciparum_panel, falciparum_panel[,1:2]==vivax_panel[,1:2])

harmony_panel <- na.omit(harmony_panel)

setwd("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")

vivax_daf <- vac69a.cytof::read_full("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")




setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only")
#setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/debarcoded/")


md <- read.csv("vac63c_metadata.csv", header = T)

prim_files <- subset(md, md$n_infection=="First", select = "file_name")

vac63_flowset <- flowCore::read.flowSet(prim_files[,1])

md <- subset(md, md$file_name %in% prim_files[,1])
md <- md[order(md$timepoint),]



#
panel <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_panel.csv", header=T)

falciparum_daf <- prepData(vac63_flowset, panel, md, md_cols =
                         list(file = "file_name", id = "sample_id", factors = c("volunteer", "timepoint", "n_infection", "batch")))

slim_falciparum_daf <- filterSCE(falciparum_daf, )
