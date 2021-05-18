#### using FlowSOM to gate out all non T cells ####



library(CATALYST)
#library(diffcyt)
# library(dplyr)
# library(tidyr)
#library(ggplot2)

`%notin%` <- Negate(`%in%`)


setwd("~/PhD/cytof/vac69b/whole_blood_single_cells/")

fcs <- list.files(path = "~/PhD/cytof/vac69b/whole_blood_single_cells/", pattern = "fcs")
#fcs <- subset(fcs, !grepl(pattern = "C45", fcs))
#fcs <- subset(fcs, grepl(pattern = "Baseline", fcs))
#fcs <- subset(fcs, !grepl(pattern = "DoD", fcs))
fcs <- subset(fcs, !grepl(pattern = "ctrl", fcs))

#vac69b_flowset <- flowCore::read.flowSet(fcs)


# 
# md <- data.frame("file_name"=fcs,
#             "volunteer"= substr(fcs, 1,3),
#             "timepoint" = substr(fcs, nchar(fcs)-6, nchar(fcs)-4))
# 
# md$timepoint <- gsub("ine", "baseline", md$timepoint)
# md$sample_id <- paste(md$volunteer, md$timepoint, sep="_")
# md$batch <- ifelse(md$volunteer %in% c("v11", "v21"), "batch_1", "batch_2")
# 


md <- read.csv("~/PhD/cytof/vac69b/T_cells_only/metadata.csv", header=T)
#split clustering to accomodate your laptops lousy 8gb of ram

md <- subset(md, md$timepoint %in% c("baseline", "dod", "ep6", "c56") & md$volunteer %in% c("v11", "v21"))

#panel <- read.csv("~/PhD/cytof/vac69b/T_cells_only/vac69b_panel.csv", header=T)
panel <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_panel.csv")
sce <- prepData(x=md$file_name, panel, md, md_cols =
                  list(file = "file_name", id = "sample_id", factors = c("batch", "volunteer", "timepoint")))


non_t_cell_markers <- c("CD45", "CD3",  "CD14", "CD16")


refined_markers <- c("CD4",
                     "CD8",
                     "Vd2",
                     "Va72",
                     "CXCR5",
                     "CD38",
                     "CD69",
                     "HLADR",
                     "ICOS",
                     "CD28",
                     "PD1",
                     #"TIM3",
                     "CD95",
                     "BCL2",
                     "CD27",
                     "Perforin",
                     "GZB",
                     "TCRgd",
                     "Tbet",
                     "Eomes",
                     #"RORgt",
                     #"GATA3",
                     "CTLA4",
                     "Ki67",
                     "CD127",
                     "CD56",
                     #"CD16",
                     "CD161",
                     "CD49d",
                     "CD25",
                     "FoxP3",
                     "CD39",
                     "CX3CR1",
                     "CD57",
                     "CD45RA",
                     "CD45RO",
                     "CCR7")


all_markers <- c(non_t_cell_markers, refined_markers)



set.seed(123);sce <- CATALYST::cluster(sce, features = refined_markers, xdim = 10, ydim = 10, maxK = 40)

ctrl_whole_blod_cluster_phenotype_noh <- plotExprHeatmap(x = sce, by = "cluster", row_clust = FALSE, col_clust = FALSE, k = "meta40", bars = TRUE, features = all_markers)

#ctrl_whole_blod_cluster_phenotype <- plotExprHeatmap(x = sce, by = "cluster", row_clust = TRUE, col_clust = FALSE, k = "meta40", bars = TRUE, features = all_markers)

# pdf("./figures/ctrl_whole_blod_cluster_phenotype_no_hierarchy.pdf", height = 8, width = 9)
# ctrl_whole_blod_cluster_phenotype_noh
# dev.off()

pdf("./figures/vac69b_whole_blood_cluster_phenotype_v11_v21.pdf", height = 8, width = 9)
ctrl_whole_blod_cluster_phenotype_noh
dev.off()


#v11, v21, baseline, dod, ep6, c56
T_cells <- filterSCE(sce, k = "meta40", cluster_id %notin% c(1:2, 5, 8, 16, 20:21, 26:30, 34:38))

#v527, baseline, dod, ep6, c56
#T_cells <- filterSCE(sce, k = "meta40", cluster_id %notin% c(1, 9:10, 18, 20:21, 27:28, 31:40))

#c56
#T_cells <- filterSCE(sce, k = "meta40", cluster_id %notin% c(1:6, 9:13, 17:23, 37))

T_cells_fs <- sce2fcs(T_cells, split_by="sample_id", keep_cd=TRUE, assay="counts")

flowCore::write.flowSet(T_cells_fs, outdir = "~/PhD/cytof/vac69b/T_cells_only/")


#### compensation ####


#reusing the spillmat from the vac63c dataset, as no antibodies were changed
spillMat <- read.csv("~/PhD/cytof/beads_and_tuning_reports/aug2019//aug2019_spillmat_final.csv", header = T, row.names = 1)

setwd("~/PhD/cytof/vac69b/T_cells_only")

md <- read.csv("metadata.csv")
md <- subset(md, md$timepoint %in% c("baseline", "dod", "ep6", "c56"))


# the internal panel has lost the Gaussian channels and acquired the cluster_id channel-
# we need to make the internal panel match with the one we pass to CATALYST so let's quickly fix that
#internal_panel <- colnames(flowCore::read.FCS(md$file_name[1]))

panel <- read.csv("vac69b_panel.csv", header=T)

sce <- prepData(x=md$file_name, panel, md, md_cols =
                  list(file = "file_name", id = "sample_id", factors = c("batch", "volunteer", "timepoint")))





# compCytof(x=batch1, y=spillMat[,-55], out_path="./comped/", method="nnls")
# system.time(CATALYST::compCytof(x=batch2, y=spillMat[,-55], out_path="./comped/", method="nnls"))
comped_sce <- CATALYST::compCytof(x=sce,
                                sm=spillMat[,-55],
                                method="nnls")


comped_T_cells_fs <- sce2fcs(comped_sce, split_by="sample_id", keep_cd=TRUE, assay="counts")

flowCore::write.flowSet(comped_T_cells_fs, outdir = "~/PhD/cytof/vac69b/T_cells_only/comped/")



