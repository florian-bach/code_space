library(CATALYST)
#library(diffcyt)
# library(dplyr)
# library(tidyr)
#library(ggplot2)

`%notin%` <- Negate(`%in%`)


setwd("~/PhD/cytof/vac69b/T_cells_only/comped/recomped")
# 
# fcs <- list.files(path = "~/PhD/cytof/vac69b/T_cells_only", pattern = "fcs")
# 
# #vac69b_flowset <- flowCore::read.flowSet(fcs)
# 
# 
# 
# md <- data.frame("file_name"=fcs,
#                  "volunteer"= substr(fcs, 1,3),
#                  "timepoint" = substr(fcs, nchar(fcs)-6, nchar(fcs)-4))
# 
# md$timepoint <- gsub("ine", "baseline", md$timepoint)
# md$sample_id <- paste(md$volunteer, md$timepoint, sep="_")
# md$batch <- ifelse(md$volunteer %in% c("v11", "v21"), "batch_1", "batch_2")
# 
# write.csv(md, "metadata.csv", row.names = F)

#split clustering to accomodate your laptops lousy 8gb of ram

md <- read.csv("~/PhD/cytof/vac69b/T_cells_only/metadata.csv")
md <- subset(md, md$timepoint %in% c("baseline", "dod", "ep6", "c56"))


# the internal panel has lost the Gaussian channels and acquired the cluster_id channel-
# we need to make the internal panel match with the one we pass to CATALYST so let's quickly fix that
#internal_panel <- colnames(flowCore::read.FCS(md$file_name[1]))

panel <- read.csv("~/PhD/cytof/vac69b/T_cells_only/vac69b_panel.csv", header=T)
panel[63, 1]="icd"
# # this neatly graps all the channels associated with metals i.e. everything but gaussian & time
# panel <- subset(panel, grepl("Di", panel$fcs_colname))
# panel <- rbind(panel, c("cluster_id", "cluster_id", "none"))
# 
# write.csv(panel, "vac69b_panel.csv", row.names = F)

sce <- prepData(x=md$file_name, panel, md, md_cols =
                  list(file = "file_name", id = "sample_id", factors = c("batch", "volunteer", "timepoint")))


non_t_cell_markers <- c("CD45", "CD3",  "CD14", "CD16")

refined_markers <- c("CD4",
                     "CD8",
                     "Vd2",
                     "Va72",
                     #"CXCR5",
                     "CD38",
                     #"CD69",
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
                     #"TCRgd",
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



set.seed(123);sce <- CATALYST::cluster(sce, features = refined_markers, xdim = 10, ydim = 10, maxK = 45)

#t_cell_phenoh <- plotExprHeatmap(x = sce, by = "cluster", row_clust = TRUE, col_clust = FALSE, k = "meta50", bars = TRUE, features = all_markers)
t_cell_pheno <- plotExprHeatmap(x = sce, by = "cluster", row_clust = FALSE, col_clust = FALSE, k = "meta45", bars = TRUE, features = all_markers)

#ctrl_whole_blod_cluster_phenotype <- plotExprHeatmap(x = sce, by = "cluster", row_clust = TRUE, col_clust = FALSE, k = "meta40", bars = TRUE, features = all_markers)

# pdf("./figures/t_cell_phenoh.pdf", height = 8, width = 9)
# t_cell_phenoh
# dev.off()

pdf("./figures/t_cell_pheno.pdf", height = 8, width = 9)
t_cell_pheno
dev.off()


meta50_table <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/prim_ter_50_merge.csv")
merged_daf <- CATALYST::mergeClusters(sce, k = "meta50", table = meta45_table, id = "flo_merge", overwrite = TRUE)




test <- c("CD38", "FoxP3")
plotScatter(sce, chs = c("Tb159Di", "Gd160Di"))






#### gate again to get rid of doublets??? ####


library(openCyto)
library(ggcyto)
# library(dplyr)
# library(tidyr)
#library(ggplot2)
# 
`%notin%` <- Negate(`%in%`)
# 
# 

setwd("~/PhD/cytof/vac69b/T_cells_only/comped/")


md <- read.csv("~/PhD/cytof/vac69b/T_cells_only/metadata.csv")
md <- subset(md, md$timepoint %in% c("baseline", "dod", "ep6", "c56"))


vac69b_flowset <- flowCore::read.flowSet(md$file_name)

vac69b_gs <- GatingSet(vac69b_flowset)


#define DNA channels
panel <- read.csv("~/PhD/cytof/vac69b/T_cells_only/vac69b_panel.csv", header=T)
dna <- grep("^Ir", panel$fcs_colname, value = TRUE)



gs_add_gating_method(vac69b_gs,
                     alias = "cells",
                     pop = "+", parent = "root",
                     dims = "Ir191Di,Ce140Di",
                     gating_method = "flowClust.2d",
                     gating_args = "K=1,quantile=0.93,target=c(5,5)")



gs_add_gating_method(vac69b_gs,
                     alias = "singlets",
                     pop = "+", parent = "cells",
                     dims = paste(dna, collapse = ","),
                     gating_method = "flowClust.2d",
                     gating_args = "K=1,quantile=0.93,target=c(5,5)")

df <- gs_pop_get_stats(vac69b_gs,
                       type = "percent",
                       nodes = c("cells", "singlets"))

df


all <- ggcyto(vac69b_gs, aes_string(dna[1], dna[2]))+
  geom_hex(bins = 100)+
  geom_gate("singlets")



fs <- gs_pop_get_data(vac69b_gs, "/cells/singlets") # get data from ’GatingSet’
write.flowSet(x = fs, outdir = "./recomped")
