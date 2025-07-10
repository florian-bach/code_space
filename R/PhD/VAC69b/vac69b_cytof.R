library(CATALYST)
#library(diffcyt)
# library(dplyr)
# library(tidyr)
#library(ggplot2)

`%notin%` <- Negate(`%in%`)


setwd("~/PhD/cytof/vac69b/")

fcs <- list.files(path = "~/PhD/cytof/vac69b/", pattern = "fcs")
#fcs <- subset(fcs, !grepl(pattern = "C45", fcs))
#fcs <- subset(fcs, grepl(pattern = "Baseline", fcs))
#fcs <- subset(fcs, !grepl(pattern = "DoD", fcs))
fcs <- subset(fcs, !grepl(pattern = "ctrl", fcs))

#vac69b_flowset <- flowCore::read.flowSet(fcs)



md <- data.frame("file_name"=fcs,
            "volunteer"= substr(fcs, 1,3),
            "timepoint" = substr(fcs, nchar(fcs)-6, nchar(fcs)-4))

md$timepoint <- gsub("ine", "baseline", md$timepoint)
md$sample_id <- paste(md$volunteer, md$timepoint, sep="_")
md$batch <- ifelse(md$volunteer %in% c("v11", "v21"), "batch_1", "batch_2")


#split clustering to accomodate your laptops lousy 8gb of ram

md <- subset(md, md$timepoint == "c56")

panel <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_panel.csv", header=T)

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

pdf("./figures/vac69b_whole_blood_cluster_phenotype_ep3_ep6.pdf", height = 8, width = 9)
ctrl_whole_blod_cluster_phenotype_noh
dev.off()


#baseline, dod, ep1
#T_cells <- filterSCE(sce, k = "meta40", cluster_id %notin% c(2, 18:22, 26, 28, 30:40))

#ep3, ep6
T_cells <- filterSCE(sce, k = "meta40", cluster_id %notin% c(1:10, 12, 17, 18, 25, 32, 39))

#c56
#T_cells <- filterSCE(sce, k = "meta40", cluster_id %notin% c())

T_cells_fs <- sce2fcs(T_cells, split_by="sample_id", keep_cd=TRUE, assay="counts")

flowCore::write.flowSet(T_cells_fs, outdir = "./T_cells_only")
