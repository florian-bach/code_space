library(CATALYST)
library(diffcyt)
# library(dplyr)
# library(tidyr)
#library(ggplot2)


setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/debarcoded/")

fcs <- list.files(pattern = "fcs")
#fcs <- subset(fcs, !grepl(pattern = "C45", fcs))
#fcs <- subset(fcs, grepl(pattern = "Baseline", fcs))
#fcs <- subset(fcs, !grepl(pattern = "DoD", fcs))
fcs <- subset(fcs, grepl(pattern = "ctrl", fcs))

vac63_flowset <- flowCore::read.flowSet(fcs)



md <- cbind("file_name"=fcs, "batch"=c("batch_1", "batch_2", "batch_3"), "sample_id"=c("batch_1", "batch_2", "batch_3"))

panel <- read.csv("vac63c_panel.csv", header=T)

sce <- prepData(vac63_flowset, panel, md, md_cols =
                  list(file = "file_name", id = "sample_id", factors = "batch"))


non_t_cell_markers <- c("CD45", "CD3",  "CD14", "CD16")

all_markers <- c(non_t_cell_markers, refined_markers)



set.seed(123);sce <- CATALYST::cluster(sce, features = refined_markers, xdim = 10, ydim = 10, maxK = 50)



plotExprHeatmap(x = sce, by = "cluster", row_clust = FALSE, k = "meta45", bars = TRUE, features = all_markers)

PBMC <- filterSCE(sce, k = "meta45", cluster_id %notin% c(4,5,43))
PBMC <- CATALYST::cluster(PBMC, features = all_markers, xdim = 10, ydim = 10, maxK = 50)

vac63c_control_pbmc_cluster_heatmap <- plotExprHeatmap(x = PBMC, by = "cluster", row_clust = TRUE, k = "meta45", bars = TRUE, features = c(all_markers, "GATA3", "RORgt"))

pdf("vac63c_control_pbmc_cluster_heatmap.pdf", height=11, width=11)
ComplexHeatmap::draw(vac63c_control_pbmc_cluster_heatmap)
dev.off()

monocytes <- filterSCE(PBMC, k = "meta45", cluster_id %in% c(1:10, 15, 22, 39))
CATALYST::plotClusterExprs(monocytes, features=non_t_cell_markers, k="meta45")

pdf("scatter_try.pdf", height=8, width=8)
plotScatter(monocytes, chs = c("CD14", "HLADR"), color_by = NULL, label = "both", facet_by = "cluster_id")
dev.off()
