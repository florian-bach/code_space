library(CATALYST)
library(diffcyt)
# library(dplyr)
# library(tidyr)
#library(ggplot2)

`%notin%` <- Negate(`%in%`)


setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/debarcoded/")

fcs <- list.files(path = "~/PhD/cytof/vac63c/normalised_renamed_comped/debarcoded/whole_blood_single_cells/", pattern = "fcs")
#fcs <- subset(fcs, !grepl(pattern = "C45", fcs))
#fcs <- subset(fcs, grepl(pattern = "Baseline", fcs))
#fcs <- subset(fcs, !grepl(pattern = "DoD", fcs))
fcs <- subset(fcs, grepl(pattern = "ctrl", fcs))

vac63_flowset <- flowCore::read.flowSet(fcs)



md <- cbind("file_name"=fcs, "batch"=c("batch_1", "batch_2", "batch_3"), "sample_id"=c("batch_1", "batch_2", "batch_3"))

panel <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_panel.csv", header=T)

sce <- prepData(vac63_flowset, panel, md, md_cols =
                  list(file = "file_name", id = "sample_id", factors = "batch"))


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



set.seed(123);sce <- CATALYST::cluster(sce, features = refined_markers, xdim = 10, ydim = 10, maxK = 50)

ctrl_whole_blod_cluster_phenotype_noh <- plotExprHeatmap(x = sce, by = "cluster", row_clust = FALSE, col_clust = FALSE, k = "meta40", bars = TRUE, features = all_markers)

ctrl_whole_blod_cluster_phenotype <- plotExprHeatmap(x = sce, by = "cluster", row_clust = TRUE, col_clust = FALSE, k = "meta40", bars = TRUE, features = all_markers)

pdf("./figures/ctrl_whole_blod_cluster_phenotype_no_hierarchy.pdf", height = 8, width = 9)
ctrl_whole_blod_cluster_phenotype_noh
dev.off()

pdf("./figures/ctrl_whole_blod_cluster_phenotype.pdf", height = 8, width = 9)
ctrl_whole_blod_cluster_phenotype
dev.off()



PBMC <- filterSCE(sce, k = "meta40", cluster_id %notin% c(1, 2, 6, 9, 11, 40))
PBMC <- CATALYST::cluster(PBMC, features = refined_markers, xdim = 10, ydim = 10, maxK = 40)

# weird <- filterSCE(PBMC, k = "meta40", cluster_id %in% c(39:40))
# not_weird <- filterSCE(PBMC, k = "meta40", cluster_id %notin% c(39:40))

PBMC <- filterSCE(PBMC, k = "meta40", cluster_id %notin% c(39: 40))
PBMC <- CATALYST::cluster(PBMC, features = refined_markers, xdim = 10, ydim = 10, maxK = 40)


vac63c_control_pbmc_cluster_heatmap <- plotExprHeatmap(x = PBMC, by = "cluster", row_clust = FALSE, col_clust = FALSE,  k = "meta30", bars = TRUE, features = c(all_markers, "GATA3", "RORgt", "140Ce"))

pdf("./figures/vac63c_control_pbmc_cluster_heatmap.pdf", height = 8, width = 9)
vac63c_control_pbmc_cluster_heatmap
dev.off()





pdf("./figures/PBMC_CXCR5_CD3_scatter.pdf")
plotScatter(PBMC, chs=c("CXCR5", "CD3"), color_by=NULL)
dev.off()

pdf("./figures/PBMC_CXCR5_HLADR_scatter.pdf")
plotScatter(PBMC, chs=c("CXCR5", "HLADR"), color_by=NULL)
dev.off()

pdf("./figures/PBMC_CXCR5_CD45RA_scatter.pdf")
plotScatter(PBMC, chs=c("CXCR5", "CD45RA"), color_by=NULL)
dev.off()

pdf("./figures/PBMC_CD3_CD45RA_scatter.pdf")
plotScatter(PBMC, chs=c("CD14", "CXCR5"), color_by=NULL)
dev.off()

meta30_table <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/debarcoded/supporting_files/meta30_PBMC_merge.csv")

PBMC <- CATALYST::mergeClusters(PBMC, k = "meta30", table = meta30_table, id = "flo_meta30")

system.time(PBMC <- scater::runUMAP(merged_PBMC,
                              subset_row=refined_markers,
                              exprs_values = "exprs",
                              scale=T))

# plotDR(x = PBMC, dr = "UMAP", color_by = "CD3")
# ggplot2::ggsave("./figures/PBMC_UMAP_CD3.png")
# plotDR(x = PBMC, dr = "UMAP", color_by = "CD4")
# ggplot2::ggsave("./figures/PBMC_UMAP_CD4.png")
# plotDR(x = PBMC, dr = "UMAP", color_by = "Va72")
# ggplot2::ggsave("./figures/PBMC_UMAP_Va72.png")
# plotDR(x = PBMC, dr = "UMAP", color_by = "TCRgd")
# ggplot2::ggsave("./figures/PBMC_UMAP_TCRgd.png")
# plotDR(x = PBMC, dr = "UMAP", color_by = "CXCR5")
# ggplot2::ggsave("./figures/PBMC_UMAP_CXCR5.png")
# plotDR(x = PBMC, dr = "UMAP", color_by = "CD14")
# ggplot2::ggsave("./figures/PBMC_UMAP_CD14.png")
# plotDR(x = PBMC, dr = "UMAP", color_by = "CD16")
# ggplot2::ggsave("./figures/PBMC_UMAP_CD16.png")
# plotDR(x = PBMC, dr = "UMAP", color_by = "CD56")
# ggplot2::ggsave("./figures/PBMC_UMAP_CD56.png")
#   
# plotDR(x = PBMC, dr = "UMAP", color_by = "flo_meta30")
# ggplot2::ggsave("./figures/PBMC_UMAP_flo_meta30.png", height = 10, width = 10)
#   
# plotDR(x = PBMC, dr = "UMAP", color_by = "ICOS")
# ggplot2::ggsave("./figures/PBMC_UMAP_ICOS.png")
# plotDR(x = PBMC, dr = "UMAP", color_by = "CD27")
# ggplot2::ggsave("./figures/PBMC_UMAP_CD27.png")
# lotDR(x = PBMC, dr = "UMAP", color_by = "PD1")
# ggplot2::ggsave("./figures/PBMC_UMAP_PD1.png")
# plotDR(x = PBMC, dr = "UMAP", color_by = "CD38")
# ggplot2::ggsave("./figures/PBMC_UMAP_CD38.png")

  

# \ooming in on those potensh RORgt+ cells, but i'm not convinced by this signal
# try <- filterSCE(PBMC, k = "meta45", cluster_id %in% c(22,2, 5))
# plotClusterExprs(try, "meta45", features = refined_markers)
# plotScatter(try, chs=c("RORgt", "CD16"), color_by=NULL, facet_by = "batch")

monocytes <- filterSCE(PBMC, k = "meta45", cluster_id %in% c(1:10, 15, 22, 39))
CATALYST::plotClusterExprs(monocytes, features=non_t_cell_markers, k="meta45")

pdf("scatter_try.pdf", height=8, width=8)
plotScatter(monocytes, chs = c("CD14", "HLADR"), color_by = NULL, label = "both", facet_by = "cluster_id")
dev.off()
