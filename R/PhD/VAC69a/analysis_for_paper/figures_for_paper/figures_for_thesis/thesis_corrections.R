library(vac69a.cytof)
library(ggplot2)

setwd("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")

smol_daf <- vac69a.cytof::read_full("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")
# #smol_daf <- read_small("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/", proportional=F, event_number=3000)
# #
# plotClusterHeatmap(smol_daf, hm2 = NULL,
#                    k = "meta45",
#                    # m = "flo_merge",
#                    cluster_anno = TRUE,
#                    draw_freqs = TRUE,
#                    scale=T
# )
#

#
#smol_daf <- filterSCE(smol_daf, timepoint!="C10")

refined_markers <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/refined_markers.csv", stringsAsFactors = F)

all_markers <- c(refined_markers$refined_markers, "CXCR5", "CD56", "CD49d")


reclust <- CATALYST::cluster(smol_daf, features = all_markers, maxK = 50, xdim = 10, ydim = 10)

inferno <- colorspace::sequential_hcl("inferno", n=8)

plotty <- CATALYST::plotExprHeatmap(reclust,
                                    # features = all_markers[-match(c("CD56", "CD49d"), all_markers)],
                                    features = all_markers,
                                    by = "cluster",
                                    row_clust = FALSE,
                                    col_clust = FALSE,
                                    hm_pal=inferno,
                                    k= "som100",
                                    assay="exprs")

pdf(file = "~/PhD/figures_for_thesis/corrections/test.pdf", height=11, width=8)
plotty
dev.off()

# umap projections ####
set.seed(1234);smol_daf <- scater::runUMAP(smol_daf,
                                             subset_row=all_markers,
                                             exprs_values = "exprs",
                                             scale=T)

#smol_daf <- CATALYST::filterSCE(smol_daf, timepoint!="C10")

big_table <- prep_sce_for_ggplot(smol_daf)


down_big_table <- big_table[seq(1,nrow(big_table), by=3), ]


supp_theme <- theme(axis.title = element_text(size = 6),
                    legend.title = element_text(size = 6),
                    legend.text = element_text(size=6))




CD4_plot<- flo_umap(down_big_table, "CD4")+supp_theme+ggtitle("CD4")
CD8_plot<- flo_umap(down_big_table, "CD8")+supp_theme+ggtitle("CD8")
Vd2_plot <- flo_umap(down_big_table, "Vd2")+supp_theme+ggtitle(expression(paste(gamma,delta)))
Va72_plot <- flo_umap(down_big_table, "Va72")+supp_theme+ggtitle(expression(paste("V",alpha,"7.2")))
CCR7_plot <- flo_umap(down_big_table, "CCR7")+supp_theme+ggtitle("CCR7")
CD45RO_plot <- flo_umap(down_big_table, "CD45RO")+supp_theme+ggtitle("CD45RO")

CD49d_plot<- flo_umap(down_big_table, "CD49d")+supp_theme+ggtitle("CD49d")
CXCR5_plot<- flo_umap(down_big_table, "CXCR5")+supp_theme+ggtitle("CXCR5")
CD56_plot<- flo_umap(down_big_table, "CD56")+supp_theme+ggtitle("CD56")


CXCR5_plot<- flo_umap(down_big_table, "CXCR5")+supp_theme+ggtitle("CXCR5")
ICOS_plot<- flo_umap(down_big_table, "ICOS")+supp_theme+ggtitle("ICOS")
PD1_plot<- flo_umap(down_big_table, "PD1")+supp_theme+ggtitle("PD1")


ash_plot_big <- cowplot::plot_grid(CD4_plot, CD8_plot, Vd2_plot, Va72_plot, CCR7_plot, CD45RO_plot, CD49d_plot, CD56_plot, CXCR5_plot, ncol=3)
ggsave("~/PhD/figures_for_thesis/corrections/ash_umap_big.png", ash_plot_big, height = 9, width=8)
