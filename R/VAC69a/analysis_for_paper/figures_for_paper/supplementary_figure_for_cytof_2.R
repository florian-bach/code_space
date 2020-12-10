# Panel A: Heatmap of all clusters' median expression of all markers ####
library(ggplot2)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)

inferno <- colorspace::sequential_hcl("inferno", n=8)
col_inferno <- circlize::colorRamp2(seq(0,1, by=1/{length(inferno)-1}), inferno)



ms_mean <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/cluster_medians_all.csv", header=T, stringsAsFactors = F)
ms_mean <- select(ms_mean, cluster_id, Marker, mean_expression)


all_mat <- as.matrix(tidyr::spread(ms_mean, cluster_id, mean_expression))

rownames(all_mat) <- all_mat[,1]
all_mat2 <- all_mat[,2:ncol(all_mat)]

scaled_all_mat <- apply(apply(all_mat2, c(1,2), as.numeric), MARGIN = 1, function(x)scales::rescale(x, to=c(0, 1)))

clustering_markers <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/refined_markers.csv", stringsAsFactors = F)


reordered_scaled_mat <- scaled_all_mat[,match(clustering_markers$refined_markers, colnames(scaled_all_mat))]



all_cluster_heatmap <- Heatmap(matrix = reordered_scaled_mat,
                               cluster_rows = T,
                               show_row_dend = FALSE,
                               show_heatmap_legend = TRUE,
                               name = "Median Marker Expression",
                               cluster_columns = FALSE,
                               column_names_gp = gpar(fontsize = 8),
                               row_names_gp = gpar(fontsize = 8),
                               row_names_side = "left",
                               col = col_inferno,
                               rect_gp = gpar(col = "white"),
                               #top_annotation = combo_top_anno,
                               column_names_rot = 45,
                               heatmap_legend_param = list(col = col_inferno,
                                                           #legend_position = "bottom",
                                                           at=c(0,0.5,1),
                                                           title = "Normalised Marker Expression",
                                                           legend_direction = "vertical",
                                                           title_position = "leftcenter-rot",
                                                           legend_height = unit(5.2, "cm"),
                                                           legend_width = unit(0.8, "cm"),
                                                           border = FALSE)
                               #heatmap_legend_param = list(col = col_fun4, title = "Normalised Frequency", title_position = "topleft"),
                               # width = unit(16, "cm"),
                               # height = unit(16*34/28, "cm"
                                             )




pdf("/home/flobuntu/PhD/cytof/vac69a/final_figures_for_paper/supp_marker_expression_heatmap_all_clusters.pdf", width=8, height=7)
draw(all_cluster_heatmap,
     heatmap_legend_side = "right"
)
dev.off()

# Panel B: Series of UMAP projections showing expression of all 10 activation markers

#read in umap projection data
big_table <- data.table::fread("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_cells_with_UMAP_and_flo_merge_cluster.csv")

# restrict data to T6 and downsample by 33% to make it look nice
short_big_table_t6 <- subset(big_table, big_table$timepoint=="T6")
short_big_table_t6 <- short_big_table_t6[seq(1,nrow(short_big_table_t6), by=3), ]

supp_theme <- theme(axis.title = element_text(size = 6),
                    legend.title = element_text(size = 6),
                    legend.text = element_text(size=6),
                    plot.title =  element_text(size=7))

CD27_plot <- flo_umap(short_big_table_t6, "CD27")+supp_theme
                                                                                                       

Ki67_plot <- flo_umap(short_big_table_t6, "Ki67")+supp_theme
                                                                                                       
CD28_plot <- flo_umap(short_big_table_t6, "CD28")+supp_theme
                                                                                                       
PD1_plot <- flo_umap(short_big_table_t6, "PD1")+supp_theme
                                                                                                     
CTLA4_plot <- flo_umap(short_big_table_t6, "CTLA4")+supp_theme
                                                                                                         

Tbet_plot <- flo_umap(short_big_table_t6, "Tbet")+supp_theme
                                                                                                       
Perforin_plot <- flo_umap(short_big_table_t6, "Perforin")+supp_theme
                                                                                                               

GZB_plot <- flo_umap(short_big_table_t6, "GZB")+supp_theme
                                                                                                     


HLADR_plot <- flo_umap(short_big_table_t6, "HLA-DR")+supp_theme
                                                                                                          
ICOS_plot <- flo_umap(short_big_table_t6, "ICOS")+supp_theme
                                                                                                       



supp_fig_2_umaps <- plot_grid(Ki67_plot,  Tbet_plot, CTLA4_plot, PD1_plot,
                         GZB_plot, Perforin_plot, HLADR_plot, CD27_plot, ICOS_plot,  ncol=3, align = "h", axis="trbl")


ggsave("/home/flobuntu/PhD/cytof/vac69a/final_figures_for_paper/supp_activation_markers_whole_umap.pdf", supp_fig_2_umaps, height=3.5, width=5.25)

