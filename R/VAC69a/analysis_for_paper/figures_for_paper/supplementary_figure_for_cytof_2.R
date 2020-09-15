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
                               show_heatmap_legend = FALSE,
                               name = "Median Marker Expression",
                               cluster_columns = FALSE,
                               row_names_side = "left",
                               col = col_inferno,
                               rect_gp = gpar(col = "white"),
                               #top_annotation = combo_top_anno,
                               column_names_rot = 45,
                               #heatmap_legend_param = list(col = col_fun4, title = "Normalised Frequency", title_position = "topleft"),
                               width = unit(16, "cm"),
                               height = unit(16*34/28, "cm")
)




png("/home/flobuntu/PhD/cytof/vac69a/final_figures_for_paper/supp_marker_expression_heatmap_all_clusters.png", width=11, height=9, units = "in", res=400)
draw(all_cluster_heatmap,
)
dev.off()

# Panel B: Series of UMAP projections showing expression of all 10 activation markers

#read in umap projection data
big_table <- data.table::fread("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_cells_with_UMAP_and_flo_merge_cluster.csv")

# restrict data to T6 and downsample by 33% to make it look nice
short_big_table_t6 <- subset(big_table, big_table$timepoint=="T6")
short_big_table_t6 <- short_big_table_t6[seq(1,nrow(short_big_table_t6), by=3), ]



CD27_plot <- flo_umap(short_big_table_t6, "CD27")+theme(axis.title = element_blank())
                                                                                                       

Ki67_plot <- flo_umap(short_big_table_t6, "Ki67")+theme(axis.title = element_blank())
                                                                                                       
CD28_plot <- flo_umap(short_big_table_t6, "CD28")+theme(axis.title = element_blank())
                                                                                                       
PD1_plot <- flo_umap(short_big_table_t6, "PD1")+theme(axis.title = element_blank())
                                                                                                     
CTLA4_plot <- flo_umap(short_big_table_t6, "CTLA4")+theme(axis.title = element_blank())
                                                                                                         

Tbet_plot <- flo_umap(short_big_table_t6, "Tbet")+theme(axis.title = element_blank())
                                                                                                       
Perforin_plot <- flo_umap(short_big_table_t6, "Perforin")+theme(axis.title = element_blank())
                                                                                                               

GZB_plot <- flo_umap(short_big_table_t6, "GZB")+theme(axis.title = element_blank())
                                                                                                     


HLADR_plot <- flo_umap(short_big_table_t6, "HLA-DR")+theme(axis.title = element_blank())
                                                                                                          
ICOS_plot <- flo_umap(short_big_table_t6, "ICOS")+theme(axis.title = element_blank())
                                                                                                       


UMAP_theme <- theme_minimal()+theme(
  panel.grid.minor = element_blank(),
  legend.position = "none",
  axis.text = element_blank()
)


{color_103_scheme <- c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                       "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                       "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                       "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                       "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                       "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                       "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                       "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                       "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                       "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                       "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                       "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                       "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")}

colcsv <- read.csv("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cluster_palette.csv", header=T, stringsAsFactors = F)
col_pal <- colcsv$x
names(col_pal) <- colcsv$X


expanded_cluster_palette <- c(col_pal, color_103_scheme[22:46])
names(expanded_cluster_palette)[11:length(expanded_cluster_palette)] <- unique(big_table$flo_label)[-match(unique(big_table$significant), unique(big_table$flo_label), nomatch = 0)]

# all cluster colours
t6_all_clusters_umap <- ggplot(short_big_table_t6, aes(x=UMAP1, y=UMAP2))+
  geom_point(aes(color=flo_label), shape=".", alpha=0.4)+
  scale_color_manual(values = expanded_cluster_palette)+
  theme_minimal()+
  ggtitle("Cluster ID")+
  UMAP_theme+
  guides(colour = guide_legend(override.aes = list(size = 2, shape=16, alpha=1), ncol = 3),
         alpha= "none")+
  theme(legend.position = "none",
        legend.text = element_text(size=6),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust=0.5))+
  coord_cartesian(xlim=c(-13, 10),
                  ylim=c(-11.2, 11.3))
                 


supp_fig_2_umaps <- plot_grid(t6_all_clusters_umap, Ki67_plot,  Tbet_plot, CTLA4_plot, PD1_plot,
                         GZB_plot, Perforin_plot, HLADR_plot, CD27_plot, ICOS_plot,  ncol=5)


ggsave("/home/flobuntu/PhD/cytof/vac69a/final_figures_for_paper/supp_activation_markers_whole_umap.png", supp_fig_2_umaps, height=4.8, width=9.6, units = "in")

