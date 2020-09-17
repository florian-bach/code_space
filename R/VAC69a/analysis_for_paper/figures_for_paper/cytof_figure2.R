#Panel A: Heatmap of Marker expression of significant clusters only ####



reordered_sig_scaled_mat <- as.matrix(read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/sig_cluster_medians_heatmap_t6.csv", header = T, row.names = 1))

sig_clusters <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/sig_t6_clusters.csv", header = TRUE, stringsAsFactors = FALSE, col.names = F)
sig_clusters <- sig_clusters[-9,]


# Right annotation ###

all_t6_data <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_t6_data.csv", stringsAsFactors = FALSE, header = T)

t6_map_data <- all_t6_data%>%
  group_by(cluster_id, timepoint) %>%
  mutate(mean_freq=mean(frequency)) %>%
  ungroup() %>%
  select(cluster_id, mean_freq, timepoint)

t6_map_data <- filter(t6_map_data, timepoint=="T6")
t6_map_data <- t6_map_data[!duplicated(t6_map_data), ]
t6_map_data <- filter(t6_map_data, cluster_id %in% sig_clusters)
t6_map_data <- t6_map_data[order(t6_map_data$mean_freq, decreasing = T),]

rereordered_sig_scaled_mat <- reordered_sig_scaled_mat[match(t6_map_data$cluster_id, rownames(reordered_sig_scaled_mat)),]



cd3_right_anno <- rowAnnotation(gap = unit(2, "mm"),
                                "Mean Frequency at T6" = anno_barplot(t6_map_data$mean_freq, which="row", axis = TRUE, ylim = c(0, 4)),
                                width = unit(4, "cm")
)



colcsv <- read.csv("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cluster_palette.csv", header=T, stringsAsFactors = F)
col_pal <- colcsv$x
names(col_pal) <- colcsv$X

colz <- unname(col_pal[match(rownames(rereordered_sig_scaled_mat), names(col_pal))])
breakz <- names(col_pal[match(rownames(rereordered_sig_scaled_mat), names(col_pal))])

cd3_right_anno_var <- rowAnnotation(gap = unit(2, "mm"),
                                    # "Mean Frequency at T6" = anno_boxplot(slim_wide_t6_map_data, which="row", axis = TRUE, gp=gpar(fill=c("#228833","#AA3377", "#AA3377", "#66CCEE", "#AA3377", "#AA3377", "#AA3377", "#4477AA", "#4477AA" ))),
                                    "Mean Percentage\nof CD3+ at T6" = anno_barplot(t6_map_data$mean_freq, which="row", axis = TRUE, ylim = c(0, 4), gp=gpar(fill=colz)),
                                    
                                    width = unit(4, "cm"),
                                    annotation_legend_param = list("Mean Percentage\nof CD3+ at T6" = list(title ="Lineage",
                                                                                                           at = breakz),
                                                                   legend_gp = gpar(fill = colz),
                                                                   title_position = "topleft")
)




draw(cd3_right_anno_var)


box_lgd <- Legend(labels =  breakz,
                  title = "Cluster_ID",
                  type = "grid",
                  legend_gp = gpar(fill = colz)
)


inferno <- colorspace::sequential_hcl("inferno", n=8)
col_inferno <- circlize::colorRamp2(seq(0,1, by=1/{length(inferno)-1}), inferno)



median_cluster_heat <- Heatmap(matrix = rereordered_sig_scaled_mat,
                               cluster_rows = FALSE,
                               name = "Median Marker Expression",
                               cluster_columns = FALSE,
                               row_names_side = "left",
                               col = col_inferno,
                               rect_gp = gpar(col = "white"),
                               #top_annotation = combo_top_anno,
                               right_annotation = cd3_right_anno_var,
                               show_heatmap_legend = FALSE,
                               column_names_rot = 60,
                               #heatmap_legend_param = list(col = col_fun4, title = "Normalised Frequency", title_position = "topleft"),
                               width = unit(16, "cm"),
                               height = unit(16*9/28, "cm")
)


png("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/diffcyt/edgeR/sig_cluster_t6_phenotype_heat_var.png", width=13, height=4, units = "in", res=400)
draw(median_cluster_heat,
     #annotation_legend_list = list(box_lgd),
     merge_legends = FALSE,
     #padding = unit(c(200, 200, 200, 200), "mm")
)
dev.off() 


#Panel B: Pie Chart of significant cluster size and lineage ####
library(gridBase)
library(grid)
library(circlize)
library(ComplexHeatmap)
library(dplyr)



all_t6_data <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_t6_data.csv", stringsAsFactors = FALSE, header = T)

# calculate average frequencies of each cluster at each timepoint, select only T6 data

t6_map_data <- all_t6_data%>%
  group_by(cluster_id, timepoint) %>%
  mutate(mean_freq=mean(frequency)) %>%
  ungroup() %>%
  select(cluster_id, mean_freq, timepoint) %>%
  filter(timepoint=="T6")

t6_map_data <- t6_map_data[!duplicated(t6_map_data), ]

#restrict data to significant clusters
sig_clusters <- scan("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/sig_t6_clusters.csv", what ="character", sep = ",", skip = 1)
sig_clusters <- sig_clusters[-9]
t6_map_data <- filter(t6_map_data, cluster_id %in% sig_clusters)
t6_map_data <- t6_map_data[order(t6_map_data$mean_freq, decreasing = T),]

# PIE CONSTRUCTION USING CIRCLIZE
pie_data <- data.frame(t6_map_data)


# manually fix lineage- this could be fixed by reading in the lineage containign file on line 108 but who cares
pie_data$lineage <- factor(c("gd", "CD4", "CD4", "MAIT", "CD4", "CD4", "CD4", "CD8", "CD8"), levels=c("CD4", "CD8", "MAIT", "gd"))
pie_data <- pie_data[order(pie_data$lineage),]


# manula(ish) mapping of colours to cluster names & lineages
colcsv <- read.csv("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cluster_palette.csv", header=T, stringsAsFactors = F)
col_pal <- colcsv$x
names(col_pal) <- colcsv$X
ordered_col_pal <- col_pal[match(pie_data$cluster_id, names(col_pal))]

#actual pie construction
circlize_plot = function() {
  
  # set the thickness of the pie and the starting point as 12 o'clock; set divisions of pie according to cluster name, thickness
  # of slice is set to cluster frequency in percent
  circos.par("track.height" = 0.6, start.degree = 90)
  circos.initialize(factors = pie_data$cluster_id, xlim=c(0,1), sector.width = pie_data$mean_freq)
  
  o.cell.padding = circos.par("cell.padding")
  
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.2, bg.border = NA, 
                         cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4]))
  
  # actually initialise the coloured slizes
  circos.track(factors = pie_data$cluster_id,
               ylim = c(0, 1),
               x=pie_data$mean_freq,
               bg.col = ordered_col_pal)
  
  
  #outer circle labelling the clusters with lineage: no color, just text
  
  highlight.sector(sector.index = grep("CD4", pie_data$cluster_id, value = T), track.index = 1, 
                   border = "black", text = "CD4", col = "white")
  highlight.sector(sector.index = grep("CD8", pie_data$cluster_id, value = T), track.index = 1, 
                   border = "black", text = "CD8", col = "white")
  highlight.sector(sector.index = grep("MAIT", pie_data$cluster_id, value = T), track.index = 1, 
                   border = "black", text = "MAIT", col = "white")
  highlight.sector(sector.index = grep("gamma delta", pie_data$cluster_id, value = T), track.index = 1, 
                   border = "black", text = "Vd2", col = "white", )
  
  
  # # Just color, no text
  # highlight.sector(sector.index = grep("CD4", pie_data$cluster_id, value = T), track.index = 1, 
  #                  col = pie_palette[1], border = NA)
  # highlight.sector(sector.index = grep("CD8", pie_data$cluster_id, value = T), track.index = 1, 
  #                  col = pie_palette[2], border = NA)
  # highlight.sector(sector.index = grep("MAIT", pie_data$cluster_id, value = T), track.index = 1, 
  #                  col = pie_palette[3], border = NA)
  # highlight.sector(sector.index = grep("Vd2", pie_data$cluster_id, value = T), track.index = 1, 
  #                  col = pie_palette[4], border = NA)
  
  circos.clear()
  
}

circlize_plot()  


# add legend for pie chart
lvl_data <- pie_data
lvl_data$lineage <- factor(lvl_data$lineage, levels = c("CD4", "CD8", "MAIT", "gd"))
lvls <- lvl_data$cluster_id[c(seq(1,7), 9, 8 )]

lgd_cluster= Legend(at =lvls, type = "grid", 
                    legend_gp = gpar(fill = unname(ordered_col_pal)[c(seq(1,7), 9, 8 )]),
                    title_position = "topleft",
                    labels_gp = gpar(fontsize = 10),
                    title = "Cluster ID")



#write out pie chart with legend

png("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cluster_activation_pie.png", width=7, height=4, units = "in", res=400)

plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circlize_plot()
par(xpd = NA)
upViewport()
draw(lgd_cluster, x = circle_size, just = "left")
#title("Average Size of Significantly Upregulated Clusters at T6", 0,6)

dev.off()



#Panel C: stacked Piechart of activated CD4 cluster sizes relative to size of memroy pool ####

# we calculate significant cluster size as a fraction of the CD4 memory pool, to do this we first assemble a data frame containing
# all non-Treg CD4 clusters, then remove the naive one. For each volunteer, the remaining clusters are then summed and each cluster
# is then divided by this sum


all_t6_data <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_t6_data.csv", stringsAsFactors = FALSE, header = T)
all_t6_data <- subset(all_t6_data, all_t6_data$timepoint == "T6")
cd4_t6_data <- subset(all_t6_data, grepl("CD4 ", all_t6_data$cluster_id))
cd4_memory_t6_data <- subset(cd4_t6_data, !grepl("Naïve", cd4_t6_data$cluster_id))


cd4_memory_t6_summary <- cd4_memory_t6_data %>%
  group_by(volunteer) %>%
  mutate("CD4_memory_percentage"=sum(frequency)) %>%
  select(volunteer, CD4_memory_percentage)

cd4_memory_t6_summary <- cd4_memory_t6_summary[!duplicated(cd4_memory_t6_summary), ]

cd4_bar_data <- cd4_summary %>%
  group_by(volunteer) %>%
  mutate("cd4_freq" = frequency/cd4_memory_t6_summary$CD4_memory_percentage[match(volunteer, cd4_memory_t6_summary$volunteer)])

cd4_bar_data <- filter(cd4_bar_data, cluster_id %in% sig_clusters)


stacked_bar_levels <- c("activated CD27- cytotoxic CD4 EM",
                        "activated CD4 EM",
                        "activated CTLA4+ CD4 EM",
                        "activated CD4 CM")

cd4_memory_activation_stacked_barchart <- ggplot(cd4_bar_data, aes(x=volunteer, y=cd4_freq, fill=factor(cluster_id, levels = stacked_bar_levels)))+
  geom_bar(stat="identity", position="stack")+
  ggtitle("Activation of CD4 Memory Cells at T6\n")+
  scale_fill_manual(values=col_pal)+
  scale_y_continuous(name = "Percentage of Activated\nCD4 memory T cells", labels=scales::percent_format(accuracy = 1))+
  theme_minimal()+
  xlab("Volunteer")+
  guides(fill=guide_legend(ncol=2))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        legend.title = element_blank(),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        legend.position = "bottom",
        legend.justification = "center",
        panel.grid.minor.y = element_blank(),
        strip.placement = "outside")


ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cd4_memory_activation_stacked_barchart.png", cd4_memory_activation_stacked_barchart, width=6, height = 6)

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cd4_memory_activation_stacked_barchart_no_leg.png", cd4_memory_activation_stacked_barchart, width=5, height = 5)



#Panel D: Show marker expression on T6 CD4 bump ####

library(vac69a.cytof)
library(ggplot2)
library(cowplot)

sig_clusters <- scan("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/sig_t6_clusters.csv", what ="character", sep = ",", skip = 1)
sig_clusters <- sig_clusters[-9]


big_table <- data.table::fread("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_cells_with_UMAP_and_flo_merge_cluster.csv")
# colnames(big_table) <- gsub("HLADR", "HLA-DR", colnames(big_table))
# big_table <- data.table::fwrite(x = big_table, "~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_cells_with_UMAP_and_flo_merge_cluster.csv")

big_table$significant <- ifelse(big_table$flo_label %in% sig_clusters, big_table$flo_label, "black")
big_table$alpha <- ifelse(big_table$flo_label %in% sig_clusters, 1, 0.5)



sig_cd4_clusters <- grep("CD4", sig_clusters, value=T)
shorter_big_table_t6 <- subset(big_table, big_table$timepoint=="T6")
shorter_big_table_t6 <- subset(shorter_big_table_t6, shorter_big_table_t6$flo_label %in% sig_clusters)


CD27_plot <- flo_umap(shorter_big_table_t6, "CD27")+theme(axis.title = element_blank())+coord_cartesian(xlim = c(0,4),
                                                                                                        ylim = c(2,4))

Ki67_plot <- flo_umap(shorter_big_table_t6, "Ki67")+theme(axis.title = element_blank())+coord_cartesian(xlim = c(0,4),
                                                                                                        ylim = c(2,4))
CD28_plot <- flo_umap(shorter_big_table_t6, "CD28")+theme(axis.title = element_blank())+coord_cartesian(xlim = c(0,4),
                                                                                                        ylim = c(2,4))
PD1_plot <- flo_umap(shorter_big_table_t6, "PD1")+theme(axis.title = element_blank())+coord_cartesian(xlim = c(0,4),
                                                                                                      ylim = c(2,4))
CTLA4_plot <- flo_umap(shorter_big_table_t6, "CTLA4")+theme(axis.title = element_blank())+coord_cartesian(xlim = c(0,4),
                                                                                                          ylim = c(2,4))

Tbet_plot <- flo_umap(shorter_big_table_t6, "Tbet")+theme(axis.title = element_blank())+coord_cartesian(xlim = c(0,4),
                                                                                                        ylim = c(2,4))
Perforin_plot <- flo_umap(shorter_big_table_t6, "Perforin")+theme(axis.title = element_blank())+coord_cartesian(xlim = c(0,4),
                                                                                                                ylim = c(2,4))

GZB_plot <- flo_umap(shorter_big_table_t6, "GZB")+theme(axis.title = element_blank())+coord_cartesian(xlim = c(0,4),
                                                                                                      ylim = c(2,4))


HLADR_plot <- flo_umap(shorter_big_table_t6, "HLA-DR")+theme(axis.title = element_blank())+coord_cartesian(xlim = c(0,4),
                                                                                                          ylim = c(2,4))
ICOS_plot <- flo_umap(shorter_big_table_t6, "ICOS")+theme(axis.title = element_blank())+coord_cartesian(xlim = c(0,4),
                                                                                                        ylim = c(2,4))


UMAP_theme <- theme_minimal()+theme(
  panel.grid.minor = element_blank(),
  legend.position = "none",
  axis.text = element_blank()
)



fig2_sig_clusters_umap <- ggplot(shorter_big_table_t6, aes(x=UMAP1, y=UMAP2))+
  geom_point(aes(color=significant), shape=".")+
  scale_color_manual(values = col_pal)+
  theme_minimal()+
  ggtitle("Cluster ID")+
  UMAP_theme+
  guides(colour = guide_legend(override.aes = list(size = 4)),
         alpha= "none")+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust=0.5))+
  coord_cartesian(xlim = c(0,4),
                  ylim = c(2,4))


short_big_table_t6 <- subset(big_table, big_table$timepoint=="T6")
short_big_table_t6 <- short_big_table_t6[seq(1,nrow(short_big_table_t6), by=3), ]


zoom_plot <-ggplot(short_big_table_t6, aes(x=UMAP1, y=UMAP2))+
  geom_point(aes(color=significant,  alpha=alpha), shape=".")+
  theme_minimal()+
  scale_color_manual(values=col_pal)+
  geom_rect(xmin=0, xmax=4, ymin=2, ymax=4, color="red", fill=NA, size=0.5)+
  UMAP_theme+
  guides(colour = guide_legend(override.aes = list(size = 4)),
         alpha= "none")+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(2.5,0.2,2.5,0.2), "cm"),
        plot.title = element_text(hjust=0.5))+
  coord_cartesian(xlim=c(-13, 10),
                  ylim=c(-11.2, 11.3))




fig_2_umaps <- plot_grid(fig2_sig_clusters_umap, Ki67_plot,  Tbet_plot, CTLA4_plot, PD1_plot,
                         GZB_plot, Perforin_plot, HLADR_plot, CD27_plot, ICOS_plot,  ncol=5)

fig_2_umaps_var <- plot_grid(zoom_plot, fig_2_umaps, rel_widths = c(1,5))
#fig_2_umaps_var <- plot_grid(zoom_plot, fig_2_umaps, rel_widths = c(1,4x))

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/fig_2_umaps_var3.png", fig_2_umaps_var, height=4, width=8, units = "in")

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/fig_2_umaps.png", fig_2_umaps, height=4, width=8, units = "in")



