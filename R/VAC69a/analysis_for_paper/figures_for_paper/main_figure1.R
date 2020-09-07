# Panel A: UMAP Density Plot ####

library(ggplot2)

big_table <- data.table::fread("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_cells_with_UMAP_and_flo_merge_cluster.csv")
sig_clusters <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/sig_t6_clusters.csv", header = TRUE, stringsAsFactors = FALSE)

#get rid of resting Vd cluster
sig_clusters <- sig_clusters[-9,2]

down_big_table <- big_table[seq(1,nrow(big_table), by=3), ]


# DOPE CONTOUR PLOT #### n=13 is your friend :*
hex_through_time <- ggplot(down_big_table, aes(x=UMAP1, y=UMAP2))+
  stat_density_2d(aes(fill = ..density..), geom = 'raster', contour = FALSE, n = 1500)+
  stat_density_2d(contour = TRUE, bins=13, color="white", size=0.05)+
  scale_x_continuous(limits=c(-13.5, 10))+
  scale_y_continuous(limits=c(-11.2, 11.3))+
  theme_minimal()+
  facet_wrap(~timepoint, ncol = 4)+
  UMAP_theme+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size=11),
        axis.title = element_blank(),
        axis.text = element_blank())+
  scale_fill_gradientn(colors=inferno_white)+
  coord_cartesian(xlim=c(-13.5, 10),
                  ylim=c(-11.2, 11.3))
#viridis::scale_fill_viridis(option="B"))


system.time(ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/proportional_raster_umap13_005_no_axis_title.png", hex_through_time, height=4, width=16))


ggsave(paste("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/proportional_contour_umap_14", "_bins",".png", sep=''), hex_through_time, height=6, width=9)

panelA_B  <- plot_grid(hex_through_time, gate_labeled_gg, ncol=2, rel_widths=c(3,1.5), labels = c("A", "B"))
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/panelA_B.png", panelA_B, height=6, width=11.25)

# Panel B: Black UMAP projection with red, labeled gates####
library(ggcyto)
library(CytoML)
library(ggplot2)

fcs_name <- "/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/post_umap/concat/all_t6.fcs"

lin_gates <- cytobank_to_gatingset("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/CytExp_295996_Gates_simplified.xml", FCS=fcs_name)

gate_name <- substr(gs_get_pop_paths(lin_gates)[2:length(gs_get_pop_paths(lin_gates))], 2, 30)
gate_name <- gsub(" ", "\n", gate_name)
gate_name[7] <- "paste(gamma, delta)"

x_coord <- c(-1.1, -2.5,   5, -10,   -9,    7,    -4.8, -10.6)
y_coord <- c(5,    -2.4, 7.5,   7,    -8,  -5.8,   8.3, -4.5)

gate_label_positions <- data.frame(gate_name, x_coord, y_coord)

umap_gate_labeled <- ggcyto(lin_gates, aes(x=UMAP1, y=UMAP2))+
    #geom_hex(bins=228)+
    geom_point(color="black", alpha=0.5, shape=".")+
    geom_gate(c(gs_get_pop_paths(lin_gates)[2:length(gs_get_pop_paths(lin_gates))]))+
    geom_text(data=gate_label_positions[7,], aes(x=x_coord, y=y_coord, label=gate_name, size=20), lineheight = 0.7, parse = T)+
    geom_text(data=gate_label_positions[-7,], aes(x=x_coord, y=y_coord, label=gate_name, size=20), lineheight = 0.7)+
    UMAP_theme+
    ggtitle("Major T cell\nlineages")+
    theme(axis.title = element_blank(), 
          strip.text = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(hjust=0.5))


gate_labeled_gg <- as.ggplot(umap_gate_labeled)

gate_labeled_gg <- gate_labeled_gg+
  coord_cartesian(xlim=c(-13, 10), ylim=c(-11.2, 11.3))

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/umap_gate_labeled.png", umap_gate_labeled, height=4, width=4)



# Panel C: stacked barchart of T cell activation, colored by lineage, split by volunteer; with lineage pies ####

#barchart

library(vac69a.cytof)

summary <- data.table::fread("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/activation_barchart_data")


lineage_palette <- c("#AA3377", "#EE6677", "#4477AA", "#66CCEE", "#228833", "#FFA500", "#BBBBBB")
names(lineage_palette) <-c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "Resting")



activation_stacked_barchart <- ggplot(summary, aes(x=volunteer, y=frequency/100, fill=lineage))+
  geom_bar(stat="identity", position="stack")+
  #geom_text(aes(y=(total_cd3/100)+0.01, label=paste0(round(total_cd3, digits = 1), "%", sep='')))+
  theme_minimal()+
  facet_wrap(~timepoint, strip.position = "bottom", ncol=4)+
  ggtitle("Overall T cell activation")+
  scale_fill_manual(values=lineage_palette, labels=c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "Resting"))+
  scale_y_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1))+
  #ylim(0,25)+
  #geom_text(aes(label=cluster_id), position = position_stack(vjust = .5))+
  theme(#legend.position = "none",
    plot.title = element_text(hjust=0.5, size=11),
    strip.text = element_text(hjust=0.5, size=10, face = "bold"),
    axis.title.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none",
    strip.placement = "outside")


ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/activation_stacked_barchart.png", activation_stacked_barchart, height=4, width=8)

# lineage pies #

#get cluster freqs and cluster lineage freqs
barchart_data <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/cluster_and_lineage_stats_for_barcharts.csv", header=T, stringsAsFactors = FALSE)

#keep only activated clusters at T6
activated_clusters <- grep("activated", unique(barchart_data$cluster_id), value = T)
barchart_data <- subset(barchart_data, barchart_data$cluster_id %in% activated_clusters)
barchart_data <- subset(barchart_data, barchart_data$timepoint=="T6")

# get average lineage activation
barchart_data <- barchart_data %>%
  group_by(timepoint, lineage) %>%
  mutate("lin_activation" = sum(lin_freq)/6) %>%
  ungroup()


# construct ggplotable dataframe; now this table contains each average for all volunteers, so get rid of duplicates, then append
# the df with 1 -the same values so that the activated and nonactivated percentages are in the same column

barchart_plot_data <- data.frame("lineage" = barchart_data$lineage, "lin_activation" = barchart_data$lin_activation)
barchart_plot_data <- barchart_plot_data[!duplicated(barchart_plot_data), ]
barchart_plot_data <- rbind(barchart_plot_data, data.frame("lineage"=barchart_plot_data$lineage, "lin_activation"=1-barchart_plot_data$lin_activation))

# add activated/resting column, make it so that the greek letters render in ggplot
barchart_plot_data$state <- c(as.character(barchart_plot_data$lineage[1:6]), rep("Resting", 6))
barchart_plot_data$state <- gsub("gd", expression(paste(gamma, delta)), barchart_plot_data$state)
barchart_plot_data$lineage <- gsub("gd", expression(paste(gamma, delta)), barchart_plot_data$lineage)


# we have to rename the colour palette because otherwise the strip.text won't parse properly
lineage_palette <- c("#AA3377", "#EE6677", "#4477AA", "#66CCEE", "#228833", "#FFA500", "#BBBBBB")
names(lineage_palette) <-c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "Resting")

pie_lineage_palette <- lineage_palette
names(pie_lineage_palette)<- c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "Resting")

#order lineage and state columns for plotting
barchart_plot_data$lineage <- factor(barchart_plot_data$lineage, levels=c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "Resting"))
barchart_plot_data$state <- factor(barchart_plot_data$state, levels=rev(c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "Resting")))

lineage_activation_pies <- ggplot(barchart_plot_data)+
  geom_bar(stat="identity", position = "stack", width=1, aes(x="", y=lin_activation, fill=state))+
  theme_minimal()+
  coord_polar(theta = "y", start = 0)+
  facet_wrap(~lineage, labeller = label_parsed, ncol=2, strip.position = "left")+
  scale_fill_manual(values=pie_lineage_palette)+
  guides()+
  ggtitle("Average Proportion of Activated Cells\nin Major T cell Lineages at T6")+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size=11),
        strip.text.y.left = element_text(hjust=0.5, size=14, face = "bold", angle = 0),
        legend.position = "none")

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/lineage_activation_pies.png", lineage_activation_pies, height=3, width=4.5)


# put the bar and pie charts together with cowlplot
panel_f <- cowplot::plot_grid(activation_stacked_barchart, lineage_activation_pies, ncol = 2, rel_widths = c(4,1), rel_heights = c(1,1), axis = "bt", align = "hv")

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/panel_f.png", panel_f, height=4, width=17)




# Panel D: UMAP projections coloured by CD38, Bcl2, all clusters, sig clusters ####

#read data
big_table <- data.table::fread("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_cells_with_UMAP_and_flo_merge_cluster.csv")

# ATTENTION: for some reason in the writing of this there were some spaces at the end of some names so if the cluster coloring
# doesn't work check that the "significant" column contains all the significant clusters!
sig_clusters <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/sig_t6_clusters.csv", header = TRUE, stringsAsFactors = FALSE)

#get rid of resting Vd cluster
sig_clusters <- sig_clusters[-9,]

# create column that colours cell either black or cluster-specifically, adjust alpha so that colored cells are completely opaque
# while black cells are partly transparent
big_table$significant <- ifelse(big_table$flo_label %in% sig_clusters, big_table$flo_label, "black")
big_table$alpha <- ifelse(big_table$flo_label %in% sig_clusters, 1, 0.5)

# restrict data to T6 and downsample by 33% to make it look nice
short_big_table_t6 <- subset(big_table, big_table$timepoint=="T6")
short_big_table_t6 <- short_big_table_t6[seq(1,nrow(short_big_table_t6), by=3), ]

##define color schemes
colcsv <- read.csv("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cluster_palette.csv", header=T, stringsAsFactors = F)
col_pal <- colcsv$x
names(col_pal) <- colcsv$X


#color_103_scheme
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

expanded_cluster_palette <- c(col_pal, color_103_scheme[22:46])
names(expanded_cluster_palette)[11:length(expanded_cluster_palette)] <- unique(big_table$flo_label)[-match(unique(big_table$significant), unique(big_table$flo_label), nomatch = 0)]

# all cluster colours
t6_all_clusters_umap <- ggplot(short_big_table_t6, aes(x=UMAP1, y=UMAP2))+
  geom_point(aes(color=flo_label), shape=".", alpha=0.4)+
  scale_color_manual(values = expanded_cluster_palette)+
  theme_minimal()+
  ggtitle("All Cells at T6 coloured\nby Cluster ID")+
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

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/t6_all_clusters_umap_without_legend.png",  t6_all_clusters_umap, height=4, width=4)



#sig cluster colours only

t6_sig_clusters_umap <- ggplot(short_big_table_t6, aes(x=UMAP1, y=UMAP2))+
  geom_point(aes(color=significant, alpha=alpha), shape=".")+
  scale_color_manual(values = col_pal)+
  theme_minimal()+
  ggtitle("Differentially Abundant\nClusters at T6")+
  UMAP_theme+
  guides(colour = guide_legend(override.aes = list(size = 4)),
         alpha= "none")+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust=0.5))+
  coord_cartesian(xlim=c(-13, 10),
                  ylim=c(-11.2, 11.3))

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/t6_sig_clusters_umap_var.png",  t6_sig_clusters_umap, height=5, width=5)


# UMAP coloured by CD38, Bcl2
cd38_plot <- flo_umap(short_big_table_t6, "CD38")+theme(axis.title = element_blank())
bcl2_plot <- flo_umap(short_big_table_t6, "BCL2")+theme(axis.title = element_blank())


horizontal_d_f_panel <- cowplot::plot_grid(cd38_plot, bcl2_plot, t6_all_clusters_umap, t6_sig_clusters_umap, ncol=4)
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/horizontal_d_f_panel.png",  horizontal_d_f_panel, height=4, width=16)


# Panel E: Differential Abundance Heatmap ####
library(colorspace)
library(scales)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)

time_col <- colorspace::sequential_hcl(5, palette = "Purple Yellow")


#cluster id and frequencies
all_t6_data <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_t6_data.csv")
#meta data
cd <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/experiment_info.csv")
#fold change and p_adj
fold_change <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/edger_t6_FC_padj.csv")

#order decreasing so top10 and bottom 10 are extreme in terms of fold change
fold_change <- fold_change[order(fold_change$logFC, decreasing = T),]

#get the correct order of cluster ids by fold change
cluster_levels <- fold_change$cluster_id


#calculate z score after asinh transform of square
t6_map_data <- all_t6_data%>%
  group_by(cluster_id) %>%
  #group_by(volunteer) %>%
  mutate(trans_freq=asin(sqrt(frequency/100))) %>%
  mutate(max_freq=max(frequency)) %>%
  mutate(trans_norm_freq=scale(trans_freq, center = TRUE, scale = TRUE)) %>%
  ungroup()



# construct matrix to be plotted in heatmap: clusters along rows, samples along columns, values=trans_norm_freq
hm_matrix <- as.data.frame(t6_map_data %>%
                             select(cluster_id, sample_id, trans_norm_freq) %>%
                             tidyr::pivot_wider(names_from = sample_id, values_from = trans_norm_freq))

rownames(hm_matrix) <- hm_matrix$cluster_id
# this orders the matrix according to fold change, as cluster_levels is..
hm_matrix <- hm_matrix[match(cluster_levels, hm_matrix$cluster_id),]

#matrix has to be numeric only
hm_matrix <- as.matrix(select(hm_matrix, -cluster_id))


#truncate matrix values so that extreme values don't pale everything elese
hm_matrix <- ifelse(hm_matrix>2.5, 2.5, hm_matrix)
hm_matrix <- ifelse(hm_matrix < -2.5, -2.5, hm_matrix)
extreme <- max(abs(range(hm_matrix)))

#bunch of differen color gradients
col_fun4 <- circlize::colorRamp2(c(-range(hm_matrix)[2], 0, range(hm_matrix)[2]), c("#0859C6", "black", "#FFA500"))

#get 12 most up/down clusters
top_mat <- head(hm_matrix, n=12)
bot_mat <- tail(hm_matrix, n=12)
combo_matrix <- rbind(top_mat, bot_mat)

# make seperate array for p_adj, significantce, fold change etc for heatmap annotations
# right annotation
p_adj <- c(head(fold_change$p_adj, n=nrow(top_mat)), tail(fold_change$p_adj, n=nrow(bot_mat)))
log2_fc <- c(head(fold_change$logFC, n=nrow(top_mat)), tail(fold_change$logFC, n=nrow(bot_mat)))
significant <- ifelse((p_adj<0.05 & abs(log2_fc)>1), "yes", "no")
sig <- c("yes"="darkgreen", "no"="lightgrey")

#top annotation
Volunteer <- c("V02" = "#FB9A99","V03" = "#E31A1C","V05" = "#A6CEE3", "V06" = "#1F78B4", "V07" = "#B2DF8A", "V09" = "#33A02C")
Timepoint <- c("Baseline"=time_col[4], "C10"=time_col[3], "DoD"=time_col[2], "T6"=time_col[1])
Batch <- c("batch_1"="lightgrey", "batch_2"="darkgrey")



combo_right_anno <-  rowAnnotation(gap = unit(2, "mm"),
                                   #space = rep("a", 24),
                                   #"log2FC" = anno_lines(log2_fc, which="row", axis = TRUE, ylim = c(-6, 6), axis_param = list(at=seq(-6, 6, by=2)), smooth = FALSE, add_points = TRUE),
                                   #"significant"=significant,
                                   "log2FC" = anno_barplot(log2_fc, which="row", axis = TRUE, ylim = c(-2, 6), axis_param = list(at=seq(-2, 6, by=2)),gp = gpar(fill = col_fun4(log2_fc))),
                                   
                                   #"p_adj"= anno_text(scales::scientific(p_adj, digits = 2), which="row"),
                                   width = unit(3, "cm") # width of the line graph
                                   #simple_anno_size = unit(2, "mm"), # width of the significance bar
)


combo_left_anno_var <-  rowAnnotation(gap = unit(5, "mm"),
                                      #annotation_name_gp = gpar(angle=45),
                                      show_annotation_name = FALSE,
                                      "significant"=significant,
                                      simple_anno_size = unit(2.5, "mm"), # width of the significance bar
                                      col=list("significant" = c("yes"="darkgreen", "no"="lightgrey")),
                                      annotation_legend_param = list(significant = list(title ="Significant",
                                                                                        at = rev(names(sig)),
                                                                                        #title_gp=gpar(angle=45),
                                                                                        legend_gp = gpar(fill = unname(sig)),
                                                                                        title_position = "topleft")
                                      )
                                      
)



combo_top_anno <- HeatmapAnnotation(gap = unit(2, "mm"), annotation_name_side = "left",
                                    Volunteer = rep(levels(cd$volunteer), 4),
                                    Timepoint = rep(levels(cd$timepoint), each=6),
                                    Batch = rep(rep(levels(cd$batch), each=3), 4),
                                    col=list(Batch = c("batch_1"="lightgrey", "batch_2"="darkgrey"),
                                             Timepoint = c("Baseline"=time_col[4], "C10"=time_col[3], "DoD"=time_col[2], "T6"=time_col[1]),
                                             Volunteer = c("V02" = "#FB9A99","V03" = "#E31A1C","V05" = "#A6CEE3", "V06" = "#1F78B4", "V07" = "#B2DF8A", "V09" = "#33A02C")
                                             
                                    ),
                                    simple_anno_size = unit(2.5, "mm"),
                                    annotation_legend_param = list(
                                      Volunteer = list(title = "Volunteer", at = names(Volunteer), legend_gp = gpar(fill = unname(Volunteer)), title_position = "topleft"),
                                      Timepoint = list(title ="Timepoint",at = names(Timepoint), legend_gp = gpar(fill = unname(Timepoint)), title_position = "topleft"),
                                      Batch = list(title = "Batch", at = names(Batch), legend_gp = gpar(fill = unname(Batch)), title_position = "topleft")
                                    )
)


#ht_opt("simple_anno_size" = unit(2, "mm"))

combo_map <- Heatmap(matrix = combo_matrix,
                     cluster_rows = FALSE,
                     name = "Normalised Frequency",
                     cluster_columns = FALSE,
                     row_names_side = "left",
                     col = col_fun4,
                     row_split = factor(rep(c("up", "down"), each = 12), levels = c("up", "down")),
                     rect_gp = gpar(col = "white"),
                     row_title = c("",""),
                     top_annotation = combo_top_anno,
                     right_annotation = combo_right_anno,
                     left_annotation = combo_left_anno_var,
                     show_heatmap_legend = TRUE,
                     column_names_rot = 45,
                     heatmap_legend_param = list(col = col_fun4, title = "Normalised Frequency", title_position = "topleft"),
                     width = unit(16, "cm"),
                     height = unit(16, "cm")
)

draw(combo_map,
     merge_legends = TRUE,
     #padding = unit(c(2, 20, 2, 2), "mm")
)


png("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/diffcyt/edgeR/improved_freq_hm_no_right.png", width=14, height=10, units = "in", res=400)
draw(combo_map,
     merge_legends = TRUE,
     #padding = unit(c(2, 20, 2, 2), "mm")
)
dev.off()
