  library(ggplot2)
  library(cowplot)
  library(vac69a.cytof)
  library(CATALYST)
  library(scales)

  
  #bigass color palette
  
  color_103_scheme <- c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
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
                        "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")
  
  inferno_mega_lite <- c("#000004", "#8A2267", "#EF802B", "#FFEC89", "#FCFFA4")
  sunset_white <- rev(c("#7D1D67", "#A52175", "#CB2F7A", "#ED4572", "#FA716C", "#FF9772", "#FFB985", "#FFD99F", "#FFFFFF"))
  
  inferno_white <- c("#FFFFFF", colorspace::sequential_hcl("inferno", n=8))
  
  UMAP_theme <- theme_minimal()+theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text = element_blank()
  )
  
  
  
  
  #functions, palettes etc. ####
  # 
  `%!in%` = Negate(`%in%`)
  # 3000 is the magic number
  #smol_daf <- read_small("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/", proportional=T, event_number=3000)
  
  #   ONLY RUN THIS WHEN CLUSTERING RESULTS ARE NEEDED, THE UMAP PROJECTION TAKES FOREVER (20+ MIN)
  smol_daf <- read_full("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")
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
  # umap_markers <- subset(refined_markers[,1], refined_markers[,1] %!in% c("CLA",  "CX3CR1", "CD95", "FoxP3", "Tbet", "Ki67", "HLADR", "ICOS", "PD1", "GZB"))
  
  # umap projections ####
  set.seed(1234);smol_daf <- scater::runUMAP(smol_daf,
                               subset_row=refined_markers[,1],
                               exprs_values = "exprs",
                               scale=T)
  
  #smol_daf <- CATALYST::filterSCE(smol_daf, timepoint!="C10")
  
  big_table <- prep_sce_for_ggplot(smol_daf)
  
  (test1 <- flo_umap(big_table, "CD38"))
  (test2 <- flo_umap(big_table, "ICOS"))
  
  
  
#data.table::fwrite(big_table, "~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_cells_with_UMAP_and_flo_merge_cluster.csv")
  
  #write.csv(big_table, "equal_small_vac69a_umap.csv")
  #write.csv(big_table, "proportional_small_vac69a_umap.csv")
  
 
ggsave(paste("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/proportional_contour_umap_14", "_bins",".png", sep=''), hex_through_time, height=6, width=9)

panelA_B  <- plot_grid(hex_through_time, gate_labeled_gg, ncol=2, rel_widths=c(3,1.5), labels = c("A", "B"))
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/panelA_B.png", panelA_B, height=6, width=11.25)

#ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/umap_through_time.png", smol_time12, width=15, height=5.54)
#ggsave("/home/flo/PhD/cytof/vac69a/figures_for_paper/umap_through_time.png", smol_time12, width=15, height=5.54)


# (big_time12 <- ggplot(big_table, aes(x=UMAP1, y=UMAP2))+
#   stat_density_2d(aes(fill = after_stat(nlevel)), geom="polygon", bins=25)+
#   scale_fill_gradientn(colors = inferno_mega_lite)+
#   xlim(c(-11, 14))+
#   ylim(c(-9, 10))+
#   theme_minimal()+
#   facet_grid(timepoint~volunteer, switch = "y")+
#   UMAP_theme+
#   theme(strip.text.y.left = element_text(size=14, angle = 0),
#         strip.text.x = element_text(size=14),
#         strip.placement = "outside")
    
# loop that makes the umap_through_time figures for each marker used in the clustering   
for (i in refined_markers[,1]){
  plt <- flo_umap(big_table, i, "timepoint")
  ggsave(paste("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/umap_through_time/", i, "_through_time.png", sep=''), plt, height=6, width=9)
  }      


# for (i in refined_markers[,1]){
#   plt <- flo_umap(down_big_table, i, "timepoint")
#   ggsave(paste("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/umap_through_time/downsampled_75_percent/", i, "_through_time.png", sep=''), plt, height=6, width=9)
# }     


# scatterplots colored by clusters ####

# add flo_merge cluster codes to daf
#coarse_merging_table <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/coarse_merge.csv", stringsAsFactors = F)
meta_45_mod_table <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/merging_tables/modified_ccp_meta45_merge.csv", header = T, stringsAsFactors = F)
merging_table1 <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/merging_tables/merging_table_april2020.csv", header=T, stringsAsFactors = F)

#get rid of spaces at beginning of string

merging_table1$new_cluster <- ifelse(substr(merging_table1$new_cluster, 1, 1)==" ", substr(merging_table1$new_cluster, 2, nchar(merging_table1$new_cluster)), merging_table1$new_cluster)
merging_table1$new_cluster <- ifelse(substr(merging_table1$new_cluster, 1, 1)==" ", substr(merging_table1$new_cluster, 2, nchar(merging_table1$new_cluster)), merging_table1$new_cluster)

merging_table1$new_cluster <- factor(merging_table1$new_cluster)

smol_merged_daf<- mergeClusters(smol_daf, k = "som100", table = meta_45_mod_table, id = "mod_meta45")
smol_merged_daf<- mergeClusters(smol_merged_daf, k = "mod_meta45", table = merging_table1, id = "flo_merge")


# turn into ggplottable object and add coumns for flo_merge name and whether it should be black or coloured

#big_table <- prep_sce_for_ggplot(sce = smol_merged_daf)
#big_table$flo_label <- as.character(cluster_ids(smol_merged_daf, k = "flo_merge"))

#data.table::fwrite(big_table, "~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_cells_with_UMAP_and_flo_merge_cluster.csv")
# making umap projections colored by cluster identity ####

big_table <- data.table::fread("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_cells_with_UMAP_and_flo_merge_cluster.csv")

# ATTENTION: for some reason in the writing of this there were some spaces at the end of some names so if the cluster coloring
# doesn't work check that the significant column contains all the significant clusters!
sig_clusters <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/sig_t6_clusters.csv", header = TRUE, stringsAsFactors = FALSE)

#get rid of resting Vd cluster
sig_clusters <- sig_clusters[-9,]


big_table$significant <- ifelse(big_table$flo_label %in% sig_clusters[,2], big_table$flo_label, "black")
big_table$alpha <- ifelse(big_table$flo_label %in% sig_clusters[,2], 1, 0.5)



#define the colour palette
cluster_palette <- color_103_scheme[2:(length(unique(big_table$significant))+1)]
cluster_palette[5] <- "#40fd14"
#cluster_palette[2] <- "#006FA6"
cluster_palette[6] <- "#0000A6"
cluster_palette[8] <- "#FFA500"
cluster_palette[9] <- "red"
cluster_palette[1] <- "#000000"
cluster_palette[10] <- "#4dabf5" 
names(cluster_palette)[1] <- "black"
names(cluster_palette)[2:length(cluster_palette)] <- unique(big_table$significant)[-match("black", unique(big_table$significant))]
cluster_palette[5:9] <- rev(cluster_palette[5:9])

# write.csv(cluster_palette,"/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cluster_palette.csv", )
# colcsv <- read.csv("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cluster_palette.csv", header=T, stringsAsFactors = F)
# 
# col_pal <- abvc$x
# names(col_pal) <- abvc$X

# restrict data to T6 and downsample by 33% to make it look nice
short_big_table_t6 <- subset(big_table, big_table$timepoint=="T6")
short_big_table_t6 <- short_big_table_t6[seq(1,nrow(short_big_table_t6), by=3), ]

t6_sig_clusters_umap <- ggplot(short_big_table_t6, aes(x=UMAP1, y=UMAP2))+
  geom_point(aes(color=significant, alpha=alpha), shape=".")+
  #UMAP_theme+
  scale_color_manual(values = cluster_palette)+
  theme_minimal()+
  ggtitle("Differentially Abundant\nClusters at T6")+
  UMAP_theme+
  guides(colour = guide_legend(override.aes = list(size = 4)),
         alpha= "none")+
  # scale_y_continuous(limits = c(-11.2, 11.3))+
  # scale_x_continuous(limits=c(-13, 10))+
  #guides(alpha = guide_legend(override.aes = list(size = 10)))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust=0.5))+
  coord_cartesian(xlim=c(-13, 10),
                  ylim=c(-11.2, 11.3))

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/t6_sig_clusters_umap_var.png",  t6_sig_clusters_umap, height=5, width=5)

cd38_plot <- flo_umap(short_big_table_t6, "CD38")+theme(axis.title = element_blank())
bcl2_plot <- flo_umap(short_big_table_t6, "BCL2")+theme(axis.title = element_blank())



panels_d_e <- plot_grid(t6_sig_clusters_umap, cd38_plot, bcl2_plot, nrow=1)

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/panels_d_e.png",  panels_d_e, height=4, width=12)

new_panels_d_e <- plot_grid(gate_labeled_gg, t6_sig_clusters_umap, cd38_plot, bcl2_plot, nrow=2)
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/new_panels_d_e.png",  new_panels_d_e, height=8, width=8)


new_panels_d_ever2 <- plot_grid(empty_gates_plot,  cd38_plot, bcl2_plot, nrow=1)
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/new_panels_d_e_ver2.png",  new_panels_d_ever2, height=4, width=12)


vertical_d_e <- plot_grid(gate_labeled_gg, t6_sig_clusters_umap, ncol = 1)
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/vertical_d_e.png", vertical_d_e, height = 8, width = 4)


#all clusters colored in
expanded_cluster_palette <- c(cluster_palette, color_103_scheme[22:46])

names(expanded_cluster_palette)[11:length(expanded_cluster_palette)] <- unique(big_table$flo_label)[-match(unique(big_table$significant), unique(big_table$flo_label), nomatch = 0)]

t6_all_clusters_umap <- ggplot(short_big_table_t6, aes(x=UMAP1, y=UMAP2))+
    #geom_point(aes(color=flo_label, alpha=scales::rescale(alpha, to=c(0.1, 1))), shape=".")+
    geom_point(aes(color=flo_label), shape=".", alpha=0.4)+
    #UMAP_theme+
    scale_color_manual(values = expanded_cluster_palette)+
    theme_minimal()+
    ggtitle("All Cells at T6 coloured\nby Cluster ID")+
    UMAP_theme+
    guides(colour = guide_legend(override.aes = list(size = 2, shape=16, alpha=1), ncol = 3),
           alpha= "none")+
    # scale_y_continuous(limits = c(-11.2, 11.3))+
    # scale_x_continuous(limits=c(-13, 10))+
    #guides(alpha = guide_legend(override.aes = list(size = 10)))+
    theme(legend.position = "none",
          legend.text = element_text(size=6),
          axis.title = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(hjust=0.5),
          #legend.key.width=unit(4,"cm"),
          )+
    coord_cartesian(xlim=c(-13, 10),
                    ylim=c(-11.2, 11.3))

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/t6_all_clusters_umap_without_legend.png",  t6_all_clusters_umap, height=4, width=4)

#ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/t6_all_clusters_umap_with_legend.png",  t6_all_clusters_umap, height=8, width=6)

#first: include my merge to DAF

#get rid of spaces at beginning of string
library(CATALYST)
meta_45_mod_table <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/merging_tables/modified_ccp_meta45_merge.csv", header = T, stringsAsFactors = F)
merging_table1 <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/merging_tables/merging_table_april2020.csv", header=T, stringsAsFactors = F)
vis_merge <- read.csv("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/figures_for_phil/merge_for_vis.csv", header=T, stringsAsFactors = F)

#get rid of spaces at beginning of string

merging_table1$new_cluster <- ifelse(substr(merging_table1$new_cluster, 1, 1)==" ", substr(merging_table1$new_cluster, 2, nchar(merging_table1$new_cluster)), merging_table1$new_cluster)
merging_table1$new_cluster <- ifelse(substr(merging_table1$new_cluster, 1, 1)==" ", substr(merging_table1$new_cluster, 2, nchar(merging_table1$new_cluster)), merging_table1$new_cluster)

merging_table1$new_cluster <- factor(merging_table1$new_cluster)

merged_daf<- mergeClusters(daf, k = "som100", table = meta_45_mod_table, id = "mod_meta45")
merged_daf<- mergeClusters(merged_daf, k = "mod_meta45", table = merging_table1, id = "flo_merge")



vis_merge$new_cluster <- factor(vis_merge$new_cluster)
merged_daf <- mergeClusters(merged_daf, k = "meta45", table = vis_merge, id = "vis_merge")





# make a dictionary for assigning colors to cluster numbers for meta45
color_dic <- color_103_scheme[1:45]
names(color_dic) <- c(1:45)

# make a dictionary for assigning colors to cluster numbers from vis_merge, while preserving names from
# flo_merge
ordered_meta_dic <- merging_table1[order(merging_table1$old_cluster),]
final_dic <- ordered_meta_dic[!duplicated(ordered_meta_dic$new_cluster),]



flumap <- plotDR(merged_daf, color_by = "vis_merge")+
  facet_wrap("timepoint")+
  scale_color_manual(values = color_dic, labels=as.character(final_dic$new_cluster))+
  theme(legend.position = "right", 
        legend.title = element_blank())


meta45_map <- plotDR(merged_daf, color_by = "meta45")+
  facet_wrap("timepoint")+
  scale_color_manual(values = color_dic)+
  theme(legend.position = "right",
        legend.title = element_blank())

meta_leg <-get_legend(meta45_map)
leg <- get_legend(flumap)  

(right_col <- cowplot::plot_grid(meta_leg, leg, ncol=1, axis= "rl", rel_widths=c(2,1)))
(left_col <- cowplot::plot_grid(meta45_map, flumap, ncol=1, axis = "trbl"))
both <- cowplot::plot_grid(left_col, leg, ncol=1, rel_heights = c(5,1), align="h", axis = "l")
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/figures_for_phil/both_umaps.png", top_row, width = 12, height=12)

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/figures_for_phil/flo_merge_umap.png", flumap, width = 13, height=6)
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/figures_for_phil/meta45_map.png", meta45_map, width = 12, height=6)


flo_merge_map <- plotDR(merged_daf, color_by = "flo_merge")+
  facet_wrap("timepoint")+
  #scale_color_manual(values=color_dic, labels=as.character(merging_table1$new_cluster))+
  scale_color_manual(values=unique(color_dic))+
  theme(legend.position = "bottom", 
        legend.title = element_text())

# plot with points and contour lines

# smol_time+
#   stat_density_2d(aes(fill = stat(nlevel)), geom="polygon", bins=11)+
#   stat_density_2d(size=0.13, color="black", bins=11)+
#   #stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon", bins=22, size=0.3, color="black", contour=T)+
#   # scale_fill_viridis(option="A")+
#   scale_color_viridis(option="A", direction=1)+
#   #scale_fill_viridis(option="B", direction=1)+
#   scale_fill_gradient(low = "white", high="black")+
#   geom_point(shape=".")+
#   xlim(c(-14, 11))+
#   ylim(c(-10, 10))+
#   theme_minimal()+
#   UMAP_theme






# make a bunch of contour plots with lapply ####



list_of_data <- split(big_table, big_table$timepoint)  

list_of_plots <- lapply(list_of_data, function(x){
  ggplot(x, aes(x=UMAP1, y=UMAP2))+
    #stat_density_2d(aes(fill = ..density..), geom = 'raster', contour = FALSE, n = 1500)+
    #stat_density_2d(contour = TRUE, bins=14, color="white", size=0.1)+
    xlim(c(-13, 10))+
    ylim(c(-9.5, 13))+
    ggtitle(unique(x$timepoint))+
    theme_minimal()+
    #facet_wrap(~timepoint)+
    UMAP_theme+
    theme(panel.grid.major = element_blank(),
          plot.title = element_text(hjust=0.5, size=14),
          plot.title.position = "bottom")+
    scale_fill_gradientn(colors=inferno_white)
})
#scale_fill_viridis(option = "A")

# don't try to render this plot in rstudio- jsut write out and look at it with
# photo viewer- way faster!

plot_grid(plotlist=list_of_plots, label)

dir_name <- "/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/"

system.time(
  lapply(list_of_plots, function(x){
    ggsave(
      paste0(dir_name, "proportional_contour_umap_", unique(x$data$timepoint), ".png", sep=''),
      x, height=6, width=9)
  })
)

ggsave(
  paste(dir_name, "labeled_contour_umap_through_time", ".png", sep=''),
  plot_grid(plotlist=list_of_plots, labels="AUTO", nrow = 1), width=5.8, height=3)


# a bunch of multiplanel plots for related markers ####

cd4_plot <- flo_umap(big_table, "CD4")
cd8_plot <- flo_umap(big_table, "CD8")
vd2_plot <- flo_umap(big_table, "Vd2")
va72_plot <- flo_umap(big_table, "Va72")


lineage_plot <- plot_grid(cd4_plot, cd8_plot, vd2_plot, va72_plot, ncol=2)
# ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/lineage_plot_ring.png", lineage_plot)
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/lineage_plot.png", lineage_plot)

cd38_plot <- flo_umap(big_table, "CD38")
hladr_plot <- flo_umap(big_table, "HLADR")
bcl2_plot <- flo_umap(big_table, "BCL2")
cd27_plot <- flo_umap(big_table, "CD27")


# cd38_bcl2_plot <-   plot_grid(cd38_plot, bcl2_plot, ncol=1)
# activation_plot  <- plot_grid(cd38_plot, bcl2_plot,  hladr_plot, cd27_plot, ncol=2)
# # ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/activation_plot.png", activation_plot)
# ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/activation_plot.png", activation_plot)
# ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cd38_bcl2_plot.png", cd38_bcl2_plot)


cd45ro_plot <- flo_umap(big_table, "CD45RO")
ccr7_plot <- flo_umap(big_table, "CCR7")
cd57_plot <- flo_umap(big_table, "CD57")
cd45ra_plot <- flo_umap(big_table, "CD45RA")

memory_plot <- plot_grid(cd45ro_plot, ccr7_plot, cd57_plot, cd45ra_plot, ncol=2)
# ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/memory_plot.png", memory_plot)
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/memory_plot.png", memory_plot)




cd25_plot <- flo_umap(big_table, "CD25")
cd127_plot <- flo_umap(big_table, "CD127")
foxp3_plot <- flo_umap(big_table, "FoxP3")
cd39_plot <- flo_umap(big_table, "CD39")

treg_plot <- cowplot::plot_grid(cd25_plot, cd127_plot, foxp3_plot, cd39_plot, ncol=2)
# ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/treg_plot.png", treg_plot, width=8, height=5.54)
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/treg_plot.png", treg_plot)






cd38_plot_through_time <- flo_umap(big_table, "CD38", "timepoint")
hladr_plot_through_time <- flo_umap(big_table, "HLADR", "timepoint")
bcl2_plot_through_time <- flo_umap(big_table, "BCL2", "timepoint")
cd27_plot_through_time <- flo_umap(big_table, "CD27", "timepoint")

activation_plots_through_time  <- list(cd38_plot_through_time, bcl2_plot_through_time,  hladr_plot_through_time, cd27_plot_through_time)
ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/activation_plot.png", activation_plot)


lapply(activation_plots_through_time, function(x){
  ggsave(paste("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/", x$labels$title ,"_through_time.png", sep=''), x,
         width = 8, height=4.5)
})



GZB_plot_through_time <- flo_umap(big_table, "GZB", "timepoint")
perforin_plot_through_time <- flo_umap(big_table, "Perforin", "timepoint")
Tbet_plot_through_time <- flo_umap(big_table, "Tbet", "timepoint")
ki67_plot_through_time <- flo_umap(big_table, "Ki67", "timepoint")

intracellular_plots_through_time  <- list(GZB_plot_through_time, perforin_plot_through_time,  Tbet_plot_through_time, ki67_plot_through_time)
# ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/activation_plot.png", activation_plot)


lapply(intracellular_plots_through_time, function(x){
  ggsave(paste("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/", x$labels$title ,"_through_time.png", sep=''), x,
         width = 8, height=4.5)
})


ctla4_plot_through_time <- flo_umap(big_table, "CTLA4", "timepoint")
pd1_plot_through_time <- flo_umap(big_table, "PD1", "timepoint")
CD28_plot_through_time <- flo_umap(big_table, "CD28", "timepoint")
CD95_plot_through_time <- flo_umap(big_table, "CD95", "timepoint")

exhaustion_plots_through_time  <- list(ctla4_plot_through_time, pd1_plot_through_time,  CD28_plot_through_time, CD95_plot_through_time)


lapply(exhaustion_plots_through_time, function(x){
  ggsave(paste("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/", x$labels$title ,"_through_time.png", sep=''), x,
         width = 8, height=4.5)
})




# loop to make a bunch of beautiful contour plots ####

n_bins <- c(8, 9, 11, 13)
n_bins <- 13
#system.time(for (i in n_bins){
  
plt <- ggplot(down_big_table, aes(x=UMAP1, y=UMAP2))+
  stat_density_2d(aes(fill = ..density..), geom = 'raster', contour = FALSE, n = 1500)+
  stat_density_2d(contour = TRUE, bins=13, color="white", size=0.05)+
  scale_x_continuous(limits=c(-14, 10), breaks = seq(-10, 10, by=5))+
  scale_y_continuous(limits=c(-11.2, 11.3), breaks = seq(-10, 10, by=5))+
  theme_minimal()+
  facet_wrap(~timepoint, ncol = 4)+
  UMAP_theme+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size=14),
        axis.title = element_text())+
  scale_fill_gradientn(colors=inferno_white)

# ggsave(paste0("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/countour_with_", i ,"_bins.png", sep=''), plt,
#        height = 4, width = 16)
# 
# })
