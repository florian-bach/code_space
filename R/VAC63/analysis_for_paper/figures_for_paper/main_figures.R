#stacked bar charts & lollipop bars####

library(dplyr)
library(tidyr)
library(ggplot2)


#colour business

colcsv <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_colours.csv")

col_pal <- colcsv$colour
names(col_pal) <- colcsv$cluster_id

lineage_palette <- c("#AA3377", "#EE6677", "#4477AA", "#66CCEE", "#228833", "#FFA500", "pink")
names(lineage_palette) <-c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")

prim_sig_clusters <- read.delim("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/prim_ter_sig_clusters.csv", header=FALSE)$V1





stacked_bar_data <-  read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_all_cluster_freqs.csv")


stacked_bar_data$volunteer <- factor(stacked_bar_data$volunteer, levels=c("v313", "v315", "v320", "v306", "v301", "v308", "v305", "v304", "v310"))

stacked_bar_data$cluster_id <- gsub("activated CD8lo Effector", "activated DN", stacked_bar_data$cluster_id)


activated_barchart_data <- subset(stacked_bar_data, stacked_bar_data$cluster_id %in% prim_sig_clusters)

activated_barchart_data$lineage <- substr(activated_barchart_data$cluster_id, nchar(activated_barchart_data$cluster_id)-2, nchar(activated_barchart_data$cluster_id))
#[1] "CD4" " gd" "AIT" "NKT" "CD8" " DN" "reg"

lin_replacement <- setNames(c("CD4", "gd", "MAIT", "NKT", "CD8", "DN", "Treg"), unique(activated_barchart_data$lineage))
activated_barchart_data$lineage <- stringr::str_replace_all(activated_barchart_data$lineage, lin_replacement)


summary <- activated_barchart_data %>%
  group_by(timepoint, cluster_id) %>%
  mutate(mean=base::mean(frequency), sd=stats::sd(frequency)) %>%
  ungroup()%>%
  mutate(timepointf=factor(timepoint, levels=c("Baseline", "DoD", "T6", "C45")))

#tail(summary)

summary$lineage <- factor(summary$lineage, levels=rev(c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")))

prim_summary_t6 <- subset(summary_t6, summary_t6$n_infection=="first" & summary_t6$cluster_id %in% sig_prim_t6_clusters)



cluster_order <-  prim_summary_t6 %>%
  group_by(cluster_id) %>%
  mutate("sum"=sum(frequency)) %>%
  arrange(desc(sum)) %>%
  select(cluster_id)

cluster_order <- rev(unique(cluster_order$cluster_id))



summary$n_infection <- ifelse(summary$volunteer %in% c("v313", "v315", "v320"), "first", "third")


summary <- filter(summary, cluster_id %in% sig_clusters)

activation_stacked_barchart_cluster <- ggplot(summary, aes(x=volunteer, y=frequency/100, fill=factor(cluster_id, levels=c(cluster_order))))+
  geom_bar(stat="identity", position="stack")+
  #geom_text(aes(y=(total_cd3/100)+0.01, label=paste0(round(total_cd3, digits = 1), "%", sep='')))+
  theme_minimal()+
  facet_wrap(~timepointf, strip.position = "bottom", ncol=4)+
  ggtitle("Overall T cell activation")+
  scale_fill_manual(values=col_pal)+
  # scale_fill_manual(values=lineage_palette, labels=rev(c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "NKT")))+
  scale_y_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid.minor.y = element_blank(),
        strip.placement = "outside")

ggsave("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/figures_for_paper/activation_stack_cluster_id.pdf", activation_stacked_barchart_cluster, height=4, width=8)


activation_stacked_barchart_lineage <- ggplot(summary, aes(x=volunteer, y=frequency/100, fill=lineage))+
  geom_bar(stat="identity", position="stack")+
  #geom_text(aes(y=(total_cd3/100)+0.01, label=paste0(round(total_cd3, digits = 1), "%", sep='')))+
  theme_minimal()+
  facet_wrap(~timepointf, strip.position = "bottom", ncol=4)+
  ggtitle("Overall T cell activation")+
  scale_fill_manual(values=lineage_palette)+
  # scale_fill_manual(values=lineage_palette, labels=rev(c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "NKT")))+
  scale_y_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid.minor.y = element_blank(),
        strip.placement = "outside")

ggsave("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/figures_for_paper/activation_stack_lineage.pdf", activation_stacked_barchart_lineage, height=4, width=8)






#prepare data for lollipop plots: only T6; calculate average between first and thirds
sig_clusters <- read.delim("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/prim_ter_sig_clusters.csv", header=FALSE)$V1

summary_t6 <- subset(summary, summary$timepoint=="T6")
summary_t6$scaled_freq <- ifelse(summary_t6$n_infection=="third",summary_t6$frequency/6, summary_t6$frequency/3)
sig_summary_t6 <- subset(summary_t6, summary_t6$cluster_id %in% sig_clusters)


sig_ter_t6_clusters <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/edgeR/ter_t6_df_edger.csv", header = TRUE, stringsAsFactors = FALSE)
sig_prim_t6_clusters <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/edgeR/prim_t6_df_edger.csv", header = TRUE, stringsAsFactors = FALSE)

sig_ter_t6_clusters <- subset(sig_ter_t6_clusters, sig_ter_t6_clusters$p_adj<0.05 & sig_ter_t6_clusters$logFC>1)$cluster_id
sig_prim_t6_clusters <- subset(sig_prim_t6_clusters, sig_prim_t6_clusters$p_adj<0.05 & sig_prim_t6_clusters$logFC>1)$cluster_id

#write.table(unique(c(sig_prim_t6_clusters,sig_ter_t6_clusters)), "prim_ter_sig_clusters.txt", sep="\t", col.names = FALSE, row.names = FALSE)

#make dataframe for first infection
prim_summary_t6 <- subset(summary_t6, summary_t6$n_infection=="first" & summary_t6$cluster_id %in% sig_prim_t6_clusters)


#make dataframe for third infection and add channel that only contains color hex codes for sig clusters, otherwise paste blank
ter_summary_t6 <- subset(summary_t6, summary_t6$n_infection=="third" & summary_t6$cluster_id %in% sig_prim_t6_clusters)
ter_summary_t6$fill <- ifelse(ter_summary_t6$cluster_id %in% sig_ter_t6_clusters, ter_summary_t6$cluster_id, "blank")


cluster_order <- prim_summary_t6 %>%
  group_by(cluster_id) %>%
  mutate("sum"=sum(frequency)) %>%
  arrange(desc(sum)) %>%
  select(cluster_id)

cluster_order <- rev(unique(cluster_order$cluster_id))

write.csv(cluster_order, "cluster_order.csv", row.names = FALSE)



prim_sig_activation_lolli <- ggplot()+
  geom_bar(data=prim_summary_t6,aes(y=factor(cluster_id, levels=cluster_order), x=scaled_freq/100, fill=factor(cluster_id, levels=cluster_order)), stat="identity", position="stack", width=1)+
  theme_minimal()+
  ggtitle("differentially UP at T6 first")+
  scale_fill_manual(values=col_pal)+
  scale_x_reverse(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1), limits = c(0.05, 0), breaks=seq(1,5)/100, expand = c(0,0))+
  guides(fill=guide_legend(reverse = TRUE, nrow = 5,keyheight = unit(2, "mm"), keywidth = unit(4, "mm")))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.placement = "outside",
        plot.margin = unit(c(1,0,1,5), "mm"),
        legend.position = "bottom",
        legend.text = element_text(size=6),
        legend.title = element_blank())


cluster_leg <- cowplot::get_legend(prim_sig_activation_lolli)

prim_sig_activation_lolli <- prim_sig_activation_lolli+theme(legend.position = "none")


ter_sig_activation_lolli <-   ggplot()+
  geom_bar(data=ter_summary_t6,aes(y=factor(cluster_id, levels=cluster_order), x=scaled_freq/100, fill=factor(fill, levels=cluster_order)), stat="identity", position="stack", width=1)+
  theme_minimal()+
  ggtitle("differentially UP at T6 third")+
  scale_fill_manual(values=col_pal)+
  scale_x_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1), limits = c(0,0.05), breaks=seq(1,5)/100 ,expand = c(0,0))+
  guides(fill=guide_legend(reverse = TRUE, nrow = 5,keyheight = unit(2, "mm"), keywidth = unit(4, "mm")))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.placement = "outside",
        plot.margin = unit(c(1,5,1,0), "mm"),
        legend.position = "none",
        legend.text = element_text(size=6),
        legend.title = element_blank())



combo_lolli <- cowplot::plot_grid(prim_sig_activation_lolli, ter_sig_activation_lolli)

combo_lolli <- cowplot::plot_grid(combo_lolli, cluster_leg, ncol=1, rel_heights = c(5,1), rel_widths = c(1,1,2),align="v", axis = "tbrl")
ggsave("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/figures_for_paper/all_activation_lollipop.pdf", combo_lolli, height=4, width=7)



# pie charts ####

library(gridBase)
library(grid)
library(circlize)
library(ComplexHeatmap)
library(dplyr)


#colour stuff

colcsv <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_colours.csv")

col_pal <- colcsv$colour
names(col_pal) <- colcsv$cluster_id

lineage_palette <- c("#AA3377", "#EE6677", "#4477AA", "#66CCEE", "#228833", "#FFA500", "pink")
names(lineage_palette) <-c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")



stacked_bar_data <-  read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_all_cluster_freqs.csv")
stacked_bar_data$cluster_id <- gsub("activated CD8lo Effector", "activated DN", stacked_bar_data$cluster_id)


prim_sig_clusters <- read.delim("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/prim_ter_sig_clusters.csv", header=FALSE)$V1

stacked_bar_data$lineage <- substr(stacked_bar_data$cluster_id, nchar(stacked_bar_data$cluster_id)-2, nchar(stacked_bar_data$cluster_id))
#[1] "CD4" " gd" "AIT" "NKT" "CD8" " DN" "reg"

lin_replacement <- setNames(c("CD4", "gd", "MAIT", "NKT", "CD8", "DN", "Treg", "gd"), unique(stacked_bar_data$lineage))
stacked_bar_data$lineage <- stringr::str_replace_all(stacked_bar_data$lineage, lin_replacement)

stacked_bar_data$n_infection <- ifelse(stacked_bar_data$volunteer %in% c("v313", "v315", "v320"), "first", "third")


prim_map_data <- stacked_bar_data %>%
  filter(timepoint=="T6", n_infection=="first") %>%
  filter(cluster_id %in% prim_sig_clusters) %>%
  group_by(cluster_id, timepoint) %>%
  mutate(mean_freq=mean(frequency)) %>%
  ungroup() %>%
  select(cluster_id, n_infection, mean_freq, lineage)

prim_summary <- prim_map_data[!duplicated(prim_map_data), ]


prim_summary$lineage <- factor(prim_summary$lineage, levels=rev(c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")))




prim_pie_data <- prim_summary
ordered_col_pal <- col_pal[match(prim_pie_data$cluster_id, names(col_pal))]


prim_pie_data <- prim_pie_data[!duplicated(prim_pie_data), ]
prim_pie_data <- prim_pie_data[order(prim_pie_data$lineage, decreasing = T),]


prim_circlize_plot = function() {
  
  circos.par("track.height" = 0.6, start.degree = 90)
  circos.par(cell.padding = c(0.02, 0, 0.02, 0))
  
  circos.initialize(factors = prim_pie_data$cluster_id, xlim = c(0,1), sector.width = prim_pie_data$mean_freq)
  
  o.cell.padding = circos.par("cell.padding")
  
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.2, bg.border = NA,
                         cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4]))
  
  circos.track(factors = prim_pie_data$cluster_id,
               ylim = c(0, 1),
               x=prim_pie_data$mean_freq,
               bg.col = col_pal[match(prim_pie_data$cluster_id, names(col_pal))])
  #circos.genomicLabels()
  
  
  #no color, just text
  # c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")
  
  highlight.sector(sector.index = prim_pie_data$cluster_id[prim_pie_data$lineage=="CD4"], track.index = 1,
                   border = "black", text = "CD4", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = prim_pie_data$cluster_id[prim_pie_data$lineage=="Treg"], track.index = 1,
                   border = "black", text = "Treg", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = prim_pie_data$cluster_id[prim_pie_data$lineage=="CD8"], track.index = 1,
                   border = "black", text = "CD8", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = prim_pie_data$cluster_id[prim_pie_data$lineage=="MAIT"], track.index = 1,
                   border = "black", text = "MAIT", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = prim_pie_data$cluster_id[prim_pie_data$lineage=="gd"], track.index = 1,
                   border = "black", text = "gd", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = prim_pie_data$cluster_id[prim_pie_data$lineage=="DN"], track.index = 1,
                   border = "black", text = "DN", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = prim_pie_data$cluster_id[prim_pie_data$lineage=="NKT"], track.index = 1,
                   border = "black", text = "NKT", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  circos.clear()
  
}

prim_lvl_data <- prim_pie_data$cluster_id


prim_leg_col_pal <- subset(col_pal, names(col_pal) %in% prim_lvl_data)

# discrete
prim_lgd_cluster= Legend(at = prim_lvl_data, type = "grid", 
                         legend_gp = gpar(fill = prim_leg_col_pal[prim_lvl_data], size=0.5),
                         title_position = "topleft",
                         labels_gp = gpar(fontsize = 5),
                         title = "Cluster ID")
# 
# 
# png("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/ter_cluster_activation_pie.png", width=7, height=4, units = "in", res=400)
# 
# plot.new()
# circle_size = unit(1, "snpc") # snpc unit gives you a square region
# pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
#                       just = c("left", "center")))
# par(omi = gridBase::gridOMI(), new = TRUE)
# circlize_plot()
# par(xpd = NA)
# upViewport()
# draw(lgd_cluster, x = circle_size, just = "left")
# #title("Average Size of Significantly Upregulated Clusters at T6", 0,6)
# 
# dev.off()



# tertiaries #


ter_t6_df <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/edgeR/ter_t6_df_edger.csv")

ter_sig_clusters <- ter_t6_df %>%
  filter(p_adj <= 0.05, logFC >=1)%>%
  select(cluster_id)



ter_map_data <- stacked_bar_data %>%
  filter(timepoint=="T6", n_infection=="third") %>%
  filter(cluster_id %in% ter_sig_clusters$cluster_id) %>%
  group_by(cluster_id, timepoint, n_infection) %>%
  mutate(mean_freq=mean(frequency)) %>%
  ungroup() %>%
  select(cluster_id, n_infection, mean_freq, lineage)

ter_summary <- ter_map_data[!duplicated(ter_map_data), ]




ter_summary$lineage <- factor(ter_summary$lineage, levels=rev(c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")))




ter_pie_data <- filter(ter_summary, n_infection=="third")

ordered_col_pal <- col_pal[match(ter_pie_data$cluster_id, names(col_pal))]



ter_pie_data <-ter_pie_data[!duplicated(ter_pie_data), ]

ter_pie_data <- ter_pie_data[order(ter_pie_data$lineage, decreasing = T),]


#circos.par(cell.padding = c(0, 0))


ter_circlize_plot = function() {
  
  
  circos.par("track.height" = 0.6, start.degree = 90)
  circos.par(cell.padding = c(0.02, 0, 0.02, 0))
  
  circos.initialize(factors = ter_pie_data$cluster_id, xlim = c(0,1), sector.width = ter_pie_data$mean_freq)
  
  o.cell.padding = circos.par("cell.padding")
  
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.2, bg.border = NA,
                         cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4]))
  
  circos.track(factors = ter_pie_data$cluster_id,
               ylim = c(0, 1),
               x=ter_pie_data$mean_freq,
               bg.col = ordered_col_pal[match(ter_pie_data$cluster_id, names(ordered_col_pal))])

  
  highlight.sector(sector.index = ter_pie_data$cluster_id[ter_pie_data$lineage=="CD4"], track.index = 1,
                   border = "black", text = "CD4", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = ter_pie_data$cluster_id[ter_pie_data$lineage=="CD8"], track.index = 1,
                   border = "black", text = "CD8", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = ter_pie_data$cluster_id[ter_pie_data$lineage=="MAIT"], track.index = 1,
                   border = "black", text = "MAIT", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = ter_pie_data$cluster_id[ter_pie_data$lineage=="gd"], track.index = 1,
                   border = "black", text = "gd", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = ter_pie_data$cluster_id[ter_pie_data$lineage=="DN"], track.index = 1,
                   border = "black", text = "DN", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = ter_pie_data$cluster_id[ter_pie_data$lineage=="NKT"], track.index = 1,
                   border = "black", text = "NKT", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  circos.clear()
  
}




pdf("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/figures_for_paper/combo_cluster_activation_pie.pdf", width=8, height=8, onefile = FALSE)

plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridBase::gridOMI(), mfrow = c(1, 3))
prim_circlize_plot()
title("First Infection",  line = -15)
ter_circlize_plot()
title("Third Infection", line = -15)
upViewport()
draw(prim_lgd_cluster, x = circle_size, just = "right")
dev.off()


# spicy UMAPs ####

library(ggplot2)
library(cowplot)
library(vac69a.cytof)
library(scales)
library(dplyr)


big_table <- data.table::fread("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_all_Tcells_with_UMAP.csv")


t6_edger <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/edgeR/t6_edgeR.csv", header = TRUE, stringsAsFactors = FALSE)
sig_t6_clusters <- subset(t6_edger, t6_edger$p_adj<0.05 & abs(t6_edger$logFC)>1)$cluster_id

sig_ter_t6_clusters <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/edgeR/ter_t6_df_edger.csv", header = TRUE, stringsAsFactors = FALSE)
sig_prim_t6_clusters <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/edgeR/prim_t6_df_edger.csv", header = TRUE, stringsAsFactors = FALSE)

sig_ter_t6_clusters <- subset(sig_ter_t6_clusters, sig_ter_t6_clusters$p_adj<0.05 & sig_ter_t6_clusters$logFC>1)$cluster_id
sig_prim_t6_clusters <- subset(sig_prim_t6_clusters, sig_prim_t6_clusters$p_adj<0.05 & sig_prim_t6_clusters$logFC>1)$cluster_id

#get rid of resting Vd cluster

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

big_table$timepoint <- factor(big_table$timepoint, levels=c("Baseline", "DoD","T6","C45"))

big_table$flo_label <- gsub("activated RA-RO- DN", "activated DN", big_table$flo_label)


big_table$prim_significant <- ifelse(big_table$flo_label %in% sig_prim_t6_clusters, big_table$flo_label, "black")
big_table$prim_alpha <- ifelse(big_table$flo_label %in% sig_prim_t6_clusters, 1, 0.6)

big_table$ter_significant <- ifelse(big_table$flo_label %in% sig_ter_t6_clusters, big_table$flo_label, "black")
big_table$ter_alpha <- ifelse(big_table$flo_label %in% sig_ter_t6_clusters, 1, 0.6)



#reduce table size by 66% 
short_big_table <- big_table[seq(1,nrow(big_table), by=3),] 



#add a column whether this cluster was significant


cluster_names <- c(unique(short_big_table$flo_label))
cols <- rev(c(color_103_scheme[40:88]))
names(cols)=cluster_names


#read in pretty divergent colours from vivax paper
colcsv <- read.csv("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cluster_palette.csv", header=T, stringsAsFactors = F)

#find all the activated clusters, replace the colours of the first nine with the pretty vivax colours minus the activated
#gamma delta pink (to avoid duplication)

vivax_colours <- colcsv$x

vivax_colours <- vivax_colours[-1]

cols[grep("activated", names(cols))[1:9]] <- vivax_colours

#light grey to
cols[grep("activated", names(cols))][15] <- "#FF913F"
#almost black to
cols[grep("activated", names(cols))][13] <- "#3B5DFF"
#light lime
cols[grep("activated", names(cols))][14] <-"#A4E804"
#dark blue
cols[grep("activated", names(cols))][10] <-"#FFFF00"
#dark yellow
cols[grep("activated", names(cols))][5] <-"#D16100"




#change pink, swap with green
cols[grep("activated ICOS- EM CD4", names(cols))] <- "#fa0064"
cols[grep("activated EM Treg", names(cols))] <- "#008941"

#change cd8 orange to something a bit redder
#cols[grep("activated HLADR+ EM CD8", names(cols))] <- "#ff623f"

# v dark blue for naive cd8s
cols[grep("Naive CD8", names(cols))] <-"#010048"
# change dn to nicer yellow #ffc900
cols[grep("activated DN", names(cols))] <- "#ffc900"

#shuffle some gd colours around
cols[grep("gd", names(cols))][6] <- "#800000"
cols[grep("gd", names(cols))][3] <- "#ff623f"

#shuffle some cd4 colours around
cols[grep("activated ICOS- EM CD4", names(cols))] <- "#004d4d"

cols[30] <- "#6e0280"

names(cols)=cluster_names






cols <- c(cols, "black"="black")

names(cols) <- gsub("activated RA-RO- DN", "activated DN", names(cols))

cluster_colors <- data.frame("cluster_id"=names(cols), "colour"=unname(cols))
cluster_colors$cluster_id <- gsub("activated RA-RO- DN", "activated DN", cluster_colors$cluster_id)

write.csv(cluster_colors, "~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_colours.csv", row.names =FALSE)




shorter_big_table_t6 <- big_table %>%
  filter(timepoint=="T6") %>%
  group_by(n_infection) %>%
  sample_n(9*10^4) %>%
  #sample_n(0.1*10^4) %>%
  ungroup()


prim_t6_umap_data <- dplyr::filter(shorter_big_table_t6, n_infection=="First")

sig_prim_colour_t_cells_umap <- ggplot(prim_t6_umap_data, aes(x=UMAP2, y=UMAP1, color=prim_significant, alpha=prim_alpha))+
  geom_point(shape=".")+
  scale_colour_manual(values = cols)+
  UMAP_theme+
  ggtitle("significantly up at T6\nfirst infection (18)")+
  theme(
    #panel.spacing.x = unit(10, "mm"),
    plot.title = element_text(size=13, hjust=0.5),
    axis.title.y = element_blank(),
    axis.title = element_text(size=10))
cowplot::ggsave2("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/figures_for_paper/sig_prim_colour_t_cells_umap.pdf", sig_prim_colour_t_cells_umap, width=3, height=3)


ter_t6_umap_data <- dplyr::filter(shorter_big_table_t6, n_infection=="Third")

sig_ter_colour_t_cells_umap <- ggplot(ter_t6_umap_data, aes(x=UMAP2, y=UMAP1, color=ter_significant, alpha=ter_alpha))+
  geom_point(shape=".")+
  scale_colour_manual(values = cols)+
  ggtitle("significantly up at T6\nthird infection (12)")+
  UMAP_theme+
  theme(
    #legend.position="right",
    #panel.spacing.x = unit(10, "mm"),
    plot.title = element_text(size=13, hjust=0.5),
    axis.title = element_blank())
cowplot::ggsave2("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/figures_for_paper/sig_ter_colour_t_cells_umap.pdf", sig_ter_colour_t_cells_umap, width=3, height=3)



all_colour_t_cells_umap <- ggplot(shorter_big_table_t6, aes(x=UMAP2, y=UMAP1, color=flo_label))+
  geom_point(shape=".")+
  scale_colour_manual(values = cols)+
  UMAP_theme+
  ggtitle("All Clusters at T6\n")+
  theme(
    plot.title = element_text(size=13, hjust=0.5),
    axis.title = element_text(size=10),
    axis.title.x = element_blank())

cowplot::ggsave2("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/figures_for_paper/all_colour_t_cells_umap.pdf", all_colour_t_cells_umap, width=3, height=3)


# combo_plot <- plot_grid(all_colour_t_cells_umap, sig_prim_colour_t_cells_umap, sig_ter_colour_t_cells_umap, nrow=1, align = "hv")
# ggsave("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/combo_plot.png", combo_plot, width=9, height=3)


# phenotypic heatmaps ####

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


scaled_mat <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_medians_prim_ter_T6.csv", row.names = 1)
reordered_scaled_mat <- scaled_mat[,match(refined_markers, colnames(scaled_mat))]
rownames(reordered_scaled_mat) <- gsub("activated RA-RO- DN", "activated DN", rownames(reordered_scaled_mat))

ter_t6_df <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/edgeR/ter_t6_df_edger.csv")

ter_sig_clusters <- ter_t6_df %>%
  filter(p_adj <= 0.05, logFC >=1)%>%
  select(cluster_id)


#subset heatmpa to be only all sig clusters
prim_sig_clusters <- read.delim("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/prim_ter_sig_clusters.csv", header=FALSE)$V1
reordered_scaled_mat <- subset(reordered_scaled_mat, rownames(reordered_scaled_mat) %in% prim_sig_clusters)



colcsv <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_colours.csv")
cluster_order <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_order.csv")$x

sig_colcsv <- subset(colcsv, colcsv$cluster_id %in% prim_sig_clusters)


col_pal <- sig_colcsv$colour
names(col_pal) <- sig_colcsv$cluster_id


both <- ifelse(rownames(reordered_scaled_mat) %in% ter_sig_clusters$cluster_id, "both", "first only")


combo_left_anno_var <-  rowAnnotation(annotation_name_gp = gpar(fontsize=10),
                                      annotation_name_rot = 45,
                                      gap = unit(1.5, "mm"),
                                      "cluster_id"=sig_colcsv$cluster_id,
                                      "n_infection"=both,
                                      show_legend = c(FALSE, TRUE),
                                      show_annotation_name = TRUE,
                                      simple_anno_size = unit(8, "mm"), # width of the significance bar
                                      col=list("n_infection" = c("both"="darkgrey", "first only"="#36454f"),
                                               "cluster_id" = col_pal),
                                      annotation_legend_param = list(n_infection = list(title ="n_infection",
                                                                                        at = c("both", "first only"),
                                                                                        title_gp=gpar(angle=45),
                                                                                        legend_gp = gpar(fill = c("darkgrey","#36454f")),
                                                                                        title_position = "topleft")
                                                                     # cluster_id = list(title="cluster_id",
                                                                     #                   at = sig_colcsv$cluster_id,
                                                                     #                   legend_gp = gpar(fill = unname(col_pal)))
                                      )
                                      
)





inferno <- colorspace::sequential_hcl("inferno", n=8)
col_inferno <- circlize::colorRamp2(seq(0,1, by=1/{length(inferno)-1}), inferno)

reordered_scaled_mat <- subset(reordered_scaled_mat, rownames(reordered_scaled_mat) %in% prim_sig_clusters)

rereordered_scaled_mat <- reordered_scaled_mat[rev(cluster_order),]


# ht_opt("ROW_ANNO_PADDING" = unit(3,"mm"))
# ht_opt(RESET = TRUE)

sig_cluster_heatmap <- Heatmap(matrix = as.matrix(rereordered_scaled_mat),
                               #cluster_rows = FALSE,
                               cluster_rows = TRUE,
                               show_row_dend = TRUE,
                               show_heatmap_legend = TRUE,
                               cluster_columns = FALSE,
                               row_names_side = "left",
                               row_dend_side = "right",
                               col = col_inferno,
                               rect_gp = gpar(col = "white"),
                               left_annotation = combo_left_anno_var,
                               column_names_rot = 45,
                               heatmap_legend_param = list(#legend_position = "top",
                                                           col=col_inferno,
                                                           title = "Normalised Marker Expression",
                                                           #legend_direction = "horizontal",
                                                           title_position = "leftcenter-rot",
                                                           legend_height = unit(6.2, "cm"),
                                                           border = FALSE)
)

pdf("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/figures_for_paper/sig_cluster_heatmap_var.pdf", width=16, height=8)
draw(sig_cluster_heatmap, padding = unit(c(2, 25, 2, 15), "mm"))
dev.off()


big_scaled_mat <- scaled_mat[,match(refined_markers, colnames(scaled_mat))]
rownames(big_scaled_mat) <- gsub("activated RA-RO- DN", "activated DN",rownames(big_scaled_mat))

all_cluster_heatmap <- Heatmap(matrix = as.matrix(big_scaled_mat),
                               #cluster_rows = FALSE,
                               cluster_rows = TRUE,
                               show_row_dend = TRUE,
                               show_heatmap_legend = TRUE,
                               cluster_columns = FALSE,
                               row_names_side = "left",
                               row_dend_side = "right",
                               col = col_inferno,
                               rect_gp = gpar(col = "white"),
                               column_names_rot = 45,
                               heatmap_legend_param = list(#legend_position = "top",
                                 col=col_inferno,
                                 title = "Normalised Marker Expression",
                                 #legend_direction = "horizontal",
                                 title_position = "leftcenter-rot",
                                 legend_height = unit(6.2, "cm"),
                                 border = FALSE)
)

pdf("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/figures_for_paper/supp_all_cluster_heatmap.pdf", width=16, height=18)
draw(all_cluster_heatmap, padding = unit(c(2, 25, 2, 15), "mm"))
dev.off()