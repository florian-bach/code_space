# open diffuerential_abundance.R and run until diffcyt constructs edgeR comparisons for each timepoint
# then proceed

# all_cluster_freqs <- diffcyt_boxplot(da_t6, merged_daf, counts=T, FDR=1, logFC = 0)
# 
# 
# 
# stacked_bar_data <-  all_cluster_freqs$data
# 
# #write.csv(t6_map_data, "/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/barchart_data.csv")
# 
# t6_map_data <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/barchart_data.csv")
#
# activated_clusters <- subset(merging_table1$new_cluster, grepl("activated", merging_table1$new_cluster))
# activated_clusters <- unique(activated_clusters)
#
# #write.csv(sig_t6_clusters, "/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/sig_t6_clusters.csv")
#
# activated_barchart_data <- subset(stacked_bar_data, stacked_bar_data$cluster_id %in% activated_clusters)
#
# # get lineage gates
# lineage_merge <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/merging_tables/lineage_merge.csv", header=T, stringsAsFactors = F)
# merged_daf<- mergeClusters(merged_daf, k = "meta45", table = lineage_merge, id = "lineage")
#
# cluster_dic <- cbind(merging_table1, "lineage"=lineage_merge$new_cluster)
#
# #categorise significant clusters into their respective lineages
# activated_barchart_data$lineage <- cluster_dic$lineage[match(activated_barchart_data$cluster_id, cluster_dic$new_cluster)]
#
#
#
# summary <- activated_barchart_data %>%
#   group_by(timepoint, cluster_id) %>%
#   mutate(mean=base::mean(frequency), sd=stats::sd(frequency))
#
# tail(summary)
#
# summary$lineage <- factor(summary$lineage, levels=c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "Resting"))
#

#taken from khroma bright 6

#data.table::fwrite(summary,  "~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/activation_barchart_data")

library(ggplot2)
library(dplyr)
library(cowplot)

colcsv <- read.csv("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cluster_palette.csv", header=T, stringsAsFactors = F)

col_pal <- colcsv$x
names(col_pal) <- colcsv$X


# Stacked Barchart Figure 1E ####

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
# 
# 
# 
# 
lineage_merge <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/merging_tables/lineage_merge.csv", header=T, stringsAsFactors = F)

merged_daf<- mergeClusters(merged_daf, k = "meta45", table = lineage_merge, id = "lineage")

cluster_dic <- cbind(merging_table1, "lineage"=lineage_merge$new_cluster)


barchart_data <- all_cluster_freqs$data

barchart_data$lineage <- cluster_dic$lineage[match(barchart_data$cluster_id, cluster_dic$new_cluster)]

barchart_data <- barchart_data %>%
  group_by(timepoint, lineage) %>%
  mutate("lin_max_count" = sum(count)) %>%
  ungroup()

barchart_data <- barchart_data %>%
  group_by(timepoint, cluster_id) %>%
  mutate("cluster_max_count" = sum(count)) %>%
  ungroup()

barchart_data$lin_freq <- barchart_data$cluster_max_count/barchart_data$lin_max_count

# 
write.csv(barchart_data, "/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/cluster_and_lineage_stats_for_barcharts.csv")


# PIES Figure 1F ####


barchart_data <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/cluster_and_lineage_stats_for_barcharts.csv", header=T, stringsAsFactors = FALSE)



activated_clusters <- grep("activated", unique(barchart_data$cluster_id), value = T)
barchart_data <- subset(barchart_data, barchart_data$cluster_id %in% activated_clusters)


barchart_data <- subset(barchart_data, barchart_data$timepoint=="T6")


barchart_data <- barchart_data %>%
  group_by(timepoint, lineage) %>%
  mutate("lin_activation" = sum(lin_freq)/6) %>%
  ungroup()




barchart_plot_data <- data.frame("lineage" = barchart_data$lineage, "lin_activation" = barchart_data$lin_activation)

barchart_plot_data <- barchart_plot_data[!duplicated(barchart_plot_data), ]

barchart_plot_data <- rbind(barchart_plot_data, data.frame("lineage"=barchart_plot_data$lineage, "lin_activation"=1-barchart_plot_data$lin_activation))

barchart_plot_data$state <- c(as.character(barchart_plot_data$lineage[1:6]), rep("Resting", 6))

barchart_plot_data$state <- gsub("gd", expression(paste(gamma, delta)), barchart_plot_data$state)
barchart_plot_data$lineage <- gsub("gd", expression(paste(gamma, delta)), barchart_plot_data$lineage)


# we have to rename the colour palette because otherwise the strip.text won't parse properly
pie_lineage_palette <- lineage_palette
names(pie_lineage_palette)<- c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "Resting")

barchart_plot_data$lineage <- factor(barchart_plot_data$lineage, levels=c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "Resting"))
barchart_plot_data$state <- factor(barchart_plot_data$state, levels=rev(c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "Resting")))

lineage_activation_pies <- ggplot(barchart_plot_data)+
  #geom_bar(stat="identity", fill="white", x="", y=1, color="black")+
  geom_bar(stat="identity", position = "stack", width=1, aes(x="", y=lin_activation, fill=state))+
  theme_minimal()+
  coord_polar(theta = "y", start = 0)+
  facet_wrap(~lineage, labeller = label_parsed, ncol=2, strip.position = "left")+
  scale_fill_manual(values=pie_lineage_palette)+
  # scale_fill_manual(values=lineage_palette, breaks = levels(barchart_plot_data$lineage))+
  guides()+
  ggtitle("Average Proportion of Activated Cells\nin Major T cell Lineages at T6")+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size=11),
        strip.text.y.left = element_text(hjust=0.5, size=14, face = "bold", angle = 0),
        legend.position = "none")

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/lineage_activation_pies.png", lineage_activation_pies, height=3, width=4.5)



panel_f <- plot_grid(activation_stacked_barchart, lineage_activation_pies, ncol = 2, rel_widths = c(4,1), rel_heights = c(1,1), axis = "bt", align = "hv")

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/panel_f.png", panel_f, height=4, width=17)




# FIGURE 2B ####

library(gridBase)
library(grid)
library(circlize)
library(ComplexHeatmap)
library(dplyr)

inferno <- colorspace::sequential_hcl("inferno", n=8)
col_inferno <-colorRamp2(seq(0,1, by=1/{length(inferno)-1}), inferno)

reordered_sig_scaled_mat <- as.matrix(read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/cluster_medians_heatmap_t6.csv", header = T, row.names = 1))

sig_clusters <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/sig_t6_clusters.csv", header = TRUE, stringsAsFactors = FALSE)
sig_clusters <- sig_clusters <- sig_clusters[-9,2]

all_t6_data <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_t6_data.csv", stringsAsFactors = FALSE, header = T)

t6_map_data <- all_t6_data%>%
  group_by(cluster_id, timepoint) %>%
  mutate(mean_freq=mean(frequency)) %>%
  ungroup() %>%
  select(cluster_id, mean_freq, timepoint)

t6_map_data <- filter(t6_map_data, timepoint=="T6")
t6_map_data <- t6_map_data[!duplicated(t6_map_data), ]


t6_map_data$cluster_id <- ifelse(substr(t6_map_data$cluster_id, 1, 1)==" ", substr(t6_map_data$cluster_id, 2, nchar(t6_map_data$cluster_id)), t6_map_data$cluster_id)
t6_map_data$cluster_id <- ifelse(substr(t6_map_data$cluster_id, 1, 1)==" ", substr(t6_map_data$cluster_id, 2, nchar(t6_map_data$cluster_id)), t6_map_data$cluster_id)

t6_map_data$cluster_id <- ifelse(substr(t6_map_data$cluster_id, nchar(t6_map_data$cluster_id), nchar(t6_map_data$cluster_id))==" ", substr(t6_map_data$cluster_id, 1, nchar(t6_map_data$cluster_id)-1), t6_map_data$cluster_id)
t6_map_data$cluster_id <- ifelse(substr(t6_map_data$cluster_id, nchar(t6_map_data$cluster_id), nchar(t6_map_data$cluster_id))==" ", substr(t6_map_data$cluster_id, 1, nchar(t6_map_data$cluster_id)-1), t6_map_data$cluster_id)



t6_map_data <- filter(t6_map_data, cluster_id %in% sig_clusters)
t6_map_data <- t6_map_data[order(t6_map_data$mean_freq, decreasing = T),]


pie_data <- data.frame(t6_map_data)


pie_data$lineage <- factor(c("gd", "CD4", "CD4", "MAIT", "CD4", "CD4", "CD4", "CD8", "CD8"), levels=c("CD4", "CD8", "MAIT", "gd"))
pie_data <- pie_data[order(pie_data$lineage),]

#pie_data$lineage <- gsub("gd", expression(paste(gamma, delta)), pie_data$lineage)

colcsv <- read.csv("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cluster_palette.csv", header=T, stringsAsFactors = F)

col_pal <- colcsv$x
names(col_pal) <- colcsv$X

ordered_col_pal <- col_pal[match(pie_data$cluster_id, names(col_pal))]

circlize_plot = function() {
  
  circos.par("track.height" = 0.6, start.degree = 90)
  circos.initialize(factors = pie_data$cluster_id, xlim=c(0,1), sector.width = pie_data$mean_freq)
  
  o.cell.padding = circos.par("cell.padding")
  
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.2, bg.border = NA, 
                         cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4]))
  
  circos.track(factors = pie_data$cluster_id,
               ylim = c(0, 1),
               x=pie_data$mean_freq,
               bg.col = ordered_col_pal)
  
  
  #no color, just text
  
  highlight.sector(sector.index = grep("CD4", pie_data$cluster_id, value = T), track.index = 1, 
                   border = "black", text = "CD4", col = "white")
  highlight.sector(sector.index = grep("CD8", pie_data$cluster_id, value = T), track.index = 1, 
                   border = "black", text = "CD8", col = "white")
  highlight.sector(sector.index = grep("MAIT", pie_data$cluster_id, value = T), track.index = 1, 
                   border = "black", text = "MAIT", col = "white")
  highlight.sector(sector.index = grep("Vd2", pie_data$cluster_id, value = T), track.index = 1, 
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

lvl_data <- pie_data
  # lvl_data$lineage <- factor(lvl_data$lineage, levels=c(expression(paste(gamma, delta)), "MAIT", "CD4", "CD8"))
  lvl_data$lineage <- factor(lvl_data$lineage, levels = c("CD4", "CD8", "MAIT", "gd"))
  
  # lvls <- lvl_data$cluster_id[c(8, 9, seq(1,7))]
  lvls <- lvl_data$cluster_id[c(seq(1,7), 9, 8 )]
  
  
  # discrete
  lgd_cluster= Legend(at =lvls, type = "grid", 
                      legend_gp = gpar(fill = unname(ordered_col_pal)[c(seq(1,7), 9, 8 )]),
                      title_position = "topleft",
                      labels_gp = gpar(fontsize = 10),
                      title = "Cluster ID")
  # discrete
  # lgd_lineage = Legend(at = c("CD4", "CD8", "MAIT", "gd"), type = "grid", 
  #                      legend_gp = gpar(fill = unname(pie_palette)[1:4]),
  #                      title_position = "topleft",
  #                      labels_gp = gpar(fontsize = 10, parse = TRUE),
  #                      title = "Lineage", )
  
  #lgd_list <- packLegend(lgd_cluster, lgd_lineage)
  
  
  circlize_plot()
  
  draw(lgd_list)
  
  lgd_lineage@grob$children
  
    
    
      
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
  


# figure 2 D memory activation barchart ####


summary <- data.table::fread("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/activation_barchart_data")

cd4_summary <- subset(summary, grepl("CD4", summary$cluster_id))

lineage_palette <- c("#AA3377", "#EE6677", "#4477AA", "#66CCEE", "#228833", "#FFA500", "#BBBBBB")
names(lineage_palette) <-c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "Resting")

cd4_summary <- cd4_summary[cd4_summary$timepoint=="T6",]




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


#stacked_bar_levels <- levels(factor(cd4_bar_data$cluster_id))[c(3,1,4,2,5)]
stacked_bar_levels <- lvls[1:5] # defined on line 287

cd4_memory_activation_stacked_barchart <- ggplot(cd4_bar_data, aes(x=volunteer, y=cd4_freq, fill=factor(cluster_id, levels = stacked_bar_levels)))+
  geom_bar(stat="identity", position="stack")+
  #facet_wrap(~volunteer, strip.position = "bottom", ncol=3)+
  ggtitle("Activation of CD4 Memory Cells at T6\n")+
  scale_fill_manual(values=col_pal)+
  scale_y_continuous(name = "Percentage of Activated\nCD4 memory T cells", labels=scales::percent_format(accuracy = 1))+
  theme_minimal()+
  xlab("Volunteer")+
  guides(fill=guide_legend(ncol=2))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        legend.title = element_blank(),
    strip.text = element_text(hjust=0.5, size=10, face = "bold"),
    #axis.title.x = element_blank(),
    #axis.text.x = element_blank(),
    legend.position = "bottom",
    legend.justification = "center",
    panel.grid.minor.y = element_blank(),
    strip.placement = "outside")


ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cd4_memory_activation_stacked_barchart.png", cd4_memory_activation_stacked_barchart, width=6, height = 6)

      ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cd4_memory_activation_stacked_barchart_no_leg.png", cd4_memory_activation_stacked_barchart, width=5, height = 5)

  # old code you probably won't need again; this sums the activated counts and divides by the sums of lineages, rather thandisaplying averages


# pie_data <- barchart_data
# 
# t6_pie_data <- subset(pie_data, pie_data$timepoint=="T6")
# 
# t6_lin_counts <- t6_pie_data %>%
#   group_by(volunteer, lineage) %>%
#   mutate("sum_of_lineage"=sum(count)) %>%
#   select(volunteer, lineage, sum_of_lineage) %>%
#   ungroup()
# 
# t6_lin_counts <- t6_lin_counts[!duplicated(t6_lin_counts), ]
# 
# 
# 
# 
# t6_act_counts <- subset(t6_pie_data, t6_pie_data$cluster_id %in% activated_clusters)
# 
# t6_act_counts <- t6_act_counts %>%
#   group_by(volunteer, lineage) %>%
#   mutate("sum_of_activated"=sum(count)) %>%
#   select(volunteer, lineage, sum_of_activated) %>%
#   ungroup()
# 
# t6_act_counts <- t6_act_counts[!duplicated(t6_act_counts), ]
# 
# t6_act_counts <- t6_act_counts[order(t6_act_counts$volunteer, t6_act_counts$lineage),]
# t6_lin_counts <- t6_lin_counts[order(t6_lin_counts$volunteer, t6_lin_counts$lineage),]
# 
# pie_data <- cbind(t6_lin_counts, t6_act_counts)  
# pie_data[,1:2] <- NULL
# pie_data$ratio <- pie_data$sum_of_activated/pie_data$sum_of_lineage
# 
# final_pie <- pie_data %>%
#   group_by(lineage) %>%
#   mutate("mean_lin_act"= mean(ratio))
# 
# pie_final_pie <- select(final_pie, lineage, mean_lin_act)
# pie_final_pie$lineage <- gsub("gd", expression(paste(gamma, delta)), pie_final_pie$lineage)
# 
# 
# lineage_activation_pies <- ggplot(pie_final_pie)+
#   geom_rect(aes(xmin=0, xmax=1, ymin=0, ymax=mean_lin_act, fill=lineage))+
#   geom_rect(aes(xmin=0, xmax=1, ymin=mean_lin_act, ymax=1), fill="grey")+
#   facet_wrap(~lineage, labeller = label_parsed, ncol=2)+
#   coord_polar(theta = "y", start = 0)+
#   theme_minimal()+
#   scale_fill_manual(values=pie_lineage_palette)+
#   ggtitle("Average Proportion of Activation by Lineage, \nacross all samples at T6")+
#   theme(axis.title = element_blank(),
#         axis.text = element_blank(),
#         panel.grid = element_blank(),
#         plot.title = element_text(hjust=0.5, size=15),
#         strip.text = element_text(hjust=0.5, size=18, face = "bold"),
#         legend.position = "none")
# 
# 
# ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/lineage_activation_pies.png", lineage_activation_pies, height=3, width=4.5)



