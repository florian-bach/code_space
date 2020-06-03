# open diffuerential_abundance.R and run until diffcyt constructs edgeR comparisons for each timepoint
# then proceed

all_cluster_freqs <- diffcyt_boxplot(da_dod, merged_daf, counts=F, FDR=1, logFC = 0)


stacked_bar_data <-  all_cluster_freqs$data

#write.csv(t6_map_data, "/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/barchart_data.csv")

t6_map_data <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/barchart_data.csv")

activated_clusters <- subset(merging_table1$new_cluster, grepl("activated", merging_table1$new_cluster))
activated_clusters <- unique(activated_clusters)

#write.csv(sig_t6_clusters, "/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/sig_t6_clusters.csv")

activated_barchart_data <- subset(stacked_bar_data, stacked_bar_data$cluster_id %in% activated_clusters)

# get lineage gates
lineage_merge <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/merging_tables/lineage_merge.csv", header=T, stringsAsFactors = F)
merged_daf<- mergeClusters(merged_daf, k = "meta45", table = lineage_merge, id = "lineage")

cluster_dic <- cbind(merging_table1, "lineage"=lineage_merge$new_cluster)

#categorise significant clusters into their respective lineages
activated_barchart_data$lineage <- cluster_dic$lineage[match(activated_barchart_data$cluster_id, cluster_dic$new_cluster)]



summary <- activated_barchart_data %>%
  group_by(timepoint, cluster_id) %>%
  mutate(mean=base::mean(frequency), sd=stats::sd(frequency))

tail(summary)

summary$lineage <- factor(summary$lineage, levels=c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "Resting"))


#taken from khroma bright 6
lineage_palette <- c("#AA3377", "#EE6677", "#4477AA", "#66CCEE", "#228833", "#CCBB44", "#BBBBBB")
names(lineage_palette) <-c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "Resting")

activation_stacked_barchart <- ggplot(summary, aes(x=volunteer, y=frequency/100, fill=lineage))+
  geom_bar(stat="identity", position="stack")+
  #geom_text(aes(y=(total_cd3/100)+0.01, label=paste0(round(total_cd3, digits = 1), "%", sep='')))+
  theme_minimal()+
  facet_wrap(~timepoint, strip.position = "bottom")+
  ggtitle("Overall T cell activation")+
  scale_fill_manual(values=lineage_palette, labels=c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "Resting"))+
  scale_y_continuous(name = "Percentage of CD3+ T cells\n\n", labels=scales::percent_format(accuracy = 1))+
  #ylim(0,25)+
  #geom_text(aes(label=cluster_id), position = position_stack(vjust = .5))+
  theme(#legend.position = "none",
    plot.title = element_text(hjust=0.5, size=15),
    strip.text = element_text(hjust=0.5, size=12, face = "bold"),
    axis.title.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none",
    strip.placement = "outside")


ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/activation_stacked_barchart.png", activation_stacked_barchart, height=4, width=4)




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


write.csv(barchart_data, "/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/cluster_and_lineage_stats_for_barcharts.csv")

barchart_data <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/cluster_and_lineage_stats_for_barcharts.csv", header=T, stringsAsFactors = FALSE)


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
  facet_wrap(~lineage, labeller = label_parsed, ncol=2)+
  scale_fill_manual(values=pie_lineage_palette)+
  # scale_fill_manual(values=lineage_palette, breaks = levels(barchart_plot_data$lineage))+
  guides()+
  ggtitle("Average Proportion of Activation by Lineage,\nacross all samples at T6")+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5, size=15),
        strip.text = element_text(hjust=0.5, size=18, face = "bold"),
        legend.position = "none")

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/lineage_activation_pies.png", lineage_activation_pies, height=3, width=4.5)



panel_f <- plot_grid(activation_stacked_barchart, lineage_activation_pies, ncol = 2, rel_widths = c(1.8,1))

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/panel_f.png", panel_f, height=6, width=9)




pie_data <- barchart_data

t6_pie_data <- subset(pie_data, pie_data$timepoint=="T6")

t6_lin_counts <- t6_pie_data %>%
  group_by(volunteer, lineage) %>%
  mutate("sum_of_lineage"=sum(count)) %>%
  select(volunteer, lineage, sum_of_lineage) %>%
  ungroup()

t6_lin_counts <- t6_lin_counts[!duplicated(t6_lin_counts), ]
  
  
  
  
t6_act_counts <- subset(t6_pie_data, t6_pie_data$cluster_id %in% activated_clusters)

t6_act_counts <- t6_act_counts %>%
  group_by(volunteer, lineage) %>%
  mutate("sum_of_activated"=sum(count)) %>%
  select(volunteer, lineage, sum_of_activated) %>%
  ungroup()
  
t6_act_counts <- t6_act_counts[!duplicated(t6_act_counts), ]

t6_act_counts <- t6_act_counts[order(t6_act_counts$volunteer, t6_act_counts$lineage),]
t6_lin_counts <- t6_lin_counts[order(t6_lin_counts$volunteer, t6_lin_counts$lineage),]

pie_data <- cbind(t6_lin_counts, t6_act_counts)  
pie_data[,1:2] <- NULL
pie_data$ratio <- pie_data$sum_of_activated/pie_data$sum_of_lineage

final_pie <- pie_data %>%
  group_by(lineage) %>%
  mutate("mean_lin_act"= mean(ratio))

pie_final_pie <- select(final_pie, lineage, mean_lin_act)
pie_final_pie$lineage <- gsub("gd", expression(paste(gamma, delta)), pie_final_pie$lineage)


lineage_activation_pies <- ggplot(pie_final_pie)+
  geom_rect(aes(xmin=0, xmax=1, ymin=0, ymax=mean_lin_act, fill=lineage))+
  geom_rect(aes(xmin=0, xmax=1, ymin=mean_lin_act, ymax=1), fill="grey")+
  facet_wrap(~lineage, labeller = label_parsed, ncol=2)+
  coord_polar(theta = "y", start = 0)+
  theme_minimal()+
  scale_fill_manual(values=pie_lineage_palette)+
  ggtitle("Average Proportion of Activation by Lineage, \nacross all samples at T6")+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5, size=15),
        strip.text = element_text(hjust=0.5, size=18, face = "bold"),
        legend.position = "none")
  
  
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/lineage_activation_pies.png", lineage_activation_pies, height=3, width=4.5)



ggplot(pie_data, aes(x="", y=ratio, fill=lineage))+
  geom_point(aes(shape=volunteer))+
  geom_boxplot()+
  facet_wrap(~lineage, scales="free")+
  theme_minimal()
