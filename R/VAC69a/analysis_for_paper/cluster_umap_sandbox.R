cluster_palette <- color_103_scheme[2:(length(unique(big_table$flo_label))+1)]
cluster_palette[8] <- "#FFA500"
cluster_palette[9] <- "red"
cluster_palette[1] <- "#000000"
names(cluster_palette)[1] <- "black"
names(cluster_palette)[2:length(unique(big_table$significant))] <- unique(big_table$significant)[-match("black", unique(big_table$significant))]



# restrict data to T6 and downsample by 33% to make it look nice
short_big_table_t6 <- subset(big_table, big_table$timepoint=="T6")
short_big_table_t6 <- short_big_table_t6[seq(1,nrow(short_big_table_t6), by=3), ]

(t6_all_clusters_umap <- ggplot(short_big_table_t6, aes(x=UMAP1, y=UMAP2))+
  #geom_point(aes(color=flo_label, alpha=scales::rescale(alpha, to=c(0.1, 1))), shape=".")+
  geom_point(aes(color=flo_label), shape=".", alpha=0.4)+
    #UMAP_theme+
  scale_color_manual(values = expanded_cluster_palette)+
  theme_minimal()+
  ggtitle("Differentially Abundant\nClusters at T6")+
  UMAP_theme+
  guides(colour = guide_legend(override.aes = list(size = 4, shape=16, alpha=1)),
         alpha= "none")+
  # scale_y_continuous(limits = c(-11.2, 11.3))+
  # scale_x_continuous(limits=c(-13, 10))+
  #guides(alpha = guide_legend(override.aes = list(size = 10)))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust=0.5))+
  coord_cartesian(xlim=c(-13, 10),
                  ylim=c(-11.2, 11.3))
)




