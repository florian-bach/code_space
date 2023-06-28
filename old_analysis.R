library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(MASS)

set.seed(12345)

inferno <- colorspace::sequential_hcl("inferno", n=8)

volunteer_colors <- c("#AA3377", "#EE6677", "pink", "#4477AA", "#66CCEE", "#228833")
names(volunteer_colors) <-  rev(c("v305", "v308", "v310", "v320", "v313", "v315"))


# harmony with cluster 12 removed ####

big_data <- SeuratDisk::LoadH5Seurat("~/postdoc/edinburgh/scRNAseq/seurat_objecitons/final_activated_subset.h5Seurat")

big_data_plot <- DimPlot(big_data, label=TRUE)
ggsave("~/postdoc/edinburgh/scRNAseq/oldno_harmony_pbmc_dot.png", big_data_plot, width=8, height=8, bg="white")

acti_comparers <- c("CD38", "HLA-DRA", "TBX21", "FOXP3", "IKZF2", "TCF7", "CCL5", "KLRB1", "GZMK", "GZMB", "GNLY", "STAT2")
pbmc_list <- c("LEF1", "TCF7", "CCR7", "IL7R", "ITGB1", "GATA3", "LGALS1", "AQP3", "COTL1", "CCL5", "GZMK", "GZMB", "NKG7", "CD8A", "CD4", "CD320", "FOXP3", "TRAV1-2", "MKI67", "TRDC")
til_list <- c("CD3E", "CD8A", "LEF1", "TCF7", "CCR7", "IL7R", "AQP3", "ZNF683", "CX3CR1", "CCL5", "GZMK", "GZMB", "IFNG", "LAG3", "FOXP3", "TRAV1-2", "TRDC", "S100A8", "CD79A", "KRT19")

p1 <- DotPlot(big_data, features = pbmc_list) + RotatedAxis()
p2 <- DotPlot(big_data, features = til_list) + RotatedAxis()

ggplot2::ggsave("~/postdoc/edinburgh/scRNAseq/oldno_harmony_pbmc_dot.png", p1, width=8, height=8, bg="white")
ggplot2::ggsave("~/postdoc/edinburgh/scRNAseq/oldno_harmony_til_dot.png", p2, width=8, height=8, bg="white")


big_comp <- FeaturePlot(big_data, features=acti_comparers, order=TRUE, split.by = "timepoint")&scale_color_gradientn(colors=inferno)
big_comp2 <- FeaturePlot(big_data, features=acti_comparers, order=TRUE, split.by = "time_n_infection")&scale_color_gradientn(colors=inferno)

ggsave("~/postdoc/edinburgh/scRNAseq/oldbig_big_data_harmony_coomparison.png", big_comp, width=8, height = 4*length(acti_comparers), dpi=444)
ggsave("~/postdoc/edinburgh/scRNAseq/oldbig_big_data_harmony_coomparison2.png", big_comp2, width=16, height = 4*length(acti_comparers), dpi=444)


cluster_counts <- data.frame(table(big_data@active.ident, big_data@meta.data$sample_ID))
colnames(cluster_counts) <- c("Cluster_ID", "Sample_ID", "Count")
cluster_counts$Timepoint <- substr(cluster_counts$Sample_ID, 6, 14)
cluster_counts$Volunteer <- substr(cluster_counts$Sample_ID, 1, 4)
cluster_counts$N_Infection <- ifelse(cluster_counts$Volunteer %in% c("v313", "v315", "v320"), "First", "Third")


cluster_counts <- cluster_counts %>%
  group_by(Sample_ID) %>%
  mutate("Percentage"=Count/sum(Count))


cluster_percentages <- ggplot(cluster_counts, aes(x=Cluster_ID, y=Percentage, fill=Timepoint))+
  geom_boxplot()+
  geom_point(position = position_dodge(width = 0.75), aes(color=Volunteer, group=Timepoint))+
  scale_y_continuous(labels = scales::label_percent())+
  scale_fill_manual(values=c("#AA3377", "#4477AA"))+
  #scale_color_manual(values = volunteer_colors)+
  facet_wrap(~N_Infection, ncol=1)+
  theme_minimal()+
  theme(strip.text = element_text(hjust = 0))

ggsave("~/postdoc/edinburgh/scRNAseq/oldbig_big_data_harmony_cluster_percentages.png", width=8, height=5, bg="white")
