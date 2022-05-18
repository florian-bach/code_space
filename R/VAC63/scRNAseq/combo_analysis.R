library(Seurat)
# 
first <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/demultiplexed_primary.h5Seurat")
third <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/demultiplexed_tertiary.h5Seurat")
# 
# first10k <- first[, sample(colnames(first), size =10000, replace=F)]
# 
#   SeuratDisk::SaveH5Seurat(first10k, "~/postdoc/scRNAseq/Seurat_Objects/10k_demultiplexed_primary.h5Seurat")
# third10k <- third[, sample(colnames(third), size =10000, replace=F)]
#   SeuratDisk::SaveH5Seurat(third10k, "~/postdoc/scRNAseq/Seurat_Objects/10k_demultiplexed_tertiary.h5Seurat")

# first <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/10k_demultiplexed_primary.h5Seurat")
# third <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/10k_demultiplexed_tertiary.h5Seurat")
#   

combo <- merge(first, third, add.cell.ids = c("First", "Third"))

combo[["percent.mt"]] <- PercentageFeatureSet(combo, pattern = "^MT-")

#filter this bad boy
combo <- subset(combo, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# pre-processing: transformations & feature selection ####

combo <- NormalizeData(combo, normalization.method = "LogNormalize", scale.factor = 10000)

#find variable genes across cells
combo <- FindVariableFeatures(combo, selection.method = "vst", nfeatures = 1500)
#top10 <- head(VariableFeatures(combo), 30)

# plot variable features with and without labels
# plot1 <- VariableFeaturePlot(combo)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# scaling data
all.genes <- rownames(combo)
combo <- ScaleData(combo)

# dimensionality reduction ####
#  Examine and visualize PCA results a few different ways
combo <- RunPCA(combo, features = VariableFeatures(object = combo))

# clustering ####


combo <- FindNeighbors(combo, dims = 1:20)
combo <- FindClusters(combo, resolution = 1.5, random.seed = 1234, algorithm = 1)

# find markers distinguishing clusters
# combo.markers <- FindAllMarkers(combo, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)



combo <- AddMetaData(
  object = combo,
  metadata = substr(colnames(combo), 1,5),
  col.name = "n_infection"
)

combo <- AddMetaData(
  object = combo,
  metadata = paste(combo$timepoint, combo$n_infection, sep="_"),
  col.name = "time_n_infection"
)

combo$time_n_infection <- factor(combo$time_n_infection, levels = c("Baseline_First", "T6_First", "Baseline_Third", "T6_Third"))

#dimensionality reduction ####
combo <- RunUMAP(combo, dims = 1:30, seed.use = 1234, n.neighbors = 50); DimPlot(combo, reduction = "umap", label=TRUE)


cluster_umap <- DimPlot(combo, reduction = "umap", label=TRUE)

ggplot2::ggsave("~/postdoc/scRNAseq/figures_for_Ange/combo_cluster_umap.png", cluster_umap, height=4, width=6, bg="white")



combo_matrix <- data.frame(t(combo@assays$RNA@scale.data))

# randomly subsample 5555 cells to speed up plotting
# set.seed(1234); small_combo_matrix <- combo_matrix[floor(runif(9999, min=1, max=nrow(combo_matrix))),]


phate_object <- phateR::phate(combo_matrix, knn = 45, decay = 70, t=7)
#small_pbmc[["phate"]] <- CreateDimReducObject(embeddings = phate_object$embedding, key = "phate_", assay = DefaultAssay(small_pbmc))

library(ggplot2)

feature1 <- "SELL"
feature2 <- "CCR7"
feature3 <- "IL7R"
feature4 <- "FOXP3"


feature1 <- "CD38"
feature2 <- "HLA.DRB1"
feature3 <- "GZMH"
feature4 <- "CCL5"

(plot1 <- ggplot(phate_object)+
    geom_point(aes(PHATE1, PHATE2, color=combo_matrix[,feature1]), alpha=0.5, size=0.4)+
    labs(color=feature1)+ 
    ggtitle(feature1)+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_blank())+
    viridis::scale_color_viridis(option="B"))

plot2 <- ggplot(phate_object)+
  geom_point(aes(PHATE1, PHATE2, color=combo_matrix[,feature2]), alpha=0.5, size=0.4)+
  ggtitle(feature2)+
  labs(color=feature2)+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_blank())+
  viridis::scale_color_viridis(option="B")


(plot3 <- ggplot(phate_object)+
    geom_point(aes(PHATE1, PHATE2, color=combo_matrix[,feature3]), alpha=0.5, size=0.4)+
    ggtitle(feature3)+
    labs(color=feature3)+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_blank())+
    viridis::scale_color_viridis(option="B"))


plot4 <- ggplot(phate_object)+
  geom_point(aes(PHATE1, PHATE2, color=combo_matrix[,feature4]), alpha=0.5, size=0.4)+
  labs(color=feature4)+
  ggtitle(feature4)+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_blank())+
  viridis::scale_color_viridis(option="B")

# activation_phate <- cowplot::plot_grid(plot1, plot2, plot3, plot4, nrow = 2)
# ggsave("~/postdoc/scRNAseq/combo_analysis_figures/activation_phate.png", activation_phate, height=6, width=8, bg="white")


# memory_phate <- cowplot::plot_grid(plot1, plot2, plot3, plot4, nrow = 2)
# ggsave("~/postdoc/scRNAseq/combo_analysis_figures/memory_phate.png", memory_phate, height=6, width=8, bg="white")


combo[["phate"]] <- CreateDimReducObject(embeddings = phate_object$embedding, key = "phate_", assay = DefaultAssay(combo))

cluster_phate <- DimPlot(combo, reduction = "phate")
ggsave("~/postdoc/scRNAseq/combo_analysis_figures/cluster_phate.png", cluster_phate, height=3, width=5, bg="white")


top20_per_cluster <- data.frame(combo.markers %>%
                                  group_by(cluster) %>%
                                  slice_min(n = 20, order_by = p_val_adj))

View(top20_per_cluster)


c14_vs_11  <- FindMarkers(combo, ident.1 = 14, ident.2=12, min.pct = 0.25)
head(c14_vs_11, n = 30)


FeaturePlot(combo, order = TRUE, c("NKG7", "CCL5"), split.by = "time_n_infection")&viridis::scale_color_viridis(option = "B")






# frequency lollipop ####
library(dplyr)
library(ggplot2)


cluster_counts <- data.frame(table(combo@active.ident, combo@meta.data$sample_ID))
colnames(cluster_counts) <- c("Cluster_ID", "Sample_ID", "Count")
cluster_counts$Timepoint <- substr(cluster_counts$Sample_ID, 6, 14)
cluster_counts$Volunteer <- substr(cluster_counts$Sample_ID, 1, 4)
cluster_counts$N_Infection <- ifelse(cluster_counts$Volunteer %in% c("v313", "v315", "v320"), "First", "Third")



cluster_counts <- cluster_counts %>%
  group_by(Sample_ID) %>%
  mutate("Percentage"=Count/sum(Count))


baseline <- subset(cluster_counts, Timepoint=="Baseline" & N_Infection =="First")
T6 <- subset(cluster_counts, Timepoint=="T6" & N_Infection =="First")


prim_lolli <- ggplot()+
  geom_bar(data=baseline, aes(y=factor(Cluster_ID, levels = unique(baseline$Cluster_ID)), x=Percentage/3, fill=factor(Cluster_ID, levels=unique(baseline$Cluster_ID))), stat="identity", position="stack", width=1)+
  theme_minimal()+
  ggtitle("Baseline")+
  scale_x_reverse(name = "Percentage of CD4+ T cells\n", labels=scales::percent_format(accuracy = 1), expand = c(0,0),limits=c(0.3, 0))+
  guides(fill=guide_legend(label.position = "top", reverse = FALSE, nrow = 1,keyheight = unit(2, "mm"), keywidth = unit(4, "mm")))+
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


cluster_leg <- cowplot::get_legend(prim_lolli)

prim_lolli <- prim_lolli+theme(legend.position = "none")


ter_lolli <-   ggplot()+
  geom_bar(data=T6, aes(y=factor(Cluster_ID, levels = unique(baseline$Cluster_ID)), x=Percentage/3, fill=factor(Cluster_ID, levels=unique(baseline$Cluster_ID))), stat="identity", position="stack", width=1)+
  theme_minimal()+
  ggtitle("T6")+
  scale_x_continuous(name = "Percentage of CD4+ T cells\n", labels=scales::percent_format(accuracy = 1), expand = c(0,0), limits = c(0,0.3))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.placement = "outside",
        plot.margin = unit(c(1,0,1,5), "mm"),
        legend.position = "none",
        legend.text = element_text(size=6),
        legend.title = element_blank())


combo_lolli <- cowplot::plot_grid(prim_lolli, ter_lolli)
combo_lolli <- cowplot::plot_grid(combo_lolli, cluster_leg, ncol=1, rel_heights = c(5,1), rel_widths = c(1,1,2),align="v", axis = "tbrl")

# title <- cowplot::ggdraw()+cowplot::draw_label("First Infection", hjust = 0.5) 
title <- cowplot::ggdraw()+cowplot::draw_label("Third Infection", hjust = 0.5) 
combo_lolli <- cowplot::plot_grid(title, combo_lolli, ncol=1, rel_heights = c(0.1,1))

#ggsave(filename = "~/postdoc/scRNAseq/combo_analysis_figures/first_freq_lollipop.png", combo_lolli, width = 8, height = 5, bg="white", dpi=444)
#ggsave(filename = "~/postdoc/scRNAseq/combo_analysis_figures/third_freq_lollipop.png", combo_lolli, width = 8, height = 5, bg="white", dpi=444)


# pseudotime ####

library(monocle3)
combo.cds <- SeuratWrappers::as.cell_data_set(combo)

combo.cds <- monocle3::cluster_cells(cds=combo.cds, reduction_method = "UMAP", k=10, random_seed = 1234)
combo.cds <- monocle3::learn_graph(combo.cds, use_partition = TRUE, verbose=TRUE)

combo.cds <- monocle3::order_cells(combo.cds, reduction_method = "UMAP")


monocle3::plot_cells(
  cds = combo.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

# add pseudotime to seurat object
combo <- AddMetaData(
  object = combo,
  metadata = combo.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime"
)

#SeuratDisk::SaveH5Seurat(combo, "~/postdoc/scRNAseq/Seurat_Objects/20k_L003_hashed.h5Seurat", overwrite = TRUE)

#SeuratDisk::SaveH5Seurat(combo, "~/postdoc/scRNAseq/Seurat_Objects/pseudotime_dimred_processed_combo20k.H5Seurat")

FeaturePlot(combo, c("pseudotime"), split.by = "n_infection")&viridis::scale_color_viridis(option = "B")


# jackstraw ####

combo <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/pseudotime_dimred_processed_combo20k.H5Seurat")

# The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for each PC
# with a uniform distribution (dashed line). ‘Significant’ PCs will show a strong enrichment of features with
# low p-values (solid curve above the dashed line). In this case it appears that there is a sharp drop-off in
# significance after the first 10-12 PCs.

# test the dimensionality of the dataset
combo <- JackStraw(combo, num.replicate = 100)
combo <- ScoreJackStraw(combo, dims = 1:20)
JackStrawPlot(combo, dims = 1:20)
