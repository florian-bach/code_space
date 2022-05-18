library(Seurat)

# first_data <- Read10X(data.dir = "~/postdoc/scRNAseq/primary/sample_feature_bc_matrix/")
# primary <- CreateSeuratObject(counts = first_data, project = "First", min.cells = 3, min.features = 200)
# smaller_primary <- primary[, sample(colnames(primary), size =10000, replace=F)]
# SeuratDisk::SaveH5Seurat(combo_object, "~/postdoc/scRNAseq/Seurat_Objects/20k_combo_first_third.h5Seurat", overwrite = TRUE)
# 
# third_data <- Read10X(data.dir = "~/postdoc/scRNAseq/tertiary/sample_feature_bc_matrix/")
# tertiary <- CreateSeuratObject(counts = third_data, project = "Third", min.cells = 3, min.features = 200)
# smaller_tertiary <- tertiary[, sample(colnames(tertiary), size =10000, replace=F)]
# SeuratDisk::SaveH5Seurat(smaller_tertiary, "~/postdoc/scRNAseq/Seurat_Objects/10k_tertiary.h5Seurat", overwrite = TRUE)
# 




# smaller_primary <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/10k_primary.h5Seurat")
# smaller_tertiary <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/10k_tertiary.h5Seurat")
# 
# combo_object <- merge(smaller_primary, smaller_tertiary, add.cell.ids = c("First", "Third"), project = "Combo_Pool")
# 
# SeuratDisk::SaveH5Seurat(combo_object, "~/postdoc/scRNAseq/Seurat_Objects/20k_comboy.h5Seurat", overwrite = TRUE)

combo_object <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/20k_comboy.h5Seurat")

combo_object[["percent.mt"]] <- PercentageFeatureSet(combo_object, pattern = "^MT-")


combo_object <- subset(vac63c, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# pre-processing: transformations & feature selection ####

combo_object <- NormalizeData(combo_object, normalization.method = "LogNormalize", scale.factor = 10000)

#find variable genes across cells
combo_object <- FindVariableFeatures(combo_object, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(combo_object)
combo_object <- ScaleData(combo_object, features = all.genes)

# dimensionality reduction ####
#  Examine and visualize PCA results a few different ways
combo_object <- RunPCA(combo_object, features = VariableFeatures(object = combo_object))


combo_object <- FindNeighbors(combo_object, dims = 1:20)
combo_object <- FindClusters(combo_object, resolution = 1, random.seed = 1234, algorithm = 1)


combo_object <- RunUMAP(combo_object, dims = 1:30, seed.use = 1234)

cluster_umap <- DimPlot(combo_object, split.by = "orig.ident")
ggplot2::ggsave("~/postdoc/scRNAseq/combo_analysis_figures/cluster_umap.png", cluster_umap, width=8, height=4)

activation_plot <- FeaturePlot(combo_object, c("NKG7", "CCL5", "GZMB", "CD38"), split.by = "orig.ident")
activation_plot2 <- FeaturePlot(combo_object, c("SELL", "TNF", "GZMH", "GNLY"), split.by = "orig.ident")
cluster0_plot <- FeaturePlot(combo_object, c("YPEL5", "SBDS", "DDIT4", "CSKMT"), split.by = "orig.ident")


ggplot2::ggsave("~/postdoc/scRNAseq/combo_analysis_figures/activation_plot.png", activation_plot, width=8, height = 14)
ggplot2::ggsave("~/postdoc/scRNAseq/combo_analysis_figures/activation_plot2.png", activation_plot2, width=8, height = 14)
ggplot2::ggsave("~/postdoc/scRNAseq/combo_analysis_figures/cluster0_plot.png", cluster0_plot, width=8, height = 14)


q <- table(combo_object@active.ident, combo_object@meta.data$orig.ident)
cluster_counts <- data.frame("First"=cluster_counts[,1], "Third"=cluster_counts[,2])
cluster_counts$cluster_ID <- rownames(cluster_counts)



library(ggplot2)

prim_lolli <- ggplot()+
  geom_bar(data=cluster_counts,aes(y=factor(cluster_ID, levels=rownames(cluster_counts)), x=First/10000, fill=factor(cluster_ID, levels=rownames(cluster_counts))), stat="identity", position="stack", width=1)+
  theme_minimal()+
  ggtitle("First Infection")+
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
  geom_bar(data=cluster_counts,aes(y=factor(cluster_ID, levels=rownames(cluster_counts)), x=Third/10000, fill=factor(cluster_ID, levels=rownames(cluster_counts))), stat="identity", position="stack", width=1)+
  theme_minimal()+
  ggtitle("Third Infection")+
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
ggsave("~/postdoc/scRNAseq/combo_analysis_figures/combo_lollipop.png", combo_lolli, height=4, width=12, bg="white")

library(dplyr)

pbmc.markers <- FindAllMarkers(combo_object, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

top20_per_cluster <- data.frame(pbmc.markers %>%
                                  group_by(cluster) %>%
                                  slice_min(n = 30, order_by = p_val_adj))

View(top20_per_cluster)

cluster5.markers <-  FindMarkers(combo_object, ident.1 = 14, ident.2 = c(12,13), min.pct = 0.25)
head(cluster5.markers, n = 30)

VlnPlot(combo_object, features = c("TNF", "LTB", "GNLY"))
