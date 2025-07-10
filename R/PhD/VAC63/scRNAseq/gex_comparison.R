library(Seurat)

gex1_data <- Read10X(data.dir = "~/postdoc/scRNAseq/gex_comparison/gex1/sample_feature_bc_matrix/")
gex1 <- CreateSeuratObject(counts = gex1_data, project = "gex1", min.cells = 3, min.features = 200)
gex1 <- gex1[, sample(colnames(gex1), size =6666, replace=F)]
SeuratDisk::SaveH5Seurat(gex1, "~/postdoc/scRNAseq/Seurat_Objects/gex1_6666.h5Seurat", overwrite = TRUE)

gex2_data <- Read10X(data.dir = "~/postdoc/scRNAseq/primary/sample_feature_bc_matrix/")
gex2 <- CreateSeuratObject(counts = gex2_data, project = "gex2", min.cells = 3, min.features = 200)
gex2 <- gex2[, sample(colnames(gex2), size =6666, replace=F)]
SeuratDisk::SaveH5Seurat(gex2, "~/postdoc/scRNAseq/Seurat_Objects/gex2_6666.h5Seurat", overwrite = TRUE)


gex3_data <- Read10X(data.dir = "~/postdoc/scRNAseq/gex_comparison/gex3/sample_feature_bc_matrix/")
gex3 <- CreateSeuratObject(counts = gex3_data, project = "gex3", min.cells = 3, min.features = 200)
gex3 <- gex3[, sample(colnames(gex3), size =6666, replace=F)]
SeuratDisk::SaveH5Seurat(gex3, "~/postdoc/scRNAseq/Seurat_Objects/gex3_6666.h5Seurat", overwrite = TRUE)

gex4_data <- Read10X(data.dir = "~/postdoc/scRNAseq/gex_comparison/gex3/sample_feature_bc_matrix/")
gex4 <- CreateSeuratObject(counts = gex4_data, project = "gex4", min.cells = 3, min.features = 200)
gex4 <- gex4[, sample(colnames(gex4), size =6666, replace=F)]
SeuratDisk::SaveH5Seurat(gex4, "~/postdoc/scRNAseq/Seurat_Objects/gex4_6666.h5Seurat", overwrite = TRUE)

gex1 <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/gex1_6666.h5Seurat")
gex2 <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/gex2_6666.h5Seurat")
gex3 <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/gex3_6666.h5Seurat")
gex4 <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/gex4_6666.h5Seurat")

#quad_combo_object <- merge(gex1, gex2, gex3, gex4, add.cell.ids = c("gex1", "gex2", "gex3", "gex4"), project = "Combo_Pool")


tertiary_combo <- merge(gex1, c(gex2, gex3, gex4), add.cell.ids = c("gex1", "gex2", "gex3", "gex4"), project = "Combo_Pool")


tertiary_combo[["percent.mt"]] <- PercentageFeatureSet(tertiary_combo, pattern = "^MT-")


tertiary_combo <- subset(tertiary_combo, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


tertiary_combo <- NormalizeData(tertiary_combo, normalization.method = "LogNormalize", scale.factor = 10000)
  
#find variable genes across cells
tertiary_combo <- FindVariableFeatures(tertiary_combo, selection.method = "vst", nfeatures = 5000)
  
all.genes <- rownames(tertiary_combo)
tertiary_combo <- ScaleData(tertiary_combo, features = all.genes)

# dimensionality reduction ####
#  Examine and visualize PCA results a few different ways
tertiary_combo <- RunPCA(tertiary_combo, features = VariableFeatures(object = tertiary_combo))


tertiary_combo <- FindNeighbors(tertiary_combo, dims = 1:40)
tertiary_combo <- FindClusters(tertiary_combo, resolution = 1.5, random.seed = 1234, algorithm = 1)


tertiary_combo <- RunUMAP(tertiary_combo, dims = 1:40, seed.use = 1234)

cluster_umap <- DimPlot(tertiary_combo, split.by = "orig.ident", label = TRUE)
DimPlot(tertiary_combo, label = TRUE)

FeaturePlot(tertiary_combo, c("CS38"), split.by = "orig.ident")

variable_markers <- FindAllMarkers(tertiary_combo, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

top30_per_cluster <- data.frame(variable_markers %>%
                                    group_by(cluster) %>%
                                    slice_min(n = 30, order_by = p_val_adj))
  


cluster_counts <- table(tertiary_combo@active.ident, tertiary_combo@meta.data$orig.ident)
cluster_counts <- data.frame("gex1"=cluster_counts[,1],"gex2"=cluster_counts[,2],"gex3"=cluster_counts[,3], "gex4"=cluster_counts[,4])
cluster_counts$cluster_ID <- rownames(cluster_counts)


long_cluster_counts <- tidyr::pivot_longer(cluster_counts, colnames(cluster_counts)[1:4], "Library")

VlnPlot(tertiary_combo, features = c("CCL5", "GNLY", "NKG7"), split.by = "orig.ident", stack = TRUE, cols = c("orange", "darkorange", "blue", "darkblue"))

ggplot(long_cluster_counts, aes(x=factor(cluster_ID, levels = as.character(0:19)), y=value/6666, fill=Library))+
  geom_bar(stat="identity", position="dodge")+
  scale_y_continuous(labels = scales::label_percent())+
  scale_fill_manual(values = c("orange", "darkorange", "blue", "darkblue"))+
  theme_minimal()



library(ggplot2)

prim_lolli <- ggplot()+
  geom_bar(data=cluster_counts,aes(y=factor(cluster_ID, levels=rownames(cluster_counts)), x=gex3/6666, fill=factor(cluster_ID, levels=rownames(cluster_counts))), stat="identity", position="stack", width=1)+
  theme_minimal()+
  ggtitle("gex3")+
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
  geom_bar(data=cluster_counts,aes(y=factor(cluster_ID, levels=rownames(cluster_counts)), x=gex4/6666, fill=factor(cluster_ID, levels=rownames(cluster_counts))), stat="identity", position="stack", width=1)+
  theme_minimal()+
  ggtitle("gex4")+
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

VlnPlot(tertiary_combo, features = c("NCAM1"), split.by = "orig.ident")

