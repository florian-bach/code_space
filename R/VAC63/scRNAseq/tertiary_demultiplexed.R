library(Seurat)

# debarcoding ####

# read primary data
pbmc.umis <- Read10X("~/postdoc/scRNAseq/tertiary/sample_feature_bc_matrix/")
# pbmc.umis <- Read10X("~/postdoc/scRNAseq/tertiary/sample_feature_bc_matrix/")
# remove random dash 1 from colnames
colnames(pbmc.umis) <- gsub("-1", "", colnames(pbmc.umis))

# For generating a hashtag count matrix from FASTQ files, please refer to
# https://github.com/Hoohm/CITE-seq-Count.  Load in the HTO count matrix
pbmc.htos <- Read10X("~/postdoc/scRNAseq/CITE_seq_count_output/trimmed_tertiary_whitelist/", gene.column = 1)
# pbmc.htos <- Read10X("~/postdoc/scRNAseq/CITE_seq_count_output/tertiary/umi_count/", gene.column = 1)
# pbmc.htos <- Read10X("~/postdoc/scRNAseq/CITE_seq_count_output/tertiary/umi_count/", gene.column = 1)

# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.
joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))

# Subset RNA and HTO counts by joint cell barcodes
pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])

# Confirm that the HTO have the correct names
#rownames(pbmc.htos) <- c("v313_Baseline", "v315_Baseline", "v320_Baseline", "v313_T6", "v315_T6", "v320_T6", "unmapped")
# rownames(pbmc.htos) <- c("v305_Baseline", "v308_Baseline", "v310_Baseline", "v305_T6", "v308_T6", "v310_T6", "unmapped")

# Setup Seurat object
pbmc.hashtag <- CreateSeuratObject(counts = pbmc.umis)

# Normalize RNA data with log normalization
pbmc.hashtag <- NormalizeData(pbmc.hashtag)
# Find and scale variable features
pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = "mean.var.plot")
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))


# Add HTO data as a new assay independent from RNA
pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = pbmc.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR", margin = 2)


pbmc.hashtag <- HTODemux(pbmc.hashtag, init = 7, assay = "HTO", positive.quantile = 0.999999999999, verbose = TRUE);table(pbmc.hashtag$HTO_classification.global)

Idents(pbmc.hashtag) <- "HTO_maxID"
Idents(pbmc.hashtag) <- "HTO_classification.global"

DefaultAssay(pbmc.hashtag) <- "HTO"
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = rownames(pbmc.hashtag),
                          verbose = FALSE)
pbmc.hashtag <- RunPCA(pbmc.hashtag, features = rownames(pbmc.hashtag), approx = FALSE)
pbmc.hashtag <- RunTSNE(pbmc.hashtag, dims = 1:7, perplexity = 100)

oligo_tsne <- FeaturePlot(pbmc.hashtag, features = gsub("_", "-", rownames(pbmc.htos)[1:6]), ncol = 3, label.size = 0.5)
ggplot2::ggsave("~/postdoc/scRNAseq/figures_for_Ange/tertiary_oligo_tsne.png", oligo_tsne, height = 7, width = 14, bg="white")


pbmc.hashtag <- FindNeighbors(pbmc.hashtag, reduction = "pca", dims = 1:6)
pbmc.hashtag <- FindClusters(pbmc.hashtag, resolution = 0.3, verbose = FALSE)
table(pbmc.hashtag@active.ident)

oligo_clusters <- DimPlot(pbmc.hashtag)
ggplot2::ggsave("~/postdoc/scRNAseq/figures_for_Ange/tertiary_oligo_clusters.png", oligo_clusters, height = 4, width = 6, bg="white")

DimPlot(pbmc.hashtag, label = TRUE)

FeaturePlot(pbmc.hashtag, features = "SELL")

pbmc.hashtag <- subset(pbmc.hashtag, seurat_clusters %in% c(0:5))


new.cluster.ids <- c("v308_Baseline", "v310_T6", "v305_T6", "v308_T6", "v310_Baseline", "v305_Baseline")
names(new.cluster.ids) <- levels(pbmc.hashtag)
pbmc.hashtag <- RenameIdents(pbmc.hashtag, new.cluster.ids)


clean_oligo_clusters <- DimPlot(pbmc.hashtag, label = TRUE)
ggplot2::ggsave("~/postdoc/scRNAseq/figures_for_Ange/tertiary_clean_oligo_clusters.png", clean_oligo_clusters, height = 4, width = 6, bg="white")

DefaultAssay(pbmc.hashtag) <- "RNA"


sample_ID <- unname(Idents(object = pbmc.hashtag))
volunteer <- substr(sample_ID, 1,4)
timepoint <- substr(sample_ID, 6, 14)

pbmc.hashtag <- AddMetaData(
  object = pbmc.hashtag,
  metadata = sample_ID,
  col.name = 'sample_ID'
)

pbmc.hashtag <- AddMetaData(
  object = pbmc.hashtag,
  metadata = volunteer,
  col.name = 'volunteer'
)
pbmc.hashtag <- AddMetaData(
  object = pbmc.hashtag,
  metadata = timepoint,
  col.name = 'timepoint'
)






pbmc.hashtag[["percent.mt"]] <- PercentageFeatureSet(pbmc.hashtag, pattern = "^MT-")

head(pbmc@meta.data, 5)

#plot QC
#VlnPlot(vac63c, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#plot QC relationship
# plot1 <- FeatureScatter(vac63c, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(vac63c, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2

#filter this bad boy
pbmc <- subset(pbmc.hashtag, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
SeuratDisk::SaveH5Seurat(pbmc, filename = "~/postdoc/scRNAseq/Seurat_Objects/demultiplexed_tertiary.h5Seurat", overwrite = TRUE)

# pre-processing: transformations & feature selection ####
library(Seurat)

pbmc <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/demultiplexed_primary.h5Seurat")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#find variable genes across cells
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# scaling data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# dimensionality reduction ####
#  Examine and visualize PCA results a few different ways
SeuratDisk::SaveH5Seurat(pbmc, filename = "~/postdoc/scRNAseq/Seurat_Objects/processed_demultiplexed_tertiary.h5Seurat", overwrite = TRUE)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
SeuratDisk::SaveH5Seurat(pbmc, filename = "~/postdoc/scRNAseq/Seurat_Objects/dimred_processed_demultiplexed_tertiary.h5Seurat", overwrite = TRUE)

pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 1, random.seed = 1234, algorithm = 1)


# find markers distinguishing clusters
# pbmc.markers <- FindAllMarkers(pbmc, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
# 
# top20_per_cluster <- data.frame(pbmc.markers %>%
#                                   group_by(cluster) %>%
#                                   slice_min(n = 20, order_by = p_val_adj))
# 
# View(top20_per_cluster)
# 

pbmc <- RunUMAP(pbmc, dims = 1:30, seed.use = 1234)

DimPlot(pbmc, reduction = "umap", split.by = "volunteer")


pbmc@meta.data$sample_ID <- factor(pbmc@meta.data$sample_ID, levels = c("v305_Baseline", "v308_Baseline", "v310_Baseline",  "v305_T6", "v308_T6","v310_T6"))
FeaturePlot(pbmc, c("CCL5", "GZMK"), split.by = "sample_ID")
FeaturePlot(pbmc, c("SELL", "CCR7"), split.by = "sample_ID")




cluster_counts <- data.frame(table(pbmc@active.ident, pbmc@meta.data$sample_ID))
colnames(cluster_counts) <- c("Cluster_ID", "Sample_ID", "Count")
cluster_counts$Timepoint <- substr(cluster_counts$Sample_ID, 6, 14)
cluster_counts$Volunteer <- substr(cluster_counts$Sample_ID, 1, 4)

cluster_counts <- cluster_counts %>%
  group_by(Sample_ID) %>%
  mutate("Percentage"=Count/sum(Count))


baseline <- subset(cluster_counts, Timepoint=="Baseline")
T6 <- subset(cluster_counts, Timepoint=="T6")

library(ggplot2)

prim_lolli <- ggplot()+
  geom_bar(data=baseline, aes(y=factor(Cluster_ID, levels = unique(baseline$Cluster_ID)), x=Percentage/3, fill=factor(Cluster_ID, levels=unique(baseline$Cluster_ID))), stat="identity", position="stack", width=1)+
  theme_minimal()+
  ggtitle("Baseline")+
  scale_x_reverse(name = "Percentage of CD4+ T cells\n", labels=scales::percent_format(accuracy = 1), expand = c(0,0),limits=c(0.2, 0))+
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
  scale_x_continuous(name = "Percentage of CD4+ T cells\n", labels=scales::percent_format(accuracy = 1), expand = c(0,0), limits = c(0,0.2))+
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


cluster13.markers <- FindMarkers(pbmc, ident.1 = 13, min.pct = 0.25)
cluster9.markers <- FindMarkers(pbmc, ident.1 = 9, min.pct = 0.25)

VlnPlot(pbmc, features = c("RPS4Y1", "EIF1AY"), split.by = "volunteer")
VlnPlot(pbmc, features = c("DDX3Y", "KDM5D"), split.by = "volunteer")
