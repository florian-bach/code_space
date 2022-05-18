# preprocessing ####
library(dplyr)
library(Seurat)
library(patchwork)

# test_data <- Read10X(data.dir = "~/postdoc/scRNAseq/cell_ranger_adventure_7/sample_feature_bc_matrix/")
test_data <- Read10X(data.dir = "~/postdoc/scRNAseq/tertiary/sample_feature_bc_matrix/")
vac63c <- CreateSeuratObject(counts = test_data, project = "vac63c_cd4_tcells", min.cells = 3, min.features = 200)
# 

# Quality Control ####
# Sort out cells with high mitochondrial reads (they are dying); then get rid of cells with less than 200 or more than 2500 reads

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
vac63c[["percent.mt"]] <- PercentageFeatureSet(vac63c, pattern = "^MT-")

head(vac63c@meta.data, 5)

#plot QC
#VlnPlot(vac63c, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#plot QC relationship
# plot1 <- FeatureScatter(vac63c, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(vac63c, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2

#filter this bad boy
pbmc <- subset(vac63c, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# pre-processing: transformations & feature selection ####

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#find variable genes across cells
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#top10 <- head(VariableFeatures(pbmc), 30)

# plot variable features with and without labels
# plot1 <- VariableFeaturePlot(pbmc)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# scaling data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# dimensionality reduction ####
#  Examine and visualize PCA results a few different ways
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# clustering ####

#print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

#VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
# 
# DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(pbmc, dims = 1:9, cells = 500, balanced = TRUE)



# The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for each PC
# with a uniform distribution (dashed line). ‘Significant’ PCs will show a strong enrichment of features with
# low p-values (solid curve above the dashed line). In this case it appears that there is a sharp drop-off in
# significance after the first 10-12 PCs.

# test the dimensionality of the dataset
# pbmc <- JackStraw(pbmc, num.replicate = 100)
# pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
# JackStrawPlot(pbmc, dims = 1:20)


# pbmc <- FindNeighbors(pbmc, dims = 1:20)
# pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 1, random.seed = 1234, algorithm = 1)

# find markers distinguishing clusters
pbmc.markers <- FindAllMarkers(pbmc, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

top20_per_cluster <- data.frame(pbmc.markers %>%
  group_by(cluster) %>%
  slice_min(n = 20, order_by = p_val_adj))

View(top20_per_cluster)

 
# top20_per_cluster <- data.frame(pbmc.markers %>%
#                                  group_by(cluster) %>%
#                                  slice_max(n = 20, order_by = avg_log2FC))


# visualising markers across clusters ####
# cluster_7_markers <- pbmc.markers %>%
#   filter(cluster==7) %>%
#   slice_max(n = 20, order_by = avg_log2FC)
#   
#visualise individual gene expression across clusters
# by cytof, the average of naive cells across baseline and T6 is ~54.5%
VlnPlot(pbmc, features = c("RGS16", "EGR1", "EGR2", "EGR3"), stack = TRUE)
VlnPlot(pbmc, features = c("CD160", "IKZF2", "CD200", "TNFRSF9"), stack = TRUE)

VlnPlot(pbmc, features = c("CCR7", "SELL", "FOS", "JUN", "LMNA"), stack = TRUE)
VlnPlot(pbmc, features = c("TXNIP", "IL32", "LTB", "ODC1", "KLRB1"), stack = TRUE)



# cluster 3, 10 --> EM (high  LMNA, low  SELL)
# cluster 6 --> naive (high SELL, low  jun/fos)
# cluster 9 --> CM (high SELL, high LMNA)
# cluster 8, 10, activated
# 0, 1, 2 naive? (low in CSRNP1, jun/fos)
# 11 translation_repressed
# 5 memory
# 4 slightly activated (IL32, KLRB1, B2M up)
# 7 KLRB1+IL7R+ memory

cluster5.markers <- FindMarkers(pbmc, ident.1 = 7, min.pct = 0.25)
head(cluster5.markers, n = 30)

cluster7.markers <- FindMarkers(pbmc, ident.1 = 7, min.pct = 0.25)
head(cluster7.markers, n = 30)

cluster4.markers <- FindMarkers(pbmc, ident.1 = 4, ident.2 = 5, min.pct = 0.25)
head(cluster4.markers, n = 30)



# naive_markers <- c("CD3D", "CD4", "CCR7", "SELL", "LEF1", "PTPRC", "TIMP1", "TNFRSF4", "ATM", "KLF2", "ITGA6")
# VlnPlot(pbmc, features = naive_markers)
# 
# treg_markers <- c("FOXP3", "IL2RA", "CCR10", "SELPLG", "IL7R", "ENTPD1")
# VlnPlot(pbmc, features = treg_markers)
# 
# activation_markers <- c("CD38", "BCL2", "CD28", "PDCD1", "CTLA4", "GNLY", "PRF1", "GZMH", "LAG3", "TIGIT")
# VlnPlot(pbmc, features = activation_markers)
# 
# chemokine_receptors <- c("CXCR3", "CXCR4", "CXCR5", "CXCR6", "CCR4", "CCR5", "CCR6", "CCR7", "CX3CR1")
# VlnPlot(pbmc, features = chemokine_receptors)
# 
# naive_central_memory <- c("LEF1", "ATM", "SELL", "KLF2", "ITGA6")
# VlnPlot(pbmc, features = naive_central_memory)
# 
# cytotoxic_markers <- c("CD3", "NK67", "GZMH", "GRNLY", "KLRG1", "KLRB1")
# VlnPlot(pbmc, features = cytotoxic_markers)
# 
# cytokines <- c("IL2", "IFNG", "CCL5", "CCL4", "TNF", "IL21", "IL10", "IL15", "IL33", "CSF1", "CSF3")
# VlnPlot(pbmc, features = cytokines)
# 
# memory_markers <-c("CD27", "CD28", "CD44", "SELL", "IL7R")
# VlnPlot(pbmc, features = memory_markers)

# find markers for every cluster compared to all remaining cells, report only the positive
# onesVlnPlot(pbmc, features = c("MS4A1", "CD79A"))

cluster5.markers <- FindMarkers(pbmc, ident.1 = "translation_repressed", ident.2 = "effector memory2?", min.pct = 0.25)
head(cluster5.markers, n = 30)


# dimensionality reduction & pseudotime ####
#pbmc <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/L003_hashed.h5Seurat")
pbmc <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/20k_L003_hashed.h5Seurat")

#smaller_pbmc <- pbmc[, sample(colnames(pbmc), size =20000, replace=F)]
SeuratDisk::SaveH5Seurat(pbmc, "~/postdoc/scRNAseq/Seurat_Objects/vac63c_primary.h5Seurat", overwrite = TRUE)


pbmc <- RunUMAP(pbmc, dims = 1:30, seed.use = 1234)

DimPlot(pbmc, reduction = "umap")

FeaturePlot(pbmc, c("NKG7", "CCL5", "GZMB", "GZMK"))
FeaturePlot(pbmc, c("SELL","CCR7", "JUNB", "LMNA"))

# pseudotime
library(monocle3)
pbmc.cds <- SeuratWrappers::as.cell_data_set(pbmc)

pbmc.cds <- monocle3::cluster_cells(cds=pbmc.cds, reduction_method = "UMAP", k=10, random_seed = 1234)
pbmc.cds <- monocle3::learn_graph(pbmc.cds, use_partition = TRUE, verbose=TRUE)

pbmc.cds <- monocle3::order_cells(pbmc.cds, reduction_method = "UMAP")


monocle3::plot_cells(
  cds = pbmc.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

# add pseudotime to seurat object
pbmc <- AddMetaData(
  object = pbmc,
  metadata = pbmc.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime"
)

#SeuratDisk::SaveH5Seurat(pbmc, "~/postdoc/scRNAseq/Seurat_Objects/20k_L003_hashed.h5Seurat", overwrite = TRUE)


FeaturePlot(pbmc, c("pseudotime", "SELL", "LMNA"))&viridis::scale_color_viridis(option = "B")


# RNA velocity ####

library(velocyto.R)
pbmc <- SCTransform(object = pbmc, assay = "spliced")

pbmc <- filter.genes.by.cluster.expression(pbmc,cluster.label,min.max.cluster.average = 0.2)

r <- pagoda2::Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T)

cluster.label <- r$clusters$PCA$multilevel # take the cluster factor that was calculated by p2
cell.colors <- pagoda2:::fac2col(cluster.label)


nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(nmat)))

pbmc <- SeuratWrappers::RunVelocity(object = pbmc, deltaT = 1, kCells = 25, fit.quantile = 0.02)


ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
names(x = ident.colors) <- levels(x = bm)
cell.colors <- ident.colors[Idents(object = bm)]
names(x = cell.colors) <- colnames(x = bm)
show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), vel = Tool(object = bm, 
                                                                                             slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)


# UMAP
pbmc <- RunUMAP(pbmc, dims = 1:30, seed.use = 1234)

FeaturePlot(pbmc, features = c("RGS16", "EGR1", "EGR2", "EGR3"))
FeaturePlot(pbmc, features = c("PRDM1", "IKZF2", "XCL1", "TNFRSF9"))



FeaturePlot(pbmc, features = c("CCL5", "GZMA", "SELL", "TNF", "EGR1"))

FeatureScatter(pbmc, feature1 = "CCL5", feature2 = "GZMH")

FeatureScatter(pbmc, feature1 = "NKG7", feature2 = "GZMH")

FeatureScatter(pbmc, feature1 = "JUN", feature2 = "FOS")

FeatureScatter(pbmc, feature1 = "TNF", feature2 = "TGFB1")



test_vector <- c(paste("IL", seq(1,7), sep=''),
                 paste("IL", letters[5:11], sep = ''))


all_interleukins <- grep('IL[0-9]', rownames(pbmc@assays$RNA), value = TRUE)
all_interleukins_plot <- VlnPlot(pbmc, features = all_interleukins, ncol = 9)
ggplot2::ggsave("~/postdoc/scRNAseq/exploratory_plots/all_interleukins_plot.png", all_interleukins_plot, width=20, height=20)

all_chemokines <- c(grep('CCL[0-9]', rownames(pbmc@assays$RNA), value = TRUE),
                    grep('CXCL[0-9]', rownames(pbmc@assays$RNA), value = TRUE))

all_chemokines_plot <- VlnPlot(pbmc, features = all_chemokines, ncol = 5)
ggplot2::ggsave("~/postdoc/scRNAseq/exploratory_plots/all_chemokines_plot.png", all_chemokines_plot, width=12, height=6)


new.cluster.ids <- c("central memory?", "effector memory?", "effector memory2?", "activated+", "naive", "intermediate activated", "cytotoxic effectors", 
                     "translation_repressed")

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

FeaturePlot(pbmc, features = "EEF1A1")




  write.csv(pbmc@meta.data,"~/postdoc/scRNAseq/exploratory_plots/metadata/sliver_seurat_metadata.csv")
