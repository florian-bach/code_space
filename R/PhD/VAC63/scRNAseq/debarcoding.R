library(Seurat)


# read primary data
pbmc.umis <- Read10X("~/postdoc/scRNAseq/primary/sample_feature_bc_matrix/")
# pbmc.umis <- Read10X("~/postdoc/scRNAseq/tertiary/sample_feature_bc_matrix/")
# remove random dash 1 from colnames
colnames(pbmc.umis) <- gsub("-1", "", colnames(pbmc.umis))

# For generating a hashtag count matrix from FASTQ files, please refer to
# https://github.com/Hoohm/CITE-seq-Count.  Load in the HTO count matrix
pbmc.htos <- Read10X("~/postdoc/scRNAseq/CITE_seq_count_output/very_trimmed_primary_white/", gene.column = 1)
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

# If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
# clustering function for large applications You can also play with additional parameters (see
# documentation for HTODemux()) to adjust the threshold for classification Here we are using
# the default settings
# pbmc.hashtag <- MULTIseqDemux(pbmc.hashtag, assay = "HTO", autoThresh = FALSE, verbose = TRUE, maxiter = 7,
#                               quantile = 0.9999999999);table(pbmc.hashtag$MULTI_ID)
# 

pbmc.hashtag <- HTODemux(pbmc.hashtag, init = 7, assay = "HTO", positive.quantile = 0.999999999999, verbose = TRUE);table(pbmc.hashtag$HTO_classification.global)

table(pbmc.hashtag$HTO_classification)


# 
# # Group cells based on the max HTO signal
Idents(pbmc.hashtag) <- "HTO_maxID"
ridge_plot <- RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]]), ncol = 3)
ggplot2::ggsave("~/postdoc/scRNAseq/CITE_seq_count_output/tertiary/figures/debarcoding_ridge_plot_short_Read.png", ridge_plot, width=16, height=8, bg = "white")
# 
# 
# hashtag_scatter <- FeatureScatter(pbmc.hashtag, feature1 = "v308-Baseline", feature2 = "v310-Baseline")
# ggplot2::ggsave("~/postdoc/scRNAseq/CITE_seq_count_output/tertiary/figures/v308_Baseline_v310_Baseline.png", hashtag_scatter, height=4, width=8, units = "in", bg = "white")
# 
# 
# # compare number of UMIs for singlets doublets and negative cells
Idents(pbmc.hashtag) <- "HTO_classification.global"
# VlnPlot(pbmc.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
# 

# First, we will remove negative cells from the object

# Calculate a tSNE embedding of the HTO data
DefaultAssay(pbmc.hashtag) <- "HTO"
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = rownames(pbmc.hashtag),
                                 verbose = FALSE)
pbmc.hashtag <- RunPCA(pbmc.hashtag, features = rownames(pbmc.hashtag), approx = FALSE)
pbmc.hashtag <- RunTSNE(pbmc.hashtag, dims = 1:7, perplexity = 100)

FeaturePlot(pbmc.hashtag, features = gsub("_", "-", rownames(pbmc.htos)[1:6]))

pbmc.hashtag <- FindNeighbors(pbmc.hashtag, reduction = "pca", dims = 1:6)
pbmc.hashtag <- FindClusters(pbmc.hashtag, resolution = 0.3, verbose = FALSE)
table(pbmc.hashtag@active.ident)

DimPlot(pbmc.hashtag, label = TRUE)

FeaturePlot(pbmc.hashtag, features = "SELL")

pbmc.hashtag <- subset(pbmc.hashtag, seurat_clusters %in% c(0:5))


new.cluster.ids <- c("v320_Baseline", "v320_T6", "v315_Baseline", "v313_T6", "v315_T6", "v313_Baseline")
names(new.cluster.ids) <- levels(pbmc.hashtag)
pbmc.hashtag <- RenameIdents(pbmc.hashtag, new.cluster.ids)


DimPlot(pbmc.hashtag, label = TRUE)


DefaultAssay(pbmc.hashtag) <- "RNA"
pbmc.hashtag[["percent.mt"]] <- PercentageFeatureSet(pbmc.hashtag, pattern = "^MT-")

SeuratDisk::SaveH5Seurat(pbmc.hashtag, filename = "~/postdoc/scRNAseq/Seurat_Objects/demultiplexed_primary.h5Seurat")

#filter this bad boy
pbmc <- subset(pbmc.hashtag, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# pre-processing: transformations & feature selection ####
all.genes <- rownames(pbmc)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")

#find variable genes across cells
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, "RNA")

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), assay = "RNA")

pbmc <- RunUMAP(pbmc, dims = 1:30, seed.use = 1234)

pbmc[["Timepoint"]] <- PercentageFeatureSet(vac63c, pattern = "^MT-")

sample_ID <- unname(Idents(object = pbmc))
volunteer <- substr(sample_ID, 1,4)
timepoint <- substr(sample_ID, 6, 14)

pbmc <- AddMetaData(
  object = pbmc,
  metadata = sample_ID,
  col.name = 'sample_ID'
)

pbmc <- AddMetaData(
  object = pbmc,
  metadata = volunteer,
  col.name = 'volunteer'
)
pbmc <- AddMetaData(
  object = pbmc,
  metadata = timepoint,
  col.name = 'timepoint'
)


DimPlot(pbmc, group.by = "timepoint")


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




cluster_counts <- data.frame("gex1"=cluster_counts[,1],"gex2"=cluster_counts[,2],"gex3"=cluster_counts[,3], "gex4"=cluster_counts[,4])
cluster_counts$cluster_ID <- rownames(cluster_counts)



pbmc.hashtag <- RunTSNE(pbmc.singlet, reduction = "pca", dims = 1:10)




HTOHeatmap(pbmc.hashtag, assay = "HTO", ncells = 5000)




# Extract the singlets
pbmc.singlet <- subset(pbmc.hashtag, idents = "Singlet")

# Select the top 1000 most variable features
pbmc.singlet <- FindVariableFeatures(pbmc.singlet, selection.method = "mean.var.plot")

# Scaling RNA data, we only scale the variable features here for efficiency
pbmc.singlet <- ScaleData(pbmc.singlet, features = VariableFeatures(pbmc.singlet))

# Run PCA
pbmc.singlet <- RunPCA(pbmc.singlet, features = VariableFeatures(pbmc.singlet))

# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
pbmc.singlet <- FindNeighbors(pbmc.singlet, reduction = "pca", dims = 1:10)
pbmc.singlet <- FindClusters(pbmc.singlet, resolution = 0.6, verbose = FALSE)
pbmc.singlet <- RunTSNE(pbmc.singlet, reduction = "pca", dims = 1:10)

# Projecting singlet identities on TSNE visualization
DimPlot(pbmc.singlet, group.by = "HTO_classification")
