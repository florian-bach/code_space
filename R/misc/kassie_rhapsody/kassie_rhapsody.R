library(Seurat)

# p258a <- readRDS("~/postdoc/stanford/misc/kassie_rhapsody/Files from Seven Bridges-selected/p258a_CKDL240012182-1A_225CC7LT4_L2_1_Seurat.rds")
# p258b <- readRDS("~/postdoc/stanford/misc/kassie_rhapsody/Files from Seven Bridges-selected/p258b_CKDL240012182-1A_225CC7LT4_L2_1_Seurat.rds")
# combo <- merge(p258a, y=p258b, add.cell.ids=c("p258a", "p258b"), project="bd_rhapsody_spleen")
# saveRDS(combo, "~/postdoc/stanford/misc/kassie_rhapsody/combo.rds")

combo <- readRDS("~/postdoc/stanford/misc/kassie_rhapsody/combo.rds")
combo$seq_batch <- ifelse(grepl("*p258a_*", colnames(combo)), "p258a", "p258b")
  
combo[["percent.mt"]] <- PercentageFeatureSet(combo, pattern = "^MT-")

VlnPlot(combo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# this step loses 18244 cells, ~19%
combo <- subset(combo, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

oligos <- grep("pAbO", rownames(combo))
combo[["oligos"]] <- CreateAssayObject(data = combo[oligos,])

combo <- FindVariableFeatures(combo, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(combo), 100)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(combo)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(combo)
combo <- ScaleData(combo, features = VariableFeatures(object = combo))
combo <- RunPCA(combo, features = VariableFeatures(object = combo))

combo <- FindNeighbors(combo, dims = 1:10)
combo <- FindClusters(combo, resolution = 0.5)

combo <- RunUMAP(combo, dims = 1:10)
DimPlot(combo, reduction = "umap")