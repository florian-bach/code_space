library(Seurat)

# mtx <- list.files(path = "/labs/prasj/BIG_Flo/boyle_reanalysis", pattern = "matrix.mtx.gz", full.names = TRUE)
# cells <- list.files(path = "/labs/prasj/BIG_Flo/boyle_reanalysis", pattern = "barcodes.tsv.gz", full.names = TRUE)
# features <- list.files(path = "/labs/prasj/BIG_Flo/boyle_reanalysis", pattern = "features.tsv.gz", full.names = TRUE)

mtx <- list.files(path = "~/Downloads/GSE217930_RAW/unzipped/", pattern = "matrix.mtx.gz", full.names = TRUE)
cells <- list.files(path = "~/Downloads/GSE217930_RAW/unzipped/", pattern = "barcodes.tsv.gz", full.names = TRUE)
features <- list.files(path = "~/Downloads/GSE217930_RAW/unzipped/", pattern = "features.tsv.gz", full.names = TRUE)

pbmc_counts <- ReadMtx(mtx = mtx[1], cells = cells[1], features = features[1])
pbmc_seurat <- CreateSeuratObject(counts = pbmc_counts)

for(i in 2:5){
  counts <- ReadMtx(mtx = mtx[i], cells = cells[i], features = features[i])
  seurat_object <- CreateSeuratObject(counts = counts)
  pbmc_seurat <- merge(pbmc_seurat, seurat_object)
}

saveRDS(pbmc_seurat, "/labs/prasj/BIG_Flo/boyle_reanalysis/raw_seurat_object")

pbmc_seurat[["percent.mt"]] <- PercentageFeatureSet(pbmc_seurat, pattern = "^MT-")
#no max number of genes??
pbmc_seurat <- subset(pbmc_seurat, subset = nFeature_RNA > 220 & percent.mt < 20)

pbmc_seurat <- NormalizeData(pbmc_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_seurat <- FindVariableFeatures(pbmc_seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc_seurat)
pbmc_seurat <- ScaleData(pbmc_seurat, features = all.genes)
pbmc_seurat <- RunPCA(pbmc_seurat, features = VariableFeatures(object = pbmc_seurat))

pbmc_seurat <- FindNeighbors(pbmc_seurat, dims = 1:10)
pbmc_seurat <- FindClusters(pbmc_seurat, resolution = 0.5, method=4)
pbmc_seurat <- RunUMAP(pbmc_seurat, dims = 1:10)
DimPlot(pbmc_seurat, reduction = "umap")



table(
  substr(
    names(pbmc_seurat$orig.ident),
    nchar(pbmc_seurat$orig.ident)-3,
    nchar(pbmc_seurat$orig.ident)
    )
  )
