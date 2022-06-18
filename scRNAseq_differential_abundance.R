combo <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/pseudotime_dimred_processed_combo20k.H5Seurat")
# combo <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/no_jackstraw_demulitplexerd_all_processed.h5Seurat")
# combo <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/activated_subset.h5Seurat")

combo.sce <- Seurat::as.SingleCellExperiment(combo)

# sample_id timepoint   batch volunteer n_cells
# 1   v02_Baseline  Baseline batch_1       v02   42045
# 2        v02_C10       C10 batch_1       v02   28507
# 3  v02_Diagnosis Diagnosis batch_1       v02   14037
# 4         v02_T6        T6 batch_1       v02   33604

# 
# (Intercept) timepointC10 timepointDiagnosis timepointT6 volunteerv03 volunteerv05 volunteerv06
# 1            1            0                  0           0            0            0            0
# 2            1            1                  0           0            0            0            0
# 3            1            0                  1           0            0            0            0
# 4            1            0                  0           1            0            0            0
# 5            1            0                  0           0            1            0            0
sum_by <- c("seurat_clusters", "sample_ID")


summed <- scuttle::aggregateAcrossCells(combo.sce, id=SummarizedExperiment::colData(combo.sce)[,sum_by])
colnames(summed) <- apply(SummarizedExperiment::colData(summed)[,sum_by], 1, function(x)paste(x, collapse="_"))

metadata <- data.frame("column_name"=colnames(SummarizedExperiment::assay(summed, "counts")))
metadata$cluster_ID <- substr(metadata$column_name, 1, 2)
metadata$cluster_ID <- gsub("_", "", metadata$cluster_ID)

metadata$timepoint <- ifelse(grepl("T6", metadata$column_name), "T6", "Baseline")

metadata$volunteer <- substr(metadata$column_name, regexpr("v", metadata$column_name), regexpr("v", metadata$column_name)+3)

metadata$n_infection <- ifelse(metadata$volunteer %in% c("v313", "v315", "v320"), "First", "Third")

metadata$sample_type <- paste(metadata$timepoint, metadata$n_infection, sep="_")
metadata$sample_id <- paste(metadata$volunteer, metadata$timepoint, sep="_")

metadata <- metadata[,-c(1:2)]
metadata <- metadata[!duplicated(metadata),]
rownames(metadata) <- metadata$sample_id

library(SingleCellExperiment)

metadata(combo.sce)$experiment_info <-metadata
metadata(combo.sce)$cluster_codes <- data.frame("seurat_clusters"=combo.sce$seurat_clusters)

combo.sce[["cluster_id"]] <- combo.sce$seurat_clusters
combo.sce[["sample_id"]] <- combo.sce$sample_ID

#exprs(combo.sce) <- matrix(assay(combo.sce, "scaledata"))

diffcyt::calcCounts(combo.sce)




design <- diffcyt::createDesignMatrix(metadata, c("timepoint", "n_infection"))

FDR_cutoff <- 0.05

t6_first_contrast <- createContrast(c(0, 1, 0))
t6_third_contrast <- createContrast(c(0, 1, 1))
base_third_contrast <- createContrast(c(0, 0, 1))



t6_first <- diffcyt(combo.sce,
                  design = design,
                  experiment_info = metadata,
                  contrast = t6_first_contrast,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "seurat_clusters",
                  verbose = T)


t6_third <- diffcyt(combo.sce,
                  design = design,
                  contrast = t6_third_contrast,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "flo_merge",
                  verbose = T)

base_third <- diffcyt(combo.sce,
                  design = design,
                  contrast = base_third_contrast,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "flo_merge",
                  verbose = T)



# milo ####



combo <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/pseudotime_dimred_processed_combo20k.H5Seurat")

library(Seurat)

activated_clusters <- subset(combo,  seurat_clusters %in% c(6:9, 10:12))
combo <- NULL
gc()
activated_clusters <- FindVariableFeatures(activated_clusters, selection.method = "vst", nfeatures = 1000)

# all.genes <- rownames(combo)
# combo <- ScaleData(combo)


activated_clusters <- RunPCA(activated_clusters, features = VariableFeatures(object = activated_clusters))

activated_clusters <- FindNeighbors(activated_clusters, dims = 1:20)
activated_clusters <- FindClusters(activated_clusters, resolution = 2, random.seed = 1234, algorithm = 1)


combo.sce <- Seurat::as.SingleCellExperiment(activated_clusters)

sum_by <- c("seurat_clusters", "sample_ID")


#write.csv(metadata, "~/postdoc/scRNAseq/metadata/sce_metadata.csv", row.names = TRUE, col.names = TRUE)

library(miloR)



combo_milo <- Milo(combo.sce)
combo_milo <- buildGraph(combo_milo, k = 30, d = 30, reduced.dim = "PCA")
combo_milo <- makeNhoods(combo_milo, prop = 0.1, k = 30, d=30, refined = TRUE, reduced_dims = "PCA")

# need this for colData & owData
library(SingleCellExperiment)

metadata <- read.csv("~/postdoc/scRNAseq/metadata/sce_metadata.csv", row.names = 1)

combo_milo <- countCells(combo_milo, meta.data = as.data.frame(colData(combo_milo)), sample="sample_ID")

combo_milo <- calcNhoodDistance(combo_milo, d=30, reduced.dim = "PCA")

da_results <- testNhoods(combo_milo, design = ~ timepoint + volunteer, design.df = metadata)

combo_milo <- buildNhoodGraph(combo_milo)

library(patchwork)
library(ggplot2)

## Plot single-cell UMAP
umap_pl <- scater::plotReducedDim(combo_milo, dimred = "UMAP", colour_by="timepoint",
                          text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(combo_milo, da_results, layout="UMAP",alpha=0.1)

umap_pl + nh_graph_pl +
  plot_layout(guides="collect")


# add cluster id to network

da_results <- annotateNhoods(combo_milo, da_results, coldata_col = "seurat_clusters")

plotDAbeeswarm(da_results, group.by = "seurat_clusters")


# establish key marker genes for each neighbourhood
combo_milo <- scuttle::logNormCounts(combo_milo)

da_results$NhoodGroup <- as.numeric(da_results$SpatialFDR < 0.1 & da_results$logFC < 0)
da_nhood_markers <- findNhoodGroupMarkers(combo_milo, da_results, subset.row = Seurat::VariableFeatures(activated_clusters), 
                                          aggregate.samples = TRUE, sample_col = "sample_ID")


## Warning: Zero sample variances detected, have been offset away from zero

## Warning: Zero sample variances detected, have been offset away from zero

head(da_nhood_markers)
