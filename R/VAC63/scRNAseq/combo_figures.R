library(Seurat)


combo <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/pseudotime_dimred_processed_combo20k.H5Seurat")

# overview figures ####

cluster_umap <- DimPlot(combo, reduction = "umap", label=TRUE)
ggplot2::ggsave("~/postdoc/scRNAseq/combo_analysis_figures/combo_cluster_umap.png", cluster_umap, height=4, width=6, bg="white")


memory_umap <- FeaturePlot(combo, features = c("SELL", "CCR7", "IL7R"), ncol=3)
activation_umap <- FeaturePlot(combo, features = c("GZMB", "GZMH", "GNLY", "NKG7", "CCL5", "CCL4"), ncol = 3)
tcr_umap <- FeaturePlot(combo, order = TRUE, features = c("TRAV38-1", "TRBV11-3"))

ggplot2::ggsave("~/postdoc/scRNAseq/combo_analysis_figures/memory_umap.png", memory_umap, height=4, width=12, bg="white")
ggplot2::ggsave("~/postdoc/scRNAseq/combo_analysis_figures/activation_umap.png", activation_umap, height=8, width=12, bg="white")
ggplot2::ggsave("~/postdoc/scRNAseq/combo_analysis_figures/tcr_umap.png", tcr_umap, height=4, width=8, bg="white")




 

n_infection_cluster_umap <- DimPlot(combo, reduction = "umap", split.by = "n_infection")
ggplot2::ggsave("~/postdoc/scRNAseq/combo_analysis_figures/n_infection_combo_cluster_umap.png", n_infection_cluster_umap, height=5, width=12, bg="white")

third <- subset(combo, subset = n_infection=="Third")

third$sample_ID <- factor(third$sample_ID, levels=c("v305_Baseline", "v308_Baseline", "v310_Baseline", "v305_T6", "v308_T6", "v310_T6"))

third_by_volutneer <- DimPlot(third, reduction = "umap", split.by = "sample_ID", ncol = 3)
ggplot2::ggsave("~/postdoc/scRNAseq/combo_analysis_figures/third_by_volutneer_cluster_umap.png", third_by_volutneer, height=10, width=12, bg="white")


# subset on activated guys ####
activated_clusters <- subset(combo,  seurat_clusters %in% c(7, 9, 10, 11, 14))

activated_clusters <- RunPCA(activated_clusters, features = VariableFeatures(object = activated_clusters))

activated_clusters <- FindNeighbors(activated_clusters, dims = 1:20)
activated_clusters <- FindClusters(activated_clusters, resolution = 1, random.seed = 1234, algorithm = 1)

# find markers distinguishing clusters
activated_clusters.markers <- FindAllMarkers(activated_clusters, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)


top20_per_cluster <- data.frame(activated_clusters.markers %>%
                                  group_by(cluster) %>%
                                  slice_min(n = 20, order_by = p_val_adj))

View(top20_per_cluster)


activated_clusters <- RunUMAP(activated_clusters, dims = 1:30, seed.use = 1234)
activated_cluster_umap <- DimPlot(activated_clusters, reduction = "umap", label=TRUE)

ggsave("~/postdoc/scRNAseq/combo_analysis_figures/activated_cluster_umap.png", activated_cluster_umap, height=4, width=6, bg="white", dpi=444)


DimPlot(activated_clusters, reduction = "umap", label=TRUE, split.by = "time_n_infection")
treg_plot <- FeaturePlot(activated_clusters, features = c("IL2RB","IL7R", "FOXP3",  "CD38"), ncol=4)

FeaturePlot()

cluster_counts <- data.frame(table(activated_clusters@active.ident, activated_clusters@meta.data$sample_ID))
colnames(cluster_counts) <- c("Cluster_ID", "Sample_ID", "Count")
cluster_counts$Timepoint <- substr(cluster_counts$Sample_ID, 6, 14)
cluster_counts$Volunteer <- substr(cluster_counts$Sample_ID, 1, 4)
cluster_counts$N_Infection <- ifelse(cluster_counts$Volunteer %in% c("v313", "v315", "v320"), "First", "Third")

activated_cluster_percentages <- ggplot(cluster_counts, aes(x=Cluster_ID, y=Percentage, fill=Timepoint))+
  geom_boxplot()+
  #geom_point(position = position_dodge(width = 0.75))+
  scale_y_continuous(labels = scales::label_percent())+
  scale_fill_manual(values=c("#AA3377", "#4477AA"))+
  facet_wrap(~N_Infection)+
  theme_minimal()

ggsave("~/postdoc/scRNAseq/combo_analysis_figures/activated_cluster_percentages.png", activated_cluster_percentages, height=3, width=7, bg="white", dpi=444)


cluster_5_7_markers <- FindMarkers(activated_clusters, ident.1 = 7, ident.2 = 5, min.pct = 0.25)
head(cluster_5_7_markers, n = 30)

big_diff_plot <- FeaturePlot(activated_clusters, features = rownames(cluster_5_7_markers)[1:36], order = TRUE, ncol=6)
ggsave("~/postdoc/scRNAseq/combo_analysis_figures/big_diff_plot.png", big_diff_plot, height=16, width=21, bg="white", dpi=444)


small_diff_plot <- FeaturePlot(activated_clusters, features = c("IL2RA", "IL7R", "FOXP3", "CD38", "NKG7", "CCL5", "GZMK", "KLRG1", "KLRB1"), order = TRUE, ncol=3)
ggsave("~/postdoc/scRNAseq/combo_analysis_figures/small_diff_plot.png", small_diff_plot, height=8, width=12, bg="white", dpi=444)

small_diff_plot2 <- FeaturePlot(activated_clusters, features = c("TNF", "IFNG", "CCR6", "CCR7", "CXCR3", "CXCR4"), order = TRUE, ncol=3)
ggsave("~/postdoc/scRNAseq/combo_analysis_figures/small_diff_plot2.png", small_diff_plot2, height=8, width=14, bg="white", dpi=444)


# pseudobulk


combo <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/pseudotime_dimred_processed_combo20k.H5Seurat")


first <- subset(combo, subset = n_infection=="First")
third <- subset(combo, subset = n_infection=="Third")

library(Libra)

first_DE <- run_de(first, cell_type_col = "seurat_clusters", label_col = c("timepoint"), replicate_col="volunteer", de_method = "edgeR", de_type="QLF")
third_DE <- run_de(combo, cell_type_col = "seurat_clusters", label_col = c("timepoint"), replicate_col="volunteer", de_method = "edgeR", de_type="QLF")

sig_first_DE <- subset(first_DE, p_val_adj<0.05)
sig_third_DE <- subset(third_DE, p_val_adj<0.05)


combo <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/pseudotime_dimred_processed_combo20k.H5Seurat")
#combo <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/no_jackstraw_demulitplexerd_all_processed.h5Seurat")

combo.sce <- Seurat::as.SingleCellExperiment(combo)

sum_by <- c("seurat_clusters", "timepoint")

summed <- scuttle::aggregateAcrossCells(combo.sce, id=SummarizedExperiment::colData(combo.sce)[,sum_by])


colnames(summed) <- apply(SummarizedExperiment::colData(summed)[,sum_by], 1, function(x)paste(x, collapse="_"))


metadata <- data.frame("column_name"=colnames(SummarizedExperiment::assay(summed, "counts")))
metadata$cluster_ID <- substr(metadata$column_name, 1, 2)
metadata$cluster_ID <- gsub("_", "", metadata$cluster_ID)

metadata$timepoint <- ifelse(grepl("T6", metadata$column_name), "T6", "Baseline")
metadata$n_infection <- ifelse(grepl("First", metadata$column_name), "First", "Third")
metadata$sample_type <- paste(metadata$timepoint, metadata$n_infection, sep="_")

raw_summed <- SummarizedExperiment::assay(summed, "counts")

design <- diffcyt::createDesignMatrix(metadata, c("timepoint", "n_infection"))


out <- scran::pseudoBulkDGE(raw_summed, 
                     label = summed$seurat_clusters,
                     # vector or factor of length equal to ncol(x), specifying the experimental condition for each column
                     col.data = metadata,
                     condition = summed$timepoint,
                     # A formula to be used to construct a design matrix from variables in col.data
                     design = ~timepoint,
                     contrast = matrix(c(0,1))
)

metadata(out[[1]])$design


c14_vs_11  <- Seurat::FindMarkers(combo, ident.1 = 6, min.pct = 0.25)
head(c14_vs_11, n = 30)
