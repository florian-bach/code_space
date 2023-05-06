library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
# 
# first <- LoadH5Seurat("~/Documents/vac63c_scrna_seq/demultiplexed_primary.h5Seurat")
# third <- LoadH5Seurat("~/Documents/vac63c_scrna_seq/demultiplexed_tertiary.h5Seurat")
# 
# 
# combo <- merge(first, third, add.cell.ids = c("First", "Third"))
# 
# combo[["percent.mt"]] <- PercentageFeatureSet(combo, pattern = "^MT-")
# 
# 
# combo <- NormalizeData(combo, normalization.method = "LogNormalize", scale.factor = 10000)
# combo <- FindVariableFeatures(combo, selection.method = "vst", nfeatures = 1000)
# all.genes <- rownames(combo)
# combo <- ScaleData(combo, features = all.genes)
# combo <- RunPCA(combo, features = VariableFeatures(object = combo))
# 
# 
# SeuratDisk::SaveH5Seurat(combo, "~/Documents/vac63c_scrna_seq/no_jackstraw_demulitplexerd_all_processed.h5Seurat")

# data processing ####
combo <- SeuratDisk::LoadH5Seurat("~/Documents/vac63c_scrna_seq/demulitplexerd_all_processed.h5Seurat")
combo <- subset(combo, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)




combo <- AddMetaData(
  object = combo,
  metadata = substr(colnames(combo), 1,5),
  col.name = "n_infection"
)

combo <- AddMetaData(
  object = combo,
  metadata = factor(paste(combo$timepoint, combo$n_infection, sep="_"), levels=c("Baseline_First", "T6_First", "Baseline_Third", "T6_Third")),
  col.name = "time_n_infection"
)




combo <- JackStraw(combo, num.replicate = 100)
combo <- ScoreJackStraw(combo, dims = 1:20)

JackStrawPlot(combo, dims = 1:20)

combo <- FindNeighbors(combo, dims = 1:12)
combo <- FindClusters(combo, resolution = 0.8)
DimPlot(combo, label=TRUE)

cluster_layout_plot <- DimPlot(combo, label=TRUE)+(FeaturePlot(combo, c("CCL5", "CCR7"), order=TRUE)&viridis::scale_color_viridis(option="A"))+plot_layout(widths = c(1, 2))
ggsave("~/Documents/vac63c_scrna_seq/figures/final_cluster_layout_plot.png", cluster_layout_plot, height = 4, width=11, bg="white", dpi=444)

cluster_layout_plot2 <- DimPlot(combo, label=TRUE)+(FeaturePlot(combo, c("GNLY", "HLA-DRB1"), order=TRUE)&viridis::scale_color_viridis(option="A"))+plot_layout(widths = c(1, 2))
ggsave("~/Documents/vac63c_scrna_seq/figures/final_cluster_layout_plot2.png", cluster_layout_plot2, height = 4, width=11, bg="white", dpi=444)

cluster_layout_plot3 <- DimPlot(combo, label=TRUE)+(FeaturePlot(combo, c("FOXP3", "IL7R"), order=TRUE)&viridis::scale_color_viridis(option="A"))+plot_layout(widths = c(1, 2))
ggsave("~/Documents/vac63c_scrna_seq/figures/final_cluster_layout_plot3.png", cluster_layout_plot3, height = 4, width=11, bg="white", dpi=444)


weird_cytotoxic_plot <- FeaturePlot(combo, c("GZMB", "TRAV38-1", "TRBV11-3"), order=TRUE, ncol=3)&viridis::scale_color_viridis(option="A")
ggsave("~/Documents/vac63c_scrna_seq/figures/weird_cytotoxic_plot.png", weird_cytotoxic_plot, height = 4, width=11, bg="white", dpi=444)


cluster_counts <- data.frame(table(combo@active.ident, combo@meta.data$sample_ID))
colnames(cluster_counts) <- c("Cluster_ID", "Sample_ID", "Count")
cluster_counts$Timepoint <- substr(cluster_counts$Sample_ID, 6, 14)
cluster_counts$Volunteer <- substr(cluster_counts$Sample_ID, 1, 4)
cluster_counts$N_Infection <- ifelse(cluster_counts$Volunteer %in% c("v313", "v315", "v320"), "First", "Third")


cluster_counts <- cluster_counts %>%
  group_by(Sample_ID) %>%
  mutate("Percentage"=Count/sum(Count))

# write.csv(cluster_counts, "~/Documents/vac63c_scrna_seq/final_cluster_counts.csv", row.names = FALSE)

all_cluster_counts <- read.csv("~/Documents/vac63c_scrna_seq/final_cluster_counts.csv")

overall_percentages <- ggplot(all_cluster_counts, aes(x=factor(Cluster_ID), y=Percentage, fill=Timepoint))+
  geom_boxplot()+
  #geom_point(position = position_dodge(width = 0.75))+
  scale_y_continuous(labels = scales::label_percent())+
  scale_fill_manual(values=c("#AA3377", "#4477AA"))+
  facet_wrap(~N_Infection)+
  xlab("Cluster_ID")+
  theme_minimal() 

ggsave("~/Documents/vac63c_scrna_seq/figures/overall_percentages.png", overall_percentages, height=3, width=6, dpi=444, bg="white")

# remove cytotoxic cluster, zoom in on memory cells ####
activated_subset <- subset(combo, seurat_clusters %in% c(3, 6, 7, 9, 10))

activated_subset <- RunPCA(activated_subset, features = VariableFeatures(object = activated_subset, selection.method = "vst", nfeatures = 1000))

activated_subset <- JackStraw(activated_subset, num.replicate = 100)
activated_subset <- ScoreJackStraw(activated_subset, dims = 1:20)

JackStrawPlot(activated_subset, dims = 1:20)



activated_subset <- RunUMAP(activated_subset, dims = 1:12)

activated_subset <- FindNeighbors(activated_subset, dims = 1:12)
#0.8
activated_subset <- FindClusters(activated_subset, resolution = 0.8)
DimPlot(activated_subset, label=TRUE)
FeaturePlot(activated_subset, c("CD38", "GNLY", "HLA-DRB1"))



DimPlot(activated_subset, label=TRUE)

FeaturePlot(activated_subset, features=c("CD38", "CCL5", "TNF", "IFNG",  "GZMK", "NKG7"))

# SeuratDisk::SaveH5Seurat(activated_subset, "~/Documents/vac63c_scrna_seq/final_activated_subset.h5Seurat")
# SeuratDisk::SaveH5Seurat(combo, "~/Documents/vac63c_scrna_seq/final_all_combo.h5Seurat")


activated_cluster_layout_plot <- DimPlot(activated_subset, label=TRUE)+(FeaturePlot(activated_subset, c("IL7R", "KLRB1"), order=TRUE)&viridis::scale_color_viridis(option="A"))+plot_layout(widths = c(1, 2))
ggsave("~/Documents/vac63c_scrna_seq/figures/activated_cluster_layout_plot.png", activated_cluster_layout_plot, height = 4, width=11, bg="white", dpi=444)

# activated_cluster_layout_plot2 <- DimPlot(activated_subset, label=TRUE)+(FeaturePlot(activated_subset, c("TCF7", "KLRB1"), order=TRUE)&viridis::scale_color_viridis(option="A"))+plot_layout(widths = c(1, 2))
# ggsave("~/Documents/vac63c_scrna_seq/figures/activated_cluster_layout_plot.png", activated_cluster_layout_plot2, height = 4, width=11, bg="white", dpi=444)
# 

activated_cluster_effector_plot <- FeaturePlot(activated_subset, c("CCL5", "GNLY", "HLA-DRB1", "GZMK"), ncol = 4, order=TRUE)&viridis::scale_color_viridis(option="A")
ggsave("~/Documents/vac63c_scrna_seq/figures/activated_cluster_effector_plot.png", activated_cluster_effector_plot, height = 4, width=11, bg="white", dpi=444)

activated_cluster_effector_plot2 <- FeaturePlot(activated_subset, c("CD38", "TCF7", "PDCD1", "GZMK"), ncol = 4, order=TRUE)&viridis::scale_color_viridis(option="A")
ggsave("~/Documents/vac63c_scrna_seq/figures/activated_cluster_effector_plot2.png", activated_cluster_effector_plot2, height = 4, width=11, bg="white", dpi=444)


cluster_counts <- data.frame(table(activated_subset@active.ident, activated_subset@meta.data$sample_ID))
colnames(cluster_counts) <- c("Cluster_ID", "Sample_ID", "Count")
cluster_counts$Timepoint <- substr(cluster_counts$Sample_ID, 6, 14)
cluster_counts$Volunteer <- substr(cluster_counts$Sample_ID, 1, 4)
cluster_counts$N_Infection <- ifelse(cluster_counts$Volunteer %in% c("v313", "v315", "v320"), "First", "Third")


cluster_counts <- cluster_counts %>%
  group_by(Sample_ID) %>%
  mutate("Percentage"=Count/sum(Count))


write.csv(cluster_counts, "~/Documents/vac63c_scrna_seq/final_activated_cluster_counts.csv", row.names = FALSE)
activated_cluster_counts <- read.csv("~/Documents/vac63c_scrna_seq/final_activated_cluster_counts.csv")

activated_cluster_percentages <- ggplot(activated_cluster_counts, aes(x=factor(Cluster_ID), y=Percentage, fill=Timepoint))+
  geom_boxplot()+
  #geom_point(position = position_dodge(width = 0.75))+
  scale_y_continuous(labels = scales::label_percent())+
  scale_fill_manual(values=c("#AA3377", "#4477AA"))+
  facet_wrap(~N_Infection)+
  theme_minimal()

ggsave("~/Documents/vac63c_scrna_seq/figures/activated_cluster_percentages.png", activated_cluster_percentages, height=3, width=6, dpi=444, bg="white")

activated_subset_markers <- FindAllMarkers(activated_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

activated_subset_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)


# differential gene expression ####


library(Libra)


first <- subset(activated_subset, subset = n_infection=="First")
third <- subset(activated_subset, subset = n_infection=="Third")


first_DE <- run_de(first, cell_type_col = "seurat_clusters", label_col = c("timepoint"), replicate_col="volunteer", de_method = "edgeR", de_type = "LRT")
sig_first_DE <- subset(first_DE, p_val_adj<0.1)
dim(sig_first_DE)
length(unique(sig_first_DE$gene))
#lrt: 155 genes, 169 gene-cluster_combos


`%notin%` <- Negate(`%in%`)
third_DE <- run_de(third, cell_type_col = "seurat_clusters", label_col = c("timepoint"), replicate_col="volunteer", de_method = "edgeR", de_type="LRT")
sig_third_DE <- subset(third_DE, p_val_adj<0.1)
dim(sig_third_DE)
length(unique(sig_third_DE$gene))
#lrt: 65 genes, 78 gene-cluster_combos

sig_first_de_big_plot <- FeaturePlot(activated_subset, features = c("GNLY", "HOPX", "BCL2A1", "COMT", "GLG1", "UBE2A"), split.by = "time_n_infection", order=TRUE)&viridis::scale_color_viridis(option="A", direction=-1)
ggsave("~/Documents/vac63c_scrna_seq/figures/sig_first_de_big_plot.png", sig_first_de_big_plot, height=13, width=13, bg="white", dpi=444)

sig_first_de_big_plot2 <- FeaturePlot(activated_subset, features = c("PSMC1", "GBP1", "LINC01871", "REST", "STAT1", "TM2D1"), split.by = "time_n_infection", order=TRUE)&viridis::scale_color_viridis(option="A", direction=-1)
ggsave("~/Documents/vac63c_scrna_seq/figures/sig_first_de_big_plot2.png", sig_first_de_big_plot2, height=13, width=13, bg="white", dpi=444)

sig_first_isgs <- FeaturePlot(activated_subset, features = c("STAT1", "NFKB1", "GBP1", "GBP2", "GBP5", "IFNG", "TBX21"), split.by = "time_n_infection", order=TRUE)&viridis::scale_color_viridis(option="A", direction=-1)
ggsave("~/Documents/vac63c_scrna_seq/figures/sig_first_isgs.png", sig_first_isgs, height=13, width=13, bg="white", dpi=444)

exhaustion_genes <- FeaturePlot(activated_subset, features = c("TCF7", "TOX", "TOX2", "HAVCR2", "PDCD1", "LAG3", "TIGIT"), split.by = "time_n_infection", order=TRUE)&viridis::scale_color_viridis(option="A", direction=-1)
ggsave("~/Documents/vac63c_scrna_seq/figures/exhaustion.png", exhaustion_genes, height=13, width=13, bg="white", dpi=444)

cytotoxicity <- FeaturePlot(activated_subset, features = c("GNLY", "NKG7", "PRF1", "GZMH", "GZMK"), split.by = "time_n_infection", order=TRUE)&viridis::scale_color_viridis(option="A", direction=-1)
ggsave("~/Documents/vac63c_scrna_seq/figures/cytotoxicity.png", cytotoxicity, height=13, width=13, bg="white", dpi=444)


cyto_chemo <- FeaturePlot(activated_subset, features = c("CCL5", "CCL4", "TNF", "IFNG", "IL1B", "IL21"), split.by = "time_n_infection", order=TRUE)&viridis::scale_color_viridis(option="A", direction=-1)
ggsave("~/Documents/vac63c_scrna_seq/figures/cyto_chemo.png", cyto_chemo, height=13, width=13, bg="white", dpi=444)

diff_tfs <- FeaturePlot(activated_subset, features = c("MKI67", "RORC", "EGR1", "EGR2", "EOMES", "IKZF2"), split.by = "time_n_infection", order=TRUE)&viridis::scale_color_viridis(option="A", direction=-1)
ggsave("~/Documents/vac63c_scrna_seq/figures/diff_tfs.png", diff_tfs, height=13, width=13, bg="white", dpi=444)

chemo_recep <- FeaturePlot(activated_subset, features = c("CXCR4", "CXCR3", "CXCR5", "CCR4", "CCR5", "CX3CR1"), split.by = "time_n_infection", order=TRUE)&viridis::scale_color_viridis(option="A", direction=-1)
ggsave("~/Documents/vac63c_scrna_seq/figures/chemo_recep.png", chemo_recep, height=13, width=13, bg="white", dpi=444)



sig_third_de_big_plot <- FeaturePlot(activated_subset, features = c("CCL4", "DERA", "SNAI3", "TRIM47", "CWC22", "NKG7"), split.by = "time_n_infection", order=TRUE)&viridis::scale_color_viridis(option="A", direction=-1)
ggsave("~/Documents/vac63c_scrna_seq/figures/sig_third_de_big_plot.png", sig_third_de_big_plot, height=13, width=13, bg="white", dpi=444)

# write.csv(first_DE, "~/Documents/vac63c_scrna_seq/first_DE.csv", row.names = FALSE)
# write.csv(third_DE, "~/Documents/vac63c_scrna_seq/third_DE.csv", row.names = FALSE)
# 
# write.table(unique(sig_first_DE$gene), file = "~/Documents/vac63c_scrna_seq/libra_first_t6.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(unique(sig_third_DE$gene), file = "~/Documents/vac63c_scrna_seq/libra_third_t6.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
#


# milo ####


# milo ####


#write.csv(metadata, "~/postdoc/scRNAseq/metadata/sce_metadata.csv", row.names = TRUE, col.names = TRUE)

library(miloR)
library(patchwork)
library(ggplot2)
library(SingleCellExperiment)

combo.sce <- Seurat::as.SingleCellExperiment(activated_subset)

combo_milo <- Milo(combo.sce)
combo_milo <- buildGraph(combo_milo, k = 30, d = 30, reduced.dim = "PCA")
combo_milo <- makeNhoods(combo_milo, prop = 0.1, k = 30, d=30, refined = TRUE, reduced_dims = "PCA")

# need this for colData & owData

metadata <- read.csv("~/Documents/vac63c_scrna_seq/sce_metadata.csv", row.names = 1)

combo_milo <- countCells(combo_milo, meta.data = as.data.frame(colData(combo_milo)), sample="sample_ID")

combo_milo <- calcNhoodDistance(combo_milo, d=30, reduced.dim = "PCA")

da_results <- testNhoods(combo_milo, design = ~ timepoint + volunteer, design.df = metadata)

combo_milo <- buildNhoodGraph(combo_milo)



## Plot single-cell UMAP
umap_pl <- scater::plotReducedDim(combo_milo, dimred = "UMAP", colour_by="timepoint",
                                  text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(combo_milo, da_results, layout="UMAP",alpha=0.1)

milo_umap <- umap_pl + nh_graph_pl +
  plot_layout(guides="collect")

ggsave("~/Documents/vac63c_scrna_seq/figures/milo_umap.png", milo_umap, height=8, width=8, bg="white", dpi=444)

# add cluster id to network

da_results <- annotateNhoods(combo_milo, da_results, coldata_col = "seurat_clusters")

milo_bee <- plotDAbeeswarm(da_results, group.by = "seurat_clusters")
ggsave("~/Documents/vac63c_scrna_seq/figures/milo_bee.png", milo_bee, height=8, width=8, bg="white", dpi=444)


# establish key marker genes for each neighbourhood
combo_milo <- scuttle::logNormCounts(combo_milo)

da_results$NhoodGroup <- as.numeric(da_results$SpatialFDR < 0.1 & da_results$logFC < 0)
da_nhood_markers <- findNhoodGroupMarkers(combo_milo, da_results, subset.row = Seurat::VariableFeatures(activated_subset), 
                                          aggregate.samples = TRUE, sample_col = "sample_ID")


## Warning: Zero sample variances detected, have been offset away from zero

## Warning: Zero sample variances detected, have been offset away from zero

write.csv(da_results, "~/Documents/vac63c_scrna_seq/milo_de_results.csv", row.names = FALSE)


# increase granularity

granular_activation <- activated_subset

granular_activation <- FindNeighbors(granular_activation, dims = 1:12)
#0.8
granular_activation <- FindClusters(granular_activation, resolution = 2.5)
DimPlot(granular_activation, label=TRUE)





cluster_counts <- data.frame(table(granular_activation@active.ident, granular_activation@meta.data$sample_ID))
colnames(cluster_counts) <- c("Cluster_ID", "Sample_ID", "Count")
cluster_counts$Timepoint <- substr(cluster_counts$Sample_ID, 6, 14)
cluster_counts$Volunteer <- substr(cluster_counts$Sample_ID, 1, 4)
cluster_counts$N_Infection <- ifelse(cluster_counts$Volunteer %in% c("v313", "v315", "v320"), "First", "Third")


cluster_counts <- cluster_counts %>%
  group_by(Sample_ID) %>%
  mutate("Percentage"=Count/sum(Count))


write.csv(cluster_counts, "~/Documents/vac63c_scrna_seq/granular_activated_cluster_counts.csv", row.names = FALSE)
activated_cluster_counts <- read.csv("~/Documents/vac63c_scrna_seq/granular_activated_cluster_counts.csv")

activated_cluster_percentages <- ggplot(activated_cluster_counts, aes(x=factor(Cluster_ID), y=Percentage, fill=Timepoint))+
  geom_boxplot()+
  #geom_point(position = position_dodge(width = 0.75))+
  scale_y_continuous(labels = scales::label_percent())+
  scale_fill_manual(values=c("#AA3377", "#4477AA"))+
  facet_wrap(~N_Infection)+
  xlab("Cluster_ID")+
  theme_minimal()

ggsave("~/Documents/vac63c_scrna_seq/figures/granular_activated_cluster_percentages.png", activated_cluster_percentages, height=3, width=12, dpi=444, bg="white")


# pseudotime ####

combo <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/final_activated_subset.h5Seurat")
library(monocle3)
combo.cds <- SeuratWrappers::as.cell_data_set(combo)
#k=8 is nice
combo.cds <- monocle3::cluster_cells(cds=combo.cds, reduction_method = "UMAP", k=8, random_seed = 1234)
combo.cds <- monocle3::learn_graph(combo.cds, use_partition = TRUE, verbose=TRUE)

combo.cds <- monocle3::order_cells(combo.cds, reduction_method = "UMAP")


monocle3::plot_cells(group_cells_by = "cluster",
                     cds = combo.cds, 
                     color_cells_by = "pseudotime",
                     show_trajectory_graph = TRUE
)


ggplot2::ggsave("~/postdoc/scRNAseq/final_analysis_and_figures/figures/pseudotime_k8.png")
# granulated clustering ####


activated_clusters <- Seurat::FindClusters(activated_clusters, resolution = 2.2)
granular_cluster_umap <- Seurat::DimPlot(activated_clusters, label=TRUE)
ggsave("~/postdoc/scRNAseq/final_analysis_and_figures/figures/granular_cluster_umap.png", granular_cluster_umap, height=4, width=6, dpi=444, bg="white")


cluster_counts <- data.frame(table(activated_clusters@active.ident, activated_clusters@meta.data$sample_ID))
colnames(cluster_counts) <- c("Cluster_ID", "Sample_ID", "Count")
cluster_counts$Timepoint <- substr(cluster_counts$Sample_ID, 6, 14)
cluster_counts$Volunteer <- substr(cluster_counts$Sample_ID, 1, 4)
cluster_counts$N_Infection <- ifelse(cluster_counts$Volunteer %in% c("v313", "v315", "v320"), "First", "Third")


cluster_counts <- cluster_counts %>%
  group_by(Sample_ID) %>%
  mutate("Percentage"=Count/sum(Count))


write.csv(cluster_counts, "~/postdoc/scRNAseq/final_analysis_and_figures/granular_activated_cluster_counts.csv", row.names = FALSE)

# differential gene expresion ####


first <- subset(combo, subset = n_infection=="First")
third <- subset(combo, subset = n_infection=="Third")


library(Libra)

first_DE <- run_de(first, cell_type_col = "seurat_clusters", label_col = c("timepoint"), replicate_col="volunteer", de_method = "edgeR", de_type = "LRT", min_cells=10)
sig_first_DE <- subset(first_DE, p_val_adj<0.1& abs(avg_logFC) > log2(1.5))
dim(sig_first_DE)
# lrt 169 gene-clsuter combinations; 155 genes

`%notin%` <- Negate(`%in%`)
third_DE <- run_de(third, cell_type_col = "seurat_clusters", label_col = c("timepoint"), replicate_col="volunteer", de_method = "edgeR", de_type="LRT")
sig_third_DE <- subset(third_DE, p_val_adj<0.1 & abs(avg_logFC) > log2(1.5))
dim(sig_third_DE)
# lrt 78 gene-cluster combinations, 65 genes


  