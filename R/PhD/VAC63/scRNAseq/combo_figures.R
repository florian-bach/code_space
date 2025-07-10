library(Seurat)


#combo <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/pseudotime_dimred_processed_combo20k.H5Seurat")
combo <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/final_analysis_and_figures/final_activated_subset.h5Seurat")

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
combo <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/pseudotime_dimred_processed_combo20k.H5Seurat")

library(Seurat)

activated_clusters <- subset(combo,  seurat_clusters %in% c(6:9, 10:12, 14))
combo <- NULL
activated_clusters <- RunPCA(activated_clusters, features = VariableFeatures(object = activated_clusters))

activated_clusters <- FindNeighbors(activated_clusters, dims = 1:20)
activated_clusters <- FindClusters(activated_clusters, resolution = 2, random.seed = 1234, algorithm = 1)

# find markers distinguishing clusters
# activated_clusters.markers <- FindAllMarkers(activated_clusters, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
# 
# 
# top20_per_cluster <- data.frame(activated_clusters.markers %>%
#                                   group_by(cluster) %>%
#                                   slice_min(n = 20, order_by = p_val_adj))
# 
# View(top20_per_cluster)
# 

activated_clusters <- RunUMAP(activated_clusters, dims = 1:30, seed.use = 1234)
activated_cluster_umap <- DimPlot(activated_clusters, reduction = "umap", label=TRUE)

ggsave("~/postdoc/scRNAseq/combo_analysis_figures/activated_cluster_umap.png", activated_cluster_umap, height=4, width=6, bg="white", dpi=444)


library(ggplot2)



DimPlot(activated_clusters, reduction = "umap", label=TRUE, split.by = "time_n_infection")
treg_plot <- FeaturePlot(activated_clusters, features = c("IL2RB","IL7R", "FOXP3",  "CD38"), ncol=4)



cluster_counts <- data.frame(table(combo@active.ident, combo@meta.data$sample_ID))
colnames(cluster_counts) <- c("Cluster_ID", "Sample_ID", "Count")
cluster_counts$Timepoint <- substr(cluster_counts$Sample_ID, 6, 14)
cluster_counts$Volunteer <- substr(cluster_counts$Sample_ID, 1, 4)
cluster_counts$N_Infection <- ifelse(cluster_counts$Volunteer %in% c("v313", "v315", "v320"), "First", "Third")




cluster_counts %>%
  group_by(Cluster_ID, Timepoint, N_Infection) %>%
  summarise(Count)


cluster_counts <- cluster_counts %>%
  group_by(Sample_ID) %>%
  mutate("Percentage"=Count/sum(Count))



activated_cluster_percentages <- ggplot(cluster_counts, aes(x=Cluster_ID, y=Count, fill=Timepoint))+
  geom_boxplot()+
  #geom_point(position = position_dodge(width = 0.75))+
  #scale_y_continuous(labels = scales::label_percent())+
  scale_fill_manual(values=c("#AA3377", "#4477AA"))+
  facet_wrap(~N_Infection, ncol=1)+
  theme_minimal()+
  theme(strip.text = element_text(hjust = 0))

ggsave("~/postdoc/scRNAseq/combo_analysis_figures/activated_cluster_percentages.png", activated_cluster_percentages, height=3, width=7, bg="white", dpi=444)


cluster_5_7_markers <- FindMarkers(activated_clusters, ident.1 = 7, ident.2 = 5, min.pct = 0.25)
head(cluster_5_7_markers, n = 30)

big_diff_plot <- FeaturePlot(activated_clusters, features = rownames(cluster_5_7_markers)[1:36], order = TRUE, ncol=6)
ggsave("~/postdoc/scRNAseq/combo_analysis_figures/big_diff_plot.png", big_diff_plot, height=16, width=21, bg="white", dpi=444)


small_diff_plot <- FeaturePlot(activated_clusters, features = c("IL2RA", "IL7R", "FOXP3", "CD38", "NKG7", "CCL5", "GZMK", "KLRG1", "KLRB1"), order = TRUE, ncol=3)
ggsave("~/postdoc/scRNAseq/combo_analysis_figures/small_diff_plot.png", small_diff_plot, height=8, width=12, bg="white", dpi=444)

small_diff_plot2 <- FeaturePlot(activated_clusters, features = c("TNF", "IFNG", "CCR6", "CCR7", "CXCR3", "CXCR4"), order = TRUE, ncol=3)
ggsave("~/postdoc/scRNAseq/combo_analysis_figures/small_diff_plot2.png", small_diff_plot2, height=8, width=14, bg="white", dpi=444)


# pseudobulk ####


combo <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/pseudotime_dimred_processed_combo20k.H5Seurat")


# first <- subset(combo, subset = n_infection=="First")
# third <- subset(combo, subset = n_infection=="Third")


activated_clusters <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/final_analysis_and_figures/final_activated_subset.h5Seurat")


first <- subset(activated_clusters, subset = n_infection=="First")
third <- subset(activated_clusters, subset = n_infection=="Third")


library(Libra)

first_DE <- run_de(first, cell_type_col = "seurat_clusters", label_col = c("timepoint"), replicate_col="volunteer", de_method = "edgeR", de_type = "LRT")
sig_first_DE <- subset(first_DE, p_val_adj<0.1)
sig_first_DE$avg_logFC <- sig_first_DE$avg_logFC*-1
write.csv(sig_first_DE, "~/postdoc/scRNAseq/final_analysis_and_figures/sig_first_DE.csv", quote = FALSE, row.names = FALSE)

#qlf: 94 genes, 99 gene-cluster_combos
#lrt: 182 genes, 190 gene-cluster_combos


`%notin%` <- Negate(`%in%`)
third_DE <- run_de(third, cell_type_col = "seurat_clusters", label_col = c("timepoint"), replicate_col="volunteer", de_method = "edgeR", de_type="LRT")
sig_third_DE <- subset(third_DE, p_val_adj<0.1)
sig_third_DE$avg_logFC <- sig_third_DE$avg_logFC*-1
write.csv(sig_third_DE, "~/postdoc/scRNAseq/final_analysis_and_figures/sig_third_DE.csv", quote = FALSE, row.names = FALSE)

#qlf: 11 genes, 13 9gene-cluster_combos
#lrt: 27 genes, 30 gene-cluster_combos


big_plot <- Seurat::FeaturePlot(activated_clusters, c("CD38", "KLRB1", "GNLY", "HLA-DRB1", "GZMK"), split.by="time_n_infection", order=TRUE)&viridis::scale_color_viridis(option="A", direction=-1)
ggplot2::ggsave("~/postdoc/scRNAseq/combo_analysis_figures/activation_n_infection_plot.png", height=16, width=18, bg="white", dpi=444)

# 301
base_DE <- Libra::run_de(subset(combo, subset = timepoint =="Baseline"), cell_type_col = "seurat_clusters", label_col = c("n_infection"), replicate_col="volunteer", de_method = "edgeR", de_type="LRT")
sig_base_DE <- subset(base_DE, p_val_adj<0.1 & abs(avg_logFC)>log2(1.5))



library(SingleCellExperiment)
activated_clusters <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/final_activated_subset.h5Seurat")
# activated_clusters <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/no_jackstraw_demulitplexerd_all_processed.h5Seurat")
# combo <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/activated_subset.h5Seurat")

#activated_clusters <- subset(combo,  seurat_clusters %in% c(6:9, 10:12))

all.genes <- rownames(activated_clusters)


combo.sce <- Seurat::as.SingleCellExperiment(activated_clusters)


sum_by <- c("sample_ID")


summed <- scuttle::aggregateAcrossCells(combo.sce, id=SummarizedExperiment::colData(combo.sce)[,sum_by])
colnames(summed) <- apply(SummarizedExperiment::colData(summed)[,sum_by], 1, function(x)paste(x, collapse="_"))



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
# 
# metadata <- data.frame("column_name"=colnames(SummarizedExperiment::assay(summed, "counts")))
# metadata$cluster_ID <- substr(metadata$column_name, 1, 2)
# metadata$cluster_ID <- gsub("_", "", metadata$cluster_ID)
# 
# metadata$timepoint <- ifelse(grepl("T6", metadata$column_name), "T6", "Baseline")
# 
# metadata$volunteer <- substr(metadata$column_name, regexpr("v", metadata$column_name), regexpr("v", metadata$column_name)+3)
# 
# metadata$n_infection <- ifelse(metadata$volunteer %in% c("v313", "v315", "v320"), "First", "Third")
# 
# metadata$sample_type <- paste(metadata$timepoint, metadata$n_infection, sep="_")
# metadata$sample_ID <- paste(metadata$volunteer, metadata$timepoint, sep="_")
metadata <- read.csv("~/postdoc/scRNAseq/metadata/sce_metadata.csv", row.names = 1)

design <- createDesignMatrix(ei, c("timepoint", "volunteer"))

# raw_summed <- SummarizedExperiment::assay(summed, "counts")

# design <- diffcyt::createDesignMatrix(metadata, c("timepoint", "volunteer"))


out <- scran::pseudoBulkDGE(summed, 
                     label = summed$cell_type,
                     # vector or factor of length equal to ncol(x), specifying the experimental condition for each column
                     col.data = metadata,
                     condition = summed$timepoint,
                     # A formula to be used to construct a design matrix from variables in col.data
                     design = ~timepoint+volunteer,
                     # sig at t6 third infection
                     contrast = matrix(c(0,1,0,0,0,0,0))
)

big_table <- data.frame(matrix(ncol=7))
colnames(big_table) <- c(colnames(out[[1]]), "cluster_ID", "gene")

for(i in 1:length(unique(summed$seurat_clusters))){
  new_entry <- data.frame(out[[i]])
  new_entry$cluster_ID <- paste("cluster_", i-1, sep='')
  new_entry$gene <- rownames(new_entry)
  big_table <- rbind(big_table, new_entry)
}

big_table <- data.frame(out[[1]])
big_table <- na.omit(big_table)

sig_t6 <- subset(big_table, FDR<0.10)


big_plot <- Seurat::FeaturePlot(combo, unique(sig_t6[order(sig_t6$FDR),7])[1:16], split.by = "timepoint", order=TRUE, ncol = 4)&viridis::scale_color_viridis(option="B")

ggplot2::ggsave(filename = "~/postdoc/scRNAseq/combo_analysis_figures/scran_DE_t6.png", big_plot, width = 20, height = 40, bg="white")
# write.csv(big_table, "~/postdoc/scRNAseq/differential_gene_expression/scran_base_t6_third.csv", row.names = FALSE)

write.table(unique(sig_first_DE$gene), file = "~/postdoc/scRNAseq/differential_gene_expression/libra_first_t6.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

out2 <- scran::pseudoBulkDGE(summed, 
                            label = summed$seurat_clusters,
                            # vector or factor of length equal to ncol(x), specifying the experimental condition for each column
                            col.data = metadata,
                            condition = summed$timepoint,
                            # A formula to be used to construct a design matrix from variables in col.data
                            design = ~timepoint+n_infection,
                            # sig at t6 in first infection
                            contrast = matrix(c(0,1,0))
)

big_table2 <- data.frame(matrix(ncol=7))
colnames(big_table2) <- c(colnames(out2[[1]]), "cluster_ID", "gene")

# write.csv(big_table2, "~/postdoc/scRNAseq/differential_gene_expression/scran_base_t6_first.csv", row.names = FALSE)

for(i in 1:length(unique(summed$seurat_clusters))){
  new_entry <- data.frame(out2[[i]])
  new_entry$cluster_ID <- paste("cluster_", i-1, sep='')
  new_entry$gene <- rownames(new_entry)
  big_table2 <- rbind(big_table2, new_entry)
}

sig_first <- subset(big_table2, FDR<0.1 & abs(logFC)>log2(1.5))







out3 <- scran::pseudoBulkDGE(summed, 
                             label = summed$seurat_clusters,
                             # vector or factor of length equal to ncol(x), specifying the experimental condition for each column
                             col.data = metadata,
                             condition = summed$timepoint,
                             # A formula to be used to construct a design matrix from variables in col.data
                             design = ~timepoint+n_infection,
                             # sig in first infection
                             # intercept, t6, third
                             contrast = matrix(c(0,0,1))
)

big_table3 <- data.frame(matrix(ncol=7))
colnames(big_table3) <- c(colnames(out3[[1]]), "cluster_ID", "gene")

for(i in 1:length(unique(summed$seurat_clusters))){
  new_entry <- data.frame(out3[[i]])
  new_entry$cluster_ID <- paste("cluster_", i-1, sep='')
  new_entry$gene <- rownames(new_entry)
  big_table3 <- rbind(big_table3, new_entry)
}

sig_base <- subset(big_table3, FDR<0.01 & abs(logFC)>log2(1.5))

# 1,0,-1
top_sig_base <- rbind(dplyr::top_n(sig_base, 10, wt = logFC),
                      dplyr::top_n(sig_base, -10, wt = logFC))

#1,0,0
top_sig_base2 <- rbind(dplyr::top_n(sig_base, 10, wt = logFC),
                      dplyr::top_n(sig_base, -10, wt = logFC))


# drilling down on activation 

cd38_umap <- FeaturePlot(activated_clusters, order = TRUE, features=c("CD38"))&viridis::scale_color_viridis(option="B")
ggsave("~/postdoc/scRNAseq/combo_analysis_figures/cd38_umap.png", cd38_umap, height=4, width=6, bg="white", dpi=444 )


big_pheno_umap <- FeaturePlot(combo, order=T, features=c("LEF1", "TCF7", "CD27", "CD28", "PDCD1", "HAVCR2", "LAG3", "TOX2", "NR4A1"))&viridis::scale_color_viridis(option="B")
ggsave("~/postdoc/scRNAseq/combo_analysis_figures/big_pheno_umap.png", big_pheno_umap, height=8, width=12, bg="white", dpi=444 )

cytotoxic_umap <- FeaturePlot(combo, order=T, features=c("KLRG1", "KLRB1", "PRF1", "GZMB", "GZMK", "GZMH", "FAS", "FOXP3", "IL2RA"))&viridis::scale_color_viridis(option="B")
ggsave("~/postdoc/scRNAseq/combo_analysis_figures/cytotoxic_umap.png", cytotoxic_umap, height=8, width=12, bg="white", dpi=444 )

functional_umap <- FeaturePlot(combo, order=T, features=c("CCL5", "CCL4", "IFNG", "TNF", "IL32", "GZMH", "FAS", "FOXP3", "IL2RA"))&viridis::scale_color_viridis(option="B")
ggsave("~/postdoc/scRNAseq/combo_analysis_figures/functional_umap.png", functional_umap, height=8, width=12, bg="white", dpi=444 )

functional_umap2 <- FeaturePlot(combo, order=T, features=c("TIGIT", "HLA-DRA", "HLA-DRB1", "CST7", "CSF1","APOBEC3H"))&viridis::scale_color_viridis(option="B")
ggsave("~/postdoc/scRNAseq/combo_analysis_figures/functional_umap2.png", functional_umap, height=8, width=12, bg="white", dpi=444 )



activated_clusters.markers <- FindAllMarkers(combo, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)



top10_per_cluster <- data.frame(activated_clusters.markers %>%
                                  group_by(cluster) %>%
                                  slice_min(n = 10, order_by = p_val_adj))




# pseudobulk on whole sample ####

library(SingleCellExperiment)
activated_clusters <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/final_activated_subset.h5Seurat")
# activated_clusters <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/no_jackstraw_demulitplexerd_all_processed.h5Seurat")
# combo <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/activated_subset.h5Seurat")

#activated_clusters <- subset(combo,  seurat_clusters %in% c(6:9, 10:12))

all.genes <- rownames(activated_clusters)


combo.sce <- Seurat::as.SingleCellExperiment(activated_clusters)


sum_by <- c("sample_ID")


summed <- scuttle::aggregateAcrossCells(combo.sce, id=SummarizedExperiment::colData(combo.sce)[,sum_by])
colnames(summed) <- apply(SummarizedExperiment::colData(summed)[,sum_by], 1, function(x)paste(x, collapse="_"))



metadata <- read.csv("~/postdoc/scRNAseq/metadata/sce_metadata.csv", row.names = 1)

design <- createDesignMatrix(ei, c("timepoint", "volunteer"))

summed$cell_type <- "CD4"

first_summed <- CATALYST::filterSCE(summed, n_infection=="First")

out <- scran::pseudoBulkDGE(first_summed, 
                            label = first_summed$cell_type,
                            # vector or factor of length equal to ncol(x), specifying the experimental condition for each column
                            col.data = subset(metadata, n_infection=="First"),
                            condition = first_summed$timepoint,
                            # A formula to be used to construct a design matrix from variables in col.data
                            design = ~timepoint+volunteer,
                            # sig at t6 third infection
                            contrast = matrix(c(0,1,0,0))
)

first_table <- data.frame(out[[1]])
first_table <- na.omit(first_table)

sig_first <- subset(first_table, FDR<0.1)

sig_first_up <- subset(sig_first, logFC>0)
sig_first_down <- subset(sig_first, logFC<0)


write.table(rownames(sig_first_up), "~/postdoc/scRNAseq/final_analysis_and_figures/scran_de_first_up_genes.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(rownames(sig_first_down), "~/postdoc/scRNAseq/final_analysis_and_figures/scran_de_first_down_genes.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)



third_summed <- CATALYST::filterSCE(summed, n_infection=="Third")


out2 <- scran::pseudoBulkDGE(third_summed, 
                            label = third_summed$cell_type,
                            # vector or factor of length equal to ncol(x), specifying the experimental condition for each column
                            col.data = subset(metadata, n_infection=="Third"),
                            condition = third_summed$timepoint,
                            # A formula to be used to construct a design matrix from variables in col.data
                            design = ~timepoint+volunteer,
                            # sig at t6 third infection
                            contrast = matrix(c(0,1,0,0))
                            #coef="timepointT6"
)
metadata(out2[[1]])$design

third_table <- data.frame(out2[[1]])
third_table <- na.omit(third_table)

sig_third <- subset(third_table, FDR<0.20)

# wiebke RNAseq

first_memory_bulk <- data.table::fread("~/PhD/RNAseq/vac63c/04May2021_007_Memory_Post_Primary-Memory_Pre_Primary_genelist.csv")
third_memory_bulk <- data.table::fread("~/PhD/RNAseq/vac63c/04May2021_009_Memory_Post_Tertiary-Memory_Pre_Tertiary_genelist.csv")


sig_first_memory_bulk <- subset(first_memory_bulk, abs(log2FoldChange)>log2(2)&padj<0.1)
sig_third_memory_bulk <- subset(third_memory_bulk, abs(log2FoldChange)>log2(2)&padj<0.1)


