
BiocManager::install(c("edgeR", "DESeq2", "limma"))

if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("neurorestore/Libra")

remotes::install_github("mojaveazure/seurat-disk")

big_data <-SeuratDisk::LoadH5Seurat("~/postdoc/edinburgh/scRNAseq/seurat_objecitons/no_jackstraw_demulitplexerd_all_processed.h5Seurat")
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
