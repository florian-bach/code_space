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


#glms ####


cluster_counts <- read.csv("~/postdoc/scRNAseq/final_analysis_and_figures/granular_activated_cluster_counts.csv")


# granuluar first ####
cluster_counts <- read.csv("~/postdoc/scRNAseq/final_analysis_and_figures/granular_activated_cluster_counts.csv")

granular_first_cluster_counts <- subset(cluster_counts, N_Infection=="First")

granular_first_n_cells <- granular_first_cluster_counts %>%
  group_by(Sample_ID)%>%
  summarise("sum"=sum(Count))

granular_first_list_of_dfs_for_glm <- split(granular_first_cluster_counts, granular_first_cluster_counts$Cluster_ID)



granular_first_list_of_models <- lapply(granular_first_list_of_dfs_for_glm, function(x) glm(Percentage~Timepoint+Volunteer, family = "binomial", weights = granular_first_n_cells$sum, data=x))
#list_of_models <- lapply(list_of_dfs_for_glm, function(x) nlme::lme(concentration~timepoint, random=~1|Volunteer, data=x))
#list_of_models <- lapply(list_of_dfs_for_glm, function(x) lme4::lmer(concentration~timepoint+(1|Volunteer), data=x))

t6_contrast <- t(matrix(c(0,1,0,0)))

# c45_contrast <- t(diffcyt::createContrast(c(0,1,0,0,0,0,0,0)))
##
granular_first_list_of_tests <- lapply(granular_first_list_of_models, function(x) multcomp::glht(x, t6_contrast))
granular_first_list_of_pvalues <- lapply(granular_first_list_of_tests, function(x) data.frame("p_raw"=summary(x)$test$pvalues, "coef"=summary(x)$test$coefficients))

granular_first_df_first_ <- do.call(rbind, granular_first_list_of_pvalues)

granular_first_df_first_$Cluster_ID <- names(granular_first_list_of_pvalues)
#
granular_first_df_first_$p_adj <- p.adjust(granular_first_df_first_$p_raw, method = "fdr")


#log(2) = 0.6931472
subset(granular_first_df_first_, p_adj<=0.05 & abs(coef)>0.6931472)

sig_granular_first <- subset(granular_first_df_first_, p_adj<=0.05 & abs(coef)>0.6931472)
sig_granular_first <- data.frame("cluster"=sig_granular_first$Cluster_ID, "raw_FC"=exp(sig_granular_first$coef), "p_adj"=sig_granular_first$p_adj)

# granular third ####

cluster_counts <- read.csv("~/postdoc/scRNAseq/final_analysis_and_figures/granular_activated_cluster_counts.csv")

granular_third_cluster_counts <- subset(cluster_counts, N_Infection=="Third")

granular_third_n_cells <- granular_third_cluster_counts %>%
  group_by(Sample_ID)%>%
  summarise("sum"=sum(Count))

granular_third_list_of_dfs_for_glm <- split(granular_third_cluster_counts, granular_third_cluster_counts$Cluster_ID)



granular_third_list_of_models <- lapply(granular_third_list_of_dfs_for_glm, function(x) glm(Percentage~Timepoint+Volunteer, family = "binomial", weights = granular_third_n_cells$sum, data=x))
#list_of_models <- lapply(list_of_dfs_for_glm, function(x) nlme::lme(concentration~timepoint, random=~1|Volunteer, data=x))
#list_of_models <- lapply(list_of_dfs_for_glm, function(x) lme4::lmer(concentration~timepoint+(1|Volunteer), data=x))

t6_contrast <- t(matrix(c(0,1,0,0)))

# c45_contrast <- t(diffcyt::createContrast(c(0,1,0,0,0,0,0,0)))
##
granular_third_list_of_tests <- lapply(granular_third_list_of_models, function(x) multcomp::glht(x, t6_contrast))
granular_third_list_of_pvalues <- lapply(granular_third_list_of_tests, function(x) data.frame("p_raw"=summary(x)$test$pvalues, "coef"=summary(x)$test$coefficients))

granular_third_df_first_ <- do.call(rbind, granular_third_list_of_pvalues)

granular_third_df_first_$Cluster_ID <- names(granular_third_list_of_pvalues)
#
granular_third_df_first_$p_adj <- p.adjust(granular_third_df_first_$p_raw, method = "fdr")


#log(2) = 0.6931472
subset(granular_third_df_first_, p_adj<=0.05 & abs(coef)>0.6931472)

sig_granular_third <- subset(granular_third_df_first_, p_adj<=0.05 & abs(coef)>0.6931472)
sig_granular_third <- data.frame("cluster"=sig_granular_third$Cluster_ID, "raw_FC"=exp(sig_granular_third$coef), "p_adj"=sig_granular_third$p_adj)

# p_raw       coef Cluster_ID        p_adj
# 4  0.000000e+00 -1.0135963          4 0.000000e+00
# 8  4.764177e-09 -0.7272765          8 4.128954e-08
# 10 1.068725e-07 -0.8859015         10 6.946711e-07
# 23 2.580773e-05  1.1334626         23 8.387511e-05
# 24 2.606005e-02  2.3119246         24 4.517075e-02

# normal first ####

cluster_counts <- read.csv("~/postdoc/scRNAseq/final_analysis_and_figures/final_activated_cluster_counts.csv")

first_cluster_counts <- subset(cluster_counts, N_Infection=="First")

first_n_cells <- first_cluster_counts %>%
  group_by(Sample_ID)%>%
  summarise("sum"=sum(Count))

first_list_of_dfs_for_glm <- split(first_cluster_counts, first_cluster_counts$Cluster_ID)



first_list_of_models <- lapply(first_list_of_dfs_for_glm, function(x) glm(Percentage~Timepoint+Volunteer, family = "binomial", weights = first_n_cells$sum, data=x))
#list_of_models <- lapply(list_of_dfs_for_glm, function(x) nlme::lme(concentration~timepoint, random=~1|Volunteer, data=x))
#list_of_models <- lapply(list_of_dfs_for_glm, function(x) lme4::lmer(concentration~timepoint+(1|Volunteer), data=x))

t6_contrast <- t(matrix(c(0,1,0,0)))

# c45_contrast <- t(diffcyt::createContrast(c(0,1,0,0,0,0,0,0)))
##
first_list_of_tests <- lapply(first_list_of_models, function(x) multcomp::glht(x, t6_contrast))
first_list_of_pvalues <- lapply(first_list_of_tests, function(x) data.frame("p_raw"=summary(x)$test$pvalues, "coef"=summary(x)$test$coefficients))

first_df_first_ <- do.call(rbind, first_list_of_pvalues)

first_df_first_$Cluster_ID <- names(first_list_of_pvalues)
#
first_df_first_$p_adj <- p.adjust(first_df_first_$p_raw, method = "fdr")


#log(2) = 0.6931472
subset(first_df_first_, p_adj<=0.05 & abs(coef)>0.6931472)

sig_first <- subset(first_df_first_, p_adj<=0.05 & abs(coef)>0.6931472)
sig_first <- data.frame("cluster"=sig_first$Cluster_ID, "raw_FC"=exp(sig_first$coef), "p_adj"=sig_first$p_adj)

# p_raw       coef Cluster_ID        p_adj
# 4 3.891074e-10 -0.8166622          4 1.426727e-09
# 5 2.745214e-07 -0.9382784          5 7.549339e-07
# 6 0.000000e+00 -2.2167921          6 0.000000e+00
# 9 0.000000e+00  2.9567040          9 0.000000e+00


# normal third ####

cluster_counts <- read.csv("~/postdoc/scRNAseq/final_analysis_and_figures/final_activated_cluster_counts.csv")

third_cluster_counts <- subset(cluster_counts, N_Infection=="Third")

third_n_cells <- third_cluster_counts %>%
  group_by(Sample_ID)%>%
  summarise("sum"=sum(Count))

third_list_of_dfs_for_glm <- split(third_cluster_counts, third_cluster_counts$Cluster_ID)



third_list_of_models <- lapply(third_list_of_dfs_for_glm, function(x) glm(Percentage~Timepoint+Volunteer, family = "binomial", weights = third_n_cells$sum, data=x))
#list_of_models <- lapply(list_of_dfs_for_glm, function(x) nlme::lme(concentration~timepoint, random=~1|Volunteer, data=x))
#list_of_models <- lapply(list_of_dfs_for_glm, function(x) lme4::lmer(concentration~timepoint+(1|Volunteer), data=x))

t6_contrast <- t(matrix(c(0,1,0,0)))

# c45_contrast <- t(diffcyt::createContrast(c(0,1,0,0,0,0,0,0)))
##
third_list_of_tests <- lapply(third_list_of_models, function(x) multcomp::glht(x, t6_contrast))
third_list_of_pvalues <- lapply(third_list_of_tests, function(x) data.frame("p_raw"=summary(x)$test$pvalues, "coef"=summary(x)$test$coefficients))

third_df_first_ <- do.call(rbind, third_list_of_pvalues)

third_df_first_$Cluster_ID <- names(third_list_of_pvalues)
#
third_df_first_$p_adj <- p.adjust(third_df_first_$p_raw, method = "fdr")


#log(2) = 0.6931472
subset(third_df_first_, p_adj<=0.05 & abs(coef)>0.6931472)

sig_third <- subset(third_df_first_, p_adj<=0.05 & abs(coef)>0.6931472)
sig_third <- data.frame("cluster"=sig_third$Cluster_ID, "raw_FC"=exp(sig_third$coef), "p_adj"=sig_third$p_adj)

# p_raw       coef Cluster_ID        p_adj
# 2 0.000000e+00 -0.9743647          2 0.000000e+00
# 7 6.575629e-12  0.9002093          7 3.616596e-11
# 9 1.056613e-03  0.7766124          9 2.324548e-03
