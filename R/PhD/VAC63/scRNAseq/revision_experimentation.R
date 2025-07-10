# granular activation stuff ####
activated_subset <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/final_analysis_and_figures/final_activated_subset.h5Seurat")

library(Seurat)
library(dplyr)

granular_activation <- subset(activated_subset, seurat_clusters %in% c(2,4:6,9))

#granular_activation <- RunUMAP(granular_activation, dims = 1:12)


granular_activation <- FindNeighbors(granular_activation, dims = 1:12)
#0.8
granular_activation <- FindClusters(granular_activation,
                                    resolution = 0.8,
                                    algorithm = 1,
                                    );DimPlot(granular_activation, label=TRUE)

FeaturePlot(granular_activation, features = c("CD38", "HLA-DRB1", "KLRB1", "TCF7"), order=TRUE)&viridis::scale_color_viridis(option="A")

combo.markers <- FindAllMarkers(granular_activation, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top20_per_cluster <- data.frame(combo.markers %>%
                                  filter(cluster %in% c(5,6)) %>%
                                  group_by(cluster) %>%
                                  slice_min(n = 30, order_by = p_val_adj))

write.csv(top20_per_cluster, "~/postdoc/scRNAseq/activated_subgroup_top_genes.csv", row.names = FALSE, quote = FALSE)

cluster_counts <- data.frame(table(granular_activation@active.ident, granular_activation@meta.data$sample_ID))
colnames(cluster_counts) <- c("Cluster_ID", "Sample_ID", "Count")
cluster_counts$Timepoint <- substr(cluster_counts$Sample_ID, 6, 14)
cluster_counts$Volunteer <- substr(cluster_counts$Sample_ID, 1, 4)
cluster_counts$N_Infection <- ifelse(cluster_counts$Volunteer %in% c("v313", "v315", "v320"), "First", "Third")


cluster_counts <- cluster_counts %>%
  group_by(Sample_ID) %>%
  mutate("Percentage"=Count/sum(Count))

cd38hi <- filter(cluster_counts, Cluster_ID %in% c(5,6))


ggplot(cd38hi, aes(x=factor(Cluster_ID), y=Percentage, fill=Timepoint))+
  geom_boxplot()+
  geom_point(position = position_dodge(width = 0.75))+#+
  scale_y_continuous(labels = scales::label_percent())+
  scale_fill_manual(values=c("#AA3377", "#4477AA"))+
  facet_wrap(~N_Infection)+
  theme_minimal()

ggplot(cd38hi, aes(x=factor(Cluster_ID), y=Percentage, fill=Timepoint, label=Count))+
  geom_boxplot()+
  geom_point(position = position_dodge(width = 0.75))+
  #geom_text(position = position_dodge(width = 0.9))+
  ggrepel::geom_label_repel()+
  scale_y_continuous(labels = scales::label_percent())+
  scale_fill_manual(values=c("#AA3377", "#4477AA"))+
  facet_wrap(~N_Infection)+
  theme_minimal()




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



# milo stuff ####


activated_subset <- SeuratDisk::LoadH5Seurat("~/postdoc/edinburgh/scRNAseq/seurat_objecitons/final_activated_subset.h5Seurat")

library(miloR)
library(patchwork)
library(ggplot2)
library(SingleCellExperiment)

combo.sce <- Seurat::as.SingleCellExperiment(activated_subset)

combo_milo <- Milo(combo.sce)
combo_milo <- buildGraph(combo_milo, k = 30, d = 30, reduced.dim = "PCA")
combo_milo <- makeNhoods(combo_milo, prop = 0.1, k = 30, d=30, refined = TRUE, reduced_dims = "PCA")

# need this for colData & owData

#metadata <- read.csv("~/Documents/vac63c_scrna_seq/sce_metadata.csv", row.names = 1)
metadata <- read.csv("~/postdoc/edinburgh/scRNAseq/sce_metadata.csv", row.names = 1)

combo_milo <- countCells(combo_milo, meta.data = as.data.frame(colData(combo_milo)), sample="sample_ID")

combo_milo <- calcNhoodDistance(combo_milo, d=30, reduced.dim = "PCA")





da_results <- testNhoods(combo_milo, design = ~ n_infection, design.df = metadata)
combo_milo <- buildNhoodGraph(combo_milo)
## Plot single-cell UMAP
umap_pl <- scater::plotReducedDim(combo_milo, dimred = "UMAP", colour_by="n_infection",
                                  text_size = 3, point_size=0.5) +
  guides(fill="none")
## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(combo_milo, da_results, layout="UMAP",alpha=0.1)

(milo_umap <- umap_pl + nh_graph_pl +
  plot_layout(guides="collect"))

ggsave("~/postdoc/scRNAseq/figures/n_infection_only_milo_umap.png", milo_umap, height=8, width=12, bg="white", dpi=444)




da_results <- testNhoods(combo_milo, design = ~ n_infection + timepoint, design.df = metadata)
combo_milo <- buildNhoodGraph(combo_milo)
## Plot single-cell UMAP
umap_pl <- scater::plotReducedDim(combo_milo, dimred = "UMAP", colour_by="n_infection",
                                  text_size = 3, point_size=0.5) +
  guides(fill="none")
## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(combo_milo, da_results, layout="UMAP",alpha=0.1)

milo_umap <- umap_pl + nh_graph_pl +
  plot_layout(guides="collect")

ggsave("~/postdoc/scRNAseq/figures/n_infection_plus_timpoeint_milo_umap.png", milo_umap, height=8, width=12, bg="white", dpi=444)





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

ggsave("~/postdoc/scRNAseq/figures/timepoint_plus_volunteer_milo_umap.png", milo_umap, height=8, width=12, bg="white", dpi=444)






da_results <- testNhoods(combo_milo, design = ~ timepoint + n_infection, design.df = metadata)
combo_milo <- buildNhoodGraph(combo_milo)
## Plot single-cell UMAP
umap_pl <- scater::plotReducedDim(combo_milo, dimred = "UMAP", colour_by="timepoint",
                                  text_size = 3, point_size=0.5) +
  guides(fill="none")
## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(combo_milo, da_results, layout="UMAP",alpha=0.1)

milo_umap <- umap_pl + nh_graph_pl +
  plot_layout(guides="collect")

ggsave("~/postdoc/scRNAseq/figures/timepoint_plus_n_infection_milo_umap.png", milo_umap, height=8, width=12, bg="white", dpi=444)




# deg stuff ####

activated_subset <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/final_analysis_and_figures/final_activated_subset.h5Seurat")

library(Libra)
library(Seurat)

# first <- subset(activated_subset, subset = n_infection=="First")
# third <- subset(activated_subset, subset = n_infection=="Third")
# 



n_infection_de <- run_de(activated_subset, cell_type_col = "seurat_clusters", label_col = c("n_infection"), replicate_col="timepoint", de_method = "edgeR", de_type = "LRT", min_cells=10)
sig_n_infection_de <- subset(n_infection_de, p_val_adj<0.1& abs(avg_logFC) > log2(1.5))
dim(sig_n_infection_de)
write.csv(sig_n_infection_de, "~/postdoc/scRNAseq/figures/sig_n_infection_DE.csv", row.names = FALSE)
# lrt 169 gene-clsuter combinations; 155 genes


# misc stuff ####

umap_pl <- scater::plotReducedDim(combo_milo, dimred = "UMAP", colour_by="volunteer",
                                  text_size = 3, point_size=1) +
  guides(fill="none")+
  scale_color_manual(values = c("brown1", "brown3", "brown4", "cornflowerblue", "blue", "darkblue"))+
  facet_wrap(~colour_by, nrow = 2)


ggsave("~/postdoc/scRNAseq/figures/vol_umap.png", umap_pl, height=8, width=12, bg="white", dpi=444)

