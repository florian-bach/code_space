## extra script, just for stats ###

library(CATALYST)
library(diffcyt)
library(vac69a.cytof)
library(SummarizedExperiment)
library(SingleCellExperiment)

#daf <- read_small("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/", proportional = T, event_number = 1000)
daf <- read_full("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")
#
merging_table1 <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/merging_table_april2020.csv", header=T, stringsAsFactors = F)

#get rid of spaces at beginning of string
merging_table1$new_cluster <- ifelse(substr(merging_table1$new_cluster, 1, 1)==" ", substr(merging_table1$new_cluster, 2, nchar(merging_table1$new_cluster)), merging_table1$new_cluster)
merging_table1$new_cluster <- ifelse(substr(merging_table1$new_cluster, 1, 1)==" ", substr(merging_table1$new_cluster, 2, nchar(merging_table1$new_cluster)), merging_table1$new_cluster)

merging_table1$new_cluster <- factor(merging_table1$new_cluster)

merged_daf<- mergeClusters(daf, k = "meta45", table = merging_table1, id = "flo_merge")


#merged_daf <- daf

ei <- metadata(merged_daf)$experiment_info
#ei$timepoint <- factor(ei$timepoint, levels=c("C10", "Baseline", "DoD", "T6"))

design <- createDesignMatrix(ei, c("timepoint", "volunteer"))
# 
# design <- model.matrix(~ei$time+ei$time:ei$volunteer)
# batch_design <- createDesignMatrix(ei, c("timepoint", "t"))

FDR_cutoff <- 0.05

pairwise_contrast_t6 <- createContrast(c(c(0, 0, 0, 1), rep(0,5)))
pairwise_contrast_dod <- createContrast(c(c(0, 0, 1, 0), rep(0,5)))
pairwise_contrast_c10 <- createContrast(c(c(0, 1, 0, 0), rep(0,5)))

da_c10 <- diffcyt(merged_daf,
                  design = design,
                  #contrast = contrast_c10,
                  contrast = pairwise_contrast_c10,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "flo_merge",
                  verbose = T)


da_dod <- diffcyt(merged_daf,
                  design = design,
                  #contrast = contrast_dod,
                  contrast = pairwise_contrast_dod,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "flo_merge",
                  verbose = T)


da_t6 <- diffcyt(merged_daf,
                 design = design,
                 #contrast = contrast_t6,
                 contrast = pairwise_contrast_t6,
                 analysis_type = "DA",
                 method_DA = "diffcyt-DA-edgeR",
                 clustering_to_use = "flo_merge",
                 verbose = T)

table(rowData(da_c10$res)$p_adj < FDR_cutoff)
# FALSE
#    42
table(rowData(da_dod$res)$p_adj < FDR_cutoff)
# # FALSEFALSE  TRUE
# 34     8
table(rowData(da_t6$res)$p_adj < FDR_cutoff)
# # FALSEFALSE  TRUE
# 28     14


# with C+10 being the intercept term, c(-1, 1) contrasts can be made between baseline and dod/t6
# the significant clusters at dod are 10, 12, 14, 37, 3, 6, 43, 45, (n=8) p_adjust range from 2.5e-6 to 2.3e-2
# the significant clusters at t6 are 37, 15, 43, 10, 36, 12, 42, 14, 6, 8, 3, 29, 18, 19 (n=14), p+adjust range from 5.2e-52 to 4.4e-2

# these are exactly the same clusters, when baseline is the intercept term and we set it to 0 for the contrasts
# this is good news!!!

# using a design matrix with timepoint and batch reduces our dod clusters to 0 and the t6 clusters to 6 (37, 43, 36, 15, 42, 12)
# removing the cluster term returns 0 significant clusters at dod and 5 at t6 (37, 32, 36, 15, 42) so there's a small effect;
# all the mismatched clusters between those models are detected in the full design matrix
plotDiffHeatmap(merged_daf, da_c10, th = FDR_cutoff, normalize = TRUE, hm1 = F, top_n = 25)
plotDiffHeatmap(merged_daf, da_dod, th = FDR_cutoff, normalize = TRUE, hm1 = F, top_n = 25)
plotDiffHeatmap(merged_daf, da_t6, th = FDR_cutoff, normalize = TRUE, hm1 = F, top_n = 25)




# glmms #

#this one works, don't change
da_formula1 <- createFormula(ei, cols_fixed = "timepoint",
                             cols_random = "volunteer")

# this one you're allowe to play with
da_formula2 <- createFormula(ei, cols_fixed = c("timepoint", "volunteer"),
                             cols_random = c("sample_id"))



glm_contrast_c10 <- createContrast(c(0, 1, 0, 0))
glm_contrast_dod <- createContrast(c(0, 0, 1, 0))
glm_contrast_t6 <- createContrast(c(0, 0, 0, 1))

glm_c10 <- diffcyt(merged_daf,
                      formula = da_formula2,
                      contrast = glm_contrast_c10,
                      analysis_type = "DA",
                      method_DA = "diffcyt-DA-GLMM",
                      clustering_to_use = "meta45",
                      verbose = T)

glm_dod <- diffcyt(merged_daf,
                      formula = da_formula2,
                      contrast = glm_contrast_dod,
                      analysis_type = "DA",
                      method_DA = "diffcyt-DA-GLMM",
                      clustering_to_use = "meta45",
                      verbose = T)

glm_t6 <- diffcyt(merged_daf,
                     contrast = glm_contrast_t6,
                     formula = da_formula2,
                     analysis_type = "DA",
                     method_DA = "diffcyt-DA-GLMM",
                     clustering_to_use = "meta45",
                     verbose = T)



table(rowData(glm_c10$res)$p_adj < FDR_cutoff)
# FALSE  TRUE 
# 18    27
table(rowData(glm_dod$res)$p_adj < FDR_cutoff)
# FALSE  TRUE 
# 9    36  
table(rowData(glm_t6$res)$p_adj < FDR_cutoff)
# FALSE  TRUE 
# 3    42


plotDiffHeatmap(merged_daf, glm_c10, th = FDR_cutoff, normalize = TRUE, hm1 = F, top_n = 25)
plotDiffHeatmap(merged_daf, glm_dod, th = FDR_cutoff, normalize = TRUE, hm1 = F, top_n = 25)
plotDiffHeatmap(merged_daf, glm_t6, th = FDR_cutoff, normalize = TRUE, hm1 = F, top_n = 25)



glm_c10_sig <- topTable(glm_c10, top_n = 5)$cluster_id
glm_dod_sig <- topTable(glm_dod, top_n = 18)$cluster_id
glm_t6_sig <- topTable(glm_t6, top_n = 24)$cluster_id


edger_dod_sig <- topTable(da_dod, top_n = 8)$cluster_id
edger_t6_sig <- topTable(da_t6, top_n = 14)$cluster_id

all(edger_dod_sig %in% glm_dod_sig) # TRUE
all(edger_t6_sig %in% glm_t6_sig) # TRUE

glm_c10_sig %in% glm_dod_sig





ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/all_clusters_counts.png", diffcyt_boxplot(glm_t6, merged_daf, counts=T, FDR=0.5), height = 12, width=12)# works
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/all_clusters_freqs.png", diffcyt_boxplot(glm_t6, merged_daf, counts=F, FDR=0.5), height = 12, width=12)# works

diffcyt_boxplot(glm_c10, merged_daf, FDR=0.01)# works
diffcyt_boxplot(glm_c10, merged_daf, counts=T, FDR=0.01)# works



