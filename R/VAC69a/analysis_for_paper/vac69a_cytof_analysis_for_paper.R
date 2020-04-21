library(diffcyt)
library(tidyr)
library(dplyr)
library(cowplot)
library(SummarizedExperiment)
library(SingleCellExperiment)

library(vac69a.cytof)
library(CATALYST)

#a couple of extra functions defined here ####
`%!in%` = Negate(`%in%`)

inferno <- colorspace::sequential_hcl("inferno", n=8)
inferno_white <- c("#FFFFFF", colorspace::sequential_hcl("inferno", n=8))
inferno_mega_lite <- c("#000004", "#8A2267", "#EF802B", "#FFEC89", "#FCFFA4")

daf <- read_full("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")

plotClusterHeatmap(daf, hm2=NULL,
                   k = "meta45",
                   #m = "flo_merge",
                   cluster_anno = TRUE,
                   draw_freqs = TRUE,
                   scale = TRUE, 
                   palette=inferno
)


# start <- Sys.time(); daf <- runUMAP(daf, exprs_values = "exprs", feature_set=refined_markers); end <- Sys.time()
# (duration <- round(end - start)) # ~20-25 min


#print(c("UMAP took ", duration, " minutes to run on 861161 cells"), sep='')

# start <- Sys.time(); daf <- runTSNE(daf, exprs_values = "exprs", feature_set=refined_markers); end <- Sys.time()
# duration <- round(end - start)
# print(c("tSNE took ", duration, " minutes to run on 861161 cells"), sep='')
# 
# 

# print(c("together, UMAP and tSNE took ", duration, " minutes to run on 10,000 cells"), sep='')




# daf100_delta <- metadata(daf100)$delta_area

#### get rid of weird super positive events ####


# removing a cluster messes with the differential expression calculation for some reason- maybe because
# it's a factor and removing the cluster doens't change the entry in the dictionary?



# cluster_ids <- seq(1,45)
# bad_clusters <- 24
# cluster_ids <- as.character(cluster_ids[-bad_clusters])
# 
# daffy <- filterSCE(merged_daf, k = "meta45", cluster_id %!in% list("trash"))





# th35 <- filterSCE(daf, k = "meta35", cluster_id %in% paste(c(16, 20, 21)))
# th40 <- filterSCE(daf, k = "meta40", cluster_id %in% paste(c(17, 23, 22, 18)))
# 
# th35 <- runUMAP(th35, exprs_values = "exprs", feature_set=refined_markers)
# th40 <- runUMAP(th40, exprs_values = "exprs", feature_set=refined_markers)
# 
# th35_plot <- plotDR(th35, color_by = "meta35")
# th40_plot <-plotDR(th40, color_by = "meta40")
# 
# plotClusterExprs(th35, k = "meta35", features = "type")
# plotClusterExprs(th40, k = "meta40", features = "type")

####  the mergening 

#32 clusters, mostly delineated along CD28, CD27, CD57, after lineage, memory and activation markers
merging_table1 <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/merging_table_april2020.csv", header=T, stringsAsFactors = F)

#get rid of spaces at beginning of string
merging_table1$new_cluster <- ifelse(substr(merging_table1$new_cluster, 1, 1)==" ", substr(merging_table1$new_cluster, 2, nchar(merging_table1$new_cluster)), merging_table1$new_cluster)
merging_table1$new_cluster <- ifelse(substr(merging_table1$new_cluster, 1, 1)==" ", substr(merging_table1$new_cluster, 2, nchar(merging_table1$new_cluster)), merging_table1$new_cluster)

merging_table1$new_cluster <- factor(merging_table1$new_cluster)

merged_daf<- mergeClusters(daf, k = "meta45", table = merging_table1, id = "flo_merge")

flo_32_cluster_heatmap <- plotClusterHeatmap(merged_daf, hm2=NULL,
                   k = "meta45",
                   #m = "flo_merge",
                   cluster_anno = FALSE,
                   draw_freqs = TRUE,
                   scale = TRUE, 
                   palette=inferno
)

ggsave("flo_32_cluster_heatmap.png", flo_32_cluster_heatmap)

### diffcyt ####z
ei <- metadata(merged_daf)$experiment_info

design <- createDesignMatrix(ei, c("timepoint", "volunteer"))

levels(ei$timepoint)

###If a design matrix has been used, the entries of contrast correspond to the columns of the design
#matrix and the length of contrast equals the number of columns in the design matrix. If a model formula
#has been used, the entries correspond to the levels of the fixed effect terms;
#and the length equals the number of levels of the fixed effect terms.



FDR_cutoff <- 0.05

# edgeR models with all timepoints####

contrast_baseline <- createContrast(c(-1, rep(0, 8)))
contrast_c10 <- createContrast(c(c(0, 1), rep(0,7)))
contrast_dod<- createContrast(c(c(0, 0, 1), rep(0,6)))
contrast_t6 <- createContrast(c(c(0, 0, 0, 1), rep(0,5)))


da_baseline <- diffcyt(merged_daf,
                       design = design,
                       contrast = contrast_baseline,
                       analysis_type = "DA",
                       method_DA = "diffcyt-DA-edgeR",
                       clustering_to_use = "flo_merge",
                       verbose = T)

da_c10 <- diffcyt(merged_daf,
                  design = design,
                  contrast = contrast_c10,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "flo_merge",
                  verbose = T)

da_dod <- diffcyt(merged_daf,
                  design = design,
                  contrast = contrast_dod,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "flo_merge",
                  verbose = T)

da_t6 <- diffcyt(merged_daf,
                 design = design,
                 contrast = contrast_t6,
                 analysis_type = "DA",
                 method_DA = "diffcyt-DA-edgeR",
                 clustering_to_use = "flo_merge",
                 verbose = T)

# results  ###
table(rowData(da_baseline$res)$p_adj < FDR_cutoff)


table(rowData(da_c10$res)$p_adj < FDR_cutoff)
# FALSE 
# 38 
table(rowData(da_dod$res)$p_adj < FDR_cutoff)
# FALSE  TRUE 
# 24    8
table(rowData(da_t6$res)$p_adj < FDR_cutoff)
# FALSE  TRUE 
# 23    15 

plotDiffHeatmap(merged_daf, da_baseline, th = FDR_cutoff, normalize = TRUE, hm1 = F)
plotDiffHeatmap(merged_daf, da_c10, th = FDR_cutoff, normalize = TRUE, hm1 = F)
plotDiffHeatmap(merged_daf, da_dod, th = FDR_cutoff, normalize = TRUE, hm1 = F)
plotDiffHeatmap(merged_daf, da_t6, th = FDR_cutoff, normalize = TRUE, hm1 = F)



#### edger with two timepoints for pairwise comparisons #### 

# first we filter the original sce to exclude unnecessary timepoints #
base_c10 <- filterSCE(merged_daf, timepoint %in% c("Baseline", "C10"))
base_dod <- filterSCE(merged_daf, timepoint %in% c("Baseline", "DoD"))
base_t6 <- filterSCE(merged_daf, timepoint %in% c("Baseline", "T6"))


# then, for each dataset, create design matrix & change metadata to change the factor levels which until now
# still include the timepoints and sample_ids that aren't present anymore; if you don't do the latter step the
# diffcyt function will fail because of an internal sanity check that attempts to match sample_ids between the
# object that contains cluster medians and the sce being tested
ei_c10 <- metadata(base_c10)$experiment_info
ei_c10$timepoint <- factor(ei_c10$timepoint, levels=c("Baseline", "C10"))
ei_c10$sample_id <- factor(as.character(ei_c10$sample_id))
c10_design <- createDesignMatrix(ei_c10, c("timepoint", "volunteer"))
metadata(base_c10)$experiment_info <- ei_c10

ei_dod <- metadata(base_dod)$experiment_info
ei_dod$timepoint <- factor(ei_dod$timepoint, levels=c("Baseline", "DoD"))
ei_dod$sample_id <- factor(as.character(ei_dod$sample_id))
dod_design <- createDesignMatrix(ei_dod, c("timepoint", "volunteer"))
metadata(base_dod)$experiment_info <- ei_dod

ei_t6 <- metadata(base_t6)$experiment_info
ei_t6$timepoint <- factor(ei_t6$timepoint, levels=c("Baseline", "T6"))
ei_t6$sample_id <- factor(as.character(ei_t6$sample_id))
t6_design <- createDesignMatrix(ei_t6, c("timepoint", "volunteer"))
metadata(base_t6)$experiment_info <- ei_t6

levels(ei_c10$timepoint)

# 
# contrast_dod<- createContrast(c(c(0, 0, 1), rep(0,6)))
# contrast_t6 <- createContrast(c(c(1, 0, 0, -1), rep(0,5)))

pairwise_contrast <- createContrast(c(c(0, 1), rep(0,5)))

pair_base_c10 <- diffcyt(base_c10,
                  design = c10_design,
                  contrast = pairwise_contrast,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "flo_merge",
                  verbose = T)

pair_base_dod <- diffcyt(base_dod,
                  design = dod_design,
                  contrast = pairwise_contrast,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "flo_merge",
                  verbose = T)

pair_base_t6 <- diffcyt(base_t6,
                 design = t6_design,
                 contrast = pairwise_contrast,
                 analysis_type = "DA",
                 method_DA = "diffcyt-DA-edgeR",
                 clustering_to_use = "flo_merge",
                 verbose = T)


table(rowData(pair_base_c10$res)$p_adj < FDR_cutoff)
# FALSE  TRUE 
# 36     2
table(rowData(pair_base_dod$res)$p_adj < FDR_cutoff)
# FALSE  TRUE 
# 34     4
table(rowData(pair_base_t6$res)$p_adj < FDR_cutoff)
# FALSE  TRUE 
# 20    18


plotDiffHeatmap(base_c10, pair_base_c10, th = FDR_cutoff, normalize = TRUE, hm1 = T)
plotDiffHeatmap(base_dod, pair_base_dod, th = FDR_cutoff, normalize = TRUE, hm1 = F)
plotDiffHeatmap(base_t6, pair_base_t6, th = FDR_cutoff, normalize = TRUE, hm1 = F)




# glmm models ####

contrast_baseline <- createContrast(c(1,0,0,0))
contrast_c10 <- createContrast(c(0,1,0,0))
contrast_dod <- createContrast(c(0,0,1,0))
contrast_t6 <- createContrast(c(0,0,0,1))

da_formula2 <- createFormula(ei,
                             cols_fixed = c("timepoint"),
                             cols_random ="sample_id")

da_baseline_vol <- diffcyt(merged_daf,
                           #design = design,
                           formula = da_formula2,
                           contrast = contrast_baseline,
                           analysis_type = "DA",
                           method_DA = "diffcyt-DA-GLMM",
                           clustering_to_use = "flo_merge",
                           verbose = T)

da_c10_vol <- diffcyt(merged_daf,
                      #design = design,
                      formula = da_formula2,
                      contrast = contrast_c10,
                      analysis_type = "DA",
                      method_DA = "diffcyt-DA-GLMM",
                      clustering_to_use = "flo_merge",
                      verbose = T)

da_dod_vol <- diffcyt(merged_daf,
                      #design = design,
                      formula = da_formula2,
                      contrast = contrast_dod,
                      analysis_type = "DA",
                      method_DA = "diffcyt-DA-GLMM",
                      clustering_to_use = "flo_merge",
                      verbose = T)

da_t6_vol <- diffcyt(merged_daf,
                     #design = design,
                     formula = c,
                     contrast = contrast_t6,
                     analysis_type = "DA",
                     method_DA = "diffcyt-DA-GLMM",
                     clustering_to_use = "flo_merge",
                     verbose = T)

### results using glmm, <y ~ timepoint + (1 | sample_id)>

table(rowData(da_dod_vol$res)$p_adj < FDR_cutoff)
# # FALSE
# # 35
table(rowData(da_t6_vol$res)$p_adj < FDR_cutoff)
# # FALSE  TRUE
# # 23
9



# it only picks up on conserved stuff, but that makes for an easier story, 8/9 clusters are up..

countz <- calcCounts(daf)

testDA_GLMM(
  d_counts,
  formula,
  contrast,
  min_cells = 3,
  min_samples = NULL,
  normalize = FALSE,
  norm_factors = "TMM"
)



### results using glmm, <y ~ timepoint + (1| sample_ID) +(1 | volunteer)>
table(rowData(da_baseline_vol$res)$p_adj < FDR_cutoff)
# FALSE 
# 35 

table(rowData(da_c10_vol$res)$p_adj < FDR_cutoff)
# FALSE   
# 14

table(rowData(da_dod_vol$res)$p_adj < FDR_cutoff)
# FALSE  TRUE 
# 16    19

table(rowData(da_t6_vol$res)$p_adj < FDR_cutoff)
# FALSE  TRUE 
# 23     9 

plotDiffHeatmap(daf, da_baseline_vol, th = FDR_cutoff, normalize = TRUE, hm1 = T)
plotDiffHeatmap(daf, da_c10_vol, th = FDR_cutoff, normalize = TRUE, hm1 = T)
plotDiffHeatmap(merged_daf, da_dod_vol, th = FDR_cutoff, normalize = TRUE, hm1 = F)
plotDiffHeatmap(merged_daf, da_t6_vol, th = FDR_cutoff, normalize = TRUE, hm1 = F)

# adjusting for individual differences with a random effect includes some spurious looking results, but 
# also includes some more that i feel are missing if there's no attempt to control for individual identity
# adding a fixed effect somehow breaks it?? might be worth running lme4 oldschool to check what's
# going on, maybe speak to some IEB people about this again...


###  make boxplots of cluster counts/frequencies ####
# topTable(da_t6_vol, show_counts = T)
dod_vol <- data.frame(topTable(pair_base_dod, _freq=T, show_counts = T))
up_dod <-  dplyr::filter(dod_vol, dod_vol$p_adj < FDR_cutoff)

long_up_dod <- gather(up_dod, sample_id, count, colnames(up_dod)[4:ncol(up_dod)])
long_up_dod$volunteer <- stringr::str_match(long_up_dod$sample_id, "V[0-9]*")[, 1]
long_up_dod$timepoint <- substr(long_up_dod$sample_id, 12,nchar(long_up_dod$sample_id))

long_up_dod$sample_id <- gsub("counts_", "", long_up_dod$sample_id)
counts <- n_cells(daf)
long_up_dod$frequency <- long_up_dod$count / counts[long_up_dod$sample_id] *100

ggplot(long_up_dod, aes(x=factor(long_up_dod$timepoint), y=long_up_dod$count))+
  geom_boxplot(aes(fill=long_up_dod$timepoint))+
  geom_point(aes(shape=long_up_dod$volunteer))+
  facet_wrap(~long_up_dod$cluster_id, scales = "free", ncol=5)+
  theme_minimal()+
  theme(axis.title = element_blank(),
        legend.title = element_blank())


t6_vol <- data.frame(topTable(pair_base_t6, all=T, show_counts = T))
up_t6 <-  dplyr::filter(t6_vol, t6_vol$p_adj < FDR_cutoff)

long_up_t6 <- gather(up_t6, sample_id, count, colnames(up_t6)[4:ncol(up_t6)])
long_up_t6$volunteer <- stringr::str_match(long_up_t6$sample_id, "V[0-9]*")[, 1]
long_up_t6$timepoint <- substr(long_up_t6$sample_id, 12,nchar(long_up_t6$sample_id))

long_up_t6$sample_id <- gsub("counts_", "", long_up_t6$sample_id)
counts <- n_cells(daf)
long_up_t6$frequency <- long_up_t6$count / counts[long_up_t6$sample_id] *100

ggplot(long_up_t6, aes(x=factor(long_up_t6$timepoint), y=long_up_t6$count))+
  geom_boxplot(aes(fill=long_up_t6$timepoint))+
  geom_point(aes(shape=long_up_t6$volunteer))+
  facet_wrap(~long_up_t6$cluster_id, scales = "free", ncol=5)+
  theme_minimal()+
  theme(axis.title = element_blank(),
        legend.title = element_blank())


# cluster x cluster correlation matrices
all_frequencies <- data.frame(topTable(da_baseline, all=T, show_counts = T))

all_frequencies <- gather(all_frequencies, sample_id, count, colnames(all_frequencies)[4:ncol(all_frequencies)])
all_frequencies$volunteer <- stringr::str_match(all_frequencies$sample_id, "V[0-9]*")[, 1]
all_frequencies$timepoint <- substr(all_frequencies$sample_id, 12,nchar(all_frequencies$sample_id))

all_frequencies$sample_id <- gsub("counts_", "", all_frequencies$sample_id)
counts <- n_cells(daf)
all_frequencies$frequency <- all_frequencies$count / counts[all_frequencies$sample_id] *100

baseline_freq_matrix <- all_frequencies %>%
  dplyr::filter(timepoint=="Baseline") %>%
  select(cluster_id, volunteer, timepoint, frequency)


# making correlation heatmaps

baseline_freq_matrix <- spread(baseline_freq_matrix, cluster_id, frequency)
baseline_spearman <- cor(baseline_freq_matrix[,3:ncol(baseline_freq_matrix)], method = "p")

unit <- dist(baseline_spearman, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
units <- hclust(unit)

specific_order <- colnames(baseline_spearman)[units$order]

#check.names=FALSE here makes sure that the +/- symbols parse and spaces aren't dots
baseline_spearman_df  <- data.frame(baseline_spearman, check.names = FALSE)
baseline_spearman_df$cluster_id_x <- rownames(baseline_spearman_df)

long_baseline_spearman <- gather(baseline_spearman_df, cluster_id_y, ro, colnames(baseline_spearman_df)[1:ncol(baseline_spearman_df)-1])

# hclust_levels <- colnames(baseline_spearman)[baselin_hclust$order]

corr_matrix_theme <-
  theme(axis.title = element_blank(),
        # axis.text.x = element_text(angle = 60, hjust = 1),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.position = "none")



ggplot(long_baseline_spearman, aes_(x=factor(long_baseline_spearman$cluster_id_x, levels = colnames(baseline_spearman)[units$order]), y=factor(long_baseline_spearman$cluster_id_y, levels=colnames(baseline_spearman)[units$order])))+
    geom_tile(aes(fill=long_baseline_spearman$ro))+
    viridis::scale_fill_viridis(option="A")+
    ggtitle("Baseline")+
    labs(fill = expression(paste("Pearson ", rho)))+
    corr_matrix_theme


# [1] "counts_V02_Baseline" "counts_V02_Baseline" "counts_V02_Baseline" "counts_V02_Baseline"
# [5] "counts_V02_Baseline" "counts_V02_Baseline"

## lm models for differential marker expression####

ei <- metadata(merged_daf)$experiment_info

ds_formula1 <- createFormula(ei, cols_fixed = "timepoint",
                             cols_random = "sample_id")

ds_formula2 <- createFormula(ei, cols_fixed = "timepoint",
                             cols_random = c("sample_id", "volunteer"))

design <- createDesignMatrix(ei, c("timepoint", "volunteer"))


lm_contrast_dod <- createContrast(c(0,0,1,0))
lm_contrast_t6 <- createContrast(c(0,0,0,1))

ds_dod_1 <- diffcyt(merged_daf,                                            
                   formula = ds_formula1, contrast = contrast_dod, design= design,                    
                   analysis_type = "DS",            
                   clustering_to_use = "flo_merge", verbose = T)               

table(rowData(ds_dod_1$res)$p_adj < FDR_cutoff)

ds_t6_1 <- diffcyt(merged_daf,                                            
                    formula = ds_formula1, contrast = lm_contrast_t6,                    
                    analysis_type = "DS", method_DS = "diffcyt-DS-LMM",            
                    clustering_to_use = "flo_merge", verbose = FALSE)               
table(rowData(ds_t6_1$res)$p_adj < FDR_cutoff)         



ds_dod_2 <- diffcyt(merged_daf,                                            
                    formula = ds_formula2, contrast = lm_contrast_dod,                    
                    analysis_type = "DS", method_DS = "diffcyt-DS-LMM",            
                    clustering_to_use = "flo_merge", verbose = FALSE)               
table(rowData(ds_dod_2$res)$p_adj < FDR_cutoff)

ds_t6_2 <- diffcyt(merged_daf,                                            
                   formula = ds_formula2, contrast = lm_contrast_t6,                    
                   analysis_type = "DS", method_DS = "diffcyt-DS-LMM",            
                   clustering_to_use = "flo_merge", verbose = FALSE)               
table(rowData(ds_t6_2$res)$p_adj < FDR_cutoff)         

        

