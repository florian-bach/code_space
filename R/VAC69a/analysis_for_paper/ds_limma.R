library(diffcyt)
library(tidyr)
library(dplyr)
library(cowplot)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(ggplot2)
library(vac69a.cytof)
library(CATALYST)

### read in metadata etc. ####
#path_to_directory <- "~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/"
#
# character.check <- is.character(path_to_direcory)
# space.check <- length(path_to_direcory) == 1
#
# if(character.check==FALSE) stop("path_to_direcory isn't a string")
# if(space.check==FALSE) stop("path_to_direcory contains spaces")



### CATALYST ####
#construct daFrame #
## md has to have particular properties: file_name=NAMED LIST (chr), ID and everything else in factors

#daf <- read_small("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/", proportional = T, event_number = 3000)
daf <- read_full("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")

coarse_table <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/merging_tables/most_coarse_merge.csv", header=T, stringsAsFactors = F)



#get rid of spaces at beginning of string

coarse_table$new_cluster <- factor(coarse_table$new_cluster)

merged_daf<- mergeClusters(daf, k = "mod_meta45", table = coarse_table, id = "coarse_merge")


plotClusterHeatmap(merged_daf, hm2=NULL,
                   k = "mod_meta45",
                   m = "coarse_merge",
                   cluster_anno = FALSE,
                   draw_freqs = TRUE,
                   scale = TRUE
)



ei <- metadata(merged_daf)$experiment_info

design <- createDesignMatrix(ei, c("timepoint", "volunteer"))

levels(ei$timepoint)

###If a design matrix has been used, the entries of contrast correspond to the columns of the design
#matrix and the length of contrast equals the number of columns in the design matrix. If a model formula
#has been used, the entries correspond to the levels of the fixed effect terms;
#and the length equals the number of levels of the fixed effect terms.

FDR_cutoff <- 0.05

# edgeR models with all timepoints####

pairwise_contrast_t6 <- createContrast(c(c(0, 0, 0, 1), rep(0,5)))
pairwise_contrast_dod <- createContrast(c(c(0, 0, 1, 0), rep(0,5)))
pairwise_contrast_c10 <- createContrast(c(c(0, 1, 0, 0), rep(0,5)))

# copntrasts for models excluding volunteer
# pairwise_contrast_t6 <- createContrast(c(c(0, 0, 0, 1), rep(0,1)))
# pairwise_contrast_dod <- createContrast(c(c(0, 0, 1, 0), rep(0,1)))
# pairwise_contrast_c10 <- createContrast(c(c(0, 1, 0, 0), rep(0,1)))


states <- names(marker_classes(merged_daf))
states <- states[-c(1:21,
                    match(c("CCR7", "CD127", "CD20", "CXCR5", "CD103", "TCRgd", "TIM3", "IntegrinB7", "CD56", "CD3", "CD49d", "CD45RA", "CD4", "Vd2", "Va72", "CD161", "FoxP3", "CD45RO"), states),
                    53, 54, 56:length(states)
                    )
                 ]

# logic <- names(marker_classes(merged_daf)) %in% states

metadata(merged_daf)$id_state_markers <- states

ds_c10 <- diffcyt(merged_daf,                                            
                  design = design, contrast = pairwise_contrast_c10,                    
                  analysis_type = "DS", method_DS = "diffcyt-DS-limma",       
                  clustering_to_use = "coarse_merge", verbose = TRUE)               


ds_dod <- diffcyt(merged_daf,                                            
                   design = design, contrast = pairwise_contrast_dod,                    
                   analysis_type = "DS", method_DS = "diffcyt-DS-limma",       
                   clustering_to_use = "coarse_merge", verbose = TRUE)               


ds_t6 <- diffcyt(merged_daf,                                            
                   design = design, contrast = pairwise_contrast_t6,                    
                   analysis_type = "DS",  method_DS =  "diffcyt-DS-limma",         
                   clustering_to_use = "coarse_merge", verbose = TRUE)



topTable(ds_dod, top_n = 5, order_by = "cluster_id",
         show_meds = TRUE, format_vals = TRUE, digits = 3)
topTable(ds_t6, top_n = 55, order_by = "cluster_id",
         show_meds = TRUE, format_vals = TRUE, digits = 3)


res_DA_c10 <- topTable(ds_c10, all = TRUE, show_logFC = T)
table(res_DA_c10$p_adj <= 0.05)

res_DA_dod <- topTable(ds_dod, all = TRUE, show_logFC = T)
table(res_DA_dod$p_adj <= 0.05)

res_DA_t6 <- topTable(ds_t6, all = TRUE, show_logFC = T)
table(res_DA_t6$p_adj <= 0.05)

sigs_t6 <- subset(res_DA_t6, p_adj<=0.05)
sigs_dod <- subset(res_DA_dod, p_adj<=0.05)

spider_sigs_t6 <- data.frame(sigs_t6$cluster_id, sigs_t6$marker_id, sigs_t6$logFC)
spider_sigs_t6$direction <- ifelse(sigs_t6$logFC>=0, "up", "down")

trimmed_spider_sigs_t6 <- subset(spider_sigs_t6, 2^spider_sigs_t6$sigs_t6.logFC>= 1.1 |  2^spider_sigs_t6$sigs_t6.logFC<=0.9)

#trimmed_spider_sigs_t6$sigs_t6.cluster_id <- gsub("gamma delta", expression(alpha) , trimmed_spider_sigs_t6$sigs_t6.cluster_id)

(ds_limma_t6_heatmap <- ggplot(trimmed_spider_sigs_t6, aes(x= sigs_t6.cluster_id, y=sigs_t6.marker_id, label=round(2^sigs_t6.logFC, digits = 2)))+
  #geom_tile(aes(fill=rescale(sigs_t6.logFC, to=c(-5, 5))))+
  geom_tile(aes(fill=sigs_t6.logFC*5))+
  geom_text(parse = TRUE)+
  scale_fill_gradient2(low = "blue", high="red", midpoint = 0)+
  scale_x_discrete(labels=c("gamma delta"=expression(paste(gamma, delta))))+
  # scale_fill_gradientn(colours = c("blue", "white", "red"),
  #                      values = c(min(spider_sigs_t6$sigs_t6.logFC), 0, 0.000001, max(spider_sigs_t6$sigs_t6.logFC)),
  #                      guide = guide_colourbar(nbin = 1000)
  #                       )+
  theme_minimal()+
  ggtitle("Fold Change in Intensity of Marker Expression\nper Cluster at T6 relative to Baseline")+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=14),
        axis.text.x = element_text(hjust=1, angle = 60),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))
)


ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/diffcyt/ds_limma/timepoint_batch_t6_heatmap.png", ds_limma_t6_heatmap, width = 7, height = 7)

plotDiffHeatmap(merged_daf, ds_dod, top_n = 5, order = TRUE,    
                th = FDR_cutoff, normalize = TRUE, hm1 = FALSE)   


plotDiffHeatmap(merged_daf, ds_t6, top_n = 55, order = TRUE,    
                th = FDR_cutoff, normalize = TRUE, hm1 = FALSE)   


table(sigs_t6$cluster_id)
table(sigs_dod$cluster_id)


t6_table <- 

# table(sigs_dod$marker_id)[order(table(sigs_dod$marker_id), decreasing = T)]
# Tbet  CD38 HLADR   PD1  CD25  BCL2   GZB  ICOS  CD27 CTLA4  CD28  Ki67  CD95 
# 6     3     2     2     2     1     1     1     1     1     1     1     0
# table(sigs_t6$marker_id)[order(table(sigs_t6$marker_id), decreasing = T)]
# 
# CD38  CD95  Tbet CTLA4  BCL2  Ki67  CD27 HLADR  ICOS   PD1   GZB  CD28  CD25 
# 25    22    22    18    13    13    12    10     9     8     7     7     6 



# table(sigs_dod$cluster_id)
# 
# activated  CD4 CM                 activated  CD8 EM 
# 0                                 0 
# activated  MAIT                 activated  Treg EM 
# 4                                 0 
# activated  Vd2+          activated HLADR+ DN TEMRA 
# 2                                 0 
# activated PD1+ CD4 EM activated PD1+CD27-HLADR+  CD4 EM 
# 0                                 2 
# activated PD1+HLADR+ CD4 EM       activated PD1+HLADR+ CD8 EM 
# 0                                 0 
# CD27- CD4 EM               CD27- CD8 Effectors 
# 0                                 1 
# CD27- CD8 TEMRA                  CD27- Vd1+ TEMRA 
# 1                                 2 
# CD27-CD28-  CD4 EM                CD28- CD8 Effector 
# 0                                 1 
# CD4 CM                            CD4 EM 
# 0                                 1 
# CD4 Naïve                       CD45RA+ DN  
# 1                                 0 
# CD57+CD27-CD28-  CD4 EM                     CD8 Effectors 
# 0                                 0 
# CD8 Naive                         CD8 TEMRA 
# 2                                 1 
# CLA+ CD4 EM                       PD1+ CD8 EM 
# 0                                 2 
# PD1+CD57+ CD4 EM      PD1+CD57+HLADR+ CD8 Effector 
# 0                                 0 
# resting  MAIT                   resting  Treg EM 
# 0                                 1 
# resting  Vd2+         resting HLADR+CLA+ Treg EM 
# 1                                 0 



#table(sigs_t6$cluster_id)
# activated  CD4 CM                 activated  CD8 EM 
# 9                                 8 
# activated  MAIT                 activated  Treg EM 
# 0                                 5 
# activated  Vd2+          activated HLADR+ DN TEMRA 
# 4                                 1 
# activated PD1+ CD4 EM activated PD1+CD27-HLADR+  CD4 EM 
# 7                                 5 
# activated PD1+HLADR+ CD4 EM       activated PD1+HLADR+ CD8 EM 
# 6                                 7 
# CD27- CD4 EM               CD27- CD8 Effectors 
# 6                                 8 
# CD27- CD8 TEMRA                  CD27- Vd1+ TEMRA 
# 4                                 6 
# CD27-CD28-  CD4 EM                CD28- CD8 Effector 
# 3                                 4 
# CD4 CM                            CD4 EM 
# 2                                 7 
# CD4 Naïve                       CD45RA+ DN  
# 3                                 1 
# CD57+CD27-CD28-  CD4 EM                     CD8 Effectors 
# 4                                 8 
# CD8 Naive                         CD8 TEMRA 
# 2                                 6 
# CLA+ CD4 EM                       PD1+ CD8 EM 
# 3                                 6 
# PD1+CD57+ CD4 EM      PD1+CD57+HLADR+ CD8 Effector 
# 6                                 9 
# resting  MAIT                   resting  Treg EM 
# 8                                 8 
# resting  Vd2+         resting HLADR+CLA+ Treg EM 
# 9                                 7 
