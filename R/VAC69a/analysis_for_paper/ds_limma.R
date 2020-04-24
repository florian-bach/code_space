library(diffcyt)
library(tidyr)
library(dplyr)
library(cowplot)
library(SummarizedExperiment)
library(SingleCellExperiment)

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

working_directory <- "~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/"

md <- read.csv(paste(working_directory, "meta_data.csv", sep=''), header=T, stringsAsFactors = F)
md$timepoint <- factor(md$timepoint, levels = c("Baseline", "C8", "C10", "C12", "DoD", "T6"))
md <- md[with(md, order(md$volunteer, md$timepoint)),]
md$file_name <- paste(working_directory, md$file_name, sep='')

#read in panel
panel <- read.csv(paste(working_directory, "VAC69_PANEL.CSV", sep=''), header = T, stringsAsFactors = F)
colnames(panel)[2] <- "marker_name"

### read in flowfiles using flowCore
vac69a <- flowCore::read.flowSet(md$file_name)


#sample.int takes a sample of the specified size from the elements of x using either with or without replacement.
# set.seed(1234); smaller_vac69a <- fsApply(vac69a, function(ff) {
#   idx <- sample.int(nrow(ff), min(sampling_ceiling, nrow(ff)))
#   ff[idx,]  # alt. ff[order(idx),]
# })

proportional=TRUE
event_number=3500
if(isTRUE(proportional))
  assign("downsample", function(fs, event_number){
    flowCore::fsApply(fs, function(ff){
      idx <- sample.int(nrow(ff), nrow(ff)/min(flowCore::fsApply(fs, nrow))*event_number)
      ff[idx,]
    })
  })

if(xor(isFALSE(proportional), is.null(proportional)))
  assign("downsample", function(fs, event_number){
    flowCore::fsApply(fs, function(ff){
      idx <- sample.int(nrow(ff), min(event_number, nrow(ff)))
      ff[idx,]
    })
  })


set.seed(1234); smaller_vac69a <- downsample(vac69a, event_number)



### CATALYST ####
#construct daFrame #
## md has to have particular properties: file_name=NAMED LIST (chr), ID and everything else in factors
sce <- CATALYST::prepData(smaller_vac69a, panel, md,
                          
                          md_cols = list(file = "file_name",
                                         id = "sample_id",
                                         factors = c("timepoint", "batch", "volunteer")
                          ),
                          
                          panel_cols = list(channel = "fcs_colname",
                                            antigen = "marker_name",
                                            class = "marker_class"
                          )
)

typs <- type_markers(sce)
typs <- typs[-match(c("CD20", "CXCR5", "CD103", "TCRgd", "IntegrinB7", "CD56", "CD3", "CD49d"), typs)]

# refined_markers <- c("CD4",
#                      "CD8",
#                      "Vd2",
#                      "Va72",
#                      "CD38",
#                      "HLADR",
#                      "ICOS",
#                      "CD28",
#                      "PD1",
#                      #"TIM3",
#                      "CD95",
#                      "BCL2",
#                      "CD27",
#                      "Perforin",
#                      "GZB",
#                      "CX3CR1",
#                      "Tbet",
#                      "CTLA4",
#                      "Ki67",
#                      "CD127",
#                      #"IntegrinB7",
#                      #"CD56",
#                      #"CD16",
#                      "CD161",
#                      #"CD49d",
#                      #"CD103",
#                      "CD25",
#                      "FoxP3",
#                      "CD39",
#                      "CLA",
#                      #"CXCR5",
#                      "CD57",
#                      "CD45RA",
#                      "CD45RO",
#                      "CCR7")

#write.csv(data.frame(refined_markers), paste(working_directory, "refined_markers.csv", sep = ''),  row.names = F)
refined_markers <- read.csv(paste(working_directory, "refined_markers.csv", sep = ''), stringsAsFactors = F)

# clustering ####
set.seed(123);sce <- CATALYST::cluster(sce, features = refined_markers[,1], xdim = 10, ydim = 10, maxK = 45)

daf <- read_full("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")



merging_table1 <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/merging_table_april2020.csv", header=T, stringsAsFactors = F)

#get rid of spaces at beginning of string
merging_table1$new_cluster <- ifelse(substr(merging_table1$new_cluster, 1, 1)==" ", substr(merging_table1$new_cluster, 2, nchar(merging_table1$new_cluster)), merging_table1$new_cluster)
merging_table1$new_cluster <- ifelse(substr(merging_table1$new_cluster, 1, 1)==" ", substr(merging_table1$new_cluster, 2, nchar(merging_table1$new_cluster)), merging_table1$new_cluster)

merging_table1$new_cluster <- factor(merging_table1$new_cluster)

merged_daf<- mergeClusters(daf, k = "meta45", table = merging_table1, id = "flo_merge")


plotClusterHeatmap(merged_daf, hm2=NULL,
                   k = "flo_merge",
                   #m = "flo_merge",
                   cluster_anno = FALSE,
                   draw_freqs = TRUE,
                   scale = TRUE, 
                   palette=inferno
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

contrast_baseline <- createContrast(c(0, rep(1, 3), rep(0, 5)))
contrast_c10 <- createContrast(c(c(0, 1), rep(0,7)))
contrast_dod<- createContrast(c(c(0, 0, 1), rep(0,6)))
contrast_t6 <- createContrast(c(c(0, 0, 0, 1), rep(0,5)))


# da_baseline <- diffcyt(merged_daf,
#                        design = design,
#                        contrast = contrast_baseline,
#                        analysis_type = "DA",
#                        method_DA = "diffcyt-DA-edgeR",
#                        clustering_to_use = "meta40",
#                        verbose = T)
# 
# da_c10 <- diffcyt(merged_daf,
#                   design = design,
#                   contrast = contrast_c10,
#                   analysis_type = "DA",
#                   method_DA = "diffcyt-DA-edgeR",
#                   clustering_to_use = "meta40",
#                   verbose = T)
# 
# da_dod <- diffcyt(merged_daf,
#                   design = design,
#                   contrast = C,
#                   analysis_type = "DA",
#                   method_DA = "diffcyt-DA-edgeR",
#                   clustering_to_use = "meta40",
#                   verbose = T)
# 
# da_t6 <- diffcyt(merged_daf,
#                  design = design,
#                  contrast = contrast_t6,
#                  analysis_type = "DA",
#                  method_DA = "diffcyt-DA-edgeR",
#                  clustering_to_use = "meta40",
#                  verbose = T)
# 
# # results  ###
# table(rowData(da_baseline$res)$p_adj < FDR_cutoff)
# 
# 
# table(rowData(da_c10$res)$p_adj < FDR_cutoff)
# # FALSE 
# # 38 
# table(rowData(da_dod$res)$p_adj < FDR_cutoff)
# # FALSE  TRUE 
# # 24    8
# table(rowData(da_t6$res)$p_adj < FDR_cutoff)
# # FALSE  TRUE 
# # 23    15 
# 
# plotDiffHeatmap(merged_daf, da_baseline, th = FDR_cutoff, normalize = TRUE, hm1 = T)
# plotDiffHeatmap(merged_daf, da_c10, th = FDR_cutoff, normalize = TRUE, hm1 = T)
# plotDiffHeatmap(merged_daf, da_dod, th = FDR_cutoff, normalize = TRUE, hm1 = T)
# plotDiffHeatmap(merged_daf, da_t6, th = FDR_cutoff, normalize = TRUE, hm1 = T)


states <- names(marker_classes(daf))
states <- states[-c(1:21,
                    match(c("CCR7", "CD127", "CD20", "CXCR5", "CD103", "TCRgd", "TIM3", "IntegrinB7", "CD56", "CD3", "CD49d", "CD45RA", "CD4", "Vd2", "Va72", "CD161", "FoxP3", "CD45RO"), states),
                    53:length(states)
                    )
                 ]

logic <- names(marker_classes(daf)) %in% states

ds_formula2 <- createFormula(ei, cols_fixed = "timepoint", cols_random = "sample_id") 

contrast_dod<- createContrast(c(c(0, 0, 1), rep(0,6)))
contrast_t6 <- createContrast(c(c(0, 0, 0, 1), rep(0,5)))

ds_dod <- diffcyt(merged_daf,                                            
                   design = design, contrast = contrast_dod,                    
                   analysis_type = "DS", method_DS = "diffcyt-DS-limma",       
                   clustering_to_use = "flo_merge", verbose = TRUE, markers_to_test = logic)               


ds_t6 <- diffcyt(merged_daf,                                            
                   design = design, contrast = contrast_t6,                    
                   analysis_type = "DS",  method_DS =  "diffcyt-DS-limma",         
                   clustering_to_use = "flo_merge", verbose = TRUE, markers_to_test = logic)



topTable(ds_dod, top_n = 5, order_by = "cluster_id",
         show_meds = TRUE, format_vals = TRUE, digits = 3)
topTable(ds_t6, top_n = 5, order_by = "cluster_id",
         show_meds = TRUE, format_vals = TRUE, digits = 3)

res_DA_dod <- topTable(ds_dod, all = TRUE)
table(res_DA_dod$p_adj <= 0.05)

res_DA_t6 <- topTable(ds_t6, all = TRUE)
table(res_DA_t6$p_adj <= 0.05)

sigs_t6 <- subset(res_DA_t6, p_adj<=0.05)
sigs_dod <- subset(res_DA_dod, p_adj<=0.05)


plotDiffHeatmap(merged_daf, ds_dod, top_n = 50, order = TRUE,    
                th = FDR_cutoff, normalize = TRUE, hm1 = FALSE)   


table(sigs_t6$cluster_id)
table(sigs_dod$cluster_id)

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
