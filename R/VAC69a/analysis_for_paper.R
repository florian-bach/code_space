library(readxl) 
library(CATALYST)
library(flowCore)
library(ggplot2)
library(RColorBrewer)
library(SingleCellExperiment)
library(CytoNorm)
library(diffcyt)
library(tidyr)
library(dplyr)
library(cowplot)
library(scater)
library(purrr)


#extra functions defined here ####
`%!in%` = Negate(`%in%`)

# color palettes defined here ####
#this is how you extract colors from an existing plot object... the order is slightly
# bongled up though...

# (plsm <- ggplot(iris, aes(x=iris$Sepal.Length, y=iris$Sepal.Width, color=iris$Petal.Width))+
#    geom_point()+
#    scale_color_gradientn(colors = inferno_mega_lite))
# 
# g <- ggplot_build(plsm)
# paste(unique(g$data[[1]]["colour"]))


myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
red_palette <- c(myPalette(100), rep(myPalette(100)[100], 200))

sc <- scale_colour_gradientn(colours = red_palette, limits=c(0, 8))

magma <- c("#FCFFB2","#FCDF96","#FBC17D","#FBA368","#FA8657","#F66B4D","#ED504A",
           "#E03B50","#C92D59","#B02363","#981D69","#81176D","#6B116F","#57096E","#43006A",
           "#300060","#1E0848","#110B2D","#080616","#000005")

magma_lite <- c("#FCFFB2","#FCDF96","#FBC17D","#FBA368","#FA8657","#F66B4D","#ED504A",
                "#E03B50","#C92D59","#B02363","#981D69","#81176D","#6B116F","#57096E","#43006A",
                "#300060","#1E0848","#110B2D")


plasma <- c("#2B078E", "#4F049B", "#3E0595", "#0D0887", "#5E02A2", "#6E01A7", "#D45470",
            "#DD6066", "#CB497A","#E56D5E", "#A92593", "#B6308C", "#F0894E", "#C03D83", "#EB7B56",
            "#F0F921", "#F69644", "#FBB434", "#FAC631", "#FBA339", "#F4E828","#F7D72D")

viridis <- c("#45155E", "#452F73", "#462369", "#440154", "#433B7E", "#414687", "#299A87",
             "#25A485", "#2B9089", "#36AD7F", "#30728D", "#2B7C8D", "#5DBE6B", "#2C868B", "#4CB575", "#FDE725",
             "#6BC760", "#93D54C", "#B0DA45", "#78CF54", "#E4E332", "#CBDF3C")

inferno <- c("#16071C", "#2C0E42", "#200D2E", "#000004", "#380D57", "#460B68", "#C84449",
             "#D74D3F", "#B93B53", "#E15D37", "#8A2267", "#9A2964", "#EF802B", "#AA325B",
             "#E86F32", "#FCFFA4", "#F59121", "#FEB431", "#FFC751", "#FBA210", "#FFEC89",
             "#FFDA6D")
inferno_lite <- c("#000004", "#380D57", "#460B68","#8A2267", "#9A2964",
                  "#AA325B","#B93B53","#C84449", "#D74D3F", "#E15D37",
                  "#E86F32","#EF802B", "#F59121", "#FBA210", "#FEB431", "#FFC751",
                  "#FFDA6D", "#FFEC89", "#FCFFA4")
#scales::show_col(inferno_lite)

inferno_mega_lite <-inferno_lite[seq(1,19, by=3)]
inferno_mega_lite <- c("#000004", "#C84449", "#E15D37", "#EF802B", "#FFC751", "#FCFFA4")
smooth_inferno_lite <- colorRampPalette(inferno_lite)


### read in metadata etc. ####
#setwd("C:/Users/Florian/PhD/cytof/vac63c/t_cells/fcs")
setwd("/Users/s1249052/PhD/cytof/vac69a/T_cells_only/fcs")

md <- read.csv("meta_data.csv", header=T) 


## md has to have particular properties: file_name=NAMED LIST (chr), ID and everything else in factors
## reorder so that it's neat 
md$timepoint <- factor(md$timepoint, levels = c("Baseline", "C8", "C10", "C12", "DoD", "T6"))
md <- md[
  with(md, order(md$volunteer, md$timepoint)),
  ]
md$file_name <- as.character(md$file_name)


#change number of volunteers/timepoints here
md <- dplyr::filter(md, timepoint %in% c("Baseline", "C10" , "DoD", "T6"))


# select whose files to import
fcs_files <- grep("fcs", list.files(), value = T)

primaries <- c()

#change number of volunteers here
timepoints <- c("C-1", "C_10", "DoD", "T_6")

for(i in timepoints){
  one_moment <- grep(i, fcs_files, value=T)
  primaries <- c(timepoints, one_moment)
}

### read in flowfiles using flowCore
vac69a <- read.flowSet(md$file_name)


### define panel using first fcs file in flowset, then clean up

# panel <- data.frame("fcs_colname"=colnames(vac69a[[1]]), "antigen"=c("Time", markernames(vac69a[[1]])))
# panel$antigen <- gsub("_", "", panel$antigen)
# 
# panel$antigen <- ifelse(nchar(panel$antigen)>6,
#                         substr(panel$antigen, 6, nchar(panel$antigen)),
#                         panel$antigen
# )
# 
# panel$antigen[3] <- "CD45"
# 
# panel$antigen <- gsub("-", "", panel$antigen)
# panel$antigen <- gsub(".", "", panel$antigen, fixed=T)
# 
# write.csv(panel, "VAC69_PANEL.CSV")
panel <- read.csv("VAC69_PANEL.CSV", header=T)
####

### CATALYST ####
#construct daFrame #
## md has to have particular properties: file_name=NAMED LIST (chr), ID and everything else in factors
daf <- prepData(vac69a, panel, md, md_cols =
                  list(file = "file_name", id = "sample_id", factors = c("timepoint", "batch", "volunteer")))


# this gets rid of barcoding and quality control channels so clustering is easier
proper_channels <- colnames(daf)[c(3,13:14,22:56,62,64)]
# get rid of CD45, CD14, CD20, CD3
t_cell_channels <- proper_channels[c(-1, -2, -10, -32)]

t_cell_channels <- c("CD4",
                     "CD8",
                     "Vd2",
                     "Va72",
                     "CD38",
                     "HLADR",
                     "ICOS",
                     "CD28",
                     "PD1",
                     "TIM3",
                     "CD95",
                     "BCL2",
                     "CD27",
                     "Perforin",
                     "GZB",
                     "CX3CR1",
                     "Tbet",
                     "CTLA4",
                     "Ki67",
                     "CD127",
                     "IntegrinB7",
                     "CD56",
                     "CD16",
                     "CD161",
                     "CD49d",
                     "CD103",
                     "CD25",
                     "FoxP3",
                     "CD39",
                     "CLA",
                     "CXCR5",
                     "CD57",
                     "CD45RA",
                     "CD45RO",
                     "CCR7")

refined_markers <- c("CD4",
                     "CD8",
                     "Vd2",
                     "Va72",
                     "CD38",
                     "HLADR",
                     "ICOS",
                     "CD28",
                     "PD1",
                     #"TIM3",
                     "CD95",
                     "BCL2",
                     "CD27",
                     "Perforin",
                     "GZB",
                     #"CX3CR1",
                     "Tbet",
                     "CTLA4",
                     "Ki67",
                     "CD127",
                     #"IntegrinB7",
                     #"CD56",
                     #"CD16",
                     "CD161",
                     #"CD49d",
                     #"CD103",
                     "CD25",
                     "FoxP3",
                     "CD39",
                     #"CLA",
                     #"CXCR5",
                     "CD57",
                     "CD45RA",
                     "CD45RO",
                     "CCR7")

# clustering ####
set.seed(1234);daf <- cluster(daf, features = refined_markers, xdim = 10, ydim = 10, maxK = 45, seed = 1234)




#start <- Sys.time(); daf <- runUMAP(daf, exprs_values = "exprs", feature_set=refined_markers); end <- Sys.time()
#duration <- round(end - start) # ~12-23 min
#print(c("UMAP took ", duration, " minutes to run on 861161 cells"), sep='')

# start <- Sys.time(); daf <- runTSNE(daf, exprs_values = "exprs", feature_set=refined_markers); end <- Sys.time()
# duration <- round(end - start)
# print(c("tSNE took ", duration, " minutes to run on 861161 cells"), sep='')
# 
# 

# print(c("together, UMAP and tSNE took ", duration, " minutes to run on 10,000 cells"), sep='')
# 
# th35_plot <- plotDR(daf, "UMAP", color_by = "meta35")
# tfh_plot <- plotDR(daf, "UMAP", color_by = "CD4")


plotClusterHeatmap(daf, hm2 = NULL,
                   #m = "meta35",
                   k = "meta45",
                   cluster_anno = TRUE,
                   draw_freqs = TRUE,
                   scale=T, 
                   palette=inferno_lite
                   )



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

####  merge clusters 


merging_table1 <- data.frame(read_excel("flo_cluster_merging1.xlsx"))

merging_table1$new_cluster <- factor(merging_table1$new_cluster)


merged_daf<- mergeClusters(daf, k = "meta45", table = merging_table1, id = "flo_merge")



plotClusterHeatmap(merged_daf, hm2 = NULL,
                   k = "flo_merge",
                   cluster_anno = TRUE,
                   draw_freqs = TRUE,
                   scale=T,
                   palette=inferno_lite)


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

contrast_baseline <- createContrast(c(1, rep(0, 8)))
contrast_c10 <- createContrast(c(c(0, 1), rep(0,7)))
contrast_dod<- createContrast(c(c(0, 0, 1), rep(0,6)))
contrast_t6 <- createContrast(c(c(1, 0, 0, -1), rep(0,5)))


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
# 23    15
table(rowData(da_t6$res)$p_adj < FDR_cutoff)
# FALSE  TRUE 
# 23    15 

plotDiffHeatmap(merged_daf, da_baseline, th = FDR_cutoff, normalize = TRUE, hm1 = T)
plotDiffHeatmap(merged_daf, da_c10, th = FDR_cutoff, normalize = TRUE, hm1 = T)
plotDiffHeatmap(merged_daf, da_dod, th = FDR_cutoff, normalize = TRUE, hm1 = F)
plotDiffHeatmap(merged_daf, da_t6, th = FDR_cutoff, normalize = TRUE, hm1 = F)


# edger with two timepoints for pairwise comparisons#
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
                             cols_fixed = "timepoint",
                             cols_random =c("volunteer","sample_id"))

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
                     formula = da_formula2,
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
# # 29     9



# it only picks up on conserved stuff, but that makes for an easier story, 8/9 clusters are up..


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
# 15     20 

plotDiffHeatmap(daf, da_baseline_vol, th = FDR_cutoff, normalize = TRUE, hm1 = T)
plotDiffHeatmap(daf, da_c10_vol, th = FDR_cutoff, normalize = TRUE, hm1 = T)
plotDiffHeatmap(merged_daf, da_dod_vol, th = FDR_cutoff, normalize = TRUE, hm1 = F)
plotDiffHeatmap(merged_daf, da_t6_vol, th = FDR_cutoff, normalize = TRUE, hm1 = F)

# adjusting for individual differences with a random effect includes some spurious looking results, but 
# also includes some more that i feel are missing if there's no attempt to control for individual identity
# adding a fixed effect somehow breaks it?? might be worth running lme4 oldschool to check what's
# going on, maybe speak to some IEB people about this again...



# topTable(da_t6_vol, show_counts = T)
dod_vol <- data.frame(topTable(pair_base_dod, all=T, show_counts = T))
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


meta_up_t6 <- filterSCE(daf, k = "meta35", cluster_id %in% up_t6)
CATALYST::plotMDS(daf, color_by = "timepoint") 

# [1] "counts_V02_Baseline" "counts_V02_Baseline" "counts_V02_Baseline" "counts_V02_Baseline"
# [5] "counts_V02_Baseline" "counts_V02_Baseline"













# plotClusterExprs(daf100, k = "meta15", markers = "state")
plotAbundances(daf, k = "meta15", by = "cluster_id", shape="volunteer")

merging_table1 <- data.frame(read_excel("flo_cluster_merging1.xlsx"))
merging_table1$original_cluster <- as.character(merging_table1$original_cluster)
merging_table1$new_cluster <- factor(merging_table1$new_cluster, levels = c("CD4+_Naive",
                                                                            "Treg",
                                                                            "CD4+_Effectors_CD57-ICOS+",
                                                                            "CD4+_Effectors_CD57intGZBint",
                                                                            "CD4+_Effectors_CD57+GZB+", 
                                                                            "CD4+_Effectors_CD27-HLADR+",
                                                                            "CD8+_Naive",
                                                                            "CD8+_EM_PD1+",
                                                                            "CD8+EffectorsCD57+Tbet+",
                                                                            "CD8+_Effectors_CD57-CD45RA+CD45RO-",
                                                                            "MAIT",
                                                                            "gd_CD161-CD28-",
                                                                            "gd_CD161+_CD28+",
                                                                            "DN",
                                                                            "trash"
))

merging_table1 <- data.frame(read_excel("flo_cluster_merging1.xlsx"))
daf100 <- mergeClusters(daf100, k = "meta15", table = merging_table1, id = "merging1")

plotDR(daf100, "UMAP", color_by="CD4")
####


### differential analysis using pairwise comparisons implemented in diffcyt ####
ei <- metadata(daf100)$experiment_info

da_formula1 <- createFormula(ei, cols_fixed = "timepoint",
                             cols_random = c("sample_id"))

(da_formula2 <- createFormula(ei,
                              cols_fixed = c("timepoint", "volunteer"), 
                              cols_random = "sample_id"))


#try som100; run colnames(metadata(daf100)$cluster_codes) to see everything

#design <- createDesignMatrix(ei, c("timepoint", "sample_id", "volunteer"))

levels(ei$timepoint)
# [1] "Baseline" "C10"      "DoD"      "T6"   

# try putting baserline level as -1 rather than 0...
contrast_baseline <- createContrast(c(1, rep(0, 3)))
contrast_c10 <- createContrast(c(0,1,0,0))
contrast_dod <- createContrast(c(0,0,1,0))
contrast_t6 <- createContrast(c(rep(0, 3), 1))
# nrow(contrast) == ncol(design)
# data.frame(parameters = colnames(design), contrast)

da_baseline <- diffcyt(daf100,
                       formula = da_formula1,
                       contrast = contrast_baseline,
                       analysis_type = "DA",
                       method_DA = "diffcyt-DA-GLMM",
                       clustering_to_use = "som100",
                       verbose = F)


da_c10 <- diffcyt(daf100,
                  formula = da_formula1,
                  #design = design,
                  contrast = contrast_c10,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-GLMM",
                  clustering_to_use = "som100",
                  verbose = F)

da_dod <- diffcyt(daf100,
                  formula = da_formula1,
                  #design = design,
                  contrast = contrast_dod,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-GLMM",
                  clustering_to_use = "som100",
                  verbose = F)


da_t6 <- diffcyt(daf100,
                 formula = da_formula1,
                 #design = design,
                 contrast = contrast_t6,
                 analysis_type = "DA",
                 method_DA = "diffcyt-DA-GLMM",
                 clustering_to_use = "som100",
                 verbose = F)


da_baseline_vol <- diffcyt(daf100,
                           formula = da_formula2,
                           contrast = contrast_baseline,
                           analysis_type = "DA",
                           method_DA = "diffcyt-DA-GLMM",
                           clustering_to_use = "som100",
                           verbose = F)


da_c10_vol <- diffcyt(daf100,
                      formula = da_formula2,
                      contrast = contrast_c10,
                      analysis_type = "DA",
                      method_DA = "diffcyt-DA-GLMM",
                      clustering_to_use = "som100",
                      verbose = F)

da_dod_vol <- diffcyt(daf100,
                      formula = da_formula2,
                      contrast = contrast_dod,
                      analysis_type = "DA",
                      method_DA = "diffcyt-DA-GLMM",
                      clustering_to_use = "som100",
                      verbose = F)


da_t6_vol <- diffcyt(daf100,
                     formula = da_formula2,
                     contrast = contrast_t6,
                     analysis_type = "DA",
                     method_DA = "diffcyt-DA-GLMM",
                     clustering_to_use = "som100",
                     verbose = F)


FDR_cutoff <- 0.05

plotDiffHeatmap(daf100, da_baseline, th = FDR_cutoff, normalize = TRUE, hm1 = T)
plotDiffHeatmap(daf100, da_c10, th = FDR_cutoff, normalize = TRUE, hm1 = T)
plotDiffHeatmap(daf100, da_dod, th = FDR_cutoff, normalize = TRUE, hm1 = T)
plotDiffHeatmap(daf100, da_t6, th = FDR_cutoff, normalize = TRUE, hm1 = T)



#### metacluster
table(rowData(da_baseline$res)$p_adj < FDR_cutoff)
# FALSE  TRUE 
# 1    14 
table(rowData(da_c10$res)$p_adj < FDR_cutoff)
# FALSE 
# 15 
table(rowData(da_dod$res)$p_adj < FDR_cutoff)
# FALSE FALSE  TRUE 
# 15     94     6
table(rowData(da_t6$res)$p_adj < FDR_cutoff)
# FALSE  TRUE FALSE  TRUE 
# 13     2     88    12

####

# export data ####

### wrapper function extracClusters doesn't work somehow? using the source code for it's easy enough to
### replicate and it seems to work just fin

cluster_ids <- cluster_ids(daf100, "meta15")
clusters <- c("2", "14")
sig_t6_daf <- daf100[cluster_ids %in% clusters,]


#### reclustering on t6 upregulated cells ####
sig_t6_daf <- runDR(sig_t6_daf, "UMAP", rows_to_use = 2000)

bcl2_plot <- plotDR(sig_t6_daf, "UMAP", color_by="BCL2", facet="timepoint")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.title = element_blank())+
  sc

cd38_plot <- plotDR(sig_t6_daf, "UMAP", color_by="CD38", facet="timepoint")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.title = element_blank())+
  sc

cluster_plot <- plotDR(sig_t6_daf, "UMAP", color_by="som100", facet="timepoint")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.title = element_blank())

plotAbundances(sig_t6_daf, k = "som100", by = "cluster_id", shape="volunteer")

plotClusterHeatmap(sig_t6_daf,
                   hm2 = 'abundances',
                   #k = "meta15",
                   k = "som100", 
                   m = "meta15",
                   cluster_anno = T                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               ,
                   draw_freqs = F)
#split_by='timepoint')


panel <- data.frame("fcs_colname"=colnames(sig_t6[[1]]), "antigen"=c("Time", markernames(sig_t6[[1]])))
panel$antigen <- gsub("_", "", panel$antigen)

panel$antigen <- ifelse(nchar(panel$antigen)>6,
                        substr(panel$antigen, 6, nchar(panel$antigen)),
                        panel$antigen
)

panel$antigen[3] <- "CD45"

panel$antigen <- gsub("-", "", panel$antigen)
panel$antigen <- gsub(".", "", panel$antigen, fixed=T)










plotMedExprs(daf100, k = "meta15", facet = "antigen", shape_by = "volunteer", group_by="timepoint")




dafplotDiffHeatmap(daf100, da_baseline_vol, th = FDR_cutoff, normalize = TRUE, hm1 = T)
plotDiffHeatmap(daf100, da_c10_vol, th = FDR_cutoff, normalize = TRUE, hm1 = T)
plotDiffHeatmap(daf100, da_dod_vol, th = FDR_cutoff, normalize = TRUE, hm1 = T, top_n = 55)
plotDiffHeatmap(daf100, da_t6_vol, th = FDR_cutoff, normalize = TRUE, hm1 = T, top_n = 55)

table(rowData(da_baseline_vol$res)$p_adj < FDR_cutoff)
table(rowData(da_c10_vol$res)$p_adj < FDR_cutoff)
table(rowData(da_dod_vol$res)$p_adj < FDR_cutoff)
table(rowData(da_t6_vol$res)$p_adj < FDR_cutoff)


##












#plotClusterHeatmap(daf196, hm2 = NULL, k = "meta25", m = NULL, cluster_anno = TRUE, draw_freqs = TRUE) 

start_time <- Sys.time()

daf100 <- runDR(daf100, "TSNE", rows_to_use = 1000)
daf100 <- runDR(daf100, "UMAP", rows_to_use = 5000)

daf196 <- runDR(daf196, "TSNE", rows_to_use = 1000)
daf196 <- runDR(daf196, "UMAP", rows_to_use = 5000)

end_time <- Sys.time() #8 min

p1 <- plotDR(daf100, "UMAP", color_by="meta15", facet=c("timepoint"))+ theme(legend.position = "none")

p2 <- plotDR(daf196, "UMAP", color_by="meta25", facet=c("timepoint"))+ theme(legend.position = "none")


, nrow = 2, rel_heights = c(5, 5, 2))



myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
red_palette <- c(myPalette(100), rep(myPalette(100)[100], 200))

sc <- scale_colour_gradientn(colours = red_palette, limits=c(0, 8))


(cd38_plot <-  plotDR(daf100, "UMAP", color_by="meta12", facet="timepoint")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none", 
          axis.title = element_blank())#+sc
)



cd38_plot$layers[[1]]$aes_params$size <- 0.1; cd38_plot




density_plot <- plotDR(daf100, "UMAP", color_by="CCR7", facet=c("volunteer", "timepoint"))+
  stat_density2d(bins=100, size=0.11, colour="maroon")

density_plot$layers[[1]] <- NULL
(density_plot <- density_plot+xlim(c(-15, 8))+ylim(c(-10,12))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.title = element_blank()))





#dot plot shenanigans####


bcl2_plot <- plotDR(daf100, "UMAP", color_by="BCL2", facet="timepoint")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.title = element_blank())+
  sc

cd38_plot <- plotDR(daf100, "UMAP", color_by="CD38", facet="timepoint")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.title = element_blank())+
  sc

cluster_plot <- plotDR(daf100, "UMAP", color_by="meta15", facet="timepoint")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.title = element_blank())

bcl2_plot$layers[[1]]$aes_params$size <- 0.1
cd38_plot$layers[[1]]$aes_params$size <- 0.1
cluster_plot$layers[[1]]$aes_params$size <- 0.1

#triple_big <- plot_grid(cd38_plot, bcl2_plot, cluster_plot, nrow = 3)
triple_big <- gridExtra::grid.arrange(cd38_plot, bcl2_plot, cluster_plot, nrow = 3)

ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/big_dot_plot.pdf", triple_big)


bcl2_plot$layers[[1]]$aes_params$size <- 0.01
cd38_plot$layers[[1]]$aes_params$size <- 0.01
cluster_plot$layers[[1]]$aes_params$size <- 0.01

triple_small <- gridExtra::grid.arrange(cd38_plot, bcl2_plot, cluster_plot, nrow = 3)


ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/small_dot_plot.pdf", triple_small)


ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/big_cd38_plot.pdf", cd38_plot)

#end####
