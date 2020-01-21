library(readxl) 
library(CATALYST)
library(flowCore)
library(ggplot2)
library(RColorBrewer)
library(SingleCellExperiment)
library(cowplot)
library(gridExtra)
library(diffcyt)

`%!in%` = Negate(`%in%`)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
red_palette <- c(myPalette(100), rep(myPalette(100)[100], 200))

sc <- scale_colour_gradientn(colours = red_palette, limits=c(0, 8))

#start_time <- Sys.time()
###


# in order to create a daFrame, the single cell experiment, files, metadata and panel information need to be assembled and passed to
# CATALYST


### read in metadata etc. ####
setwd("C:/Users/Florian/PhD/cytof/vac63c/t_cells/fcs")
#setwd("/Users/s1249052/PhD/cytof/vac69a/T_cells_only/fcs")

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
timepoints <- c("C-1", "C_12", "DoD", "T_6")

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
daf <- daFrame(vac69a, panel, md, md_cols =
                 list(file = "file_name", id = "sample_id", factors = c("timepoint", "batch", "volunteer")))

# this gets rid of barcoding and quality control channels so clustering is easier
proper_channels <- colnames(daf)[c(3,13:14,22:56,62,64)]
# get rid of CD45, CD14, CD20, CD3
t_cell_channels <- proper_channels[c(-1, -2, -10, -32)]

marker_levels <- c("CD4",
                   "CD8",
                   "CX3CR1",
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

t_cell_channels <- marker_levels


### metaclustering merged central memory and naive cells...
daf100 <- cluster(daf, cols_to_use = t_cell_channels, xdim = 10, ydim = 10, maxK = 15, seed = 1234)

# daf100_delta <- metadata(daf100)$delta_area

#### get rid of weird super positive events ####
cluster_ids <- cluster_ids(daf100, "meta15")
clusters <- "15"
daf100 <- daf100[cluster_ids %!in% clusters,]####

plotClusterHeatmap(daf100,
                   hm2 = 'TIM3',
                   #k = "meta15",
                   k = "meta15", 
                   #m = "meta15",
                   cluster_anno = T                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               ,
                   draw_freqs = T,
                   scale=T)
                   
#split_by='timepoint')
               

plotClusterExprs(daf100, k = "meta15", markers = "state")
plotAbundances(daf100, k = "meta15", by = "cluster_id", shape="volunteer")

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
