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

setwd("/home/florian/PhD/cytof/vac63c/t_cells/experiment_279034_files")
# setwd("C:/Users/Florian/PhD/cytof/vac63c/t_cells/experiment_279034_files")
# setwd("/Users/s1249052/PhD/cytof/vac63c/fcs_files/T\ cells/")

#### how the metadata file was made ####
md <- read.delim("experiment_279034_annotations.tsv")
md <- dplyr::select(md, FCS.Filename, Timepoints, Individuals)
colnames(md) <- c("file_name", "timepoint", "volunteer")

md$timepoint <- gsub("C-1", "Baseline", md$timepoint)
md$timepoint <- gsub("+", "", md$timepoint, fixed = T)

batch <- stringr::str_match(md$file_name, "b[0-9]*")[, 1]
batch <- ifelse(is.na(batch)==T, "b3", batch)

sample_id <- paste(md$timepoint, "_v", md$volunteer, sep='')

#md <- read_excel("meta_data.xlsx") 
vac63c_metadata <- data.frame("file_name"=md$file_name, "sample_id"=sample_id, "timepoint"=md$timepoint, "batch"=batch, "volunteer"=paste("v", md$volunteer, sep=''))
write.csv(vac63c_metadata, "vac63c_metadata.csv") 

#### let's go ####
md <- read.csv("vac63c_metadata.csv")

#change number of volunteers here

md <- dplyr::filter(md, volunteer %in% c("v315", "v313", "v320"))

md$timepoint <- factor(md$timepoint, levels = c("Baseline", "DoD", "T6", "C45"))
md$sample_id <- factor(md$sample_id, levels = md$sample_id[order(md$timepoint)])

# select whose files to import

#fcs_files <- grep("fcs", list.files("/Users/s1249052/PhD/cytof/vac63c/fcs_files/T\ cells/"), value = T)
#fcs_files <- grep("fcs", list.files("C:/Users/Florian/PhD/cytof/vac63c/t_cells/experiment_279034_files"), value = T)
fcs_files <- grep("fcs", list.files("/home/florian/PhD/cytof/vac63c/t_cells/experiment_279034_files"), value = T)


primaries <- c()

#change number of volunteers here
people <- c("315", "313", "320")

for(i in people){
  one_person <- grep(i, fcs_files, value=T)
  primaries <- c(primaries, one_person)
}

vac63c_primaries <- read.flowSet(primaries)



## antigen colnames CANNOT start with a digit, it won't parse

#### how the panel was made ####
try <- vac63c_primaries[[1]]
panel <- data.frame("fcs_colname"=colnames(try), "antigen"=c("Time", markernames(try)))
panel$antigen <- gsub("_", "", panel$antigen)

panel$antigen <- ifelse(nchar(panel$antigen)>6,
                                       substr(panel$antigen, 6, nchar(panel$antigen)),
                                       panel$antigen
                                       )
panel$antigen[3] <- "CD45"

write.csv(panel, "vac63c_panel.csv")
#### ####

panel <- read.csv("vac63c_panel.csv")
# spot check that all panel columns are in the flowSet object
# all(panel$fcs_colname %in% colnames(vac63c_primaries))
# all(colnames(vac63c_primaries) %in% panel$fcs_colname)
# # 
# all(md$file_name %in% sampleNames(vac63c_primaries))
# all(sampleNames(vac63c_primaries) %in% md$file_name)




## md has to have particular properties: file_name=NAMED LIST (chr), ID AND THEN EVERYTHING ELSE IN FACTORS
daf <- daFrame(vac63c_primaries, panel, md, md_cols =
                 list(file = "file_name", id = "sample_id", factors = c("timepoint", "batch", "volunteer")))


# this gets rid of barcoding and quality control channels so clustering is easier
proper_channels <- colnames(daf)[c(3,14:15,23:57,63,65)]
t_cell_channels <- proper_channels[c(-1, -2, -32)]

marker_levels <- c("CD4",
                   "CD8",
                   "Vd2",
                   "Va72",
                   "CD38",
                   "CD69",
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
                   "Tbet",
                   "Eomes",
                   "RORgt",
                   "GATA3",
                   "CTLA4",
                   "Ki67",
                   "CD127",
                   "CD56",
                   "CD16",
                   "CD161",
                   "CD49d",
                   "CD25",
                   "FoxP3",
                   "CD39",
                   "CXCR5",
                   "CX3CR1",
                   "CD57",
                   "CD45RA",
                   "CD45RO",
                   "CCR7")

t_cell_channels <- marker_levels


##
daf <- cluster(daf, cols_to_use = t_cell_channels, xdim = 10, ydim = 10, maxK = 25, seed = 1234)

##metadata(daf)$delta_area

plotMDS(daf, color_by = "timepoint")

plotExprHeatmap(daf, bin_anno = TRUE, row_anno = TRUE)


# run FlowSOM

#plotNRS(daf, color_by="timepoint")

# density plot, no contours
# vd2_plot <- plotDR(daf, "UMAP", color_by="Vd2", facet=c("volunteer","timepoint"))+
#   stat_density2d(aes(fill = ..density..), geom = 'tile', contour=F, bins=100, size=0.2)
# 
# vd2_plot$layers[[1]] <- NULL
# vd2_plot+xlim(c(-7, 15))+ylim(c(-7, 15))+
#   scale_fill_distiller(palette = 'RdPu', trans="reverse")



plotClusterHeatmap(daf, hm2 = "abundances",
                   k = "meta25",
                   m = NULL,
                   cluster_anno = TRUE,
                   draw_freqs = TRUE,
                   scale=F) 




ei <- metadata(daf)$experiment_info

da_formula1 <- createFormula(ei, cols_fixed = "timepoint",
                             cols_random = c("sample_id"))

(da_formula2 <- createFormula(ei,
                              cols_fixed = c("timepoint", "volunteer"), 
                              cols_random = "sample_id"))


#try som100; run colnames(metadata(daf100)$cluster_codes) to see everything

#design <- createDesignMatrix(ei, c("timepoint", "sample_id", "volunteer"))

levels(ei$timepoint)
# [1] "Baseline" "DoD"      "T6"       "C45"  

# try putting baserline level as -1 rather than 0...
contrast_baseline <- createContrast(c(1, rep(0, 3)))
contrast_dod<- createContrast(c(0,1,0,0))
contrast_t6 <- createContrast(c(0,0,1,0))
contrast_c45 <- createContrast(c(rep(0, 3), 1))
# nrow(contrast) == ncol(design)
# data.frame(parameters = colnames(design), contrast)

da_baseline <- diffcyt(daf,
                       formula = da_formula1,
                       contrast = contrast_baseline,
                       analysis_type = "DA",
                       method_DA = "diffcyt-DA-GLMM",
                       clustering_to_use = "meta25",
                       verbose = T)


da_dod <- diffcyt(daf,
                  formula = da_formula1,
                  #design = design,
                  contrast = contrast_dod,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-GLMM",
                  clustering_to_use = "meta25",
                  verbose = F)

da_t6 <- diffcyt(daf,
                  formula = da_formula1,
                  #design = design,
                  contrast = contrast_t6,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-GLMM",
                  clustering_to_use = "meta25",
                  verbose = F)


da_c45 <- diffcyt(daf,
                 formula = da_formula1,
                 #design = design,
                 contrast = contrast_c45,
                 analysis_type = "DA",
                 method_DA = "diffcyt-DA-GLMM",
                 clustering_to_use = "meta25",
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

plotDiffHeatmap(daf, da_baseline, th = FDR_cutoff, normalize = TRUE, hm1 = T, draw_freq=T)
plotDiffHeatmap(daf, da_dod, th = FDR_cutoff, normalize = TRUE, hm1 = T)
plotDiffHeatmap(daf, da_t6, th = FDR_cutoff, normalize = TRUE, hm1 = T, top_n=50)
plotDiffHeatmap(daf, da_c45, th = FDR_cutoff, normalize = TRUE, hm1 = T, top_n=50)


table(rowData(da_baseline$res)$p_adj < FDR_cutoff)
table(rowData(da_dod$res)$p_adj < FDR_cutoff)
table(rowData(da_t6$res)$p_adj < FDR_cutoff)
table(rowData(da_c45$res)$p_adj < FDR_cutoff)

# <3
dod_tbl <- data.frame(topTable(da_dod, show_props = TRUE, format_vals = TRUE, digits = 2))

sig_dod <- dplyr::select(dod_tbl, colnames(dod_tbl)[c(1:6,grep("DoD", colnames(dod_tbl), fixed=T))])
sig_dod_clusters <- subset(rowData(da_dod$res), rowData(da_dod$res)$p_adj < FDR_cutoff)$cluster_id 

sig_dod <- filter(sig_dod, sig_dod$cluster_id %in% sig_dod_clusters)
long_sig_dod <- gather(sig_dod, sample, frequency, colnames(sig_dod)[4:length(colnames(sig_dod))])
long_sig_dod$volunteer <- substr(long_sig_dod$sample, nchar(long_sig_dod$sample)-3, nchar(long_sig_dod$sample))
long_sig_dod$timepoint <- ifelse(grepl("DoD", long_sig_dod$sample)==T, "DoD", "Baseline")

ggplot(long_sig_dod, aes(x=timepoint, y=frequency, group=sample, fill=volunteer))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_minimal()



t6_tbl <- data.frame(topTable(da_t6, show_props = TRUE, format_vals = TRUE, digits = 2))

sig_t6 <- dplyr::select(t6_tbl, colnames(t6_tbl)[c(1:6,grep("T6", colnames(t6_tbl), fixed=T))])
sig_t6_clusters <- subset(rowData(da_t6$res), rowData(da_t6$res)$p_adj < FDR_cutoff)$cluster_id 


sig_t6 <- filter(sig_t6, sig_t6$cluster_id %in% sig_t6_clusters)
long_sig_t6 <- gather(sig_t6, sample, frequency, colnames(sig_t6)[4:length(colnames(sig_t6))])
long_sig_t6$volunteer <- substr(long_sig_t6$sample, nchar(long_sig_t6$sample)-3, nchar(long_sig_t6$sample))
long_sig_t6$timepoint <- ifelse(grepl("T6", long_sig_t6$sample)==T, "T6", "Baseline")


ggplot(long_sig_t6, aes(x=timepoint, y=frequency, group=sample, fill=volunteer))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_minimal()+
  facet_wrap(~long_sig_t6$cluster_id, scales="free")
  
diffcyt::  
  









sigt6 <- subset(rowData(da_t6$res), rowData(da_t6$res)$p_adj < FDR_cutoff) 

sigt6$cluster_id

sig_clusters <- daf[cluster_id %in% as.character(sigt6$cluster_id), "meta25")

cluster_ids <- cluster_ids(daf, "meta25")
sig_clusters <- as.character(sigt6$cluster_id)
sig_t6_daf <- daf[cluster_ids %in% sig_clusters,]




daf <- runDR(daf, "UMAP", rows_to_use = 5000)


table(rowData(da_dod_vol$res)$p_adj < FDR_cutoff)



u <- CATALYST::filterSCE(daf, k = "meta25",
               cluster_id %in% sigt6$cluster_id)

cowplot::plot_grid(
  plotDR(daf, "UMAP", color_by = "BCL2",  facet="timepoint"),
  plotDR(sig_t6_daf, "UMAP", color_by = "BCL2",  facet="timepoint"), nrow = 2)



myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
red_palette <- c(myPalette(100), rep(myPalette(100)[100], 200))

sc <- scale_colour_gradientn(colours = red_palette, limits=c(0, 8))









# density plot, with contours
density_plot <- plotDR(daffy, "UMAP", color_by="CCR7", facet=c("timepoint"))+
  stat_density2d(bins=100, size=0.11, colour="maroon")

density_plot$layers[[1]] <- NULL
density_plot <- density_plot+xlim(c(-15, 8))+ylim(c(-8,15))+
theme(panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(), 
      axis.title = element_blank())


cd38_plot <- plotDR(daf, "UMAP", color_by="CD38", facet="timepoint")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.title = element_blank())+
  sc


cd38_plot$layers[[1]]$aes_params$size <- 1

bcl2_plot <- plotDR(daf, "UMAP", color_by="BCL2", facet="timepoint")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.title = element_blank())+
  sc

ggsave("dens_cd38_bcl2.png", cowplot::plot_grid(density_plot, cd38_plot, bcl2_plot, nrow=3), dpi=400, width=6.4, height = 6.4)



cxcr5_plot <- plotDR(daf, "UMAP", color_by="CXCR5", facet="timepoint")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  sc
  

pd1_plot <- plotDR(daf, "UMAP", color_by="PD1", facet="timepoint")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+sc


icos_plot <- plotDR(daf, "UMAP", color_by="ICOS", facet="timepoint")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  sc

cowplot::plot_grid(cxcr5_plot, pd1_plot, icos_plot, nrow=3)



#### cytoNorm ####

devtools::install_github('saeyslab/CytoNorm')
remotes::install_github("r-lib/rlang")
