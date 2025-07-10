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

#setwd("/home/florian/PhD/cytof/vac63c/t_cells/experiment_279034_files")
# setwd("C:/Users/Florian/PhD/cytof/vac63c/t_cells/experiment_279034_files")
setwd("/Users/s1249052/PhD/cytof/vac63c/fcs_files/experiment_279034_files")

#### how the metadata file was made ####
# md <- read.delim("experiment_279034_annotations.tsv")
# md <- dplyr::select(md, FCS.Filename, Timepoints, Individuals)
# colnames(md) <- c("file_name", "timepoint", "volunteer")
# 
# md$timepoint <- gsub("C-1", "Baseline", md$timepoint)
# md$timepoint <- gsub("+", "", md$timepoint, fixed = T)
# 
# batch <- stringr::str_match(md$file_name, "b[0-9]*")[, 1]
# batch <- ifelse(is.na(batch)==T, "b3", batch)
# 
# sample_id <- paste(md$timepoint, "_v", md$volunteer, sep='')
# 
# #md <- read_excel("meta_data.xlsx") 
# vac63c_metadata <- data.frame("file_name"=md$file_name, "sample_id"=sample_id, "timepoint"=md$timepoint, "batch"=batch, "volunteer"=paste("v", md$volunteer, sep=''))
# write.csv(vac63c_metadata, "vac63c_metadata.csv") 

#### let's go ####
md <- read.csv("vac63c_metadata.csv")

#change number of volunteers here

md <- dplyr::filter(md, volunteer %in% c("v315", "v313", "v320"))

md$timepoint <- factor(md$timepoint, levels = c("Baseline", "DoD", "T6", "C45"))
md$sample_id <- factor(md$sample_id, levels = md$sample_id[order(md$timepoint)])

# select whose files to import

fcs_files <- grep("fcs", list.files("/Users/s1249052/PhD/cytof/vac63c/fcs_files/experiment_279034_files"), value = T)
#fcs_files <- grep("fcs", list.files("C:/Users/Florian/PhD/cytof/vac63c/t_cells/experiment_279034_files"), value = T)
#fcs_files <- grep("fcs", list.files("/home/florian/PhD/cytof/vac63c/t_cells/experiment_279034_files"), value = T)


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
# try <- vac63c_primaries[[1]]
# panel <- data.frame("fcs_colname"=colnames(try), "antigen"=c("Time", markernames(try)))
# panel$antigen <- gsub("_", "", panel$antigen)
# 
# panel$antigen <- ifelse(nchar(panel$antigen)>6,
#                                        substr(panel$antigen, 6, nchar(panel$antigen)),
#                                        panel$antigen
#                                        )
# panel$antigen[3] <- "CD45"
# 
# write.csv(panel, "vac63c_panel.csv")
#### make sce object ####

panel <- read.csv("vac63c_panel.csv")
# spot check that all panel columns are in the flowSet object
# all(panel$fcs_colname %in% colnames(vac63c_primaries))
# all(colnames(vac63c_primaries) %in% panel$fcs_colname)
# # 
# all(md$file_name %in% sampleNames(vac63c_primaries))
# all(sampleNames(vac63c_primaries) %in% md$file_name)




## md has to have particular properties: file_name=NAMED LIST (chr), ID AND THEN EVERYTHING ELSE IN FACTORS
daf <- prepData(vac63c_primaries, panel, md, md_cols =
                 list(file = "file_name", id = "sample_id", factors = c("timepoint", "batch", "volunteer")))


# this gets rid of barcoding and quality control channels so clustering is easier
proper_channels <- rownames(daf)[c(3,14:15,23:57,63,65)]
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


#  FLowSOM clustering ####
start <- Sys.time()
daf <- cluster(daf, features = t_cell_channels, xdim = 10, ydim = 10, maxK = 35, seed = 1234)

# UMAP ####
# set.seed(1234)
# daf <- runUMAP(daf, exprs_values = "exprs", feature_set=t_cell_channels)
# plotDR(daf, color_by="meta25")+facet_wrap("timepoint")
# end <- Sys.time()
# print(paste0("clustering & UMAP projecting took ", round(start-end, digits=1), " minutes"))
# 11.6 minutes

metadata(daf)$delta_area

# plotMDS(daf, color_by = "timepoint")
# 
# plotExprHeatmap(daf, bin_anno = TRUE, row_anno = TRUE)


# run FlowSOM

#plotNRS(daf, color_by="timepoint")

# density plot, no contours
# vd2_plot <- plotDR(daf, "UMAP", color_by="Vd2", facet=c("volunteer","timepoint"))+
#   stat_density2d(aes(fill = ..density..), geom = 'tile', contour=F, bins=100, size=0.2)
# 
# vd2_plot$layers[[1]] <- NULL
# vd2_plot+xlim(c(-7, 15))+ylim(c(-7, 15))+
#   scale_fill_distiller(palette = 'RdPu', trans="reverse")



plotClusterHeatmap(daf,draw_freqs = T,
                   k = "meta35",
                   m = NULL,
                   cluster_anno = TRUE,
                  
                   scale=T) 




plotAbundances(daf, k = "meta12", by = "sample_id", group_by = "timepoint")

#### diffcyt ####

ei <- metadata(daf)$experiment_info

da_formula1 <- createFormula(ei, cols_fixed = "timepoint")

da_formula2 <- createFormula(ei,
                              cols_fixed = "timepoint",
                              cols_random = "volunteer")


#try som100; run colnames(metadata(daf100)$cluster_codes) to see everything

design <- createDesignMatrix(ei, c("timepoint", "volunteer"))
data.frame(parameters = colnames(design), contrast_t6)

levels(ei$timepoint)
# [1] "Baseline" "DoD"      "T6"       "C45"  

# try putting baserline level as -1 rather than 0...
contrast_baseline <- createContrast(c(1, rep(0, 3)))
contrast_dod<- createContrast(c(0,1,0,0,0,0))
contrast_t6 <- createContrast(c(0,0,1,0,0,0))
contrast_c45 <- createContrast(c(rep(0, 3), 1, 0, 0))

data.frame(parameters = colnames(design), contrast_c45)


#design <- createDesignMatrix(ei(daf), cols_design = "timepoint")

# nrow(contrast) == ncol(design)
# data.frame(parameters = colnames(design), contrast)


#### models without volunteer as effect ####


da_baseline <- diffcyt(daf,
                       formula = da_formula1,
                       experiment_info = ei,
                       marker_info = panel,
                       contrast = contrast_baseline,
                       analysis_type = "DA",
                       method_DA = "diffcyt-DA-GLMM",
                       clustering_to_use = "meta35",
                       verbose = T)


da_dod <- diffcyt(daf,
                  design = design,
                  contrast = contrast_dod,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "meta35",
                  verbose = F)

da_t6 <- diffcyt(daf,
                  formula = da_formula1,
                  design = design,
                  contrast = contrast_t6,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "meta35",
                  verbose = F)


da_c45 <- diffcyt(daf,
                 design = design,
                 contrast = contrast_c45,
                 analysis_type = "DA",
                 method_DA = "diffcyt-DA-edgeR",
                 clustering_to_use = "meta35",
                 verbose = F)

#### models with volunteer as effect ####
da_baseline_vol <- diffcyt(daf,
                           formula = da_formula2,
                           contrast = contrast_baseline,
                           analysis_type = "DA",
                           method_DA = "diffcyt-DA-GLMM",
                           clustering_to_use = "meta25",
                           verbose = F)

da_dod_vol <- diffcyt(daf,
                      formula = da_formula2,
                      contrast = contrast_dod,
                      analysis_type = "DA",
                      method_DA = "diffcyt-DA-GLMM",
                      clustering_to_use = "meta35",
                      verbose = F)


da_t6_vol <- diffcyt(daf,
                     formula = da_formula2,
                     contrast = contrast_t6,
                     analysis_type = "DA",
                     method_DA = "diffcyt-DA-GLMM",
                     clustering_to_use = "meta35",
                     verbose = F)

da_c45_vol <- diffcyt(daf,
                      formula = da_formula2,
                      contrast = contrast_c45,
                      analysis_type = "DA",
                      method_DA = "diffcyt-DA-GLMM",
                      clustering_to_use = "meta25",
                      verbose = F)


#### diffcyt heatmaps ####
FDR_cutoff <- 0.05

plotDiffHeatmap(daf, da_baseline, th = FDR_cutoff, normalize = TRUE, hm1 = T)
plotDiffHeatmap(daf, da_dod, th = FDR_cutoff, normalize = TRUE, hm1 = T)
#mostly mait cells going away, nothing up
plotDiffHeatmap(daf, da_t6, th = FDR_cutoff, normalize = TRUE, hm1 = T, all=T)
# # UPREGULATED AT T6: 15, 24, 9, 12, 18, 7, 22, 19, 16
plotDiffHeatmap(daf, da_c45, th = FDR_cutoff, normalize = TRUE, hm1 = T, top_n=50)
# 
plotDiffHeatmap(daf, da_baseline_vol, th = FDR_cutoff, normalize = TRUE, hm1 = T)
plotDiffHeatmap(daf, da_dod_vol, th = FDR_cutoff, normalize = TRUE, hm1 = T, all = T)
plotDiffHeatmap(daf, da_t6_vol, th = FDR_cutoff, normalize = TRUE, hm1 = T, all = T)
plotDiffHeatmap(daf, da_c45, th = FDR_cutoff, normalize = TRUE, hm1 = T, top_n=50)

#upregulated when including volunteer: 7, 24, 18, 13, 9, 21
#
# 
table(rowData(da_baseline$res)$p_adj < FDR_cutoff)
table(rowData(da_dod$res)$p_adj < FDR_cutoff)
table(rowData(da_t6$res)$p_adj < FDR_cutoff)
table(rowData(da_c45$res)$p_adj < FDR_cutoff)


table(rowData(da_baseline_vol$res)$p_adj < FDR_cutoff)
table(rowData(da_dod_vol$res)$p_adj < FDR_cutoff)
table(rowData(da_t6_vol$res)$p_adj < FDR_cutoff)
table(rowData(da_c45_vol$res)$p_adj < FDR_cutoff)


#### FILTERING ON UPREGULATED STUFF
sig_t6 <- subset(rowData(da_t6_vol$res), rowData(da_t6_vol$res)$p_adj < FDR_cutoff)
up_t6 <- paste(c(10, 31, 17, 21, 22, 18, 27, 14, 29, 8, 32, 28, 23, 26))

meta_up_t6 <- filterSCE(daf, k = "meta35", cluster_id %in% up_t6)

#### color scheme business ####
#pastel scheme from https://graphicdesign.stackexchange.com/questions/3682/where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d
meta35_color_scheme <- c(
  "#4A3B53","#0AA6D8","#FFB500","#6A3A4C","#456648","#B79762","#3B5DFF","#0CBD66","#4FC601",
  "#D157A0","#34362D","#885578","#8FB0FF","#006FA6","#BA0900","#7A4900","#FFF69F","#FFDBE5",
  "#EEC3FF","#B4A8BD","#3A2465","#7A87A1","#1CE6FF","#6F0062","#0086ED","#6B002C","#C0B9B2",
  "#04F757","#000000","#FF2F80","#FFAA92","#A30059","#FF34FF","#013349","#C2FFED")

names(meta35_color_scheme) <- as.character(seq(1, length(meta35_color_scheme)))

#maximally different 64 colours from https://graphicdesign.stackexchange.com/questions/3682/where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d
# hard_scheme <- c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
#                  "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
#                  "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
#                  "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
#                  "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
#                  "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
#                  "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
#                  "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
#                  "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
#                  "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
#                  "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
#                  "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
#                  "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")
# 
# names(hard_scheme) <- as.character(seq(1, length(hard_scheme)))

#random sampling
# scale1 <- sample(hard_scheme, size = 35);names(scale1) <- as.character(seq(1, length(scale1)))
# scale2 <- sample(hard_scheme, size = 35);names(scale2) <- as.character(seq(1, length(scale2)))
# scale3 <- sample(hard_scheme, size = 35);names(scale3) <- as.character(seq(1, length(scale3)))
# scale4 <- sample(hard_scheme, size = 35);names(scale4) <- as.character(seq(1, length(scale4)))
# scale5 <- sample(hard_scheme, size = 35);names(scale5) <- as.character(seq(1, length(scale5)))
# 
# scales::show_col(scale1)
# 
# scale2 <- sample(hard_scheme, size = 35);names(scale2) <- as.character(seq(1, length(scale2)))
# scales::show_col(scale2)
# 
# scales::show_col(scale3)
# scales::show_col(scale4)
# scales::show_col(scale5)

# plot_grid(
#   plotDR(daf, color_by = "meta25")+scale_color_manual(values=scale1),
#   plotDR(daf, color_by = "meta25")+scale_color_manual(values=scale2),
#   plotDR(daf, color_by = "meta25")+scale_color_manual(values=scale3),
#   plotDR(daf, color_by = "meta25")+scale_color_manual(values=scale4),
#   plotDR(daf, color_by = "meta25")+scale_color_manual(values=scale5), labels = c("scale1", "scale2", "scale3", "scale4", "scale5"))


# #### UMAP of filtered populations ####
# 
# plot_grid(
#   plotDR(daf, color_by = "meta25")+facet_wrap("timepoint")+scale_color_manual(values=meta25_color_scheme),
#   plotDR(meta_up_t6, color_by = "meta25")+facet_wrap("timepoint")+scale_color_manual(values=meta25_color_scheme)
#   )
# 
# plot_grid(
# plotDR(daf, color_by = "CD161")+facet_wrap("timepoint"),
# plotDR(meta_up_t6, color_by = "CD161")+facet_wrap("timepoint")
# )

# heatmaps of upregulated clusters ####

# extract lower level functions from catalyst in order to make heatmaps of median marker

# change yhr cd$cluster_id line to access different metaclustering/som codes

# expression per cluster ####
shplit_cells <- function(x, by) {
  stopifnot(is.character(by), by %in% colnames(colData(x)))
  cd <- data.frame(colData(x))
  cd$cluster_id <- cluster_ids(x, "meta35")
  dt <- data.table::data.table(cd, i = seq_len(ncol(x)))
  dt_split <- split(dt, by = by, sorted = TRUE, flatten = FALSE)
  map_depth(dt_split, length(by), "i")
}


ahgg <- function(x, by, fun = c("median", "mean", "sum")) {
  fun <- switch(match.arg(fun), 
                median = rowMedians, mean = rowMeans, sum = rowSums)
  cs <- shplit_cells(x, by)
  pb <- map_depth(cs, -1, function(i) {
    if (length(i) == 0) return(numeric(nrow(x)))
    fun(assay(x, "exprs")[, i, drop = FALSE])
  })
  map_depth(pb, -2, function(u) as.matrix(data.frame(
    u, row.names = rownames(x), check.names = FALSE)))
}


# make maxtrix of median expression ####

ms <- ahgg(meta_up_t6[type_markers(meta_up_t6), ], by = c("cluster_id", "sample_id"))
ms <- lapply(ms, reshape::melt, varnames = c("antigen", "sample_id"))
ms <- bind_rows(ms, .id = "cluster_id")


# transform each channel according to its min/max ####
scaled_ms <- ms %>%
  group_by(antigen) %>%
  group_by(sample_id) %>%
  dplyr::mutate(value=scales::rescale(value, to=c(0,1))) %>%
  ungroup()

# select channels to plot, select representative heatmap ####
my_palette <- c("#D53E4F","#D96459","#F2AE72","#588C73","#1A9CC7")

heatmap_channels <- c("CD4", "CD8", "Vd2", "Va72", "CD38", "CD69", "HLADR", "ICOS", "CD28", "PD1", 
                      "CD95", "BCL2", "CD27", "Perforin", "GZB", "Tbet", "Eomes", "CTLA4", "Ki67", "CD127", "CD56", "CD16",
                      "CD161", "CD49d", "CD25", "FoxP3", "CD39", "CXCR5", "CX3CR1", "CD57", 
                      "CD45RA", "CD45RO", "CCR7")

scaled_ms <- dplyr::filter(scaled_ms, antigen %in% heatmap_channels)
# scaled_ms$timepoint <- substr(scaled_ms$sample_id, 1, 1)

# scaled_ms <- scaled_ms %>%
#   group_by(timepoint, cluster_id, antigen) %>%
#   dplyr::mutate(median_expression=median(value))

scaled_ms <- dplyr::filter(scaled_ms, sample_id == "T6_v315")

specific_cluster_levels <- c("10", "17", "21", "22", "18", "27", "14", "8", "29","28", "32","23","31",  "26")

# plot ####
(up_t6_heatmap <- ggplot(data = scaled_ms, aes_(x = factor(scaled_ms$cluster_id, levels = specific_cluster_levels), y = factor(scaled_ms$antigen, levels = rev(t_cell_channels)), group=scaled_ms$sample_id))+
  geom_tile(aes(fill=scaled_ms$value))+
  scale_fill_gradientn(colors=rev(my_palette), "white")+
  xlab("Metacluster")+
  #ggtitle("Heatmap of Clusters Upregulated Post-Treatment in Naive Volunteers")
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line.y.left = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position="none", 
        legend.title = element_blank(),
        #plot.title = element_text(size = 45, hjust = 0.5),
        plot.margin = unit(c(1,0,1,0), "cm"))
)

ggsave("/Users/s1249052/PhD/presentations/NexGenImmunology_2020/primary_up_t6_heatmap.png", up_t6_heatmap, width=8.4, height=7.4)

# extract cluster frequencies #### <3
# dod_tbl <- data.frame(topTable(da_dod, show_props = TRUE, format_vals = TRUE, digits = 2))
# 
# sig_dod <- dplyr::select(dod_tbl, colnames(dod_tbl)[c(1:6,grep("DoD", colnames(dod_tbl), fixed=T))])
# sig_dod_clusters <- subset(rowData(da_dod$res), rowData(da_dod$res)$p_adj < FDR_cutoff)$cluster_id 
# 
# sig_dod <- filter(sig_dod, sig_dod$cluster_id %in% sig_dod_clusters)
# long_sig_dod <- gather(sig_dod, sample, frequency, colnames(sig_dod)[4:length(colnames(sig_dod))])
# long_sig_dod$volunteer <- substr(long_sig_dod$sample, nchar(long_sig_dod$sample)-3, nchar(long_sig_dod$sample))
# long_sig_dod$timepoint <- ifelse(grepl("DoD", long_sig_dod$sample)==T, "DoD", "Baseline")



t6_tbl <- data.frame(topTable(da_t6, show_props = TRUE, format_vals = TRUE, digits = 2))

sig_t6 <- dplyr::select(t6_tbl, colnames(t6_tbl)[c(1:6,grep("T6", colnames(t6_tbl), fixed=T))])
sig_t6_clusters <- subset(rowData(da_t6$res), rowData(da_t6$res)$p_adj < FDR_cutoff)$cluster_id 
sig_t6_clusters_up <- up_t6

sig_t6 <- dplyr::filter(sig_t6, sig_t6$cluster_id %in% sig_t6_clusters_up)
long_sig_t6 <- gather(sig_t6, sample, frequency, colnames(sig_t6)[4:length(colnames(sig_t6))])
long_sig_t6$volunteer <- substr(long_sig_t6$sample, nchar(long_sig_t6$sample)-3, nchar(long_sig_t6$sample))
long_sig_t6$timepoint <- ifelse(grepl("T6", long_sig_t6$sample)==T, "T6", "Baseline")

# figures of cluster frequencies ####


# ggplot(long_sig_dod, aes(x=timepoint, y=frequency, group=sample, fill=volunteer))+
#   geom_bar(stat = "identity", position = "dodge")+
#   theme_minimal()


subset_colors <- c("#41B3A3", "#C38D9E", "#E27D60", "#F64C72", "#E8A87C")
names(subset_colors) <- c("CD4+", "MAIT+", "CD8+", "VD2+", "DN")

#  !!!this file has to be sorted in excel to have cluster id in ascending order!!!
cluster_annotations <- read_excel("/Users/s1249052/PhD/presentations/NexGenImmunology_2020/flo_primary_cluster_merging1.xlsx")

cluster_subset <- stringr::str_to_upper(cluster_annotations$new_cluster)
#cluster_subset[grep("VD2", cluster_subset)] <- "vd2"
names(cluster_subset) <- seq(1,35)
long_sig_t6$cluster_kind <- cluster_subset[long_sig_t6$cluster_id]


cluster_names=paste("Metacluster ", specific_cluster_levels, sep='')
names(cluster_names) <- specific_cluster_levels

cluster_labeller <- function(variable,value){
  return(cluster_names[value])
}

long_sig_t6$cluster_id <- factor(long_sig_t6$cluster_id, levels=specific_cluster_levels)
long_sig_t6$cluster_kind <- factor(long_sig_t6$cluster_kind, levels=c("CD4+", "VD2+", "CD8+", "MAIT+", "DN"))
  
prim_boxplot <- ggplot(long_sig_t6, aes(x=timepoint, y=frequency))+
  geom_boxplot(aes(fill=cluster_kind))+
  geom_point(color="black", size=1.8, aes(shape=volunteer))+
  theme_minimal()+
  facet_wrap(~cluster_id, scales="free",labeller = cluster_labeller)+
  #facet_wrap(~cluster_id, scales="free")+
  scale_fill_manual(values=subset_colors, labels=c("CD4+", expression(paste(gamma, delta)), "CD8+", "MAIT",  "DN"), guide=guide_legend(label.hjust=0.5))+
  labs(fill = "", shape="Volunteer")+
  xlab("")+
  ylab("Percentage of CD3+ T cells")+
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18),
        strip.text = element_text(size=14))
  
ggsave("/Users/s1249052/PhD/presentations/NexGenImmunology_2020/prim_boxplot.png", prim_boxplot, width = 9, height=7)  
  



# hot pie stuff should go here ####

prime_pie_data <- select(long_sig_t6, cluster_id, volunteer, timepoint, frequency)
prime_pie_data <- prime_pie_data[order(prime_pie_data$cluster_id),]

sum_pi <- prime_pie_data %>%
     group_by(timepoint, cluster_id) %>%
     dplyr::mutate(cluster_sum=sum(frequency)) %>%
     ungroup() %>%
     select(-frequency)%>%
     slice(seq(from=1, to=nrow(.), by=3))

sum_pi <- sum_pi %>%
     group_by(timepoint) %>%
     dplyr::mutate(ymax=cumsum(cluster_sum)) %>%
     dplyr::arrange(timepoint)
     


sum_pi$ymin <- c(0, sum_pi$ymax[1:nrow(sum_pi)-1])

sum_pi <- sum_pi %>% 
  group_by(cluster_id) %>%
  dplyr::mutate(fold_change=cluster_sum[2]/cluster_sum[1])

sum_pi$fold_change <- ifelse(sum_pi$fold_change>100, 100, sum_pi$fold_change)

pie_theme <- theme(aspect.ratio=1,
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(), 
      legend.position = "none")

base_sum_pi <- dplyr::filter(sum_pi, timepoint=="Baseline")
t6_sum_pi <- dplyr::filter(sum_pi, timepoint=="T6")

base_pie <- ggplot()+
  geom_rect(data=base_sum_pi, aes(xmin=1, xmax=2, ymin=base_sum_pi$ymin/3, ymax=base_sum_pi$ymax/3, color=base_sum_pi$cluster_id, fill=base_sum_pi$cluster_id))+
  coord_polar(theta="y")+
  scale_fill_manual(values=rev(meta25_color_scheme))+
  scale_color_manual(values=rev(meta25_color_scheme))+
  xlim(0,3)+
  theme_void()+
  pie_theme


(t6_pie <- ggplot()+
  geom_rect(data=t6_sum_pi, aes(xmin=log(t6_sum_pi$fold_change, base = 8), xmax=log(t6_sum_pi$fold_change, base = 8)+1, ymin=t6_sum_pi$ymin/3, ymax=t6_sum_pi$ymax/3, color=t6_sum_pi$cluster_id, fill=t6_sum_pi$cluster_id))+
  coord_polar(theta="y")+
  scale_fill_manual(values=rev(meta25_color_scheme))+
  scale_color_manual(values=rev(meta25_color_scheme))+
  #xlim(0,3)+
  theme_void()+
  pie_theme)

plot_grid(base_pie, t6_pie, rel_widths = c(1,8))



ggplot()+
  geom_rect(data=sum_pi, aes(xmin=1, xmax=2, ymin=sum_pi$ymin/3, ymax=sum_pi$ymax/3, fill=sum_pi$cluster_id, group=sum_pi$cluster_id))+
  coord_polar(theta="y")+
  scale_fill_manual(values=rev(meta25_color_scheme))+
  xlim(0,3)+
  theme_minimal()


ggplot(pi, aes(x=timepoint, y=frequency/max, fill=cluster_id, group=volunteer))+
  geom_bar(stat="identity")+
  coord_polar(theta="y")+
  facet_wrap(~timepoint)
  theme_minimal()


ggplot(pi, aes(x=timepoint, y=frequency/max, fill=cluster_id, group=volunteer))+
  geom_bar(stat="identity")+
  coord_polar(theta="y")+
  theme_minimal()


sigt6 <- subset(rowData(da_t6$res), rowData(da_t6$res)$p_adj < FDR_cutoff) 

sigt6$cluster_id

sig_clusters <- daf[cluster_id %in% as.character(sigt6$cluster_id), "meta25")

cluster_ids <- cluster_ids(daf, "meta25")
sig_clusters <- as.character(sigt6$cluster_id)
sig_t6_daf <- daf[cluster_ids %in% sig_clusters,]


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
