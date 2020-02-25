library(Rphenograph)
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
library(Rtsne)


BiocManager::install("cytofkit")

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
scales::show_col(inferno_lite)

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
set.seed(1234);daf <- cluster(daf, features = refined_markers, xdim = 10, ydim = 10, maxK = 45, seed = 1234)

#daf <- filterSCE(daf, volunteer=="V03")

big_table <- t(data.frame(assays(daf)$exprs))

slim_table <- data.frame(big_table[,colnames(big_table) %in% refined_markers])
dim(slim_table)


# check that this is arcsinh transformed (cofactor=5) before running Phenograph!! this is default when
# processing with CATALYST beforehand

start_time <- Sys.time()
pheno <- Rphenograph(slim_table) # 33 clusters
end_time <- Sys.time()

end_time-start_time # 2h7min


UMAP_theme <- theme_minimal()+theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position = "none",
  axis.title = element_blank()
)

slim_umap <- data.frame(reducedDim(daf, "UMAP"))
colnames(slim_umap) <- c("UMAP1", "UMAP2")

slim_table <- cbind(slim_table, slim_umap)

slim_table$pheno_cluster <- factor(membership(pheno[[2]]))

slim_table$meta40 <- cluster_ids(daf, k="meta40")
slim_table$meta32 <- cluster_ids(daf, k="meta32")

# create dataset with only ~40,000 cells to make plotting faster
slimmed <- slim_table[seq(1, nrow(slim_table), by = 20),]



# plot making area ####

# color schemes

color_35_scheme <- c(
  "#4A3B53","#0AA6D8","#FFB500","#6A3A4C","#456648","#B79762","#3B5DFF","#0CBD66","#4FC601",
  "#D157A0","#34362D","#885578","#8FB0FF","#006FA6","#BA0900","#7A4900","#FFF69F","#FFDBE5",
  "#EEC3FF","#B4A8BD","#3A2465","#7A87A1","#1CE6FF","#6F0062","#0086ED","#6B002C","#C0B9B2",
  "#04F757","#000000","#FF2F80","#FFAA92","#A30059","#FF34FF","#013349","#C2FFED")

color_103_scheme <- c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                      "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                      "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                      "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                      "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                      "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                      "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                      "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                      "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                      "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                      "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                      "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                      "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")


# plots to compare clustering #

pheno_graph <- ggplot(slimmed, aes(x=slimmed$UMAP1, y=slimmed$UMAP2, color=slimmed$pheno_cluster))+
  geom_point(shape=1)+UMAP_theme+ggtitle("pheno_cluster, m=32")+scale_color_manual(values = color_35_scheme)

meta32_graph <- ggplot(slimmed, aes(x=slimmed$UMAP1, y=slimmed$UMAP2, color=slimmed$meta32))+
  geom_point(shape=1)+UMAP_theme+ggtitle("flowSOM, m=32")+scale_color_manual(values = color_35_scheme)

meta40_graph <- ggplot(slimmed, aes(x=slimmed$UMAP1, y=slimmed$UMAP2, color=slimmed$meta40))+
  geom_point(shape=1)+UMAP_theme+ggtitle("flowSOM, m=40")+scale_color_manual(values = color_103_scheme)


three_plots <- plot_grid(pheno_graph, meta32_graph, meta40_graph, ncol=3)
ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/phenograph_vs_flowsom.png", three_plots, height = 2.85, width=11.295)


# plots with important lineage markers

cd4 <- ggplot(slimmed, aes(x=slimmed$UMAP1, y=slimmed$UMAP2, color=slimmed$CD4))+
  geom_point(shape=1)+UMAP_theme+ggtitle("CD4")+
  scale_color_gradientn(colors=inferno_lite)

cd8 <- ggplot(slimmed, aes(x=slimmed$UMAP1, y=slimmed$UMAP2, color=slimmed$CD8))+
  geom_point(shape=1)+UMAP_theme+ggtitle("CD8")+
  scale_color_gradientn(colors=inferno_lite)

vd2 <- ggplot(slimmed, aes(x=slimmed$UMAP1, y=slimmed$UMAP2, color=slimmed$Vd2))+
  geom_point(shape=1)+UMAP_theme+ggtitle("Vd2")+
  scale_color_gradientn(colors=inferno_lite)

va7 <- ggplot(slimmed, aes(x=slimmed$UMAP1, y=slimmed$UMAP2, color=slimmed$Va72))+
  geom_point(shape=1)+UMAP_theme+ggtitle("Va7.2")+
  scale_color_gradientn(colors=inferno_lite)


lin_plots <- plot_grid(cd4, cd8, vd2, va7, ncol=2)
ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/lineage_umap.png", lin_plots, height = 5.71, width=7.53)


# plots with important memory/differentiation markers

cd27 <- ggplot(slimmed, aes(x=slimmed$UMAP1, y=slimmed$UMAP2, color=slimmed$CD27))+
  geom_point(shape=1)+UMAP_theme+ggtitle("CD27")+
  scale_color_gradientn(colors=inferno_lite)

cd57 <- ggplot(slimmed, aes(x=slimmed$UMAP1, y=slimmed$UMAP2, color=slimmed$CD28))+
  geom_point(shape=1)+UMAP_theme+ggtitle("CD57")+
  scale_color_gradientn(colors=inferno_lite)

cd45ro <- ggplot(slimmed, aes(x=slimmed$UMAP1, y=slimmed$UMAP2, color=slimmed$CD45RO))+
  geom_point(shape=1)+UMAP_theme+ggtitle("CD45RO")+
  scale_color_gradientn(colors=inferno_lite)

ccr7 <- ggplot(slimmed, aes(x=slimmed$UMAP1, y=slimmed$UMAP2, color=slimmed$CCR7))+
  geom_point(shape=1)+UMAP_theme+ggtitle("CCR7")+
  scale_color_gradientn(colors=inferno_lite)

cd25 <- ggplot(slimmed, aes(x=slimmed$UMAP1, y=slimmed$UMAP2, color=slimmed$CD25))+
  geom_point(shape=1)+UMAP_theme+ggtitle("CD25")+
  scale_color_gradientn(colors=inferno_lite)

cd127 <- ggplot(slimmed, aes(x=slimmed$UMAP1, y=slimmed$UMAP2, color=slimmed$CD127))+
  geom_point(shape=1)+UMAP_theme+ggtitle("CD127")+
  scale_color_gradientn(colors=inferno_lite)


type_plots <- plot_grid(cd45ro, ccr7, cd25, cd127, cd27, cd57, ncol = 3)
ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/type_umap.png", type_plots, height = 5.71, width=11.295)


#plots with important activation markers

cd38 <- ggplot(slimmed, aes(x=slimmed$UMAP1, y=slimmed$UMAP2, color=slimmed$CD38))+
  geom_point(shape=1)+UMAP_theme+ggtitle("CD38")+
  scale_color_gradientn(colors=inferno_lite)

bcl2 <- ggplot(slimmed, aes(x=slimmed$UMAP1, y=slimmed$UMAP2, color=slimmed$BCL2))+
  geom_point(shape=1)+UMAP_theme+ggtitle("Bcl-2")+
  scale_color_gradientn(colors=inferno_lite)

hladr <- ggplot(slimmed, aes(x=slimmed$UMAP1, y=slimmed$UMAP2, color=slimmed$HLADR))+
  geom_point(shape=1)+UMAP_theme+ggtitle("HLA-DR")+
  scale_color_gradientn(colors=inferno_lite)


state_plots <- plot_grid(cd38, bcl2, hladr, ncol = 3)
ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/state_umap.png", state_plots, height = 2.85, width=11.295)

# quality control plots of batches, timepoints and volunteers;

time_point_plot <- plotDR(daf, color_by = "volunteer")+facet_wrap("timepoint")+UMAP_theme+scale_size_manual(values=c(2,1.5,1,0.5))
time_point_plot$layers[[1]]$aes_params$shape <- 1
#time_point_plot$layers[[1]]$aes_params$alpha <- 0.3

batch_plot <- plotDR(daf, color_by = "meta35")+facet_wrap("batch")+UMAP_theme+scale_size_manual(values=c(1, 0.5))+
  scale_color_manual(values=color_35_scheme)

batch_plot$layers[[1]]$aes_params$shape <- 1
#batch_plot$layers[[1]]$aes_params$layers[[1]]$aes_params$alpha <- 0.3

ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/time_point_plot.png", time_point_plot, height = 5.71, width=7.53)
ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/batch_plot.png", batch_plot, height = 2.855, width=7.53)
