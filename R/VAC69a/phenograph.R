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

daf <- filterSCE(daf, volunteer=="V03")

big_table <- t(data.frame(assays(daf)$exprs))

slim_table <- data.frame(big_table[,colnames(big_table) %in% refined_markers])
dim(slim_table)


# check that this is arcsinh transformed (cofactor=5) before running Phenograph!! #

start_time <- Sys.time()
pheno <- Rphenograph(slim_table)


slim_table$pheno_cluster <- factor(membership(pheno[[2]]))

slim_tbl <- dplyr::select(slim_table, -pheno_cluster)
slim_sne <- Rtsne(slim_tbl, perplexity=30)
end_time <- Sys.time()

slim_table <- cbind(slim_table, slim_sne[2])

tsne_theme <- theme_void()+theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position = "none", 
  axis.title = element_blank()
)



pheno_graph <- ggplot(slim_table, aes(x=slim_table$Y.1, y=slim_table$Y.2, color=slim_table$pheno_cluster))+
  geom_point()+tsne_theme+ggtitle("pheno_cluster")+
  scale_color_gradientn(colors=inferno_lite)

cd4 <- ggplot(slim_table, aes(x=slim_table$Y.1, y=slim_table$Y.2, color=slim_table$CD4))+
  geom_point()+tsne_theme+ggtitle("CD4")+
  scale_color_gradientn(colors=inferno_lite)

cd8 <- ggplot(slim_table, aes(x=slim_table$Y.1, y=slim_table$Y.2, color=slim_table$CD8))+
  geom_point()+tsne_theme+ggtitle("CD8")+
  scale_color_gradientn(colors=inferno_lite)

vd2 <- ggplot(slim_table, aes(x=slim_table$Y.1, y=slim_table$Y.2, color=slim_table$Vd2))+
  geom_point()+tsne_theme+ggtitle("Vd2")+
  scale_color_gradientn(colors=inferno_lite)

va7 <- ggplot(slim_table, aes(x=slim_table$Y.1, y=slim_table$Y.2, color=slim_table$Va72))+
  geom_point()+tsne_theme+ggtitle("Va7.2")+
  scale_color_gradientn(colors=inferno_lite)

cd38 <- ggplot(slim_table, aes(x=slim_table$Y.1, y=slim_table$Y.2, color=slim_table$CD38))+
  geom_point()+tsne_theme+ggtitle("CD38")+
  scale_color_gradientn(colors=inferno_lite)

bcl2 <- ggplot(slim_table, aes(x=slim_table$Y.1, y=slim_table$Y.2, color=slim_table$BCL2))+
  geom_point()+tsne_theme+ggtitle("Bcl-2")+
  scale_color_gradientn(colors=inferno_lite)

pheno_overview <- plot_grid(cd4, cd8, vd2, va7, cd38, bcl2, ncol = 2)

ggsave(pheno_graph, "/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/pheno_graph.png")
ggsave(pheno_overview, "/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/pheno_overview.png")


