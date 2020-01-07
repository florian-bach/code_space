library(readxl) 
library(CATALYST)
library(flowCore)
library(ggplot2)
library(RColorBrewer)
library(SingleCellExperiment)
library(cowplot)

`%!in%` = Negate(`%in%`)
start_time <- Sys.time()
###


# in order to create a daFrame, the single cell experiment, files, metadata and panel information need to be assembled and passed to
# CATALYST


### read in metadata
setwd("/Users/s1249052/PhD/cytof/vac69a/T_cells_only/fcs")

md <- read.csv("meta_data.csv", header=T) 


## md has to have particular properties: file_name=NAMED LIST (chr), ID and everything else in factors
## reorder so that it's neat 
md$timepoint <- factor(md$timepoint, levels = c("Baseline", "C8", "C10", "C12", "DoD", "T6"))
md <- md[
        with(md, order(md$volunteer, md$timepoint)),
         ]
md$file_name <- as.character(md$file_name)


#change number of volunteers here
md <- dplyr::filter(md, volunteer %in% c("V315", "V313", "V320"))


md$timepoint <- factor(md$timepoint, levels = c("Baseline", "DoD", "T+6", "C+45"))
md$sample_id <- factor(md$sample_id, levels = md$sample_id[order(md$timepoint)])

# select whose files to import

fcs_files <- grep("fcs", list.files(), value = T)

primaries <- c()

#change number of volunteers here
timepoints <- c("C-1", "C_12", "DoD", "T_6")

for(i in timepoints){
  one_moment <- grep(i, fcs_files, value=T)
  primaries <- c(timepoints, one_moment)
}

vac63c_primaries <- read.flowSet(primaries)











### read in flowfiles using flowCore
vac69a <- read.flowSet(md$file_name)


### define panel using first fcs file in flowset, then clean up

panel <- data.frame("fcs_colname"=colnames(vac69a[[1]]), "antigen"=c("Time", markernames(vac69a[[1]])))
panel$antigen <- gsub("_", "", panel$antigen)

panel$antigen <- ifelse(nchar(panel$antigen)>6,
                        substr(panel$antigen, 6, nchar(panel$antigen)),
                        panel$antigen
)

panel$antigen[3] <- "CD45"

panel$antigen <- gsub("-", "", panel$antigen)
panel$antigen <- gsub(".", "", panel$antigen, fixed=T)




### construct daFrame
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



daf100 <- cluster(daf, cols_to_use = t_cell_channels, xdim = 10, ydim = 10, maxK = 15, seed = 1234)
daf196 <- cluster(daf, cols_to_use = t_cell_channels, xdim = 14, ydim = 14, maxK = 25, seed = 1234)

daf100_delta <- metadata(daf100)$delta_area
daf196_delta <- metadata(daf196)$delta_area
daf100_delta
daf196_delta

end_time <- Sys.time()

duration = end_time - start_time
#duration = ~5min


plotClusterHeatmap(daf100, hm2 = NULL, k = "meta15", m = NULL, cluster_anno = TRUE, draw_freqs = TRUE) 
plotClusterHeatmap(daf196, hm2 = NULL, k = "meta25", m = NULL, cluster_anno = TRUE, draw_freqs = TRUE) 

start_time <- Sys.time()

daf100 <- runDR(daf100, "TSNE", rows_to_use = 500)
daf100 <- runDR(daf100, "UMAP", rows_to_use = 5000)

daf196 <- runDR(daf196, "TSNE", rows_to_use = 500)
daf196 <- runDR(daf196, "UMAP", rows_to_use = 5000)

end_time <- Sys.time() #8 min

p1 <- plotDR(daf100, "UMAP", color_by="meta15", facet=c("timepoint"))+ theme(legend.position = "none")

p2 <- plotDR(daf196, "UMAP", color_by="meta25", facet=c("timepoint"))+ theme(legend.position = "none")

lgd <- get_legend(p2)

plot_grid(p1, p2, lgd, nrow = 1, rel_widths = c(5, 5, 2))



myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
red_palette <- c(myPalette(100), rep(myPalette(100)[100], 200))

sc <- scale_colour_gradientn(colours = red_palette, limits=c(0, 8))


cd38_plot <- plotDR(daf100, "UMAP", color_by="CD38", facet="timepoint")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.title = element_blank())+
  sc

lgd <- get_legend(p2)


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
