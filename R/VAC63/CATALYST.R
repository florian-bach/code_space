library(readxl) 
library(CATALYST)
library(flowCore)
library(ggplot2)
library(RColorBrewer)

setwd("/Users/s1249052/PhD/cytof/vac63c/fcs_files/T\ cells/")


md <- read_excel("meta_data.xlsx") 

#cchange number of volunteers here
md <- dplyr::filter(md, volunteer %in% c("V315", "V313", "V320"))


md$timepoint <- factor(md$timepoint, levels = c("Baseline", "DoD", "T+6", "C+45"))
md$sample_id <- factor(md$sample_id, levels = md$sample_id[order(md$timepoint)])

# select whose files to import

fcs_files <- grep("fcs", list.files("/Users/s1249052/PhD/cytof/vac63c/fcs_files/T\ cells/"), value = T)

primaries <- c()

#change number of volunteers here
people <- c("315", "313", "320")

for(i in people){
  one_person <- grep(i, fcs_files, value=T)
  primaries <- c(primaries, one_person)
}

vac63c_primaries <- read.flowSet(primaries)



## antigen colnames CANNOT start with a digit, otherwise it won't parse

try <- vac63c_primaries[[1]]
panel <- data.frame("fcs_colname"=colnames(try), "antigen"=c("Time", markernames(try)))
panel$antigen <- gsub("_", "", panel$antigen)

panel$antigen <- ifelse(nchar(panel$antigen)>6,
                                       substr(panel$antigen, 6, nchar(panel$antigen)),
                                       panel$antigen
                                       )
panel$antigen[3] <- "CD45"


# spot check that all panel columns are in the flowSet object
# all(panel$fcs_colname %in% colnames(vac63c_primaries))
# all(colnames(vac63c_primaries) %in% panel$fcs_colname)
# # 
# all(md$file_name %in% sampleNames(vac63c_primaries))
# all(sampleNames(vac63c_primaries) %in% md$file_name)




## THIS ABSOLUTELY HAS TO BE SPECIFIED!!! NAMED LIST MUST ONLY HAVE FILE, ID AND THEN EVERYTHING ELSE IN FACTORS
daf <- daFrame(vac63c_primaries, panel, md, md_cols =
                 list(file = "file_name", id = "sample_id", factors = c("timepoint", "batch", "volunteer")))


# this gets rid of barcoding and quality control channels so clustering is easier
proper_channels <- colnames(daf)[c(3,14:15,23:57,63,65)]
t_cell_channels <- proper_channels[c(-1, -2, -32)]

plotMDS(daf, color_by = "timepoint")

plotExprHeatmap(daf, bin_anno = TRUE, row_anno = TRUE)

set.seed(1234)


# run FlowSOM
daf <- cluster(daf, cols_to_use = t_cell_channels, xdim = 10, ydim = 10, maxK = 20, seed = 1234)
daf <- runDR(daf, "UMAP", rows_to_use = 5000, cols_to_use=t_cell_channels, overwrite=T)


# density plot, no contours
# vd2_plot <- plotDR(daf, "UMAP", color_by="Vd2", facet=c("volunteer","timepoint"))+
#   stat_density2d(aes(fill = ..density..), geom = 'tile', contour=F, bins=100, size=0.2)
# 
# vd2_plot$layers[[1]] <- NULL
# vd2_plot+xlim(c(-7, 15))+ylim(c(-7, 15))+
#   scale_fill_distiller(palette = 'RdPu', trans="reverse")


myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
red_palette <- c(myPalette(100), rep(myPalette(100)[100], 100))

sc <- scale_colour_gradientn(colours = red_palette, limits=c(0, 8))



# density plot, with contours
density_plot <- plotDR(daf, "UMAP", color_by="CCR7", facet=c("timepoint"))+
  stat_density2d(bins=100, size=0.18, colour="maroon")

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


bcl2_plot <- plotDR(daf, "UMAP", color_by="BCL2", facet="timepoint")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.title = element_blank())+
  sc

ggsave("dens_cd38_bcl2.png", cowplot::plot_grid(density_plot, cd38_plot, bcl2_plot, nrow=3), dpi=400, width=6.4, height = 6.4)



plotDR(daf, "UMAP", color_by="Vd2", facet="timepoint")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  sc




