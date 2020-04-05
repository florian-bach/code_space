library(ggcyto)
library(flowCore)
library(CATALYST)
library(ggplot2)
library(viridis)
library(SingleCellExperiment)
library(cowplot)
library(SummarizedExperiment)

#functions, palettes etc. ####

`%!in%` = Negate(`%in%`)

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

inferno_mega_lite <- c("#000004", "#8A2267", "#EF802B", "#FFEC89", "#FCFFA4")

UMAP_theme <- theme_minimal()+theme(
  panel.grid.major = element_blank(),
  legend.position = "none",
  axis.title = element_blank(),
  plot.title = element_text(hjust = 0.5)
)


flo_umap <- function(df, color_by, facet_by=NULL){
  
  if(is.null(facet_by))
    assign("facet", ~"")
    assign("facet_title", element_blank())
  
  data <- dplyr::select(df, UMAP1, UMAP2, color_by, facet_by)
  colnames(data)[3] <- "color"
  

  data$color <- scales::rescale(data$color, to=c(0,5))
  

  (plt <- ggplot(data, aes(x=UMAP1, y=UMAP2, color=color))+
      geom_point(size=0.7)+
      scale_color_gradientn(colors = inferno_mega_lite)+
      #facet_wrap(facet_two~facet_one)
      UMAP_theme+
      facet_wrap(facet)+
      ggtitle(color_by)+
      theme(strip.text = facet_title)
      
    )
  
  
}

# read in data ####

setwd("C:/Users/bachf/PhD/cytof/vac69a/T_cells_only/fcs")
#setwd("/home/flobuntu/PhD/cytof/vac69a/T_cells_only/fcs")
#setwd("/Users/s1249052/PhD/cytof/vac69a/T_cells_only/fcs")

RNGversion("3.6.3")

md <- read.csv("new_files.csv", header=T) 

md$timepoint <- factor(md$timepoint, levels = c("Baseline", "C8", "C10", "C12", "DoD", "T6"))
md <- md[
  with(md, order(md$volunteer, md$timepoint)),
  ]
md$file_name <- as.character(md$file_name)


#change number of volunteers/timepoints here
md <- dplyr::filter(md, timepoint %in% c("Baseline", "C10" , "DoD", "T6"))

#import panel 
panel <- read.csv("VAC69_PANEL.CSV", header=T, stringsAsFactors = F)
colnames(panel)[2] <- "marker_name"

# select whose files to import
fcs_file
s <- grep("fcs", list.files(), value = T)

primaries <- c()

timepoints <- c("C-1", "C_10", "DoD", "T_6")

for(i in timepoints){
  one_moment <- grep(i, fcs_files, value=T)
  primaries <- c(timepoints, one_moment)
}

### read in flowfiles using flowCore
vac69a <- read.flowSet(md$file_name)

#downsampling to equal size
sampling_ceiling <- 3000
# Being reproducible is a plus
set.seed(1234)

# sample.int takes a sample of the specified size from the elements of x using either with or without replacement.
# smaller_vac69a <- fsApply(vac69a, function(ff) {
#   idx <- sample.int(nrow(ff), min(sampling_ceiling, nrow(ff)))
#   ff[idx,]  # alt. ff[order(idx),]
# })

downsample <- function(fs, floor){fsApply(fs, function(ff) {
  idx <- sample.int(nrow(ff), nrow(ff)/min(fsApply(fs, nrow))*floor)
  ff[idx,]  # alt. ff[order(idx),]
})}

smaller_vac69a <- downsample(vac69a, 2000)

smol_daf <- prepData(smaller_vac69a, panel, md, md_cols =
                  list(file = "file_name", id = "sample_id", factors = c("timepoint", "batch", "volunteer")),
                panel_cols = list(channel = "fcs_colname", antigen = "marker_name", class =
                                    "marker_class"))

smol_daf <- cluster(smol_daf, features = refined_markers, xdim = 10, ydim = 10, maxK = 45, seed = 1234)


# plotClusterHeatmap(smol_daf, hm2 = NULL,
#                    k = "meta45",
#                    # m = "flo_merge",
#                    cluster_anno = TRUE,
#                    draw_freqs = TRUE,
#                    scale=T
# )
# 


smol_daf <- filterSCE(smol_daf, timepoint!="C10")
set.seed(1234)

# umap projections ####
smol_daf <- scater::runUMAP(smol_daf,
                             subset_row=refined_markers,
                             exprs_values = "exprs",
                             scale=T)


# 
# not_scaled7 <- plotDR(not_scaled_daf, color_by = "CCR7")+facet_wrap("timepoint")
# not_scaled45 <- plotDR(not_scaled_daf, color_by = "CD45RA")+facet_wrap("timepoint")
# 
# cata <- runDR(smol_daf, "UMAP", features = refined_markers, assay = "exprs")
# plotDR(cata, color_by="CD4")
# 
# cowplot::plot_grid(not_scaled7, scaled7,not_scaled45, scaled45, ncol=2)
# 
# # pca 50 / NULL ?
# 
# start <- Sys.time(); smol_daf <- runDR(smol_daf, cells=2500, dr="UMAP", features=refined_markers, assay = "exprs"); end <- Sys.time()
# (duration <- round(end - start)) #2min
# 
# plotDR(smol_daf, "UMAP", color_by = "CD4")



big_table <- data.frame(t(data.frame(assays(smol_daf)$exprs)))
big_table <- data.frame(cbind(big_table, colData(smol_daf)))

slim_umap <- data.frame(reducedDim(smol_daf, "UMAP"))
colnames(slim_umap) <- c("UMAP1", "UMAP2")

big_table <- data.frame(cbind(big_table, slim_umap), stringsAsFactors = F)


(smol_time12 <- ggplot(big_table, aes(x=UMAP1, y=UMAP2))+
  stat_density_2d(aes(fill = after_stat(nlevel)), geom="polygon", bins=20)+
  scale_fill_gradientn(colors = inferno_mega_lite)+
  xlim(c(-14, 12))+
  ylim(c(-10, 11))+
  theme_minimal()+
  facet_wrap(~timepoint)+
  UMAP_theme+
  theme(strip.text = element_text(size=14)))

setwd("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper")

ggsave("proportional_contour_umap.png", smol_time12, height=6, width=9)


#ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/umap_through_time.png", smol_time12, width=15, height=5.54)
#ggsave("/home/flo/PhD/cytof/vac69a/figures_for_paper/umap_through_time.png", smol_time12, width=15, height=5.54)


big_time12 <- ggplot(big_table, aes(x=UMAP1, y=UMAP2))+
  stat_density_2d(aes(fill = after_stat(nlevel)), geom="polygon", bins=25)+
  scale_fill_gradientn(colors = inferno_mega_lite)+
  xlim(c(-11, 14))+
  ylim(c(-9, 10))+
  theme_minimal()+
  facet_grid(timepoint~volunteer, switch = "y")+
  UMAP_theme+
  theme(strip.text.y.left = element_text(size=14, angle = 0),
        strip.text.x = element_text(size=14),
        strip.placement = "outside")



ggsave("big_time12.png", big_time12, height=12, width=18)



plotDR(smol_daf, color_by="CD45RO")+facet_grid(timepoint~volunteer)



cd4_plot <- flo_umap(big_table, "CD4")
cd8_plot <- flo_umap(big_table, "CD8")
vd2_plot <- flo_umap(big_table, "Vd2")
va72_plot <- flo_umap(big_table, "Va72")


scaled_lineage_plot <- plot_grid(cd4_plot, cd8_plot, vd2_plot, va72_plot, ncol=2)
ggsave("scaled_lineage_plot.png", scaled_lineage_plot)
# ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/lineage_plot_ring.png", lineage_plot)
ggsave("/home/flo/PhD/cytof/vac69a/figures_for_paper/lineage_plot.png", lineage_plot)


cd38_plot <- flo_umap(big_table, "CD38")
hladr_plot <- flo_umap(big_table, "HLADR")
bcl2_plot <- flo_umap(big_table, "BCL2")
cd27_plot <- flo_umap(big_table, "CD27")

activation_plot  <- plot_grid(cd38_plot, bcl2_plot,  hladr_plot, cd27_plot, ncol=2)
# ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/activation_plot.png", activation_plot)
ggsave("/home/flo/PhD/cytof/vac69a/figures_for_paper/activation_plot.png", activation_plot)


cd45ro_plot <- flo_umap(big_table, "CD45RO")
ccr7_plot <- flo_umap(big_table, "CCR7")
cd57_plot <- flo_umap(big_table, "CD57")
cd28_plot <- flo_umap(big_table, "CD28")

memory_plot <- plot_grid(cd45ro_plot, ccr7_plot, cd57_plot, cd28_plot, ncol=2)
# ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/memory_plot.png", memory_plot)
ggsave("/home/flo/PhD/cytof/vac69a/figures_for_paper/memory_plot.png", memory_plot)



  
cd38_time_plot <- flo_umap(big_table, "CD38")
hladr_time_plot <- flo_umap(big_table, "HLADR")
bcl2_time_plot <- flo_umap(big_table, "BCL2")
cd27_time_plot <- flo_umap(big_table, "CD27")

activation_time_plot <- plot_grid(cd38_time_plot, hladr_time_plot, bcl2_time_plot, cd27_time_plot, ncol=2)
# ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/activation_time_plot.png", activation_time_plot, width=8, height=5.54)
ggsave("/home/flo/PhD/cytof/vac69a/figures_for_paper/activation_time_plot.png", activation_time_plot, width=8, height=5.54)



cd25_plot <- flo_umap(big_table, "CD25")
cd127_plot <- flo_umap(big_table, "CD127")
foxp3_plot <- flo_umap(big_table, "FoxP3")
cd39_plot <- flo_umap(big_table, "CD39")

treg_plot <- plot_grid(cd25_plot, cd127_plot, foxp3_plot, cd39_plot, ncol=2)
# ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/treg_plot.png", treg_plot, width=8, height=5.54)
ggsave("/home/flo/PhD/cytof/vac69a/figures_for_paper/treg_plot.png", treg_plot, width=8, height=5.54)


# plot with points and contour lines

# smol_time+
#   stat_density_2d(aes(fill = stat(nlevel)), geom="polygon", bins=11)+
#   stat_density_2d(size=0.13, color="black", bins=11)+
#   #stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon", bins=22, size=0.3, color="black", contour=T)+
#   # scale_fill_viridis(option="A")+
#   scale_color_viridis(option="A", direction=1)+
#   #scale_fill_viridis(option="B", direction=1)+
#   scale_fill_gradient(low = "white", high="black")+
#   geom_point(shape=".")+
#   xlim(c(-14, 11))+
#   ylim(c(-10, 10))+
#   theme_minimal()+
#   UMAP_theme




# regular scatter plots ####
big_table <- data.frame(t(data.frame(assays(smol_daf)$exprs)))
big_table <- data.frame(cbind(big_table, colData(smol_daf)))

ggplot(big_table, aes(x=CD3, y=CD56))+
  geom_hex(bins=128)+
  # scale_x_flowCore_fasinh()+
  # scale_y_flowCore_fasinh()+
  scale_fill_viridis(option="A")

ggplot(big_table, aes(x=CD3, y=CD16))+
  geom_hex(bins=128)+
  # scale_x_flowCore_fasinh()+
  # scale_y_flowCore_fasinh()+
  scale_fill_viridis(option="A")

ggplot(big_table, aes(x=CD4, y=CD8))+
  geom_point(shape=".")+
  #geom_hex(bins=128)+
  facet_wrap(~timepoint)+
  scale_x_flowCore_fasinh()+
  scale_y_flowCore_fasinh()+
  scale_color_viridis(option="A")+
  scale_fill_viridis(option="A")+
  theme_minimal()


# p <- ggcyto(smaller_vac69a, aes(x=CD3, y=CD56))
#   
# #   geom_hex(bins=8000)+
# #   geom_point(shape='.')+
# #   scale_x_flowCore_fasinh()+
# #   scale_y_flowCore_fasinh()+
# #   theme_minimal()
# 
# p+geom_point(shape='.')+
#   scale_x_flowCore_fasinh()+
#   scale_y_flowCore_fasinh()

 