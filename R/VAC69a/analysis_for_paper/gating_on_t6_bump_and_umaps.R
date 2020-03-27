library(ggcyto)
library(flowCore)
library(CATALYST)
library(SingleCellExperiment)
library(cowplot)


#functions, palettes etc. ####

`%!in%` = Negate(`%in%`)

cd4_t6_markers <- c(#"CD4",
                     #"CD8",
                     #"Vd2",
                     #"Va72",
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


inferno_lite <- c("#000004", "#380D57", "#460B68","#8A2267", "#9A2964",
                  "#AA325B","#B93B53","#C84449", "#D74D3F", "#E15D37",
                  "#E86F32","#EF802B", "#F59121", "#FBA210", "#FEB431", "#FFC751",
                  "#FFDA6D", "#FFEC89", "#FCFFA4")

inferno_lite <- colorRampPalette(inferno_lite)


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
  

# GATING ####

# read in data
setwd("C:/Users/bachf/PhD/cytof/vac69a/T_cells_only/fcs/post_umap")

flz <- list.files(pattern = "*.fcs")

vac69a <- read.flowSet(flz)


sampling_ceiling <- 3000
# Being reproducible is a plus
set.seed(1234)

# sample.int takes a sample of the specified size from the elements of x using either with or without replacement.
smaller_vac69a <- fsApply(vac69a, function(ff) {
  idx <- sample.int(nrow(ff), min(sampling_ceiling, nrow(ff)))
  ff[idx,]  # alt. ff[order(idx),]
})

#gate that mf

t6_gate <-rectangleGate("UMAP2"=c(1, 2.85), "UMAP1"=c(-2.4, -0.5))

result <- Subset(smaller_vac69a, flowCore::filter(smaller_vac69a, t6_gate))


ggcyto(vac69a[c(7,8,15,16, 23, 24)], aes(x=UMAP1, y=UMAP2))+
  stat_density_2d(aes(fill = after_stat(nlevel)), geom="polygon", bins=25)+
  xlim(c(-11, 14))+
  ylim(c(-9, 10))+
  #geom_hex(bins = 190)+
  #geom_point(shape=".")+
  #geom_density()+
  geom_gate(t6_gate)+
  geom_stats(type = "percent", adjust =1.85)+
  theme_minimal()

t6_bump <- Subset(vac69a, flowCore::filter(vac69a, t6_gate))

#import into CATALYST ####

setwd("C:/Users/bachf/PhD/cytof/vac69a/T_cells_only/fcs")

md <- read.csv("new_files.csv", header=T) 

md$timepoint <- factor(md$timepoint, levels = c("Baseline", "C10", "DoD", "T6"))

md$file_name <- as.character(md$file_name)

md$file_name <- paste(substr(md$file_name, 0, nchar(md$file_name)-4),"_post_umap.fcs", sep='')

#import panel 
panel <- read.csv("VAC69_PANEL.CSV", header=T, stringsAsFactors = F)
colnames(panel)[2] <- "marker_name"

#the order of entries in the md file matter, be careful
t6_daf <- prepData(t6_bump, panel, md, md_cols =
                       list(file = "file_name", id = "sample_id", factors = c("timepoint", "batch", "volunteer")),
                       panel_cols = list(channel = "fcs_colname", antigen = "marker_name", class =
                                         "marker_class"))

t6_daf <- cluster(t6_daf, features = cd4_t6_markers, xdim = 10, ydim = 10, maxK = 45, seed = 1234)


# plotClusterHeatmap(smol_daf, hm2 = NULL,
#                    k = "meta45",
#                    # m = "flo_merge",
#                    cluster_anno = TRUE,
#                    draw_freqs = TRUE,
#                    scale=T
# )
# 


t6_daf <- filterSCE(t6_daf, timepoint!="C10")
set.seed(1234)

# umap projections ####
t6_daf <- scater::runUMAP(t6_daf,
                            subset_row=cd4_t6_markers,
                            exprs_values = "exprs",
                            scale=T)

plotDR(t6_daf, color_by = "meta10")+facet_grid(volunteer~timepoint)
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



big_table <- data.frame(t(data.frame(assays(t6_daf)$exprs)))
big_table <- data.frame(cbind(big_table, colData(t6_daf)))

slim_umap <- data.frame(reducedDim(t6_daf, "UMAP"))
colnames(slim_umap) <- c("UMAP1", "UMAP2")

big_table <- data.frame(cbind(big_table, slim_umap), stringsAsFactors = F)


smol_time12 <- ggplot(big_table, aes(x=UMAP1, y=UMAP2))+
  stat_density_2d(aes(fill = after_stat(nlevel)), geom="polygon", bins=25)+
  scale_fill_gradientn(colors = inferno_mega_lite)+
  xlim(c(-11, 14))+
  ylim(c(-9, 10))+
  theme_minimal()+
  facet_wrap(~timepoint)+
  UMAP_theme+
  theme(strip.text = element_text(size=14))




cd38_plot <- flo_umap(big_table, "CD38")
hladr_plot <- flo_umap(big_table, "HLADR")
bcl2_plot <- flo_umap(big_table, "BCL2")
cd27_plot <- flo_umap(big_table, "CD27")

activation_plot  <- plot_grid(cd38_plot, bcl2_plot,  hladr_plot, cd27_plot, ncol=2)

ggsave("C:/Users/bachf/PhD/cytof/vac69a/figures_for_paper/umap/t6_activation_plot.png", activation_plot)
# ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/activation_plot.png", activation_plot)
# ggsave("/home/flo/PhD/cytof/vac69a/figures_for_paper/activation_plot.png", activation_plot)


cd45ro_plot <- flo_umap(big_table, "CD45RO")
ccr7_plot <- flo_umap(big_table, "CCR7")
cd57_plot <- flo_umap(big_table, "CD57")
cd45ra_plot <- flo_umap(big_table, "CD45RA")

memory_plot <- plot_grid(cd45ro_plot, ccr7_plot, cd45ra_plot, cd57_plot, ncol=2)
ggsave("C:/Users/bachf/PhD/cytof/vac69a/figures_for_paper/umap/t6_memory_plot.png", memory_plot)
# ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/memory_plot.png", memory_plot)
# ggsave("/home/flo/PhD/cytof/vac69a/figures_for_paper/memory_plot.png", memory_plot)



cd25_plot <- flo_umap(big_table, "CD25")
cd127_plot <- flo_umap(big_table, "CD127")
foxp3_plot <- flo_umap(big_table, "FoxP3")
cd39_plot <- flo_umap(big_table, "CD39")

treg_plot <- plot_grid(cd25_plot, cd127_plot, foxp3_plot, cd39_plot, ncol=2)
ggsave("C:/Users/bachf/PhD/cytof/vac69a/figures_for_paper/umap/t6_treg_plot.png", treg_plot)
# ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/treg_plot.png", treg_plot, width=8, height=5.54)
# ggsave("/home/flo/PhD/cytof/vac69a/figures_for_paper/treg_plot.png", treg_plot, width=8, height=5.54)




pd1_plot <- flo_umap(big_table, "PD1")
ctla4_plot <- flo_umap(big_table, "CTLA4")
icos_plot <- flo_umap(big_table, "ICOS")
cd161_plot <- flo_umap(big_table, "CD161")

t6_bump_plot <- plot_grid(pd1_plot, ctla4_plot, icos_plot, cd161_plot, ncol=2)
ggsave("C:/Users/bachf/PhD/cytof/vac69a/figures_for_paper/umap/t6_bump_plot.png", t6_bump_plot)


ki67_plot <- flo_umap(big_table, "Ki67")
tbet_plot <- flo_umap(big_table, "Tbet")
foxp3_plot <- flo_umap(big_table, "FoxP3")
gzb_plot <- flo_umap(big_table, "GZB")

t6_bump_2_plot <- plot_grid(ki67_plot, tbet_plot, foxp3_plot, gzb_plot, ncol=2)
ggsave("C:/Users/bachf/PhD/cytof/vac69a/figures_for_paper/umap/t6_bump_2_plot.png", t6_bump_2_plot)









plotClusterHeatmap(t6_daf, hm2 = NULL,
                   k = "meta20",
                   # m = "flo_merge",
                   cluster_anno = TRUE,
                   draw_freqs = TRUE,
                   scale=T
)
