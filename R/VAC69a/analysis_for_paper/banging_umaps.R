library(ggcyto)
library(flowCore)
library(CATALYST)
library(ggplot2)
library(viridis)



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
      geom_point(shape=".")+
      scale_color_gradientn(colors = inferno_mega_lite)+
      #facet_wrap(facet_two~facet_one)
      UMAP_theme+
      facet_wrap(facet)+
      ggtitle(color_by)+
      theme(strip.text = facet_title)
      
    )
  
  
}


#downsampling to equal size
sampling_ceiling <- 3000
# Being reproducible is a plus
set.seed(1234)

# sample.int takes a sample of the specified size from the elements of x using either with or without replacement.
smaller_vac69a <- fsApply(vac69a, function(ff) {
  idx <- sample.int(nrow(ff), min(sampling_ceiling, nrow(ff)))
  ff[idx,]  # alt. ff[order(idx),]
})



p <- ggcyto(smaller_vac69a, aes(x=CD3, y=CD56))
  
#   geom_hex(bins=8000)+
#   geom_point(shape='.')+
#   scale_x_flowCore_fasinh()+
#   scale_y_flowCore_fasinh()+
#   theme_minimal()

p+geom_point(shape='.')+
  scale_x_flowCore_fasinh()+
  scale_y_flowCore_fasinh()


smol_daf <- prepData(smaller_vac69a, panel, md, md_cols =
                  list(file = "file_name", id = "sample_id", factors = c("timepoint", "batch", "volunteer")),
                panel_cols = list(channel = "fcs_colname", antigen = "marker_name", class =
                                    "marker_class"))

smol_daf <- cluster(smol_daf, features = refined_markers, xdim = 10, ydim = 10, maxK = 45, seed = 1234)

smol_daf <- filterSCE(smol_daf, timepoint!="C10")


plotClusterHeatmap(smol_daf, hm2 = NULL,
                   k = "meta45",
                   # m = "flo_merge",
                   cluster_anno = TRUE,
                   draw_freqs = TRUE,
                   scale=T
)








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


# umap projections ####
start <- Sys.time(); smol_daf <- runUMAP(smol_daf,
                                         exprs_values = "exprs",
                                         feature_set=refined_markers)
                                        ; end <- Sys.time()

(duration <- round(end - start)) #2min
plotDR(smol_daf, color_by = "CD4")+facet_wrap("timepoint")

# pca 50 / NULL ?

start <- Sys.time(); smol_daf <- runDR(smol_daf, dr="UMAP", features=refined_markers); end <- Sys.time()
(duration <- round(end - start)) #2min

plotDR(smol_daf, color_by = "CD4")+facet_wrap("timepoint")



big_table <- data.frame(t(data.frame(assays(smol_daf)$exprs)))
big_table <- data.frame(cbind(big_table, colData(smol_daf)))

slim_umap <- data.frame(reducedDim(smol_daf, "UMAP"))
colnames(slim_umap) <- c("UMAP1", "UMAP2")

big_table <- data.frame(cbind(big_table, slim_umap), stringsAsFactors = F)


smol_time12 <- ggplot(big_table, aes(x=UMAP1, y=UMAP2))+
  stat_density_2d(aes(fill = after_stat(nlevel)), geom="polygon", bins=12)+
  scale_fill_gradientn(colors = inferno_mega_lite)+
  xlim(c(-12, 12))+
  ylim(c(-10, 10))+
  theme_minimal()+
  facet_wrap(~timepoint)+
  UMAP_theme+
  theme(strip.text = element_text(size=14))


#ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/umap_through_time.png", smol_time12, width=15, height=5.54)
ggsave("/home/flo/PhD/cytof/vac69a/figures_for_paper/umap_through_time.png", smol_time12, width=15, height=5.54)








cd4_plot <- flo_umap(big_table, "CD4")
cd8_plot <- flo_umap(big_table, "CD8")
vd2_plot <- flo_umap(big_table, "Vd2")
va72_plot <- flo_umap(big_table, "Va72")


lineage_plot <- plot_grid(cd4_plot, cd8_plot, vd2_plot, va72_plot, ncol=2)
# ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/lineage_plot_ring.png", lineage_plot)
ggsave("/home/flo/PhD/cytof/vac69a/figures_for_paper/lineage_plot.png", lineage_plot)


cd38_plot <- flo_umap(big_table, "CD38")
hladr_plot <- flo_umap(big_table, "HLADR")
bcl2_plot <- flo_umap(big_table, "BCL2")
cd27_plot <- flo_umap(big_table, "CD27")

activation_plot <- plot_grid(cd38_plot, hladr_plot, bcl2_plot, cd27_plot, ncol=2)
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

 