library(ggcyto)



#downsampling to equal size
sampling_ceiling <- 3500
# Being reproducible is a plus
set.seed(1234)



# sample.int takes a sample of the specified size from the elements of x using either with or without replacement.
smaller_vac69a <- fsApply(vac69a, function(ff) {
  idx <- sample.int(nrow(ff), min(sampling.ceiling, nrow(ff)))
  ff[idx,]  # alt. ff[order(idx),]
})



p <- ggcyto(smaller_vac69a, aes(x=`BCL-2`, y=CD38))
  #geom_hex(bins=8000)+
  #geom_point(shape='.')+
  #stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
  # geom_density2d(colour="black", bins=50, contour =T, size=0.1)+
  # geom_density2d(bins=50, contour=F)+
  stat_density_2d(aes(fill = after_stat(nlevel)), geom = "polygon")+
  #geom_density2d(colour = "black")+
  #stat_density2d(aes(fill = ..density..), geom = "tile")+
  #scale_fill_viridis(option="A", begin = 0, end = 1)+
  scale_fill_viridis(option="A")+
  scale_x_flowCore_fasinh()+
  scale_y_flowCore_fasinh()+
  theme_minimal()


smol_daf <- prepData(smaller_vac69a, panel, md, md_cols =
                  list(file = "file_name", id = "sample_id", factors = c("timepoint", "batch", "volunteer")),
                panel_cols = list(channel = "fcs_colname", antigen = "marker_name", class =
                                    "marker_class"))

smol_daf <- cluster(smol_daf, features = refined_markers, xdim = 10, ydim = 10, maxK = 45, seed = 1234)

smol_daf <- filterSCE(smol_daf, timepoint!="C10")


UMAP_theme <- theme_minimal()+theme(
  panel.grid.major = element_blank(),
  legend.position = "none",
  axis.title = element_blank(),
  plot.title = element_text(hjust = 0.5)
)




start <- Sys.time(); smol_daf <- runUMAP(smol_daf, exprs_values = "exprs", feature_set=refined_markers); end <- Sys.time()
(duration <- round(end - start)) #2min


smol_time <- plotDR(smol_daf, color_by = "CD16")+facet_wrap(~timepoint, ncol = 3)
smol_time$layers[[1]] <- NULL

inferno_mega_lite <- c("#000004", "#8A2267", "#EF802B", "#FFEC89", "#FCFFA4")

smol_time12<-smol_time+
  stat_density_2d(aes(fill = after_stat(nlevel)), geom="polygon", bins=12)+
  scale_fill_gradientn(colors = inferno_mega_lite)+
  xlim(c(-11, 13))+
  ylim(c(-10, 10))+
  theme_minimal()+
  UMAP_theme


plot_grid(smol_time10, smol_time12, smol_time_auto, ncol=1)

flo_umap <- function(sce, color_by, ncol){
  plt <- plotDR(sce, color_by = color_by)
  plt$layers[[1]] <- NULL
  
  plt+
    geom_point(shape="o", size=1.9)+
    theme_minimal()+
    scale_color_viridis(option="A", begin = 0.1)+
    UMAP_theme+
    ggtitle(color_by)
}

cd4_plot <- flo_umap(smol_daf, "CD4")
cd8_plot <- flo_umap(smol_daf, "CD8")
vd2_plot <- flo_umap(smol_daf, "Vd2")
va72_plot <- flo_umap(smol_daf, "Va72")


lineage_plot <- plot_grid(cd4_plot, cd8_plot, vd2_plot, va72_plot, ncol=2)
ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/lineage_plot.png", lineage_plot)

cd38_plot <- flo_umap(smol_daf, "CD38")
hladr_plot <- flo_umap(smol_daf, "HLADR")
bcl2_plot <- flo_umap(smol_daf, "BCL2")
cd27_plot <- flo_umap(smol_daf, "CD27")

activation_plot <- plot_grid(cd38_plot, hladr_plot, bcl2_plot, cd27_plot, ncol=2)
ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/activation_plot.png", activation_plot)

cd45ro_plot <- flo_umap(smol_daf, "CD45RO")
ccr7_plot <- flo_umap(smol_daf, "CCR7")
cd57_plot <- flo_umap(smol_daf, "CD57")
cd127_plot <- flo_umap(smol_daf, "CD127")

memory_plot <- plot_grid(cd45ro_plot, ccr7_plot, cd57_plot, cd127_plot, ncol=2)
ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/memory_plot.png", memory_plot)



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

 