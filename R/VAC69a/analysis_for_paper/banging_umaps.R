library(ggcyto)

#downsampling

sampling_ceiling <- 5000
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
  stat_density_2d(aes(fill = after_stat(nlevel)), geom = "polygon"))+
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



UMAP_theme <- theme_minimal()+theme(
  panel.grid.major = element_blank(),
  legend.position = "none",
  axis.title = element_blank()
)



start <- Sys.time(); smol_daf <- runUMAP(smol_daf, exprs_values = "exprs", feature_set=refined_markers); end <- Sys.time()
(duration <- round(end - start)) #2min

smol_time <- plotDR(smol_daf, color_by = "CD38")+facet_wrap(~timepoint, ncol = 2)

smol_time$layers[[1]] <- NULL


smol_time+
     stat_density_2d(aes(fill = stat(nlevel)), geom="polygon", bins=10)+
     #stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon", bins=12, size=0.3, color="black", contour=T)+
     # scale_fill_viridis(option="A")+
     scale_color_viridis(option="B", direction=1)+
     #scale_fill_viridis(option="A", direction=1)+
     scale_fill_gradient(low = "lavenderblush1", high="black")+
     geom_point(shape=".")+
     # xlim(c(-14, 11))+
     # ylim(c(-10, 10))+
     theme_minimal()+
     UMAP_theme
 