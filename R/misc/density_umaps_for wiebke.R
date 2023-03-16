# i can't remember if I put this function into wiebkeTOF but in case I didn't here it is:
# this extracts the UMAP coordinates from the SCE and adds them to a dataframe that contains the expression levels so you can
# really personalise your visualisations

prep_sce_for_ggplot <- function(sce){
  
  #extract event data, metadata and transpose
  df <- data.frame(t(data.frame(SummarizedExperiment::assays(sce)$exprs)))
  df <- data.frame(cbind(df, SingleCellExperiment::colData(sce)))
  
  #extract UMAP coordinates
  umaps <- data.frame(SingleCellExperiment::reducedDim(sce, "UMAP"))
  colnames(umaps) <- c("UMAP1", "UMAP2")
  
  # put it together
  df <- data.frame(cbind(df, umaps), stringsAsFactors = F)
  return(df)
}

big_table <- prep_sce_for_ggplot(merged_sce)

# I use this to write out the table, it's much easier to just load this rather than do all the upstream stuff again
# data.table::fwrite(big_table, "~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_cells_with_UMAP_and_flo_merge_cluster.csv")

# easily downsample the big table for plotting purposes, here we're only keeping every 3rd cell
down_big_table <- big_table[seq(1,nrow(big_table), by=3), ]



# DOPE CONTOUR PLOT
# this took my old laptop between 10-20 minutes to run properly. only one line is the problem and i commented it out
# don't run it at first to have a look at the isobars and if you think they look good, bins=13 was my friend, it may not be yours, experiment

# the n=1500 in the raster graphic command is what makes it so damn slow, every panel contains 1500^2 pixels i.e. many. start at lower numbers to have a look, but i found 1500 worked for me 
# once you do run the raster graphic, DONT attempt to look at it in RStudio. Your session may freeze and crash and it's super annoying. just savve it as png and look at it in preview

hex_through_time <- ggplot(down_big_table, aes(x=UMAP1, y=UMAP2))+
  #stat_density_2d(aes(fill = ..density..), geom = 'raster', contour = FALSE, n = 1500)+
  stat_density_2d(contour = TRUE, bins=13, color="white", size=0.05)+
  scale_x_continuous(limits=c(-13.5, 10))+
  scale_y_continuous(limits=c(-11.2, 11.3))+
  theme_minimal()+
  facet_wrap(~timepoint, ncol = 4)+
  UMAP_theme+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size=11),
        axis.title = element_blank(),
        axis.text = element_blank())+
  scale_fill_gradientn(colors=inferno_white)+
  coord_cartesian(xlim=c(-13.5, 10),
                  ylim=c(-11.2, 11.3))

# system.time will print the amount of time that elapsed to execute the command within it, you don't need it but i find it useful to know how long this step takes, you could optimise this while doing something not computer related because your processor may be very occupied
system.time(ggsave("~/graphs/plot.png", hex_through_time, height=4, width=16))

