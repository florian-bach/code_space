library(ggcyto)
library(flowCore)
library(CATALYST)
library(ggplot2)
library(viridis)
library(SingleCellExperiment)
library(cowplot)
library(SummarizedExperiment)
library(vac69a.cytof)

#functions, palettes etc. ####

`%!in%` = Negate(`%in%`)

smol_daf <- read_small("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/", proportional=T, event_number=800)
smol_daf <- read_small("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/", proportional=F, event_number=3000)

# plotClusterHeatmap(smol_daf, hm2 = NULL,
#                    k = "meta45",
#                    # m = "flo_merge",
#                    cluster_anno = TRUE,
#                    draw_freqs = TRUE,
#                    scale=T
# )
# 


smol_daf <- filterSCE(smol_daf, timepoint!="C10")

refined_markers <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/refined_markers.csv", stringsAsFactors = F)
# umap projections ####
set.seed(1234);smol_daf <- scater::runUMAP(smol_daf,
                             subset_row=refined_markers[,1],
                             exprs_values = "exprs",
                             scale=T)

big_table <- prep_sce_for_ggplot(smol_daf)

#write.csv(big_table, "equal_small_vac69a_umap.csv")
#write.csv(big_table, "proportional_small_vac69a_umap.csv")

#big_table <- data.table::fread("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/big_table.csv")
inferno_mega_lite <- c("#000004", "#8A2267", "#EF802B", "#FFEC89", "#FCFFA4")

UMAP_theme <- theme_minimal()+theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position = "none",
  axis.title = element_blank()
)

(smol_time12 <- ggplot(big_table, aes(x=UMAP1, y=UMAP2))+
  stat_density_2d(aes(fill = after_stat(level)), geom="polygon", bins=18)+
  xlim(c(-12, 12))+
  ylim(c(-12, 11))+
  theme_minimal()+
  facet_wrap(~timepoint)+
  UMAP_theme+
  theme(strip.text = element_text(size=14))+
  scale_fill_gradientn(colours = inferno_mega_lite))

#ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/proportional_contour_umap.png", smol_time12, height=6, width=9)



#ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/umap_through_time.png", smol_time12, width=15, height=5.54)
#ggsave("/home/flo/PhD/cytof/vac69a/figures_for_paper/umap_through_time.png", smol_time12, width=15, height=5.54)


# (big_time12 <- ggplot(big_table, aes(x=UMAP1, y=UMAP2))+
#   stat_density_2d(aes(fill = after_stat(nlevel)), geom="polygon", bins=25)+
#   scale_fill_gradientn(colors = inferno_mega_lite)+
#   xlim(c(-11, 14))+
#   ylim(c(-9, 10))+
#   theme_minimal()+
#   facet_grid(timepoint~volunteer, switch = "y")+
#   UMAP_theme+
#   theme(strip.text.y.left = element_text(size=14, angle = 0),
#         strip.text.x = element_text(size=14),
#         strip.placement = "outside")

  cd4_plot <- flo_umap(big_table, "CD4")
  cd8_plot <- flo_umap(big_table, "CD8")
  vd2_plot <- flo_umap(big_table, "Vd2")
  va72_plot <- flo_umap(big_table, "Va72")
  
  
(  lineage_plot <- plot_grid(cd4_plot, cd8_plot, vd2_plot, va72_plot, ncol=2))
  # ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/lineage_plot_ring.png", lineage_plot)
  ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/lineage_plot.png", lineage_plot)
  
  
  cd38_plot <- flo_umap(big_table, "CD38")
  hladr_plot <- flo_umap(big_table, "HLADR")
  bcl2_plot <- flo_umap(big_table, "BCL2")
  cd27_plot <- flo_umap(big_table, "CD27")
  
  activation_plot  <- plot_grid(cd38_plot, bcl2_plot,  hladr_plot, cd27_plot, ncol=2)
  # ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/activation_plot.png", activation_plot)
  ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/activation_plot.png", activation_plot)
  
  
  cd45ro_plot <- flo_umap(big_table, "CD45RO")
  ccr7_plot <- flo_umap(big_table, "CCR7")
  cd57_plot <- flo_umap(big_table, "CD57")
  cd28_plot <- flo_umap(big_table, "CD28")
  
  memory_plot <- plot_grid(cd45ro_plot, ccr7_plot, cd57_plot, cd28_plot, ncol=2)
  # ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/memory_plot.png", memory_plot)
  ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/memory_plot.png", memory_plot)
  
  
  
  
  cd25_plot <- flo_umap(big_table, "CD25")
  cd127_plot <- flo_umap(big_table, "CD127")
  foxp3_plot <- flo_umap(big_table, "FoxP3")
  cd39_plot <- flo_umap(big_table, "CD39")
  
  treg_plot <- plot_grid(cd25_plot, cd127_plot, foxp3_plot, cd39_plot, ncol=2)
  # ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/treg_plot.png", treg_plot, width=8, height=5.54)
  ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/treg_plot.png", treg_plot)


  
  
  
  
  cd38_plot_through_time <- flo_umap(big_table, "CD38", "timepoint")
  hladr_plot_through_time <- flo_umap(big_table, "HLADR", "timepoint")
  bcl2_plot_through_time <- flo_umap(big_table, "BCL2", "timepoint")
  cd27_plot_through_time <- flo_umap(big_table, "CD27", "timepoint")
  
  activation_plots_through_time  <- list(cd38_plot_through_time, bcl2_plot_through_time,  hladr_plot_through_time, cd27_plot_through_time)
  # ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/activation_plot.png", activation_plot)
  
  
  lapply(activation_plots_through_time, function(x){
    ggsave(paste("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/", x$labels$title ,"_through_time.png", sep=''), x,
           width = 8, height=4.5)
    })
  

  
  GZB_plot_through_time <- flo_umap(big_table, "GZB", "timepoint")
  perforin_plot_through_time <- flo_umap(big_table, "Perforin", "timepoint")
  Tbet_plot_through_time <- flo_umap(big_table, "Tbet", "timepoint")
  ki67_plot_through_time <- flo_umap(big_table, "Ki67", "timepoint")
  
  intracellular_plots_through_time  <- list(GZB_plot_through_time, perforin_plot_through_time,  Tbet_plot_through_time, ki67_plot_through_time)
  # ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/activation_plot.png", activation_plot)
  
  
  lapply(intracellular_plots_through_time, function(x){
    ggsave(paste("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/", x$labels$title ,"_through_time.png", sep=''), x,
           width = 8, height=4.5)
  })
  

  ctla4_plot_through_time <- flo_umap(big_table, "CTLA4", "timepoint")
  pd1_plot_through_time <- flo_umap(big_table, "PD1", "timepoint")
  CD28_plot_through_time <- flo_umap(big_table, "CD28", "timepoint")
  CD127_plot_through_time <- flo_umap(big_table, "CD127", "timepoint")
  
  exhaustion_plots_through_time  <- list(ctla4_plot_through_time, pd1_plot_through_time,  CD28_plot_through_time, CD127_plot_through_time)
  # ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/activation_plot.png", activation_plot)
  
  
  lapply(exhaustion_plots_through_time, function(x){
    ggsave(paste("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/", x$labels$title ,"_through_time.png", sep=''), x,
           width = 8, height=4.5)
  })
  

  
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

 