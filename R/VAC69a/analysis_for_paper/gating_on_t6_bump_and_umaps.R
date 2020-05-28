library(ggcyto)
library(flowCore)
library(CATALYST)
library(SingleCellExperiment)
library(cowplot)
library(vac69a.cytof)
library(CytoML)

#functions, palettes etc. ####

`%!in%` = Negate(`%in%`)

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



# GATING ####

# read in data
# setwd("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/post_umap/")
# 
# flz <- list.files(pattern = "*.fcs")
# 
# umap_vac69a <- read.flowSet(flz)
# 
# 
# sampling_ceiling <- 2500
# # Being reproducible is a plus
# set.seed(1234)
# 
# # sample.int takes a sample of the specified size from the elements of x using either with or without replacement.
# smaller_umap_vac69a <- fsApply(umap_vac69a, function(ff) {
#   idx <- sample.int(nrow(ff), min(sampling_ceiling, nrow(ff)))
#   ff[idx,]  # alt. ff[order(idx),]
# })

fcs_name <- "/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/post_umap/concat/all.fcs"

lin_gates <- CytoML::cytobank_to_gatingset("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/CytExp_295996_Gates_simplified.xml", FCS=fcs_name)


#flow_file <- read.FCS(filename = fcs_name)
#lin_gates <- CytoML::read.gatingML.cytobank("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/gating_ml_for_umap.xml")

gate_name <- substr(gs_get_pop_paths(lin_gates)[2:length(gs_get_pop_paths(lin_gates))], 2, 30)
gate_name <- gsub(" ", "\n", gate_name)
gate_name[7] <- "paste(gamma, delta)"
# x_coord <- c(-1.1, -2.5,   5, -10,   -9,    7,    1,   -4.8, -10.6)
# y_coord <- c(5,    -2.4, 7.5,   7,     -8, -5.8, -0.7,    8.3, -4.5)
x_coord <- c(-1.1, -2.5,   5, -10,   -9,    7,    -4.8, -10.6)
y_coord <- c(5,    -2.4, 7.5,   7,    -8,  -5.8,   8.3, -4.5)

gate_label_positions <- data.frame(gate_name, x_coord, y_coord)

(umap_gate_labeled <- ggcyto(lin_gates, aes(x=UMAP1, y=UMAP2))+
  geom_hex(bins=228)+
  geom_gate(c(gs_get_pop_paths(lin_gates)[2:length(gs_get_pop_paths(lin_gates))]))+
  geom_text(data=gate_label_positions[7,], aes(x=x_coord, y=y_coord, label=gate_name, size=20), lineheight = 0.7, parse = T)+
  geom_text(data=gate_label_positions[-7,], aes(x=x_coord, y=y_coord, label=gate_name, size=20), lineheight = 0.7)+
    UMAP_theme+
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme(axis.title = element_text(), 
        panel.grid = element_blank(),
        strip.text = element_blank(),
        plot.title = element_blank()))

gate_labeled_gg <- as.ggplot(umap_gate_labeled)

gate_labeled_gg <- gate_labeled_gg+
  scale_x_continuous(limits=c(-13, 10))+
  scale_y_continuous(limits = c(-11.2, 11.3))+

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/umap_gate_labeled.png", umap_gate_labeled, height=4, width=4)








nodes <- gs_pop_get_children(lin_gates[[1]], "root")
nodes <- nodes[-c(1, 4, 5, 6, 7, 14)]
nodes <- substr(nodes, 2, nchar(nodes))

ggcyto(smaller_umap_vac69a[[1]], aes(x=UMAP2, y=UMAP1))+
  
  #geom_hex(bins=190)+
  # geom_gate(nodes)+
  geom_gate("Naive CD4+")


     geom_stats(type = "gate_name")


#gate that mf

t6_gate <-rectangleGate("UMAP2"=c(-4.5, -2.5), "UMAP1"=c(-2.4, 0))

t6_gate_plot <- ggcyto(umap_vac69a[c(4,8,16, 20, 24)], aes(x=UMAP1, y=UMAP2))+
  #stat_density_2d(aes(fill = after_stat(nlevel)), geom="polygon", bins=25)+
  xlim(c(-11, 14))+
  ylim(c(-10, 10))+
  geom_hex(bins = 190)+
  #geom_point(shape=".")+
  #geom_density()+
  geom_gate(t6_gate)+
  geom_stats(type = "percent", adjust =1.85)+
  theme_minimal()

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/t6_gate_umap.png", t6_gate_plot, height=6, width=9)

baseline_t6_gate_plot <- ggcyto(umap_vac69a[c(4,8,16, 20, 24)-3], aes(x=UMAP1, y=UMAP2))+
  #stat_density_2d(aes(fill = after_stat(nlevel)), geom="polygon", bins=25)+
  xlim(c(-11, 14))+
  ylim(c(-10, 10))+
  geom_hex(bins = 190)+
  #geom_point(shape=".")+
  #geom_density()+
  geom_gate(t6_gate)+
  geom_stats(type = "percent", adjust =1.85)+
  theme_minimal()

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/baseline_t6_gate_plot.png", baseline_t6_gate_plot, height=6, width=9)


t6_bump <- Subset(umap_vac69a, flowCore::filter(umap_vac69a, t6_gate))


#import into CATALYST ####

setwd("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")

md <- read.csv("meta_data.csv", header=T) 

md$timepoint <- factor(md$timepoint, levels = c("Baseline", "C10", "DoD", "T6"))

md$file_name <- as.character(md$file_name)

md$file_name <- paste(substr(md$file_name, 0, nchar(md$file_name)-4),"_post_umap.fcs", sep='')

#import panel 
panel <- read.csv("VAC69_PANEL.CSV", header=T, stringsAsFactors = F)
colnames(panel)[2] <- "marker_name"

#add umap channels into panel
panel <- rbind(panel, c("UMAP1", "UMAP1", "none"), c("UMAP2", "UMAP2", "none"))

#the order of entries in the md file matter, be careful
t6_daf <- prepData(t6_bump, panel, md, md_cols =
                       list(file = "file_name", id = "sample_id", factors = c("timepoint", "batch", "volunteer")),
                       panel_cols = list(channel = "fcs_colname", antigen = "marker_name", class =
                                         "marker_class"))

t6_daf <- cluster(t6_daf, features = cd4_t6_markers, xdim = 10, ydim = 10, maxK = 45, seed = 1234)


plotClusterHeatmap(t6_daf, hm2 = NULL,
                   k = "meta7",
                   # m = "flo_merge",
                   cluster_anno = TRUE,
                   draw_freqs = TRUE,
                   scale=T,
                   palette = inferno_mega_lite)

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/t6_bump_meta7_heatmap.png", unlist(meta7_heatmap))


t6_daf <- filterSCE(t6_daf, timepoint!="C10")
set.seed(1234)

# umap projections ####
t6_daf <- scater::runUMAP(t6_daf,
                            subset_row=cd4_t6_markers,
                            exprs_values = "exprs",
                            scale=T)


t6_bump_meta7 <- plotDR(t6_daf, color_by = "meta7")+facet_grid(~timepoint)
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/t6_bump_meta7.png", t6_bump_meta7)
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


(smol_time12 <- ggplot(big_table, aes(x=UMAP1, y=UMAP2))+
  #stat_density_2d(aes(fill = after_stat(nlevel)), geom="polygon", bins=25)+
  scale_fill_gradientn(colors = inferno_mega_lite)+
  geom_point(shape=".")+
  theme_minimal()+
  facet_wrap(~timepoint)+
  UMAP_theme+
  theme(strip.text = element_text(size=14)))

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/t6_bump_through_time.png", smol_time12)


big_table <- dplyr::select(big_table, -UMAP1, -UMAP2)
colnames(big_table)[c(ncol(big_table)-1,ncol(big_table))] <- c("UMAP1", "UMAP2")

cd38_plot <- flo_umap(big_table, "CD38")
hladr_plot <- flo_umap(big_table, "HLADR")
bcl2_plot <- flo_umap(big_table, "BCL2")
cd27_plot <- flo_umap(big_table, "CD27")

activation_plot  <- plot_grid(cd38_plot, bcl2_plot,  hladr_plot, cd27_plot, ncol=2)


#ggsave("C:/Users/bachf/PhD/cytof/vac69a/figures_for_paper/umap/t6_activation_plot.png", activation_plot)
# ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/activation_plot.png", activation_plot)
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cd4_t6_activation_plot.png", activation_plot)


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
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/t6_bump_plot.png", t6_bump_plot)


ki67_plot <- flo_umap(big_table, "Ki67")
tbet_plot <- flo_umap(big_table, "Tbet")
foxp3_plot <- flo_umap(big_table, "FoxP3")
gzb_plot <- flo_umap(big_table, "GZB")

t6_bump_2_plot <- plot_grid(ki67_plot, tbet_plot, foxp3_plot, gzb_plot, ncol=2)
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/t6_bump_2_plot.png", t6_bump_2_plot)




DrawGate(flz[1], c("UMAP1, UMAP2"), gate_type=polygon, N = 1, axis = "x", adjust = 1.5)



