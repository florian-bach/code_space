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

fcs_name <- "/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/all_events/vac69V05_T+6_comped.fcs"

flow_file <- read.FCS(filename = fcs_name)

lin_gates <- CytoML::cytobank_to_gatingset("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/how_T_cells_were_made_CytExp_293788_Gates_v11.xml", FCS=fcs_name)


dna_gate_plot <- as.ggplot(ggcyto(lin_gates, aes(x=Ir191Di, y=Ce140Di))+
  geom_hex(bins=256)+
  geom_point(color="black", alpha=0.5, shape=".")+
  geom_gate(gs_get_pop_paths(lin_gates)[2])+
  stat_density_2d(contour = TRUE, bins=13, color="white", size=0.1)+
  # geom_text(data=gate_label_positions[7,], aes(x=x_coord, y=y_coord, label=gate_name, size=20), lineheight = 0.7, parse = T)+
  # geom_text(data=gate_label_positions[-7,], aes(x=x_coord, y=y_coord, label=gate_name, size=20), lineheight = 0.7)+
  ggtitle("All Cells")+
  geom_stats()+
  xlab("191Ir")+
  ylab("140Ce")+
  theme_minimal()+
  scale_x_continuous(breaks=seq(0,8,by=2))+
  theme(strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none"))


singlet_gate_plot <- as.ggplot(ggcyto(lin_gates, aes(x=Ir191Di, y=Event_length))+
  geom_point(color="black", alpha=0.5, shape=".", position = position_jitter(height = 1))+
  geom_gate(gs_get_pop_paths(lin_gates)[3])+
  stat_density_2d(contour = TRUE, bins=13, color="white", size=0.1)+
  ggtitle("Singlets")+
  geom_stats()+
    xlab("191Ir")+
    ylab("Event Length")+
    theme_minimal()+
  scale_x_continuous(breaks=seq(0,8,by=2))+
    
  theme(strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5)))

(singlet_gate_plot <- singlet_gate_plot+
  coord_cartesian(xlim = c(0,7), ylim = c(0,100)))


# single_cell_plot <- ggcyto(lin_gates, aes(x=Ir191Di, y=Pt195Di))+
#   geom_point(color="black", alpha=0.5, shape=".")+
#   geom_gate(gs_get_pop_paths(lin_gates)[4])+
#   theme_minimal()+
#   geom_stats()+
#   ggtitle("Singlets Cells")+
#   theme(strip.text = element_blank(),
#         plot.title = element_text(hjust = 0.5))+
#   coord_cartesian(xlim=c(0,8),
#                   ylim=c(0,8),
#                   expand = TRUE)
# 


white_cell_plot <- as.ggplot(ggcyto(lin_gates, aes(x=Ir193Di, y=CD45))+
  geom_point(color="black", alpha=0.5, shape=".")+
  geom_gate(gs_get_pop_paths(lin_gates)[5])+
  theme_minimal()+
  stat_density_2d(contour = TRUE, bins=13, color="white", size=0.1)+
  geom_stats()+
    xlab("193Ir")+
    ylab("CD45")+
  ggtitle("White Cells")+
  theme(strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  coord_cartesian(xlim=c(0,8),
                  ylim=c(0,8)))


t_cell_plot <- as.ggplot(ggcyto(lin_gates, aes(x=CD3, y=CD20))+
  geom_gate(gs_get_pop_paths(lin_gates)[8])+
  geom_point(color="black", alpha=0.5, shape=".")+
  stat_density_2d(contour = TRUE, bins=10, color="white", size=0.1)+
  theme_minimal()+
    xlab("CD3")+
    ylab("CD20")+
  geom_stats()+
  ggtitle("T Cells")+
  theme(strip.text = element_blank(),
  plot.title = element_text(hjust = 0.5))
  )
  

gating_strat <- cowplot::plot_grid(dna_gate_plot, singlet_gate_plot, white_cell_plot, t_cell_plot, nrow=1)

ggsave("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/gating_strategy.png", height=2.2, width=8)
