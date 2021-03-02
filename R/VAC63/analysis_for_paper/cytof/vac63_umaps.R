library(ggplot2)
library(cowplot)
library(vac69a.cytof)
library(scales)
library(dplyr)


big_table <- data.table::fread("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_all_Tcells_with_UMAP.csv")

t6_edger <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/edgeR/t6_edgeR.csv", header = TRUE, stringsAsFactors = FALSE)
sig_t6_clusters <- subset(t6_edger, t6_edger$p_adj<0.05 & abs(t6_edger$logFC)>1)$cluster_id

#get rid of resting Vd cluster

#bigass color palette

color_103_scheme <- c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                      "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                      "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                      "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                      "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                      "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                      "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                      "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                      "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                      "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                      "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                      "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                      "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")

inferno_mega_lite <- c("#000004", "#8A2267", "#EF802B", "#FFEC89", "#FCFFA4")
sunset_white <- rev(c("#7D1D67", "#A52175", "#CB2F7A", "#ED4572", "#FA716C", "#FF9772", "#FFB985", "#FFD99F", "#FFFFFF"))

# colcsv <- read.csv("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cluster_palette.csv", header=T, stringsAsFactors = F)
# 
# col_pal <- colcsv$x
# names(col_pal) <- colcsv$X


inferno_white <- c("#FFFFFF", colorspace::sequential_hcl("inferno", n=8))

UMAP_theme <- theme_minimal()+theme(
  panel.grid.minor = element_blank(),
  legend.position = "none",
  axis.text = element_blank()
)

big_table$timepoint <- factor(big_table$timepoint, levels=c("Baseline", "DoD","T6","C45"))

all_black_t_cells_umap <- ggplot(big_table, aes(x=UMAP2, y=UMAP1))+
  facet_wrap(n_infection~timepoint, ncol=4)+
  geom_point(color="black", shape="o")+
  stat_density_2d(contour = TRUE, bins=13, color="red", size=0.15)+
  UMAP_theme

ggsave("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/all_black_t_cells_umap.png", width=13, height=7)



inferno_white <- c("#FFFFFF", colorspace::sequential_hcl("inferno", n=8))


all_cd38_t_cells_umap <- ggplot(big_table, aes(x=UMAP2, y=UMAP1, color=CD38))+
  facet_wrap(n_infection~timepoint, ncol=4)+
  geom_point(shape="o")+
  scale_colour_gradientn(colors=inferno_mega_lite)+
  UMAP_theme

ggsave("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/all_cd38_t_cells_umap.png", all_cd38_t_cells_umap, width=13, height=7)

all_bcl2_t_cells_umap <- ggplot(big_table, aes(x=UMAP2, y=UMAP1, color=BCL2))+
  facet_wrap(n_infection~timepoint, ncol=4)+
  geom_point(shape="o")+
  scale_colour_gradientn(colors=inferno_mega_lite)+
  UMAP_theme

ggsave("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/all_bcl2_t_cells_umap.png", all_bcl2_t_cells_umap, width=13, height=7)



refined_markers <- c("CD4",
                     "CD8",
                     "Vd2",
                     "Va72",
                     #"CXCR5",
                     "CD38",
                     #"CD69",
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
                     #"TCRgd",
                     "Tbet",
                     "Eomes",
                     #"RORgt",
                     #"GATA3",
                     "CTLA4",
                     "Ki67",
                     "CD127",
                     "CD56",
                     #"CD16",
                     "CD161",
                     "CD49d",
                     "CD25",
                     "FoxP3",
                     "CD39",
                     "CX3CR1",
                     "CD57",
                     "CD45RA",
                     "CD45RO",
                     "CCR7")


short_big_table <- big_table[seq(1,nrow(big_table), by=3),] 


for(channel in refined_markers){
  plotplot <- flo_umap(short_big_table, color_by = channel, facet_by = c("timepoint", "n_infection"))
  ggsave(paste("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/umaps/t6_n_infection_umap_", channel, ".png", sep = ""), plotplot,  width=13, height=7)
  print(paste("I printed the ", channel, " plot for you :) "))
}



short_big_table$significant <- ifelse(short_big_table$flo_label %in% sig_t6_clusters, short_big_table$flo_label, "black")
short_big_table$alpha <- ifelse(short_big_table$flo_label %in% sig_t6_clusters, 1, 0.5)



cluster_names <- c(unique(short_big_table$flo_label))
cols <- rev(c(color_103_scheme[1:length(cluster_names)]))
names(cols)=cluster_names


#read in pretty divergent colours from vivax paper
colcsv <- read.csv("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cluster_palette.csv", header=T, stringsAsFactors = F)

#find all the activated clusters, replace the colours of the first nine with the pretty vivax colours minus the activated
#gamma delta pink (to avoid duplication)

vivax_colours <- colcsv$x
vivax_colours[4] <- "#6F0062"

vivax_colours <- vivax_colours[-c(1,3)]

cols[grep("activated", names(cols))[2:9]] <- vivax_colours 

names(cols)=cluster_names



cols <- c(cols, "black"="black")


write.csv(data.frame("cluster_id"=names(cols), "colour"=unname(cols)), "~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_colours.csv", row.names =FALSE)


short_big_table_t6 <- short_big_table %>%
  filter(timepoint=="T6") %>%
  group_by(n_infection) %>%
  # sample_n(5*10^4) %>%
  sample_n(0.1*10^4) %>%
  ungroup()




                
sig_colour_t_cells_umap <- ggplot(short_big_table_T6, aes(x=UMAP2, y=UMAP1, color=significant, alpha=alpha))+
  facet_wrap(~n_infection, ncol=2)+
  geom_point(shape="o")+
  scale_colour_manual(values = cols)+
  UMAP_theme+
  ggtitle("Differentially Abundant Clusters at T6")+
  theme(strip.text = element_text(size=15),
      panel.spacing.x = unit(10, "mm"),
      plot.title = element_text(size=17),
      axis.title = element_text(size=13))

#ggsave("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/sig_colour_t_cells_umap.png", sig_colour_t_cells_umap, width=9, height=7)



all_colour_t_cells_umap <- ggplot(short_big_table_T6, aes(x=UMAP2, y=UMAP1, color=flo_label, alpha=alpha))+
  facet_wrap(~n_infection, ncol=2)+
  geom_point(shape="o")+
  scale_colour_manual(values = cols)+
  UMAP_theme+
  ggtitle("All Clusters at T6")+
  theme(strip.text = element_text(size=15),
        panel.spacing.x = unit(10, "mm"),
        plot.title = element_text(size=17),
        axis.title = element_text(size=13))

#ggsave("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/all_colour_t_cells_umap.png", all_colour_t_cells_umap, width=9, height=7)


combo_plot <- plot_grid(all_colour_t_cells_umap, sig_colour_t_cells_umap, nrow=2)
ggsave("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/combo_plot.png", combo_plot, width=7.5, height=8)



