library(gridBase)
library(grid)
library(circlize)
library(ComplexHeatmap)
library(dplyr)


stacked_bar_data <-  read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_all_cluster_freqs.csv")

prim_sig_clusters <- read.delim("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/prim_ter_sig_clusters.csv", header=FALSE)$V1

stacked_bar_data$lineage <- substr(stacked_bar_data$cluster_id, nchar(stacked_bar_data$cluster_id)-2, nchar(stacked_bar_data$cluster_id))
#[1] "CD4" " gd" "AIT" "NKT" "CD8" " DN" "reg"

lin_replacement <- setNames(c("CD4", "gd", "MAIT", "NKT", "CD8", "DN", "Treg", "gd"), unique(stacked_bar_data$lineage))
stacked_bar_data$lineage <- stringr::str_replace_all(stacked_bar_data$lineage, lin_replacement)

stacked_bar_data$n_infection <- ifelse(stacked_bar_data$volunteer %in% c("v313", "v315", "v320"), "first", "third")


prim_map_data <- stacked_bar_data %>%
  filter(timepoint=="T6", n_infection=="first") %>%
  filter(cluster_id %in% prim_sig_clusters) %>%
  group_by(cluster_id, timepoint) %>%
  mutate(mean_freq=mean(frequency)) %>%
  ungroup() %>%
  select(cluster_id, n_infection, mean_freq, lineage)

prim_summary <- prim_map_data[!duplicated(prim_map_data), ]



#tail(summary)

prim_summary$lineage <- factor(prim_summary$lineage, levels=rev(c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")))


#taken from khroma bright 6

colcsv <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_colours.csv")

col_pal <- colcsv$colour
names(col_pal) <- colcsv$cluster_id


lineage_palette <- c("#AA3377", "#EE6677", "#4477AA", "#66CCEE", "#228833", "#FFA500", "pink")
names(lineage_palette) <-c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")



prim_pie_data <- prim_summary

ordered_col_pal <- col_pal[match(prim_pie_data$cluster_id, names(col_pal))]



prim_pie_data <- prim_pie_data[!duplicated(prim_pie_data), ]

prim_pie_data <- prim_pie_data[order(prim_pie_data$lineage, decreasing = T),]


prim_circlize_plot = function() {
  
  circos.par("track.height" = 0.6, start.degree = 90)
  circos.par(cell.padding = c(0.02, 0, 0.02, 0))
  
  circos.initialize(factors = prim_pie_data$cluster_id, xlim = c(0,1), sector.width = prim_pie_data$mean_freq)
  
  o.cell.padding = circos.par("cell.padding")
  
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.2, bg.border = NA,
                         cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4]))

  circos.track(factors = prim_pie_data$cluster_id,
               ylim = c(0, 1),
               x=prim_pie_data$mean_freq,
               bg.col = col_pal[match(prim_pie_data$cluster_id, names(col_pal))])
  #circos.genomicLabels()
  
  
  #no color, just text
  # c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")
  
  highlight.sector(sector.index = prim_pie_data$cluster_id[prim_pie_data$lineage=="CD4"], track.index = 1,
                   border = "black", text = "CD4", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = prim_pie_data$cluster_id[prim_pie_data$lineage=="Treg"], track.index = 1,
                   border = "black", text = "Treg", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = prim_pie_data$cluster_id[prim_pie_data$lineage=="CD8"], track.index = 1,
                   border = "black", text = "CD8", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = prim_pie_data$cluster_id[prim_pie_data$lineage=="MAIT"], track.index = 1,
                   border = "black", text = "MAIT", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = prim_pie_data$cluster_id[prim_pie_data$lineage=="gd"], track.index = 1,
                   border = "black", text = "gd", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = prim_pie_data$cluster_id[prim_pie_data$lineage=="DN"], track.index = 1,
                   border = "black", text = "DN", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = prim_pie_data$cluster_id[prim_pie_data$lineage=="NKT"], track.index = 1,
                   border = "black", text = "NKT", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  circos.clear()
  
}

prim_lvl_data <- prim_pie_data$cluster_id

prim_leg_col_pal <- subset(ordered_col_pal, names(ordered_col_pal) %in% prim_lvl_data)
# discrete
prim_lgd_cluster= Legend(at = prim_lvl_data, type = "grid", 
                    legend_gp = gpar(fill = unname(prim_leg_col_pal), size=0.5),
                    title_position = "topleft",
                    labels_gp = gpar(fontsize = 5),
                    title = "Cluster ID")
# discrete
# lgd_lineage = Legend(at = c("CD4", "CD8", "MAIT", "gd"), type = "grid", 
#                      legend_gp = gpar(fill = unname(pie_palette)[1:4]),
#                      title_position = "topleft",
#                      labels_gp = gpar(fontsize = 10, parse = TRUE),
#                      title = "Lineage", )

#lgd_list <- packLegend(lgd_cluster, lgd_lineage)


png("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/ter_cluster_activation_pie.png", width=7, height=4, units = "in", res=400)

plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridBase::gridOMI(), new = TRUE)
circlize_plot()
par(xpd = NA)
upViewport()
draw(lgd_cluster, x = circle_size, just = "left")
#title("Average Size of Significantly Upregulated Clusters at T6", 0,6)

dev.off()



# tertiaries ####


ter_t6_df <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/edgeR/ter_t6_df_edger.csv")

ter_sig_clusters <- ter_t6_df %>%
  filter(p_adj <= 0.05, logFC >=1)%>%
  select(cluster_id)



ter_map_data <- stacked_bar_data %>%
  filter(timepoint=="T6", n_infection=="third") %>%
  filter(cluster_id %in% ter_sig_clusters$cluster_id) %>%
  group_by(cluster_id, timepoint, n_infection) %>%
  mutate(mean_freq=mean(frequency)) %>%
  ungroup() %>%
  select(cluster_id, n_infection, mean_freq, lineage)

ter_summary <- ter_map_data[!duplicated(ter_map_data), ]



#tail(summary)

ter_summary$lineage <- factor(ter_summary$lineage, levels=rev(c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")))


#taken from khroma bright 6

colcsv <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_colours.csv")

col_pal <- colcsv$colour
names(col_pal) <- colcsv$cluster_id


lineage_palette <- c("#AA3377", "#EE6677", "#4477AA", "#66CCEE", "#228833", "#FFA500", "pink")
names(lineage_palette) <-c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")



ter_pie_data <- filter(ter_summary, n_infection=="third")

ordered_col_pal <- col_pal[match(ter_pie_data$cluster_id, names(col_pal))]



ter_pie_data <-ter_pie_data[!duplicated(ter_pie_data), ]

ter_pie_data <- ter_pie_data[order(ter_pie_data$lineage, decreasing = T),]


circos.par(cell.padding = c(0, 0))

sqrt(4.038931)

ter_circlize_plot = function() {
  
  
  # circos.par("canvas.xlim" = c(-2, 2), "canvas.ylim" = c(-2, 2))
  
  circos.par("track.height" = 0.6, start.degree = 90)
  circos.par(cell.padding = c(0.02, 0, 0.02, 0))
  
  circos.initialize(factors = ter_pie_data$cluster_id, xlim = c(0,1), sector.width = ter_pie_data$mean_freq)
  
  o.cell.padding = circos.par("cell.padding")
  
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.2, bg.border = NA,
                         cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4]))
  
  circos.track(factors = ter_pie_data$cluster_id,
               ylim = c(0, 1),
               x=ter_pie_data$mean_freq,
               bg.col = ordered_col_pal[match(ter_pie_data$cluster_id, names(ordered_col_pal))])
  #circos.genomicLabels()
  
  
  #no color, just text
  # c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")
  
  highlight.sector(sector.index = ter_pie_data$cluster_id[ter_pie_data$lineage=="CD4"], track.index = 1,
                   border = "black", text = "CD4", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = ter_pie_data$cluster_id[ter_pie_data$lineage=="CD8"], track.index = 1,
                   border = "black", text = "CD8", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = ter_pie_data$cluster_id[ter_pie_data$lineage=="MAIT"], track.index = 1,
                   border = "black", text = "MAIT", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = ter_pie_data$cluster_id[ter_pie_data$lineage=="gd"], track.index = 1,
                   border = "black", text = "gd", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = ter_pie_data$cluster_id[ter_pie_data$lineage=="DN"], track.index = 1,
                   border = "black", text = "DN", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  highlight.sector(sector.index = ter_pie_data$cluster_id[ter_pie_data$lineage=="NKT"], track.index = 1,
                   border = "black", text = "NKT", col = "white", cex=0.7, niceFacing = TRUE, facing = "clockwise")
  circos.clear()
  
}

ter_circlize_plot()

ter_lvl_data <- ter_pie_data$cluster_id

ter_leg_col_pal <- subset(ordered_col_pal, names(ordered_col_pal) %in% ter_lvl_data)
# discrete
ter_lgd_cluster= Legend(at = ter_lvl_data, type = "grid", 
                         legend_gp = gpar(fill = unname(ter_leg_col_pal)),
                         title_position = "topleft",
                         labels_gp = gpar(fontsize = 5),
                         title = "Cluster ID")


png("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/ter_cluster_activation_pie.png", width=7, height=4, units = "in", res=400)

plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridBase::gridOMI(), new = TRUE)
ter_circlize_plot()
par(xpd = NA)
upViewport()
draw(ter_lgd_cluster, x = circle_size, just = "left")
#title("Average Size of Significantly Upregulated Clusters at T6", 0,6)

dev.off()



  
  
  pdf("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/combo_cluster_activation_pie.pdf", width=8, height=8, onefile = FALSE)
  
  plot.new()
  circle_size = unit(1, "snpc") # snpc unit gives you a square region
  pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                        just = c("left", "center")))
  par(omi = gridBase::gridOMI(), mfrow = c(1, 3))
  prim_circlize_plot()
  ter_circlize_plot()
  #par(new = FALSE)
  upViewport()
  draw(prim_lgd_cluster, x = circle_size, just = "right")
  dev.off()

