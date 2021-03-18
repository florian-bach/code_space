# preamble ####

library(cowplot)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(gridBase)
library(grid)
library(dplyr)
library(tidyr)
library(ggplot2)

`%notin%` <- Negate(`%in%`)

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


# first import some stuff that so that we don't have to keep defining everything, even if we want to just run a chunk
stacked_bar_data <-  read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_all_cluster_freqs.csv")
stacked_bar_data$volunteer <- factor(stacked_bar_data$volunteer, levels=c("v313", "v315", "v320", "v306", "v301", "v308", "v305", "v304", "v310"))
stacked_bar_data$cluster_id <- gsub("activated CD8lo Effector", "activated DN", stacked_bar_data$cluster_id)
stacked_bar_data$timepoint <- gsub("DoD", "Diagnosis", stacked_bar_data$timepoint)
stacked_bar_data$timepoint <- factor(stacked_bar_data$timepoint, levels=c("Baseline", "Diagnosis", "T6", "C45"))

# add column for lineage, relying on cluster_id names
stacked_bar_data$lineage <- substr(stacked_bar_data$cluster_id, nchar(stacked_bar_data$cluster_id)-2, nchar(stacked_bar_data$cluster_id))
#[1] "CD4" " gd" "AIT" "NKT" "CD8" " DN" "reg"

lin_replacement <- setNames(c("CD4", "gd", "MAIT", "NKT", "CD8", "DN", "Treg", "gd"), unique(stacked_bar_data$lineage))
stacked_bar_data$lineage <- stringr::str_replace_all(stacked_bar_data$lineage, lin_replacement)
stacked_bar_data$lineage <- factor(stacked_bar_data$lineage, levels=rev(c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")))



#all significant clusters at T6 in primaries (all in tertiaries are contained within as well)
ter_sig_clusters <- scan("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/ter_sig_t6_clusters.txt", what = "")
prim_sig_clusters <- scan("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/prim_sig_t6_clusters.txt", what = "")

cluster_order <- scan("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/cluster_order.csv", what="")

#colour business
colcsv <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_colours.csv")
col_pal <- colcsv$colour
names(col_pal) <- colcsv$cluster_id

lineage_palette <- c("#AA3377", "#EE6677", "#4477AA", "#66CCEE", "#228833", "#FFA500", "pink")
names(lineage_palette) <-c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")

#vol_pal <- colorspace::qualitative_hcl(n=9, palette = "dark3")[c(1,9,8,2:7)]
vol_pal <-c("#FFC800",
            "#EE000C",
            "#FF8000",
            "#0080FF",
            "#009B95",
            "#E54787",
            "#8000FF",
            "#66CCFF",
            "#4B0055")

names(vol_pal) <- c("v313", "v315", "v320", "v306", "v301", "v308", "v305", "v304", "v310")



para_vol_pal <- vol_pal
names(para_vol_pal ) <- c("313", "315", "320", "1039", "1040", "1061", "1068", "1075", "6032")

lymph_vol_pal <- para_vol_pal
names(lymph_vol_pal) <- paste("v", names(para_vol_pal), sep="")


time_col <- colorspace::sequential_hcl(5, palette = "Purple Yellow")

inferno_white <- c("#FFFFFF", colorspace::sequential_hcl("inferno", n=8))


UMAP_theme <- theme_minimal()+theme(
  panel.grid.minor = element_blank(),
  legend.position = "none",
  axis.text = element_blank()
)





# barcharts & lollipops ####


# only display DA clusters
activated_barchart_data <- subset(stacked_bar_data, stacked_bar_data$cluster_id %in% prim_sig_clusters)

# caluclate mean frequencies for barcharts & lollipops; scaled_freq=divide freq by number of volunteers so that
# the sum of each group's volunteeers' freqs together equal arthmetic mean
summary <- activated_barchart_data %>%
  group_by(timepoint, cluster_id) %>%
  mutate(mean=base::mean(frequency), sd=stats::sd(frequency)) %>%
  ungroup()%>%
  mutate(timepointf = factor(timepoint, levels=c("Baseline", "Diagnosis", "T6", "C45"))) %>%
  mutate(n_infection = ifelse(.$volunteer %in% c("v313", "v315", "v320"), "first", "third")) %>%
  mutate(scaled_freq = ifelse(.$n_infection=="third",.$frequency/6, .$frequency/3))




summary_t6 <- filter(summary, timepoint=="T6")


# only display first infection
prim_summary_t6 <- subset(summary_t6, summary_t6$n_infection=="first" & summary_t6$cluster_id %in% prim_sig_clusters)

#determine cluster order based on size across all samples
cluster_order <-  prim_summary_t6 %>%
  group_by(cluster_id) %>%
  mutate("sum"=sum(frequency)) %>%
  arrange(desc(sum)) %>%
  select(cluster_id)

cluster_order <- rev(unique(cluster_order$cluster_id))



activation_stacked_barchart_cluster <- ggplot(summary, aes(x=volunteer, y=frequency/100, fill=factor(cluster_id, levels=c(cluster_order))))+
  geom_bar(stat="identity", position="stack")+
  theme_minimal()+
  facet_wrap(~timepointf, strip.position = "bottom", ncol=4)+
  ggtitle("Overall T cell activation")+
  scale_fill_manual(values=col_pal)+
  scale_y_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid.minor.y = element_blank(),
        strip.placement = "outside")

ggsave("/home/flobuntu/PhD/figures_for_thesis/chapter_03/vac63c_activation_stack_cluster_id.pdf", activation_stacked_barchart_cluster, height=4, width=8)


cd4_data <- subset(summary, grepl("CD4", summary$cluster_id))

cd4_activation_stacked_barchart <- ggplot(cd4_data, aes(x=volunteer, y=frequency/100, fill=factor(cluster_id, levels=c(cluster_order))))+
  geom_bar(stat="identity", position="stack")+
  theme_minimal()+
  facet_wrap(~timepointf, strip.position = "bottom", ncol=4)+
  ggtitle("Overall T cell activation")+
  scale_fill_manual(values=col_pal)+
  scale_y_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid.minor.y = element_blank(),
        strip.placement = "outside")

ggsave("/home/flobuntu/PhD/figures_for_thesis/chapter_03/cd4_activation_stack.pdf", cd4_activation_stacked_barchart, height=4, width=8)




all_activated_data <- subset(stacked_bar_data, grepl("*activated*", stacked_bar_data$cluster_id))


all_summary <- all_activated_data %>%
  mutate(timepointf = factor(timepoint, levels=c("Baseline", "Diagnosis", "T6", "C45"))) %>%
  mutate(n_infection = ifelse(.$volunteer %in% c("v313", "v315", "v320"), "first", "third")) %>%
  mutate(scaled_freq = ifelse(.$n_infection=="third",.$frequency/6, .$frequency/3))

all_activation_stacked_barchart_lineage <- ggplot(all_summary, aes(x=volunteer, y=frequency/100, fill=lineage))+
  geom_bar(stat="identity", position="stack")+
  theme_minimal()+
  facet_wrap(~timepointf, strip.position = "bottom", ncol=4)+
  scale_fill_manual(values=lineage_palette)+
  scale_y_continuous(name = "Percentage of CD3+ T cells activated\n", labels=scales::percent_format(accuracy = 1))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=13, face = "bold"),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_text(angle = 45, hjust=1, size=13),
        panel.grid.minor.y = element_blank(),
        strip.placement = "outside")

ggsave("/home/flobuntu/PhD/figures_for_thesis/chapter_03/vac63c_all_activation_stack_lineage.pdf", all_activation_stacked_barchart_lineage, height=4, width=8)



var_cluster_order <- unique(summary$cluster_id[order(summary$lineage)])

activation_stacked_barchart_cluster_var <- ggplot(summary, aes(x=volunteer, y=frequency/100, fill=factor(cluster_id, levels = var_cluster_order)))+
  geom_bar(stat="identity", position="stack")+
  theme_minimal()+
  facet_wrap(~timepointf, strip.position = "bottom", ncol=4)+
  scale_fill_manual(values=col_pal)+
  scale_y_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1))+
  guides(fill=guide_legend(override.aes = list(size=0.5), ncol = 1, keywidth = 0.5))+
  theme(plot.title = element_text(hjust=0.5, size=14),
        strip.text = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        legend.position="none",
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.placement = "outside")



all_activation_stacked_barchart_lineage_var <- all_activation_stacked_barchart_lineage+theme(legend.position = "bottom",
                                                                                             legend.direction = "horizontal",
                                                                                             legend.title = element_blank(),
                                                                                             axis.title = element_blank())+guides(fill=guide_legend(override.aes = list(size=0.5), nrow=1))


combo_cluster_lineage_plot <- plot_grid(activation_stacked_barchart_cluster_var, all_activation_stacked_barchart_lineage_var, ncol=1, rel_heights = c(1,1.2), align = "v", axis="ltbr")

acti.grob <- textGrob("Percentage of CD3+ T cells Activated", 
                      gp=grid::gpar(fontsize=14), rot = 90)

#add to plot
combo_cluster_lineage_plot <- gridExtra::grid.arrange(gridExtra::arrangeGrob(combo_cluster_lineage_plot, left = acti.grob))



ggsave("~/PhD/figures_for_thesis/chapter_03/combo_cluster_lineage_plot.pdf", combo_cluster_lineage_plot, height=9, width = 9)







#make dataframe for third infection and add channel that only contains color hex codes for sig clusters, otherwise paste blank
ter_summary_t6 <- subset(summary_t6, summary_t6$n_infection=="third" & summary_t6$cluster_id %in% prim_sig_clusters)
ter_summary_t6$fill <- ifelse(ter_summary_t6$cluster_id %in% ter_sig_clusters, ter_summary_t6$cluster_id, "blank")


prim_sig_activation_lolli <- ggplot()+
  geom_bar(data=prim_summary_t6,aes(y=factor(cluster_id, levels=cluster_order), x=scaled_freq/100, fill=factor(cluster_id, levels=cluster_order)), stat="identity", position="stack", width=1)+
  theme_minimal()+
  ggtitle("differentially UP at T6 first")+
  scale_fill_manual(values=col_pal)+
  scale_x_reverse(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1), limits = c(0.05, 0), breaks=seq(1,5)/100, expand = c(0,0))+
  guides(fill=guide_legend(reverse = TRUE, nrow = 5,keyheight = unit(2, "mm"), keywidth = unit(4, "mm")))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.placement = "outside",
        plot.margin = unit(c(1,0,1,5), "mm"),
        legend.position = "bottom",
        legend.text = element_text(size=6),
        legend.title = element_blank())


cluster_leg <- cowplot::get_legend(prim_sig_activation_lolli)

prim_sig_activation_lolli <- prim_sig_activation_lolli+theme(legend.position = "none")


ter_sig_activation_lolli <-   ggplot()+
  geom_bar(data=ter_summary_t6,aes(y=factor(cluster_id, levels=cluster_order), x=scaled_freq/100, fill=factor(fill, levels=cluster_order)), stat="identity", position="stack", width=1)+
  theme_minimal()+
  ggtitle("differentially UP at T6 third")+
  scale_fill_manual(values=col_pal)+
  scale_x_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1), limits = c(0,0.05), breaks=seq(1,5)/100 ,expand = c(0,0))+
  guides(fill=guide_legend(reverse = TRUE, nrow = 5,keyheight = unit(2, "mm"), keywidth = unit(4, "mm")))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.placement = "outside",
        plot.margin = unit(c(1,5,1,0), "mm"),
        legend.position = "none",
        legend.text = element_text(size=6),
        legend.title = element_blank())



combo_lolli <- cowplot::plot_grid(prim_sig_activation_lolli, ter_sig_activation_lolli)

combo_lolli <- cowplot::plot_grid(combo_lolli, cluster_leg, ncol=1, rel_heights = c(5,1), rel_widths = c(1,1,2),align="v", axis = "tbrl")
ggsave("/home/flobuntu/PhD/figures_for_thesis/chapter_03/vac63c_sig_activation_lollipop.pdf", combo_lolli, height=4, width=7)


# vivax vs falciparum barcharts ####

falci_activation_data <- all_summary

vivax_activation_data <- data.table::fread("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/activation_barchart_data")
vivax_activation_data$volunteer <- gsub("V", "v", vivax_activation_data$volunteer, fixed=T)

vivax_activation_data$timepoint <- gsub("DoD", "Diagnosis", vivax_activation_data$timepoint, fixed=T)

harmonised_colnames <- subset(colnames(vivax_activation_data), colnames(vivax_activation_data) %in% colnames(falci_activation_data))


falci_activation_data <- falci_activation_data %>%
  select(all_of(harmonised_colnames)) %>%
  filter(timepoint %in% c("Baseline", "Diagnosis", "T6"))%>%
  filter(lineage %notin% c("NKT")) %>%
  mutate(species="P. falciparum")

vivax_activation_data <- vivax_activation_data %>%
  mutate(species="P. vivax") %>%
  filter(timepoint %in% c("Baseline", "Diagnosis", "T6"))%>%
  filter(lineage %notin% c("NKT")) %>%
  select(all_of(colnames(falci_activation_data))) 
  

vivax_falci_summary <- rbind(vivax_activation_data, falci_activation_data)

vivax_falci_summary$species <- factor(vivax_falci_summary$species, levels=c("P. vivax", "P. falciparum"))

(vivax_falci_activation_stacked <- ggplot(vivax_falci_summary, aes(x=volunteer, y=frequency/100, fill=lineage))+
  geom_bar(stat="identity", position="stack")+
  theme_minimal()+
  facet_grid(timepoint~species, scales = "free_x")+
  scale_fill_manual(values=lineage_palette)+
  scale_y_continuous(name = "Percentage of CD3+ T cells activated\n", labels=scales::percent_format(accuracy = 1))+
  guides(fill=guide_legend(nrow = 1, title="Lineage"))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        legend.position="bottom",
        legend.direction = "horizontal",
        strip.text.x = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid.minor.y = element_blank(),
        strip.placement = "outside"))

ggsave("/home/flobuntu/PhD/figures_for_thesis/chapter_03/vivax_falci_activation_stack.pdf", vivax_falci_activation_stacked, height=7, width=6)


res <- vivax_falci_summary %>%
  filter(timepoint=="T6")%>%
  group_by(volunteer, lineage)%>%
  summarise("perc"=sum(frequency))

res %>%
  group_by(species) %>%
  summarise(mean(perc))



  # pie charts ####

#make the pie chart for primaries

prim_map_data <- prim_summary_t6 %>%
  group_by(cluster_id, timepoint) %>%
  mutate(mean_freq=mean(frequency)) %>%
  ungroup() %>%
  select(cluster_id, n_infection, mean_freq, lineage)

prim_pie_data <- prim_map_data[!duplicated(prim_map_data), ]
prim_pie_data <- prim_pie_data[order(prim_pie_data$lineage, decreasing = T),]

#make the colour palette so that it's in the right order
ordered_col_pal <- col_pal[match(prim_pie_data$cluster_id, names(col_pal))]




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

# this retrieves the cluster_ids in the same order in which they were plotted so that the order in the legend is the
# same; then we make the legend palette based on that

prim_lvl_data <- prim_pie_data$cluster_id

prim_leg_col_pal <- subset(col_pal, names(col_pal) %in% prim_lvl_data)

# make the legend
prim_lgd_cluster= Legend(at = prim_lvl_data, type = "grid", 
                         legend_gp = gpar(fill = prim_leg_col_pal[prim_lvl_data], size=0.5),
                         title_position = "topleft",
                         labels_gp = gpar(fontsize = 5),
                         title = "Cluster ID")


# make the pie plot for tertiaries #


ter_map_data <- stacked_bar_data %>%
  mutate(n_infection = ifelse(.$volunteer %in% c("v313", "v315", "v320"), "first", "third")) %>%
  filter(timepoint=="T6", n_infection=="third") %>%
  filter(cluster_id %in% ter_sig_clusters) %>%
  group_by(cluster_id, timepoint) %>%
  mutate(mean_freq=mean(frequency)) %>%
  ungroup() %>%
  select(cluster_id, n_infection, mean_freq, lineage)

ter_pie_data <- ter_map_data[!duplicated(ter_map_data), ]
ter_pie_data <- ter_pie_data[order(ter_pie_data$lineage, decreasing = T),]


#circos.par(cell.padding = c(0, 0))


ter_circlize_plot = function() {
  
  
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

#save a plot with both pies and the legend from the first one
pdf("/home/flobuntu/PhD/figures_for_thesis/chapter_03/vac63c_cluster_activation_pies.pdf", width=8, height=8, onefile = FALSE)

plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridBase::gridOMI(), mfrow = c(1, 3))
prim_circlize_plot()
title("First Infection",  line = -15)
ter_circlize_plot()
title("Third Infection", line = -15)
upViewport()
draw(prim_lgd_cluster, x = circle_size, just = "right")
dev.off()




# all frequency plots ####

equal_breaks <- function(n = 3, s = 0.05, ...){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1+2*s)
    seq(min(x)+d, max(x)-d, length=n)
  }
}

all_freq_data <- data.frame(stacked_bar_data, "n_infection" = ifelse(stacked_bar_data$volunteer %in% c("v313", "v315", "v320"), "first", "third")) 


all_freq_data$cluster_id <- gsub(pattern = "activated CD27-ICOS-HLADR+ EM CD4", replacement = "activated CD27- ICOS-HLADR+ EM CD4",
                                 fixed = TRUE,
                                 all_freq_data$cluster_id) 

gg_prim_sig_clusters <- gsub(pattern = "activated CD27-ICOS-HLADR+ EM CD4", replacement = "activated CD27- ICOS-HLADR+ EM CD4",
                             fixed = TRUE, prim_sig_clusters)

all_cluster_freqs <- ggplot(all_freq_data, aes(x=timepoint, y=frequency/100))+
  geom_boxplot(aes(fill=n_infection), outlier.shape = NA)+
  geom_point(aes(group=n_infection, colour=volunteer), position = position_dodge(width = 0.75), size=0.5)+
  facet_wrap(~cluster_id, ncol=7, scales = "free_y", labeller=label_wrap_gen(width = 12))+
  theme_minimal()+
  scale_x_discrete(name = "Timepoint")+
  scale_y_continuous(name = "Fraction of CD3+ T cells", label=scales::label_percent(), breaks=equal_breaks(n=4, s=0.05))+
  scale_fill_manual(values = c("Baseline"=time_col[4], "first"=time_col[2], "third"=time_col[1], "C45"=time_col[5]))+
  guides(color=guide_legend(nrow = 1, title = NULL, order = 1, keywidth = 0.5, override.aes = list(size = 2)),
         fill=guide_legend(title=NULL, order = 2))+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=5),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        #legend.box = "vertical",
        legend.spacing.x = unit(1, "mm"),
        strip.text = element_text(size=5.5),
        plot.margin = margin(c(2,2,2,2)))


ggsave("~/PhD/figures_for_thesis/chapter_03/vac63c_all_cluster_freqs.pdf", all_cluster_freqs, width=8.3, height=11.7)






all_cluster_counts <- ggplot(all_freq_data, aes(x=timepoint, y=count+1))+
  geom_boxplot(aes(fill=n_infection), outlier.shape = NA)+
  geom_point(aes(group=n_infection, colour=volunteer), position = position_dodge(width = 0.75), size=0.5)+
  facet_wrap(~cluster_id, ncol=7, labeller=label_wrap_gen(width = 12))+
  theme_minimal()+
  scale_x_discrete(name = "Timepoint")+
  scale_y_log10(name = "Number of Cells per Sample",)+
  scale_fill_manual(values = c("Baseline"=time_col[4], "first"=time_col[2], "third"=time_col[1], "C45"=time_col[5]))+
  guides(color=guide_legend(nrow = 1, title = NULL, order = 1, keywidth = 0.5, override.aes = list(size = 2)),
         fill=guide_legend(title=NULL, order = 2))+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=5),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        #legend.box = "vertical",
        legend.spacing.x = unit(1, "mm"),
        strip.text = element_text(size=5.5),
        plot.margin = margin(c(2,2,2,2)))


ggsave("~/PhD/figures_for_thesis/chapter_03/vac63c_all_cluster_log_counts.pdf", all_cluster_counts, width=8.3, height=11.7)



sig_freq_data <- filter(all_freq_data, cluster_id %in% gg_prim_sig_clusters)

sig_freq_data$fill <- ifelse(sig_freq_data$n_infection=="third", ifelse(sig_freq_data$cluster_id %in% ter_sig_clusters, sig_freq_data$n_infection, NA), sig_freq_data$n_infection)

sig_cluster_freqs <- ggplot(sig_freq_data, aes(x=factor(timepoint), y=frequency/100))+
  geom_boxplot(aes(fill=fill), outlier.shape = NA)+
  geom_point(aes(group=n_infection, colour=volunteer), position = position_dodge(width=0.75), size=0.5)+
  facet_wrap(~cluster_id, ncol=6, scales = "free", labeller=label_wrap_gen(width = 12))+
  theme_minimal()+
  scale_x_discrete(name = "Timepoint")+
  scale_y_continuous(name = "log(%) of CD3+ T cells", label=scales::label_percent(accuracy = 0.01), trans = "log10")+
  scale_fill_manual(values = c("first"=time_col[2], "third"=time_col[1]), breaks =  c("first", "third"))+
  guides(color=guide_legend(nrow = 1, title = NULL, order = 1, keywidth = 0.5, override.aes = list(size = 2)),
         fill=guide_legend(title=NULL, order = 2))+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=6),
        axis.text.y = element_text(size=6),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        #legend.box = "vertical",
        legend.spacing.x = unit(1, "mm"),
        strip.text = element_text(size=7),
        plot.margin = margin(c(2,2,2,2)))


ggsave("~/PhD/figures_for_thesis/chapter_03/vac63c_sig_cluster_freqs.pdf", sig_cluster_freqs, width=8.3, height=6)



# phenotypic heatmaps ####

#read in table of cluster medians, reorder channels, fix names
scaled_mat <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_medians_prim_ter_T6.csv", row.names = 1)
reordered_scaled_mat <- scaled_mat[,match(refined_markers, colnames(scaled_mat))]
rownames(reordered_scaled_mat) <- gsub("activated RA-RO- DN", "activated DN", rownames(reordered_scaled_mat))


#subset heatmap to be only all sig clusters
reordered_scaled_mat <- subset(reordered_scaled_mat, rownames(reordered_scaled_mat) %in% prim_sig_clusters)

#make vector that corresponds to whether the clulsters in the heatmap were differentially up in first only or also third
both <- ifelse(rownames(reordered_scaled_mat) %in% ter_sig_clusters, "both", "first only")


combo_left_anno_var <-  rowAnnotation(annotation_name_gp = gpar(fontsize=10),
                                      annotation_name_rot = 45,
                                      gap = unit(1.5, "mm"),
                                      "cluster_id"=prim_sig_clusters,
                                      "n_infection"=both,
                                      show_legend = c(FALSE, TRUE),
                                      show_annotation_name = TRUE,
                                      simple_anno_size = unit(8, "mm"), # width of the significance bar
                                      col=list("n_infection" = c("both"="darkgrey", "first only"="#36454f"),
                                               "cluster_id" = col_pal),
                                      annotation_legend_param = list(n_infection = list(title ="n_infection",
                                                                                        at = c("first only", "both"),
                                                                                        title_gp=gpar(angle=45),
                                                                                        legend_gp = gpar(fill = c("darkgrey","#36454f")),
                                                                                        title_position = "topleft")
                                      )
                                      
)


#make colour palette for heatmap
inferno <- colorspace::sequential_hcl("inferno", n=8)
col_inferno <- circlize::colorRamp2(seq(0,1, by=1/{length(inferno)-1}), inferno)

reordered_scaled_mat <- subset(reordered_scaled_mat, rownames(reordered_scaled_mat) %in% prim_sig_clusters)


sig_cluster_heatmap <- Heatmap(matrix = as.matrix(reordered_scaled_mat),
                               cluster_rows = TRUE,
                               show_heatmap_legend = TRUE,
                               cluster_columns = FALSE,
                               row_names_side = "left",
                               row_dend_side = "right",
                               col = col_inferno,
                               rect_gp = gpar(col = "white"),
                               left_annotation = combo_left_anno_var,
                               column_names_rot = 45,
                               heatmap_legend_param = list(col=col_inferno,
                                                           title = "Normalised Marker Expression",
                                                           title_position = "leftcenter-rot",
                                                           legend_height = unit(6.2, "cm"),
                                                           border = FALSE)
)

pdf("/home/flobuntu/PhD/figures_for_thesis/chapter_03/vac63c_sig_cluster_heatmap.pdf", width=16, height=8)
draw(sig_cluster_heatmap, padding = unit(c(2, 25, 2, 15), "mm"))
dev.off()



big_scaled_mat <- scaled_mat[,match(refined_markers, colnames(scaled_mat))]
rownames(big_scaled_mat) <- gsub("activated RA-RO- DN", "activated DN",rownames(big_scaled_mat))

all_cluster_heatmap <- Heatmap(matrix = as.matrix(big_scaled_mat),
                               cluster_rows = TRUE,
                               show_row_dend = TRUE,
                               show_heatmap_legend = TRUE,
                               cluster_columns = FALSE,
                               row_names_side = "left",
                               row_dend_side = "right",
                               col = col_inferno,
                               rect_gp = gpar(col = "white"),
                               column_names_rot = 45,
                               heatmap_legend_param = list(col=col_inferno,
                                                           title = "Normalised Marker Expression",
                                                           title_position = "leftcenter-rot",
                                                           legend_height = unit(6.2, "cm"),
                                                           border = FALSE)
)

pdf("//home/flobuntu/PhD/figures_for_thesis/chapter_03/vac63c_all_cluster_heatmap.pdf", width=16, height=18)
draw(all_cluster_heatmap, padding = unit(c(2, 25, 2, 15), "mm"))
dev.off()

# spicy umaps ####

#read in UMAP coordinates
big_table <- data.table::fread("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_all_Tcells_with_UMAP.csv")

# quick edits
big_table$timepoint <- gsub("DoD", "Diagnosis", big_table$timepoint)

big_table$timepoint <- factor(big_table$timepoint, levels=c("Baseline", "Diagnosis","T6","C45"))
big_table$flo_label <- gsub("activated RA-RO- DN", "activated DN", big_table$flo_label)


#create columns that that mark cells that are in significant cluster; if they are significant they are completely opaque;
#if not they're a bit transparent
big_table$prim_significant <- ifelse(big_table$flo_label %in% prim_sig_clusters, big_table$flo_label, "black")
big_table$prim_alpha <- ifelse(big_table$flo_label %in% prim_sig_clusters, 1, 0.6)

big_table$ter_significant <- ifelse(big_table$flo_label %in% ter_sig_clusters, big_table$flo_label, "black")
big_table$ter_alpha <- ifelse(big_table$flo_label %in% ter_sig_clusters, 1, 0.6)

shorter_big_table_t6 <- big_table %>%
  filter(timepoint=="T6") %>%
  group_by(n_infection) %>%
  sample_n(9*10^4) %>%
  ungroup()


prim_t6_umap_data <- filter(shorter_big_table_t6, n_infection=="First")
ter_t6_umap_data <- filter(shorter_big_table_t6, n_infection=="Third")


sig_prim_colour_t_cells_umap <- ggplot(prim_t6_umap_data, aes(x=UMAP2, y=UMAP1, color=prim_significant, alpha=prim_alpha))+
  geom_point(shape=".")+
  scale_colour_manual(values = col_pal)+
  UMAP_theme+
  ggtitle("significantly up at T6\nfirst infection (18)")+
  theme(
    plot.title = element_text(size=13, hjust=0.5),
    axis.title.y = element_blank(),
    axis.title = element_text(size=10))
#cowplot::ggsave2("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/figures_for_paper/sig_prim_colour_t_cells_umap.pdf", sig_prim_colour_t_cells_umap, width=3, height=3)



sig_ter_colour_t_cells_umap <- ggplot(ter_t6_umap_data, aes(x=UMAP2, y=UMAP1, color=ter_significant, alpha=ter_alpha))+
  geom_point(shape=".")+
  scale_colour_manual(values = col_pal)+
  ggtitle("significantly up at T6\nthird infection (12)")+
  UMAP_theme+
  theme(
    plot.title = element_text(size=13, hjust=0.5),
    axis.title = element_blank())
#cowplot::ggsave2("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/figures_for_paper/sig_ter_colour_t_cells_umap.pdf", sig_ter_colour_t_cells_umap, width=3, height=3)



all_colour_t_cells_umap <- ggplot(shorter_big_table_t6, aes(x=UMAP2, y=UMAP1, color=flo_label))+
  geom_point(shape=".")+
  scale_colour_manual(values = col_pal)+
  UMAP_theme+
  ggtitle("All Clusters at T6\n")+
  theme(
    plot.title = element_text(size=13, hjust=0.5),
    axis.title = element_text(size=10),
    axis.title.x = element_blank())

cowplot::ggsave2("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/figures_for_paper/all_colour_t_cells_umap.pdf", all_colour_t_cells_umap, width=3, height=3)


combo_plot <- cowplot::plot_grid(all_colour_t_cells_umap, sig_prim_colour_t_cells_umap, sig_ter_colour_t_cells_umap, nrow=1, align = "hv")
ggsave("~/PhD/figures_for_thesis/chapter_03/vac63c_sig_cluster_umap.png", combo_plot, width=9, height=3)


  
  
  down_big_table_t6 <- big_table %>%
    group_by(n_infection, timepoint) %>%
    sample_n(9*10^4) %>%
    ungroup()
  
  
  hex_through_time <- ggplot(down_big_table_t6, aes(x=UMAP2, y=UMAP1))+
    stat_density_2d(aes(fill = ..density..), geom = 'raster', contour = FALSE, n = 1500)+
    stat_density_2d(contour = TRUE, bins=13, color="white", size=0.05)+
    # geom_point(color="black", alpha=0.5)+
    # stat_density_2d(contour = TRUE, bins=13, color="red", size=0.05)+
    theme_minimal()+
    scale_y_continuous(limits=c(-15, 12))+
    scale_x_continuous(limits=c(-15, 12))+
    facet_grid(n_infection~timepoint, switch = "y")+
    theme(legend.position = "none",
          strip.placement = "outside",
          panel.spacing = unit(5, "mm"),
          panel.grid = element_blank(),
          strip.text = element_text(size=12))+
    scale_fill_gradientn(colors=c(inferno_white[c(1, 2, 2:9,9, 9)]))
    
  
system.time(ggsave("/home/flobuntu/PhD/figures_for_thesis/chapter_03/density_umap_var2.png", hex_through_time, height=6, width=12))
  

    
# phenotypic UMAPs & gate_labeled_gg  ####   
    
library(vac69a.cytof)

shortest_big_table_t6 <- big_table %>%
     group_by(n_infection) %>%
     sample_n(2*10^4) %>%
     ungroup()
    


supp_theme <- theme(axis.title = element_text(size = 6),
                    legend.title = element_text(size = 6),
                    legend.text = element_text(size=6))


cd4_plot<- vac63c_flo_umap(shortest_big_table_t6, "CD4", facet_by = "n_infection")+supp_theme+ggtitle("CD4")
cd8_plot<- vac63c_flo_umap(shortest_big_table_t6, "CD8", facet_by = "n_infection")+supp_theme+ggtitle("CD8")
vd2_plot<- vac63c_flo_umap(shortest_big_table_t6, "Vd2", facet_by = "n_infection")+supp_theme+ggtitle("Vd2")
va72_plot<- vac63c_flo_umap(shortest_big_table_t6, "Va72", facet_by = "n_infection")+supp_theme+ggtitle("Va72")
foxp3_plot<- vac63c_flo_umap(shortest_big_table_t6, "FoxP3", facet_by = "n_infection")+supp_theme+ggtitle("FoxP3")
cd27_plot<- vac63c_flo_umap(shortest_big_table_t6, "CD27", facet_by = "n_infection")+supp_theme+ggtitle("CD27")
ccr7_plot<- vac63c_flo_umap(shortest_big_table_t6, "CCR7", facet_by = "n_infection")+supp_theme+ggtitle("CCR7")
cd45ra_plot<- vac63c_flo_umap(shortest_big_table_t6, "CD45RA", facet_by = "n_infection")+supp_theme+ggtitle("CD45RA")
cd38_plot<- vac63c_flo_umap(shortest_big_table_t6, "CD38", facet_by = "n_infection")+supp_theme+ggtitle("CD38")
bcl2_plot<- vac63c_flo_umap(shortest_big_table_t6, "BCL2", facet_by = "n_infection")+supp_theme+ggtitle("Bcl2")

vac63c_pheno_umap <- cowplot::plot_grid(cd4_plot, cd8_plot, vd2_plot,
                                        va72_plot, foxp3_plot,  ccr7_plot,
                                        cd45ra_plot, cd27_plot, cd38_plot, bcl2_plot, align = "hv", axis = "tblr", ncol = 2)
 
ggsave("/home/flobuntu/PhD/figures_for_thesis/chapter_03/vac63c_pheno_umap.png", vac63c_pheno_umap, height = 10, width=8)      



gate_label <- ggplot(shortest_big_table_t6, aes(x=UMAP2, y=UMAP1))+
  geom_point(shape=".", color = "black", alpha=0.8)+
  scale_colour_manual()+
  UMAP_theme+
  ggtitle("Major T Cell Lineages\n")+
  theme(
    plot.title = element_text(size=13, hjust=0.5),
    axis.title = element_text(size=10))

cowplot::ggsave2("/home/flobuntu/PhD/figures_for_thesis/chapter_03/gate_unlabeled.png", gate_label, width=4, height=4)


# ds_limma   ####
# 

num_wide_scaled_ms <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/sample_wise_cluster_expression_matrix.csv", header=T, row.names = 1)

# ter_dod_ms <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/ds_limma/ter_dod_ms.csv", header=T, row.names = 1)
# prim_t6_ms <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/ds_limma/prim_t6_ms.csv", header=T, row.names = 1)
# ter_t6_ms <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/ds_limma/ter_t6_ms.csv", header=T, row.names = 1)
# 
#
num_wide_scaled_ms <- ter_dod_ms



rownames(num_wide_scaled_ms) <- gsub("CD8 Effector", "Effector CD8", rownames(num_wide_scaled_ms))
rownames(num_wide_scaled_ms) <- gsub("CD8 CM", "CM CD8", rownames(num_wide_scaled_ms))

colnames(num_wide_scaled_ms) <- gsub("DoD", "Diagnosis",colnames(num_wide_scaled_ms), fixed = TRUE)


limma_col_split <- factor(rep(c("Baseline", "Diagnosis", "T6", "C45"), each=9), levels=c(c("Baseline", "Diagnosis", "T6", "C45")))

open_bracket_positions <- sapply(rownames(num_wide_scaled_ms), function(x) regexpr(pattern = "*\\(", text=x))
limma_row_split_lineage <- substr(rownames(num_wide_scaled_ms),open_bracket_positions-5, open_bracket_positions-2)
sig_limma_markers <- unique(substr(rownames(num_wide_scaled_ms),open_bracket_positions+1, nchar(rownames(num_wide_scaled_ms))-1))

limma_row_split_lineage <- limma_row_split_lineage %>%
  gsub(" ", "", ., fixed = TRUE)


num_wide_scaled_ms <- num_wide_scaled_ms[,paste(rep(c("v313", "v315", "v320", "v306", "v301", "v308", "v305", "v304", "v310"),times=4), rep(c("Baseline", "Diagnosis", "T6", "C45"), each=9), sep="_")]


colnames(num_wide_scaled_ms) <- gsub("_", " ", colnames(num_wide_scaled_ms))

num_wide_scaled_ms <- ifelse(num_wide_scaled_ms>abs(min(num_wide_scaled_ms)), abs(min(num_wide_scaled_ms)), as.matrix(num_wide_scaled_ms))

col_fun_ds_limma <- circlize::colorRamp2(c(min(num_wide_scaled_ms), 0, max(num_wide_scaled_ms)), c("#0859C6", "black", "#FFA500"))
#

top_anno <-  columnAnnotation(annotation_name_gp = gpar(fontsize=10),
                                      # annotation_name_rot = 45,
                                      gap = unit(1.5, "mm"),
                                      "volunteer"=rep(c("v313", "v315", "v320", "v306", "v301", "v308", "v305", "v304", "v310"),times=4),
                                      "n_infection"=rep(c("first", "first", "first", "third", "third", "third", "third", "third", "third"),times=4),
                                      show_legend = c(FALSE, FALSE),
                                      show_annotation_name = TRUE,
                                      simple_anno_size = unit(4, "mm"), # width of the significance bar
                                      col=list("n_infection" = c("third"="darkgrey", "first"="#36454f"),
                                               "volunteer" = vol_pal),
                                      annotation_legend_param = list(n_infection = list(title ="n_infection",
                                                                                        at = c("first", "third"),
                                                                                        #title_gp=gpar(angle=45),
                                                                                        legend_gp = gpar(fill = c("#36454f", "darkgrey")),
                                                                                        title_position = "topleft")
                                      )

)



  limma_leg = Legend(at = c("v313", "v315", "v320", "v306", "v301", "v308", "v305", "v304", "v310"),
                     type = "grid",
                     legend_gp = gpar(fill = vol_pal, size=1),
                     title_position = "topleft",
                     direction = "horizontal",
                     nrow = 1,
                     labels_gp = gpar(fontsize =9),
                     title = "Volunteer")

  limma_leg2 = Legend(at = c("first", "third"),
                     type = "grid",
                     legend_gp = gpar(fill = c("#36454f", "darkgrey")),
                     title_position = "topleft",
                     direction = "horizontal",
                     nrow = 1,
                     labels_gp = gpar(fontsize =9),
                     title = "N_infection")


  median_cluster_heat <- Heatmap(matrix = num_wide_scaled_ms,
                                 cluster_rows = TRUE,
                                 column_order = 1:36,
                                 name = "Scaled Marker Expression",
                                 cluster_columns = FALSE,
                                 show_row_dend = TRUE,
                                 row_dend_side = "right",
                                 row_names_side = "left",
                                 col = col_fun_ds_limma,
                                 column_names_gp = gpar(fontsize = 11),
                                 column_split = limma_col_split,
                                 top_annotation = top_anno,
                                 split = limma_row_split_lineage,
                                 rect_gp = gpar(col = "white"),
                                 show_parent_dend_line = FALSE,
                                 show_heatmap_legend = TRUE,
                                 column_names_rot = 45,
                                 heatmap_legend_param = list(legend_position = "top",
                                                             col=col_fun_ds_limma,
                                                             title = "Normalised Marker Expression",
                                                             legend_direction = "horizontal",
                                                             title_position = "topcenter",
                                                             legend_width = unit(6.2, "cm"),
                                                             border = FALSE)
                                 #width = unit(16, "cm"),
                                 #height = unit(16*9/28, "cm")
  )


  # pdf("~/PhD/figures_for_thesis/chapter_03/cross_sectional_t6_limma.pdf", width = 10, height=12)
  # draw(median_cluster_heat,
  #      annotation_legend_list = list(limma_leg, limma_leg2),
  #      merge_legends = TRUE, heatmap_legend_side = "bottom")
  # dev.off()


  pdf("~/PhD/figures_for_thesis/chapter_03/ter_dod_limma.pdf", width = 10, height=5)
  draw(median_cluster_heat,
       annotation_legend_list = list(limma_leg, limma_leg2),
       merge_legends = TRUE, heatmap_legend_side = "bottom")
  dev.off()

  # pdf("~/PhD/figures_for_thesis/chapter_03/prim_t6_limma.pdf", width = 10, height=12)
  # draw(median_cluster_heat,
  #      annotation_legend_list = list(limma_leg, limma_leg2),
  #      merge_legends = TRUE, heatmap_legend_side = "bottom")
  # dev.off()
  
  # pdf("~/PhD/figures_for_thesis/chapter_03/ter_t6_limma.pdf", width = 10, height=6)
  # draw(median_cluster_heat,
  #      annotation_legend_list = list(limma_leg, limma_leg2),
  #      merge_legends = TRUE, heatmap_legend_side = "bottom")
  # dev.off()
  

# differential abundance modelling comparison ####


pql_glmm_first_t6 <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/glmmPQL_first_t6.csv", header=T, row.names = 1)
pql_glmm_first_t6$model <- "MASS::glmmPQL"
pql_glmm_first_t6$n_infection <- "first"

pql_glmm_third_t6 <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/glmmPQL_third_t6.csv", header=T, row.names = 1)
pql_glmm_third_t6$model <- "MASS::glmmPQL"
pql_glmm_third_t6$n_infection <- "third"

glmer_first_t6 <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/glmer_first_t6.csv", header=T, row.names = 1)
glmer_first_t6$model <- "lme4::glmer"
glmer_first_t6$n_infection <- "first"

glmer_third_t6 <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/glmer_third_t6.csv", header=T, row.names = 1)
glmer_third_t6$model <- "lme4::glmer"
glmer_third_t6$n_infection <-  "third"

nbzimm_first_t6 <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/nbzimm_first_t6.csv", header=T, row.names = 1)
nbzimm_first_t6$model <- "NBZIMM::glmm.nb"
nbzimm_first_t6$n_infection <- "first"

nbzimm_third_t6 <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/nbzimm_third_t6.csv", header=T, row.names = 1)
nbzimm_third_t6$model <- "NBZIMM::glmm.nb"
nbzimm_third_t6$n_infection <-  "third"

edger_first_t6 <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/edgeR/prim_t6_df_edger.csv", header=T)
edger_first_t6$model <- "edgeR::glmFit" 
edger_first_t6$n_infection <- "first"
edger_first_t6 <- edger_first_t6[,colnames(pql_glmm_third_t6)]

edger_third_t6 <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/edgeR/ter_t6_df_edger.csv", header=T)
edger_third_t6$model <- "edgeR::glmFit"
edger_third_t6$n_infection <- "third"
edger_third_t6 <- edger_third_t6[,colnames(pql_glmm_third_t6)]

big_model_df <- rbind(pql_glmm_first_t6, pql_glmm_third_t6, glmer_first_t6, glmer_third_t6, nbzimm_first_t6, nbzimm_third_t6, edger_first_t6, edger_third_t6)


performance <- big_model_df %>%
  mutate(direction=ifelse(logFC>0, "up", "down"))%>%
  group_by(model, n_infection, direction) %>%
  filter(abs(logFC)>1, p_adj<0.05)%>%
  tally(name = "significant clusters")

write.csv(performance, "/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/glmm_performance.csv", row.names = F)

big_model_df$model <- factor(big_model_df$model, levels=c("lme4::glmer", "NBZIMM::glmm.nb", "MASS::glmmPQL", "edgeR::glmFit"))

cytof_volcano_plot <- ggplot(big_model_df, aes(x=logFC, y=-log10(p_adj)))+
  geom_rect(xmin=1, xmax=10, ymax=-log10(0.05), ymin=20, fill="darkgoldenrod1", alpha=0.01)+
  geom_rect(xmin=-10, xmax=-1, ymax=-log10(0.05), ymin=20, fill="blue", alpha=0.01)+
  geom_point(aes(color=cluster_id), size=0.7)+
  theme_minimal()+
  # geom_vline(xintercept = 1, linetype="dashed")+
  # geom_vline(xintercept = -1, linetype="dashed")+
  # geom_hline(yintercept = -log10(0.05), linetype="dashed")+
  # geom_hline(yintercept = log10(0.05), linetype="dashed")+
  scale_x_continuous(limits=c(-10,10))+
  scale_y_continuous(expand = c(0, 1.1))+
  scale_color_manual(values = col_pal)+
  facet_grid(n_infection~model, switch = "y")+
  theme(legend.position = "none",
        strip.placement = "outside",
        axis.text = element_text(size=7),
        panel.spacing = unit(5, "mm"),
        strip.text.y = element_text(size=11),
        strip.text.x = element_text(size=9))

ggsave("~/PhD/figures_for_thesis/chapter_03/cytof_volcano_plot.pdf", cytof_volcano_plot, width=8, height=5)

# parasitaemia ####

vac63c_parasitaemia <- read.csv("~/PhD/clinical_data/vac63c/VAC063_parasitaemias_all.csv")

long_vac63c_parasitaemia <- gather(vac63c_parasitaemia, vol_id, parasitaemia, colnames(vac63c_parasitaemia)[2:ncol(vac63c_parasitaemia)])

long_vac63c_parasitaemia$Volunteer <- substr(long_vac63c_parasitaemia$vol_id, 1, 5)
long_vac63c_parasitaemia$N_infection <- ifelse(grepl("First", long_vac63c_parasitaemia$vol_id)==T, "First",
                                               ifelse(grepl("Second", long_vac63c_parasitaemia$vol_id)==T, "Second", "Third"))


long_vac63c_parasitaemia$Volunteer <- gsub("X", "", long_vac63c_parasitaemia$Volunteer)
long_vac63c_parasitaemia$Volunteer <- gsub("_", "", long_vac63c_parasitaemia$Volunteer)

long_vac63c_parasitaemia <- long_vac63c_parasitaemia[!is.na(long_vac63c_parasitaemia$parasitaemia),]

long_vac63c_parasitaemia <- long_vac63c_parasitaemia %>%
  filter(Volunteer %in% c("313", "315", "320", "1039", "1040", "1061", "1068", "1075", "6032")) %>%
  filter(N_infection %in% c("First", "Third"))



# para_vol_pal <- c("#8000FF",
#                   "#E54787",
#                   "#0080FF",
#                   "#FF8000",
#                   "#009B95",
#                   "#EE000C",
#                   "#FFC800",
#                   "#66CCFF",
#                   "#4B0055")
# 


vac63c_indie_paras <- ggplot(data=long_vac63c_parasitaemia, aes(x=Timepoint, y=parasitaemia, group=vol_id))+
  geom_line(aes(color=Volunteer))+
  geom_point(aes(fill=Volunteer), shape=21, size=2, stroke=0.1, color="black")+
  scale_y_log10()+
  theme_minimal()+
  scale_colour_manual(values=para_vol_pal)+
  ylab("Parasites / mL")+
  xlab("Days Post Infection")+
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        axis.title.x = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=11))

indie_para_leg <- get_legend(vac63c_indie_paras)

vac63c_indie_paras <- vac63c_indie_paras+theme(legend.position = "none")

vac63c_group_paras <- ggplot(data=long_vac63c_parasitaemia, aes(x=Timepoint, y=parasitaemia, group=vol_id))+
  geom_point(aes(color=N_infection))+
  geom_point(aes(fill=Volunteer), shape=21, alpha=0, size = 2, stroke=0.1, color="black")+
  geom_line(aes(color=N_infection))+
  scale_y_log10()+
  theme_minimal()+
  ylab("Parasites / mL")+
  xlab("Days Post Infection")+
  labs(color="Infection")+
  scale_color_manual(values = c("First"=time_col[2], "Third"=time_col[1]))+
  scale_fill_manual(values=para_vol_pal)+
  guides(fill=guide_legend(override.aes = list(alpha=1)))+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        axis.title.x = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=11))

# 
# group_para_leg <- get_legend(vac63c_group_paras)
# vac63c_group_paras <- vac63c_group_paras+theme(legend.position = "none")
# 
# 
# 
# 
# paras_lgd <- plot_grid(indie_para_leg, group_para_leg, ncol=1, align="hv", axis="tblr")
# vac63c_group_paras <- vac63c_group_paras+theme(legend.position = "none")      


# make common x axis title
days.grob <- textGrob("Days Post Infection", 
                            gp=grid::gpar(fontsize=14))

#add to plot
vac63c_parasites_figure <- plot_grid(vac63c_indie_paras, vac63c_group_paras)
vac63c_parasites_figure <- gridExtra::grid.arrange(gridExtra::arrangeGrob(vac63c_parasites_figure, bottom = days.grob))

ggsave("/home/flobuntu/PhD/figures_for_thesis/chapter_03/vac63c_parasitaemia.pdf", vac63c_parasites_figure, width=8, height=4)

# indie variation: aitchison ####

cytof_data <- stacked_bar_data

# aitchison can't deal with 0s :/ so we have to get rid of those clusters or impute them but the matrix is too sparse..
zero_clusters <- unique(subset(stacked_bar_data$cluster_id, stacked_bar_data$frequency==0))

cytof_data$trans_freq=asin(sqrt(cytof_data$frequency/100))


wide_cytof <- data.frame(cytof_data %>%
                           select(cluster_id, sample_id, frequency) %>%
                           filter(cluster_id %notin% zero_clusters)%>%
                           spread(sample_id, frequency)
)

wide_cytof <- as.matrix(select(wide_cytof, -cluster_id))
class(wide_cytof) <- "double"


cytof_mds <- data.frame(cmdscale(robCompositions::aDist(t(wide_cytof))))

# only activated clusters- individual variation is reduced overall, but still doesn't converge much thorugh time,
# v05 & v09 cluster sperate from the other volutneers


colnames(cytof_mds) <- c("MDS1", "MDS2")

cytof_mds$sample_ID <- rownames(cytof_mds)
cytof_mds$timepoint <- substr(cytof_mds$sample_ID, 6, nchar(cytof_mds$sample_ID))
cytof_mds$timepoint <- gsub("DoD", "Diagnosis", cytof_mds$timepoint)
cytof_mds$timepoint <- factor(cytof_mds$timepoint, levels=c("Baseline", "Diagnosis", "T6", "C45"))

cytof_mds$volunteer <- substr(cytof_mds$sample_ID, 1, 4)
cytof_mds$volunteer <- factor(cytof_mds$volunteer, levels=c("v313", "v315", "v320", "v306", "v301", "v308", "v305", "v304", "v310"))
cytof_mds$n_infection <- ifelse(cytof_mds$volunteer %in% c("v313", "v315", "v320"), "First", "Third") 



(indie_aitchison_cytof <- ggplot(cytof_mds, aes(x=MDS1, y=MDS2))+
  geom_point(aes(colour=n_infection), alpha=0)+
  geom_point(aes(shape=timepoint, fill=volunteer))+
  #ggrepel::geom_label_repel(aes_string(label = "sample_id"), show.legend = FALSE)+ 
  theme_minimal()+
  scale_shape_manual(values = c("Baseline"=21, "C10"=24, "Diagnosis"=22, "T6"=23, "C45"=25))+
  scale_fill_manual(values = vol_pal)+
  scale_color_manual(values = c("First"=time_col[2], "Third"=time_col[1]))+
  guides(fill=guide_legend(title="Volunteer", order = 1, override.aes = list(size = 3, fill = unname(vol_pal), shape=21)),
          colour=guide_legend(title="N_infection", order = 3, override.aes = list(size = 3, alpha=1)),
          shape=guide_legend(title="Timepoint", order = 2, override.aes = list(size = 3)))+
  theme(legend.title = element_text(size=8),
        axis.title = element_blank(),
        legend.position = "right"))


indie_leg <- get_legend(indie_aitchison_cytof)
indie_aitchison_cytof <- indie_aitchison_cytof+theme(legend.position = "none")

n_infection_aitchison_cytof <- ggplot(cytof_mds, aes(x=MDS1, y=MDS2))+
  geom_point(aes(shape=timepoint, color=n_infection, fill=n_infection))+
  theme_minimal()+
  scale_color_manual(values = c("First"=time_col[2], "Third"=time_col[1]))+
  scale_shape_manual(values = c("Baseline"=21, "C10"=24, "Diagnosis"=22, "T6"=23, "C45"=25))+
  scale_fill_manual(values = c("First"=time_col[2], "Third"=time_col[1]))+
  guides(color = guide_legend(title="N_infection", override.aes = list(size = 1)))+
  theme(legend.title = element_text(size=8),
        legend.position = "none",
        axis.title = element_blank())


#add to plot
vac63c_indie_var <- plot_grid(indie_aitchison_cytof, n_infection_aitchison_cytof, nrow = 1, axis = "tblr", align="hv")

#ggsave("~/PhD/figures_for_thesis/chapter_03/Aitchisons_CyTOF.png", vac63c_indie_var2, width=8, height=5)   



# indie variation: median marker expression ####

# how the sausage was made

# 
# library(CATALYST)
# library(diffcyt)
# 
# setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only")
# #setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/debarcoded/")
# 
# fcs <- list.files(pattern = "fcs")
# #fcs <- subset(fcs, !grepl(pattern = "C45", fcs))
# #fcs <- subset(fcs, grepl(pattern = "Baseline", fcs))
# #fcs <- subset(fcs, !grepl(pattern = "DoD", fcs))
# fcs <- subset(fcs, !grepl(pattern = "ctrl", fcs))
# fcs <- subset(fcs, !grepl(pattern = "307", fcs))
# fcs <- subset(fcs, !grepl(pattern = "302", fcs))
# 
# 
# vac63_flowset <- flowCore::read.flowSet(fcs)
# 
# 
# 
# md <- read.csv("vac63c_metadata.csv", header = T)
# 
# md <- subset(md, md$file_name %in% fcs)
# md <- md[order(md$timepoint),]
# 
# 
# #
# panel <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_panel.csv", header=T)
# 
# merged_daf <- prepData(vac63_flowset, panel, md, md_cols =
#                          list(file = "file_name", id = "sample_id", factors = c("volunteer", "timepoint", "n_infection", "batch")))
# 
# # all_markers <- c("CD45", "CD3", "CD3", "CD14", "CD16", refined_markers)
# 
# shplit_cells <- function(x, by) {
#   stopifnot(is.character(by), by %in% colnames(colData(x)))
#   cd <- data.frame(colData(x))
#   cd$cluster_id <- cluster_ids(x, "coarse_merge")
#   dt <- data.table::data.table(cd, i = seq_len(ncol(x)))
#   dt_split <- split(dt, by = by, sorted = TRUE, flatten = FALSE)
#   purrr::map_depth(dt_split, length(by), "i")
# }
# 
# 
# ahgg <- function(x, by, fun = c("median", "mean", "sum")) {
#   fun <- switch(match.arg(fun), 
#                 median = rowMedians, mean = rowMeans, sum = rowSums)
#   cs <- shplit_cells(x, by)
#   pb <- purrr::map_depth(cs, -1, function(i) {
#     if (length(i) == 0) return(numeric(nrow(x)))
#     fun(assay(x, "exprs")[, i, drop = FALSE])
#   })
#   purrr::map_depth(pb, -2, function(u) as.matrix(data.frame(
#     u, row.names = rownames(x), check.names = FALSE)))
# }
# 
# 
# 
# ms <- ahgg(merged_daf[sig_limma_markers, ], by = c("cluster_id", "sample_id"))
# ms <- lapply(ms, reshape2::melt, varnames = c("antigen", "sample_id"))
# ms <- bind_rows(ms, .id = "cluster_id")
# 
# ms$timepoint <- substr(ms$sample_id, 5, nchar(as.character(ms$sample_id)))
# 
# write.csv(ms, "~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/median_expression_on_each_cluster.csv", row.names = F)

ms <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/median_expression_on_each_cluster.csv", header=T)


slimmer_ms <- ms %>%
  mutate("contrast"= paste(antigen, " on ", cluster_id))%>%
  select(contrast, sample_id, value) %>%
  pivot_wider(names_from = sample_id, values_from = value)

slimmer_ms <- data.frame(slimmer_ms)
rownames(slimmer_ms) <- slimmer_ms$contrast

slimmer_ms$contrast <- NULL

mds <- limma::plotMDS(slimmer_ms)
df <- data.frame(MDS1 = mds$x, MDS2 = mds$y)
md <- S4Vectors::metadata(merged_daf)$experiment_info
m <- match(rownames(df), md$sample_id)
df <- data.frame(df, md[m, ])




# do it by sample across all cells
# 
# 
# 
# 
# 
# cs_by_s <- split(seq_len(ncol(merged_daf)), merged_daf$sample_id)
# es <- as.matrix(SummarizedExperiment::assay(merged_daf, "exprs"))
# ms <- vapply(cs_by_s, function(cs) Biobase::rowMedians(es[, cs, drop = FALSE]), 
#              numeric(nrow(merged_daf)))
# rownames(ms) <- rownames(merged_daf)
# 
# #write.csv(ms, "~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/median_expression")
# # state markers, type markers or both
# ms <- subset(ms, rownames(ms)%in%sig_limma_markers)
# 
# mds <- limma::plotMDS(ms, plot = FALSE)
# df <- data.frame(MDS1 = mds$x, MDS2 = mds$y)
# md <- S4Vectors::metadata(merged_daf)$experiment_info
# m <- match(rownames(df), md$sample_id)
# df <- data.frame(df, md[m, ])
df$timepoint <- gsub("DoD", "Diagnosis", df$timepoint)
df$timepoint <- factor(df$timepoint, levels=c("Baseline", "Diagnosis", "T6", "C45"))
df$volunteer <- factor(df$volunteer, levels=c("v313", "v315", "v320", "v306", "v301", "v308", "v305", "v304", "v310"))
#df2 <- filter(df, timepoint %in% c("Baseline", "DoD"))

state_markers_mds_indie <- ggplot(df, aes(x=MDS1, y=MDS2))+
     geom_point(aes(colour=n_infection), alpha=0)+
     geom_point(aes(shape=timepoint, fill=volunteer))+
     theme_minimal()+
     scale_shape_manual(values = c("Baseline"=21, "C10"=24, "Diagnosis"=22, "T6"=23, "C45"=25))+
     scale_fill_manual(values = vol_pal)+
     scale_color_manual(values = c("First"=time_col[2], "Third"=time_col[1]))+
     guides(fill=guide_legend(title="Volunteer", order = 1, override.aes = list(size = 3, fill = unname(vol_pal), shape=21)),
            colour=guide_legend(title="N_infection", order = 3, override.aes = list(size = 3, alpha=1)),
            shape=guide_legend(title="Timepoint", order = 2, override.aes = list(size = 3)))+
     theme(legend.title = element_text(size=8),
           axis.title = element_blank(),
           legend.position = "none")


n_infection_state_markers_mds <- ggplot(df, aes(x=MDS1, y=MDS2))+
  geom_point(aes(shape=timepoint, color=n_infection, fill=n_infection))+
  theme_minimal()+
  scale_color_manual(values = c("First"=time_col[2], "Third"=time_col[1]))+
  scale_shape_manual(values = c("Baseline"=21, "C10"=24, "Diagnosis"=22, "T6"=23, "C45"=25))+
  scale_fill_manual(values = c("First"=time_col[2], "Third"=time_col[1]))+
  guides(color = guide_legend(title="N_infection", override.aes = list(size = 1)))+
  theme(legend.title = element_text(size=8),
        legend.position = "none",
        axis.title = element_blank())

marker_mds <- plot_grid(state_markers_mds_indie, n_infection_state_markers_mds, align="hv", axis = "tblr")
marker.grob <- textGrob("Marker Expression", gp=grid::gpar(fontsize=12))
vac63c_aitch_mds <- grid.arrange(arrangeGrob(marker_mds , top = marker.grob ))

vac63c_indie_var <- plot_grid(indie_aitchison_cytof, n_infection_aitchison_cytof, axis = "tblr", align="hv")
comp.grob <- textGrob("Composition", gp=grid::gpar(fontsize=12))
vac63c_marker_mds <- grid.arrange(arrangeGrob(vac63c_indie_var , top = comp.grob ))

mds.grob <-  textGrob("MDS1", gp=grid::gpar(fontsize=12))
mds.grob2 <-  textGrob("MDS2", gp=grid::gpar(fontsize=12), rot=90)

falci_indie_var <- plot_grid(vac63c_marker_mds, vac63c_aitch_mds,
                             align="hv", axis = "tblr", nrow=2)
falci_indie_var <- grid.arrange(arrangeGrob(falci_indie_var , bottom = mds.grob, left= mds.grob2))

falci_indie_var2 <- plot_grid(falci_indie_var, indie_leg, rel_widths = c(6,1))

ggsave("~/PhD/figures_for_thesis/chapter_03/vac63c_indie_var.png", falci_indie_var2, width=7, height=6)


#lymphocytes ####


vac63c_lymph <- read.csv("~/PhD/clinical_data/vac63c/VAC063_haem_all_sequenced_WBC_real_percent.csv")

vac63c_lymph <- dplyr::filter(vac63c_lymph, Leukocytes=="Lymphocytes")

vac63c_lymph <- select(vac63c_lymph, Volunteer_code, trial_number, timepoint, N_infection, Leukocytes, cell_counts)
vac63c_lymph <- dplyr::filter(vac63c_lymph, timepoint %in% c("C-1", "Diagnosis", "D+6"))
vac63c_lymph <- dplyr::filter(vac63c_lymph, N_infection %in% c("First", "Diagnosis", "Third"))
vac63c_lymph$timepoint <- gsub("D+6", "T6", vac63c_lymph$timepoint, fixed=T)
vac63c_lymph$timepoint <- gsub("C-1", "Baseline", vac63c_lymph$timepoint, fixed=T)

vac63c_lymph$Volunteer_code <- gsub("V", "v", vac63c_lymph$Volunteer_code, fixed=T)
vac63c_lymph<- filter(vac63c_lymph, vac63c_lymph$Volunteer_code %in% names(lymph_vol_pal))

indie_lymph_vac63c <- ggplot(vac63c_lymph, aes(x=factor(timepoint, levels=c("Baseline", "Diagnosis", "T6")), y=cell_counts*1000, group=trial_number))+
  geom_point(aes(color=Volunteer_code))+
  geom_line(aes(color=Volunteer_code))+
  theme_minimal()+
  scale_color_manual(values=lymph_vol_pal)+
  xlab("Timepoint")+
  ylab(expression(Lymphocytes~"/"~mu*L~blood))+
  scale_y_continuous(label=scales::comma)+
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14),
        axis.text.x = element_text(size=12, angle=45, hjust=1),
        plot.margin = unit(c(2,2,2,2), "mm"))


# group_lymph_vac63c <- ggplot(vac63c_lymph, aes(x=factor(timepoint, levels=c("C-1", "Diagnosis", "D+6")), y=cell_counts, group=trial_number))+
#   geom_point(aes(color=vac63c_lymph$N_infection))+
#   geom_line(aes(color=vac63c_lymph$N_infection))+
#   theme_minimal()+
#   xlab("Timepoint")+
#   ylab(bquote('Lymphocytes ('*10^6~cells~'/ mL)'))+
#   labs(color="Infection")
#   


group_lymph_vac63c_box <- ggplot(vac63c_lymph, aes(x=factor(timepoint, levels=c("Baseline", "Diagnosis", "T6")), y=cell_counts*1000, fill=N_infection))+
  geom_point(aes(colour=Volunteer_code), alpha=0)+
  geom_boxplot()+
  theme_minimal()+
  xlab("Timepoint")+
  ylab(bquote('Lymphocytes ('*10^6~cells~'/ mL)'))+
  labs(fill="Infection")+
  scale_y_continuous(label=scales::comma)+
  scale_color_manual(values=lymph_vol_pal)+
  scale_fill_manual(values=c("First"=time_col[2], "Third"=time_col[1]))+
  guides(colour=guide_legend(title="", override.aes = list(alpha=1)))+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size=12, angle=45, hjust=1),
        plot.margin = unit(c(2,2,2,2), "mm"))
 



# lgd <- get_legend(group_lymph_vac63c_box)
# group_lymph_vac63c_box <- group_lymph_vac63c_box + theme(legend.position = "none")

# make common x axis title
timepoint.grob <- textGrob("Timepoint", 
                           gp=gpar(fontsize=14))

#add to plot
vac63c_lymphocytes_figure <- plot_grid(indie_lymph_vac63c, group_lymph_vac63c_box, align="h", axis="tblr")
#vac63c_lymphocytes_figure <- grid.arrange(arrangeGrob(vac63c_lymphocytes_figure, bottom = timepoint.grob))


ggsave("~/PhD/figures_for_thesis/chapter_03/vac63c_lymphocytes_figure.png", vac63c_lymphocytes_figure, width=7, height=3.15)



