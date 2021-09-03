# Panel A: Heatmap of all clusters' median expression of all markers ####
library(ggplot2)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)

inferno <- colorspace::sequential_hcl("inferno", n=8)
col_inferno <- circlize::colorRamp2(seq(0,1, by=1/{length(inferno)-1}), inferno)



ms_mean <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/cluster_medians_all.csv", header=T, stringsAsFactors = F)
ms_mean <- select(ms_mean, cluster_id, Marker, mean_expression)


all_mat <- as.matrix(tidyr::spread(ms_mean, cluster_id, mean_expression))

rownames(all_mat) <- all_mat[,1]
all_mat2 <- all_mat[,2:ncol(all_mat)]

scaled_all_mat <- apply(apply(all_mat2, c(1,2), as.numeric), MARGIN = 1, function(x)scales::rescale(x, to=c(0, 1)))

clustering_markers <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/refined_markers.csv", stringsAsFactors = F)


reordered_scaled_mat <- scaled_all_mat[,match(clustering_markers$refined_markers, colnames(scaled_all_mat))]



all_cluster_heatmap <- Heatmap(matrix = reordered_scaled_mat,
                               cluster_rows = T,
                               show_row_dend = FALSE,
                               show_heatmap_legend = TRUE,
                               name = "Median Marker Expression",
                               cluster_columns = FALSE,
                               column_names_gp = gpar(fontsize = 8),
                               row_names_gp = gpar(fontsize = 8),
                               row_names_side = "left",
                               col = col_inferno,
                               rect_gp = gpar(col = "white"),
                               #top_annotation = combo_top_anno,
                               column_names_rot = 45,
                               heatmap_legend_param = list(col = col_inferno,
                                                           #legend_position = "bottom",
                                                           at=c(0,0.5,1),
                                                           title = "Normalised Marker Expression",
                                                           legend_direction = "vertical",
                                                           title_position = "leftcenter-rot",
                                                           legend_height = unit(5.2, "cm"),
                                                           legend_width = unit(0.8, "cm"),
                                                           border = FALSE)
                               #heatmap_legend_param = list(col = col_fun4, title = "Normalised Frequency", title_position = "topleft"),
                               # width = unit(16, "cm"),
                               # height = unit(16*34/28, "cm"
                                             )




pdf("/home/flobuntu/PhD/cytof/vac69a/final_figures_for_paper/supp_marker_expression_heatmap_all_clusters.pdf", width=8, height=7)
draw(all_cluster_heatmap,
     heatmap_legend_side = "right"
)
dev.off()

# Panel B: Series of UMAP projections showing expression of all 10 activation markers

#read in umap projection data
big_table <- data.table::fread("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_cells_with_UMAP_and_flo_merge_cluster.csv")

# restrict data to T6 and downsample by 33% to make it look nice
short_big_table_t6 <- subset(big_table, big_table$timepoint=="T6")
short_big_table_t6 <- short_big_table_t6[seq(1,nrow(short_big_table_t6), by=3), ]

supp_theme <- theme(axis.title = element_text(size = 6),
                    legend.title = element_text(size = 6),
                    legend.text = element_text(size=6),
                    plot.title =  element_text(size=7))

CD27_plot <- flo_umap(short_big_table_t6, "CD27")+supp_theme
                                                                                                       

Ki67_plot <- flo_umap(short_big_table_t6, "Ki67")+supp_theme
                                                                                                       
CD28_plot <- flo_umap(short_big_table_t6, "CD28")+supp_theme
                                                                                                       
PD1_plot <- flo_umap(short_big_table_t6, "PD1")+supp_theme
                                                                                                     
CTLA4_plot <- flo_umap(short_big_table_t6, "CTLA4")+supp_theme
                                                                                                         

Tbet_plot <- flo_umap(short_big_table_t6, "Tbet")+supp_theme
                                                                                                       
Perforin_plot <- flo_umap(short_big_table_t6, "Perforin")+supp_theme
                                                                                                               

GZB_plot <- flo_umap(short_big_table_t6, "GZB")+supp_theme
                                                                                                     


HLADR_plot <- flo_umap(short_big_table_t6, "HLA-DR")+supp_theme
                                                                                                          
ICOS_plot <- flo_umap(short_big_table_t6, "ICOS")+supp_theme
                                                                                                       



supp_fig_2_umaps <- plot_grid(Ki67_plot,  Tbet_plot, CTLA4_plot, PD1_plot,
                         GZB_plot, Perforin_plot, HLADR_plot, CD27_plot, ICOS_plot,  ncol=3, align = "h", axis="trbl")


ggsave("/home/flobuntu/PhD/cytof/vac69a/final_figures_for_paper/supp_activation_markers_whole_umap.pdf", supp_fig_2_umaps, height=3.5, width=5.25)




#### vivax falciparum comparison ####

library(dplyr)
library(tidyr)
library(ggplot2)



lineage_palette <- c("#AA3377", "#EE6677", "#4477AA", "#66CCEE", "#228833", "#FFA500", "pink")
names(lineage_palette) <-c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")

#vac63c data

# first import some stuff that so that we don't have to keep defining everything, even if we want to just run a chunk
stacked_bar_data <-  read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_all_cluster_freqs.csv")
stacked_bar_data$volunteer <- factor(stacked_bar_data$volunteer, levels=c("v313", "v315", "v320", "v306", "v301", "v308", "v305", "v304", "v310"))
stacked_bar_data$cluster_id <- gsub("activated CD8lo Effector", "activated DN", stacked_bar_data$cluster_id)

# add column for lineage, relying on cluster_id names
stacked_bar_data$lineage <- substr(stacked_bar_data$cluster_id, nchar(stacked_bar_data$cluster_id)-2, nchar(stacked_bar_data$cluster_id))
#[1] "CD4" " gd" "AIT" "NKT" "CD8" " DN" "reg"

lin_replacement <- setNames(c("CD4", "gd", "MAIT", "NKT", "CD8", "DN", "Treg", "gd"), unique(stacked_bar_data$lineage))
stacked_bar_data$lineage <- stringr::str_replace_all(stacked_bar_data$lineage, lin_replacement)
stacked_bar_data$lineage <- factor(stacked_bar_data$lineage, levels=rev(c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")))

vac63c_data <- subset(stacked_bar_data, stacked_bar_data$volunteer %in% c("v313", "v315", "v320"))
vac63c_data$species <- "P. falciparum"

vac69a_data <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_cluster_freqs.csv", header = T)[,-1]

vac69a_data$species <- "P. vivax"

vac69a_data$lineage <- ifelse(grepl("CD4", vac69a_data$cluster_id), "CD4", NA)
vac69a_data$lineage <- ifelse(grepl("CD8", vac69a_data$cluster_id), "CD8", vac69a_data$lineage)
vac69a_data$lineage <- ifelse(grepl("Treg", vac69a_data$cluster_id), "Treg", vac69a_data$lineage)
vac69a_data$lineage <- ifelse(grepl("MAIT", vac69a_data$cluster_id), "MAIT", vac69a_data$lineage)
vac69a_data$lineage <- ifelse(grepl("gamma delta", vac69a_data$cluster_id), "gd", vac69a_data$lineage)
vac69a_data$lineage <- ifelse(grepl("DN", vac69a_data$cluster_id), "DN", vac69a_data$lineage)
vac69a_data$lineage <- ifelse(grepl("DP", vac69a_data$cluster_id), "DP", vac69a_data$lineage)

vac69a_data$volunteer <- gsub("V", "v", vac69a_data$volunteer, fixed=T)
vac69a_data$timepoint <- gsub("DoD", "Diagnosis", vac69a_data$timepoint, fixed=T)

vac69a_data$cluster_id <- gsub("naive", "Naive", vac69a_data$cluster_id, fixed=T)

all_data <- rbind(vac69a_data, vac63c_data)

# all cd3 activation ####

activated_vivax <- subset(vac69a_data, grepl("activated", vac69a_data$cluster_id))
activated_vivax$species <- "P. vivax"

activated_falciparum <- subset(vac63c_data, grepl("activated", vac63c_data$cluster_id))
activated_falciparum$species <- "P. falciparum"

all_activated <- rbind(activated_falciparum, activated_vivax)

all_activated <- filter(all_activated, timepoint=="T6")


activation_stacked_barchart <- ggplot(all_activated, aes(x=volunteer, y=frequency/100, fill=lineage))+
  geom_bar(stat="identity", position="stack")+
  #geom_text(aes(y=(total_cd3/100)+0.01, label=paste0(round(total_cd3, digits = 1), "%", sep='')))+
  theme_minimal()+
  #scale_x_discrete(drop=T)+
  facet_grid(timepoint~species, drop = T, scales = "free_x", space = "free")+
  ggtitle("Overall T cell Activation")+
  scale_fill_manual(values=lineage_palette)+#, labels=c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "Resting"))+
  scale_y_continuous(name = "Fraction of T cells activated", labels=scales::percent_format(accuracy = 1))+
  guides(fill=guide_legend(title="Lineage"))+
  theme(
    plot.title = element_text(hjust=0.5, size=12),
    strip.text = element_text(),
    strip.background = element_blank(),
    strip.text.x = element_text(face = "italic"), 
    axis.title.x = element_blank(),
    axis.text.x = element_text(hjust=0.5, angle=45),
    panel.spacing.x = unit(0.8,"lines"),
    axis.title.y = element_text(size=10),
    panel.grid.minor.y = element_blank())

# lineage_leg <- cowplot::get_legend(activation_stacked_barchart)
# 
# activation_stacked_barchart <- activation_stacked_barchart+theme(legend.position = "none")

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/vivax_falciparum_all_activation.pdf", activation_stacked_barchart, height = 5, width=5.5)





# non naive cd4 t cell activation ####

cd4_memory <- all_data %>%
  filter(timepoint=="T6", lineage=="CD4")

cd4_memory <- subset(cd4_memory, !grepl("Naive", cd4_memory$cluster_id))

cd4_memory_size <- cd4_memory %>%
  group_by(volunteer) %>%
  mutate("CD4_memory_percentage"=sum(frequency)) %>%
  select(volunteer, CD4_memory_percentage)%>%
  ungroup()

cd4_memory_size <- cd4_memory_size[!duplicated(cd4_memory_size), ]

cd4_bar_data <- cd4_memory %>%
  group_by(volunteer) %>%
  mutate("cd4_freq" = frequency/cd4_memory_size$CD4_memory_percentage[match(volunteer, cd4_memory_size$volunteer)])

cd4_bar_data <- subset(cd4_bar_data, grepl("activated",cd4_bar_data$cluster_id))




cd4_acti <- cd4_bar_data %>%
  group_by(volunteer, species) %>%
  summarise(sum=sum(cd4_freq))

wilcox.test(cd4_acti$sum~cd4_acti$species)
#p-value = 0.04762

(cd4_bar_plot <- ggplot(cd4_bar_data, aes(x=volunteer, y=cd4_freq, fill=lineage))+
    geom_bar(stat="identity", position="stack")+
    #geom_text(aes(y=(total_cd3/100)+0.01, label=paste0(round(total_cd3, digits = 1), "%", sep='')))+
    theme_minimal()+
    #scale_x_discrete(drop=T)+
    facet_grid(timepoint~species, drop = T, scales = "free_x", space = "free")+
    ggtitle("Non-Naive CD4 T cell Activation")+
    scale_fill_manual(values=lineage_palette, labels=c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "Resting"))+
    scale_y_continuous(name = "Fraction of Memory CD4 T cells Activated", labels=scales::percent_format(accuracy = 1))+
    theme(
      plot.title = element_text(hjust=0.5, size=12),
      strip.text = element_text(),
      strip.background = element_blank(),
      strip.text.x = element_text(face = "italic"), 
      axis.title.x = element_blank(),
      axis.text.x = element_text(hjust=0.5, angle=45),
      panel.spacing.x = unit(0.8,"lines"),
      axis.title.y = element_text(size=10),
      panel.grid.minor.y = element_blank(),
      legend.position = "none"))

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/vivax_falciparum_cd4_memory_activation.pdf", cd4_bar_plot, height = 5, width=5)





# Treg activation ####


all_data <- all_data %>%
  group_by(lineage, volunteer, timepoint) %>%
  mutate("sum"=sum(frequency)) %>%
  mutate("scaled_freq"=frequency/sum)


all_data$pie_fill <- ifelse(grepl("*activated*", all_data$cluster_id), as.character(all_data$lineage), "resting")

treg_memory <- all_data %>%
  filter(timepoint=="T6", lineage=="Treg")

treg_acti <- treg_memory %>%
  filter(pie_fill=="Treg") %>%
  group_by(volunteer, species, pie_fill) %>%
  summarise(sum=sum(scaled_freq))

wilcox.test(treg_acti$sum~treg_acti$species)



phil_palette <- c(lineage_palette, "resting"="grey")

treg_acti_plot <- ggplot(treg_memory, aes(x=volunteer, y=scaled_freq, fill=factor(pie_fill, levels=rev(c("Treg", "resting")))))+
  geom_bar(stat="identity", position="stack")+
  theme_minimal()+
  facet_grid(timepoint~species, scales="free_x", space="free")+
  #coord_polar(theta = "y")+
  scale_fill_manual(values=phil_palette)+
  scale_y_continuous(name = "Percentage of Treg Activated", labels=scales::percent_format(accuracy = 1))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text.x = element_text(hjust=0.5, size=10, face = "italic"),
        axis.title.x = element_blank(),
        legend.position="bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid.minor.y = element_blank(),
        strip.placement = "outside")



ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/vivax_falciparum_treg_activation.pdf", treg_acti_plot, height = 5, width=5)
