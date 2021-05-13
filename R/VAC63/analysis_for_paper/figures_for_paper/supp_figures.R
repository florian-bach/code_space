# preamble ####

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


all_cells_acti_summary <- stacked_bar_data %>%
  group_by(lineage, volunteer, timepoint) %>%
  mutate("sum"=sum(frequency)) %>%
  mutate("scaled_freq"=frequency/sum)


all_cells_acti_summary$pie_fill <- ifelse(grepl("*activated*", all_cells_acti_summary$cluster_id), as.character(all_cells_acti_summary$lineage), "resting")




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






#gd mait activation stack ####


all_gd_data <- subset(stacked_bar_data, grepl("gd", stacked_bar_data$cluster_id, fixed=T))
all_mait_data <- subset(stacked_bar_data, grepl("MAIT", stacked_bar_data$cluster_id, fixed=T))


View(all_gd_data %>%
       group_by(volunteer, timepoint) %>%
       summarise(sum(frequency)))
# 
# gd_activation_stacked_barchart <- ggplot(all_mait_data, aes(x=volunteer, y=frequency/100, fill=factor(cluster_id)))+
#   geom_bar(stat="identity", position="stack")+
#   theme_minimal()+
#   facet_wrap(~timepoint, strip.position = "bottom", ncol=4)+
#   scale_fill_manual(values=col_pal)+
#   scale_y_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1))+
#   theme(plot.title = element_text(hjust=0.5, size=11),
#         strip.text = element_text(hjust=0.5, size=10, face = "bold"),
#         axis.title.x = element_blank(),
#         legend.position="bottom",
#         legend.direction = "horizontal",
#         legend.title = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust=1),
#         panel.grid.minor.y = element_blank(),
#         strip.placement = "outside")
# 
# ggsave("/home/flobuntu/PhD/figures_for_thesis/chapter_03/gd_freq_plot.pdf", gd_activation_stacked_barchart, height=4, width=8)
# 
# 

all_mamma_delta_data <- rbind(all_gd_data, all_mait_data)

all_mamma_delta_data_summary <- all_mamma_delta_data %>%
  group_by(lineage, volunteer, timepoint) %>%
  mutate("sum"=sum(frequency)) %>%
  mutate("scaled_freq"=frequency/sum)



activated_all_mamma_delta_data <- subset(all_mamma_delta_data, grepl("*activated*", all_mamma_delta_data$cluster_id))


activated_all_mamma_delta_plot <- ggplot(activated_all_mamma_delta_data, aes(x=volunteer, y=frequency/100, fill=factor(lineage)))+
  geom_bar(stat="identity", position="stack")+
  theme_minimal()+
  facet_wrap(~timepoint, strip.position = "bottom", ncol=4)+
  scale_fill_manual(values=lineage_palette)+
  scale_y_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        legend.position="bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid.minor.y = element_blank(),
        strip.placement = "outside")

ggsave("/home/flobuntu/PhD/figures_for_thesis/chapter_03/gd_mait_stack.pdf", activated_all_mamma_delta_plot, height=4, width=9)




all_mamma_delta_data_summary$pie_fill <- ifelse(grepl("*activated*", all_mamma_delta_data$cluster_id), as.character(all_mamma_delta_data$lineage), "resting")

phil_palette <- c(lineage_palette, "resting"="grey")

mait_gd_acti_plot <- ggplot(all_mamma_delta_data_summary, aes(x=volunteer, y=scaled_freq, fill=factor(pie_fill, levels=rev(c("gd", "MAIT", "resting")))))+
  geom_bar(stat="identity", position="stack")+
  theme_minimal()+
  facet_grid(lineage~timepoint)+
  #coord_polar(theta = "y")+
  scale_fill_manual(values=phil_palette)+
  scale_y_continuous(name = "Percentage of Lineage Activated", labels=scales::percent_format(accuracy = 1))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        legend.position="bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid.minor.y = element_blank(),
        strip.placement = "outside")


#ggsave("/home/flobuntu/PhD/figures_for_thesis/chapter_03/mait_gd_acti_plot.pdf", mait_gd_acti_plot, height=8, width=9)
ggsave("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/figures_for_paper/mait_gd_acti_plot.pdf", mait_gd_acti_plot, height=8, width=9)






# ds_limma   ####
# 

num_wide_scaled_ms <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/sample_wise_cluster_expression_matrix.csv", header=T, row.names = 1)

ter_dod_ms <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/ds_limma/ter_dod_ms.csv", header=T, row.names = 1)
prim_t6_ms <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/ds_limma/prim_t6_ms.csv", header=T, row.names = 1)
ter_t6_ms <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/ds_limma/ter_t6_ms.csv", header=T, row.names = 1)
# 
#
# slot ds limma dataset in here! 
num_wide_scaled_ms <- prim_t6_ms

#label each df
prim_t6_ms$n_infection <- "first"
ter_t6_ms$n_infection <- "third"

#add contrast as column
prim_t6_ms$rowname <- rownames(prim_t6_ms)
ter_t6_ms$rowname <- rownames(ter_t6_ms)

#combine
combo_t6_ms <- rbind(prim_t6_ms, ter_t6_ms)

#identify shared features
combo_t6_ms$n_infection <- ifelse(duplicated(combo_t6_ms$rowname), "both", combo_t6_ms$n_infection)
boths <- combo_t6_ms$rowname[combo_t6_ms$n_infection=="both"]

# remove duplicates
num_wide_scaled_ms <- combo_t6_ms[!(combo_t6_ms$rowname %in% boths & combo_t6_ms$n_infection !="both"),-38]

limma_n_infection <- num_wide_scaled_ms$n_infection


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


combo_left_anno_var <-  rowAnnotation(annotation_name_gp = gpar(fontsize=10),
                                      annotation_name_rot = 45,
                                      gap = unit(1.5, "mm"),
                                      "Significant In"=limma_n_infection,
                                      show_legend = FALSE,
                                      show_annotation_name = TRUE,
                                      simple_anno_size = unit(4, "mm"), # width of the significance bar
                                      col=list("Significant In" = c("first"="#36454f","third"="darkgrey", "both"="#c227ba")),
                                      annotation_legend_param = list("Significant In" = list(title ="Significant In",
                                                                                        at = c("first", "both", "third"),
                                                                                        title_gp=gpar(angle=45),
                                                                                        legend_gp = gpar(fill = c("#36454f", "#darkgrey", "#c227ba")),
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

limma_leg3 = Legend(at = c("first", "both", "third"),
                    type = "grid",
                    legend_gp = gpar(fill = c("#36454f", "#c227ba", "darkgrey")),
                    title_position = "topleft",
                    direction = "horizontal",
                    nrow = 1,
                    labels_gp = gpar(fontsize =9),
                    title = "Significant In")

rownames(num_wide_scaled_ms) <- gsub(")1", ")", rownames(num_wide_scaled_ms))

median_cluster_heat <- Heatmap(matrix = num_wide_scaled_ms,
                               cluster_rows = TRUE,
                               column_order = 1:36,
                               name = "Scaled Marker Expression",
                               cluster_columns = FALSE,
                               show_row_dend = TRUE,
                               row_dend_side = "right",
                               row_names_side = "left",
                               col = col_fun_ds_limma,
                               left_annotation = combo_left_anno_var,
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


pdf("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/figures_for_paper/combo_prim_ter_t6_limma.pdf", width = 10, height=12)
draw(median_cluster_heat,
     annotation_legend_list = list(limma_leg, limma_leg2, limma_leg3),
     merge_legends = TRUE, heatmap_legend_side = "bottom",)
dev.off()
