#Panel A: Heatmap of Marker expression of significant clusters only ####



reordered_sig_scaled_mat <- as.matrix(read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/sig_cluster_medians_heatmap_t6.csv", header = T, row.names = 1))

sig_clusters <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/sig_t6_clusters.csv", header = TRUE, stringsAsFactors = FALSE, col.names = F)
sig_clusters <- sig_clusters[-9,]


# Right annotation ###

all_t6_data <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_t6_data.csv", stringsAsFactors = FALSE, header = T)

t6_map_data <- all_t6_data%>%
  group_by(cluster_id, timepoint) %>%
  mutate(mean_freq=mean(frequency)) %>%
  ungroup() %>%
  select(cluster_id, mean_freq, timepoint)

t6_map_data <- filter(t6_map_data, timepoint=="T6")
t6_map_data <- t6_map_data[!duplicated(t6_map_data), ]
t6_map_data <- filter(t6_map_data, cluster_id %in% sig_clusters)
t6_map_data <- t6_map_data[order(t6_map_data$mean_freq, decreasing = T),]

rereordered_sig_scaled_mat <- reordered_sig_scaled_mat[match(t6_map_data$cluster_id, rownames(reordered_sig_scaled_mat)),]



cd3_right_anno <- rowAnnotation(gap = unit(2, "mm"),
                                "Mean Frequency at T6" = anno_barplot(t6_map_data$mean_freq, which="row", axis = TRUE, ylim = c(0, 4)),
                                width = unit(4, "cm")
)



colcsv <- read.csv("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cluster_palette.csv", header=T, stringsAsFactors = F)
col_pal <- colcsv$x
names(col_pal) <- colcsv$X

colz <- unname(col_pal[match(rownames(rereordered_sig_scaled_mat), names(col_pal))])
breakz <- names(col_pal[match(rownames(rereordered_sig_scaled_mat), names(col_pal))])

cd3_right_anno_var <- rowAnnotation(gap = unit(2, "mm"),
                                    # "Mean Frequency at T6" = anno_boxplot(slim_wide_t6_map_data, which="row", axis = TRUE, gp=gpar(fill=c("#228833","#AA3377", "#AA3377", "#66CCEE", "#AA3377", "#AA3377", "#AA3377", "#4477AA", "#4477AA" ))),
                                    "Mean Percentage\nof CD3+ at T6" = anno_barplot(t6_map_data$mean_freq, which="row", axis = TRUE, ylim = c(0, 4), gp=gpar(fill=colz)),
                                    
                                    width = unit(4, "cm"),
                                    annotation_legend_param = list("Mean Percentage\nof CD3+ at T6" = list(title ="Lineage",
                                                                                                           at = breakz),
                                                                   legend_gp = gpar(fill = colz),
                                                                   title_position = "topleft")
)




draw(cd3_right_anno_var)


box_lgd <- Legend(labels =  breakz,
                  title = "Cluster_ID",
                  type = "grid",
                  legend_gp = gpar(fill = colz)
)


inferno <- colorspace::sequential_hcl("inferno", n=8)
col_inferno <- circlize::colorRamp2(seq(0,1, by=1/{length(inferno)-1}), inferno)



median_cluster_heat <- Heatmap(matrix = rereordered_sig_scaled_mat,
                               cluster_rows = FALSE,
                               name = "Median Marker Expression",
                               cluster_columns = FALSE,
                               row_names_side = "left",
                               col = col_inferno,
                               rect_gp = gpar(col = "white"),
                               #top_annotation = combo_top_anno,
                               right_annotation = cd3_right_anno_var,
                               show_heatmap_legend = FALSE,
                               column_names_rot = 60,
                               #heatmap_legend_param = list(col = col_fun4, title = "Normalised Frequency", title_position = "topleft"),
                               width = unit(16, "cm"),
                               height = unit(16*9/28, "cm")
)


png("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/diffcyt/edgeR/sig_cluster_t6_phenotype_heat_var.png", width=13, height=4, units = "in", res=400)
draw(median_cluster_heat,
     #annotation_legend_list = list(box_lgd),
     merge_legends = FALSE,
     #padding = unit(c(200, 200, 200, 200), "mm")
)
dev.off() 


#Panel B: Pie Chart of significant cluster size and lineage ####




#Panel C: stacked Piechart of activated CD4 cluster sizes relative to size of memroy pool ####
#Panel D: Show marker expression on T6 CD4 bump ####
