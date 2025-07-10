library(colorspace)
library(scales)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
#get experiment info
# batch <- rep(rep(levels(cd$batch), each=3), 4)
# volunteer <- rep(levels(cd$volunteer), 4)
# timepoint <- rep(levels(cd$timepoint), each=6)

#names(volunteer) <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")

#define colors for timepoints
time_col <- colorspace::sequential_hcl(5, palette = "Purple Yellow")


#cluster id and frequencies
all_t6_data <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_t6_data.csv")
#meta data
cd <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/experiment_info.csv")
#fold change and p_adj
fold_change <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/edger_t6_FC_padj.csv")

#order decreasing so top10 and bottom 10 are extreme in terms of fold change
fold_change <- fold_change[order(fold_change$logFC, decreasing = T),]

#get the correct order of cluster ids by fold change
cluster_levels <- fold_change$cluster_id


#calculate z score after asinh transform of square
t6_map_data <- all_t6_data%>%
  group_by(cluster_id) %>%
  #group_by(volunteer) %>%
  mutate(trans_freq=asin(sqrt(frequency/100))) %>%
  mutate(max_freq=max(frequency)) %>%
  mutate(trans_norm_freq=scale(trans_freq, center = TRUE, scale = TRUE)) %>%
  ungroup()



# construct matrix to be plotted in heatmap: clusters along rows, samples along columns, values=trans_norm_freq
hm_matrix <- as.data.frame(t6_map_data %>%
                             select(cluster_id, sample_id, trans_norm_freq) %>%
                             tidyr::pivot_wider(names_from = sample_id, values_from = trans_norm_freq))

rownames(hm_matrix) <- hm_matrix$cluster_id
# this orders the matrix according to fold change, as cluster_levels is..
hm_matrix <- hm_matrix[match(cluster_levels, hm_matrix$cluster_id),]

#matrix has to be numeric only
hm_matrix <- as.matrix(select(hm_matrix, -cluster_id))


#truncate matrix values so that extreme values don't pale everything elese
hm_matrix <- ifelse(hm_matrix>2.5, 2.5, hm_matrix)
hm_matrix <- ifelse(hm_matrix < -2.5, -2.5, hm_matrix)
extreme <- max(abs(range(hm_matrix)))

#bunch of differen color gradients
    col_fun4 <- circlize::colorRamp2(c(-range(hm_matrix)[2], 0, range(hm_matrix)[2]), c("#0859C6", "black", "#FFA500"))

#get 12 most up/down clusters
top_mat <- head(hm_matrix, n=12)
bot_mat <- tail(hm_matrix, n=12)
combo_matrix <- rbind(top_mat, bot_mat)

# make seperate array for p_adj, significantce, fold change etc for heatmap annotations
# right annotation
p_adj <- c(head(fold_change$p_adj, n=nrow(top_mat)), tail(fold_change$p_adj, n=nrow(bot_mat)))
log2_fc <- c(head(fold_change$logFC, n=nrow(top_mat)), tail(fold_change$logFC, n=nrow(bot_mat)))
significant <- ifelse((p_adj<0.05 & abs(log2_fc)>1), "yes", "no")
sig <- c("yes"="darkgreen", "no"="lightgrey")

#top annotation
Volunteer <- c("V02" = "#FB9A99","V03" = "#E31A1C","V05" = "#A6CEE3", "V06" = "#1F78B4", "V07" = "#B2DF8A", "V09" = "#33A02C")
Timepoint <- c("Baseline"=time_col[4], "C10"=time_col[3], "DoD"=time_col[2], "T6"=time_col[1])
Batch <- c("batch_1"="lightgrey", "batch_2"="darkgrey")



combo_right_anno <-  rowAnnotation(gap = unit(2, "mm"),
                                  #space = rep("a", 24),
                                  #"log2FC" = anno_lines(log2_fc, which="row", axis = TRUE, ylim = c(-6, 6), axis_param = list(at=seq(-6, 6, by=2)), smooth = FALSE, add_points = TRUE),
                                  #"significant"=significant,
                                  "log2FC" = anno_barplot(log2_fc, which="row", axis = TRUE, ylim = c(-2, 6), axis_param = list(at=seq(-2, 6, by=2)),gp = gpar(fill = col_fun4(log2_fc))),
                                  
                                  #"p_adj"= anno_text(scales::scientific(p_adj, digits = 2), which="row"),
                                  width = unit(3, "cm") # width of the line graph
                                  #simple_anno_size = unit(2, "mm"), # width of the significance bar
                                  )
                                  




# combo_right_anno_var <-  rowAnnotation(gap = unit(5, "mm"),
#                                    #space = rep("a", 24),
#                                    #"log2FC" = anno_lines(log2_fc, which="row", axis = TRUE, ylim = c(-6, 6), axis_param = list(at=seq(-6, 6, by=2)), smooth = FALSE, add_points = TRUE),
#                                    "significant"=significant,
#                                    #"log2FC" = anno_lines(log2_fc, which="row", axis = TRUE, ylim = c(-4, 6), axis_param = list(at=seq(-4, 6, by=2)), add_points = TRUE, smooth = FALSE),
#                                    
#                                    #"p_adj"= anno_text(scales::scientific(p_adj, digits = 2), which="row"),
#                                    width = unit(3, "cm"), # width of the line graph
#                                    simple_anno_size = unit(2, "mm"), # width of the significance bar
#                                    col=list(significant = c("yes"="darkgreen", "no"="lightgrey"),
#                                             space=c("a"="white")),
#                                    annotation_legend_param = list(significant = list(title ="Significant",at = rev(names(sig)), legend_gp = gpar(fill = unname(sig)), title_position = "topleft")
#                                    )
#                                    
# )



combo_left_anno_var <-  rowAnnotation(gap = unit(5, "mm"),
                                      #annotation_name_gp = gpar(angle=45),
                                      show_annotation_name = FALSE,
                                      "significant"=significant,
                                       simple_anno_size = unit(2.5, "mm"), # width of the significance bar
                                       col=list("significant" = c("yes"="darkgreen", "no"="lightgrey")),
                                       annotation_legend_param = list(significant = list(title ="Significant",
                                                                                         at = rev(names(sig)),
                                                                                         #title_gp=gpar(angle=45),
                                                                                         legend_gp = gpar(fill = unname(sig)),
                                                                                         title_position = "topleft")
                                       )
                                       
)



combo_top_anno <- HeatmapAnnotation(gap = unit(2, "mm"), annotation_name_side = "left",
  Volunteer = rep(levels(cd$volunteer), 4),
  Timepoint = rep(levels(cd$timepoint), each=6),
  Batch = rep(rep(levels(cd$batch), each=3), 4),
  col=list(Batch = c("batch_1"="lightgrey", "batch_2"="darkgrey"),
           Timepoint = c("Baseline"=time_col[4], "C10"=time_col[3], "DoD"=time_col[2], "T6"=time_col[1]),
           Volunteer = c("V02" = "#FB9A99","V03" = "#E31A1C","V05" = "#A6CEE3", "V06" = "#1F78B4", "V07" = "#B2DF8A", "V09" = "#33A02C")
           
  ),
  simple_anno_size = unit(2.5, "mm"),
  annotation_legend_param = list(
  Volunteer = list(title = "Volunteer", at = names(Volunteer), legend_gp = gpar(fill = unname(Volunteer)), title_position = "topleft"),
  Timepoint = list(title ="Timepoint",at = names(Timepoint), legend_gp = gpar(fill = unname(Timepoint)), title_position = "topleft"),
  Batch = list(title = "Batch", at = names(Batch), legend_gp = gpar(fill = unname(Batch)), title_position = "topleft")
  )
)


#ht_opt("simple_anno_size" = unit(2, "mm"))

combo_map <- Heatmap(matrix = combo_matrix,
        cluster_rows = FALSE,
        name = "Normalised Frequency",
        cluster_columns = FALSE,
        row_names_side = "left",
        col = col_fun4,
        row_split = factor(rep(c("up", "down"), each = 12), levels = c("up", "down")),
        rect_gp = gpar(col = "white"),
        row_title = c("",""),
        top_annotation = combo_top_anno,
        right_annotation = combo_right_anno,
        left_annotation = combo_left_anno_var,
        show_heatmap_legend = TRUE,
        column_names_rot = 45,
        heatmap_legend_param = list(col = col_fun4, title = "Normalised Frequency", title_position = "topleft"),
        width = unit(16, "cm"),
        height = unit(16, "cm")
)

draw(combo_map,
     merge_legends = TRUE,
     #padding = unit(c(2, 20, 2, 2), "mm")
     )


  png("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/diffcyt/edgeR/improved_freq_hm_no_right.png", width=14, height=10, units = "in", res=400)
draw(combo_map,
     merge_legends = TRUE,
     #padding = unit(c(2, 20, 2, 2), "mm")
)
dev.off()



# old designs you'll proabbaly not go back too ####



# manual control of legens



# freq_lgd <- Legend(col = col_fun4, title = "Normalised Frequency", title_position = "topleft")
# 
# vol_lgd <- Legend(title = "Volunteer", at = names(Volunteer), legend_gp = gpar(fill = unname(Volunteer)), title_position = "topleft")
# 
# time_lgd <- Legend(title ="Timepoint",at = names(Timepoint), legend_gp = gpar(fill = unname(Timepoint)), title_position = "topleft")
# 
# batch_lgd <- Legend(title ="Batch",at = names(Batch), legend_gp = gpar(fill = unname(Batch)), title_position = "topleft")
# 
# sig_lgd <- Legend(title ="Significant",at = names(sig), legend_gp = gpar(fill = unname(sig)), title_position = "topleft")

#make heatmap annotation object hard coding each color for each property
top_anno <- HeatmapAnnotation(
  Volunteer = rep(levels(cd$volunteer), 4),
  Timepoint = rep(levels(cd$timepoint), each=6),
  Batch = rep(rep(levels(cd$batch), each=3), 4),
  col=list(Batch = c("batch_1"="lightgrey", "batch_2"="darkgrey"),
           Timepoint = c("Baseline"=time_col[4], "C10"=time_col[3], "DoD"=time_col[2], "T6"=time_col[1]),
           Volunteer = c("V02" = "#FB9A99","V03" = "#E31A1C","V05" = "#A6CEE3", "V06" = "#1F78B4", "V07" = "#B2DF8A", "V09" = "#33A02C")
           
  ),simple_anno_size = unit(3, "mm")
)





freq_hm <- Heatmap(matrix = hm_matrix,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   row_names_side = "left",
                   col = col_fun4,
                   rect_gp = gpar(col = "white"),
                   top_annotation = top_anno,
                   show_heatmap_legend = FALSE,
                   width = unit(12, "cm"),
                   height = unit(17, "cm")
)
freq_hm




top_mat <- head(hm_matrix, n=10)

top_line_anno <-  rowAnnotation("log2FC" = anno_lines(head(fold_change$logFC, n=nrow(top_mat)), which="row", axis = TRUE, smooth = TRUE, add_points = TRUE))

bot_line_anno <-  rowAnnotation("log2FC" = anno_lines(tail(fold_change$logFC, n=nrow(bot_mat)), which="row", axis = TRUE, smooth = TRUE, add_points = TRUE))



top_hm <- Heatmap(matrix = top_mat,
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  row_names_side = "left",
                  col = col_fun4,
                  rect_gp = gpar(col = "white"),
                  top_annotation = top_anno,
                  right_annotation = top_line_anno,
                  show_heatmap_legend = FALSE,
                  width = unit(12, "cm"),
                  height = unit(6, "cm")
)



bot_mat <- tail(hm_matrix, n=10)

bot_hm <- Heatmap(matrix = bot_mat,
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  row_names_side = "left",
                  col = col_fun4,
                  rect_gp = gpar(col = "white"),
                  right_annotation = bot_line_anno,
                  show_heatmap_legend = FALSE,
                  width = unit(12, "cm"),
                  height = unit(6, "cm")
)


draw(top_hm %v% bot_hm)

