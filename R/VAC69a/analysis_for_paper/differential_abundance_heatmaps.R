library(colorspace)
library(scales)
library(ComplexHeatmap)

#get experiment info
cd <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/experiment_info.csv")
# batch <- rep(rep(levels(cd$batch), each=3), 4)
# volunteer <- rep(levels(cd$volunteer), 4)
# timepoint <- rep(levels(cd$timepoint), each=6)

#names(volunteer) <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")

#define colors for timepoints
time_col <- colorspace::sequential_hcl(5, palette = "Purple Yellow")

#make heatmap annotation object hard coding each color for each property
top_anno <- HeatmapAnnotation(
  Batch = rep(rep(levels(cd$batch), each=3), 4),
  Volunteer = rep(levels(cd$volunteer), 4),
  Timepoint = rep(levels(cd$timepoint), each=6),
  col=list(Batch = c("batch_1"="lightgrey", "batch_2"="darkgrey"),
           Timepoint = c("Baseline"=time_col[4], "C10"=time_col[3], "DoD"=time_col[2], "T6"=time_col[1]),
           Volunteer = c("V02" = "#FB9A99","V03" = "#E31A1C","V05" = "#A6CEE3", "V06" = "#1F78B4", "V07" = "#B2DF8A", "V09" = "#33A02C")
           
  )
)


#get count data, p values, all that stuff from diccyt
t6_map_data <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/t6_map_data.csv", header = TRUE, stringsAsFactors = TRUE)

#make the matrix object that will be the input for ComplexHeatmap::Heatmap
hm_matrix <- as.data.frame(t6_map_data %>%
                             select(cluster_id, sample_id, trans_norm_freq) %>%
                             tidyr::pivot_wider(names_from = sample_id, values_from = trans_norm_freq))

rownames(hm_matrix) <- hm_matrix$cluster_id
hm_matrix <- as.matrix(select(hm_matrix, -cluster_id))




#truncate matrix values so that extreme values don't pale everything elese
hm_matrix <- ifelse(hm_matrix>2.5, 2.5, hm_matrix)
hm_matrix <- ifelse(hm_matrix < -2.5, -2.5, hm_matrix)
extreme <- max(abs(range(hm_matrix)))

#bunch of differen color gradients
col_fun <- circlize::colorRamp2(c(-extreme, 0, extreme), c("#002940", "white", "#FFA500"))
col_fun2 <- circlize::colorRamp2(c(-extreme, 0, extreme), c("#0C71E0", "white", "#FFA500"))
col_fun3 <- circlize::colorRamp2(c(range(hm_matrix)[1], 0, range(hm_matrix)[2]), c("#0859C6", "white", "#FFA500"))
col_fun4 <- circlize::colorRamp2(c(-range(hm_matrix)[2], 0, range(hm_matrix)[2]), c("#0859C6", "black", "#FFA500"))

freq_hm <- Heatmap(matrix = hm_matrix,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   row_names_side = "left",
                   col = col_fun4,
                   rect_gp = gpar(col = "white"),
                   top_annotation = top_anno,
                   show_heatmap_legend = FALSE,
                   width = unit(14, "cm"),
                   height = unit(8, "cm")
)
freq_hm




#save as pdf
pdf("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/diffcyt/edgeR/freq_hm.pdf", width=11, height=6)
freq_hm
dev.off()

#add raster arguments in order to save as png
freq_hm <- Heatmap(matrix = hm_matrix,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   row_names_side = "left",
                   col = col_fun4,
                   rect_gp = gpar(col = "white"),
                   top_annotation = top_anno,
                   show_heatmap_legend = TRUE,
                   width = unit(14, "cm"),
                   height = unit(8, "cm"),
                   use_raster = TRUE,
                   raster_device = "png")

#save as png
png("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/diffcyt/edgeR/freq_hm.png", width=11, height=6, units = "in", res=400)
freq_hm
dev.off()





# some unfinished legend business to make the legends horizontal
volunteer_col <- c("V02" = "#FB9A99","V03" = "#E31A1C","V05" = "#A6CEE3", "V06" = "#1F78B4", "V07" = "#B2DF8A", "V09" = "#33A02C")
batch_col <- c("batch_1"="lightgrey", "batch_2"="darkgrey")
timepoint_col = c("Baseline"=time_col[1], "C10"=time_col[2], "DoD"=time_col[3], "T6"=time_col[4])


volunteer_lgd = Legend(at = names(volunteer_col), legend_gp = gpar(fill=volunteer_col), title = "Volunteer",
                       nr = 1)

batch_lgd = Legend(at = names(batch_col), legend_gp = gpar(fill=batch_col), title = "Batch",
                   nr = 1)

timepoint_lgd = Legend(at = names(timepoint_col), legend_gp = gpar(fill=timepoint_col), title = "Timepoint",
                       nr = 1)




pd <- packLegend(volunteer_lgd, batch_lgd, timepoint_lgd, direction = "horizontal")

