library(tidyr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

# data preperation ####

setwd("~/PhD/plasma/vac69a/")


data3 <- read.csv("big_plasma_table.csv")

data3[,3:ncol(data3)] <- log2(data3[,3:ncol(data3)])




data3 <- data3 %>%
  mutate(Volunteer = gsub("00", "0", Volunteer)) %>%
  mutate(timepoint = gsub("C-1", "Baseline", timepoint)) %>%
  mutate(timepoint = gsub("+", "", timepoint, fixed = T))

data3$timepoint <- factor(data3$timepoint, levels=c("Baseline", "DoD", "T6", "C45"))
data3 <- arrange(data3, timepoint, Volunteer)

data3$Sample_ID <- paste(data3$Volunteer, "_", data3$timepoint, sep='')

#censor volutneer 03
data3 <- subset(data3, data3$Volunteer != "v03")



data4 <- data3
data4[, 1:2] <- NULL


trans_data <- data.frame("Sample_ID"=data4$Sample_ID, apply(data4[,1:ncol(data4)-1], 2, function(x) scale(x, center = TRUE, scale = TRUE)))


plasma_levels <- read.delim("analytes_sorted_by_padj.txt")


t_trans_data <- t(trans_data)
colnames(t_trans_data) <- t_trans_data[1,]

plasma_matrix <- t_trans_data[2:nrow(t_trans_data),]

class(plasma_matrix) <- "numeric"


plasma_matrix <- plasma_matrix[match(plasma_levels$Analyte, rownames(plasma_matrix)),]


plasma_matrix <- ifelse(plasma_matrix>2.5, 2.5, plasma_matrix)
plasma_matrix <- ifelse(plasma_matrix < -2.5, -2.5, plasma_matrix)
plasma_matrix <- plasma_matrix[1:22,]

# top anno ####

#define colors for timepoints
time_col <- colorspace::sequential_hcl(5, palette = "Purple Yellow")

Volunteer <- c("v02" = "#FB9A99","v03" = "#E31A1C","v05" = "#A6CEE3", "v06" = "#1F78B4", "v07" = "#B2DF8A")
Timepoint <- c("Baseline"=time_col[4], "DoD"=time_col[2], "T6"=time_col[1], "C45"=time_col[5])


combo_top_anno <- HeatmapAnnotation(gap = unit(2, "mm"), annotation_name_side = "left",
                                    Volunteer = data3$Volunteer,
                                    Timepoint = data3$timepoint,
                                    col=list(Timepoint = c("Baseline"=time_col[4], "DoD"=time_col[2], "T6"=time_col[1], "C45"=time_col[5]),
                                             Volunteer = c("v02" = "#FB9A99","v03" = "#E31A1C","v05" = "#A6CEE3", "v06" = "#1F78B4", "v07" = "#B2DF8A")),
                                    simple_anno_size = unit(2.5, "mm"),
                                    annotation_legend_param = list(
                                      Volunteer = list(title = "Volunteer", at = names(Volunteer), legend_gp = gpar(fill = unname(Volunteer)), title_position = "topleft"),
                                      Timepoint = list(title ="Timepoint",at = names(Timepoint), legend_gp = gpar(fill = unname(Timepoint)), title_position = "topleft")
                                    )
)


# left anno ####
number_of_hits <- 7

significant <-  c(rep("yes",number_of_hits), rep("no", nrow(plasma_matrix)-number_of_hits))
sig <- c("yes"="darkgreen", "no"="lightgrey")


left_anno <-  rowAnnotation(gap = unit(5, "mm"),
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




# right anno ####

#heatmap construction ####


col_fun4 <- circlize::colorRamp2(c(-2.5, 0, 2.5), c("#0859C6", "black", "#FFA500"))



combo_map <- Heatmap(matrix = plasma_matrix,
                     cluster_rows = FALSE,
                     name = "Normalised Plasma Concentration",
                     cluster_columns = FALSE,
                     row_names_side = "left",
                     col = col_fun4,
                     #row_split = factor(rep(c("up", "down"), each = 12), levels = c("up", "down")),
                     rect_gp = gpar(col = "white"),
                     #row_title = c("",""),
                     top_annotation = combo_top_anno,
                     #right_annotation = combo_right_anno,
                     left_annotation = left_anno,
                     show_heatmap_legend = TRUE,
                     column_names_rot = 45,
                     heatmap_legend_param = list(col = col_fun4, title = "Z-Score", title_position = "topleft"),
                     #width = unit(16, "cm"),
                     #height = unit(16, "cm")
)



png("./figures/plasma_zscore_heatmap_without_v03_fdr_10-e1.png", width=12, height=9, units = "in", res=400)
draw(combo_map,
     merge_legends = TRUE,
     #padding = unit(c(2, 20, 2, 2), "mm")
)
dev.off()










