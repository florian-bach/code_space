library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")
sig_base_zero_infectiontype <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/differential_abundance/sig_base_zero_infectiontype.csv")

col_fun_pearson <- circlize::colorRamp2(c(-1, 0, 1), c("#0859C6", "black", "#FFA500"))

# correlation heatmaps ####

##  sympomtatic day 0####


S_day0 <- clean_data %>%
  filter(timepoint %in% c("day0"), infectiontype=="S")%>%
  filter(targetName %in% sig_base_zero_infectiontype$targetName)%>%
  filter(!is.na(parasitedensity))%>%
  pivot_wider(names_from = targetName, values_from = concentration, id_cols = c(id, parasitedensity))

S_day0_cor <- cor(S_day0[,2:ncol(S_day0)], method = "spearman")

day0_dist <- dist(S_day0_cor, method = "euclidean", diag = FALSE, upper = FALSE)
day0_hclust <- hclust(day0_dist)

# S_day0_cor_matrix <- as.matrix(S_day0_cor)
S_day0_cor_matrix <- as.matrix(S_day0_cor[rownames(S_day0_cor)[rev(day0_hclust$order)],colnames(S_day0_cor)[rev(day0_hclust$order)]])



S_day0_cor_heatmap <- Heatmap(matrix = S_day0_cor_matrix,
                           cluster_rows = FALSE,
                           cluster_columns=FALSE,
                           # show_row_dend = FALSE,
                           # show_column_dend = TRUE,
                           show_heatmap_legend = TRUE,
                           name = "Symptomatic Day 0",
                           #cluster_columns = FALSE,
                           column_names_gp = gpar(fontsize = 6),
                           row_names_gp = gpar(fontsize = 6),
                           row_names_side = "left",
                           col = col_fun_pearson,
                           column_names_rot = 90)

png("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/s_day0_big_heat.png", height=10, width = 10, units = "in", res=444)
draw(S_day0_cor_heatmap, padding=unit(c(2,2,2,2), "mm"))
dev.off()



##  asymptomatic day0 ####

A_day0 <- clean_data %>%
  filter(timepoint %in% c("day0"), infectiontype=="A")%>%
  filter(targetName %in% sig_base_zero_infectiontype$targetName)%>%
  filter(!is.na(parasitedensity))%>%
  pivot_wider(names_from = targetName, values_from = concentration, id_cols = c(id, parasitedensity))

A_day0_cor <- cor(A_day0[,2:ncol(A_day0)], method = "spearman")

# A_day0_cor_matrix <- as.matrix(A_day0_cor)
  
A_day0_cor_matrix <- as.matrix(A_day0_cor[rownames(A_day0_cor)[rev(day0_hclust$order)],colnames(A_day0_cor)[rev(day0_hclust$order)]])


col_fun_pearson <- circlize::colorRamp2(c(min(A_day0_cor_matrix), 0, max(A_day0_cor_matrix)), c("#0859C6", "black", "#FFA500"))

A_day0_cor_heatmap <- Heatmap(matrix = A_day0_cor_matrix,
                              cluster_rows = FALSE,
                              cluster_columns=FALSE,
                              show_row_dend = FALSE,
                              # show_column_dend = TRUE,
                              # show_heatmap_legend = TRUE,
                              name = "Asymptomatic Day 0",
                              #cluster_columns = FALSE,
                              column_names_gp = gpar(fontsize = 6),
                              row_names_gp = gpar(fontsize = 6),
                              row_names_side = "left",
                              col = col_fun_pearson,
                              column_names_rot = 90)

png("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/A_day0_big_heat.png", height=10, width = 10, units = "in", res=444)
draw(A_day0_cor_heatmap, padding=unit(c(2,2,2,2), "mm"))
dev.off()





##  sympomtatic baseline ####



S_baseline <- clean_data %>%
  filter(timepoint %in% c("baseline"), infectiontype=="S")%>%
  filter(targetName %in% sig_base_zero_infectiontype$targetName)%>%
  filter(!is.na(parasitedensity))%>%
  pivot_wider(names_from = targetName, values_from = concentration, id_cols = c(id, parasitedensity))

S_baseline_cor <- cor(S_baseline[,2:ncol(S_baseline)], method = "spearman")

# S_baseline_cor_matrix <- as.matrix(S_baseline_cor)

baseline_dist <- dist(S_baseline_cor, method = "euclidean", diag = FALSE, upper = FALSE)
baseline_hclust <- hclust(baseline_dist)

S_baseline_cor_matrix <- as.matrix(S_baseline_cor[rownames(S_baseline_cor)[rev(baseline_hclust$order)],colnames(S_baseline_cor)[rev(baseline_hclust$order)]])


col_fun_pearson <- circlize::colorRamp2(c(min(S_baseline_cor_matrix), 0, max(S_baseline_cor_matrix)), c("#0859C6", "black", "#FFA500"))

S_baseline_cor_heatmap <- Heatmap(matrix = S_baseline_cor_matrix,
                              cluster_rows = FALSE,
                              cluster_columns=FALSE,
                              # show_row_dend = FALSE,
                              # show_column_dend = TRUE,
                              show_heatmap_legend = TRUE,
                              name = "Symptomatic Baseline",
                              #cluster_columns = FALSE,
                              column_names_gp = gpar(fontsize = 6),
                              row_names_gp = gpar(fontsize = 6),
                              row_names_side = "left",
                               col = col_fun_pearson,
                              column_names_rot = 90)

png("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/s_baseline_big_heat.png", height=10, width = 10, units = "in", res=444)
draw(S_baseline_cor_heatmap, padding=unit(c(2,2,2,2), "mm"))
dev.off()



##  asymptomatic baseline ####

A_baseline <- clean_data %>%
  filter(timepoint %in% c("baseline"), infectiontype=="A")%>%
  filter(targetName %in% sig_base_zero_infectiontype$targetName)%>%
  pivot_wider(names_from = targetName, values_from = concentration, id_cols = c(id))

A_baseline_cor <- cor(A_baseline[,2:ncol(A_baseline)], method = "spearman")

# A_baseline_cor_matrix <- as.matrix(A_baseline_cor)
A_baseline_cor_matrix <- as.matrix(A_baseline_cor[rownames(S_baseline_cor)[rev(baseline_hclust$order)[-70]],
                                                  colnames(S_baseline_cor)[rev(baseline_hclust$order)[-70]]])



A_baseline_cor_heatmap <- Heatmap(matrix = A_baseline_cor_matrix,
                              cluster_rows = FALSE,
                              cluster_columns=FALSE,
                              show_row_dend = FALSE,
                              # show_column_dend = TRUE,
                              # show_heatmap_legend = TRUE,
                              name = "Symtpomatic Day 0",
                              #cluster_columns = FALSE,
                              column_names_gp = gpar(fontsize = 6),
                              row_names_gp = gpar(fontsize = 6),
                              row_names_side = "left",
                              col = col_fun_pearson,
                              column_names_rot = 90)

png("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/A_baseline_big_heat.png", height=10, width = 10, units = "in", res=444)
draw(A_baseline_cor_heatmap, padding=unit(c(2,2,2,2), "mm"))
dev.off()




##  all data ####
overall <- clean_data %>%
  filter(timepoint %in% c("baseline", "day0", "day14"), infectiontype%in%c("A", "S"))%>%
  filter(targetName %in% sig_base_zero_infectiontype$targetName)%>%
  filter(!is.na(parasitedensity))%>%
  pivot_wider(names_from = targetName, values_from = concentration, id_cols = c(sample_id, parasitedensity))

overall_cor <- cor(overall[,2:ncol(overall)], method = "spearman")

day0_dist <- dist(overall_cor, method = "euclidean", diag = FALSE, upper = FALSE)
day0_hclust <- hclust(baseline_dist)

# overall_cor_matrix <- as.matrix(overall_cor)
overall_cor_matrix <- as.matrix(overall_cor[rownames(overall_cor)[rev(day0_hclust$order)],colnames(overall_cor)[rev(day0_hclust$order)]])



overall_cor_heatmap <- Heatmap(matrix = overall_cor_matrix,
                              cluster_rows = TRUE,
                              cluster_columns=TRUE,
                              show_row_dend = FALSE,
                              show_column_dend = TRUE,
                              show_heatmap_legend = TRUE,
                              name = "All data",
                              #cluster_columns = FALSE,
                              column_names_gp = gpar(fontsize = 6),
                              row_names_gp = gpar(fontsize = 6),
                              row_names_side = "left",
                              col = col_fun_pearson,
                              column_names_rot = 90)

png("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/overall_big_heat.png", height=10, width = 10, units = "in", res=444)
draw(overall_cor_heatmap, padding=unit(c(2,2,2,2), "mm"))
dev.off()



##  all symptomatic data ####
overall_S <- clean_data %>%
  filter(timepoint %in% c("baseline", "day0", "day7", "day14"), infectiontype %in% c("S"))%>%
  filter(targetName %in% sig_base_zero_infectiontype$targetName)%>%
  filter(!is.na(parasitedensity))%>%
  pivot_wider(names_from = targetName, values_from = concentration, id_cols = c(sample_id, parasitedensity))

overall_S_cor <- cor(overall_S[,2:ncol(overall_S)], method = "spearman")

overall_s_dist <- dist(overall_S_cor, method = "euclidean", diag = FALSE, upper = FALSE)
overall_s_hclust <- hclust(overall_s_dist)

# overall_S_cor_matrix <- as.matrix(overall_S_cor)
overall_S_cor_matrix <- as.matrix(overall_S_cor[rownames(overall_S_cor)[rev(overall_s_hclust$order)],colnames(overall_S_cor)[rev(overall_s_hclust$order)]])


overall_S_cor_heatmap <- Heatmap(matrix = overall_S_cor_matrix,
                               cluster_rows = FALSE,
                               cluster_columns=FALSE,
                               # show_row_dend = FALSE,
                               # show_column_dend = TRUE,
                               show_heatmap_legend = TRUE,
                               name = "Overall Symptomatic",
                               #cluster_columns = FALSE,
                               column_names_gp = gpar(fontsize = 6),
                               row_names_gp = gpar(fontsize = 6),
                               row_names_side = "left",
                               col = col_fun_pearson,
                               column_names_rot = 90)

png("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/overall_S_big_heat.png", height=10, width = 10, units = "in", res=444)
draw(overall_S_cor_heatmap, padding=unit(c(2,2,2,2), "mm"))
dev.off()

##  all asymptomatic data ####
overall_A <- clean_data %>%
  filter(timepoint %in% c("baseline", "day0", "day14"), infectiontype%in%c("A"))%>%
  filter(targetName %in% sig_base_zero_infectiontype$targetName)%>%
  filter(!is.na(parasitedensity))%>%
  pivot_wider(names_from = targetName, values_from = concentration, id_cols = c(sample_id, parasitedensity))

overall_A_cor <- cor(overall_A[,2:ncol(overall_A)], method = "spearman")
# 
# day0_dist <- dist(overall_A_cor, method = "euclidean", diag = FALSE, upper = FALSE)
# day0_hclust <- hclust(baseline_dist)

# overall_A_cor_matrix <- as.matrix(overall_A_cor)
overall_A_cor_matrix <- as.matrix(overall_A_cor[rownames(overall_A_cor)[rev(overall_s_hclust$order)],colnames(overall_A_cor)[rev(overall_s_hclust$order)]])



overall_A_cor_heatmap <- Heatmap(matrix = overall_A_cor_matrix,
                                 cluster_rows = FALSE,
                                 cluster_columns=FALSE,
                                 # show_row_dend = FALSE,
                                 # show_column_dend = TRUE,
                                 show_heatmap_legend = TRUE,
                                 name = "Overall Asymptomatic",
                                 #cluster_columns = FALSE,
                                 column_names_gp = gpar(fontsize = 6),
                                 row_names_gp = gpar(fontsize = 6),
                                 row_names_side = "left",
                                 col = col_fun_pearson,
                                 column_names_rot = 90)

png("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/overall_A_big_heat.png", height=10, width = 10, units = "in", res=444)
draw(overall_A_cor_heatmap, padding=unit(c(2,2,2,2), "mm"))
dev.off()
