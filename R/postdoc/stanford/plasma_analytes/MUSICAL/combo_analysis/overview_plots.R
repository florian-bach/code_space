library(tidyr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(patchwork)

`%notin%` <- Negate(`%in%`)

da_boxplot_theme <- theme(legend.position = "none",
                          axis.title = element_blank())

#unclean data
unclean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/unclean_musical_combo_with_metadata.csv")
unclean_data <- unclean_data %>%
  mutate(timepoint = factor(timepoint, levels=c("bad_baseline", "baseline", "day0", "day7", "day14", "day28")))

#clean data
clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")
clean_data <- clean_data %>%
  mutate(timepoint = factor(timepoint, levels=c("bad_baseline", "baseline", "day0", "day7", "day14", "day28")))

# big heatmap ####
## unclean ####
t_nulisa <-t(as.matrix(unclean_data %>%
                         # filter(mean_z_conc > -0.4)%>%
                         select(sample_id, targetName, z_conc)%>%
                         pivot_wider(names_from = targetName, values_from = z_conc)))


colnames(t_nulisa) <- t_nulisa[1,]

num_nulisa_matrix <- t_nulisa[-1,]
class(num_nulisa_matrix) <- "numeric"

short_num_nulisa_matrix <- num_nulisa_matrix

# z_col_fun <- circlize::colorRamp2(c(min(short_num_nulisa_matrix)/1.5, 0, abs(min(short_num_nulisa_matrix)/1.5)), c("#123499", "white", "#FFA500"))
z_col_fun <- circlize::colorRamp2(c(min(-10.90915)/1.5, 0, abs(min(-10.90915)/1.5)), c("#123499", "white", "#FFA500"))

# big_heatmap_split <- substr(colnames(t_nulisa), 5, nchar(colnames(t_nulisa)))
# 
# big_heatmap_split <- factor(big_heatmap_split, levels=sort(unique(big_heatmap_split))[c(1:5, 7, 6)])


big_heatmap <- Heatmap(matrix = short_num_nulisa_matrix,
                       cluster_rows = TRUE,
                       cluster_columns=TRUE,
                       show_row_dend = FALSE,
                       show_column_dend = FALSE,
                       show_heatmap_legend = TRUE,
                       # column_split = big_heatmap_split,
                       name = "z score",
                       column_names_gp = gpar(fontsize = 3),
                       row_names_gp = gpar(fontsize = 0),
                       row_names_side = "left",
                       col = z_col_fun,
                       column_names_rot = 90)



png("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/unclean_z_heat.png", width=16, height = 9, res = 444, units = "in")
draw(big_heatmap)
dev.off()

## clean data ####

t_nulisa <-t(as.matrix(clean_data %>%
                         # filter(mean_z_conc > -0.4)%>%
                         select(sample_id, targetName, z_conc)%>%
                         pivot_wider(names_from = targetName, values_from = z_conc)))


colnames(t_nulisa) <- t_nulisa[1,]

num_nulisa_matrix <- t_nulisa[-1,]
class(num_nulisa_matrix) <- "numeric"

short_num_nulisa_matrix <- num_nulisa_matrix

# z_col_fun <- circlize::colorRamp2(c(min(short_num_nulisa_matrix)/1.5, 0, abs(min(short_num_nulisa_matrix)/1.5)), c("#123499", "white", "#FFA500"))
z_col_fun <- circlize::colorRamp2(c(min(-10.90915)/1.5, 0, abs(min(-10.90915)/1.5)), c("#123499", "white", "#FFA500"))

# big_heatmap_split <- substr(colnames(t_nulisa), 5, nchar(colnames(t_nulisa)))
# 
# big_heatmap_split <- factor(big_heatmap_split, levels=sort(unique(big_heatmap_split))[c(1:5, 7, 6)])


big_heatmap <- Heatmap(matrix = short_num_nulisa_matrix,
                       cluster_rows = TRUE,
                       cluster_columns=TRUE,
                       show_row_dend = FALSE,
                       show_column_dend = FALSE,
                       show_heatmap_legend = TRUE,
                       # column_split = big_heatmap_split,
                       name = "z score",
                       column_names_gp = gpar(fontsize = 3),
                       row_names_gp = gpar(fontsize = 0),
                       row_names_side = "left",
                       col = z_col_fun,
                       column_names_rot = 90)



png("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/clean_z_heat.png", width=16, height = 9, res = 444, units = "in")
draw(big_heatmap)
dev.off()


# PCA ####
  ## unclean ####

id_columns <- c("sample_id", "id", "gender_categorical", "ageyrs", "timepoint", "infectiontype", "plate_number", "qpcr", "qpcr_cat")

unclean_wide_df2 <- unclean_data %>%
  # filter(targetName %in% most_variable$targetName)%>%
  pivot_wider(names_from = targetName, values_from = concentration, id_cols = all_of(id_columns))


rownames(unclean_wide_df2) <- unclean_wide_df2$sample_id
# each row = measurement; each column = feature
unclean_big_pca <-  prcomp(unclean_wide_df2[,(length(id_columns)+1):ncol(unclean_wide_df2)], center = T)
unclean_pca_plot_data <- as.data.frame(cbind(unclean_wide_df2, unclean_big_pca$x))

removed_samples <- unique(unclean_data$sample_id[unclean_data$sample_id %notin% clean_data$sample_id])

post_hoc <- unique(unclean_data$sample_id[unclean_data$barcode%in%c("D1RTFU", "D12FNR")])
  
pca_removal_plot <- unclean_pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  mutate("over_six"=if_else(ageyrs>6, "over 6y", "under 6y"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(., aes(x=PC1, y=PC2))+
  geom_point(aes(color=ifelse(sample_id %in% removed_samples, "removed", "kept")
                 ))+
  # geom_text(aes(label=sample_id))+
  # stat_ellipse()+
  # xlab(paste("PC1 ", data.frame(summary(unclean_big_pca)[6])[2,1]*100, "%", sep = ""))+
  # ylab(paste("PC2 ", data.frame(summary(unclean_big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=c("black", "red"))+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/pca_removal_plot.png", width=6, height = 6, dpi=444, bg="white")



pca_removal_plot2 <- unclean_pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  mutate("over_six"=if_else(ageyrs>6, "over 6y", "under 6y"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(., aes(x=PC3, y=PC4))+
  geom_point(aes(color=CTSS))+
  # stat_ellipse()+
  xlab(paste("PC3 ", data.frame(summary(unclean_big_pca)[6])[2,3]*100, "%", sep = ""))+
  ylab(paste("PC4 ", data.frame(summary(unclean_big_pca)[6])[2,4]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_viridis_c()+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/pca_removal_plot2.png", width=6, height = 6, dpi=444, bg="white")



loop_data <- unclean_pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  mutate("over_six"=if_else(ageyrs>6, "over 6y", "under 6y"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))

plot_list <- list()

for(color_by in c("gender_categorical", "infectiontype", "over_six", "plate_number", "qpcr_cat")){
  for(i in seq(1,7, by=2)){
    
    xvar <- paste("PC", i, sep="")
    yvar <- paste("PC", i+1, sep="")
    
    # print(head(loop_data[,xvar]))
    # print(head(loop_data[,yvar]))
    
    plt <- ggplot(loop_data, aes_string(x=xvar, y=yvar, color=color_by))+
      geom_point()+
      xlab(paste0("PC", i, " ", data.frame(summary(unclean_big_pca)[6])[2,i]*100, "%", sep = ""))+
      ylab(paste0("PC", i+1, " ", data.frame(summary(unclean_big_pca)[6])[2,i+1]*100, "%", sep = ""))+
      theme_minimal()+
      theme(legend.title = element_blank())+
      scale_color_manual(values=viridis::magma(length(unique(loop_data[,color_by]))+1))
    #ggrepel::geom_label_repel(aes_string(label = "sample_id"), show.legend = FALSE)
    
    plot_list[[i]] <- plt
    
  }
  big_plot <- plot_list[[1]] + plot_list[[3]] +
    plot_list[[5]] + plot_list[[7]] 
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/unclean_many_pca_", color_by, ".png", sep=""), big_plot, height = 12, width=12, dpi=444, bg="white")
  
  small_plot <-  plot_list[[3]]
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/unclean_pc3_", color_by, ".png", sep=""), small_plot, height = 6, width=6, dpi=444, bg="white")
  
}

  ## clean ####

id_columns <- c("sample_id", "id", "gender_categorical", "ageyrs", "timepoint", "infectiontype", "plate_number", "qpcr", "qpcr_cat")
wide_df2 <- clean_data %>%
  filter(targetName !="CTSS")%>%
  # filter(targetName %in% most_variable$targetName)%>%
  filter(sample_id %notin% c("X384_S.1_D1V93U", "X744_A.1_D1KT5Z", "X323_A14_D1DZ2D", "X496_NM7_D1CAYS", "X316_S.1_D12FNR", "X667_NM7_D1GBYA", "X 176 S_t14 D1FXRJ", "X164_NM0_D1WA6Q"))%>%
  pivot_wider(names_from = targetName, values_from = concentration, id_cols = all_of(id_columns))


rownames(wide_df2) <- wide_df2$sample_id
# each row = measurement; each column = feature
big_pca <-  prcomp(wide_df2[,(length(id_columns)+1):ncol(wide_df2)], center = T)
pca_plot_data <- as.data.frame(cbind(wide_df2, big_pca$x))



loadings_df <- data.frame(big_pca$rotation)
loadings_df$targetName <- rownames(loadings_df)
loadings_df$targetName <- factor(loadings_df$targetName, levels = loadings_df$targetName[order(loadings_df$PC1)])

pc_cols <- colorspace::sequential_hcl(nrow(loadings_df), palette = "Purple Yellow")
names(pc_cols) <- loadings_df$targetName[order(loadings_df$PC1)]

PC3_plot <- ggplot(loadings_df, aes(x=factor(targetName, levels = targetName[order(loadings_df$PC3)]), y=PC3, fill=targetName))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = pc_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")

PC8_plot <- ggplot(loadings_df, aes(x=factor(targetName, levels = targetName[order(loadings_df$PC8)]), y=PC8, fill=targetName))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = pc_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")


age_pca_plot <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  mutate("over_six"=if_else(ageyrs>6, "over 6y", "under 6y"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(., aes(x=PC1, y=PC2, color=over_six))+
  geom_point()+
  # stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=c("violet", "chartreuse"))+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/age_pca_plot.png", age_pca_plot, width=4, height=4, bg="white", dpi=444)


sex_pca_plot <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  mutate("over_six"=if_else(ageyrs>6, "over 6y", "under 6y"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(., aes(x=PC1, y=PC2, color=gender_categorical))+
  geom_point()+
  # stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=c("darkblue", "orange"))+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/sex_pca_plot.png", sex_pca_plot, width=4, height=4, bg="white", dpi=444)


time_pca_plot <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  mutate("over_six"=if_else(ageyrs>6, "over 6y", "under 6y"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(., aes(x=PC1, y=PC2, color=timepoint))+
  geom_point()+
  # stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=viridis::magma(5))+
  # facet_wrap(~timepoint)+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/time_pca_plot.png", time_pca_plot, width=4, height=4, bg="white", dpi=444)


plate_number_pca_plot <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  mutate("over_six"=if_else(ageyrs>6, "over 6y", "under 6y"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(., aes(x=PC1, y=PC2, fill=plate_number))+
  geom_point(shape=21)+
  # stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  # scale_fill_manual(values=colorspace::qualitative_hcl(n=5, palette = "dark2"))+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/plate_number_pca_plot.png", plate_number_pca_plot, width=4, height=4, bg="white", dpi=444)


infectiontype_pca_plot <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  mutate("over_six"=if_else(ageyrs>6, "over 6y", "under 6y"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(., aes(x=PC1, y=PC2, color=infectiontype))+
  geom_point()+
  # stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=viridis::magma(6))+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/infectiontype_pca_plot.png", infectiontype_pca_plot, width=4, height=4, bg="white", dpi=444)



qpcr_pca_plot <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  mutate(qpcr_cat=factor(qpcr_cat, levels=c("0", ">1", ">10", ">10e2", ">10e3", ">10e4", ">10e5")))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(., aes(x=PC1, y=PC2, color=qpcr_cat))+
  geom_point()+
  # stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=colorspace::sequential_hcl(n = 7, "LaJolla", rev = TRUE))+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/qpcr_pca_plot.png", qpcr_pca_plot, width=4, height=4, bg="white", dpi=444)

### PC || covariate examination####

loop_data <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  mutate("over_six"=if_else(ageyrs>6, "over 6y", "under 6y"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  mutate(qpcr_cat=factor(qpcr_cat, levels=c("0", ">1", ">10", ">10e2", ">10e3", ">10e4", ">10e5")))
  

plot_list <- list()

for(color_by in c("gender_categorical", "infectiontype", "over_six", "plate_number", "qpcr_cat")){
  for(i in seq(1,7, by=2)){
    
    xvar <- paste("PC", i, sep="")
    yvar <- paste("PC", i+1, sep="")
    
    # print(head(loop_data[,xvar]))
    # print(head(loop_data[,yvar]))
    
    plt <- ggplot(loop_data, aes_string(x=xvar, y=yvar, color=color_by))+
      geom_point()+
      xlab(paste0("PC", i, " ", data.frame(summary(big_pca)[6])[2,i]*100, "%", sep = ""))+
      ylab(paste0("PC", i+1, " ", data.frame(summary(big_pca)[6])[2,i+1]*100, "%", sep = ""))+
      theme_minimal()+
      theme(legend.title = element_blank())+
      scale_color_manual(values=viridis::magma(length(unique(loop_data[,color_by]))+1))
    #ggrepel::geom_label_repel(aes_string(label = "sample_id"), show.legend = FALSE)
    
    plot_list[[i]] <- plt
    
  }
  big_plot <- plot_list[[1]] + plot_list[[3]] +
              plot_list[[5]] + plot_list[[7]] 
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/clean_many_pca_", color_by, ".png", sep=""), big_plot, height = 12, width=12, dpi=444, bg="white")
  
  small_plot <-  plot_list[[3]]
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/clean_pc3_", color_by, ".png", sep=""), small_plot, height = 6, width=6, dpi=444, bg="white")
  
}

# plot parasite densities ####
## qpcr ####
asymp_pcr_plot <- clean_data %>%
  filter(infectiontype%in%c("A"), timepoint %notin% c("day28"), timepoint!="bad_baseline")%>%
  distinct(id, timepoint, qpcr)%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=qpcr+0.1, color=id))+
  geom_point()+
  geom_line(aes(group=id))+
  ggtitle("asymptomatic parasitemia")+
  scale_color_manual(values=viridis::magma(60))+
  scale_y_log10(limits=c(0.1, 10^6))+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/asymp_pcr_plot.png", asymp_pcr_plot, width=5, height=5, dpi=444, bg="white")


symp_pcr_plot <- clean_data %>%
  filter(infectiontype%in%c("S"), timepoint %notin% c("day28"), timepoint!="bad_baseline")%>%
  distinct(id, timepoint, qpcr)%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=qpcr+0.1, color=id))+
  geom_point()+
  geom_line(aes(group=id))+
  ggtitle("symptomatic parasitemia")+
  scale_color_manual(values=viridis::magma(60))+
  scale_y_log10(limits=c(0.1, 10^6))+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/symp_pcr_plot.png", symp_pcr_plot, width=5, height=5, dpi=444, bg="white")


## blood smears ####

asymp_bs_plot <- clean_data %>%
  filter(infectiontype%in%c("A"), timepoint %notin% c("day28"), timepoint!="bad_baseline")%>%
  distinct(id, timepoint, parasitedensity)%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=parasitedensity+0.1, color=id))+
  geom_point()+
  geom_line(aes(group=id))+
  ggtitle("asymptomatic parasitemia")+
  scale_color_manual(values=viridis::magma(60))+
  scale_y_log10(limits=c(0.1, 10^6))+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/asymp_bs_plot.png", asymp_bs_plot, width=5, height=5, dpi=444, bg="white")


symp_bs_plot <- clean_data %>%
  filter(infectiontype%in%c("S"), timepoint %notin% c("day28"), timepoint!="bad_baseline")%>%
  distinct(id, timepoint, parasitedensity)%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=parasitedensity+0.1, color=id))+
  geom_point()+
  geom_line(aes(group=id))+
  ggtitle("symptomatic parasitemia")+
  scale_color_manual(values=viridis::magma(60))+
  scale_y_log10(limits=c(0.1, 10^6))+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/symp_pcr_plot.png", symp_bs_plot, width=5, height=5, dpi=444, bg="white")

## correlating blood smears with pcr####

pardens_qpcr_cor <- clean_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"), infectiontype %in% c("A", "S"))%>%
  distinct(id, timepoint, parasitedensity, qpcr, infectiontype)%>%
  ggplot(aes(x=qpcr+0.1, y=parasitedensity+0.1))+
  geom_point(aes(color=infectiontype))+
  geom_hline(yintercept = 500, linetype="dashed")+
  ggpubr::stat_cor()+
  geom_smooth(method="lm")+
  ggtitle("correlating microscopy with qPCR")+
  scale_color_manual(values=viridis::magma(3))+
  scale_y_log10(limits=c(0.1, 10^6))+
  scale_x_log10(limits=c(0.1, 10^6))+
  facet_wrap(~timepoint)+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/pardens_qpcr_cor.png", pardens_qpcr_cor, width=5, height=5, dpi=444, bg="white")


# analytes that vary by parasitemia, but not infectiontype; 59/95
clean_data %>%
  filter(infectiontype %in% c("A", "S"))%>%
  filter(targetName %in% c("CTLA4", "IL2RA", "TNFRSF1A", "FASLG", "CCL4", "GZMA"))%>%
  ggplot(., aes(x=log_qpcr, y=concentration, color=infectiontype))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~targetName)+
  scale_color_manual(values=viridis::magma(n=3))+
  theme_minimal()

# analytes that vary by timepoint, but not parasitemia;
clean_data %>%
  filter(infectiontype %in% c("A", "S"))%>%
  filter(targetName %in% c("LAMP3", "CLEC4A", "CCL22", "CXADR", "CD276", "CD200"))%>%
  ggplot(., aes(x=log_qpcr, y=concentration, color=infectiontype))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~targetName)+
  scale_color_manual(values=viridis::magma(n=3))+
  theme_minimal()

clean_data %>%
  filter(infectiontype %in% c("A", "S"), timepoint %notin% c("bad_baseline", "day28", "day7"))%>%
  filter(targetName %in% c("LAMP3", "CLEC4A", "CCL22", "CXADR", "CD276", "CD200"))%>%
  ggplot(., aes(x=timepoint, y=concentration, fill=infectiontype))+
  geom_boxplot()+
  # geom_point()+
  # geom_smooth(method="lm")+
  facet_wrap(~targetName)+
  scale_fill_manual(values=viridis::magma(n=3))+
  theme_minimal()

clean_data %>%
  filter(infectiontype %in% c("A", "S"))%>%
  filter(targetName %in% c("IL10", "CXCL10", "IL6", "TNFSF11", "KNG1", "C1QA"))%>%
  ggplot(., aes(x=log_qpcr, y=concentration, color=infectiontype))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~targetName)+
  scale_color_manual(values=viridis::magma(n=3))+
  theme_minimal()


clean_data %>%
  filter(infectiontype %in% c("A", "S"), timepoint %notin% c("bad_baseline", "day28"))%>%
  select(-targetName)%>%
  filter(!duplicated(sample_id))%>%
  ggplot(., aes(x=timepoint, y=log_qpcr, fill=infectiontype))+
  # scale_y_log10()+
  # geom_point(position = position_dodge(width=0.75))+
  # geom_boxplot()+
  geom_violin(draw_quantiles = seq(0,1,0.25))+
  scale_fill_manual(values=viridis::magma(n=3))+
  theme_minimal()


 clean_data %>%
  filter(infectiontype %in% c("A", "S"), timepoint %notin% c("bad_baseline", "day28"))%>%
  select(-targetName)%>%
  filter(!duplicated(sample_id))%>%
  ggplot(., aes(x=timepoint, y=parasitedensity+0.1, fill=infectiontype))+
  scale_y_log10()+
  geom_violin(draw_quantiles = seq(0,1,0.25))+
  scale_fill_manual(values=viridis::magma(n=3))+
  theme_minimal()
  

# sandbox ####
# 
# library(ggstats)
# 
# 
# # plot linear regression final results ####
# for(i in 1:length(symp_only_purff$model)){
#   
#   plt <- ggcoef_table(symp_only_purff$model[[i]],
#                       show_p_values = TRUE)
#   ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/model_plots/symp_", symp_only_purff$targetName[[i]], ".png", sep=""), plt, height=4, width=8, dpi=444, bg="white")
#   
# }
# 
# # plot model selection intermediates ####
# 
# mixed_model_selection <- clean_data %>%
#   filter(timepoint!="day28", timepoint!="bad_baseline")%>%
#   group_by(targetName)%>%
#   mutate(log_qpcr=log10(qpcr+0.1))%>%
#   nest() %>%
#   mutate(simple_model=map(data,  ~lme4::lmer(concentration~timepoint+infectiontype+(1|id), data=.)), simple_AIC=map_dbl(simple_model, ~AIC(.)))%>%
#   mutate(qpcr_model=map(data,     ~lme4::lmer(concentration~timepoint+infectiontype+log_qpcr+(1|id), data=.)), qpcr_AIC=map_dbl(qpcr_model, ~AIC(.)))%>%
#   mutate(age_model=map(data,     ~lme4::lmer(concentration~timepoint+infectiontype+ageyrs+(1|id), data=.)), age_AIC=map_dbl(age_model, ~AIC(.)))%>%
#   mutate(sex_model=map(data,     ~lme4::lmer(concentration~timepoint+infectiontype+gender_categorical+(1|id), data=.)), sex_AIC=map_dbl(sex_model, ~AIC(.)))%>%
#   mutate(age_sex_model=map(data, ~lme4::lmer(concentration~timepoint+infectiontype+ageyrs+gender_categorical+(1|id), data=.)), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))
# 
# 
# 
# for(i in 1:nrow(mixed_model_selection)){
#   
#   models <- list(
#     "simple model" = mixed_model_selection$simple_model[[i]],
#     "qpcr model" = mixed_model_selection$qpcr_model[[i]],
#     "age model" = mixed_model_selection$age_model[[i]],
#     "sex model" = mixed_model_selection$sex_model[[i]],
#     "age sex model" = mixed_model_selection$age_sex_model[[i]]
#   )
#   
#   plt <- ggcoef_compare(models, type = "faceted")
#   ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/model_plots/mixed_model_selection", mixed_model_selection$targetName[[i]], ".png", sep=""), plt, height=6, width=18, dpi=444, bg="white")
#   
# }
# 
# 
# 
# 
# pca_plot_data %>%
#   filter(timepoint %in% c("day0", "day7", "day14", "baseline"), infectiontype %in% c("A", "S"))%>%
#   ggplot(., aes(x=PC1, y=PC2, color=timepoint))+
#   geom_point()+
#   stat_ellipse()+
#   xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
#   ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
#   theme_minimal()+
#   facet_wrap(~infectiontype)+
#   scale_color_manual(values=viridis::magma(5))+
#   theme(legend.title = element_blank(),
#         axis.text = element_blank())
# 
# pca_plot_data %>%
#   filter(timepoint %in% c("baseline"), infectiontype %in% c("A", "S"))%>%
#   ggplot(., aes(x=PC1, y=PC2, color=infectiontype))+
#   geom_point()+
#   stat_ellipse()+
#   xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
#   ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
#   theme_minimal()+
#   scale_color_manual(values=viridis::magma(3))+
#   theme(legend.title = element_blank(),
#         axis.text = element_blank())
# 
# 
# qpcr_cat_at_day0 <- clean_data %>%
#   filter(infectiontype %in% c("A", "S"), !is.na(qpcr_cat), timepoint %in% c("baseline", "day0", "day7", "day14"))%>%
#   group_by(infectiontype, id)%>%
#   reframe("day0_qpcr_cat"=qpcr_cat[timepoint=="day0"],
#           "day0_qpcr"=log10(qpcr+0.1)[timepoint=="day0"])%>%
#   distinct()
# 
# 
# clean_data %>%
#   filter(infectiontype%in%c("A", "S"), timepoint!="day28", timepoint!="bad_baseline", !is.na(qpcr_cat))%>%
#   left_join(., qpcr_cat_at_day0, by=c("infectiontype", "id"))%>%
#   filter(targetName %in% c("IFNG", "TNF", "CCL4"))%>%
#   # filter(targetName %in% c("IL15", "SDC1", "EPO", "TNFRSF8", "PDCD1", "IRAK4", "CXCL11""ANGPT2"))%>%
#   mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
#   mutate(qpcr_cat = factor(qpcr_cat, levels=c("0", ">1", ">10", ">10e2", ">10e3", ">10e4", ">10e5")))%>%
#   mutate("over_six"=if_else(ageyrs>6, "over 6y", "under 6y"))%>%
#   ggplot(aes(x=factor(timepoint), y=concentration, fill=day0_qpcr))+
#   ggtitle("day0_qpcr")+
#   # geom_point()+#
#   # geom_line(aes(group=id))+
#   geom_boxplot(outliers = FALSE)+
#   # geom_violin(draw_quantiles = seq(0,1,0.25))+
#   # ggtitle("upregulated on day 7 post malaria")+
#   facet_wrap(~targetName+infectiontype, scales = "free", ncol=2)+
#   scale_fill_manual(values=viridis::magma(8))+
#   theme_minimal()+
#   da_boxplot_theme+
#   theme(legend.position = "right")
# 
# ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/il10_tnf_qpcr_cat_plot.png", il10_tnf_qpcr_cat_plot, width=8, height=4, dpi=444, bg='white')
# 
# pca_plot_data %>%
#   filter(infectiontype %in% c("A", "S"))%>%
#   left_join(., qpcr_cat_at_day0, by=c("infectiontype", "id"))%>%
#   ggplot(., aes(x=PC1, y=PC2, color=day0_qpcr))+
#   geom_point()+
#   # stat_ellipse()+
#   xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
#   ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
#   theme_minimal()+
#   scale_color_manual(values=viridis::magma(7))+
#   theme(legend.title = element_blank(),
#         axis.text = element_blank())
# 
# 
# # check whether
# clean_data %>%
#   filter(infectiontype%in%c("S"), timepoint!="day28", timepoint!="bad_baseline", !is.na(qpcr_cat))%>%
#   left_join(., qpcr_cat_at_day0, by=c("infectiontype", "id"))%>%
#   # filter(targetName %in% c("IL15", "SDC1", "EPO", "TNFRSF8", "PDCD1", "IRAK4", "CXCL11""ANGPT2"))%>%
#   mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
#   mutate(qpcr_cat = factor(qpcr_cat, levels=c("0", ">1", ">10", ">10e2", ">10e3", ">10e4", ">10e5")))%>%
#   mutate("over_six"=if_else(ageyrs>6, "over 6y", "under 6y"))%>%
#   filter(targetName %in% sig_para$targetName[9:16])%>%
#   filter(timepoint %in% c("baseline"))%>%
#   ggplot(aes(x=concentration, y=day0_qpcr, color=timepoint))+
#   ggtitle("day0_qpcr")+
#   # geom_point()+
#   geom_smooth(method="lm")+
#   ggpubr::stat_cor(method = "spearman")+
#   geom_point()+#
#   # geom_line(aes(group=id))+
#   # geom_boxplot(outliers = FALSE)+
#   # geom_violin(draw_quantiles = seq(0,1,0.25))+
#   # ggtitle("upregulated on day 7 post malaria")+
#   facet_wrap(~targetName+infectiontype, scales = "free", ncol=4)+
#   scale_fill_manual(values=viridis::magma(8))+
#   theme_minimal()

