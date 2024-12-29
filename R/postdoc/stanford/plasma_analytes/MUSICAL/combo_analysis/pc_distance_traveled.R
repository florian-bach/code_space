library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)

time_cols <- list("baseline"="#E4DEBD",
                  "day0" = "#C03F3E",
                  "day7" = "#D87E1F",
                  "day14" = "#E6B85F")
clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")

combo_as_results <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/differential_abundance/combo_as_purff.csv", check.names = FALSE)

sig_targets <- combo_as_results %>%
  filter(padj<0.05, contrast %in% c("baseline S - day0 S", "baseline A - day0 A", "baseline A - day14 A"))%>%
  distinct(targetName)

controllers <- clean_data %>%
  filter(infectiontype %in% c("A", "S"), timepoint=="day14")%>%
  group_by(infectiontype)%>%
  mutate(day14_para=if_else(parasitedensity > 10, "parasitemic_day14", "no_parasites_day14"))%>%
  distinct(id, infectiontype, day14_para)
  

wide_df2 <- clean_data %>%
  filter(targetName %in% sig_targets$targetName)%>%
  pivot_wider(names_from = targetName, values_from = concentration, id_cols = c("id", "timepoint", "infectiontype", "sample_id"))%>%
  left_join(., controllers, by=c("id", "infectiontype"))

id_columns <- c("sample_id", "id_cat", "gender_categorical", "ageyrs", "timepoint", "infectiontype", "plate_number", "qpcr", "qpcr_cat", "day14_para")

rownames(wide_df2) <- wide_df2$sample_id
# each row = measurement; each column = feature
big_pca <-  prcomp(wide_df2[, unique(sig_targets$targetName)], center = T)
pca_plot_data <- as.data.frame(cbind(wide_df2, big_pca$x))

big_indie_pca <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"), infectiontype %in% c("A", "S"))%>%
  ggplot(., aes(x=PC1, y=PC2, color=timepoint))+
  geom_point()+
  facet_wrap(~id+infectiontype)+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=time_cols)+
  theme(legend.title = element_blank(),
        axis.text = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1))

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/big_indie_sig_only_pca.png", big_indie_pca, width = 10, height = 10, dpi=444, bg="white")


big_indie_pca2 <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"), infectiontype=="A")%>%
  ggplot(., aes(x=PC1, y=PC2, color=timepoint))+
  geom_point()+
  facet_wrap(~id+infectiontype)+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=time_cols)+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/big_indie_asymp_pca.png", big_indie_pca2, width = 10, height = 10, dpi=444, bg="white")


big_indie_pca3 <-

pca_plot_data2 <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"), infectiontype=="S")

  ggplot(pca_plot_data2, aes(x=PC1, y=PC2, color=timepoint, group=factor(id)))+
  geom_point()+
  geom_segment(arrow = arrow(), aes(x=PC1[timepoint=="baseline"],
                                    y=PC2[timepoint=="baseline"],
                                    xend=PC1[timepoint=="day 0"],
                                    yend=PC2[timepoint=="day 0"], group=factor(id)))+
  facet_wrap(~id+infectiontype)+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=time_cols)+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/big_indie_symp_pca.png", big_indie_pca3, width = 10, height = 10, dpi=444, bg="white")


big_pca_plot <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"), infectiontype %in% c("A", "S"))%>%
  ggplot(., aes(x=PC1, y=PC2, fill=factor(timepoint, levels=c("baseline", "day0", "day7", "day14"))))+
  geom_point(size=3, shape=21)+
  facet_wrap(~infectiontype)+
  # geom_contour(aes(color=timepoint))+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_fill_manual(values=time_cols)+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/big_pca.png", big_pca_plot, width = 8, height = 5, dpi=444, bg="white")


distance_traveled <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  select(id, timepoint, day14_para, infectiontype, PC1, PC2)%>%
  pivot_wider(names_from = timepoint, values_from = c(PC1, PC2))%>%
  mutate(base_d0_PC1=abs(PC1_day0-PC1_baseline),
         base_d0_PC2=abs(PC2_day0-PC2_baseline),
         base_d14_PC1=abs(PC1_day0-PC1_baseline),
         base_d14_PC2=abs(PC2_day0-PC2_baseline),
         base_d0_distance=sqrt(base_d0_PC1^2 + base_d0_PC2^2),
         base_d14_distance=sqrt(base_d14_PC1^2 + base_d14_PC2^2))


## pc1 & pc2 distance ####

id_levels <- distance_traveled%>%
  filter(infectiontype=="S", !is.na(base_d0_PC1))%>%
  distinct(id, base_d0_distance)%>%
  arrange(base_d0_distance)

overall_distance_lollipop <- distance_traveled%>%
  filter(infectiontype %in% c("A", "S"), !is.na(base_d0_distance), id %in% id_levels$id)%>%
  mutate(base_d0_distance=if_else(infectiontype=="A", -base_d0_distance, base_d0_distance))%>%
  arrange(base_d0_distance)%>%
  mutate(id=factor(id, levels=id_levels$id))%>%
  ggplot(., aes(x=base_d0_distance, y=id, fill=id))+
  geom_bar(stat="identity")+
  theme_minimal()+
  ggtitle("PC distance travelled from baseline to day 0")+
  facet_wrap(~infectiontype, scales="free_x")+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank())+
  scale_fill_manual(values = viridis::viridis(n=nrow(id_levels)))


ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/overall_distance_lollipop.png", overall_distance_lollipop, bg="white", height=8, width = 5, dpi=444)

# dotplot
overall_distance_dotplot <- distance_traveled %>%
  pivot_wider(names_from = infectiontype, values_from = base_d0_distance, id_cols = id)%>%
  ggplot(., aes(x=S, y=A))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_minimal()+
  ggpubr::stat_cor(method="spearman")+
  ggtitle("PC distance travelled from baseline to day 0")+
  theme(plot.title = element_text(size=13))

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/overall_distance_dotplot.png", overall_distance_dotplot, bg="white", height=4, width = 4, dpi=444)

## pc1  distance ####
id_levels <- distance_traveled%>%
  filter(infectiontype=="S", !is.na(base_d0_PC1))%>%
  distinct(id, base_d0_PC1)%>%
  arrange(abs(base_d0_PC1))

pc1_distance_lollipop <- distance_traveled%>%
  filter(infectiontype %in% c("A","S"), id %in% id_levels$id)%>%
  mutate(base_d0_PC1=if_else(infectiontype=="A", -base_d0_PC1, base_d0_PC1))%>%
  arrange(base_d0_PC1)%>%
  mutate(id=factor(id, levels=id_levels$id))%>%
  ggplot(., aes(x=base_d0_PC1, y=id, fill=id))+
  geom_bar(stat="identity")+
  theme_minimal()+
  ggtitle("PC1 distance travelled from baseline to day 0")+
  facet_wrap(~infectiontype, scales = "free_x")+
  # scale_x_continuous(limits=c(0, 35))+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank())+
  scale_fill_manual(values = viridis::viridis(n=nrow(id_levels)))


ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/pc1_distance_lollipop.png", pc1_distance_lollipop, bg="white", height=8, width = 5, dpi=444)

pc1_dotplot <- distance_traveled %>%
  pivot_wider(names_from = infectiontype, values_from = base_d0_PC1, id_cols = id)%>%
  ggplot(., aes(x=S, y=A))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_minimal()+
  ggpubr::stat_cor(method="spearman")+
  ggtitle("PC1 distance travelled from baseline to day 0")+
  theme(plot.title = element_text(size=13))

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/pc1_distance_dotplot.png", pc1_dotplot, bg="white", height=4, width = 4, dpi=444)


## pc2  distance ####
id_levels <- distance_traveled%>%
  filter(infectiontype=="S", !is.na(base_d0_PC2))%>%
  distinct(id, base_d0_PC2)%>%
  arrange(abs(base_d0_PC2))

pc2_distance_lollipop <- distance_traveled%>%
  filter(infectiontype %in% c("A","S"), id %in% id_levels$id)%>%
  mutate(base_d0_PC2=if_else(infectiontype=="A", -base_d0_PC2, base_d0_PC2))%>%
  arrange(base_d0_PC2)%>%
  mutate(id=factor(id, levels=id_levels$id))%>%
  ggplot(., aes(x=base_d0_PC2, y=id, fill=id))+
  geom_bar(stat="identity")+
  theme_minimal()+
  ggtitle("PC2 distance travelled from baseline to day 0")+
  facet_wrap(~infectiontype, scales="free_x")+
  # scale_x_continuous(limits=c(0, 35))+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank())+
  scale_fill_manual(values = viridis::viridis(n=nrow(id_levels)))

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/pc2_distance_lollipop.png", pc2_distance_lollipop, bg="white", height=8, width = 5, dpi=444)



pc2_dotplot <- distance_traveled %>%
  pivot_wider(names_from = infectiontype, values_from = base_d0_PC2, id_cols = id)%>%
  ggplot(., aes(x=S, y=A))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_minimal()+
  ggpubr::stat_cor(method="spearman")+
  ggtitle("PC2 distance travelled from baseline to day 0")+
  theme(plot.title = element_text(size=13))

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/pc1_distance_dotplot.png", pc2_dotplot, bg="white", height=4, width = 4, dpi=444)

## pc 2 distance at 14 ####

id_levels <- distance_traveled%>%
  filter(infectiontype=="A", !is.na(base_d14_PC2))%>%
  distinct(id, base_d14_PC2)%>%
  arrange(desc(abs(base_d14_PC2)))

pc2_day14_distance_lollipop <- distance_traveled%>%
  filter(infectiontype %in% c("A","S"), id %in% id_levels$id)%>%
  mutate(day0_s_day14_a=if_else(infectiontype=="A", -base_d14_distance, base_d0_distance))%>%
  arrange(base_d14_distance)%>%
  mutate(id=factor(id, levels=rev(id_levels$id)))%>%
  ggplot(., aes(x=day0_s_day14_a, y=id, fill=id, color = day14_para))+
  geom_bar(stat="identity", linewidth = 0.51)+
  theme_minimal()+
  ggtitle("distance travelled from baseline to day0 (S) day 14 (A)")+
  facet_wrap(~infectiontype, scales="free_x")+
  scale_color_manual(values = c("white", "black"))+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank())+
  scale_fill_manual(values = viridis::viridis(n=nrow(id_levels)))

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/pc2_day14_distance_lollipop.png", pc2_day14_distance_lollipop, bg="white", height=8, width = 5, dpi=444)

# 
base_d14_dotplot <- distance_traveled %>%
  mutate(base_d14_PC2=if_else(infectiontype=="A", base_d14_distance, base_d0_distance))%>%
  pivot_wider(names_from = infectiontype, values_from = base_d14_distance, id_cols = c(id))%>%
  ggplot(., aes(x=S, y=A))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_minimal()+
  ggpubr::stat_cor(method="spearman")+
  ggtitle("PC distance travelled from baseline to day 14")+
  theme(plot.title = element_text(size=13))

# distance_traveled %>%
#   mutate(base_d14_PC2=if_else(infectiontype=="A", base_d14_distance, base_d0_distance))%>%
#   pivot_wider(names_from = infectiontype, values_from = base_d14_distance, id_cols = c(id))%>%
#   ggplot(., aes(x=S, y=A))+
#   geom_point()+
#   geom_smooth(method="lm")+
#   theme_minimal()+
#   ggpubr::stat_cor(method="spearman")+
#   ggtitle("PC distance travelled from baseline to day 14")+
#   theme(plot.title = element_text(size=13))
# 
# ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/pc1_distance_dotplot.png", pc2_dotplot, bg="white", height=4, width = 4, dpi=444)





distance_traveled %>%
  filter(infectiontype %in% c("A", "S"), !is.na(base_d0_PC1))%>%
  pivot_wider(names_from = infectiontype, values_from = c(base_d0_distance), id_cols=(id))%>%
  ggplot(., aes(x=S, y=A))+
  geom_point()+
  xlab("distance traveled S")+
  ylab("distance traveled A")+
  ggpubr::stat_cor(method="spearman")+
  theme_minimal()


distance_traveled %>%
  filter(infectiontype %in% c("A", "S"), !is.na(base_d0_PC1))%>%
  pivot_wider(names_from = infectiontype, values_from = c(base_d0_PC1), id_cols=(id))%>%
  ggplot(., aes(x=S, y=A))+
  geom_point()+
  xlab("distance traveled PC1 S")+
  ylab("distance traveled PC1 A")+
  ggpubr::stat_cor(method="spearman")+
  theme_minimal()

distance_traveled %>%
  filter(infectiontype %in% c("A", "S"), !is.na(base_d0_PC2))%>%
  pivot_wider(names_from = infectiontype, values_from = c(base_d0_PC2), id_cols=(id))%>%
  ggplot(., aes(x=S, y=A))+
  geom_point()+
  xlab("distance traveled PC2 S")+
  ylab("distance traveled PC2 A")+
  ggpubr::stat_cor(method="spearman")+
  theme_minimal()

  
