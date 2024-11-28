library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")

wide_df2 <- clean_data %>%
  filter(targetName !="CTSS")%>%
  pivot_wider(names_from = targetName, values_from = concentration, id_cols = c("id", "timepoint", "infectiontype", "sample_id"))


rownames(wide_df2) <- wide_df2$sample_id
# each row = measurement; each column = feature
big_pca <-  prcomp(wide_df2[,(length(id_columns)+1):ncol(wide_df2)], center = T)
pca_plot_data <- as.data.frame(cbind(wide_df2, big_pca$x))

big_indie_pca <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  ggplot(., aes(x=PC1, y=PC2, color=timepoint))+
  geom_point()+
  facet_wrap(~id+infectiontype)+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=viridis::magma(n=5))+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/big_indie_pca.png", big_indie_pca, width = 10, height = 10, dpi=444, bg="white")


distance_traveled <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  select(id, timepoint, infectiontype, PC1, PC2)%>%
  pivot_wider(names_from = timepoint, values_from = c(PC1, PC2))%>%
  mutate(base_d0_PC1=PC1_day0-PC1_baseline,
         base_d0_PC2=PC2_day0-PC2_baseline,
         base_d14_PC1=PC1_day0-PC1_baseline,
         base_d14_PC2=PC2_day0-PC2_baseline,
         base_d0_distance=sqrt(base_d0_PC1^2 + base_d0_PC2^2),
         base_d14_distance=sqrt(base_d14_PC1^2 + base_d14_PC2^2))


id_levels <- distance_traveled%>%
  filter(infectiontype=="S", !is.na(base_d0_PC1))%>%
  distinct(id, base_d0_distance)%>%
  arrange(base_d0_distance)


lollipop_s <- distance_traveled%>%
  filter(infectiontype %in% c("S"), !is.na(base_d0_PC1))%>%
  arrange(base_d0_distance)%>%
  mutate(id=factor(id, levels=id_levels$id))%>%
  ggplot(., aes(x=base_d0_distance, y=id, fill=id))+
  geom_bar(stat="identity")+
  theme_minimal()+
  facet_wrap(~infectiontype)+
  scale_x_continuous(limits=c(0, 35))+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank())+
  scale_fill_manual(values = viridis::viridis(n=nrow(id_levels)))


lollipop_a <- distance_traveled%>%
  filter(infectiontype %in% c("A"), id %in% lollipop_s$data$id)%>%
  mutate(base_d0_PC1=ifelse(is.na(base_d0_PC1), 0,base_d0_PC1))%>% 
  # mutate(id=factor(id, levels=id_levels$id))%>%
  ggplot(., aes(x=-base_d0_distance, y=id, fill=id))+
  geom_bar(stat="identity")+
  theme_minimal()+
  facet_wrap(~infectiontype)+
  scale_x_continuous(limits=c(-35, 0))+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank())+
  scale_fill_manual(values = viridis::viridis(n=nrow(id_levels)))

lollipop_a + lollipop_s


distance_traveled %>%
  filter(infectiontype %in% c("A", "S"), !is.na(base_d0_PC1))%>%
  pivot_wider(names_from = infectiontype, values_from = c(base_d0_distance), id_cols=(id))%>%
  ggplot(., aes(x=S, y=A))+
    geom_point()+
  xlab("distance traveled S")+
  ylab("distance traveled A")+
  ggpubr::stat_cor()+
    theme_minimal()
  
