library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbiplot)
library(plotly)
library(shadowtext)

#library(rgl) # want some 3d?

#data <- read.csv("/Users/s1249052/PhD/plasma/vac69a/Vivax_plasma_analytes2_no_inequalities.csv")

plasma <- read.csv("~/PhD/plasma/vac69a/Vivax_plasma_analytes2_no_inequalities.csv", header=T, stringsAsFactors = F)
#plasma <- filter(plasma, Volunteer!="v009")

long_data <- gather(plasma, Analyte, Concentration, colnames(plasma)[3:ncol(plasma)])

long_data$Analyte <- substr(long_data$Analyte, 1, nchar(long_data$Analyte)-7)

long_data$Analyte <- gsub(".", "", long_data$Analyte, fixed = T)
long_data$Concentration <- as.numeric(long_data$Concentration)

my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")


setwd("~/PhD/plasma/vac69a/figures/")

# big figure with all analytes through time ####

all_analytes_sans_v9_plot <- ggplot(long_data, aes(x=factor(timepoint, levels=c("C-1", "DoD", "T+6", "C+45")), y=Concentration, color=Volunteer), group=Volunteer)+
  geom_point()+
  geom_line(aes(group=Volunteer))+
  facet_wrap(~ Analyte, scales = "free", ncol=8)+
  scale_y_log10()+
  theme_bw()+
  scale_color_manual(values=my_paired_palette)+
  theme(axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"))

ggsave(filename = "all_analytes_sans_v9_log_plot.png", all_analytes_sans_v9_plot, width = 16, height=9)

#viva_data <- filter(long_data, Analyte %in% c("CXCL10", "IL12p70", "IL10", "TNFRII", "IL1RA", "ICAM1"))

  
# figures of fold change ####  
  
fc_data <- long_data
#'fc_data$Sample_ID <- paste(fc_data$Volunteer, fc_data$Timepoint, sep="_")
fc_data <- pivot_wider(long_data, values_from = Concentration, names_from = timepoint)
colnames(fc_data)[3:6] <- c("Baseline", "DoD", "T6", "C45") 
fc_data$Baseline_DoD <- fc_data$DoD/fc_data$Baseline
fc_data$Baseline_T6 <- fc_data$T6/fc_data$Baseline
fc_data$Baseline_C45 <- fc_data$C45/fc_data$Baseline
fc_data <- subset(fc_data, !fc_data$Volunteer=="v009")

long_fc_data <- gather(fc_data, Comparison, Fold_Change, c("Baseline_DoD", "Baseline_T6"))

fc_data2 <-  gather(fc_data, Timepoint, Concentration, c("Baseline", "DoD", "T6", "C45"))

mean_fc <- fc_data2 %>%
  group_by(Analyte, Timepoint) %>%
  mutate("Mean_FC_DoD" =mean(Baseline_DoD)) %>%
  ungroup() %>%
  group_by(Analyte, Timepoint) %>%
  mutate("Mean_FC_T6"=mean(Baseline_T6)) %>%
  ungroup()

mean_fc2 <- select(mean_fc, Analyte, Mean_FC_DoD) 
mean_fc2 <- mean_fc2[!duplicated(mean_fc2),]

cluster_mat <- as.matrix(mean_fc2$Mean_FC_DoD)
rownames(cluster_mat) <-mean_fc2$Analyte

bas_dod_dist <- dist(cluster_mat, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
base_dod_hclust <- hclust(baseline_dist)


fc_levels <- rownames(cluster_mat)[order(cluster_mat, decreasing = F)]


#fold change at DOD and T6 split by volunteer #

fc_split_by_volunteer <- ggplot(long_fc_data, aes(x=Volunteer, y=factor(Analyte, levels=fc_levels)))+
  geom_tile(aes(fill=log2(Fold_Change)))+
  scale_fill_gradientn(name="log2FC",
                       values = scales::rescale(c(min(long_fc_data$Fold_Change), 0, max(long_fc_data$Fold_Change)), to=c(0,1)),
                       colors = c("#0859C6","black","#FFA500"))+
  facet_wrap(~Comparison)+
  theme_void()+
  theme(
    axis.text.x = element_text(angle=45, hjust = 1, vjust = 1),
    axis.text.y = element_text(),
    plot.title = element_text(hjust=0.5),
    legend.title = element_text(),
    legend.margin=margin(0,0,0,0),
    legend.position = "right",
    legend.box.margin=margin(0,0,0,0))

ggsave("fc_split_by_volunteer.png", fc_split_by_volunteer, height=6, width=6)


#mean fold change at DoD and T6 ##

mean_fc3 <- mean_fc %>%
  select(Analyte, Mean_FC_DoD, Mean_FC_T6) %>%
  gather(Comparison, Fold_Change, c( Mean_FC_DoD, Mean_FC_T6)) %>%
  mutate("log2_Fold_Change"=log2(Fold_Change))

mean_fc3 <- mean_fc3[!duplicated(mean_fc3),]


mean_fc_plot <- ggplot(mean_fc3, aes(x=Comparison, y=factor(Analyte, levels=fc_levels)))+
  geom_tile(aes(fill=log2_Fold_Change))+
  scale_fill_gradientn(name="log2FC",
                       values = scales::rescale(c(min(mean_fc3$log2_Fold_Change), 0, max(mean_fc3$log2_Fold_Change)), to=c(0,1)),
                       colors = c("#0859C6","black","#FFA500"))+
  theme_void()+
  theme(
    axis.text.x = element_text(angle=45, hjust = 1, vjust = 1),
    axis.text.y = element_text(),
    plot.title = element_text(hjust=0.5),
    legend.title = element_text(),
    legend.margin=margin(0,0,0,0),
    legend.position = "right",
    legend.box.margin=margin(0,0,0,0))

ggsave("mean_fc_plot.png", mean_fc_plot, height=6, width=6)











