library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbiplot)
library(rgl)

data <- read.csv("/Users/s1249052/PhD/plasma/vac69a/vac69a_Vivax_plasma_analytes_no_inequalities.csv")
data<- na.omit(data)
data$timepoint <- gsub("D", "DoD", data$timepoint, fixed=T)
data$timepoint <- gsub("DoDoDoD+6", "T+6", data$timepoint, fixed=T)
data$timepoint <- gsub("DoD+6", "T+6", data$timepoint, fixed=T)

long_data <- gather(data, Analyte, Concentration, colnames(data)[3:37])

long_data$Analyte <- substr(long_data$Analyte, 1, nchar(long_data$Analyte)-7)

long_data$Analyte <- gsub(".", "", long_data$Analyte, fixed = T)
long_data$Concentration <- as.numeric(long_data$Concentration)



ggplot(long_data, aes(x=factor(long_data$timepoint, levels=c("C-1", "DoD", "T+6", "C+45")), y=Concentration, color=Volunteer), group=Volunteer)+
  geom_point()+
  geom_line(aes(group=long_data$Volunteer))+
  facet_wrap(~ Analyte, scales = "free")


####     pca plot

data2 <- spread(long_data, Analyte, Concentration)


plasma_pca <- prcomp(data2[,3:37], center = T, scale. = T)

summary(plasma_pca)  

pca12 <- ggbiplot(plasma_pca, circle = T, choices=1:2, obs.scale = 1, var.scale = 1,var.axes=FALSE, groups = data$Volunteer, labels = data$timepoint)+
  theme_minimal()+

pca12


pca23 <- ggbiplot(plasma_pca, choices=2:3, obs.scale = 1, var.scale = 1,var.axes=FALSE, groups = data$Volunteer, labels = data$timepoint)+
  theme(legend.title = element_blank())
pca23
#plotly::ggplotly(pca, tooltip="timepoint")

rgl::plot3d(plasma_pca$x[,1], plasma_pca$x[,2], plasma_pca$x[,3], col=as.numeric(as.factor(data$timepoint)))

my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")


#####     handmade

ggplot()+
  # geom_point(aes(x=plasma_pca$x[,1], y=plasma_pca$x[,2], color=data$Volunteer), size=4)+
  # ggrepel::geom_label_repel(aes(x=plasma_pca$x[,1], y=plasma_pca$x[,2], color=data$Volunteer, label=data$timepoint), point.padding = 0.2, size=4)+
  geom_text(aes(x=plasma_pca$x[,1], y=plasma_pca$x[,2], color=data$Volunteer, label=data$timepoint), size=5, fontface="bold")+
  # geom_label(aes(x=plasma_pca$x[,1], y=plasma_pca$x[,2], color=data$Volunteer, label=data$timepoint), size=4)+
  xlab(paste("PC1 ", summary(plasma_pca)$importance[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", summary(plasma_pca)$importance[2,2]*100, "%", sep = ""))+
  scale_shape_manual(values=c(16:19))+
  scale_color_manual(values=my_paired_palette)+
  theme_minimal()+
  theme(legend.title = element_blank(),
        axis.text = element_text(size=14),
        axis.title = element_text(size=17))
   
###########




#cytof test


cytof <- read.csv("/Users/s1249052/PhD/cytof/better_gating/big_flowsoms/FlowSOM_all_cd4s_baseline_dod_t6_results/results/cluster_abundances.csv")
colnames(cytof)[3:20] <- paste(rep(c("V06", "V07", "V09", "V02", "V03", "V05"), each=3), rep(c("Baseline_01", "DoD_01", "DoD+6_01"), times=6), sep='_')

long_cytof <- gather(cytof, Dataset, Frequency, colnames(cytof[3:20]))

long_cytof$Volunteer <- substr(long_cytof$Datase, 1, 3)
long_cytof$Timepoint <- substr(long_cytof$Dataset, 5, nchar(long_cytof$Dataset)-3)
long_cytof$ClusterID <- paste("C_", long_cytof$ClusterID, sep='')
long_cytof <- select(long_cytof, -MetaclusterID, -Dataset)

broad_cytof <- spread(long_cytof, ClusterID, Frequency)

cytof_pca <- prcomp(broad_cytof[,3:102], center = T, scale. = T)










cytof <- read.csv("/Users/s1249052/PhD/cytof/vac63c/analysis/FlowSOM_all_cd4+_2_results/results/cluster_abundances.csv")

colnames(cytof)[3:49] <- substr(colnames(cytof)[3:49], nchar(colnames(cytof)[3:49])-10, nchar(colnames(cytof)[3:49])-4)

# get rid of control files
cytof <- cytof[,-grep("tr", colnames(cytof), fixed = T)]
colnames(cytof) <- gsub("_", "", colnames(cytof), fixed=T)
colnames(cytof) <- gsub(".", "", colnames(cytof), fixed=T)


long_cytof <- gather(cytof, Dataset, Frequency, colnames(cytof[3:46]))

long_cytof$Volunteer <- substr(long_cytof$Datase, 1, 3)
long_cytof$Timepoint <- substr(long_cytof$Dataset, 4, nchar(long_cytof$Dataset))
long_cytof$ClusterID <- paste("C_", long_cytof$ClusterID, sep='')
long_cytof <- select(long_cytof, -MetaclusterID, -Dataset)

broad_cytof <- spread(long_cytof, ClusterID, Frequency)

cytof_pca <- prcomp(broad_cytof[,3:102], center = T, scale. = T)
















ggplot()+
   # geom_point(aes(x=cytof_pca$x[,1], y=cytof_pca$x[,2]), size=4)+
  # ggrepel::geom_label_repel(aes(x=cytof_pca$x[,1], y=cytof_pca$x[,2], color=data$Volunteer, label=data$timepoint), point.padding = 0.2, size=4)+
  #geom_text(aes(x=cytof_pca$x[,1], y=cytof_pca$x[,2], color=broad_cytof$Volunteer, label=broad_cytof$Timepoint), size=5, fontface="bold")+
  geom_line(aes(x=cytof_pca$x[,1], y=cytof_pca$x[,2], group=broad_cytof$Volunteer, color=broad_cytof$Volunteer))+
  # geom_label(aes(x=cytof_pca$x[,1], y=cytof_pca$x[,2], color=data$Volunteer, label=data$timepoint), size=4)+
  xlab(paste("PC1 ", summary(cytof_pca)$importance[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", summary(cytof_pca)$importance[2,2]*100, "%", sep = ""))+
  scale_shape_manual(values=c(16:19))+
  #scale_color_manual(values=my_paired_palette)+
  theme_minimal()+
  theme(legend.title = element_blank(),
        axis.text = element_text(size=14),
        axis.title = element_text(size=17))



rgl::plot3d(cytof_pca$x[,1], cytof_pca$x[,2], cytof_pca$x[,3], col=as.numeric(as.factor(broad_cytof$Volunteer)))
