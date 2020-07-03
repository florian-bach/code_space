library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbiplot)
library(plotly)
library(shadowtext)

#library(rgl) # want some 3d?

#data <- read.csv("/Users/s1249052/PhD/plasma/vac69a/Vivax_plasma_analytes2_no_inequalities.csv")

plasma <- read.csv("~/PhD/plasma/vac69a/Vivax_plasma_analytes2_no_inequalities.csv", header=T, stringsAsFactors = F)
plasma <- filter(plasma, Volunteer!="v009")

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
  facet_wrap(~ Analyte, scales = "free", ncol=6)+
  scale_y_log10()+
  theme_bw()+
  scale_color_manual(values=my_paired_palette)+
  theme(axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"))

  ggsave(filename = "all_analytes_sans_v9_log_plot.png", all_analytes_sans_v9_plot, width = 8, height=9)

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



####     pca stuff ####


data2 <- spread(long_data, Analyte, Concentration)
split_data <- split(data2, data2$Volunteer)

my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")
names(my_paired_palette) <- names(list_of_pcas)

data3 <- data2
data3[,3:ncol(data3)] <- log10(data3[,3:ncol(data3)])

big_pca <-  prcomp(data3[,3:ncol(data3)], center = T)
big_pca2 <- cbind(data3[, 1:2], big_pca$x)



all_vols_together_vol_color <- ggplot(big_pca2, aes(x=PC1, y=PC2))+
  geom_shadowtext(aes(label=timepoint, color=Volunteer), size=5, fontface="bold")+
  scale_color_manual(values=my_paired_palette)+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size=10),
        panel.border = element_rect(color="black", fill=NA),
        axis.title = element_text(size=12),
        #axis.text = element_blank(),
        
        plot.title = element_text(size=14, hjust=0.5)
  )

ggsave("all_vols_together_vol_color.png", all_vols_together_vol_color)



arrow_pca <- subset(big_pca2, big_pca2$timepoint %in% c("C-1", "DoD"))

wide_arrow_data <- arrow_pca[, 1:4]
wide_arrow_data <- pivot_wider(wide_arrow_data, names_from = timepoint, values_from = c(PC1, PC2))



arrow_pca_plot <- ggplot(arrow_pca, aes(x=PC1, y=PC2, group=Volunteer))+
  geom_point(aes(color=Volunteer))+
  geom_segment(data=wide_arrow_data, aes(x=`PC1_C-1`,xend=PC1_DoD, y=`PC2_C-1`, yend=PC2_DoD, color=Volunteer), arrow =arrow(length = unit(0.2, "cm")) )+
  #geom_line(arrow = arrow(length = unit(0.2, "cm")))+
  scale_color_manual(values=my_paired_palette)+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  theme(
        axis.text = element_text(size=10),
        panel.border = element_rect(color="black", fill=NA),
        axis.title = element_text(size=12),
        #axis.text = element_blank(),
        
        plot.title = element_text(size=14, hjust=0.5)
  )

ggsave("arrow_pca.png", arrow_pca_plot)





  list_of_pcas2 <- lapply(split_data, function(x){
  pca <- prcomp(x[,3:ncol(x)], center = T)
  cbind(x[,1:2], pca$x)
})


list_of_perc2 <- lapply(split_data, function(x){
  pca <- prcomp(x[,3:ncol(x)], center = T)
})

top_hits2 <- lapply(list_of_perc2, function(x){
  head(
    x$rotation[order((x$rotation[,1]), decreasing=T),],
    n=10)
})





top_hits <- lapply(list_of_perc, function(x){
  head(
    x$rotation[order(x$rotation[,1], decreasing=T),],
    n=15)
})

top_hits

  head(
    list_of_perc[1]$v002$rotation[order
      (abs(list_of_perc[1]$v002$rotation[,1]), decreasing=T),],
    n=10)









geom_text(aes(x=p$x[,1], y=p$x[,2], color=data$Volunteer, label=data$timepoint), size=7, fontface="bold")
# ggplot()+
#   +   geom_text(aes(x=p$x[,1], y=p$x[,2], color=data$Volunteer, label=data$timepoint), size=7, fontface="bold")+
#   +   xlab(paste("PC1 ", summary(p)$importance[2,1]*100, "%", sep = ""))+
#   +   ylab(paste("PC2 ", summary(p)$importance[2,2]*100, "%", sep = ""))+
#   +   scale_shape_manual(values=c(16:19))+
#   +   #scale_color_manual(values=my_paired_palette)+
#   +   theme_minimal()+
#   +   theme(legend.title = element_blank(),
#             +         axis.text = element_text(size=20),
#             +         axis.title = element_text(size=25))
# 
# 




plasma_pca <- prcomp(data2[,3:ncol(data2)], center = T, scale. = T)



summary(plasma_pca)  

# View(plasma_pca$rotation)

look <- data.frame(Analyte=rownames(plasma_pca$x), plasma_pca$x)

pca12 <- ggbiplot(plasma_pca, circle = T, choices=1:2, obs.scale = 1, var.scale = 1,var.axes=FALSE, groups = data$Volunteer, labels = data$timepoint)+
  theme_minimal()+

pca12


pca23 <- ggbiplot(plasma_pca, choices=2:3, obs.scale = 1, var.scale = 1,var.axes=FALSE, groups = data$Volunteer, labels = data$timepoint)+
  theme(legend.title = element_blank())
pca23
#plotly::ggplotly(pca, tooltip="timepoint")

rgl::plot3d(x=plasma_pca$x[,1], y=plasma_pca$x[,2], z=plasma_pca$x[,3], xlab=NULL, ylab=NULL, zlab=NULL,col=as.numeric(as.factor(data$Volunteer)), main=F, sub=F type='s', box=F, axes=F)
if (!rgl.useNULL())
  play3d(spin3d(axis = c(1, 0, 1), rpm = 2000), duration = 100)



rgl::spin3d(d3d, axis = c(0, 0, 1), rpm=10,dev = rgl.cur(), subscene = par3d("listeners", dev = dev))









#####     handmade

pca_viva <-ggplot()+
  # geom_point(aes(x=plasma_pca$x[,1], y=plasma_pca$x[,2], color=data$Volunteer), size=4)+
  # ggrepel::geom_label_repel(aes(x=plasma_pca$x[,1], y=plasma_pca$x[,2], color=data$Volunteer, label=data$timepoint), point.padding = 0.2, size=4)+
  geom_text(aes(x=plasma_pca$x[,1], y=plasma_pca$x[,2], color=data$Volunteer, label=data$timepoint), size=7, fontface="bold")+
  # geom_label(aes(x=plasma_pca$x[,1], y=plasma_pca$x[,2], color=data$Volunteer, label=data$timepoint), size=4)+
  xlab(paste("PC1 ", summary(plasma_pca)$importance[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", summary(plasma_pca)$importance[2,2]*100, "%", sep = ""))+
  scale_shape_manual(values=c(16:19))+
  scale_color_manual(values=my_paired_palette)+
  theme_minimal()+
  theme(legend.title = element_blank(),
        axis.text = element_text(size=20),
        axis.title = element_text(size=25))
   
ggsave("/Users/s1249052/PhD/plasma/vac69a/pca_viva.png", pca_viva, width=12, height = 10)

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









##############################################


# yeast1 <- data.frame(spline(x=filter(fiss$strain=="wt")$minute, y=fiss$count, n=100))
# yeast2 <- data.frame(spline(x=fiss$minute, y=fiss$count, n=100))
# 
# yeasts <- list(yeast1, yeast2)
# 

continuous_data <- long_data

continuous_data$timepoint <- as.numeric(factor(long_data$timepoint, levels=c("C-1", "DoD", "T+6","C+45")))*100-100


conti_list <- split(continuous_data, continuous_data$Volunteer)
yeet <- lapply(conti_list, FUN=function(x){split(x, x$Analyte)})
yeet[[1]] <- NULL
yeetus <- yeet[[1]]

applied_yeet <- lapply(yeetus, FUN = function(x){spline(x$timepoint, x$Concentration, n=100)})  

# works
# lapply(conti_list,
#        function(x){print(
#          c(x$timepoint[1:3], x$Concentration[1:3])
#          )
#          }
#        )
#

v2 <- filter(continuous_data, Volunteer=="v002")


ggplotly(ggplot()+
  geom_point(data=v2, aes(x = timepoint, y = Concentration, color= Analyte, group = Analyte))+
  #facet_wrap(~Analyte, scales="free")+
  #scale_y_continuous(limits=c(1,10000))+
  #geom_line(data=hello, aes(x=x, y=y))
  #geom_line(data=dplyr::bind_rows(yeasts, .id="df"), aes(x=x, y=y, color=df))
  geom_line(data=dplyr::bind_rows(applied_yeet, .id="df"), aes(x=x, y=y, color=df))+
  #facet_wrap(~df, scales="free")+
  #scale_y_log10()+
  theme_minimal()+
  theme(legend.position = "none"))


fiss

# garbage code you'll probably never need again ####


# data2 <- spread(long_data, Analyte, Concentration)
# split_data <- split(data2, data2$Volunteer)
# 
# my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")
# names(my_paired_palette) <- names(list_of_pcas)
# 
# 
# transformed_list_of_analytes <- lapply(
#   split_data, function(x) cbind(x[,1:2],log10(x[,3:ncol(x)])))
# 
# 
# list_of_pcas <- lapply(transformed_list_of_analytes, function(x){
#   pca <- prcomp(x[,3:ncol(x)], center = T)
#   cbind(x[,1:2], pca$x)
# })
# 
# 
# list_of_perc <- lapply(transformed_list_of_analytes, function(x){
#   pca <- prcomp(x[,3:ncol(x)], center = T)
# })
# 
# time_col <- colorspace::sequential_hcl(5, palette = "Purple Yellow")
# time_col_scheme <- c("C-1"=time_col[4], "C+45"=time_col[5], "DoD"=time_col[2], "T+6"=time_col[1])
# 
# for(i in 1:length(list_of_pcas)){
#   
#   p <- list_of_pcas[[i]]
#   q <-list_of_perc[[i]]
#   
#   assign(paste(names(list_of_pcas[i]), "_pca_plot", sep=''),
#          
#          ggplot(p, aes(x=PC1, y=PC2, colour=Volunteer))+
#            scale_color_manual(values=my_paired_palette)+
#            geom_shadowtext(aes(label=timepoint, color=timepoint), size=5, fontface="bold")+
#            ggtitle(paste(names(list_of_pcas[i])))+
#            scale_color_manual(values=time_col_scheme)+
#            xlab(paste("PC1 ", data.frame(summary(q)[6])[2,1]*100, "%", sep = ""))+
#            ylab(paste("PC2 ", data.frame(summary(q)[6])[2,2]*100, "%", sep = ""))+
#            theme_minimal()+
#            theme(legend.position = "none",
#                  axis.text = element_text(size=10),
#                  axis.title = element_text(size=12),
#                  #axis.text = element_blank(),
#                  plot.title = element_text(size=14, hjust=0.5)
#            )
#   )
# }
# 
# 
# individual_pcas <- gridExtra::grid.arrange(v002_pca_plot, v003_pca_plot, v005_pca_plot, v006_pca_plot, v007_pca_plot)
# 
# ggsave("individual_pcas.png", individual_pcas, width=10, height=10)
# 
# all_vols_together <- do.call(rbind, list_of_pcas)
# all_perc_together <- do.call(rbind, list_of_perc)
# 
# all_vols_together_individual_color <- ggplot(all_vols_together, aes(x=PC1, y=PC2, colour=Volunteer))+
#   scale_color_manual(values=my_paired_palette)+
#   geom_shadowtext(aes(label=timepoint, color=Volunteer), size=5, fontface="bold")+
#   #scale_color_manual(values=time_col_scheme)+
#   xlab(paste("PC1 ", data.frame(summary(q)[6])[2,1]*100, "%", sep = ""))+
#   ylab(paste("PC2 ", data.frame(summary(q)[6])[2,2]*100, "%", sep = ""))+
#   theme_minimal()+
#   theme(legend.position = "none",
#         axis.text = element_text(size=10),
#         axis.title = element_text(size=12),
#         #axis.text = element_blank(),
#         plot.title = element_text(size=14, hjust=0.5)
#   )
# 
# all_vols_together_time_color <- ggplot(all_vols_together, aes(x=PC1, y=PC2))+
#   geom_shadowtext(aes(label=timepoint, color=timepoint), size=5, fontface="bold")+
#   scale_color_manual(values=time_col_scheme)+
#   xlab(paste("PC1 ", data.frame(summary(q)[6])[2,1]*100, "%", sep = ""))+
#   ylab(paste("PC2 ", data.frame(summary(q)[6])[2,2]*100, "%", sep = ""))+
#   theme_minimal()+
#   theme(legend.position = "none",
#         axis.text = element_text(size=10),
#         panel.border = element_rect(color="black", fill=NA),
#         axis.title = element_text(size=12),
#         #axis.text = element_blank(),
#         
#         plot.title = element_text(size=14, hjust=0.5)
#   )
# 
# 
# plot_grid(all_vols_together_time_color, all_vols_together_individual_color, ncol=1)




