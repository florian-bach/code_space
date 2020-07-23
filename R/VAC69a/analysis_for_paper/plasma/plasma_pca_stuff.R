####     pca stuff ####

library(tidyr)
library(dplyr)
library(ggplot2)

setwd("~/PhD/plasma/vac69a/")

long_data <- read.csv("long_plasma_data.csv", header = TRUE, stringsAsFactors = FALSE)
sig_analytes <- read.delim("analytes_sorted_by_padj.txt", header = T)

#filter for top 20 analytes
long_data <- subset(long_data, long_data$Analyte %in% sig_analytes$Analyte[1:12])

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
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,3]*100, "%", sep = ""))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size=10),
        panel.border = element_rect(color="black", fill=NA),
        axis.title = element_text(size=12),
        #axis.text = element_blank(),
        
        plot.title = element_text(size=14, hjust=0.5)
  )

ggsave("./figures/all_vols_together_top12_analytes.png", all_vols_together_vol_color)



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

ggsave("./figures/arrow_pca_n12.png", arrow_pca_plot)




pc_one <- rownames(head(abs(big_pca$rotation[order(abs(big_pca$rotation[,1]), decreasing = T),1:3]), n=20))
pc_one_data <- subset(long_data, long_data$Analyte %in% pc_one)
pc_one_data$AnalyteF <- factor(pc_one_data$Analyte, levels=pc_one)

pc_two <- rownames(head(abs(big_pca$rotation[order(abs(big_pca$rotation[,2]), decreasing = T),1:3]), n=20))
pc_two_data <- subset(long_data, long_data$Analyte %in% pc_two)
pc_two_data$AnalyteF <- factor(pc_two_data$Analyte, levels=pc_two)



pc_one_plot <- ggplot(pc_one_data, aes(x=factor(timepoint, levels=c("C-1", "DoD", "T+6", "C+45")), y=Concentration, color=Volunteer), group=Volunteer)+
  geom_point()+
  geom_line(aes(group=Volunteer))+
  facet_wrap(~ AnalyteF, scales = "free", ncol=5)+
  
  theme_bw()+
  scale_color_manual(values=my_paired_palette)+
  theme(axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"))

ggsave(filename = "./figures/pc_one_deconvoluted.png", pc_one_plot, width = 16, height=9)



pc_two_plot <- ggplot(pc_two_data, aes(x=factor(timepoint, levels=c("C-1", "DoD", "T+6", "C+45")), y=Concentration, color=Volunteer), group=Volunteer)+
  geom_point()+
  geom_line(aes(group=Volunteer))+
  facet_wrap(~ AnalyteF, scales = "free", ncol=5)+
  scale_y_log10()+
  theme_bw()+
  scale_color_manual(values=my_paired_palette)+
  theme(axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"))

ggsave(filename = "./figures/pc_two_deconvoluted.png", pc_two_plot, width = 16, height=9)







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
