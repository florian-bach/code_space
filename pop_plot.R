library(tidyr)
library(dplyr)
library(ggplot2)

data <- read.csv("~/PhD/cytof/vac69a/Vac69a_michalina_compensated_Exported_Stats 4.csv")
str(data)
colnames(data) <-c("CD4+", "Vd2+", "CD8+", "MAIT", "Tregs", "DN", "Activated", "Gate", "Timepoint", "Volunteer") 
head(data)


long_data <- gather(data, Population, Percentage, colnames(data)[1:7])
long_data$Gatef <- factor(long_data$Gate, levels=c("All T Cells", "Activated T cells"))
my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")


pie_plot <- ggplot(long_data, aes(x=factor(Timepoint, levels=c("C-1", "C+8", "C+10", "C+12", "DoD", "T+6")), y=Percentage, group=Volunteer, fill=Volunteer))+
         geom_bar(stat="identity", position="dodge")+
         facet_grid(Population~Volunteer+Gatef, scales="free")+
         scale_fill_manual(values=my_paired_palette)+
         theme_bw()+
         xlab("Timepoint")+
         theme(legend.position = "none",
               axis.text.x = element_text(angle = 60, hjust = 1, size=12))


setwd("/Users/s1249052/PhD/cytof/better_gating/double_flowsoms/figures")

ggsave("pop_plot.pdf", pop_plot)




##############       pie charts yoooooo


norm_data <- data
norm_data$Activated <- rep(norm_data$Activated[1:30], times=2)
norm_data$sum <- apply(norm_data[,1:6], 1, sum)
norm_data$sum <- norm_data$sum/100

norm_data[,1:6] <- norm_data[,1:6]/norm_data$sum
norm_data$sum <- NULL

long_norm_data <- gather(norm_data, Population, Percentage, colnames(norm_data)[1:6])

pie_data <- filter(long_norm_data, Timepoint=="T+6")
pie_data <- filter(pie_data, Population != "Activated")
pie_data <- filter(pie_data, Gate=="Activated T cells")


(pie_plot <- ggplot(pie_data, aes(x=Activated*2, y=Percentage, fill=Population, width=Activated))+
  geom_bar(stat="identity", color="black")+
  coord_polar("y")+
  facet_grid(~Volunteer)+
  theme_void()+
  theme(legend.title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
  )

vols <- c("Volunteer 02", "Volunteer 03", "Volunteer 05", "Volunteer 06", "Volunteer 07", "Volunteer 09")
freqs <- c(9.39, 19.89, 10.3, 18, 12.16, 15.1)





up_dod_dod6 <- filter(long_abun_clusters, Comparison=="dod_dod6")
up_dod_dod6 <- filter(up_dod_dod6, Fold_Change>1)
up_dod_dod6 <- filter(up_dod_dod6, Timepoint=="post")

up_dod_dod6[,c(2,4)] <- NULL
up_dod_dod6$Count <- up_dod_dod6$Count*100
up_dod_dod6$SubSet <- c("CD4", "CD4", "CD8", "Vd2", "Vd2", "CD8", "MAIT", "CD4", "CD8")

#up_dod_dod6 <- data.frame(node=up_dod_dod6$parent, parent=up_dod_dod6$ClusterID, size=up_dod_dod6$Count, colour=paste("Cluster_", up_dod_dod6$ClusterID))
up_dod_dod6 <- up_dod_dod6[order(up_dod_dod6$SubSet),]
up_dod_dod6$ClusterID <- paste("Cluster_", up_dod_dod6$ClusterID, sep='')




up_dod_dod6$ymin[2:nrow(up_dod_dod6)] <- cumsum(up_dod_dod6$Count)[1:nrow(up_dod_dod6)-1]
up_dod_dod6$ymax <- up_dod_dod6$ymin + up_dod_dod6$Count

my_palette <- colorRampPalette(c("#ffd4c4", "#ffb5e2", "#ce70e8", "#a347e5", "#6a39ff"))(n=20)


ggplot(up_dod_dod6)+
  
  #geom_rect(aes(fill=Fold_Change, ymin=ymin, ymax=ymax, xmax=2.8, xmin=2), inherit.aes = F)+
  geom_rect(aes(fill=SubSet, ymin=ymin, ymax=ymax, xmax=7, xmin=6), inherit.aes = F)+
   geom_rect(aes(fill=ClusterID, ymin=ymin, ymax=ymax, xmax=5, xmin=3), inherit.aes = F)+
   geom_text(aes(x=4.1, y=(ymin+ymax)/2, label = paste(ClusterID)),size=2.4, inherit.aes = F)+
   geom_text(aes(x=6.5, y=(ymin+ymax)/2, label = paste(SubSet)), size=2.7, inherit.aes = F)+
  #scale_fill_gradient(low = "black", high = "red")+
  theme(aspect.ratio=1,
        axis.text = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())+
  xlim(c(0, 7))+
  coord_polar("y")

  rev(c("#ffd4c4", "#ffb5e2", "#ce70e8", "#a347e5"))
