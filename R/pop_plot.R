library(tidyr)
library(dplyr)
library(ggplot2)
library(colorspace)

# laptop
# data <- read.csv("C:/Users/Florian/PhD/cytof/vac69a/Vac69a_michalina_compensated_Exported_Stats\ 4.csv")
# iMac
data <- read.csv("/Users/s1249052/PhD/cytof/vac69a/Vac69a_michalina_compensated_Exported_Stats\ 4.csv")

str(data)
colnames(data) <-c("CD4+", "Vd2+", "CD8+", "MAIT", "Tregs", "DN", "Activated", "Gate", "Timepoint", "Volunteer") 
head(data)


long_data <- gather(data, Population, Percentage, colnames(data)[1:7])
long_data$Gatef <- factor(long_data$Gate, levels=c("All T Cells", "Activated T cells"))
my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")


pop_plot <- ggplot(long_data, aes(x=factor(Timepoint, levels=c("C-1", "C+8", "C+10", "C+12", "DoD", "T+6")), y=Percentage, group=Volunteer, fill=Volunteer))+
         geom_bar(stat="identity", position="dodge")+
         facet_grid(Population~Volunteer+Gatef, scales="free")+
         scale_fill_manual(values=my_paired_palette)+
         theme_bw()+
         xlab("Timepoint")+
         theme(legend.position = "none",
               axis.text.x = element_text(angle = 60, hjust = 1, size=12),
               strip.text.x = element_text(size = 10))


setwd("/Users/s1249052/PhD/cytof/better_gating/double_flowsoms/figures")

ggsave("pop_plot.png", pop_plot, width = 14, height=8)




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


(pie_plot <- ggplot(pie_data, aes(x=Activated/2, y=Percentage, fill=Population, width=Activated))+
  geom_bar(stat="identity", color="black")+
  coord_polar("y")+
  facet_wrap(~Volunteer, ncol=3)+
  theme_void()+
  theme(legend.title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
  )

vols <- c("Volunteer 02", "Volunteer 03", "Volunteer 05", "Volunteer 06", "Volunteer 07", "Volunteer 09")
freqs <- c(9.39, 19.89, 10.3, 18, 12.16, 15.1)





up_dod_dod6 <- dplyr::filter(long_abun_clusters, Comparison=="dod_dod6")
up_dod_dod6 <- dplyr::filter(up_dod_dod6, Fold_Change>1)
up_dod_dod6 <- dplyr::filter(up_dod_dod6, Timepoint=="post")

up_dod_dod6[,c(2,4)] <- NULL
up_dod_dod6$Count <- up_dod_dod6$Count*100
up_dod_dod6$SubSet <- c("CD4", "CD4", "CD8", "Vd2", "Vd2", "CD8", "MAIT", "CD4", "CD8")

up_dod_dod6 <- up_dod_dod6[order(up_dod_dod6$SubSet),]
#up_dod_dod6$ClusterID <- paste0("Cluster_", up_dod_dod6$ClusterID)

up_dod_dod6$ymin[1]<-0
up_dod_dod6$ymin[2:nrow(up_dod_dod6)] <- cumsum(up_dod_dod6$Count)[1:nrow(up_dod_dod6)-1]
up_dod_dod6$ymax <- up_dod_dod6$ymin + up_dod_dod6$Count

my_long_palette <- c(rev(sequential_hcl(5, "Heat")),qualitative_hcl(100, "Cold"), qualitative_hcl(150, "Dynamic"))

ggplot(up_dod_dod6)+
 
  geom_rect(aes(fill=Fold_Change, ymin=up_dod_dod6$ymin, ymax=up_dod_dod6$ymax, xmax=2.8, xmin=2, colour=Fold_Change))+
  geom_rect(aes(fill=as.numeric(factor(ClusterID))*100, ymin=ymin, ymax=ymax, xmax=6, xmin=3))+
  geom_rect(aes(fill=as.numeric(factor(SubSet))*500, ymin=ymin, ymax=ymax, xmax=9, xmin=6.2))+
  scale_color_gradientn(colours=rev(sequential_hcl(5, "Heat")))+
  scale_fill_gradientn(guide = FALSE, colors=my_long_palette, values=c(scales::rescale(seq(1,50), to=c(0,0.8)),scales::rescale(c(500,900,1000,1500,1800,2000,2500), to=c(0.8,1))))+
  
  geom_text(aes(x=4.5, y=(ymin+ymax)/2, label = paste("Cluster", ClusterID, sep="\n")),size=2.4,fontface="bold")+
  geom_text(aes(x=7.5, y=(ymin+ymax)/2, label = paste(SubSet)), size=2.7, fontface="bold")+
  
  guides(color=guide_colorbar(barwidth = 2))+
  labs(color = "Fold Change")+
  theme(aspect.ratio=1,
        legend.title = element_text(),
        legend.title.align = 0.5,
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())+
  xlim(c(0, 9))+
  coord_polar("y")



# 
# 
# compute_angle = function(perc){
#   angle = -1
#   # Or even more compact, but less readable
#   if(perc < 0.5) # 1st half [90, -90]
#     angle = (180 - (perc/0.5) * 180) - 90
#   else # 2nd half [90, -90]
#     angle = (90 - ((perc - 0.5)/0.5) * 180)
#   
#   return(angle)
# }






