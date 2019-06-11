library(ggplot2)
library(colorspace)
library(gridExtra)

setwd("/Users/s1249052/PhD/cytof/better_gating/double_flowsoms/figures")


flist <- list.files(".", pattern="*.csv")

up_dod_dod6 <- data.frame()

for(i in flist){
  tmp <- read.csv(i)
  tmp$Volunteer <- substr(i, 1, 3)
  up_dod_dod6 <- rbind(up_dod_dod6, tmp)
}


#my_palette <- c(qualitative_hcl(nrow(up_dod_dod6), "Dynamic"), qualitative_hcl(5, "Dark3"))

#my_palette <- c("CD4", "CD8", "MAIT", "Vd2", "DN")

my_palette <- c(qualitative_hcl(64, "Dynamic"), qualitative_hcl(5, "Dark3"))
names(my_palette) <- seq(1,69)
names(my_palette)[65:69] <- c("CD4", "MAIT", "CD8", "Vd2", "DN")


for (i in unique(up_dod_dod6$Volunteer)){
  
  tmp_dat <- NULL
  tmp_dat <- dplyr::filter(up_dod_dod6, up_dod_dod6$Volunteer==i)
  
  tmp_dat[nrow(tmp_dat)+1, ] = list(0,0,0,0, "CD4", -1,-1, unique(tmp_dat$Volunteer))
  tmp_dat[nrow(tmp_dat)+1, ] <- list(0,0,0,0, "CD8", -1,-1, unique(tmp_dat$Volunteer))
  tmp_dat[nrow(tmp_dat)+1, ] <- list(0,0,0,0, "MAIT", -1,-1, unique(tmp_dat$Volunteer))
  tmp_dat[nrow(tmp_dat)+1, ] <- list(0,0,0,0, "Vd2", -1,-1, unique(tmp_dat$Volunteer))
  tmp_dat[nrow(tmp_dat)+1, ] <- list(0,0,0,0, "DN", -1,-1, unique(tmp_dat$Volunteer))
  
  tmp_dat$ymin <- tmp_dat$ymin/max(tmp_dat$ymax)+4
  tmp_dat$ymax <- tmp_dat$ymax/max(tmp_dat$ymax)+4
  
  assign(
    paste(unique(tmp_dat$Volunteer), "_pie_plot", sep=''), 
  
  ggplot()+
  # circles made of numbers indicating fold change
  geom_text(aes(x=7.05, y=seq(0,max(tmp_dat$ymax), by=max(tmp_dat$ymax)/nrow(tmp_dat)/8), label=paste0("4")), size=2)+
  geom_text(aes(x=8.05, y=seq(0,max(tmp_dat$ymax), by=max(tmp_dat$ymax)/nrow(tmp_dat)/8), label=paste0("8")), size=2)+
  geom_text(aes(x=9.05, y=seq(0,max(tmp_dat$ymax), by=max(tmp_dat$ymax)/nrow(tmp_dat)/8), label=paste0("16")), size=2)+
  geom_text(aes(x=10.05, y=seq(0,max(tmp_dat$ymax), by=max(tmp_dat$ymax)/nrow(tmp_dat)/8), label=paste0("32")), size=2)+

  # clusters and subset pies
  geom_rect(data=tmp_dat, aes_(fill=factor(tmp_dat$SubSet, levels=c("CD4", "CD8", "MAIT", "Vd2", "DN")), ymin=tmp_dat$ymin, ymax=tmp_dat$ymax, xmax=20, xmin=17))+
  geom_rect(data=tmp_dat, aes_(fill=factor(tmp_dat$ClusterID), ymin=tmp_dat$ymin, ymax=tmp_dat$ymax, xmax=10+log2(tmp_dat$Fold_Change), xmin=6+log2(tmp_dat$Fold_Change)))+
    
  #labelling those pies
  #geom_text(data=tmp_dat, aes(x=4.5+log2(tmp_dat$Fold_Change), y=(ymin+ymax)/2, label = paste("C_", ClusterID, sep="")),size=2.4,fontface="bold")+
  geom_text(data=tmp_dat, aes_(x=18.5, y=(tmp_dat$ymin+tmp_dat$ymax)/2, label = paste(tmp_dat$SubSet), size=1.5, fontface="bold"), show.legend = F)+
  #beautify
  scale_fill_manual(guide="none", values=my_palette)+
    theme_void()+
    theme(aspect.ratio=1,
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust=0.5, size = 20, face = "bold"))+
  coord_polar(theta = "y")+
  ggtitle(paste("Volunteer", substr(i, 2, 3)))+
  
  ylim(4,5)+
  xlim(4, 20)
  )
 }


#ggsave("v09_hot_pie.pdf", V09, height = 4, width = 4)
#grid.arrange(V02_pie_plot, V03_pie_plot, V05_pie_plot)
grid.arrange(V02_pie_plot, V03_pie_plot, V05_pie_plot, V06_pie_plot, V07_pie_plot, V09_pie_plot, ncol=2, nrow=3, widths=rep(5, 3), heights=rep(5, 2))

ggsave("all_v_hot_pie.png", grid.arrange(V02_pie_plot, V03_pie_plot, V05_pie_plot, V06_pie_plot, V07_pie_plot, V09_pie_plot, ncol=2, nrow=3, widths=rep(5, 2), heights=rep(5, 3)), width = 10, height=15)

