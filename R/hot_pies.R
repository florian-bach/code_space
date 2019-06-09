library(ggplot2)
library(colorspace)
library(gridExtra)

setwd("C:/Users/Florian/PhD/cytof/vac69a/double_flowsoms/figures")


flist <- list.files(".", pattern="*.csv")

up_dod_dod6 <- data.frame()

for(i in flist){
  tmp <- read.csv(i)
  up_dod_dod6 <- rbind(up_dod_dod6, tmp)
}

up_dod_dod6$Volunteer <- as.character(up_dod_dod6$Volunteer)

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
  geom_text(aes(x=5.05, y=seq(0,max(tmp_dat$ymax), by=max(tmp_dat$ymax)/nrow(tmp_dat)/8), label=paste0("4")), size=1.2)+
  geom_text(aes(x=6.05, y=seq(0,max(tmp_dat$ymax), by=max(tmp_dat$ymax)/nrow(tmp_dat)/8), label=paste0("8")), size=1.2)+
  geom_text(aes(x=7.05, y=seq(0,max(tmp_dat$ymax), by=max(tmp_dat$ymax)/nrow(tmp_dat)/8), label=paste0("16")), size=1.2)+
  geom_text(aes(x=8.05, y=seq(0,max(tmp_dat$ymax), by=max(tmp_dat$ymax)/nrow(tmp_dat)/8), label=paste0("32")), size=1.2)+

  # clusters and subset pies
  geom_rect(data=tmp_dat, aes_(fill=factor(tmp_dat$SubSet, levels=c("CD4", "CD8", "MAIT", "Vd2", "DN")), ymin=tmp_dat$ymin, ymax=tmp_dat$ymax, xmax=16, xmin=14))+
  geom_rect(data=tmp_dat, aes_(fill=factor(tmp_dat$ClusterID), ymin=tmp_dat$ymin, ymax=tmp_dat$ymax, xmax=6+log2(tmp_dat$Fold_Change), xmin=3+log2(tmp_dat$Fold_Change)))+
    
  #labelling those pies
  #geom_text(data=tmp_dat, aes(x=4.5+log2(tmp_dat$Fold_Change), y=(ymin+ymax)/2, label = paste("C_", ClusterID, sep="")),size=2.4,fontface="bold")+
  geom_text(data=tmp_dat, aes_(x=15, y=(tmp_dat$ymin+tmp_dat$ymax)/2, label = paste(tmp_dat$SubSet), size=2, fontface="bold"), show.legend = F)+
  #beautify
  scale_fill_manual(guide="none", values=my_palette)+
  #scale_x_continuous(limits = c(0, 18))+
  theme(aspect.ratio=1,
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())+
  coord_polar(theta = "y")+
  theme_void()+
  ylim(4,5)
  )
 }


# make 0-1 transform and then specify size with ymax somehow?


#ggsave("v09_hot_pie.pdf", V09, height = 4, width = 4)
#grid.arrange(V02_pie_plot, V03_pie_plot, V05_pie_plot)
grid.arrange(V02_pie_plot, V03_pie_plot, V05_pie_plot, V06_pie_plot, V07_pie_plot, V09_pie_plot, ncol=3, nrow=2, widths=rep(3.3, 3), heights=rep(4, 2))
ggsave("all_v_hot_pie.pdf", grid.arrange(V02_pie_plot, V03_pie_plot, V05_pie_plot, V06_pie_plot, V07_pie_plot, V09_pie_plot, ncol=3, nrow=2, widths=rep(3.3, 3), heights=rep(4, 2)), width = 10, height=8)

