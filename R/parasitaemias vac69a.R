library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(tidyr)
library(gtools)
library(directlabels)

setwd("/Users/s1249052/PhD/oxford/vac69")

data <- read.csv("vac69a_parasitaemia.csv", header=T)

# get rid of garbage timepoints that mess up graph

parasitaemias <- gather(data, Timepoint, Genomes, c("C.1", colnames(data)[12:35]) )
parasitaemias[,2:13] <- NULL

parasitaemias$Timepoint <- gsub("C.1", "C-1", parasitaemias$Timepoint) 
parasitaemias$Timepoint <- gsub("D5.5.1", "D5.5", parasitaemias$Timepoint) 
parasitaemias$Timepoint <- gsub("D", "C+", parasitaemias$Timepoint) 



parasitaemias$Volunteer <- paste("V", substr(as.character(parasitaemias$Volunteer), nchar(as.character(parasitaemias$Volunteer))-1, nchar(as.character(parasitaemias$Volunteer))), sep='')

my_palette <- c("#588C73", "#D53E4F", "#F2E394",   "#F2AE72",   "#D96459",   "#1A9CC7")
my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")



dods <- data.frame(vol=c('v02', 'v03', 'v05', 'v06', 'v07', 'v09'),
                   paras=c(4907, 8054, 11707, 6441, 9061, 10718),
                   dod=c('D15.5', 'D12.5', 'D15', 'D15', 'D15.5', 'D16.5'))
dods$inter <- match(dods$dod, unique(parasitaemias$Timepoint))



ggplot(data=parasitaemias[!is.na(parasitaemias$Genomes),], aes(x=factor(Timepoint, levels = unique(mixedsort(parasitaemias$Timepoint))), y=Genomes, group=factor(Volunteer)))+
  geom_point(aes(color=factor(Volunteer)), size=2.5)+
  geom_line(aes(color=factor(Volunteer)), size=2)+
  scale_color_manual(values=my_paired_palette)+
  theme_bw()+
  xlab("Day of Infection")+
  ylab("Genome Copies / mL")+
  scale_y_continuous(trans="log10", limits=c(1, 27000), breaks=c(10, 100, 1000, 10000, 25000))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20) ,
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1, size=18, color="black"),
        axis.title.x = element_text(size=24, color="black"),
        axis.title.y = element_text(size=24, color="black"),
        axis.text.y = element_text(size=20, color="black"))#+
  geom_segment(data=dods, aes(x=inter, xend=inter, y = 0, yend=paras), colour=my_paired_palette, position=position_jitter(), inherit.aes=F)#+
  geom_segment(data=dods, aes(x=0, xend=inter, y = paras, yend=paras), colour=my_paired_palette, inherit.aes=F)
  
              
 ggsave ("/Users/s1249052/PhD/oxford/vac69/parasitaemias_vac69.pdf", height = 8, width=10)
 ggsave ("/Users/s1249052/PhD/oxford/vac69/parasitaemias_vac69.png", height = 8, width=10)
 

ggplot(dods, aes(x=vol, y=paras, fill=vol))+
  geom_bar(stat='identity')+
  theme_bw()+
  scale_fill_manual(values=my_paired_palette)+
  ylab("Parasitaemia at DoD")+
  geom_text(aes(label=dod), position=position_dodge(width=0.9), vjust=-0.25, size=6)+
  geom_hline(yintercept=10000)+
  theme(legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=18, color="black"),
        axis.title.y = element_text(size=24, color="black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=20, color="black"))
 
ggsave ("/Users/s1249052/PhD/oxford/vac69/dod_parasitaemias_vac69.pdf", height = 8, width=10)