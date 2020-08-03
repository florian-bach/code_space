library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggthemes)
library(cowplot)
library(plyr)
library(gtools)


setwd("C:/users/FLorian/PhD/clinical data/vac68/")
clinic<-read.csv("clinvac68y.csv", header=T)

clinic$Individual<-ifelse(clinic$trial_number==6801004, "04", "08")
clinic$trial_number<-NULL
clinic$review_ae_unsol<-NULL
clinic$pyrexia_temp<-NULL

clinic4<-subset(clinic, clinic$Individual=="04")
clinic8<-subset(clinic, clinic$Individual=="08")



clinic4<-melt(clinic4, id=c("Individual", "timepoint"))
colnames(clinic4)<-c("Individual", "Timepoint", "Symptom", "Severity")

clinic8<-melt(clinic8, id=c("Individual", "timepoint"))
colnames(clinic8)<-c("Individual", "Timepoint", "Symptom", "Severity")

clinic4[is.na(clinic4)] <- 0
clinic8[is.na(clinic8)] <- 0

clinic4$Severity<-as.factor(clinic4$Severity)
clinic8$Severity<-as.factor(clinic8$Severity)

clinic4<-clinic4[order(clinic4$Timepoint),]
clinic8<-clinic8[order(clinic8$Timepoint),]

clinic4$Timepoint<-mixedsort(clinic4$Timepoint)
clinic8$Timepoint<-mixedsort(clinic8$Timepoint)

c4plot<-ggplot(data=clinic4, aes(x=clinic4$Timepoint, y=clinic4$Symptom))+
  geom_tile(aes(fill=Severity), color="white")+
  ylab("Symptom")+
  xlab("Timepoint")+
  scale_fill_manual(values= c("grey", "yellow", "orange", "red"))+
  ggtitle("Volunteer 04")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.title = element_text(size = 20, face= "bold"))+
  theme(axis.text=element_text(size=16), plot.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"))


c4plot

c8plot<-ggplot(data=clinic8, aes(x=clinic8$Timepoint, y=clinic8$Symptom))+
  geom_tile(aes(fill=Severity), color="white")+
  geom_tile(aes(fill=Severity), color="white")+
  ylab("Symptom")+
  xlab("Timepoint")+
  scale_fill_manual(values= c("grey", "yellow", "orange", "red"))+
  ggtitle("Volunteer 08")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.title = element_text(size = 20, face= "bold"))+
  theme(axis.text=element_text(size=16), plot.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"))



vac68clin<-plot_grid(c4plot, c8plot, align = "hv", cols=1, labels = "AUTO", label_size = 19)
ggsave("clinvac68.pdf", height=12, width=12)








data <- read.csv("~/PhD/clinical_data/vac69a/symptoms_vac69a.csv", header = T, stringsAsFactors = F)


long_data <- tidyr::gather(data, Symptom, Severity, colnames(data)[c(11, 13:ncol(data))])

long_data$timepoint <- substr(as.character(long_data$timepoint), stringr::str_locate(pattern="_", string=as.character(long_data$timepoint)[1])+1, nchar(as.character(long_data$timepoint)))
long_data$timepoint <- ifelse(substr(long_data$timepoint, 1, 1)=="_", substr(long_data$timepoint, 2, nchar(long_data$timepoint)), long_data$timepoint)

long_data$timepoint <- gsub("chall__", "C", long_data$timepoint)
long_data$timepoint <- gsub("___", "", long_data$timepoint)
long_data$timepoint <- gsub("_or_c_28", "", long_data$timepoint)
long_data$timepoint <- gsub("postchall_ep", "T", long_data$timepoint)
long_data$timepoint <- gsub("_am", "", long_data$timepoint)
long_data$timepoint <- gsub("_pm", "_5", long_data$timepoint)
long_data$timepoint <- gsub("diagnosis_or_c_21", "_5", long_data$timepoint)

timepoint_levels <- unique(long_data$timepoint[gtools::mixedorder(long_data$timepoint)])

symptom_heatmap <- ggplot(long_data, aes(x=factor(timepoint, levels=timepoint_levels), y=Symptom))+
  geom_tile(aes(fill=factor(Severity), width=0.92, height=0.92), color=ifelse(grepl("diagnosis", long_data$timepoint), "black", "grey"))+
  scale_fill_manual(values =  list("lightgrey", "yellow", "orange", "red"))+
  facet_wrap(~trial_number, scales="free")+
  theme_minimal()+
  guides(fill=guide_legend(title="Severity"))+
  theme(axis.text.x = element_text(hjust=1, angle=60),
        axis.title = element_blank(),
        legend.title = element_blank())

ggsave("/home/flobuntu/PhD/clinical_data/vac69a/figures/symptom_heatmap.png", symptom_heatmap, height=8, width=14)


  # 
  # first
  # 
  # combo_figure <- cowplot::plot_grid(first, second, ncol=1)

