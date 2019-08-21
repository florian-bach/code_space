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
