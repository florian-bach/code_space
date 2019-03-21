library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggthemes)
library(cowplot)
library(plyr)
library(gtools)


setwd("/Users/s1249052/PhD/clinical data")
clinic<-read.csv("Copy of VAC068 2018-07-12-1.csv", header=T)

clinic$Individual<-ifelse(clinic$trial_number==6801004, "04", "08")
clinic$trial_number<-NULL
clinic$review_ae_unsol<-NULL
clinic$pyrexia_temp<-NULL

symptoms<-clinic[,colnames(clinic)[9:24]]
symptoms$timepoint<-clinic$timepoint

symptoms<-melt(symptoms, id=c("Individual", "timepoint"))
colnames(symptoms)<-c("Individual", "Timepoint", "Symptom", "Severity")

#symptoms[is.na(symptoms)] <- 0
symptoms$Severity<-as.factor(symptoms$Severity)
symptoms<-symptoms[order(symptoms$Timepoint),]
symptoms$Timepoint<-substring(symptoms$Timepoint, 3,8)

symptoms$Timepoint<-factor(symptoms$Timepoint, levels=mixedsort(symptoms$Timepoint))

symptoms$Timepoint

# "10am" "10pm" "11am" "11pm" "12am" "12pm" "13am" "13pm" "14am" "14pm" "6pm" 
# [12] "7am"  "7pm"  "8am"  "8pm"  "9am"  "9pm"  " 1"   " 11"  " 13"  " 15"  " 2"  
# [23] " 3"   " 5"   " 7"   " 9"  

scale_x_discrete(limits=c()

ggplot(symptoms, aes(x=mixedsort(symptoms$Timepoint), y=Symptom, fill=Severity))+
  geom_tile(stat="identity")+
  facet_wrap(~ Individual)+
  ylab("Symptom")+
  xlab("Timepoint")+
  scale_fill_manual(values= c("grey", "yellow", "orange", "red"))+
  #ggtitle("Volunteer 04")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.title = element_text(size = 20, face= "bold"))+
  theme(axis.text=element_text(size=16), plot.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"))


c4plot

ggplot(data=subset(symptoms, Individual=="04"), aes(x=Timepoint, y=Symptom))+
  geom_tile(aes(fill=Severity), color="white")+
  geom_tile(aes(fill=Severity), color="white")+
  ylab("Symptom")+
  xlab("Timepoint")+
  scale_fill_manual(values= c("grey", "yellow", "orange", "red"))+
  ggtitle("Volunteer 04")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.title = element_text(size = 20, face= "bold"))+
  theme(axis.text=element_text(size=16), plot.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"))



vac68clin<-plot_grid(c4plot, c8plot, align = "hv", cols=1, labels = "AUTO", label_size = 19)
ggsave("clinvac68.pdf", height=12, width=12)
