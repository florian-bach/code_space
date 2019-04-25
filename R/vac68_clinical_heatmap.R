library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggthemes)
library(cowplot)
library(plyr)
library(gtools)
library(dplyr)
library(tidyr)

setwd("/Users/s1249052/PhD/clinical data/vac68")
clinic<-read.csv("symptoms_vac68.csv", header=T)

#clinic2<-mutate(clinic, volunteer=substring(clinic$trial_number, 6, 8))

clinic3 <- gather(clinic, symptom, severity, 10:26)

clinic4 <- clinic3 %>%
            select(volunteer=trial_number, symptom=symptom, severity=severity, timepoint=timepoint)

clinic4<-na.omit(clinic4)

#get rid of temperature & unsolicited aes as symptom
clinic4<-clinic4[!grepl("pyrexia_temp", clinic4$symptom),]
clinic4<-clinic4[!grepl("review_ae_unsol", clinic4$symptom),]

clinic4$timepoint <- gsub("am", "", clinic4$timepoint)
clinic4$timepoint <- gsub("pm", ".5", clinic4$timepoint)

### make function that capitalises things

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

###

clinic4$symptom <- sapply(clinic4$symptom, simpleCap)
clinic4$symptom <- gsub("_", " ", clinic4$symptom)
clinic4$volunteer <- gsub(68010, "V", clinic4$volunteer)
###

clinic5 <- clinic4[clinic4$timepoint %in% unique(clinic4$timepoint)[1:22],]

ggplot(clinic5, aes(x=factor(timepoint, levels=mixedsort(unique(clinic5$timepoint))), y=symptom))+
     geom_tile(aes(fill=as.factor(severity)), color="white")+
     geom_vline(xintercept=match("C+14", mixedsort(unique(clinic5$timepoint)))+0.5)+
     facet_wrap( ~ volunteer, nrow=1)+
     scale_fill_manual(values= c("grey", "yellow", "orange", "red"))+
     ggtitle("Symtpoms in VAC68")+
     labs(fill="Severity")+
     theme(axis.text.x = element_text(angle = 60, hjust = 1),
           axis.title = element_blank(),
           axis.text=element_text(size=16),
           plot.title = element_text(size = 20, face = "bold"),
           legend.title.align=0.5)
  
ggsave("clinvac68.pdf")
           