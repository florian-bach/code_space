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
data$flo_timepoint <- gsub("am ", "am", data$flo_timepoint)
data$flo_timepoint <- gsub("pm ", "pm", data$flo_timepoint)
data$flo_timepoint <- gsub("DoD ", "Diagnosis", data$flo_timepoint)

long_data <- tidyr::gather(data, Symptom, Severity, colnames(data)[c(12, 14:ncol(data))])

long_data <- mutate(long_data, Volunteer=gsub("69010", "V", long_data$trial_number))

timepoint_levels <- list("Timepoints"=unique(long_data$flo_timepoint[gtools::mixedorder(long_data$flo_timepoint)]))

timepoint_levels$Timepoints[[length(timepoint_levels$Timepoints)+1]] <- timepoint_levels$Timepoints[32]

timepoint_levels <- timepoint_levels$Timepoints[-32]

myLoc <- (which(levels(long_data$flo_timepoint) == "DoD ") +
            which(levels(long_data$flo_timepoint) == "C15 pm")) / 
  2

symptom_heatmap <- ggplot(long_data, aes(x=factor(flo_timepoint, levels=timepoint_levels), y=Symptom))+
  geom_tile(aes(fill=factor(Severity), width=0.93, height=0.93), color=ifelse(grepl("DoD", long_data$flo_timepoint), "black", "lightgrey"))+
  scale_fill_manual(values =  list("lightgrey", "yellow", "orange", "red"))+
  facet_wrap(~Volunteer, scales="free")+
  theme_minimal()+
  guides(fill=guide_legend(title="Severity"))+
  theme(axis.text.x = element_text(hjust=1, angle=60),
        axis.title = element_blank(),
        legend.title = element_blank())

ggsave("/home/flobuntu/PhD/clinical_data/vac69a/figures/symptom_heatmap.png", symptom_heatmap, height=8, width=14)

# make a figure for number of AEs per timepoint

library(dplyr)
library(tidyr)

long_data$flo_timepoint <- factor(long_data$flo_timepoint)


adverse_events <- long_data %>%
  filter(Severity > 0) %>%
  group_by(Volunteer, flo_timepoint, Severity) %>%
  summarise(ae_count = n())


colored_stack <- ggplot(adverse_events,  aes(x=factor(flo_timepoint, levels=timepoint_levels), y=ae_count/6, fill=factor(Severity, levels=paste(rev(1:3)))))+
  geom_bar(stat="identity", position = "stack")+
  scale_fill_manual(values =  list("1"="yellow", "2"="orange", "3"="red"))+
  #facet_wrap(~Volunteer)+
  ylab("Average Number of AEs per Volunteer")+
  xlab("Timepoint")+
  ggtitle("Adverse Events")+
  geom_vline(aes(xintercept = 21.5))+
  guides(fill=guide_legend(title="Severity"))+
  scale_y_continuous(limits = c(0,8), breaks = seq(0, 8, by=2))+
  theme_minimal()+
  theme(axis.text.x = element_text(hjust=1, angle=60, size=8),
        plot.title = element_text(hjust=0.5),
        panel.grid.minor = element_blank())

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/adverse_events_stacked.png", colored_stack, height=4, width=6)



fever <- subset(long_data, long_data$pyrexia_temp>37)

volunteer_colours <- list("V02" = "#FB9A99",
                          "V03" = "#E31A1C",
                          "V05" = "#A6CEE3",
                          "V06" = "#1F78B4",
                          "V07" = "#B2DF8A",
                          "V09" = "#33A02C")


volunteer_palette <- unlist(unname(volunteer_colours))
names(volunteer_palette) <- names(volunteer_colours)

fever_curves <- ggplot(fever, aes(x=factor(flo_timepoint, levels=timepoint_levels), y=pyrexia_temp, color=Volunteer, group=Volunteer))+
  scale_fill_manual(values=volunteer_palette)+
  scale_color_manual(values=volunteer_palette)+
  geom_line(aes(color=Volunteer), size=1.1)+
  geom_point(fill="white", stroke=1, shape=21)+
  ggtitle("Fever")+
  xlab("Timepoint")+
  ylab(expression(paste("Temperature ",degree,"C",sep="")))+
  scale_y_continuous(limits=c(37.5, 40), breaks = seq(37.5, 40, by=0.5) )+
  theme_minimal()+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=60, size=8))

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/fever_curves.png", fever_curves)

clinical_graphs <- cowplot::plot_grid(colored_stack, fever_curves, ncol = 2, rel_widths = c(1.6,1))

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/aes_and_fever_curves.png", clinical_graphs, height=4, width=8*4/3)

