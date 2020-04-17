library(ggplot2)
library(tidyr)
library(dplyr)
library(colorspace)

vac69b_palette <- c("#FB9A99","#A6CEE3", "#B2DF8A", "#D33F6A", "#E99A2C")
names(vac69b_palette) <- c("v2","v5","v7","v11","v21")

parasitaemia_theme <-  theme(legend.title = element_blank(),
                             legend.text = element_text(size = 20) ,
                             axis.line = element_line(colour = "black"),
                             axis.text.x = element_text(angle = 60, hjust = 1, size=18, color="black"),
                             axis.title.x = element_text(size=24, color="black"),
                             axis.title.y = element_text(size=24, color="black"),
                             axis.text.y = element_text(size=20, color="black"),
                             legend.position = "top")



data <- read.csv("~/PhD/clinical_data/vac69b/parasitaemia/vac69b_parasitaemia.csv", header=T, stringsAsFactors = F)
ctrls <- subset(data, data$Vaccination=="Control")

data <- gather(data, Timepoint, Parasitaemia, colnames(data)[5:ncol(data)])
ctrls <- gather(ctrls, Timepoint, Parasitaemia, colnames(ctrls)[5:ncol(ctrls)])


ctrls$Volunteer <- factor(ctrls$Volunteer, levels=c("v2","v5","v7","v11","v21"))

ctrl_plot <- ggplot(ctrls[!is.na(ctrls$Parasitaemia),], aes(x=factor(Timepoint, levels = unique(mixedsort(Timepoint))), y=Parasitaemia, group=Volunteer))+
  geom_point(aes(color=Volunteer), size=2.5)+
  geom_line(aes(color=Volunteer), size=2)+
  scale_color_manual(values=vac69b_palette)+
  theme_bw()+
  xlab("Day of Infection")+
  ylab("Genome Copies / mL")+
  scale_y_continuous(trans="log10", limits=c(1, 27000), breaks=c(10, 100, 1000, 10000, 25000))+
  scale_x_discrete(limits=c("Baseline", paste("D", seq(6, 23), sep='')))+
  parasitaemia_theme

ggsave("~/PhD/clinical_data/vac69b/figures/vac69b_parasitaemia.png", ctrl_plot, height = 8, width=10)


vacc_plot <- ggplot(data[!is.na(data$Parasitaemia),], aes(x=factor(Timepoint, levels = unique(mixedsort(Timepoint))), y=Parasitaemia, group=Volunteer))+
  geom_point(aes(color=factor(Vaccination)), size=2.5)+
  geom_line(aes(color=factor(Vaccination)), size=2)+
  #scale_color_manual(values=vac69b_palette)+
  theme_bw()+
  xlab("Day of Infection")+
  ylab("Genome Copies / mL")+
  scale_y_continuous(trans="log10", limits=c(1, 27000), breaks=c(10, 100, 1000, 10000, 25000))+
  scale_x_discrete(limits=c("Baseline", paste("D", seq(6, 23), sep='')))+
  parasitaemia_theme

ggsave("~/PhD/clinical_data/vac69b/figures/vac69b_parasitaemia_vaccine.png", vacc_plot, height = 8, width=10)




vac69a <- read.csv("/home/flobuntu/PhD/clinical_data/vac69a/parasitaemia/vac69a_parasitaemia.csv", stringsAsFactors = F)
vac69a$Volunteer <- paste("v", substr(vac69a$Volunteer, nchar(vac69a$Volunteer), nchar(vac69a$Volunteer)), sep='')

vac69a <- subset(vac69a, Volunteer %in% c("v2", "v5", "v7"))
vac69a<- data.frame(Volunteer=vac69a[,1], N_Infection="Primary", vac69a[,2:ncol(vac69a)])                          
colnames(vac69a)[3] <- "Baseline"

vac69a <- gather(vac69a, Timepoint, Parasitaemia, colnames(vac69a)[3:ncol(vac69a)])


combined <- rbind(vac69a, select(ctrls, colnames(vac69a)))
combined$Sample_ID <- paste(combined$Volunteer, combined$N_Infection, sep='_')
combined$Sample_ID <- factor(combined$Sample_ID, levels = unique(mixedsort(combined$Sample_ID)))

my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", vac69b_palette[c("v11","v21")])
names(my_paired_palette) <- c(levels(combined$Sample_ID))

n_infection_plot_groups <- ggplot(combined[!is.na(combined$Parasitaemia),], aes(x=factor(Timepoint, levels = unique(mixedsort(Timepoint))), y=Parasitaemia, group=Sample_ID))+
  geom_point(aes(color=factor(N_Infection)), size=2.5)+
  geom_line(aes(color=factor(N_Infection)), size=2)+
  #scale_color_manual(values=vac69b_palette)+
  theme_bw()+
  xlab("Day of Infection")+
  ylab("Genome Copies / mL")+
  scale_y_continuous(trans="log10", limits=c(1, 27000), breaks=c(10, 100, 1000, 10000, 25000))+
  scale_x_discrete(limits=c("Baseline", paste("D", seq(6, 23), sep='')))+
  parasitaemia_theme

ggsave("~/PhD/clinical_data/vac69b/figures/vac69AB_parasitaemia_n_infection_groups.png", n_infection_plot_groups, height = 8, width=10)


n_infection_plot_indie <- ggplot(combined[!is.na(combined$Parasitaemia),], aes(x=factor(Timepoint, levels = unique(mixedsort(Timepoint))), y=Parasitaemia, group=Sample_ID))+
  geom_point(aes(color=factor(Sample_ID)), size=2.5)+
  geom_line(aes(color=factor(Sample_ID)), size=2)+
  scale_color_manual(values=my_paired_palette)+
  theme_bw()+
  xlab("Day of Infection")+
  ylab("Genome Copies / mL")+
  scale_y_continuous(trans="log10", limits=c(1, 27000), breaks=c(10, 100, 1000, 10000, 25000))+
  scale_x_discrete(limits=c("Baseline", paste("D", seq(6, 23), sep='')))+
  parasitaemia_theme

ggsave("~/PhD/clinical_data/vac69b/figures/vac69AB_parasitaemia_n_infection_indie.png", n_infection_plot_indie, height = 8, width=10)

