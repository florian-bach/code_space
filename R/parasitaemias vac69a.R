library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(tidyr)
library(gtools)


setwd("/Users/s1249052/PhD/oxford/vac69")

data <- read.csv("vac69a_parasitaemia.csv", header=T)

# get rid of garbage timepoints that mess up graph

parasitaemias <- gather(data, Timepoint, Genomes, c("C.1", colnames(data)[12:35]) )
parasitaemias[,2:13] <- NULL

parasitaemias$Timepoint <- gsub("C.1", "C-1", parasitaemias$Timepoint) 


my_palette <- c("#588C73", "#D53E4F", "#F2E394",   "#F2AE72",   "#D96459",   "#1A9CC7")
my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")

ggplot(data=parasitaemias[!is.na(parasitaemias$Genomes),], aes(x=factor(Timepoint, levels = unique(mixedsort(parasitaemias$Timepoint))), y=Genomes, group=factor(Volunteer)))+
  geom_point(aes(color=factor(Volunteer)), size=2.5)+
  geom_line(aes(color=factor(Volunteer)), size=2)+
  scale_color_manual(values=my_paired_palette)+
  theme_bw()+
  xlab("Day of Infection")+
  ylab("Genome Copies / mL")+
  scale_y_continuous(trans="log10", limits=c(1, 27000), breaks=c(10, 100, 1000, 10000, 25000))+
  theme(legend.title = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1, size=12, color="black"),
        axis.title.x = element_text(size=16, color="black"),
        axis.title.y = element_text(size=16, color="black"),
        axis.text.y = element_text(size=14, color="black"))


c("#1A87C7", "#FFFFFF", "#FFFFFF")warnings()
geom_line(data=dfr[!is.na(dfr$y),])scale_color_brewer(palette="Set1")+
  scale_y_continuous(trans="log10", limits=c(10, 40000), breaks=c(10, 100, 1000, 10000, 30000))+
  scale_x_continuous(breaks=seq(6, 15, by=1))+
  theme(text = element_text(size=17, face="bold", color="black"),
        axis.text.x = element_text(size=14, color="black"),
        axis.text.y = element_text(size=14, color="black"))+
  coord_fixed(ratio=2)+
  geom_smooth(method = "loess", se=TRUE, aes=(color=factor(Volunteer)), formula = y ~ x)
 
ggsave ("C:/users/Florian/PhD/reports/figures/parasitaemias vac68.pdf", width=8, height=6)

attach(parasitaemias)

colnames(parasitaemias)<-c('Volunteer', 'Timepoint', 'Genome_Copies')

para04<-subset(parasitaemias, parasitaemias$Volunteer=='01-004')
para08<-subset(parasitaemias, parasitaemias$Volunteer=='01-008')

model1<-glm(Genome_Copies~Timepoint*Volunteer, data=parasitaemias)
AIC(model1)
model1<-glm(Genome_Copies~Timepoint+Volunteer, data=parasitaemias)
AIC(model1)#661.929

######## individual volunteers not significantly different
model1<-glm(Genome_Copies~Timepoint, data=parasitaemias)
AIC(model1)#660.1242

AImodel2<-glm(Genome_Copies~Timepoint+I(Timepoint)^2, data=parasitaemias)
AIC(model2)#661.929
#make data frames for each person for linear regression, starting at d8.5

frame04 = parasitaemias[seq(9, nrow(parasitaemias), 2), ]
frame08 = parasitaemias[seq(10, nrow(parasitaemias), 2), ]

# drop letter D from time point label

frame04[] <- lapply(frame04, gsub, pattern='D', replacement='')
frame08[] <- lapply(frame08, gsub, pattern='D', replacement='')

#make timepoints float

frame04$variable<-as.numeric(frame04$variable)
frame08$variable<-as.numeric(frame08$variable)

#make values float

frame04$value<-as.numeric(frame04$value)
frame08$value<-as.numeric(frame08$value)


#linear regression

model04<-lm(value~variable, data=frame04)
print(model04)

model04a<-glm(value~variable+I(variable^2), data=frame04)
print(model04a)
3861*x -38610

model08<-glm(value~variable, data=frame08)
print(model08)

fun04<-function(x) 3861*x -38610
fun08<-function(x) -4487 + 897.8*x
stat_function(fun = fun.1) + xlim(-5,5)

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p+
  stat_function(fun = fun04)+
  xlim(8.5, 14.5)+
  scale_y_continuous(trans="log10", limits=c(10, 40000), breaks=c(10, 100, 1000, 10000, 40000))




para04<-subset(parasitaemias, parasitaemias$Volunteer=='01-004')
time2<-para04$Timepoint^2
linear<-lm(Genome_Copies~Timepoint, data=para04)
quadratic<-lm(Genome_Copies~Timepoint+time2, data=para04)
forcelin<-lm(Genome_Copies~0+Timepoint, data=para04)
forcequad<-lm(Genome_Copies~0+Timepoint+I(Timepoint^2), data=para04)

> AIC(linear)
# [1] 348.1928
> AIC(quadratic)
# [1] 328.8093
coef(quadratic)

AIC(forcelin)
AIC(forcequad)
coef(forcequad)

# Timepoint I(Timepoint^2) 
# -2240.3213       238.9384 

fun04<-


para08<-subset(parasitaemias, parasitaemias$Volunteer=='01-008')
time2<-para08$Timepoint^2
linear<-lm(Genome_Copies~Timepoint, data=para08)
quadratic<-lm(Genome_Copies~Timepoint+time2, data=para08)
> AIC(linear)
# [1] 313.6454
> AIC(quadratic)
# [1] 295.4272
coef(quadratic)

forcelin<-lm(Genome_Copies~0+Timepoint, data=para08)
forcequad<-lm(Genome_Copies~0+Timepoint+I(Timepoint^2), data=para08)
AIC(forcelin)
AIC(forcequad)

coef(forcequad)




##########################        blablablablabla




# volunteer 04

25*x^3.25=31010
x^3.25=1240.4
1240.4^(1/3.25)
8.95 = PMR




#volunteer 08

79*x^2.75=16717
16388/79

x^2.75= 6078.909
x^2.5=
209^(1/2.75)
6.977335
207.443^(1/2.5)
8.473416



