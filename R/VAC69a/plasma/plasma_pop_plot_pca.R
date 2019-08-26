library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbiplot)


data <- read.csv("/Users/s1249052/PhD/plasma/vac69a/vac69a_Vivax_plasma_analytes_no_inequalities.csv")
data<- na.omit(data)
data$timepoint <- gsub("D", "DoD", data$timepoint, fixed=T)
data$timepoint <- gsub("DoDoDoD+6", "T+6", data$timepoint, fixed=T)
data$timepoint <- gsub("DoD+6", "T+6", data$timepoint, fixed=T)

long_data <- gather(data, Analyte, Concentration, colnames(data)[3:37])

long_data$Analyte <- substr(long_data$Analyte, 1, nchar(long_data$Analyte)-7)

long_data$Analyte <- gsub(".", "", long_data$Analyte, fixed = T)
long_data$Concentration <- as.numeric(long_data$Concentration)



ggplot(long_data, aes(x=factor(long_data$timepoint, levels=c("C-1", "DoD", "T+6", "C+45")), y=Concentration, color=Volunteer), group=Volunteer)+
  geom_point()+
  geom_line(aes(group=long_data$Volunteer))+
  facet_wrap(~ Analyte, scales = "free")


####     pca plot

data <- spread(long_data, Analyte, Concentration)


plasma_pca <- prcomp(data[,3:37], center = T, scale. = T)

summary(plasma_pca)  

pca <- ggbiplot(plasma_pca, choices=2:3, obs.scale = 1, var.scale = 1,var.axes=FALSE, groups = data$Volunteer)

plotly::ggplotly(pca+theme(legend.title = element_blank()),tiptool="timepoint")

timepoint

