library(ggplot2)
library(tidyr)
library(dplyr)


setwd("C:/Users/Florian/PhD/cytof/vac69a/clinical_data")

data <- read.csv("symptoms.csv", header=T)


colnames(data)[1] <- "Volunteer"

data$Volunteer <- paste("Volunteer", substr(data$Volunteer, 6, 7), sep='')

data$timepoint <- gsub("_am", "", data$timepoint, fixed=T)
data$timepoint <- gsub("_pm", ".5", data$timepoint, fixed=T)
data$timepoint <- substr(data$timepoint, nchar(as.character(data$timepoint))-3, nchar(as.character(data$timepoint)))
data$timepoint <- gsub("_", "", data$timepoint, fixed=T)
data$timepoint <- paste("C+", data$timepoint, sep='')

data$timepoint <- gsub("C+ep", "EP+", data$timepoint, fixed=T)
data$timepoint <- gsub("C+c", "C+", data$timepoint, fixed=T)
data$timepoint <- gsub("C+6", "T+6", data$timepoint, fixed=T)


data$pyrexia_temp <- NULL
data$X <- NULL


long_data <- gather(data, Symptom, Severity, colnames(data)[3:16]) 



for(i in unique(long_data$Volunteer)){
  ifelse(i %in% c("Volunteer  02","Volunteer  06"), assign("result", element_text(size=35)), assign("result", element_blank()))
  ifelse(i %in% c("Volunteer  05","Volunteer  09"), assign("result1", "right"), assign("result1", "left"))
  
  sub_set <- filter(long_data, Volunteer == i)

  assign(i,
         

  ggplot(sub_set, aes(x=timepoint, y=Symptom))+
  geom_tile(aes(fill=as.character(Severity)), color="white")+
  scale_fill_manual(values= c("grey", "yellow", "orange", "red"))+
  theme(panel.border = element_blank(),
          axis.text.y.left = result,
          axis.line.y.left = element_blank(),
          axis.line.y.right = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 33),
          axis.text.y.right = element_text(size = 35),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.title = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 45, hjust = 0.5),
          plot.margin = unit(c(1,0,1,0), "cm"))
  )
}


ggsave("symptom_plot.png", symptom_plot, width=15, height=10)


plot_grid(`Volunteer  02`, `Volunteer  03`, `Volunteer  05`, `Volunteer  06`, `Volunteer  07`, `Volunteer  09`)
