library(ggplot2)
library(tidyr)
library(dplyr)
library(gtools)
library(cowplot)


setwd("C:/Users/Florian/PhD/presentations/icpvr2019/")


data <- read.csv("boid_fig_cont.csv")
colnames(data)[1] <- "Patient"

#levels=c("201-500", "501-1000", "1000-2000", "2000-4000", "4000-6000" , "6000-8000", "8000-10000", ">10000") 


# set up boundaries for intervals/bins
breaks <- c(0,10,50,100,200, 500,1000,2000,4000,6000, 8000, 10000, 13000)
# specify interval/bin labels
labels <- c("<10",  "11-50", "51-100", "101-200", "201-500", "501-1000", "1001-2000", "2001-4000", "4001-6000" , "6001-8000", "8001-10000", ">10000")
# bucketing data points into bins
bins_first <- cut(data$First.Fever, breaks, include.lowest = T, right=FALSE, labels=labels)
bins_last <- cut(data$Last.Fever, breaks, include.lowest = T, right=FALSE, labels=labels)
# inspect bins
summary(bins_first)

y <- cbind(data, bins_first, bins_last)

boyd_plot <- ggplot()+
geom_point(data=y, aes_(x=y$bins_first, y=y$bins_last), size=4, color="red")+
  xlab("Parasitaemia at First Fever")+
  ylab("Parasitaemia at Last Fever")+
  scale_x_discrete(drop=FALSE)+
  scale_y_discrete(drop=FALSE)+
  geom_abline(slope=1, intercept=0, size=1.2, linetype = "dashed")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=23),
        axis.title.x = element_text(size=23),
        legend.position = "bottom",
        legend.title = element_blank())
  

ggsave("boyd_plot.png", height=8, width=8)


data <- read.csv("ciuca.csv")
colnames(data)[1] <- "Order of Infection"

ciuca_plot <- ggplot(data, aes(x=`Order of Infection`, y=Tolerant))+
  geom_bar(stat="identity", fill="#ff1493")+
  scale_y_continuous(labels = scales::percent)+
  xlab("Order of Infection")+
  ylab("No Fever Despite Patent Parasitaemia")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        legend.title = element_blank())

ggsave("ciuca_plot.png", ciuca_plot, height = 8, width=8)




data <- read.csv("collins.csv")


head(data)

data$Max_count_ul <- as.character(data$Max_count_ul)
data$Max_count_ul <- gsub(",", "", data$Max_count_ul, fixed=T)

data$Max_temp <- as.character(data$Max_temp)
data$days_above_38.3 <- as.character(data$days_above_38.3)

parasite_plot <- ggplot(data, aes(x=Infection, y=as.numeric(Max_count_ul), fill=Infection))+
  geom_boxplot()+
  ylab(expression(paste("Maximum Parasite Density /  " , mu*"L", sep = '')))+
  scale_fill_manual(values = c("#f3d250", "#ff1493"))+
  theme_bw()+
  theme(legend.title = element_blank(), 
        legend.position = "none",
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=23),
        axis.title.x = element_text(size=23))


ggsave("collins_max_parasites.png", parasite_plot, width = 5, height = 6)

max_feverplot <- ggplot(data, aes(x=Infection, y=as.numeric(Max_temp), fill=Infection))+
  geom_boxplot()+
  ylab(expression(paste("Maximum Body Temperature")))+
  scale_fill_manual(values = c("#f3d250", "#ff1493"))+
  theme_bw()+
  theme(legend.title = element_blank())

duration_feverplot <- ggplot(data, aes(x=Infection, y=as.numeric(days_above_38.3), fill=Infection))+
  geom_boxplot()+
  ylab(expression(paste("Days With Fever Over 38.3°C")))+
  scale_fill_manual(values = c("#f3d250", "#ff1493"))+
  theme_bw()+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=23),
        axis.title.x = element_text(size=23))


ggsave("collins_duration_feverplot.png", duration_feverplot, width = 5, height = 6)





