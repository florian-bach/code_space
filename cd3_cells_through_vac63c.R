library(tidyr)
library(ggplot2)
library(RColorBrewer)

remove(list = ls())

data <- read.csv("/Users/s1249052/PhD/flow_data/vac63c/cd3_through_time.csv")

colnames(data)[c(2,4)] <- c("C-1", "T+6")

long_data <- gather(data, Timepoint, Frequency, colnames(data)[2:4])

long_data$Frequency <- long_data$Frequency/100

cd3 <- ggplot(long_data, aes(x=factor(Timepoint, levels=c("C-1", "DoD", "T+6")), y=Frequency, group=Volunteer, colour=Volunteer))+
  geom_point(size=3)+
  geom_line(size=2)+
  theme_bw()+
  xlab("Timepoint")+
  ylab("Percentage of White Cells \n")+
  scale_y_continuous(limits=c(0,0.37),labels = scales::percent_format(accuracy=1))+
  scale_colour_brewer(palette="Paired", direction=-1)+
  theme(axis.text = element_text(size=16),
        legend.position = "none",
        axis.title = element_text(size=20))

ggsave("/Users/s1249052/PhD/flow_data/vac63c/cd3_through_time.pdf", cd3)
