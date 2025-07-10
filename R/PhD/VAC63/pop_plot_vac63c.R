library(tidyr)
library(dplyr)
library(ggplot2)
library(colorspace)
library(cowplot)



slimmer_data <- read.csv("/home/florian/PhD/cytof/vac63c/activation_gating_stats.csv", header=T)

colnames(slimmer_data)[3:8] <- c("T Cells", "CD4+", "gd", "CD8+", "MAIT", "DN")

slimmer_data <- slimmer_data[-grep("control", slimmer_data[,2]),]            

slimmer_data[,1] <- as.character(slimmer_data[,1])
slimmer_data[,2] <- as.character(slimmer_data[,2])

long_data <- gather(slimmer_data, Cells, Percentage, colnames(slimmer_data)[3:8])



long_data$order_vol <- factor(long_data$Individuals, levels=c("313", "315", "320", "302", "307", "301", "304", "305", "306", "308", "310"))
long_data$order_cell <- factor(long_data$Cells, levels=c("T Cells", "CD4+", "CD8+", "gd", "MAIT", "DN"))

(pop_plot <- ggplot(long_data, aes(x=factor(Timepoints, levels=c("C-1", "DoD", "T+6", "C+45")), y=Percentage, group=order_vol, fill=order_vol))+
  geom_bar(stat="identity", position="dodge")+
  facet_grid(order_cell~order_vol, scales="free")+
  scale_y_continuous(labels = scales::percent_format())+
  theme_bw()+
  xlab("Timepoint")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, size=12),
        strip.text.x = element_text(size = 10)))
