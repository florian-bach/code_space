library(ggplot2)
library(tidyr)
library(dplyr)

data <- read.csv("/Users/s1249052/PhD/cytof/vac63c/analysis/activation_gating_stats.csv")

colnames(data)

data <- data[-grep("control", data[,2]),]

data$Individuals <- as.character(data$Individuals)
data$Timepoints <- as.character(data$Timepoints)

colnames(data)[3:8] <- c("T cells", "CD4+", "gd", "CD8+", "MAIT", "DN")

long_data <- gather(data, Cells, Percentage, colnames(data)[3:8])



long_data$order_pop <- factor(long_data$Cells, levels = c("T cells", "CD4+", "CD8+", "gd", "MAIT", "DN"))
long_data$order_vol <- factor(long_data$Individuals, levels = c("313", "315", "320", "302", "307", "301", "304", "305", "306", "308", "310"))

(pop_plot <- ggplot(long_data, aes(x=factor(Timepoints, levels=c("C-1", "DoD", "T+6", "C+45")), y=Percentage, group=order_vol, fill=order_vol))+
  geom_bar(stat="identity", position="dodge")+
  facet_grid(order_pop~order_vol, scales="free")+
  theme_bw()+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  xlab("Timepoint")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, size=12),
        strip.text.x = element_text(size = 10),
        axis.text.y = element_text(size=14))
)


ggsave("pop_plot.pdf", pop_plot, width=15, height=15)
