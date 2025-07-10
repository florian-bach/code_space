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

long_data <- dplyr::filter(long_data, Individuals %in% c("313", "315", "320", "301", "304", "305", "306", "308", "310"))

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

setwd("/Users/s1249052/PhD/cytof/vac63c/figures")

ggsave("pop_plot.pdf", pop_plot, width=15, height=15)








data <- read.csv("/Users/s1249052/PhD/cytof/vac63c/analysis/pop_plot/poppy_plot.csv")

long_data <- gather(data, Cell_Type, Percentage, colnames(data)[1:5])
long_data <- dplyr::filter(long_data, Individuals %in% c("313", "315", "320", "301", "304", "305", "306", "308", "310"))
long_data$Individuals <- paste("v", long_data$Individuals, sep='')
long_data$Individuals <- factor(long_data$Individuals, levels=c("v313", "v315", "v320", "v301", "v304", "v305", "v306", "v308", "v310"))

long_data$Cell_Type <- gsub(".", "+", long_data$Cell_Type, fixed=T)
long_data$Cell_Type <- factor(long_data$Cell_Type, levels=c("CD4+", "CD8+", "MAIT", "VD2+", "DN"))

long_data$Timepoints <- gsub("C-1", "Baseline", long_data$Timepoints)
long_data$Timepoints <- gsub("T+6", "T6", long_data$Timepoints, fixed = T)
long_data$Timepoints <- gsub("DoD", "Diagnosis", long_data$Timepoints, fixed = T)

long_data <- dplyr::filter(long_data, Timepoints %in% c("Baseline", "DoD", "T6"))

subset_colors <- c("#41B3A3", "#C38D9E", "#E27D60", "#F64C72", "#E8A87C")
names(subset_colors) <- c("CD4+", "MAIT", "CD8+", "VD2+", "DN")

long_data$Percentage <- abs(long_data$Percentage)

ter_poppy_plot <- ggplot(dplyr::filter(long_data, Individuals %in% c("v301", "v304", "v305", "v306", "v308", "v310")), aes(x=factor(Timepoints, levels = c("Baseline", "Diagnosis", "T6")), y=Percentage))+
  geom_bar(stat="identity",position = "stack", aes(fill=Cell_Type))+
  facet_wrap(~Individuals)+
  scale_fill_manual(values=subset_colors, labels=c("CD4+", "CD8+", "MAIT", expression(paste(gamma, delta)), "DN"), guide=guide_legend(label.hjust=0.5))+
  xlab("")+
  ylab("")+
  scale_y_continuous(breaks = seq(0,50, by=10), limits = c(0,55))+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 9),
        legend.position = "none",
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        strip.text = element_text(size=14, color = "white"),
        strip.background = element_rect(fill="#FB027F"))


prim_poppy_plot <- ggplot(dplyr::filter(long_data, Individuals %in% c("v313", "v315", "v320")), aes(x=factor(Timepoints, levels = c("Baseline", "Diagnosis", "T6")), y=Percentage))+
  geom_bar(stat="identity", position = "stack", aes(fill=Cell_Type))+
  facet_wrap(~Individuals)+
  scale_fill_manual(values=subset_colors, labels=c("CD4+", "CD8+", "MAIT", expression(paste(gamma, delta)), "DN"), guide=guide_legend(label.hjust=0.5))+
  xlab("")+
  ylab("")+
  scale_y_continuous(breaks = seq(0,50, by=10), limits = c(0,55))+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 9),
        legend.position = "none",
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        strip.text = element_text(size=14, color = "white"),
        strip.background = element_rect(fill="#6666FF"))


yayuh <- plot_grid(prim_poppy_plot, ter_poppy_plot, ncol = 1, rel_heights = c(1,2))

ggsave("/Users/s1249052/PhD/presentations/NexGenImmunology_2020/pop_boxplot.png", yayuh, width = 7, height=7)  



