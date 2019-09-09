library(tidyr)
library(dplyr)
library(ggplot2)
library(colorspace)

# iMac
data <- read.csv("/Users/s1249052/PhD/cytof/vac69a/ics/ics_flow_experiments_060919.csv")
data[,colnames(data)[4:length(colnames(data))]] <- data[,colnames(data)[4:length(colnames(data))]]/100

data$Individuals <- paste0("0", data$Individuals, sep='')
data$Individuals <- as.character(data$Individuals)

colnames(data) <- gsub("..", " ", colnames(data), fixed=T)
colnames(data) <- gsub(".", " ", colnames(data), fixed=T)


long_data <- gather(data, Cytokine, Percentage, colnames(data)[4:length(colnames(data))])
long_data_no_t_cells <- filter(long_data, !grepl("T cells", long_data$Cytokine))





c6_data <- filter(long_data, long_data$Timepoint=="C+6")
c6_t_cells <- filter(c6_data, grepl("T cells",c6_data$Cytokine))
c6_t_cells$Cytokine <- substr(c6_t_cells$Cytokine, 1, nchar(c6_t_cells$Cytokine)-10)

dod_data <- filter(long_data, long_data$Timepoint=="DoD")
dod_t_cells <- filter(dod_data, grepl("T cells", dod_data$Cytokine))
dod_t_cells$Cytokine <- substr(dod_t_cells$Cytokine, 1, nchar(dod_t_cells$Cytokine)-10)


my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")
names(my_paired_palette) <- c("02", "03", "05", "06", "07", "09")


(c6_plot <- ggplot(c6_t_cells, aes(x=Timepoint, y=Percentage, fill=factor(Individuals)))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_grid(Conditions~Cytokine, scales = "free")+
  scale_fill_manual(values=my_paired_palette, name="Volunteer")+
  scale_y_continuous(labels = scales::percent_format(accuracy=0.1))+
  ggtitle("C+6, all conditions")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))
  )

ggsave("/Users/s1249052/PhD/cytof/vac69a/ics/figures/all_conditions_c6.png", c6_plot)

(dod_plot <- ggplot(dod_t_cells, aes(x=Timepoint, y=Percentage, fill=factor(Individuals)))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_grid(Conditions~Cytokine, scales = "free")+
  scale_fill_manual(values=my_paired_palette, name="Volunteer")+
  scale_y_continuous(labels = scales::percent_format(accuracy=0.1))+
  ggtitle("DoD, 1/5 only")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))
)

ggsave("/Users/s1249052/PhD/cytof/vac69a/ics/figures/dil_only_dod.png", dod_plot)

time_t_cells <- filter(long_data, grepl("T cells", long_data$Cytokine))
time_t_cells <- filter(time_t_cells, grepl("dil", time_t_cells$Conditions))

time_plot <- ggplot(time_t_cells, aes(x=Timepoint, y=Percentage, group=factor(Individuals), color=factor(Individuals)))+
    geom_point()+
    geom_line()+
    facet_grid(Conditions~Cytokine, scales = "free")+
    scale_colour_manual(values=my_paired_palette, name="Volunteer")+
    scale_y_continuous(labels = scales::percent_format(accuracy=0.1))+
    ggtitle("C+6 vs DoD, 1/5 only")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5))
ggsave("/Users/s1249052/PhD/cytof/vac69a/ics/figures/dil_only_c6_and_dod.png", time_plot)




#############   split by cd4 and cd8 cells    #############


# 
# 
# long_data_no_t_cells$Cell <- ifelse(grepl("CD4 neg", long_data_no_t_cells$Cytokine)==T, "CD4-", "CD4+")
# 
# 
# 
# c6_data <- filter(long_data_no_t_cells, long_data_no_t_cells$Timepoint=="C+6")
# c6_data$Cytokine <- substr(c6_data$Cytokine, 1, 4)
# c6_data$Cytokine <- gsub(" ", "", c6_data$Cytokine)
# 
# dod_data <- filter(long_data_no_t_cells, long_data_no_t_cells$Timepoint=="DoD")
# dod_data$Cytokine <- substr(dod_data$Cytokine, 1, 4)
# dod_data$Cytokine <- gsub(" ", "", dod_data$Cytokine)
# 
# my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")
# names(my_paired_palette) <- c("02", "03", "05", "06", "07", "09")
# 
# 
# 
# (c6_plot <- ggplot(c6_data, aes(x=Timepoint, y=Percentage, fill=factor(Individuals)))+
#     geom_bar(stat = "identity", position = "dodge", colour="black")+
#     facet_grid(Conditions~Cytokine+Cell, scales = "free")+
#     scale_fill_manual(values=my_paired_palette, name="Volunteer")+
#     scale_y_continuous(labels = scales::percent_format(accuracy=0.1))+
#     ggtitle("C+6, all conditions")+
#     theme_minimal()+
#     theme(plot.title = element_text(hjust = 0.5))
# )
# 
# 
# (dod_plot <- ggplot(dod_data, aes(x=Timepoint, y=Percentage, fill=Individuals))+
#     geom_bar(stat = "identity", position = "dodge")+
#     facet_grid(Conditions~Cytokine+Cell, scales = "free")+
#     scale_fill_manual(values=my_paired_palette, name="Volunteer")+
#     scale_y_continuous(labels = scales::percent_format(accuracy=0.01))+
#     ggtitle("DoD, 1/5 only")+
#     theme_minimal()+
#     theme(plot.title = element_text(hjust = 0.5))
# )
# 
# 
# time_t_cells <- filter(long_data, grepl("T cells", long_data$Cytokine))
# time_t_cells <- filter(time_t_cells, grepl("dil", time_t_cells$Conditions))
# 
# (time_plot <- ggplot(time_t_cells, aes(x=Timepoint, y=Percentage, group=factor(Individuals), color=factor(Individuals)))+
#     geom_point()+
#     geom_line()+
#     facet_grid(Conditions~Cytokine, scales = "free")+
#     scale_fill_manual(values=my_paired_palette, name="Volunteer")+
#     scale_y_continuous(labels = scales::percent_format(accuracy=0.1))+
#     ggtitle("C+6 vs DoD, 1/5 only")+
#     theme_minimal()+
#     theme(plot.title = element_text(hjust = 0.5))
# )
# 
# 
# 





¢¢¢¢¢¢¢¢