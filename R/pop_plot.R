library(tidyr)
library(dplyr)
library(ggplot2)
library(colorspace)

# laptop
data <- read.csv("C:/Users/Florian/PhD/cytof/vac69a/Vac69a_michalina_compensated_Exported_Stats\ 4.csv")
# iMac
# data <- read.csv("/Users/s1249052/PhD/cytof/vac69a/Vac69a_michalina_compensated_Exported_Stats\ 4.csv")

str(data)
colnames(data) <-c("CD4+", "Vd2+", "CD8+", "MAIT", "Tregs", "DN", "Activated", "Gate", "Timepoint", "Volunteer") 
head(data)


long_data <- gather(data, Population, Percentage, colnames(data)[1:7])
long_data$Gatef <- factor(long_data$Gate, levels=c("All T Cells", "Activated T cells"))
my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")


pop_plot <- ggplot(long_data, aes(x=factor(Timepoint, levels=c("C-1", "C+8", "C+10", "C+12", "DoD", "T+6")), y=Percentage, group=Volunteer, fill=Volunteer))+
         geom_bar(stat="identity", position="dodge")+
         facet_grid(Population~Volunteer+Gatef, scales="free")+
         scale_fill_manual(values=my_paired_palette)+
         theme_bw()+
         xlab("Timepoint")+
         theme(legend.position = "none",
               axis.text.x = element_text(angle = 60, hjust = 1, size=12),
               strip.text.x = element_text(size = 10))


setwd("/Users/s1249052/PhD/cytof/better_gating/double_flowsoms/figures")

ggsave("pop_plot.png", pop_plot, width = 14, height=8)



# convert number to percentages, get rid of superfluous columns & rows

dod6_data <- subset(data, Timepoint %in% c("C-1","DoD", "T+6"))
dod6_data[dod6_data$Gate=="Activated T cells",] <- NA
dod6_data <- na.omit(dod6_data)
dod6_data[,1:6] <- dod6_data[,1:6]/100
dod6_data[,1:6] <- dod6_data[,1:6] * dod6_data$Activated
dod6_data$Activated <- NULL
dod6_data$Tregs <- NULL

long_dod6_data <- gather(dod6_data, Population, Percentage, colnames(dod6_data)[1:5])
# 
# # create ymin and ymax columns for geom_rect() by subsetting dataframe and doing cumsum
# 
# pie_dod6_data <- long_dod6_data%>%
#   dplyr::arrange(Volunteer, Timepoint, Population) %>%
#   group_by(Timepoint, Volunteer) %>%
#   mutate(ymax=cumsum(Percentage)) %>%
#   ungroup()
# 
# pie_dod6_data$ymin <- c(0, pie_dod6_data$ymax[1:nrow(pie_dod6_data)-1])
# pie_dod6_data$Population <- gsub("+", "", pie_dod6_data$Population, fixed = T)
# 

# for (i in unique(pie_dod6_data$Volunteer)){
#   
#   tmp_dat <- NULL
#   tmp_dat <- dplyr::filter(pie_dod6_data, pie_dod6_data$Volunteer==i)
#   
#   tmp_dat$ymin <- tmp_dat$ymin/max(tmp_dat$ymax)+4
#   tmp_dat$ymax <- tmp_dat$ymax/max(tmp_dat$ymax)+4
#   
#   assign(paste(unique(tmp_dat$Volunteer), "_ratio", sep=''), max(tmp_dat$ymax)/8.4675213)
#   
#   tmp_dat$ymin <- tmp_dat$ymin/max(tmp_dat$ymax)+4
#   tmp_dat$ymax <- tmp_dat$ymax/max(tmp_dat$ymax)+4
#   
#   
#   assign(
#     paste("V", substr(unique(tmp_dat$Volunteer),11,12), "_pop_pie_plot", sep=''), 
#     
#     ggplot()+
#       # clusters and subset pies
#       geom_rect(data=tmp_dat, aes_(fill=factor(tmp_dat$Population, levels=c("CD4", "CD8", "MAIT", "Vd2", "DN")), ymin=tmp_dat$ymin, ymax=tmp_dat$ymax, xmax=6, xmin=4))+
#       
#       scale_fill_manual(values=qualitative_hcl(6, "Dark3"))+
#       facet_wrap(~Timepoint)+
#       theme(aspect.ratio=1,
#             axis.text = element_blank(),
#             axis.title = element_blank(),
#             axis.line = element_blank(),
#             axis.ticks = element_blank())+
#       coord_polar(theta = "y")+
#       theme_void()
#     )
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 

long_dod6_data <- gather(dod6_data, Population, Percentage, colnames(dod6_data)[1:5])

long_dod6_data$Population <- gsub("+", "", long_dod6_data$Population, fixed = T)
long_dod6_data$Popf <- factor(long_dod6_data$Population, levels=c("CD4", "CD8", "MAIT", "Vd2", "DN"))

test <- long_dod6_data[order(long_dod6_data$Popf),]

colour_fixing <- c("#ff1493", "#f3d250", "#90ccf4", "#D96459", "#0bb38f")
names(colour_fixing) <- c("CD4", "MAIT", "CD8", "Vd2", "DN")


(dod6_pop_plot <- ggplot(test, aes(x=factor(Timepoint, levels=c("C-1", "DoD", "T+6")), group=Volunteer, fill=Population, y=Percentage))+
    geom_bar(stat="identity", position="stack")+
    facet_wrap(~Volunteer)+
    scale_fill_manual(values=colour_fixing)+
    theme_bw()+
    xlab("Timepoint")+
    scale_y_continuous(labels=function(x) paste0(x,"%"))+
    ylab("Percenentage Activated")+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=12),
          legend.title = element_blank(),
          axis.text.y = element_text(size=12),
          strip.text = element_text(size=16, face = "bold"),
          legend.text = element_text(size=16),
          strip.background = element_blank(),
          axis.title.y = element_text(size=16, face = "bold"),
          axis.title.x = element_text(size=16, face = "bold")))



setwd("C:/Users/Florian/PhD/cytof/vac69a/double_flowsoms/figures/")

ggsave("dod6_pop_plot.png", dod6_pop_plot, width=11, height=8)


test <- subset(long_dod6_data, Volunteer=="Volunteer 02")
test2 <- subset(test, Timepoint=="DoD")

