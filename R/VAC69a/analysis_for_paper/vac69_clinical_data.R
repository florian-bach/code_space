library(ggplot2)
library(tidyr)
library(dplyr)
library(gtools)
library(cowplot)


#setwd("C:/Users/Florian/PhD/cytof/vac69a/clinical_data")
setwd("/Users/s1249052/PhD/clinical_data/vac69a")
##########   symptoms findings   #################




data <- read.csv("symptoms.csv", header=T)


colnames(data)[1] <- "Volunteer"

data$Volunteer <- paste("Volunteer", substr(data$Volunteer, 6, 7), sep=' ')

data$timepoint <- gsub("_am", "", data$timepoint, fixed=T)
data$timepoint <- gsub("_pm", ".5", data$timepoint, fixed=T)
data$timepoint <- substr(data$timepoint, nchar(as.character(data$timepoint))-3, nchar(as.character(data$timepoint)))
data$timepoint <- gsub("_", "", data$timepoint, fixed=T)
data$timepoint <- paste("C+", data$timepoint, sep='')

data$timepoint <- gsub("C+ep", "EP+", data$timepoint, fixed=T)
data$timepoint <- gsub("C+c", "C+", data$timepoint, fixed=T)
data$timepoint <- gsub("C+6", "T+6", data$timepoint, fixed=T)
data$timepoint <- gsub("C+0", "C+", data$timepoint, fixed=T)


data$pyrexia_temp <- NULL
data$X <- NULL


long_data <- gather(data, Symptom, Severity, colnames(data)[3:16]) 

specific_levels <- mixedsort(unique(data$timepoint))

better_levels <- append(specific_levels, specific_levels[34:36], after=31)

better_levels <- better_levels[-c(37:39)]

for(i in unique(long_data$Volunteer)){
 
 sub_set <- filter(long_data, Volunteer == i)

 ifelse(unique(sub_set$Volunteer) %in% c("Volunteer 02","Volunteer 06"), assign("result", element_text(size=12)), assign("result", element_blank()))
 ifelse(unique(sub_set$Volunteer) %in% c("Volunteer 05","Volunteer 09"), assign("result1", "right"), assign("result1", "left"))
 
 
 
 assign(paste("V", substr(i, 11,12), sep=''),
     

 ggplot(sub_set, aes(x=factor(timepoint, levels=better_levels), y=Symptom))+
 geom_tile(aes(fill=as.character(Severity)), color="white")+
 scale_fill_manual(values= c("grey", "yellow", "orange", "red"))+
 scale_y_discrete(position = result1)+
 ggtitle(unique(sub_set$Volunteer))+
 theme(panel.border = element_blank(),
       axis.text.y.left = result,
       axis.line.y.left = element_blank(),
       axis.line.y.right = element_blank(),
       axis.ticks.y = element_blank(),
       axis.title.y = element_blank(),
       axis.title.x = element_blank(),
       axis.text.y.right = element_text(size = 12),
       strip.background = element_blank(),
       legend.position = "none",
       plot.title = element_text(size = 16, hjust = 0.5),
       axis.text.x = element_text(angle = 60, hjust = 1, size=10),
       plot.margin = unit(c(1,0,1,0), "cm"))
 )
}




(symptom_plot <- plot_grid(V02, V03, V05, V06, V07, V09, ncol = 3, rel_widths = c(1,1,1)))

ggsave("symptom_plot.png", symptom_plot, width=19, height=10)




##########   haematological findings   #################



data <- read.csv("haem.csv")

colnames(data)[1] <- "Volunteer"

data$Volunteer <- paste("Volunteer", substr(data$Volunteer, 6, 7), sep=' ')


short <- filter(data, timepoint %in% c("_C_1", "_C28", "_T6", "_C90", "_C1_7", "_EP"))

short2 <- short[,-c(seq(5,13, by=2), 16)]

short2$timepoint <- as.character(short2$timepoint)

short2$timepoint[short2$timepoint=="_C_1"] <- "C-1"
short2$timepoint[short2$timepoint=="_C28"] <- "C+28"
short2$timepoint[short2$timepoint=="_T6"] <- "T+6"
short2$timepoint[short2$timepoint=="_C90"] <- "C+90"
short2$timepoint[short2$timepoint=="_EP"] <- "EP+1"
short2$timepoint[short2$timepoint=="_C1_7"] <- "C+7"



long_short <- gather(short2, Marker, Value, colnames(short2)[4:10])

long_short$Marker[long_short$Marker=="haemaglobin"] <- "haemoglobin"

my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")


ggplot(long_short, aes(x=factor(timepoint, levels=c("C-1", "C+7", "EP+1", "T+6", "C+28","C+90")), y=Value, group=Volunteer))+
 geom_point(aes(colour=long_short$Volunteer))+
 geom_line(aes(colour=long_short$Volunteer))+
 facet_wrap(~Marker, scales="free")+
 scale_colour_manual(values=my_paired_palette)+
 theme_bw()+
 theme(axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=16),
    strip.text = element_text(size=16, face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1, size=12),
    axis.text.y = element_text(size=12)
 )


###########       extra figure for T cells     ##########

data <- read.csv("t_cell_counts_vac69a.csv")

data$T.cell.count <- data$T.cell.count*1000
data$Volunteer <- paste("Volunteer 0", data$Volunteer, sep='')

ggplot(data, aes(x=factor(timepoint, levels=c("C-1", "DoD", "T+6")), y=T.cell.count, group=Volunteer))+
  geom_point(aes(colour=as.character(data$Volunteer)), size=3)+
  geom_line(aes(colour=as.character(data$Volunteer)), size=1.5)+
  scale_colour_manual(values=my_paired_palette)+
  theme_bw()+
  ylim(0,1750)+
  ylab(expression(paste("T Cells / ", mu*"L")))+
  theme(axis.title.x = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(size=20),
        strip.text = element_text(size=20, face = "bold"),
        axis.text.x = element_text(angle = 60, hjust = 1, size=14),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=20),
        legend.position = "bottom",
        legend.title = element_blank())

  ggsave("T_Cell_Count.png", height=8.5, width=7)




ggsave("haematological.png", width=10, height=10)

#######################     clinical chemistry    ########################



data <- read.csv("biochem.csv", header=T)

data_no_ae <- select(data, -c(colnames(data)[grep("_ae",colnames(data) ,fixed=T)]))

long_data <- gather(data_no_ae, biochem, value, colnames(data_no_ae)[13:20])

long_data <- filter(long_data, long_data$flo_timepoint!="extra")

alt <- subset(long_data, long_data$timepoint=="_T6"&long_data$biochem=="alt")

data <- read.csv("biochem.csv")

colnames(data)[1] <- "Volunteer"

data$Volunteer <- paste("Volunteer", substr(data$Volunteer, 6, 7), sep=' ')


short <- filter(data, timepoint %in% c("_C_1", "_C28", "_T6", "_C90", "_C1_7", "_EP"))

short2 <- short[,-c(seq(5,19, by=2), 20:24)]

short2$timepoint <- as.character(short2$timepoint)

short2$timepoint[short2$timepoint=="_C_1"] <- "C-1"
short2$timepoint[short2$timepoint=="_C28"] <- "C+28"
short2$timepoint[short2$timepoint=="_T6"] <- "T+6"
short2$timepoint[short2$timepoint=="_C90"] <- "C+90"
short2$timepoint[short2$timepoint=="_EP"] <- "EP+1"
short2$timepoint[short2$timepoint=="_C1_7"] <- "C+7"



long_short <- gather(short2, Marker, Value, colnames(short2)[4:11])



my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")


ggplot(long_short, aes(x=factor(timepoint, levels=c("C-1", "C+7", "EP+1", "T+6", "C+28","C+90")), y=Value, group=Volunteer))+
 geom_point(aes(colour=long_short$Volunteer))+
 geom_line(aes(colour=long_short$Volunteer))+
 facet_wrap(~Marker, scales="free")+
 scale_colour_manual(values=my_paired_palette)+
 theme_bw()+
 theme(axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=16),
    strip.text = element_text(size=16, face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1, size=12),
    axis.text.y = element_text(size=12)
 )



ggsave("biochem.png", width=10, height=10)



########       blood counts        ########


data <- read.csv("blood_counts.csv")


data$Volunteer <- paste("Volunteer", data$Volunteer, sep=' ')



short2 <- data
short2$timepoint <- as.character(short2$timepoint)

short2$timepoint[short2$timepoint=="_C_1"] <- "C-1"
short2$timepoint[short2$timepoint=="_C28"] <- "C+28"
short2$timepoint[short2$timepoint=="_T6"] <- "T+6"
short2$timepoint[short2$timepoint=="_C90"] <- "C+90"
short2$timepoint[short2$timepoint=="_EP"] <- "EP+1"
short2$timepoint[short2$timepoint=="_C1_7"] <- "C+7"


long_short <- gather(short2, Cell_Type, Count, colnames(short2)[c(2:4, 6)])


ggplot(long_short, aes(x=factor(timepoint, levels=c("C-1", "EP+1", "T+6")), y=Count, group=Volunteer))+
  geom_point(aes(colour=as.character(long_short$Volunteer)), size=3)+
  geom_line(aes(colour=as.character(long_short$Volunteer)), size=1.5)+
  scale_colour_manual(values=my_paired_palette)+
  theme_bw()+
  facet_wrap(~Cell_Type, scales="free")
  #ylab(expression(paste("T Cells / ", mu*"L")))+
  theme(axis.title.x = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(size=20),
        strip.text = element_text(size=20, face = "bold"),
        axis.text.x = element_text(angle = 60, hjust = 1, size=14),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=20),
        legend.position = "bottom",
        legend.title = element_blank())

ggsave("T_Cell_Count.png", height=8.5, width=7)




############        vac 63 lymph counts        ############


setwd("C:/Users/Florian/PhD/cytof/vac63c")

data <- read.csv("lymph_stats.csv")

colnames(data)[1:3] <- c("T_of_lymph", "B_of_Lymph", "T_of_all")
data$B_of_all <- data$B_of_Lymph/data$T_of_lymph*data$T_of_all

long_data <- gather(data, Cells, Frequency, c("T_of_lymph", "B_of_Lymph", "T_of_all", "B_of_all"))

lymph_of_all_only <- long_data[long_data$Cells%in%c("T_of_all", "B_of_all"),]

ggplot(lymph_of_all_only, aes(x=factor(Timepoint, levels=c("C-1", "C+9", "C+10", "C+14", "C+16", "D+6")), y=Frequency, group=Volunteer))+
  geom_point()+
  geom_line()+
  facet_grid(Cells~Volunteer)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=60, hjust=1))
