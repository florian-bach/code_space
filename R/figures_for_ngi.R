library(ggplot2)
library(colorspace)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)

### wassupppp ####

setwd("/Users/s1249052/PhD/presentations/NexGenImmunology_2020")


n_infection_palette <- c("#6666FF", "#ffbf00", "#FB027F")
names(n_infection_palette) <- c("First", "Second", "Third")

indie_color_scheme <- rev(c("#023fa5", "#7d87b9", "#bec1d4", "#d6bcc0", "#bb7784", "#8e063b", "#4a6fe3", "#8595e1", "#b5bbe3", "#e6afb9", "#e07b91", "#d33f6a", "#11c638", "#8dd593", "#c6dec7", "#ead3c6", "#f0b98d", "#ef9708", "#0fcfc0", "#9cded6", "#d5eae7", "#f3e1eb", "#f6c4e1", "#f79cd4",'#ff00ff'))
alan_color_scheme <- c("#ff0000", "#b00000", "#870000", "#550000", "#e4e400", "#baba00", "#878700", "#545400", "#00ff00", "#00b000", "#008700", "#005500", "#00ffff", "#00b0b0", "#008787", "#005555", "#b0b0ff", "#8484ff", "#4949ff", "#0000ff", "#ff00ff", "#b000b0", "#870087", "#550055", "#e4e4e4", "#bababa", "#878787", "#545454")
diana_color_scheme <- c("#9FA2FF", "#808080", "#8000FF", "#E54787", "#0080FF", "#4B0055", "#87DF6B", "#FF8000", "#009B95", "#EE000C", "#FFC800", "#66CCFF", "#F8A29E")

  
  # VAC63C lymphocytes ####

vac63c_lymph <- read.csv("VAC063_haem_all_sequenced_WBC_real_percent.csv")

vac63c_lymph <- dplyr::filter(vac63c_lymph, Leukocytes=="Lymphocytes")

vac63c_lymph <- select(vac63c_lymph, Volunteer_code, trial_number, timepoint, N_infection, Leukocytes, cell_counts)
vac63c_lymph <- dplyr::filter(vac63c_lymph, timepoint %in% c("C-1", "Diagnosis", "D+6"))
vac63c_lymph <- dplyr::filter(vac63c_lymph, N_infection %in% c("First", "Diagnosis", "Third"))
vac63c_lymph$timepoint <- gsub("D+6", "T6", vac63c_lymph$timepoint, fixed=T)
vac63c_lymph$timepoint <- gsub("C-1", "Baseline", vac63c_lymph$timepoint, fixed=T)

indie_lymph_vac63c <- ggplot(vac63c_lymph, aes(x=factor(timepoint, levels=c("Baseline", "Diagnosis", "T6")), y=cell_counts, group=trial_number))+
  geom_point(aes(color=vac63c_lymph$Volunteer_code))+
  geom_line(aes(color=vac63c_lymph$Volunteer_code))+
  theme_minimal()+
  scale_color_manual(values=diana_color_scheme)+
  xlab("Timepoint")+
  ylab(bquote('Lymphocytes ('*10^6~cells~'/ mL)'))+
  theme(legend.position="none",
        axis.text = element_text(size=16),
        axis.title = element_text(size=22),
        axis.title.x = element_blank())


# group_lymph_vac63c <- ggplot(vac63c_lymph, aes(x=factor(timepoint, levels=c("C-1", "Diagnosis", "D+6")), y=cell_counts, group=trial_number))+
#   geom_point(aes(color=vac63c_lymph$N_infection))+
#   geom_line(aes(color=vac63c_lymph$N_infection))+
#   theme_minimal()+
#   xlab("Timepoint")+
#   ylab(bquote('Lymphocytes ('*10^6~cells~'/ mL)'))+
#   labs(color="Infection")
#   


group_lymph_vac63c_box <- ggplot(vac63c_lymph, aes(x=factor(timepoint, levels=c("Baseline", "Diagnosis", "T6")), y=cell_counts, fill=N_infection))+
  geom_boxplot()+
  theme_minimal()+
  xlab("Timepoint")+
  ylab(bquote('Lymphocytes ('*10^6~cells~'/ mL)'))+
  labs(fill="Infection")+
  scale_fill_manual(values=n_infection_palette)+
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=22),
        axis.title.x = element_blank(),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18))




lgd <- get_legend(group_lymph_vac63c_box)
group_lymph_vac63c_box <- group_lymph_vac63c_box + theme(legend.position = "none")

# make common x axis title
timepoint.grob <- textGrob("Timepoint", 
                   gp=gpar(fontsize=22))

#add to plot
vac63c_lymphocytes_figure <- plot_grid(indie_lymph_vac63c, group_lymph_vac63c_box)
vac63c_lymphocytes_figure <- grid.arrange(arrangeGrob(vac63c_lymphocytes_figure, bottom = timepoint.grob))

vac63c_lymphocytes_figure <- plot_grid(vac63c_lymphocytes_figure, lgd, nrow = 1, rel_widths = c(10, 2))


ggsave("vac63c_lymphocytes_figure.png", vac63c_lymphocytes_figure, width=11, height=6)



# VAC63C PARASITAEMIA ####

vac63c_parasitaemia <- read.csv("VAC063_parasitaemias_all.csv", header=T)
long_vac63c_parasitaemia <- gather(vac63c_parasitaemia, vol_id, parasitaemia, colnames(vac63c_parasitaemia)[2:ncol(vac63c_parasitaemia)])

long_vac63c_parasitaemia$Volunteer <- substr(long_vac63c_parasitaemia$vol_id, 1, 5)
long_vac63c_parasitaemia$N_infection <- ifelse(grepl("First", long_vac63c_parasitaemia$vol_id)==T, "First",
                                               ifelse(grepl("Second", long_vac63c_parasitaemia$vol_id)==T, "Second", "Third"))


long_vac63c_parasitaemia$Volunteer <- gsub("X", "", long_vac63c_parasitaemia$Volunteer)
long_vac63c_parasitaemia$Volunteer <- gsub("_", "", long_vac63c_parasitaemia$Volunteer)

long_vac63c_parasitaemia <- long_vac63c_parasitaemia[!is.na(long_vac63c_parasitaemia$parasitaemia),]
long_vac63c_parasitaemia <- dplyr::filter(long_vac63c_parasitaemia, N_infection %in% c("First", "Third"))





vac63c_indie_paras <- ggplot(data=long_vac63c_parasitaemia, aes(x=Timepoint, y=parasitaemia, group=vol_id))+
  geom_point(aes(color=long_vac63c_parasitaemia$Volunteer))+
  geom_line(aes(color=long_vac63c_parasitaemia$Volunteer))+
  scale_y_log10()+
  theme_minimal()+
  scale_color_manual(values=diana_color_scheme)+
  ylab("Parasites / mL")+
  xlab("Days Post Infection")+
  theme(legend.position = "none",
        axis.text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.title = element_text(size=22))
       


vac63c_group_paras <- ggplot(data=long_vac63c_parasitaemia, aes(x=Timepoint, y=parasitaemia, group=vol_id))+
  geom_point(aes(color=long_vac63c_parasitaemia$N_infection))+
  geom_line(aes(color=long_vac63c_parasitaemia$N_infection))+
  scale_y_log10()+
  theme_minimal()+
  ylab("Parasites / mL")+
  xlab("Days Post Infection")+
  labs(color="Infection")+
  scale_color_manual(values=n_infection_palette)+
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size=20),
        axis.title = element_text(size=22),
        axis.title.x = element_blank(),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18))



paras_lgd <- get_legend(vac63c_group_paras)
vac63c_group_paras <- vac63c_group_paras+theme(legend.position = "none")      


# make common x axis title
days.grob <- grid::textGrob("Days Post Infection", 
                   gp=grid::gpar(fontsize=22))

#add to plot
vac63c_parasites_figure <- plot_grid(vac63c_indie_paras, vac63c_group_paras)
vac63c_parasites_figure <- grid.arrange(arrangeGrob(vac63c_parasites_figure, bottom = days.grob))

vac63c_parasites_figure <- plot_grid(vac63c_parasites_figure, paras_lgd, nrow = 1, rel_widths = c(10, 2))



#vac63c_parasites_figure <- plot_grid(vac63c_indie_paras, vac63c_group_paras, paras_lgd, nrow = 1, rel_widths = c(5, 5, 2))

ggsave("vac63c_parasites_figure.png", vac63c_parasites_figure, width=11, height=6)

#yellow


  