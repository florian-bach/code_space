# preamble ####
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)


vac69c_colours <-list("v07" = "#F0E442",
                      "v11" = "#D33F6A",
                      "v21" = "#E99A2C",
                      "v34" = "red", 
                      "v35" = "#6C3461",
                      "v37" = "#5E9B8A",
                      "v39" = "#FF08E8",
                      "v40" = "#61E160",
                      "v42" = "#6488EA",
                      "v44" = "darkblue")


vac69c_colours <- unlist(unname(vac69c_colours))
names(vac69c_colours) <- names(vac69c_colours)

n_infection_palette <- c(rgb(5,50,80, maxColorValue = 255),
                         rgb(250, 100, 0, maxColorValue = 255),
                         rgb(40,210,250, maxColorValue = 255))

names(n_infection_palette) <- c("First", "Second", "Third")

#show_col(vac69c_colours)



# biochem ####


biochem <- read.csv("~/PhD/clinical_data/vac69c/biochem.csv", header=TRUE, stringsAsFactors = FALSE)

biochem$volunteer <- paste("v", substr(biochem$trial_number, 6, 8), sep="")

lousy_biochem_timepoints <- unique(biochem$timepoint)
# [1] "C_RESCREEN" "C_C2"       "C_C17"      "C_C814"     "C_C28DIAG"  "C_T1"       "C_T3"      
# [8] "C_T6"       "C_EV"       "C_C56"      "SCR"        "C_C96"

good_biochem_timepoints <- c("ReScreening", "Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "EV", "C56", "Screening", "C96")

biochem_timepoint_replacement <- setNames(good_biochem_timepoints, lousy_biochem_timepoints)
biochem$timepoint <- stringr::str_replace_all(biochem$timepoint, biochem_timepoint_replacement)

biochem$n_infection <- ifelse(biochem$volunteer == "v07", "Third",
                              ifelse(biochem$volunteer %in% c("v11", "v21"), "Second", "First"))

quant_biochem <- biochem %>%
  select(volunteer, timepoint, n_infection, date_performed, sodium, potassium, urea, creatinine, bilirubin, alt, alkphos, albumin) %>%
  pivot_longer(cols = c(sodium, potassium, urea, creatinine, bilirubin, alt, alkphos, albumin), names_to = "Analyte", values_to = "Concentration")%>%
  filter(timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96"))


indie_biochem_plot <- ggplot(quant_biochem, aes(x=factor(timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96")), y=Concentration, group=volunteer))+
  scale_color_manual(values=vac69c_colours)+
  geom_line(aes(color=volunteer), size=0.9)+
  geom_point(aes(color=volunteer), fill="white", stroke=1, shape=21, size=0.9)+
  facet_wrap(~Analyte, scales="free")+
  theme_minimal()+
  xlab("Timepoint")+
  ylab("Concentration")+
  #scale_y_continuous(label=scales::comma)+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8), 
        axis.title.x = element_blank(),
        #legend.position = "none",
        strip.text = element_text(size=10))

ggsave("~/PhD/clinical_data/vac69c/figures/indie_biochem.png", indie_biochem_plot, height=6, width=8, bg="white")



group_biochem_plot <- ggplot(quant_biochem, aes(x=factor(timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96")), y=Concentration, group=volunteer))+
  scale_color_manual(values=n_infection_palette)+
  geom_line(aes(color=n_infection), size=0.9)+
  geom_point(aes(color=n_infection), fill="white", stroke=1, shape=21, size=0.9)+
  facet_wrap(~Analyte, scales="free")+
  theme_minimal()+
  xlab("Timepoint")+
  ylab("Concentration")+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8), 
        axis.title.x = element_blank(),
        strip.text = element_text(size=10))

ggsave("~/PhD/clinical_data/vac69c/figures/n_infection_biochem.png", indie_biochem_plot, height=6, width=8, bg="white")


alt_data <- filter(quant_biochem, Analyte=="alt")

indie_alt_data_plot <- ggplot(alt_data, aes(x=factor(timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96")), y=Concentration, color=volunteer, group=volunteer))+
  #scale_fill_manual(values=volunteer_palette)+
  #scale_color_manual(values=volunteer_palette)+
  geom_line(size=0.9)+
  geom_point(fill="white", stroke=1, shape=21, size=0.9)+
  theme_minimal()+
  geom_hline(yintercept = 35, color="red", linetype="dashed")+
  geom_hline(yintercept = 45, color="blue", linetype="dashed")+
  ggtitle("ALT during VAC69c")+
  xlab("Timepoint")+
  ylab("ALT (IU / L)\n")+
  scale_y_log10()+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8), 
        axis.title.x = element_blank(),
        legend.title =  element_blank(),
        strip.text = element_text(size=10))


group_alt_data_plot <- ggplot(alt_data, aes(x=factor(timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96")), y=Concentration, color=n_infection, group=volunteer))+
  scale_color_manual(values=n_infection_palette)+
  geom_line(size=0.9)+
  geom_point(fill="white", stroke=1, shape=21, size=0.9)+
  theme_minimal()+
  geom_hline(yintercept = 35, color="red", linetype="dashed")+
  geom_hline(yintercept = 45, color="blue", linetype="dashed")+
  ggtitle("ALT during VAC69c")+
  xlab("Timepoint")+
  ylab("ALT (IU / L)\n")+
  scale_y_log10()+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8), 
        axis.title.x = element_blank(),
        legend.title =  element_blank(),
        strip.text = element_text(size=10))


combo_alt_plot <- cowplot::plot_grid(indie_alt_data_plot, group_alt_data_plot, nrow = 1)

ggsave("~/PhD/clinical_data/vac69c/figures/combo_alt_plot.png", combo_alt_plot, height=3, width=8, bg="white")


ae_biochem <- biochem %>%
  select(volunteer, timepoint, date_performed, sodium_ae, potassium_ae, urea_ae, creatinine_ae, bilirubin_ae, alt_ae, alkphos_ae, albumin_ae) %>%
  pivot_longer(cols = c(sodium_ae, potassium_ae, urea_ae, creatinine_ae, bilirubin_ae, alt_ae, alkphos_ae, albumin_ae), names_to = "Analyte", values_to = "Severity")%>%
  filter(timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96"))

ae_biochem$timepoint <- factor(ae_biochem$timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96"))

ggplot(ae_biochem,  aes(x="", fill=factor(Severity)))+
  geom_bar(stat="count", position ="stack")+
  scale_fill_manual(values =  list("0"="lightgrey", "1"="#FACA0F", "2"="chocolate1", "3"="red"))+
  facet_grid(timepoint~Analyte, switch = "y")+
  coord_polar(theta = "y")+
  ggtitle("Adverse Events\n")+
  #scale_y_continuous(limits = c(0,8), breaks = seq(0, 8, by=2))+
  theme_void()+
  theme(plot.title = element_text(hjust=0.5, vjust=0, size=10),
        panel.grid.minor = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        axis.title = element_blank())+
  guides(fill=guide_legend(title="Severity",
                           override.aes = list(size = 0.1),
                           keywidth = 0.5,
                           keyheight = 0.5))

# haem ####

haem <- read.csv("~/PhD/clinical_data/vac69c/haem.csv", header=TRUE, stringsAsFactors = FALSE)

haem$volunteer <- paste("v", substr(haem$trial_number, 6, 8), sep="")

lousy_haem_timepoints <- unique(haem$timepoint)
# [1] "C_RESCREEN" "C_C2"       "C_C17"      "C_C814"     "C_C28DIAG"  "C_T1"       "C_T3"      
# [8] "C_T6"       C_C56"  "C_EV"    "SCR"        "C_C96"

good_haem_timepoints <- c("ReScreening", "Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "EV", "Screening", "C96")

haem_timepoint_replacement <- setNames(good_haem_timepoints, lousy_haem_timepoints)
haem$timepoint <- stringr::str_replace_all(haem$timepoint, haem_timepoint_replacement)

haem$n_infection <- ifelse(haem$volunteer == "v07", "Third",
                              ifelse(haem$volunteer %in% c("v11", "v21"), "Second", "First"))

quant_haem <- haem %>%
  mutate("haemoglobin"=haemaglobin)%>%
  select(volunteer, timepoint, n_infection, date_performed, haemoglobin, wbc, platelets, haematocrit, red_cell_count, neutrophils, lymphocytes, monocytes, eosinophils, basophils) %>%
  pivot_longer(cols = c(haemoglobin, wbc, platelets, haematocrit, red_cell_count, neutrophils, lymphocytes, monocytes, eosinophils, basophils), names_to = "Cell", values_to = "Frequency")%>%
  filter(timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96"))
  


indie_haem <- ggplot(quant_haem, aes(x=factor(timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96")), y=Frequency, group=volunteer))+
  scale_color_manual(values=vac69c_colours)+
  geom_line(aes(color=volunteer), size=0.9)+
  geom_point(aes(color=volunteer), fill="white", stroke=1, shape=21, size=0.9)+
  facet_wrap(~Cell, scales="free", ncol = 5)+
  theme_minimal()+
  xlab("Timepoint")+
  ylab("Frequency")+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8), 
        axis.title.x = element_blank(),
        #legend.position = "none",
        strip.text = element_text(size=10))

ggsave("~/PhD/clinical_data/vac69c/figures/indie_haem.png", indie_haem, height=5, width=12, bg="white")


lymph_data <- filter(quant_haem, Cell %in% c("platelets", "lymphocytes"))

lymph_data$n_infection <- factor(lymph_data$n_infection, levels=c("First", "Second", "Third"))

indie_thrombo_lymph <- ggplot(lymph_data, aes(x=factor(timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96")), y=Frequency, group=volunteer))+
  scale_color_manual(values=vac69c_colours)+
  geom_line(aes(color=volunteer), size=0.9)+
  geom_point(aes(color=volunteer), fill="white", stroke=1, shape=21, size=0.9)+
  theme_minimal()+
  facet_wrap(~Cell, scales="free")+
  xlab("Timepoint")+
  ylab("Frequency")+
  #scale_y_continuous(label=scales::comma)+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8), 
        axis.title.x = element_blank(),
        #legend.position = "none",
        strip.text = element_text(size=10))



group_thrombo_lymph <- ggplot(lymph_data, aes(x=factor(timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96")), y=Frequency, group=volunteer, color=n_infection))+
  scale_color_manual(values=n_infection_palette)+
  geom_line(size=0.9)+
  geom_point(fill="white", stroke=1, shape=21, size=0.9)+
  theme_minimal()+
  facet_wrap(~Cell, scales="free")+
  xlab("Timepoint")+
  ylab("Frequency")+
  #scale_y_continuous(label=scales::comma)+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8), 
        axis.title.x = element_blank(),
        #legend.position = "none",
        strip.text = element_text(size=10))


combo_haem_plot <- cowplot::plot_grid(indie_thrombo_lymph, group_thrombo_lymph, nrow=2)

ggsave("~/PhD/clinical_data/vac69c/figures/combo_haem_plot.png", combo_haem_plot, height=6, width=6, bg="white")

ae_haem <- haem %>%
  select(volunteer, timepoint, date_performed, haemaglobin_ae, wbc_ae, platelets_ae, neutrophils_ae, lymphocytes_ae) %>%
  pivot_longer(cols = c(haemaglobin_ae, wbc_ae, platelets_ae, neutrophils_ae, lymphocytes_ae), names_to = "Analyte", values_to = "Severity")%>%
  filter(timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96"))


ae_haem$timepoint <- factor(ae_haem$timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96"))

ggplot(ae_haem,  aes(x="", fill=factor(Severity)))+
  geom_bar(stat="count", position ="stack")+
  scale_fill_manual(values =  list("0"="lightgrey", "1"="#FACA0F", "2"="chocolate1", "3"="red"))+
  facet_grid(timepoint~Analyte, switch = "y")+
  coord_polar(theta = "y")+
  ggtitle("Adverse Events\n")+
  #scale_y_continuous(limits = c(0,8), breaks = seq(0, 8, by=2))+
  theme_void()+
  theme(plot.title = element_text(hjust=0.5, vjust=0, size=10),
        panel.grid.minor = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        axis.title = element_blank())+
  guides(fill=guide_legend(title="Severity",
                           override.aes = list(size = 0.1),
                           keywidth = 0.5,
                           keyheight = 0.5))


# symptoms ####