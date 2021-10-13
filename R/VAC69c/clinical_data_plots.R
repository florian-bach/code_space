library(dplyr)
library(ggplot2)
library(tidyr)


biochem <- read.csv("~/PhD/clinical_data/vac69c/biochem.csv", header=TRUE, stringsAsFactors = FALSE)

biochem$volunteer <- paste("v", substr(biochem$trial_number, 6, 8), sep="")

lousy_biochem_timepoints <- unique(biochem$timepoint)
# [1] "C_RESCREEN" "C_C2"       "C_C17"      "C_C814"     "C_C28DIAG"  "C_T1"       "C_T3"      
# [8] "C_T6"       "C_EV"       "C_C56"      "SCR"        "C_C96"

good_biochem_timepoints <- c("ReScreening", "Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "EV", "C56", "Screening", "C96")

biochem_timepoint_replacement <- setNames(good_biochem_timepoints, lousy_biochem_timepoints)
biochem$timepoint <- stringr::str_replace_all(biochem$timepoint, biochem_timepoint_replacement)

biochem$n_infection <- ifelse()

quant_biochem <- biochem %>%
  select(volunteer, timepoint, date_performed, sodium, potassium, urea, creatinine, bilirubin, alt, alkphos, albumin) %>%
  pivot_longer(cols = c(sodium, potassium, urea, creatinine, bilirubin, alt, alkphos, albumin), names_to = "Analyte", values_to = "Concentration")%>%
  filter(timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96"))


ggplot(quant_biochem, aes(x=factor(timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96")), y=Concentration, color=volunteer, group=volunteer))+
  #scale_fill_manual(values=volunteer_palette)+
  #scale_color_manual(values=volunteer_palette)+
  geom_line(size=0.9)+
  geom_point(fill="white", stroke=1, shape=21, size=0.9)+
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

alt_data <- filter(quant_biochem, Analyte=="alt")

ggplot(alt_data, aes(x=factor(timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96")), y=Concentration, color=volunteer, group=volunteer))+
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
