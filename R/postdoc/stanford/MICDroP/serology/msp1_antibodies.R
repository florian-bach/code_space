library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(emmeans)

treatment_palette <-  c("darkred", "darkblue")
names(treatment_palette) <-c("Placebo", "DP 1 year")

msp_luminex <- readxl::read_excel("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/msp_data/slim_msd_data.xlsx")

plate_map <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/final_plate_map.csv")
mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")
raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/MICDroP Data/MICDROP expanded database through April 30th 2026.dta")

plate_map <- plate_map%>%
  mutate(well=paste(plate.row, plate.column, sep=""),
         location=paste(dane.plate, well, sep=""))

slim_plate_map <- plate_map%>%
  select(location, id, date, ageinwks)%>%
  mutate("date"=as.Date(lubridate::parse_date_time(.$date, orders="%m/%d/%y")))

all_samples <- distinct(plate_map, id, date, ageinwks)%>%
  mutate("date"=as.Date(lubridate::parse_date_time(.$date, orders="%m/%d/%y")))

# add epi data ####
## add visit-level metadata ####
metadata_columns <- c("id", "date", "mstatus", "pardens", "qPCRparsdens", "fever", "febrile", "gender")

dobs <- raw_data%>%
  distinct(id, dob)%>%
  filter(!is.na(dob))

visit_metadata <- raw_data%>%
  mutate(date=as.Date(lubridate::parse_date_time(date, "%y/%m/%d")))%>%
  select(all_of(metadata_columns))%>%
  right_join(., all_samples, by=c("id", "date"))%>%
  mutate(dob2=dobs$dob[match(id, dobs$id)])%>%
  mutate(real_dob=as.Date(lubridate::parse_date_time(dob2, "%y/%m/%d")))%>%
  # mutate(ageinwks=round((date-real_dob)/7))%>%
  mutate(timepoint_num=as.numeric(ageinwks))%>%
  mutate(timepoint_num=case_when(timepoint_num==9~8,
                                 timepoint_num==53~52,
                                 timepoint_num==105~104,
                                 .default=timepoint_num))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(gender_categorical=if_else(gender==1, "male", "female"))


## static meta data (infection, birth outcomes) ####
static_metadata <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/micdrop_metadata_for_immunology.csv")

long_luminex <- msp_luminex%>%
  mutate(location=paste(Plate, Well, sep=""))%>%
  pivot_longer(cols=c(MSP1), names_to = "antigen", values_to = "MFI")%>%
  left_join(., slim_plate_map, by="location")%>%
  left_join(., visit_metadata, by=c("id", "date"))%>%
  mutate(date=as.Date(lubridate::parse_date_time(date, "%y/%m/%d")))%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  mutate(log_mfi=log10(MFI+1),
         timepointw=paste(timepoint_num, "weeks"),
         timepoint=factor(timepointw, levels=c("8 weeks", "24 weeks", "52 weeks", "104 weeks")))%>%
  filter(!is.na(log_mfi), is.finite(log_mfi))%>%
  left_join(., static_metadata, by=c("id"))%>%
  mutate(antigen=gsub("_", " ", antigen),
         id_cat=as.character(id))%>%
  # there's a couple of mispipetted samples that have NAs for timepoint
  filter(!is.na(timepoint_num))%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))


time_treatment_plot <- long_luminex %>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  filter(!is.na(MFI))%>%
  ggplot(., aes(x=treatmentarm, y=MFI, fill=treatmentarm))+
  geom_violin()+
  geom_boxplot(width=0.2, position = position_dodge(width=0.9),  outliers = F, color="grey")+
  scale_y_log10()+
  ggtitle("MSP1")+
  facet_wrap(~timepoint, nrow=1)+
  ggpubr::stat_compare_means(vjust=-0.2)+
  scale_fill_manual(values=rev(c("firebrick", "darkblue")))+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/msp_data/figures/time_treatment_plot.png", time_treatment_plot, width=8, height=5, dpi=444)
ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/Pf_antibodies/MSP1/figures/time_treatment_plot.png", time_treatment_plot, width=8, height=5, dpi=444)
ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/Pf_antibodies/MSP1/figures/time_treatment_plot.pdf", time_treatment_plot, width=8, height=5, dpi=444)



 pardens_plot <- long_luminex %>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  filter(!is.na(MFI), timepoint_num!=8)%>%
  ggplot(., aes(x=qPCRparsdens+0.01, y=MFI, color=treatmentarm))+
  geom_point()+
  scale_y_log10()+
  scale_x_log10()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman", size=5, label.y = c(1, 1.3))+
  xlab("qPCR parasite density")+
  ggtitle("parasite density is correlated with MSP1 antibodies\n(8 week timepoint excluded)")+
  # facet_wrap(~timepoint, scales="free", labeller = label_wrap_gen(width = 10))+
  # viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  scale_color_manual(values=rev(c("firebrick", "darkblue")))+
  theme_minimal()+
  theme(legend.position="bottom",
        legend.direction = "horizontal")
ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/msp_data/figures/pardens_plot.png", pardens_plot, width=5, height=5, dpi=444)

total_n_para_12_plot <- long_luminex %>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  filter(!is.na(MFI), timepoint_num==52)%>%
  ggplot(., aes(x=factor(total_n_para_12), y=MFI, fill=treatmentarm))+
  geom_boxplot()+
  scale_y_log10()+
  xlab("number of parastiemic months in the first tyear of life")+
  ylab("MFI at 52 weeks")+
  # ggtitle("parasite density is correlated with MSP1 antibodies\n(8 week timepoint excluded)")+
  facet_wrap(~treatmentarm, scales="free", labeller = label_wrap_gen(width = 10))+
  # viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  scale_fill_manual(values=rev(c("firebrick", "darkblue")))+
  theme_minimal()+
  theme(legend.position="bottom",
        legend.direction = "horizontal")
ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/msp_data/figures/total_n_para_12_plot.png", total_n_para_12_plot, width=8, height=5, dpi=444)



write.csv(long_luminex, "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/msp_data/micdrop_msd_luminex_data.csv")
