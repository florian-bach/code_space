library(tidyr)
library(dplyr)
library(xlsx)
library(ggplot2)
library(purrr)
library(data.table)
# library(mclust)
library(emmeans)

mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")
maternal_treatment_arms <- haven::read_dta("~/Library/CloudStorage/Box-Box/DP+SP study/Databases and preliminary findings/Final database used for analyses/DPSP treatment allocation_FINAL.dta")
epi_data <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/micdrop_metadata_for_immunology.csv")
nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")
msd_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/msd_results_v2.csv")
vaccines = c("Diptheria",     "Measles" ,      "Mumps",         "Pertussis",     "Polio",
             "Rotavirus" ,    "Rubella",       "Tetanus", "Pneumo.1.4.14")

long_msd <- msd_data%>%
  mutate(sample=paste(SubjectID, "_", "tp", TimePt, sep=""))%>%
  mutate(id=SubjectID, timepoint=paste(TimePt, "weeks"))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))%>%
  select(-SubjectID, -TimePt)%>%
  pivot_longer(cols=-c(sample, id, timepoint), names_to = "antigen", values_to = "titer")%>%
  mutate(Ig_class=ifelse(grepl("*IgA$", antigen), "IgA", "IgG"))%>%
  filter(!is.na(titer))


# additional data from the moment of sampling
additional_clinical_data <- nulisa_data%>%
  distinct(id, sample, date, ageinwks, gender_categorical, mstatus, qPCRparsdens)

antibodies_and_epi <- epi_data%>%
  inner_join(., long_msd, by="id")%>%
  left_join(additional_clinical_data, by="sample")
# figure 1 timeline of antibodies, maternal & birth characteristics characteristics ####
## timeline of antibodies ####
## birth characteristics ####
# "GAcomputed" : nothing
# "SGA": signal for SARS-CoV2
# "LBWdich": not really reliable but generally less than NBW


## gravidity???
#nothing
antibodies_and_epi%>%
  filter(timepoint=="8 weeks", !is.na(GAcomputed))%>%
  filter(Ig_class=="IgG")%>%
  ggplot(., aes(x=GAcomputed, y=titer))+
  geom_point()+
  geom_smooth(method='lm')+
  ggpubr::stat_cor(method = "spearman")+
  scale_y_log10()+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()

# signal for SARS-CoV2
antibodies_and_epi%>%
  filter(timepoint=="24 weeks", !is.na(SGA))%>%
  filter(Ig_class=="IgG")%>%
  ggplot(., aes(x=factor(SGA), y=titer))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ggpubr::stat_compare_means(vjust = 0.5, comparisons = list(c("0", "1")))+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()

#have only 4 so not reliable, but generally bad, especially at 24 weeks
antibodies_and_epi%>%
  filter(timepoint=="24 weeks", !is.na(LBWdich))%>%
  filter(Ig_class=="IgG")%>%
  ggplot(., aes(x=factor(LBWdich), y=titer))+
  geom_boxplot(outliers = F)+
  geom_point()+
  scale_y_log10()+
  ggpubr::stat_compare_means(vjust = 0.5, comparisons = list(c("0", "1")))+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()
## maternal characteristics ####
### rogerson / anyHP ####
# nothing at 8 weeks
antibodies_and_epi%>%
  filter(timepoint=="8 weeks", !is.na(rogerson))%>%
  filter(Ig_class=="IgG")%>%
  ggplot(., aes(x=factor(rogerson), y=titer))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ggpubr::stat_compare_means(vjust = 0.5, comparisons = list(c("0", "1"),
                                                             c("0", "4")))+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()

antibodies_and_epi%>%
  filter(timepoint=="8 weeks", !is.na(anyHP))%>%
  filter(Ig_class=="IgG")%>%
  ggplot(., aes(x=factor(anyHP), y=titer))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ylab("concentration at 8 weeks")+
  ggpubr::stat_compare_means(vjust = 0.5, comparisons = list(c("0", "1")), p.adjust.method = "BH")+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()
### maternal treatment group ####
#### # FluA and Ev71 signal
antibodies_and_epi%>%
  filter(timepoint=="8 weeks", !is.na(mom_rx))%>%
  filter(Ig_class=="IgG")%>%
  ggplot(., aes(x=factor(mom_rx), y=titer))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ylab("concentration at 8 weeks")+
  ggpubr::stat_compare_means(vjust = 0.5, comparisons = list(c("DP", "SP"),
                                                             c("DPSP", "SP")))+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()

# need to check maternal infection intensity 

#nothing
antibodies_and_epi%>%
  filter(timepoint=="8 weeks", !is.na(mom_rx))%>%
  filter(Ig_class=="IgA")%>%
  ggplot(., aes(x=factor(mom_rx), y=titer))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ggpubr::stat_compare_means(vjust = 0.5, comparisons = list(c("DP", "SP"),                                                         c("DPSP", "SP")))+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()






antibodies_and_epi%>%
  filter(timepoint=="52 weeks", !is.na(total_n_para_12))%>%
  filter(Ig_class=="IgG")%>%
  ggplot(., aes(x=factor(total_n_para_12), y=titer))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  # ggpubr::stat_compare_means(vjust = 0.5, comparisons = list(c("0", "1"),
  #                                                            c("0", "4")))+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()


antibodies_and_epi%>%
  filter(timepoint=="24 weeks", !is.na(qPCRparsdens))%>%
  filter(Ig_class=="IgG")%>%
  ggplot(., aes(x=log10(qPCRparsdens+0.001), y=titer))+
  geom_point()+
  geom_smooth(method='lm')+
  ggpubr::stat_cor(method = "spearman")+
  scale_y_log10()+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()


# anemia ####