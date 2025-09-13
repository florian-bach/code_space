library(tidyr)
library(dplyr)
library(xlsx)
library(ggplot2)
library(purrr)
library(emmeans)


nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")
mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")
maternal_treatment_arms <- haven::read_dta("~/Library/CloudStorage/Box-Box/DP+SP study/Databases and preliminary findings/Final database used for analyses/DPSP treatment allocation_FINAL.dta")
msd_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/batch_one.csv")
cmv_data <-  read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/CMV_ELISA_Uganda.csv", skip = 1)


vaccines = c("Diptheria",     "Measles" ,      "Mumps",         "Pertussis",     "Polio",
             "Rotavirus" ,    "Rubella",       "Tetanus", "Pneumo.1.4.14")

nulisa_data <- nulisa_data%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))%>%
  mutate(CMV=cmv_data$Diagnosis[match(id, cmv_data$id)],
         id_cat=paste(id),
         timepoint=factor(timepoint, levels = c("8 weeks", "24 weeks", "52 weeks")))



metadata_columns <- c("id", "anyDP", "treatmentarm",  "dob", "date", "ageinwks", "gender_categorical", "mstatus", "qPCRparsdens", "visittype", "fever", "febrile", "rogerson", "GAcomputed", "gi", "SGA", "qPCRdich", "mqPCRparsdens")


long_msd <- msd_data%>%
  mutate(sample=paste(SubjectID, "_", "tp", TimePt, sep=""))%>%
  mutate(id=SubjectID, timepoint=paste(TimePt, "weeks"))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))%>%
  select(-SubjectID, -TimePt)%>%
  pivot_longer(cols=-c(sample, id, timepoint), names_to = "antigen", values_to = "titer")

long_msd <- long_msd%>%
  mutate(CMV=cmv_data$Diagnosis[match(id, cmv_data$id)],
         id_cat=paste(id),
         timepoint=factor(timepoint, levels = c("8 weeks", "24 weeks", "52 weeks")))


cmv_elisa_purrf <- nulisa_data %>%
  filter(!is.na(CMV), mstatus==0, timepoint %in% c("8 weeks", "24 weeks", "52 weeks"))%>%
  group_by(targetName)%>%
  nest() %>%
  mutate(time_model=purrr::map(data, ~lme4::lmer(conc~timepoint*CMV+log_qpcr+gender+hbs+(1|id_cat), data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ CMV | timepoint)))%>%
  mutate(emm2=map(time_model, ~emmeans(., specs = pairwise ~ timepoint | CMV)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  mutate("8 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("24 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[4])) %>%
  mutate("52 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[7])) %>%
  pivot_longer(cols=ends_with("weeks"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

cmv_sigs <- cmv_elisa_purrf%>%
  filter(padj < 0.1)%>%
  select(targetName, contrast, padj)


nulisa_data%>%
  mutate(CMV=if_else(is.na(CMV), "UNKNOWN", CMV))%>%
  mutate(CMV=factor(CMV, levels=c("Negative", "Positive", "UNKNOWN")))%>%
  filter(timepoint %in% c("8 weeks", "24 weeks", "52 weeks"))%>%
  mutate(timepoint=factor(timepoint, levels = c("8 weeks", "24 weeks", "52 weeks")))%>%
  filter(targetName %in% cmv_sigs$targetName)%>%
  ggplot(., aes(x=timepoint, y=conc, fill=CMV))+
  geom_boxplot()+
  facet_wrap(~targetName+anyHP, scales="free")+
  scale_fill_viridis_d(option = "F")+
  theme_minimal()+
  theme(axis.title = element_blank())
  
  
cmv_msd_purrf <- long_msdanyHPcmv_msd_purrf <- long_msd%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(titer~timepoint*CMV+(1|id_cat), data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ CMV | timepoint)))%>%
  mutate(emm2=map(time_model, ~emmeans(., specs = pairwise ~ timepoint | CMV)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  mutate("8 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("24 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("52 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  pivot_longer(cols=ends_with("weeks"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))
  
cmv_msd_sigs <- cmv_msd_purrf%>%
  filter(padj<0.15)%>%
  select(antigen, contrast, padj)

long_msd%>%
  filter(!is.na(CMV), timepoint %in% c("8 weeks", "24 weeks", "52 weeks"))%>%
  filter(antigen %in% cmv_msd_sigs$antigen)%>%
  ggplot(., aes(x=timepoint, y=titer, fill=CMV))+
  geom_boxplot()+
  scale_y_continuous(trans="log10")+
  facet_wrap(~antigen)+
  theme_minimal()


