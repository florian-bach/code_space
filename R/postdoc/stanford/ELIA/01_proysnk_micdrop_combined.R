library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(emmeans)

`%notin%` <- Negate(`%in%`)

micdrop_nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks", "68 weeks")))%>%
  filter(targetName %notin% c("CTSS", "LTA|LTB", "IFNA2"))

slim_micdrop_data <- micdrop_nulisa_data%>%
  filter(mstatus==0, treatmentarm!="DP 2 years")%>%
  mutate(id=paste("micdrop_", id, sep=''))%>%
  select(id, timepoint, anyDP, treatmentarm, gender_categorical, targetName, conc, plate)%>%
  mutate(study="MICDROP")
  
prosynk_nulisa_data <- readRDS("~/Library/CloudStorage/Box-Box/Florian Bach's Externally Shareable Files/ELIA/PROSYNK_data_for_Stanford_ELIA/Data/Main/rds/NULISAseq_Data.rds")
prosynk_enrolment_data <- readRDS("~/Library/CloudStorage/Box-Box/Florian Bach's Externally Shareable Files/ELIA/PROSYNK_data_for_Stanford_ELIA/Data/Main/rds/Enrolment_Data.rds")

slim_prosynk_data <- prosynk_nulisa_data%>%
  left_join(., prosynk_enrolment_data, by="IDNumber")%>%
  mutate(id=IDNumber, timepoint=case_match(VisitID,
                                               1~"8 weeks",
                                               2~"12 weeks",
                                               3~"24 weeks",
                                               4~"52 weeks"
                                              ),
         treatmentarm=case_match(Arm,
                                 1~"Labinic",
                                 2~"Lab4b",
                                 3~"Probiotic",
                                 4~"Placebo"),
         gender_categorical=case_match(enr_Sex, 1~"male", 2~"female", 3~"undetermined"),
         conc=targetName_Value,
         anyDP=case_match(Arm,
                          1~"yes",
                          2~"yes",
                          3~"yes",
                          4~"no"), 
         plate=paste("Plate", Plate))%>%
  select(id, timepoint, treatmentarm, gender_categorical, targetName, conc, plate)%>%
  mutate(study="PROSYNK")

#combine and calculate protein & studywise z scores

clean_data <- bind_rows(slim_micdrop_data, slim_prosynk_data)%>%
  filter(targetName %in% slim_micdrop_data$targetName & targetName %in% slim_prosynk_data$targetName)%>%
  filter(timepoint %in% c("8 weeks", "24 weeks", "52 weeks"))%>%
  group_by(targetName, study)%>%
  mutate(normalized_conc=scale(conc, center = T, scale = T))%>%
  ungroup()

write.csv(clean_data, "~/Library/CloudStorage/Box-Box/Florian Bach's Externally Shareable Files/ELIA/raw_prosynk_micdrop_combined.csv", row.names = F)

