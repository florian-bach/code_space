library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(emmeans)

`%notin%` <- Negate(`%in%`)

clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks", "68 weeks")))%>%
  filter(targetName %notin% c("CTSS", "LTA|LTB", "IFNA2"))

mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")


clean_data <- clean_data%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         treatmentarm=case_match(treatmentarm,
                    1~"Placebo",
                    2~"DP 1 year",
                    3~"DP 2 years"))

wide_data <- clean_data%>%
  filter(sample!="11660_tp24_Repeat", treatmentarm!= "DP 2 years", timepoint_num %in% c(8, 52))%>%
  mutate(anymalar=if_else(total_n_malaria_12==0, 0, 1))%>%
  filter(mstatus==0)%>%
  #filter(timepoint_num==52)%>%
  group_by(targetName)%>%
  mutate(z_score=scale(conc, center = T))%>%
  pivot_wider(names_from = c(timepoint), values_from = conc, id_cols = c("id", "treatmentarm","anymalar", "targetName"))%>%
  mutate(fold_change_8_52=`52 weeks`-`8 weeks`)%>%
  pivot_wider(names_from = targetName, values_from = fold_change_8_52, id_cols=c("id", "treatmentarm","anymalar"))

write.csv(wide_data, "~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/wide_data.csv", row.names=FALSE)
