library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(emmeans)

# preamble
epi_data <- read.csv("~/postdoc/stanford/clinical_data/MUSICAL/arefin_epi/primary_samples_2025Jan9.csv")

# "baseline_date_bool", "cleaner_baseline_date_bool", "parasitemia_12months_prior_micro",
# "parasitemia_12months_after_micro", "parasitemia_12months_prior_qpcr", "parasitemia_12months_after_qpcr",
# "count_of_symp_mal_12months_before","count_of_symp_mal_12months_after"

slim_epi_data <- epi_data%>%
  mutate(id=combined_id)%>%
  filter(baseline_date_bool=="True")%>%
  mutate(infectiontype=ifelse(grepl("S", inf_type,), "S", "A"))%>%
  select(id, infectiontype, baseline_date_bool, cleaner_baseline_date_bool, parasitemia_12months_prior_micro, parasitemia_12months_after_micro,
         parasitemia_12months_prior_qpcr, parasitemia_12months_after_qpcr, count_of_symp_mal_12months_before, count_of_symp_mal_12months_after)%>%
  mutate(symp_ratio_before=count_of_symp_mal_12months_before/parasitemia_12months_prior_micro,
         symp_ratio_after=count_of_symp_mal_12months_after/parasitemia_12months_after_micro)


nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")

slim_nulisa_data <- nulisa_data %>%
  filter(infectiontype %in% c("A", "S"))%>%
  select(id, date, id_cat, ageyrs, gender_categorical, timepoint, timepoint_imm, infectiontype, targetName, concentration, new_qpcr, parasitedensity)

combo_data <- left_join(slim_nulisa_data, slim_epi_data, by=c("id", "infectiontype"))

# overview plots ####
combo_data %>%
  filter(targetName %in% c("IL10", "GZMA", "LAG3"))%>%
  filter(timepoint=="baseline")%>%
  ggplot(., aes(x=count_of_symp_mal_12months_before, y=concentration, fill=factor(count_of_symp_mal_12months_before)))+
  geom_boxplot()+
  ggtitle("higher baseline levels of IL10 and LAG3 are associated with\nmore frequent symptomatic malaria in the 12 months prior")+
  viridis::scale_fill_viridis(discrete = T)+
  theme_minimal()+
  facet_wrap(~targetName, scales = "free")+
  theme(legend.position = "none")

combo_data %>%
  filter(targetName %in% c("CRP", "IL6", "IL27", "TNF"))%>%
  filter(timepoint=="day0", infectiontype=="S")%>%
  ggplot(., aes(x=count_of_symp_mal_12months_before, y=concentration, fill=factor(count_of_symp_mal_12months_before)))+
  geom_boxplot()+
  ggtitle("peak inflammation not clearly associated with previous malaria incidence")+
  viridis::scale_fill_viridis(discrete = T)+
  theme_minimal()+
  facet_wrap(~targetName, scales="free")+
  theme(legend.position = "none")


combo_data %>%
  filter(targetName %in% c("IL10", "LILRB2", "LAG3", "CRP"))%>%
  filter(timepoint=="day0", infectiontype=="A")%>%
  ggplot(., aes(x=count_of_symp_mal_12months_before, y=concentration, fill=factor(count_of_symp_mal_12months_before)))+
  geom_boxplot()+
  ggtitle("more malaria associated with higher response at asymptomatic day0")+
  viridis::scale_fill_viridis(discrete = T)+
  theme_minimal()+
  facet_wrap(~targetName, scales="free")+
  theme(legend.position = "none")



combo_data %>%
  filter(targetName %in% c("CRP", "IL27", "TNF", "IL6"))%>%
  filter(timepoint=="day0", infectiontype=="S")%>%
  ggplot(., aes(x=count_of_symp_mal_12months_after, y=concentration, fill=factor(count_of_symp_mal_12months_after)))+
  geom_boxplot()+
  viridis::scale_fill_viridis(discrete = T)+
  ggtitle("higher inflammation at symptomatic day 0 potentially associated with subsequently lower incidence of symptomatic malaria")+
  theme_minimal()+
  facet_wrap(~targetName, scales="free")+
  theme(legend.position = "none")
  

combo_data %>%
  filter(targetName %in% c("CRP", "IL27", "TNF", "IL6", "IL10", "IL1RN"))%>%
  filter(timepoint=="day0", infectiontype=="S")%>%
  ggplot(., aes(x=parasitemia_12months_after_micro, y=concentration, fill=factor(parasitemia_12months_after_micro)))+
  geom_boxplot()+
  viridis::scale_fill_viridis(discrete = T)+
  ggtitle("higher inflammation at symptomatic day 0 potentially associated with fewer subsequent parasitemic months")+
  theme_minimal()+
  facet_wrap(~targetName, scales="free")+
  theme(legend.position = "none")


combo_data %>%
  filter(targetName %in% c("IL10", "LILRB2", "LAG3", "CRP"))%>%
  filter(timepoint=="baseline")%>%
  ggplot(., aes(x=symp_ratio_after, y=concentration))+
  geom_point()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method = "spearman")+
  viridis::scale_fill_viridis(discrete = T)+
  ggtitle("more Tr1 proteins in plasma associated with lower symptomatic ratio of subsequent infection")+
  theme_minimal()+
  facet_wrap(~targetName, scales="free")+
  theme(legend.position = "none")




# both baselines ####
baseline_as_purff <- combo_data %>%
  filter(infectiontype%in%c("A", "S"), timepoint=="baseline")%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  group_by(targetName)%>%
  nest() %>%
  #mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id_cat), data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~count_of_symp_mal_12months_after+ageyrs+gender_categorical+(1|id_cat), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  # mutate(emm=map(model, ~emmeans(., specs = ~ timepoint:infectiontype)))%>%
  mutate(emm=map(model, ~emmeans(., ~count_of_symp_mal_12months_after)))%>%

  mutate(emm_contrast=map(emm, ~contrast(., "pairwise")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise")))%>%
  
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  
  mutate("baseline A - baseline S"=map_dbl(emm_contrast_summary2, ~.$p.value[1])) %>%
  mutate("baseline A - day0 A"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("baseline A - day14 A"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  mutate("baseline S - day0 S"=map_dbl(emm_contrast_summary, ~.$p.value[7])) %>%
  mutate("baseline S - day7 S"=map_dbl(emm_contrast_summary, ~.$p.value[8])) %>%
  mutate("baseline S - day14 S"=map_dbl(emm_contrast_summary, ~.$p.value[9]))%>%
  mutate("coef baseline A - baseline S"=map_dbl(emm_contrast_summary2, ~.$estimate[1])) %>%
  mutate("coef baseline A - day0 A"=map_dbl(emm_contrast_summary, ~.$estimate[1])) %>%
  mutate("coef baseline A - day14 A"=map_dbl(emm_contrast_summary, ~.$estimate[3])) %>%
  mutate("coef baseline S - day0 S"=map_dbl(emm_contrast_summary, ~.$estimate[7])) %>%
  mutate("coef baseline S - day7 S"=map_dbl(emm_contrast_summary, ~.$estimate[8])) %>%
  mutate("coef baseline S - day14 S"=map_dbl(emm_contrast_summary, ~.$estimate[9]))%>%
  pivot_longer(cols=starts_with("baseline"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))%>%
  ungroup()%>%
  pivot_longer(cols=starts_with("coef"), names_to = "coef_name", values_to = "coef")%>%
  rowwise()%>%
  filter(grepl(contrast, coef_name))



# epi data for Reuben ####

data_for_reuben <- combo_data%>%
  select(-targetName, -concentration, -baseline_date_bool, -cleaner_baseline_date_bool)%>%
  distinct()

write.csv(data_for_reuben, "~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/Cohort data/cohort_data_for_reuben.csv")
