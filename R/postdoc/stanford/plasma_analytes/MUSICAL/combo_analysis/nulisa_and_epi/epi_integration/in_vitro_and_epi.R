library(tidyr)
library(dplyr)
library(ggplot2)

`%notin%` <- Negate(`%in%`)

epi_data <- read.csv("~/postdoc/stanford/clinical_data/MUSICAL/arefin_epi/primary_samples_2025Jan9.csv")

slim_epi_data <- epi_data%>%
  mutate(id=combined_id)%>%
  filter(baseline_date_bool=="True")%>%
  mutate(infectiontype=ifelse(grepl("S", inf_type,), "S", "A"))%>%
  select(id, infectiontype, baseline_date_bool, cleaner_baseline_date_bool, parasitemia_12months_prior_micro, parasitemia_12months_after_micro,
         parasitemia_12months_prior_qpcr, parasitemia_12months_after_qpcr, count_of_symp_mal_12months_before, count_of_symp_mal_12months_after)%>%
  mutate(symp_ratio_before=count_of_symp_mal_12months_before/parasitemia_12months_prior_micro,
         symp_ratio_after=count_of_symp_mal_12months_after/parasitemia_12months_after_micro)


cell_count_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/df_jason_analysis.csv")

slim_cell_count_data <- cell_count_data%>%
  pivot_longer(cols = ends_with("Frequency"), names_to = "gate", values_to = "freq")%>%
  mutate(timepoint_imm=timepoint)%>%
  mutate("timepoint"=case_when(timepoint_imm==-2 & id %notin% c(176, 363, 577) ~"bad_baseline",
                               timepoint_imm==-2 & id %in% c(176, 363, 577) ~"baseline",
                               timepoint_imm==-1~"baseline",
                               timepoint_imm==0~"day0",
                               timepoint_imm==7~"day7",
                               timepoint_imm==14~"day14",
                               timepoint_imm==28~"day28"))%>%
  select(id, timepoint, timepoint_imm, infectiontype, gate, stim, freq)

combo_data <- left_join(slim_cell_count_data, slim_epi_data, by=c("id", "infectiontype"))


combo_data %>%
  filter(gate%in%c("Tr1_Frequency", "Tr1_IL10_Frequency", "Tr1_IFNg_Frequency"), timepoint=="baseline", stim=="iRBCs")%>%
  ggplot(., aes(x=count_of_symp_mal_12months_after, y=freq, fill=factor(count_of_symp_mal_12months_after)))+
  geom_boxplot(outliers = F)+
  viridis::scale_fill_viridis(discrete = T)+
  # ggtitle("more Tr1 proteins in plasma associated with lower symptomatic ratio of subsequent infection")+
  theme_minimal()+
  facet_wrap(~gate, scales="free")+
  theme(legend.position = "none")

combo_data %>%
  filter(gate=="Tr1_Frequency", timepoint=="baseline", stim=="unstim")%>%
  ggplot(., aes(x=count_of_symp_mal_12months_after, y=freq, fill=factor(count_of_symp_mal_12months_after)))+
  # geom_boxplot(outliers = F)+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  viridis::scale_fill_viridis(discrete = T)+
  # ggtitle("more Tr1 proteins in plasma associated with lower symptomatic ratio of subsequent infection")+
  theme_minimal()+
  theme(legend.position = "none")


combo_data %>%
  filter(gate%in%c("Tr1_Frequency", "Tr1_IL10_Frequency", "Tr1_IFNg_Frequency"), timepoint=="baseline", stim=="iRBCs")%>%
  ggplot(., aes(x=count_of_symp_mal_12months_before, y=freq, fill=factor(count_of_symp_mal_12months_before)))+
  geom_boxplot(outliers = F)+
  viridis::scale_fill_viridis(discrete = T)+
  # ggtitle("more Tr1 proteins in plasma associated with lower symptomatic ratio of subsequent infection")+
  theme_minimal()+
  theme(legend.position = "none")


#ratio
combo_data %>%
  filter(gate=="Tr1_Frequency", timepoint=="baseline", stim=="unstim")%>%
  filter(timepoint=="baseline")%>%
  ggplot(., aes(x=symp_ratio_after, y=freq))+
  geom_point()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method = "spearman")+
  # ggtitle("more Tr1 proteins in plasma associated with lower symptomatic ratio of subsequent infection")+
  theme_minimal()+
  theme(legend.position = "none")
