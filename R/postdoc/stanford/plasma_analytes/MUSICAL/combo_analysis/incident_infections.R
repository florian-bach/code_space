library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(emmeans)

`%notin%` <- Negate(`%in%`)

nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/tr1_paper/revised_baseline_clean_musical_combo_with_metadata.csv")
nulisa_data <- nulisa_data%>%
  mutate(infectiontype=substr(infectiontype, 1, 1))

# need to convert this to make infection itself incident or not; 
# this data only contains info on whether the date counts as an incident infection or not (so day0)
all_timepoints <-nulisa_data%>%
  filter(infectiontype%in% c("S", "A"))%>%
  distinct(id, date, infectiontype, timepoint)

day0s <- nulisa_data%>%
  filter(timepoint=="day0"&infectiontype=="S" | timepoint=="day0"&infectiontype=="A")%>%
  distinct(id, infectiontype, date)

incident_infections <- read.csv("~/postdoc/stanford/clinical_data/MUSICAL/PRISM Border Cohort Study all visits database_FINAL_probs.csv")
incident_infections <- incident_infections%>%
  mutate(id= cohortid)





#this is good
additional_timepoints <- incident_infections%>%
  filter(
  id==132& date=="2021-06-21"|
  id==261& date=="2022-01-11"|
  id==317& date=="2022-01-12"|
  id==410& date=="2021-02-16"|
  id==495& date=="2022-01-21")%>%
  mutate(infectiontype=c("S", "S", "S", "A","A"))
  
day0s_expanded <- rbind()
# 1   132 S                 
# 2   261 S               
# 3   268 S                
# 4   316 A                 
# 5   317 S                 
# 6   336 A                 
# 7   410 A                 
# 8   495 S                 


day0s_with_infs <- day0s%>%
  inner_join(., incident_infections, by=c('id', 'date'))%>%
  bind_rows(., additional_timepoints)

nulisa_data_with_infs <- nulisa_data %>%
  left_join(., day0s_with_infs, by=c('id', 'infectiontype'))%>%
  group_by(id, infectiontype)%>%
  mutate(incident_infection=ifelse(any(incident==1), "incident", "non_incident"))

combo_as_purff <- nulisa_data_with_infs %>%
  filter(infectiontype%in%c("A", "S"), timepoint %notin% c("day28"), timepoint!="bad_baseline")%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  # mutate(coif=paste(COI, "parasites"))%>%
  group_by(targetName)%>%
  nest() %>%
  #mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id_cat), data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+COI+ageyrs+gender_categorical+(1|id_cat), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  # mutate(emm=map(model, ~emmeans(., specs = ~ timepoint:infectiontype)))%>%
  mutate(emm=map(model, ~emmeans(., specs = pairwise ~ COI | timepoint : infectiontype)))%>%

  mutate(emm_contrast=map(emm, ~contrast(., "pairwise")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise")))%>%
  
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  
  mutate("baseline A incident - A non_incident"=map_dbl(emm_contrast_summary2, ~.$p.value[2])) %>%
  mutate("day0 A incident - A non_incident"=map_dbl(emm_contrast_summary2, ~.$p.value[8])) %>%
  mutate("day14 A incident - A non_incident"=map_dbl(emm_contrast_summary2, ~.$p.value[14])) %>%
  mutate("baseline S incident - S non_incident"=map_dbl(emm_contrast_summary2, ~.$p.value[5])) %>%
  mutate("day0 S incident - S non_incident"=map_dbl(emm_contrast_summary2, ~.$p.value[11])) %>%
  mutate("day14 S incident - S non_incident"=map_dbl(emm_contrast_summary2, ~.$p.value[17])) %>%
  pivot_longer(cols=ends_with("non_incident"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))


nulisa_data_with_infs%>%
  filter(targetName%in%c("IL17C", "MMP12", "NGF", "S100A9","TNFRSF4"), infectiontype %in% c("A", "S"), timepoint!="bad_baseline")%>%
  ggplot(., aes(x=timepoint, y=concentration, fill=incident_infection))+
  geom_boxplot()+
  facet_wrap(~targetName+infectiontype)+
  theme_minimal()



day0s_with_infs%>%
  ggplot(., aes(x=COI))+
  geom_density(aes(color=infectiontype))+
  theme_minimal()

coi_plot <- day0s_with_infs%>%
  add_row(COI=c(10:13, 16:19), infectiontype="zer")%>%
  ggplot(., aes(x=factor(COI)))+
  scale_alpha_manual(values=c(1,1,0))+
  geom_histogram(aes(fill=infectiontype, alpha = infectiontype), stat="count", position = position_dodge())+
  scale_fill_manual(values=c("grey", "darkred", "yellow"), limits=c("A", "S"))+
  scale_color_manual(values=c("grey", "darkred", "yellow"), guide="none")+
  # scale_x_continuous(breaks=seq(1,20,1))+
  xlab("Number of AMA1 Genotypes")+
  ylab("")+
  scale_y_continuous(breaks=seq(0,20,2))+
  guides(alpha=guide_none())+
  theme_minimal(base_size = 16)+
  theme(legend.title = element_blank())

ggsave("~/postdoc/stanford/abstracts/keystone25/poster_figures/coi_plot.png", coi_plot, height=6, width=8, dpi=444, bg="white")
