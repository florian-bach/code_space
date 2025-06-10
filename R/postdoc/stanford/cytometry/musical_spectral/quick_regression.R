library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(emmeans)

`%notin%`=Negate(`%in%`)

musical_cluster_freqs <- read.csv("/Volumes/lab_prasj/BIG_Flo/MUSICAL/spectral/50K_live_singlets/meta35_cluster_frequencies.csv")

meta35 <- musical_cluster_freqs%>%
  group_by(cluster_id, id, timepoint, infectiontype)%>%
  summarise("count"=n())%>%
  group_by(id, timepoint, infectiontype)%>%
  mutate("n_cells"=sum(count))%>%
  mutate(freq=count/n_cells)%>%
  mutate(timepoint_imm=substr(timepoint, 2, nchar(timepoint)))%>%
  mutate("timepoint"=case_when(timepoint_imm==-2 & id %notin% c(176, 363, 577) ~"bad_baseline",
                               timepoint_imm==-2 & id %in% c(176, 363, 577) ~"baseline",
                               timepoint_imm==-1~"baseline",
                               timepoint_imm==0~"day0",
                               timepoint_imm==7~"day7",
                               timepoint_imm==14~"day14",
                               timepoint_imm==28~"day28"))%>%
  mutate(timepoint=factor(timepoint, levels=c("bad_baseline", "baseline", "day0", "day7", "day14", "day28")))%>%
  filter(!grepl("^lrs", id))
  
  
# write.csv(meta35, "~/postdoc/stanford/cytometry/spectral/MUSICAL/new_unmixed/50k_live_singlets/flowsom_meat35_cluster_freqs.csv")
  

meta35_purrf <- meta35%>%
  mutate(id_cat=as.character(id))%>%
  filter(infectiontype %in% c("A", "S"), timepoint!="bad_baseline", timepoint!="day28")%>%
  group_by(cluster_id)%>%
  nest()%>%
  mutate(model=map(data, ~lme4::lmer(freq~timepoint*infectiontype+(1|id_cat), data= .)))%>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(emm=map(model, ~emmeans(., specs = pairwise ~ timepoint | infectiontype)))%>%
  mutate(emm2=map(model, ~emmeans(., specs = pairwise ~ infectiontype | timepoint)))%>%
  
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise")))%>%
  
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  
  mutate("baseline A - baseline S"=map_dbl(emm_contrast_summary2, ~.$p.value[1])) %>%
  mutate("day0 A - day0 S"=map_dbl(emm_contrast_summary2, ~.$p.value[2])) %>%
  
  mutate("baseline A - day0 A"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("baseline A - day14 A"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  mutate("day0 A - day14 A"=map_dbl(emm_contrast_summary, ~.$p.value[5])) %>%
  mutate("baseline S - day0 S"=map_dbl(emm_contrast_summary, ~.$p.value[7])) %>%
  mutate("baseline S - day7 S"=map_dbl(emm_contrast_summary, ~.$p.value[8])) %>%
  mutate("baseline S - day14 S"=map_dbl(emm_contrast_summary, ~.$p.value[9]))%>%
  
  mutate("coef baseline A - baseline S"=map_dbl(emm_contrast_summary2, ~.$estimate[1])) %>%
  mutate("coef day0 A - day0 S"=map_dbl(emm_contrast_summary2, ~.$estimate[2])) %>%
  
  mutate("coef baseline A - day0 A"=map_dbl(emm_contrast_summary, ~.$estimate[1])) %>%
  mutate("coef baseline A - day14 A"=map_dbl(emm_contrast_summary, ~.$estimate[3])) %>%
  mutate("coef day0 A - day14 A"=map_dbl(emm_contrast_summary, ~.$estimate[5])) %>%
  mutate("coef baseline S - day0 S"=map_dbl(emm_contrast_summary, ~.$estimate[7])) %>%
  mutate("coef baseline S - day7 S"=map_dbl(emm_contrast_summary, ~.$estimate[8])) %>%
  mutate("coef baseline S - day14 S"=map_dbl(emm_contrast_summary, ~.$estimate[9]))%>%
  pivot_longer(cols=c(starts_with("baseline"),starts_with("day0")), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))%>%
  ungroup()%>%
  pivot_longer(cols=starts_with("coef"), names_to = "coef_name", values_to = "coef")%>%
  rowwise()%>%
  filter(grepl(contrast, coef_name))

meta35_purrf%>%
  group_by(contrast)%>%
  summarise("n"=sum(padj<0.05))

#nothing
sig_d0_purrf <- meta35_purrf%>%
  filter(p<0.01, contrast=="baseline A - baseline S")%>%
  distinct(cluster_id)

sig_d0_purrf <- meta35_purrf%>%
  filter(padj<0.05, contrast=="baseline S - day0 S")%>%
  distinct(cluster_id)

sig_d7_purrf <- meta35_purrf%>%
  filter(padj<0.05, contrast=="baseline S - day7 S")%>%
  distinct(cluster_id)

  
meta35%>%
  filter(infectiontype %in% c("S"), timepoint!="bad_baseline", timepoint!="day28")%>%
  filter(cluster_id%in%sig_d0_purrf$cluster_id)%>%
  ggplot(., aes(x=timepoint, y=freq, fill = infectiontype))+
  geom_boxplot(outliers = F)+
  facet_wrap(~cluster_id, scales="free")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))
  

# parc ####

parc_data <- read.csv("~/postdoc/stanford/cytometry/spectral/MUSICAL/new_unmixed/50k_live_singlets/parc_cluster_freqs.csv")

parc_freqs <- parc_data %>%
  mutate(id_cat=as.character(id))%>%
  select(-ends_with("count"))%>%
  pivot_longer(cols=ends_with(".percent_total"), names_to = "cluster_id", values_to = "freq")%>%
  mutate(cluster_id=gsub(".percent_total", "", cluster_id))%>%
  mutate(timepoint_imm=substr(timepoint, 2, nchar(timepoint)))%>%
  mutate("timepoint"=case_when(timepoint_imm==-2 & id %notin% c(176, 363, 577) ~"bad_baseline",
                               timepoint_imm==-2 & id %in% c(176, 363, 577) ~"baseline",
                               timepoint_imm==-1~"baseline",
                               timepoint_imm==0~"day0",
                               timepoint_imm==7~"day7",
                               timepoint_imm==14~"day14",
                               timepoint_imm==28~"day28"))%>%
  mutate(timepoint=factor(timepoint, levels=c("bad_baseline", "baseline", "day0", "day7", "day14", "day28")))%>%
  filter(!grepl("^lrs", id))



parc_purrrf <- parc_freqs%>%
  filter(infectiontype %in% c("A", "S"), timepoint!="bad_baseline", timepoint!="day28")%>%
  group_by(cluster_id)%>%
  nest()%>%
  mutate(model=map(data, ~lme4::lmer(freq~timepoint*infectiontype+(1|id_cat), data= .)))%>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(emm=map(model, ~emmeans(., specs = pairwise ~ timepoint | infectiontype)))%>%
  mutate(emm2=map(model, ~emmeans(., specs = pairwise ~ infectiontype | timepoint)))%>%
  
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise", adjust="none")))%>%
  
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  
  mutate("baseline A - baseline S"=map_dbl(emm_contrast_summary2, ~.$p.value[1])) %>%
  mutate("day0 A - day0 S"=map_dbl(emm_contrast_summary2, ~.$p.value[2])) %>%
  
  mutate("baseline A - day0 A"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("baseline A - day14 A"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  mutate("day0 A - day14 A"=map_dbl(emm_contrast_summary, ~.$p.value[5])) %>%
  mutate("baseline S - day0 S"=map_dbl(emm_contrast_summary, ~.$p.value[7])) %>%
  mutate("baseline S - day7 S"=map_dbl(emm_contrast_summary, ~.$p.value[8])) %>%
  mutate("baseline S - day14 S"=map_dbl(emm_contrast_summary, ~.$p.value[9]))%>%
  
  mutate("coef baseline A - baseline S"=map_dbl(emm_contrast_summary2, ~.$estimate[1])) %>%
  mutate("coef day0 A - day0 S"=map_dbl(emm_contrast_summary2, ~.$estimate[2])) %>%
  
  mutate("coef baseline A - day0 A"=map_dbl(emm_contrast_summary, ~.$estimate[1])) %>%
  mutate("coef baseline A - day14 A"=map_dbl(emm_contrast_summary, ~.$estimate[3])) %>%
  mutate("coef day0 A - day14 A"=map_dbl(emm_contrast_summary, ~.$estimate[5])) %>%
  mutate("coef baseline S - day0 S"=map_dbl(emm_contrast_summary, ~.$estimate[7])) %>%
  mutate("coef baseline S - day7 S"=map_dbl(emm_contrast_summary, ~.$estimate[8])) %>%
  mutate("coef baseline S - day14 S"=map_dbl(emm_contrast_summary, ~.$estimate[9]))%>%
  
  pivot_longer(cols=c(starts_with("baseline"),starts_with("day0")), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))%>%
  ungroup()%>%
  pivot_longer(cols=starts_with("coef"), names_to = "coef_name", values_to = "coef")%>%
  rowwise()%>%
  filter(grepl(contrast, coef_name))

parc_purrrf%>%
  group_by(contrast)%>%
  summarise("n"=sum(padj<0.05))

base_diffs <- parc_purrrf%>%
  filter(contrast=="day0 A - day0 S")%>%
  filter(padj<0.1)

day7_diffs <- parc_purrrf%>%
  filter(contrast=="baseline S - day7 S")%>%
  filter(padj<0.05)

parc_freqs%>%
  filter(infectiontype %in% c("S"), timepoint%in%c("baseline", "day0", "day7"))%>%
  filter(cluster_id%in%day7_diffs$cluster_id)%>%
  ggplot(., aes(x=timepoint, y=freq, fill = infectiontype))+
  geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = F)+
  facet_wrap(~cluster_id, scales="free")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))


