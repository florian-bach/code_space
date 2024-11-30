library(emmeans)
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)



`%notin%` <- Negate(`%in%`)

da_boxplot_theme <- theme(legend.position = "none",
                          axis.title = element_blank())

# read clean data ####
clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")

clean_data <- clean_data %>%
  mutate(timepoint = factor(timepoint, levels=c("bad_baseline", "baseline", "day0", "day7", "day14", "day28")))

# emmeans ####
combo_as_purff <- clean_data %>%
  filter(infectiontype%in%c("A", "S"), timepoint %notin% c("day7", "day28"), timepoint!="bad_baseline")%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  group_by(targetName)%>%
  nest() %>%
  #mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  # mutate(emm=map(model, ~emmeans(., specs = ~ timepoint:infectiontype)))%>%
  mutate(emm=map(model, ~emmeans(., specs = pairwise ~ timepoint | infectiontype)))%>%
  mutate(emm2=map(model, ~emmeans(., specs = pairwise ~ infectiontype | timepoint)))%>%
  
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise")))%>%
  
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  
  mutate("baseline A - baseline S"=map(emm_contrast_summary2, ~.$p.value[1])) %>%
  mutate("baseline A - day0 A"=map(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("baseline A - day14 A"=map(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("baseline S - day0 S"=map(emm_contrast_summary, ~.$p.value[4])) %>%
  mutate("baseline S - day14 S"=map(emm_contrast_summary, ~.$p.value[5]))%>%
  pivot_longer(cols=starts_with("baseline"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))
  

# incorporate fold change cuttoff?
combo_as_purff%>%
  group_by(contrast)%>%
  summarise("n"=sum(padj<0.05))

sig_base_day0_s <-combo_as_purff%>%
  filter(padj<0.05 & contrast=="baseline S - day0 S")%>%
  arrange(padj)

# plot differences ####
## baseline day 0  S ####

for(i in 1:ceiling(nrow(sig_base_day0_s)/16)){
  
  plt <- clean_data %>%
    filter(infectiontype %in% c("S"), timepoint!="day28", timepoint!="bad_baseline")%>%
    filter(targetName %in% sig_base_day0_s$targetName[seq((i-1)*16+1, i*16)])%>%
    mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
    ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
    geom_line(aes(group=id), alpha=0.2)+
    geom_boxplot(outliers = FALSE)+
    facet_wrap(~targetName, scales = "free")+
    scale_fill_manual(values=viridis::magma(5))+
    theme_minimal()+
    da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/model_plots/baseline_S_day0_S", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}


## baseline day 0  A ####

sig_base_day0_a <-combo_as_purff%>%
  filter(padj<0.05 & contrast=="baseline A - day0 A")%>%
  arrange(padj)

for(i in 1:ceiling(nrow(sig_base_day0_a)/16)){
  
  plt <- clean_data %>%
    filter(infectiontype %in% c("A"), timepoint!="day28", timepoint!="bad_baseline")%>%
    filter(targetName %in% sig_base_day0_s$targetName[seq((i-1)*16+1, i*16)])%>%
    mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
    ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
    geom_line(aes(group=id), alpha=0.2)+
    geom_boxplot(outliers = FALSE)+
    facet_wrap(~targetName, scales = "free")+
    scale_fill_manual(values=viridis::magma(5))+
    theme_minimal()+
    da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/model_plots/baseline_A_day0_A", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}


## baseline baseline A  S ####

sig_base_base <-combo_as_purff%>%
  filter(p<0.05 & contrast=="baseline A - baseline S")%>%
  arrange(p)



for(i in 1:ceiling(nrow(sig_base_base)/16)){
  
  plt <- clean_data %>%
    filter(infectiontype %in% c("A", "S"), timepoint=="baseline")%>%
    filter(targetName %in% sig_base_base$targetName[seq((i-1)*16+1, i*16)])%>%
    mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
    ggplot(aes(x=infectiontype, y=concentration, fill=infectiontype))+
    geom_line(aes(group=id), alpha=0.2)+
    geom_boxplot(outliers = FALSE)+
    facet_wrap(~targetName, scales = "free")+
    scale_fill_manual(values=viridis::magma(3))+
    theme_minimal()+
    da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/model_plots/baseline_AS", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}




## baseline day 14  A ####

sig_base_day14_a <-combo_as_purff%>%
  filter(padj<0.05 & contrast=="baseline A - day14 A")%>%
  arrange(padj)

for(i in 1:ceiling(nrow(sig_base_day14_a)/16)){
  
  plt <- clean_data %>%
    filter(infectiontype %in% c("A"), timepoint!="day28", timepoint!="bad_baseline")%>%
    mutate(day14_para=if_else(timepoint=="day14" & parasitedensity > 10, "parasitemic_day14", "no_parasites_day14"))%>%
    group_by(id)%>%
    mutate(class2= if_else(any(day14_para=="parasitemic_day14"), "non_controller", "controller"))%>%
    filter(targetName %in% sig_base_day14_a$targetName[seq((i-1)*16+1, i*16)])%>%
    mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
    ggplot(aes(x=factor(timepoint), y=concentration, fill=class2))+
    geom_line(aes(group=id), alpha=0.2)+
    geom_boxplot(outliers = FALSE)+
    facet_wrap(~targetName, scales = "free")+
    scale_fill_manual(values=viridis::magma(3))+
    theme_minimal()+
    theme(legend.title = element_blank(),
          legend.position = "right")
    da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/model_plots/baseline_A_day14_A", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}


sig_base_day14_a$targetName [sig_base_day14_a$targetName %notin% sig_base_day0_a$targetName]
sig_base_day0_a$targetName [sig_base_day0_a$targetName %notin% sig_base_day14_a$targetName ]




# are parasite controlllers different ####

day14_a_purff <- clean_data %>%
  filter(infectiontype%in%c("A"), timepoint %notin% c("day7", "day28"), timepoint!="bad_baseline")%>%
  mutate(day14_para=if_else(timepoint=="day14" & parasitedensity > 10, "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), "non_controller", "controller"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  group_by(targetName)%>%
  nest() %>%
  #mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~timepoint*class2+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  # mutate(emm=map(model, ~emmeans(., specs = ~ timepoint:infectiontype)))%>%
  # mutate(emm=map(model, ~emmeans(., specs = pairwise ~ timepoint | infectiontype)))%>%
  mutate(emm2=map(model, ~emmeans(., specs = pairwise ~ class2 | timepoint)))%>%
  
  # mutate(emm_contrast=map(emm, ~contrast(., "pairwise")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise", adjust="none")))%>%
  
  # mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  
  mutate(p=map_dbl(emm_contrast_summary2, ~.$p.value[1])) %>%
  # mutate("baseline A - day0 A"=map(emm_contrast_summary, ~.$p.value[1])) %>%
  # mutate("baseline A - day14 A"=map(emm_contrast_summary, ~.$p.value[2])) %>%
  # mutate("baseline S - day0 S"=map(emm_contrast_summary, ~.$p.value[4])) %>%
  # mutate("baseline S - day14 S"=map(emm_contrast_summary, ~.$p.value[5]))%>%
  # pivot_longer(cols=starts_with("baseline"), names_to = "contrast", values_to = "p")%>%
  ungroup()%>%
  mutate(padj = p.adjust(p, method="fdr"))

sig_ctrl <- day14_a_purff %>%
  filter(p<0.05)%>%
  arrange(p)

clean_data %>%
  filter(infectiontype %in% c("A"), timepoint!="day28", timepoint!="bad_baseline")%>%
  mutate(day14_para=if_else(timepoint=="day14" & parasitedensity > 10, "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), "non_controller", "controller"))%>%
  filter(targetName %in% sig_ctrl$targetName)%>%#[seq((i-1)*16+1, i*16)])%>%
  filter(timepoint=="baseline")%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=class2))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetName, scales = "free")+
  scale_fill_manual(values=viridis::magma(3))+
  theme_minimal()+
  theme(legend.title = element_blank(),
        legend.position = "right")
da_boxplot_theme