library(emmeans)
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)


`%notin%`=Negate(`%in%`)

da_boxplot_theme <- theme(legend.position = "none",
                          axis.title = element_blank())

# read clean data ####
clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")

clean_data <- clean_data %>%
  mutate(timepoint = factor(timepoint, levels=c("bad_baseline", "baseline", "day0", "day7", "day14", "day28")))

# emmeans ####
combo_as_purff <- clean_data %>%
  filter(infectiontype%in%c("A", "S"), timepoint %notin% c("day28"), timepoint!="bad_baseline")%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  group_by(targetName)%>%
  nest() %>%
  #mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id_cat), data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id_cat), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  # mutate(emm=map(model, ~emmeans(., specs = ~ timepoint:infectiontype)))%>%
  mutate(emm=map(model, ~emmeans(., specs = pairwise ~ timepoint | infectiontype)))%>%
  mutate(emm2=map(model, ~emmeans(., specs = pairwise ~ infectiontype | timepoint)))%>%
  
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
  

# incorporate fold change cuttoff?
# reducing data to only include infections with baselines less than a year apart
# does not meaningfully change results, no differences between baselines
combo_as_purff%>%
  group_by(contrast)%>%
  summarise("n"=sum(padj<0.05))

sig_base_day0_s <-combo_as_purff%>%
  filter(padj<0.05 & contrast=="baseline S - day0 S")%>%
  select(targetName, contrast, "coef baseline S - day0 S", padj)%>%
  arrange(padj)

combo_as_results <- combo_as_purff %>%
  select(targetName, contrast, coef_name, coef, padj)%>%
  arrange(padj)
 
write.csv(combo_as_results, "~/postdoc/stanford/plasma_analytes/MUSICAL/combo/differential_abundance/combo_as_purff.csv", row.names = FALSE)
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
    filter(infectiontype %in% c("A", "S"), timepoint!="day28", timepoint!="bad_baseline")%>%
    filter(targetName %in% sig_base_day0_a$targetName[seq((i-1)*16+1, i*16)])%>%
    mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
    ggplot(aes(x=factor(timepoint), y=concentration, fill=infectiontype))+
    geom_line(aes(group=id), alpha=0.2)+
    geom_boxplot(outliers = FALSE)+
    facet_wrap(~targetName, scales = "free")+
    # scale_fill_manual(values=viridis::magma(5))+
    scale_fill_manual(values = c("A"="grey", "S"="darkred"))+
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


sig_base_as <- combo_as_purff%>%
  filter(contrast=="baseline A - baseline S", padj <0.1)%>%
  distinct(targetName)

clean_data%>%
  mutate(inf_type_dich=ifelse(grepl("^S",infectiontype), "S", ifelse(grepl("^A", infectiontype), "A", infectiontype)))%>%
  filter(timepoint=="baseline", inf_type_dich%in%c("A", "S"), targetName %in% sig_base_as$targetName)%>%
  ggplot(., aes(x=inf_type_dich, y=concentration, fill=inf_type_dich))+
  geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName, scales="free")+
  scale_fill_manual(values=c("darkgrey", "darkred"))+
  theme_minimal()+
  ggtitle("baseline")+
  xlab("Asymptomatic / Symptomatic")+
  theme(legend.position = "none")


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
    ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
    geom_line(aes(group=id), alpha=0.2)+
    geom_boxplot(outliers = FALSE)+
    facet_wrap(~targetName, scales = "free")+
    scale_fill_manual(values=viridis::magma(5))+
    theme_minimal()+
    theme(legend.title = element_blank(),
          legend.position = "right")
    da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/model_plots/baseline_A_day14_A", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}


sig_base_day14_a$targetName [sig_base_day14_a$targetName %notin% sig_base_day0_a$targetName]
sig_base_day0_a$targetName [sig_base_day0_a$targetName %notin% sig_base_day14_a$targetName ]




# are parasite controlllers different ####

day14_a_purff <- nulisa_data %>%
  filter(infectiontype%in%c("A"), timepoint %notin% c("day7", "day28", "bad_baseline"))%>%
  mutate(day14_para=if_else(timepoint=="day14" & parasitedensity > 1, "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), "non_controller", "controller"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  group_by(targetName)%>%
  nest() %>%
  #mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id_cat), data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~timepoint*class2+ageyrs+gender_categorical+(1|id_cat), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(emm=map(model, ~emmeans(., specs = ~ timepoint:class2)))%>%
  # mutate(emm=map(model, ~emmeans(., specs = pairwise ~ timepoint | infectiontype)))%>%
  mutate(emm2=map(model, ~emmeans(., specs = pairwise ~ class2 | timepoint)))%>%
  
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise", adjust="none")))%>%
  
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  
  mutate("controller - non_controller baseline p"=map_dbl(emm_contrast_summary2, ~.$p.value[1]),
         "controller - non_controller day0 p"=map_dbl(emm_contrast_summary2, ~.$p.value[2]),
         "controller - non_controller day14 p"=map_dbl(emm_contrast_summary2, ~.$p.value[3]),
         "baseline controller - day14 controller p"=map_dbl(emm_contrast_summary, ~.$p.value[2]),
         "baseline non_controller - day14 non_controller p"=map_dbl(emm_contrast_summary, ~.$p.value[14]))%>%
  # mutate("baseline A - day0 A"=map(emm_contrast_summary, ~.$p.value[1])) %>%
  # mutate("baseline A - day14 A"=map(emm_contrast_summary, ~.$p.value[2])) %>%
  # mutate("baseline S - day0 S"=map(emm_contrast_summary, ~.$p.value[4])) %>%
  # mutate("baseline S - day14 S"=map(emm_contrast_summary, ~.$p.value[5]))%>%
  pivot_longer(cols=ends_with("p"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

#nothing
sig_base_control <- day14_a_purff %>%
  filter(contrast=="controller - non_controller baseline p", padj<0.05)%>%
  select(targetName, p, padj)%>%
  arrange(padj)

#nothing
sig_day0_control <- day14_a_purff %>%
  filter(contrast=="controller - non_controller day0 p", padj<0.05)%>%
  select(targetName, p, padj)%>%
  arrange(padj)

#30!
sig_day14_control <- day14_a_purff %>%
  filter(contrast=="controller - non_controller day14 p", padj<0.05)%>%
  select(targetName, p, padj)%>%
  arrange(padj)


sig_ctrl_base_14 <- day14_a_purff%>%
  filter(contrast=="baseline controller - day14 controller p", padj<0.05)%>%
  select(targetName, p, padj)%>%
  arrange(padj)

sig_non_ctrl_base_14 <- day14_a_purff%>%
  filter(contrast=="baseline non_controller - day14 non_controller p", padj<0.05)%>%
  select(targetName, p, padj)%>%
  arrange(padj)
#Anything greater than 0.001 based on the qPCR should be considered a true positive; although values below 0.1
#should not be considered accurate in terms of quantification we still believe they represent true positives.
#In other words, 0.1 is the limit of quantification and LOD is likely around 0.02 or so,
#but values below 0.1 and above 0.001 are qualitatively positive

ctrL_plot <- clean_data %>%
  filter(infectiontype %in% c("A"), timepoint!="day28", timepoint!="bad_baseline")%>%
  mutate(day14_para=if_else(timepoint=="day14" & qpcr > 1, "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), "non_controller", "controller"))%>%
  # filter(targetName %in% sig_ctrl$targetName)%>%
  filter(targetName %in% sig_day14_control$targetName)%>%#[seq((i-1)*16+1, i*16)])%>%
  filter(class2 %in% c("non_controller", "controller"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=class2))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetName, scales = "free", nrow=5)+
  scale_fill_manual(values=viridis::magma(n=3))+
  theme_minimal()+
  theme(legend.title = element_blank(),
        legend.position = "right",
        axis.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/model_plots/ctrl_plot.png", ctrL_plot, height=5.33, width=8, dpi=444, bg="white")


ctrL_plot2 <- clean_data %>%
  filter(infectiontype %in% c("A"), timepoint!="day28", timepoint!="bad_baseline")%>%
  mutate(day14_para=if_else(timepoint=="day14" & qpcr > 1, "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), "non_controller", "controller"))%>%
  # filter(targetName %in% sig_ctrl$targetName)%>%
  filter(targetName %in% c("IL15", sig_non_ctrl_base_14$targetName))%>%#[seq((i-1)*16+1, i*16)])%>%
  filter(class2 %in% c("non_controller", "controller"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=class2))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetName, scales = "free", nrow=5)+
  scale_fill_manual(values=viridis::magma(n=3))+
  theme_minimal()+
  theme(legend.title = element_blank(),
        legend.position = "right",
        axis.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/model_plots/ctrl_plot2.png", ctrL_plot2, height=5.33, width=8, dpi=444, bg="white")


# parasitemia only model ####

combo_as_purff2 <- clean_data %>%
  filter(infectiontype%in%c("A", "S"), timepoint %notin% c("day7", "day28"), timepoint!="bad_baseline")%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  group_by(targetName)%>%
  nest() %>%
  #mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id_cat), data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~infectiontype*log_qpcr+ageyrs+gender_categorical+(1|id_cat), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  # mutate(emm=map(model, ~emmeans(., specs = ~ timepoint:infectiontype)))%>%
  # mutate(emm=map(model, ~emmeans(., specs = pairwise ~ timepoint:log_qpcr)))%>%
  mutate(emm2=map(model, ~emmeans(., specs = pairwise ~ log_qpcr | infectiontype)))%>%
  
  # mutate(emm_contrast=map(emm, ~contrast(., "pairwise")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise")))%>%
  
  # mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
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


# overview_plots ####


clean_data %>%
  filter(infectiontype %in% c("A", "S"), timepoint%notin%c("day28","bad_baseline", "day7"))%>%
  mutate(day14_para=if_else(timepoint=="day14" & qpcr > 1, "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), "non_controller", "controller"))%>%
  # filter(targetName %in% sig_ctrl$targetName)%>%
  filter(targetName %in% c("VSNL1", "GFAP", "GDF2", "SPP1"))%>%#[seq((i-1)*16+1, i*16)])%>%
  filter(class2 %in% c("non_controller", "controller"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=factor(timepoint)))+
  geom_line(aes(group=id), alpha=0.2)+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetName, scales = "free", nrow=2)+
  scale_fill_manual(values=time_cols)+
  theme_minimal()+
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.title = element_blank())
da_boxplot_theme