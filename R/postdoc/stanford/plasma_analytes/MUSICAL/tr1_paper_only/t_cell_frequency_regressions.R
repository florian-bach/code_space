library(emmeans)
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)

stim_palette <- c("darkred", "darkblue", "black")
names(stim_palette) <- c("iRBCs", "PMA", "unstim")

time_cols <- list("baseline"="#E4DEBD",
                  "day0" = "#C03F3E",
                  "day7" = "#D87E1F",
                  "day14" = "#E6B85F")



long_combo <- read.csv("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/cell_freqs_and_nulisa.csv")
long_combo <-long_combo%>%
  mutate(timepoint=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))
  



individuals_with_parasites_at_day14 <- long_combo%>%
  distinct(id, infectiontype, timepoint, log_qpcr)%>%
  filter(infectiontype %in% c("A"), timepoint!="day28", timepoint!="bad_baseline")%>%
  mutate(day14_para=if_else(timepoint=="day14" & log_qpcr > 0, "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), "non_controller", "controller"))%>%
  filter(class2=="non_controller")%>%
  distinct(id, infectiontype, class2)
  
  
long_combo%>%
  filter(gate %in% c("IL10_Frequency", "IFNg_Frequency", "IL21_Frequency"), stim!="unstim")%>%
  distinct(id, infectiontype, timepoint, gate, stim, freq)%>%
  filter(infectiontype%in%c("A", "S"), !is.na(timepoint))%>%
  mutate(class2=case_when(infectiontype=="A"&id%in%individuals_with_parasites_at_day14$id~"parasites at day 14",
                          infectiontype=="A"&id%notin%individuals_with_parasites_at_day14$id~"no parasites at day 14", .default = NULL))%>%
  ggplot(., aes(x=class2, y=freq/100, fill=timepoint))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = F)+
  scale_y_continuous(labels=scales::label_percent())+
  scale_fill_manual(values=time_cols)+
  scale_color_manual(values=c("black", "darkred", "black"))+
  facet_wrap(~gate+infectiontype+stim, scales="free", nrow=2)+
  theme_minimal()


long_combo%>%
  filter(gate %in% c("Tr1_IL10_Frequency", "Tr1_IFNg_Frequency", "Tr1_IL21_Frequency"), stim!="unstim")%>%
  distinct(id, infectiontype, timepoint, gate, stim, freq)%>%
  filter(infectiontype%in%c("A", "S"), !is.na(timepoint))%>%
  mutate(class2=case_when(infectiontype=="A"&id%in%individuals_with_parasites_at_day14$id~"parasites at day 14",
                          infectiontype=="A"&id%notin%individuals_with_parasites_at_day14$id~"no parasites at day 14", .default = NULL))%>%
  ggplot(., aes(x=class2, y=freq/100, fill=timepoint, color=class2))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = F)+
  scale_y_continuous(labels=scales::label_percent())+
  scale_fill_manual(values=time_cols)+
  scale_color_manual(values=c("black", "darkred", "black"))+
  facet_wrap(~gate+infectiontype+stim, scales="free", nrow=2)+
  theme_minimal()



long_combo%>%
  filter(gate %in% c("Tregs_PD1_Frequency", "cTfh_Tr1_Frequency", "Th1_CD25_Frequency"), stim=="unstim")%>%
  distinct(id, infectiontype, timepoint, gate, stim, freq)%>%
  filter(infectiontype%in%c("A", "S"), !is.na(timepoint))%>%
  mutate(class2=case_when(infectiontype=="A"&id%in%individuals_with_parasites_at_day14$id~"parasites at day 14",
                          infectiontype=="A"&id%notin%individuals_with_parasites_at_day14$id~"no parasites at day 14", .default = NULL))%>%
  ggplot(., aes(x=class2, y=freq/100, fill=timepoint, color=class2))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = F)+
  scale_y_continuous(labels=scales::label_percent())+
  scale_fill_manual(values=time_cols)+
  scale_color_manual(values=c("black", "darkred", "black"))+
  facet_wrap(~gate+infectiontype, scales="free", ncol=2)+
  theme_minimal()

# intracellular cytokine stuff in stim ####
t_freq_purf <- long_combo%>%
  filter(gate %in% c("Tr1_IL10_Frequency", "Tr1_IFNg_Frequency", "Tr1_IL21_Frequency", "IL10_Frequency", "IFNg_Frequency", "IL21_Frequency"), stim!="unstim")%>%
  distinct(id, infectiontype, gender_categorical, ageyrs, timepoint, gate, stim, freq)%>%
  filter(infectiontype %in% c("A", "S"))%>%
  mutate(id_cat=paste("id_", id, sep=""))%>%
  group_by(gate, stim)%>%
  nest() %>%
  mutate(model=map(data, ~lme4::lmer(freq~timepoint*infectiontype+ageyrs+gender_categorical+(1|id_cat), data=.))) %>%
  # mutate(model=map(data, ~lme4::lmer(freq~timepoint*infectiontype+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(emm=map(model, ~emmeans(., specs = pairwise ~ timepoint | infectiontype)))%>%
  mutate(emm2=map(model, ~emmeans(., specs = pairwise ~ infectiontype | timepoint)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  mutate("baseline A - baseline S"=map_dbl(emm_contrast_summary2, ~.$p.value[1])) %>%
  mutate("baseline A - day0 A"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("baseline A - day14 A"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("baseline S - day0 S"=map_dbl(emm_contrast_summary, ~.$p.value[7])) %>%
  mutate("baseline S - day7 S"=map_dbl(emm_contrast_summary, ~.$p.value[9])) %>%
  mutate("baseline S - day14 S"=map_dbl(emm_contrast_summary, ~.$p.value[8]))%>%
  mutate("coef baseline A - baseline S"=map_dbl(emm_contrast_summary2, ~.$estimate[1])) %>%
  mutate("coef baseline A - day0 A"=map_dbl(emm_contrast_summary, ~.$estimate[1])) %>%
  mutate("coef baseline A - day14 A"=map_dbl(emm_contrast_summary, ~.$estimate[2])) %>%
  mutate("coef baseline S - day0 S"=map_dbl(emm_contrast_summary, ~.$estimate[7])) %>%
  mutate("coef baseline S - day7 S"=map_dbl(emm_contrast_summary, ~.$estimate[9])) %>%
  mutate("coef baseline S - day14 S"=map_dbl(emm_contrast_summary, ~.$estimate[8]))%>%
  pivot_longer(cols=starts_with("baseline"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))%>%
  ungroup()%>%
  pivot_longer(cols=starts_with("coef"), names_to = "coef_name", values_to = "coef")%>%
  rowwise()%>%
  filter(grepl(contrast, coef_name))

t_freq_purf%>%
  filter(padj<0.05)%>%
  select(gate, contrast, stim, padj)%>%
  print(n=40)

sig_freqs <- t_freq_purf%>%
  filter(padj<0.05)%>%
  mutate(infectiontype=substr(contrast, nchar(contrast)-1, nchar(contrast)))

sig_freqs_with_data <- long_combo%>%
  # filter(gate %in% c("Tr1_IL10_Frequency", "Tr1_IFNg_Frequency", "IL10_Frequency", "IFNg_Frequency"), stim!="unstim")%>%
  distinct(id, infectiontype, timepoint, gate, stim, freq)%>%
  filter(infectiontype %in% c("A", "S"))%>%
  full_join(., sig_freqs, by=c("gate", "stim", "infectiontype"))%>%
  select(id, infectiontype, timepoint, gate, stim, freq)


# phenotypic unstim
t_freq_purf2 <- long_combo%>%
  filter(gate %in% c("Th17_Frequency", "Th2_Frequency", "Th1_Frequency", "CD49b_LAG3_Frequency", "Tr1_Frequency", "cTfh_Frequency", "CXCR5_Frequency", "FOXP3_Frequency", "Tregs_Frequency"), stim=="unstim")%>%
  distinct(id, infectiontype, gender_categorical, ageyrs, timepoint, gate, stim, freq)%>%
  filter(infectiontype %in% c("A", "S"))%>%
  mutate(id_cat=paste("id_", id, sep=""))%>%
  group_by(gate, stim)%>%
  nest() %>%
  mutate(model=map(data, ~lme4::lmer(freq~timepoint*infectiontype+ageyrs+gender_categorical+(1|id_cat), data=.))) %>%
  # mutate(model=map(data, ~lme4::lmer(freq~timepoint*infectiontype+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(emm=map(model, ~emmeans(., specs = pairwise ~ timepoint | infectiontype)))%>%
  mutate(emm2=map(model, ~emmeans(., specs = pairwise ~ infectiontype | timepoint)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise", adjust="none")))%>%
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

t_freq_purf2%>%
  filter(padj<0.05)%>%
  select(gate, contrast, stim, padj)%>%
  print(n=40)
# parasite controllers functional ####

t_freq_ctrl_purf <- long_combo%>%
  filter(gate %in% c("Tr1_IL10_Frequency", "Tr1_IFNg_Frequency", "Tr1_IL21_Frequency", "IL10_Frequency", "IFNg_Frequency", "IL21_Frequency"), stim!="unstim")%>%
  distinct(id, gender_categorical, ageyrs, infectiontype, timepoint, gate, stim, freq)%>%
  filter(infectiontype %in% c("A"))%>%
  mutate(id_cat=paste("id_", id, sep=""),
         class2=ifelse(id %in% individuals_with_parasites_at_day14$id, "susceptible", "controller"))%>%
  group_by(gate, stim,)%>%
  nest() %>%
  # mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id_cat), data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(freq~timepoint*class2+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(emm=map(model, ~emmeans(., specs = pairwise ~ timepoint | class2)))%>%
  mutate(emm2=map(model, ~emmeans(., specs = pairwise ~ class2 | timepoint)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(.,"pairwise", adjust="none")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  mutate("baseline A - baseline S"=map_dbl(emm_contrast_summary2, ~.$p.value[1])) %>%
  mutate("baseline - day0 no parasites"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("baseline - day14 no parasites"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("baseline - day0 yes parasites"=map_dbl(emm_contrast_summary, ~.$p.value[4])) %>%
  mutate("baseline - day14 yes parasites"=map_dbl(emm_contrast_summary, ~.$p.value[5])) %>%
  mutate("baseline day0 - day14 yes parasites"=map_dbl(emm_contrast_summary, ~.$p.value[6])) %>%
  mutate("coef baseline - day0 no parasites"=map_dbl(emm_contrast_summary, ~.$estimate[1])) %>%
  mutate("coef baseline - day14 no parasites"=map_dbl(emm_contrast_summary, ~.$estimate[2])) %>%
  mutate("coef baseline - day0 no parasites"=map_dbl(emm_contrast_summary, ~.$estimate[4])) %>%
  mutate("coef baseline - day14 no parasites"=map_dbl(emm_contrast_summary, ~.$estimate[5])) %>%
  mutate("coef day0 - day14 yes parasites"=map_dbl(emm_contrast_summary, ~.$estimate[6])) %>%
  pivot_longer(cols=starts_with("baseline"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))%>%
  ungroup()%>%
  pivot_longer(cols=starts_with("coef"), names_to = "coef_name", values_to = "coef")%>%
  rowwise()%>%
  filter(grepl(contrast, coef_name))

t_freq_ctrl_purf%>%
  select(gate, contrast, p, padj)%>%
  arrange(p)

# parasite controllers phenotypic####
t_freq_ctrl_purf2 <- long_combo%>%
  filter(gate %in% c("Th17_Frequency", "Th2_Frequency", "Th1_Frequency", "CD49b_LAG3_Frequency", "Tr1_Frequency", "cTfh_Frequency", "CXCR5_Frequency", "FOXP3_Frequency", "Tregs_Frequency"), stim=="unstim")%>%
  distinct(id, gender_categorical, ageyrs, infectiontype, timepoint, gate, stim, freq)%>%
  filter(infectiontype %in% c("A"))%>%
  mutate(id_cat=paste("id_", id, sep=""),
         class2=ifelse(id %in% individuals_with_parasites_at_day14$id, "susceptible", "controller"))%>%
  group_by(gate, stim,)%>%
  nest() %>%
  # mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id_cat), data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(freq~timepoint*class2+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(emm=map(model, ~emmeans(., specs = pairwise ~ timepoint | class2)))%>%
  mutate(emm2=map(model, ~emmeans(., specs = pairwise ~ class2 | timepoint)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(.,"pairwise", adjust="none")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  mutate("baseline no parasites - baseline yes parasites"=map_dbl(emm_contrast_summary2, ~.$p.value[1])) %>%
  mutate("baseline - day0 no parasites"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("baseline - day14 no parasites"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("baseline - day0 yes parasites"=map_dbl(emm_contrast_summary, ~.$p.value[4])) %>%
  mutate("baseline - day14 yes parasites"=map_dbl(emm_contrast_summary, ~.$p.value[5])) %>%
  mutate("baseline day0 - day14 yes parasites"=map_dbl(emm_contrast_summary, ~.$p.value[6])) %>%
  mutate("coef baseline no parasites - baseline yes parasites"=map_dbl(emm_contrast_summary, ~.$estimate[1])) %>%
  mutate("coef baseline - day0 no parasites"=map_dbl(emm_contrast_summary, ~.$estimate[1])) %>%
  mutate("coef baseline - day14 no parasites"=map_dbl(emm_contrast_summary, ~.$estimate[2])) %>%
  mutate("coef baseline - day0 no parasites"=map_dbl(emm_contrast_summary, ~.$estimate[4])) %>%
  mutate("coef baseline - day14 no parasites"=map_dbl(emm_contrast_summary, ~.$estimate[5])) %>%
  mutate("coef day0 - day14 yes parasites"=map_dbl(emm_contrast_summary, ~.$estimate[6])) %>%
  pivot_longer(cols=starts_with("baseline"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))%>%
  ungroup()%>%
  pivot_longer(cols=starts_with("coef"), names_to = "coef_name", values_to = "coef")%>%
  rowwise()%>%
  filter(grepl(contrast, coef_name))

t_freq_ctrl_purf2%>%
  select(gate, contrast, p, padj)%>%
  arrange(p)
# figures ####


functional_s <- long_combo%>%
  filter(gate %in% c("IFNg_Frequency", "IL10_Frequency",  "Tr1_IFNg_Frequency", "Tr1_IL10_Frequency"), stim=="iRBCs")%>%
  mutate(gate=gsub("_", " ", gate))%>%
  distinct(id, infectiontype, timepoint, gate, stim, freq)%>%
  filter(infectiontype%in%c("S"), timepoint%in%c("baseline", "day7"))%>%
  # mutate(class2=case_when(infectiontype=="A"&id%in%individuals_with_parasites_at_day14$id~"parasites at day 14",
  #                         infectiontype=="A"&id%notin%individuals_with_parasites_at_day14$id~"no parasites at day 14", .default = NULL))%>%
  ggplot(., aes(x=timepoint, y=freq/100, fill=timepoint))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = F)+
  scale_y_continuous(labels=scales::label_percent())+
  scale_fill_manual(values=time_cols)+
  facet_wrap(~gate+stim, scales="free", ncol=4)+
  theme_minimal(12)+
  ylab("% of Parent")+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/abstracts/keystone25/poster_figures/functional_s.png", functional_s, width=10, height=4, bg='white')

functional_a <- long_combo%>%
  filter(gate %in% c("IFNg_Frequency", "IL10_Frequency",  "Tr1_IFNg_Frequency", "Tr1_IL10_Frequency"), stim=="iRBCs")%>%
   distinct(id, infectiontype, timepoint, gate, stim, freq)%>%
  mutate(gate=gsub("_", " ", gate))%>%
  filter(infectiontype%in%c("A"), timepoint%in%c("baseline", "day0", "day14"))%>%
  mutate(class2=case_when(infectiontype=="A"&id%in%individuals_with_parasites_at_day14$id~"parasites at day 14",
                          infectiontype=="A"&id%notin%individuals_with_parasites_at_day14$id~"no parasites at day 14", .default = NULL))%>%
  ggplot(., aes(x=timepoint, y=freq/100, fill=timepoint, color=class2))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = F)+
  scale_y_continuous(labels=scales::label_percent())+
  scale_fill_manual(values=time_cols)+
  scale_color_manual(values=c("black", "darkred"))+
  facet_wrap(~gate+stim, scales="free", ncol=4)+
  theme_minimal()+
  guides(fill=guide_none())+
  ylab("% of Parent")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/abstracts/keystone25/poster_figures/functional_a.png", functional_a, width=10, height=4, bg='white')


 functional_combo <- long_combo%>%
  filter(gate %in% c("Tr1_IFNg_Frequency", "Tr1_IL10_Frequency", "Tr1_IL21_Frequency"), stim=="iRBCs")%>%
  mutate(gate=gsub("_", " ", gate))%>%
  distinct(id, infectiontype, timepoint, gate, stim, freq)%>%
  filter(infectiontype%in%c("S", "A"), timepoint%in%c("baseline", "day7", "day14"))%>%
  mutate(class2=case_when(infectiontype=="A"&id%in%individuals_with_parasites_at_day14$id~"parasites at day 14",
                          infectiontype=="A"&id%notin%individuals_with_parasites_at_day14$id~"no parasites at day 14", .default = NULL))%>%
  mutate(timepoint=factor(timepoint, levels=c("baseline", "day7", "day14")))%>%
  ggplot(., aes(x=timepoint, y=freq/100, fill=timepoint, color=class2))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = F)+
  scale_y_continuous(labels=scales::label_percent())+
  scale_fill_manual(values=time_cols)+
  scale_color_manual(values=c("black", "darkred"), na.value = "black")+
  facet_wrap(~infectiontype+gate, scales="free", ncol=3)+
  theme_minimal(12)+
  ylab("% of Parent")+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/abstracts/keystone25/poster_figures/functional_combo.png", functional_combo, width=8, height=6, bg='white')


phenotypic_as <- long_combo%>%
  filter(gate %in% c("Th1_Frequency", "Tr1_Frequency"), stim=="unstim")%>%
  distinct(id, infectiontype, timepoint, gate, stim, freq)%>%
  filter(infectiontype%in%c("A", "S"), timepoint%in%c("baseline", "day0", "day7", "day14"))%>%
  mutate(class2=case_when(infectiontype=="A"&id%in%individuals_with_parasites_at_day14$id~"parasites at day 14",
                          infectiontype=="A"&id%notin%individuals_with_parasites_at_day14$id~"no parasites at day 14", .default = NULL))%>%
  ggplot(., aes(x=timepoint, y=freq/100, fill=timepoint, color=class2))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = F)+
  ggpubr::stat_compare_means()+
  scale_y_continuous(labels=scales::label_percent())+
  scale_fill_manual(values=time_cols)+
  scale_color_manual(values=c("black", "darkred"))+
  facet_wrap(~infectiontype+gate, scales="free")+
  theme_minimal()+
  guides(fill=guide_none())+
  ylab("% of Parent")+
  theme(legend.position = "bottom",
        legend.title = element_blank())

ggsave("~/postdoc/stanford/abstracts/keystone25/poster_figures/phenotypic_as.png", functional_a, width=5, height=20, bg='white')





long_combo%>%
  filter(stim=="unstim")%>%
  distinct(id, infectiontype, timepoint, gate, stim, freq)%>%
  filter(infectiontype%in%c("A", "S"), timepoint%in%c("baseline"))%>%
  mutate(class2=case_when(infectiontype=="A"&id%in%individuals_with_parasites_at_day14$id~"parasites at day 14",
                          infectiontype=="A"&id%notin%individuals_with_parasites_at_day14$id~"no parasites at day 14", .default = NULL))%>%
  ggplot(., aes(x=timepoint, y=freq/100, fill=infectiontype))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = F)+
  scale_y_continuous(labels=scales::label_percent())+
  # scale_fill_manual(values=time_cols)+
  scale_fill_manual(values=c("black", "darkred"))+
  facet_wrap(~gate, scales="free")+
  theme_minimal()+
  guides(fill=guide_none())+
  ylab("% of Parent")+
  theme(legend.position = "bottom",
        legend.title = element_blank())



long_combo%>%
  filter(stim=="unstim")%>%
  distinct(id, infectiontype, timepoint, gate, stim, freq)%>%
  filter(infectiontype%in%c("A"), timepoint%in%c("baseline"))%>%
  mutate(class2=case_when(infectiontype=="A"&id%in%individuals_with_parasites_at_day14$id~"parasites at day 14",
                          infectiontype=="A"&id%notin%individuals_with_parasites_at_day14$id~"no parasites at day 14", .default = NULL))%>%
  ggplot(., aes(x=timepoint, y=freq/100, fill=class2))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = F)+
  scale_y_continuous(labels=scales::label_percent())+
  # scale_fill_manual(values=time_cols)+
  scale_fill_manual(values=c("black", "darkred"))+
  facet_wrap(~gate, scales="free")+
  theme_minimal()+
  guides(fill=guide_none())+
  ylab("% of Parent")+
  theme(legend.position = "bottom",
        legend.title = element_blank())
