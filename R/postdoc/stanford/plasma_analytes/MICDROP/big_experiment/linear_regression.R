library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(emmeans)

`%notin%` <- Negate(`%in%`)

clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks", "68 weeks")))%>%
  filter(targetName %notin% c("CTSS", "LTA|LTB", "IFNA2"))


# modelz ####
## change through time ####
time_purf <- clean_data%>%
  filter(mstatus==0)%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(conc~timepoint+gender_categorical+log_qpcr+(1|id), data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ timepoint)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate("8 weeks - 24 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("8 weeks - 52 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("24 weeks - 52 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  pivot_longer(cols=ends_with("weeks"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

sig_time <- time_purf%>%
  filter(padj < 0.05)

sig_8_24 <- sig_time%>%
  filter(contrast=="8 weeks - 24 weeks")%>%
  arrange(padj)

sig_8_52 <- sig_time%>%
  filter(contrast=="8 weeks - 52 weeks")%>%
  arrange(padj)

sig_24_52 <- sig_time%>%
  filter(contrast=="24 weeks - 52 weeks")%>%
  arrange(padj)

for(i in 1:ceiling(nrow(sig_8_24)/16)){
  
  plt <- clean_data %>%
    filter(targetName %in% sig_8_24$targetName[seq((i-1)*16+1, i*16)],
           timepoint!="68 weeks")%>%
    mutate(targetNamef=factor(targetName, levels=sig_8_24$targetName[seq((i-1)*16+1, i*16)]))%>%
    ggplot(aes(x=factor(timepoint_num), y=conc, fill=timepoint))+
    geom_line(aes(group=id), alpha=0.2)+
    geom_boxplot(outliers = FALSE)+
    facet_wrap(~targetNamef, scales = "free")+
    scale_fill_manual(values=viridis::viridis(3))+
    xlab("age in weeks")+
    theme_minimal()+
    theme(legend.position = "none",
          axis.title.y = element_blank())
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/analytes_through_time/sig_8_24", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}


for(i in 1:ceiling(nrow(sig_8_52)/16)){
  
  plt <- clean_data %>%
    filter(targetName %in% sig_8_52$targetName[seq((i-1)*16+1, i*16)],
           timepoint!="68 weeks")%>%
    mutate(targetNamef=factor(targetName, levels=sig_8_52$targetName[seq((i-1)*16+1, i*16)]))%>%
    ggplot(aes(x=factor(timepoint_num), y=conc, fill=timepoint))+
    geom_line(aes(group=id), alpha=0.2)+
    geom_boxplot(outliers = FALSE)+
    facet_wrap(~targetNamef, scales = "free")+
    scale_fill_manual(values=viridis::magma(3))+
    xlab("age in weeks")+
    theme_minimal()+
    theme(legend.position = "none",
          axis.title.y = element_blank())
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/analytes_through_time/sig_8_52", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}


for(i in 1:ceiling(nrow(sig_24_52)/16)){
  
  plt <- clean_data %>%
    filter(targetName %in% sig_24_52$targetName[seq((i-1)*16+1, i*16)],
           timepoint!="68 weeks")%>%
    mutate(targetNamef=factor(targetName, levels=sig_24_52$targetName[seq((i-1)*16+1, i*16)]))%>%
    ggplot(aes(x=factor(timepoint_num), y=conc, fill=timepoint))+
    geom_line(aes(group=id), alpha=0.2)+
    geom_boxplot(outliers = FALSE)+
    facet_wrap(~targetNamef, scales = "free")+
    scale_fill_manual(values=viridis::magma(3))+
    xlab("age in weeks")+
    theme_minimal()+
    theme(legend.position = "none",
          axis.title.y = element_blank())
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/analytes_through_time/sig_24_52", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}

## change through mstatus ####
mstatus_purf <- clean_data%>%
  filter(ageinwks<60, targetName %notin% c("CTSS", "LTA|LTB", "IFNA2"))%>%
  mutate(mstatusf=factor(case_when(mstatus==0&qPCRparsdens>0.5~"asymptomatic",
                                   mstatus==0&qPCRparsdens==0 ~"uninfected",
                                   mstatus==1~"symptomatic", .default = "unsure"
                                   )))%>%
  filter(mstatusf !="unsure")%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(conc~timepoint*mstatusf+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ mstatusf | timepoint)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate("symp vs uninfected 8 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  mutate("symp vs uninfected 24 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[6])) %>%
  mutate("symp vs uninfected 52 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[9])) %>%
  mutate("asymp vs uninfected 8 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("asymp vs uninfected 24 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[5])) %>%
  mutate("asymp vs uninfected 52 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[8])) %>%
  mutate("asymp vs symp 8 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("asymp vs symp 24 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[4])) %>%
  mutate("asymp vs symp 52 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[7])) %>%
  pivot_longer(cols=ends_with("weeks"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

# mstatus_purf%>%
#   select(targetName, contrast, p, padj)%>%
#   write.csv(., "~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/mstatus_purf.csv", row.names = F)

sig_malaria <- mstatus_purf%>%
  filter(padj < 0.05, grepl("^symp", contrast))

sig_asymp <- mstatus_purf%>%
  filter(padj < 0.05, grepl("^asymp vs uninfected", contrast))

asymp_specific <- unique(sig_asymp$targetName)[unique(sig_asymp$targetName) %notin% unique(sig_malaria$targetName)]

# 
# sig_8 <- sig_mstatus%>%
#   filter(contrast=="8 weeks")%>%
#   arrange(padj)
# 
# sig_24 <- sig_mstatus%>%
#   filter(contrast=="24 weeks")%>%
#   arrange(padj)
# 
# sig_52 <- sig_mstatus%>%
#   filter(contrast=="52 weeks")%>%
#   arrange(padj)
unique_sigs <- unique(sig_mstatus$targetName)

for(i in 1:ceiling(length(unique_sigs)/16)){
  
  plt <- clean_data %>%
    filter(targetName %in% unique_sigs[seq((i-1)*16+1, i*16)],
           timepoint!="68 weeks")%>%
    mutate(targetNamef=factor(targetName, levels=unique_sigs[seq((i-1)*16+1, i*16)]))%>%
    ggplot(aes(x=factor(timepoint_num), y=conc, fill=factor(mstatus)))+
    # geom_line(aes(group=id), alpha=0.2)+
    geom_boxplot(outliers = FALSE)+
    facet_wrap(~targetNamef, scales = "free")+
    scale_fill_manual(values=c("darkgrey", "darkred"))+
    xlab("age in weeks")+
    theme_minimal()+
    theme(legend.position = "none",
          axis.title.y = element_blank())
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/mstatus_analytes", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}

age_change_mstatus_plot <- clean_data %>%
  filter(targetName %in% c("CD70", "IL18", "CTLA4", "LAG3", "GZMA", "IFNG", "CCL2", "CCL8", "IL10 "),
         timepoint!="68 weeks")%>%
  ggplot(aes(x=factor(timepoint_num), y=conc, fill=factor(mstatus)))+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetName, scales = "free", nrow=2)+
  scale_fill_manual(values=c("darkgrey", "darkred"))+
  xlab("age in weeks")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/age_change_mstatus_plot.png", age_change_mstatus_plot, height=4, width=8, dpi=444, bg="white")

#chemokine decreases = unmoved during asymp; result of lack of consumption by T cells?
age_change_asymp_plot <- clean_data %>%
  filter(targetName %in% asymp_specific,
         #c("SDC1", "LILRB2", "IL10", "VCAM1", "LAG3"),
         timepoint!="68 weeks")%>%
  mutate(mstatus=if_else(mstatus==0&qPCRparsdens>1, 0.5, mstatus))%>%
  ggplot(aes(x=factor(timepoint_num), y=conc, fill=factor(mstatus)))+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetName, scales = "free", nrow=4)+
  scale_fill_manual(values=c("darkgrey", "tan1", "darkred"))+
  xlab("age in weeks")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/age_change_mstatus_plot.png", age_change_mstatus_plot, height=4, width=8, dpi=444, bg="white")


# add n_infection etc ###

mic_drop <-  haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/Jun25/MICDSpecimenBoxJun25_withclinical.dta")

mic_drop_data <- mic_drop %>%
  filter(mstatus != 0, !is.na(mstatus), incidentmalaria==1 | is.na(incidentmalaria) & mstatus==2 | incidentmalaria==0 & mstatus==2)%>%
  group_by(id) %>%
  add_count(name="total_n_infection") %>%
  mutate(n_infection = seq(1, max(total_n_infection)),
         mstatus = case_match(mstatus,
                              0~"no malaria",
                              1~"uncomplicated",
                              2~"complicated",
                              3~"quinine for AL failure",
                              4~"Q/AS failure"))%>%
  mutate(date=as.character(date))%>%
  select(id, date, n_infection)

infections <- clean_data%>%
  filter(mstatus==1)%>%
  distinct(id, date, ageinwks)


n_infections <- left_join(infections, mic_drop_data, by=c("id", "date"))


analytes <- unique(clean_data$targetName)

nulisa_to_show <- clean_data %>%
  left_join(., n_infections, by=c("id", "date"))%>%
  mutate(n_infection=ifelse(n_infection>=4, "4+", n_infection))%>%
  filter(targetName %in% c("CTLA4", "CD274", "CD80", "PDCD1"))%>%
  mutate(targetNamef=factor(targetName, levels=c("CTLA4", "CD274", "CD80", "PDCD1")))%>%
  ggplot(aes(x=factor(n_infection), y=conc, fill=factor(mstatus)))+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetNamef, scales = "free", nrow=2)+
  scale_fill_manual(values=viridis::rocket(n = 5))+
  xlab("Order of infection")+
  ylab("Concentration in Blood")+
  theme_minimal(base_size = 16)+
  theme(legend.position = "none")

ggsave("~/Downloads/nulisa_for_malaria_day.png", nulisa_to_show, height=6, width=6, dpi = 444, bg="white")



nulisa_to_show2 <- clean_data %>%
  left_join(., n_infections, by=c("id", "date"))%>%
  mutate(n_infection=ifelse(n_infection>=4, "4+", n_infection))%>%
  filter(!is.na(n_infection))%>%
  filter(targetName %in% c("CD70", "IL18", "CTLA4", "LAG3", "GZMA", "IFNG", "CCL2", "CCL8", "IL10 "))%>%
  # mutate(targetNamef=factor(targetName, levels=c("CD70", "IL18", "CTLA4", "LAG3", "GZMA", "IFNG", "CCL2", "CCL8", "IL10 ")))%>%
  ggplot(aes(x=factor(n_infection), y=conc, fill=factor(mstatus)))+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetName, scales = "free", nrow=2)+
  scale_fill_manual(values=c("darkred"))+
  xlab("order of malaria episode")+
  ylab("concentration in blood")+
  theme_minimal(base_size = 16)+
  theme(legend.position = "none")

ggsave("~/Downloads/nulisa_for_boyd_lab.png", nulisa_to_show2, height=4, width=8, dpi = 444, bg="white")


n_infection_purf <- clean_data%>%
  left_join(., n_infections, by=c("id", "date"))%>%
  filter(mstatus==1, targetName %notin% c("CTSS", "LTA|LTB", "IFNA2"))%>%
  mutate(n_infection_cat=case_when(n_infection<4~as.character(n_infection),
                                   n_infection>3~"4+"))%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(time_model=map(data, ~lm(conc~factor(n_infection_cat), data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ n_infection_cat)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate("1 vs 2"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("1 vs 3"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("1 vs 4+"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  pivot_longer(cols=starts_with("1"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

n_infection_sigs <-n_infection_purf%>%
  filter(p<0.01)%>%
  arrange(p)


clean_data%>%
  left_join(., n_infections, by=c("id", "date"))%>%
  filter(targetName %in% grep("^IL*", n_infection_purf$targetName, value=T),
         mstatus==1)%>%
  mutate(n_infection_cat=case_when(n_infection<4~as.character(n_infection),
                                   n_infection>3~"4+"))%>%
  mutate(targetNamef=factor(targetName))%>%
  ggplot(aes(x=factor(n_infection_cat), y=conc, fill=factor(n_infection_cat)))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetNamef, scales = "free")+
  viridis::scale_fill_viridis(discrete = T)+
  xlab("number of malaria episodes")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_blank())


# IL10, IL22, CTLA4, CD70, CD274, PDCD1
mstatus_targets_of_interest <- c("IL18", "IL27", "IFNG",
                                 "IL18BP", "NTF3", "FASLG")

analytes_during_malaria_n_infection <- clean_data%>%
  left_join(., n_infections, by=c("id", "date"))%>%
  filter(targetName %in% mstatus_targets_of_interest,
         mstatus==1)%>%
  mutate(n_infection_cat=case_when(n_infection<4~as.character(n_infection),
                                   n_infection>3~"4+"))%>%
  mutate(targetNamef=factor(targetName, levels=mstatus_targets_of_interest))%>%
  ggplot(aes(x=factor(n_infection_cat), y=conc, fill=factor(n_infection_cat)))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetNamef, scales = "free")+
  viridis::scale_fill_viridis(discrete = T)+
  xlab("number of malaria episodes")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/analytes_during_malaria_n_infection.png", analytes_during_malaria_n_infection, height=4, width=8, dpi=444, bg="white")


analytes_during_malaria_age <- clean_data%>%
  filter(targetName %in% mstatus_targets_of_interest,
         , timepoint!="68 weeks")%>%
  mutate(targetNamef=factor(targetName, levels=mstatus_targets_of_interest))%>%
  ggplot(aes(x=factor(timepoint), y=conc, fill=factor(mstatus)))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetNamef, scales = "free")+
  scale_fill_manual(values=c("darkgrey", "darkred"))+
  xlab("age")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/analytes_during_malaria_age.png", analytes_during_malaria_age, height=4, width=8, dpi=444, bg="white")

## change through infection history####
### negative binomial; this is needed for para because variance = 4 and mean =1 ####
purf_52 <- clean_data%>%
  filter(mstatus==0, timepoint=="52 weeks")%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~MASS::glm.nb(total_n_malaria_12~conc+gender_categorical, data=.))) %>%
  mutate(n_para_model=map(data, ~MASS::glm.nb(total_n_para_12~conc+gender_categorical, data=.))) %>%
  # mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(total_n_malaria_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  # mutate(n_para_model=map(data, ~glmmTMB::glmmTMB(total_n_para_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  mutate(n_para_model_summary=map(n_para_model, ~summary(.)))%>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[11]))%>%
  mutate(n_para_p=map_dbl(n_para_model_summary, ~coef(.)[11]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"),
         n_para_padj=p.adjust(n_para_p, method="BH"))

purf_52%>%
  select(targetName, n_malaria_padj, n_para_padj)%>%
  write.csv(., "~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/n_infection_purf.csv", row.names = F)

sig_mala_12 <- purf_52%>%
  filter(n_malaria_padj<0.1)

sig_para_12 <- purf_52%>%
  filter(n_para_padj<0.1)

for(i in 1:ceiling(length(sig_mala_12)/16)){
  
  plt <- clean_data %>%
    filter(targetName %in% sig_mala_12$targetName,
           timepoint=="52 weeks")%>%
    mutate(targetNamef=factor(targetName, levels=sig_mala_12$targetName))%>%
    ggplot(aes(x=factor(total_n_malaria_12), y=conc, fill=factor(total_n_malaria_12)))+
    # geom_line(aes(group=id), alpha=0.2)+
    geom_boxplot(outliers = FALSE)+
    facet_wrap(~targetNamef, scales = "free")+
    viridis::scale_fill_viridis(discrete = T)+
    xlab("number of malaria episodes")+
    theme_minimal()+
    theme(legend.position = "none",
          axis.title.y = element_blank())
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/n_mala_no_pcr_analytes", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=4, width=8, dpi=444, bg="white")
  
}

for(i in 1:ceiling(length(sig_para_12)/16)){
  
  plt <- clean_data %>%
    filter(targetName %in% sig_para_12$targetName,
           timepoint=="52 weeks")%>%
    mutate(targetNamef=factor(targetName, levels=sig_para_12$targetName))%>%
    ggplot(aes(x=factor(total_n_para_12), y=conc, fill=factor(total_n_para_12)))+
    # geom_line(aes(group=id), alpha=0.2)+
    geom_boxplot(outliers = FALSE)+
    facet_wrap(~targetNamef, scales = "free", nrow=2)+
    viridis::scale_fill_viridis(discrete = T)+
    xlab("number of parasitemic months")+
    theme_minimal()+
    theme(legend.position = "none",
          axis.title.y = element_blank())
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/n_para_no_pcr_analytes", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=12, dpi=444, bg="white")
  
}
# 
# micdrop_nulisa_with_meta%>%
#   filter(mstatus==0, timepoint=="52 weeks")%>%
#   filter(targetName%in%sig_para_12$targetName)%>%
#   ggplot(., aes(x=total_n_para_12, y=conc, fill=factor(total_n_para_12)))+
#   facet_wrap(~targetName, scales="free")+
#   geom_boxplot(outliers = F)+
#   viridis::scale_fill_viridis(discrete = T)+
#   theme_minimal()+
#   theme(legend.position = "none")
# 
# micdrop_nulisa_with_meta%>%
#   filter(mstatus==0, timepoint=="52 weeks")%>%
#   filter(targetName%in%sig_mala_12$targetName)%>%
#   ggplot(., aes(x=total_n_malaria_12, y=conc, fill=factor(total_n_malaria_12)))+
#   facet_wrap(~targetName, scales="free")+
#   geom_boxplot(outliers = F)+
#   viridis::scale_fill_viridis(discrete = T)+
#   theme_minimal()+
#   theme(legend.position = "none")


### logistic; any malaria / para ####
purf_52_bino <- clean_data%>%
  filter(mstatus==0, timepoint=="52 weeks")%>%
  group_by(targetName)%>%
  mutate(any_malaria=ifelse(total_n_malaria_12==0, 0, 1))%>%
  mutate(any_para=ifelse(total_n_para_12==0, 0, 1))%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~glm(any_malaria~conc+gender_categorical, data=., family="binomial"))) %>%
  mutate(n_para_model=map(data, ~glm(any_para~conc+gender_categorical, data=., family="binomial"))) %>%
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  mutate(n_para_model_summary=map(n_para_model, ~summary(.)))%>%
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[11]))%>%
  mutate(n_para_p=map_dbl(n_para_model_summary, ~coef(.)[11]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"),
         n_para_padj=p.adjust(n_para_p, method="BH"))

sig_mala_bino_12 <- purf_52_bino%>%
  filter(n_malaria_padj<0.1)

sig_para_bino_12 <- purf_52_bino%>%
  filter(n_para_padj<0.1)



clean_data%>%
  mutate(any_malaria=ifelse(total_n_malaria_12==0, 0, 1))%>%
  filter(mstatus==0, timepoint=="52 weeks")%>%
  filter(targetName%in%sig_mala_12$targetName)%>%
  ggplot(., aes(x=any_malaria, y=conc, fill=factor(any_malaria)))+
  facet_wrap(~targetName, scales="free")+
  geom_boxplot(outliers = F)+
  viridis::scale_fill_viridis(discrete = T)+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/Downloads/any_malar.png", width=8, height=4, dpi=444)
### normal timepoint model with infection history as categorical variable ####
any_para_purf <- clean_data%>%
  filter(mstatus==0, treatmentarm!="DP 2 years")%>%
  mutate(any_para=ifelse(total_n_para_12==0, 0, 1))%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(conc~timepoint*any_para+(1|id), data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ any_para | timepoint)))%>%
  mutate(emm2=map(time_model, ~emmeans(., specs = pairwise ~ timepoint | any_para)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  mutate("8 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("24 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("52 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  pivot_longer(cols=ends_with("weeks"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))


any_para_nulisa_fc_data <- clean_data %>%
  mutate(any_para=ifelse(total_n_para_12==0, "none", "some"))%>%
  filter(treatmentarm != "DP 2 years", timepoint !=c("68 weeks"))%>%
  group_by(targetName, timepoint, any_para)%>%
  summarise("mean_conc"=mean(conc, na.rm = T))%>%
  pivot_wider(names_from = any_para, values_from = mean_conc)%>%
  mutate(fc=some-none)



any_para_nulisa_purf_fc <- any_para_purf%>%
  mutate("timepoint"=contrast)%>%
  left_join(., any_para_nulisa_fc_data, by=c("targetName", "timepoint"))%>%
  mutate(contrast=factor(contrast, levels=c("8 weeks", "24 weeks", "52 weeks")))

any_para_nulisa_purf_fc%>%
  select(targetName, contrast, fc, p, padj)%>%
  write.csv(., "~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/any_para_regression.csv")

treatment_nulisa_sigs2 <- any_para_nulisa_purf_fc %>%
  filter(padj<0.1 & abs(fc)>=0.5)%>%
  select(targetName, contrast, p, padj)


any_mala_purf <- clean_data%>%
  filter(mstatus==0, treatmentarm!="DP 2 years")%>%
  mutate(any_malaria=ifelse(total_n_malaria_12==0, 0, 1))%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(conc~timepoint*any_malaria+(1|id), data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ any_malaria | timepoint)))%>%
  mutate(emm2=map(time_model, ~emmeans(., specs = pairwise ~ timepoint | any_malaria)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  mutate("8 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("24 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("52 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  pivot_longer(cols=ends_with("weeks"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))



any_mala_nulisa_fc_data <- clean_data %>%
  mutate(any_malaria=ifelse(total_n_malaria_12==0, "none", "some"))%>%
  filter(treatmentarm != "DP 2 years", timepoint !=c("68 weeks"))%>%
  group_by(targetName, timepoint, any_malaria)%>%
  summarise("mean_conc"=mean(conc, na.rm = T))%>%
  pivot_wider(names_from = any_malaria, values_from = mean_conc)%>%
  mutate(fc=some-none)


any_mala_purf_fc <- any_mala_purf%>%
  mutate("timepoint"=contrast)%>%
  left_join(., any_mala_nulisa_fc_data, by=c("targetName", "timepoint"))%>%
  mutate(contrast=factor(contrast, levels=c("8 weeks", "24 weeks", "52 weeks")))

any_mala_purf_fc%>%
  select(targetName, contrast, fc, p, padj)%>%
  write.csv(., "~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/any_malaria_regression.csv")


any_mala_sigs <- treatment_nulisa_purf_fc2 %>%
  filter(padj<0.1 & abs(fc)>=0.5 )
### poisson; this is OK for malaria because variance = 0.28 and mean =0.22####
purf_52_poisson <- clean_data%>%
  filter(mstatus==0, timepoint=="52 weeks")%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~glm(total_n_malaria_12~conc, data=., family="poisson"))) %>%
  mutate(n_para_model=map(data, ~glm(total_n_para_12~conc, data=., family="poisson"))) %>%
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  mutate(n_para_model_summary=map(n_para_model, ~summary(.)))%>%
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[8]))%>%
  mutate(n_para_p=map_dbl(n_para_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"),
         n_para_padj=p.adjust(n_para_p, method="BH"))

sig_mala_poisson_12 <- purf_52_poisson%>%
  filter(n_malaria_padj<0.1)

sig_para_poisson_12 <- purf_52_poisson%>%
  filter(n_para_padj<0.1)

clean_data%>%
  filter(mstatus==0, timepoint=="52 weeks")%>%
  filter(targetName%in%sig_mala_poisson_12$targetName)%>%
  ggplot(., aes(x=total_n_malaria_12, y=conc, fill=factor(total_n_malaria_12)))+
  facet_wrap(~targetName, scales="free")+
  geom_boxplot(outliers = F)+
  viridis::scale_fill_viridis(discrete = T)+
  theme_minimal()+
  theme(legend.position = "none")



## treatment arm ####

mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")

clean_data <- clean_data%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         treatmentarm=case_match(treatmentarm,
                    1~"Placebo",
                    2~"DP 1 year",
                    3~"DP 2 years"))

treatment_purf <- clean_data%>%
  filter(mstatus==0, treatmentarm!="DP 2 years")%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(conc~timepoint*treatmentarm+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ treatmentarm | timepoint)))%>%
  mutate(emm2=map(time_model, ~emmeans(., specs = pairwise ~ timepoint | treatmentarm)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  mutate("8 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("24 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("52 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  pivot_longer(cols=ends_with("weeks"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

sigs <- treatment_purf%>%
  filter(padj<0.1)%>%
  filter(contrast=="52 weeks")

time_treatment_plot <- clean_data%>%
  filter(treatmentarm!="DP 2 years", timepoint=="52 weeks", mstatus==0)%>%
  filter(targetName %in% sigs$targetName)%>%
  ggplot(aes(x=factor(timepoint_num), y=conc, fill=treatmentarm))+
  geom_violin(draw_quantiles = seq(0,1,0.25), color="white")+
  ggpubr::stat_compare_means(size=2, )+
  # geom_boxplot(outliers=F)+
  # # stat_summary(aes(group=c(treatmentarm)), geom="cro", fun.y=median, colour="white")+
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  facet_wrap(~targetName, scales="free", nrow=2)+
  scale_fill_manual(values=c("darkred", "#00555A"))+
  theme_minimal()+
  theme(legend.title = element_blank(),
        axis.title = element_blank())
  

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/time_treatment_plot.png", width=8, height=4, dpi=444, bg="white")

treatment_purf%>%
  select(targetName, contrast, p, padj)%>%
  write.csv(., "~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/nulisa_treatment_regression.csv")

# sickle trait ####
sickle_purf <- clean_data%>%
  filter(mstatus==0, timepoint!="68 weeks", hbs!="HbSS")%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(conc~timepoint*hbs+(1|id), data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ hbs | timepoint)))%>%
  mutate(emm2=map(time_model, ~emmeans(., specs = pairwise ~ timepoint | hbs)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  mutate("8 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("24 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("52 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  pivot_longer(cols=ends_with("weeks"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

sigs <- sickle_purf%>%
  filter(padj<0.1)


clean_data%>%
  filter(mstatus==0, timepoint!="68 weeks", hbs!="HbSS")%>%
  filter(targetName%in%c(sigs$targetName))%>%
  ggplot(., aes(x=timepoint, y=conc, fill=hbs))+
  ggpubr::stat_compare_means()+
  geom_point(position=position_dodge(width=0.75), alpha=0.1)+
  geom_boxplot(outliers = F)+
  theme_minimal()


# hierarchical clustering ####

