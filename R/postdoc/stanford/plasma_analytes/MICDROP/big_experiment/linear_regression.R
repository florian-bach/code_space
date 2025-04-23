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

## change through with mstatus ####
mstatus_purf <- clean_data%>%
  filter(ageinwks<60, targetName %notin% c("CTSS", "LTA|LTB", "IFNA2"))%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(conc~timepoint*mstatus+(1|id), data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ mstatus | timepoint)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate("8 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("24 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("52 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  pivot_longer(cols=ends_with("weeks"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))


sig_mstatus <- mstatus_purf%>%
  filter(padj < 0.05)
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
infections <- clean_data%>%
  filter(mstatus==1)%>%
  distinct(id, date, ageinwks)


# add n_infection etc

mic_drop <-  haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/Mar25/MICDSpecimenBoxMar25_withclinical.dta")

mic_drop_data <- mic_drop %>%
  filter(mstatus != 0, !is.na(mstatus))%>%
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

n_infections <- left_join(infections, mic_drop_data, by=c("id", "date"))


analytes <- unique(clean_data$targetName)

nulisa_to_show <- clean_data %>%
  left_join(., n_infections, by=c("id", "date"))%>%
  mutate(n_infection=ifelse(n_infection>=4, "4+", n_infection))%>%
  filter(targetName %in% c("CTLA4", "CD274", "CD80", "PDCD1"),
         mstatus==1)%>%
  mutate(targetNamef=factor(targetName, levels=c("CTLA4", "CD274", "CD80", "PDCD1")))%>%
  ggplot(aes(x=factor(n_infection), y=conc, fill=factor(n_infection)))+
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
  filter(targetName %in% c("CTLA4", "CD274", "CD80", "PDCD1"))%>%
  mutate(targetNamef=factor(targetName, levels=c("CTLA4", "CD274", "CD80", "PDCD1")))%>%
  ggplot(aes(x=factor(timepoint_num), y=conc, fill=factor(mstatus)))+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetNamef, scales = "free", nrow=2)+
  scale_fill_manual(values=c("black", "darkred"))+
  xlab("age in weeks")+
  ylab("Concentration in Blood")+
  theme_minimal(base_size = 16)+
  theme(legend.position = "none")

nulisa_to_show2+nulisa_to_show


## change through infection history####
### negative binomial; this is needed for para because variance = 4 and mean =1 ####
purf_52 <- clean_data%>%
  filter(mstatus==0, timepoint=="52 weeks")%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~MASS::glm.nb(total_n_malaria_12~conc+log_qpcr+gender_categorical, data=.))) %>%
  mutate(n_para_model=map(data, ~MASS::glm.nb(total_n_para_12~conc+log_qpcr+gender_categorical, data=.))) %>%
  # mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(total_n_malaria_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  # mutate(n_para_model=map(data, ~glmmTMB::glmmTMB(total_n_para_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  mutate(n_para_model_summary=map(n_para_model, ~summary(.)))%>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[14]))%>%
  mutate(n_para_p=map_dbl(n_para_model_summary, ~coef(.)[14]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"),
         n_para_padj=p.adjust(n_para_p, method="BH"))

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


