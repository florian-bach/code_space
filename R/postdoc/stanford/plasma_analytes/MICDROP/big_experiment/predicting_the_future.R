clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks", "68 weeks")))%>%
  filter(targetName %notin% c("CTSS", "LTA|LTB", "IFNA2"))



future_purf_52 <- clean_data%>%
  filter(mstatus==0, timepoint=="52 weeks")%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~MASS::glm.nb(total_n_malaria_12_24~conc+gender_categorical+log_qpcr, data=.))) %>%
  mutate(n_para_model=map(data, ~MASS::glm.nb(total_n_para_12_24~conc+gender_categorical+log_qpcr, data=.))) %>%
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

# negative binomial means coefficient is on natural log scale.
#so coef of -0.3 -> each 1-unit increase in conc is associated with ~26.6% decrease in the expected number of infections
# higher conc at month 12 is associated with fewer cumulative infections over the first year.


sig_mala_12 <- future_purf_52%>%
  filter(n_malaria_padj<0.1)

sig_para_12 <- future_purf_52%>%
  filter(n_para_padj<0.1)



clean_data %>%
  filter(targetName %in% unique(sig_para_12$targetName),
         timepoint=="52 weeks",
         !is.na(total_n_para_12_24),
         log_qpcr<0.001
         )%>%
  mutate(pardens_dich = if_else(log_qpcr>0.001, "some", "none"))%>%
  mutate(targetNamef=factor(targetName, levels=sig_mala_12$targetName))%>%
  ggplot(aes(x=factor(total_n_para_12_24), y=conc, fill=factor(treatmentarm)))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot()+
  facet_wrap(~targetNamef+treatmentarm, scales = "free")+
  viridis::scale_fill_viridis(discrete = T)+
  xlab("number of parasitemic episodes in months 12-24")+
  facet_wrap(~targetName)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_blank())





# 12_18 only ####
## n cases ####


future_purf_52 <- clean_data%>%
  filter(mstatus==0, timepoint=="52 weeks")%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~MASS::glm.nb(total_n_malaria_12_18~conc+gender_categorical+log_qpcr, data=.))) %>%
  mutate(n_para_model=map(data, ~MASS::glm.nb(total_n_para_12_18~conc+gender_categorical+log_qpcr, data=.))) %>%
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

# negative binomial means coefficient is on natural log scale.
#so coef of -0.3 -> each 1-unit increase in conc is associated with ~26.6% decrease in the expected number of infections
# higher conc at month 12 is associated with fewer cumulative infections over the first year.


sig_mala_12 <- future_purf_52%>%
  filter(n_malaria_padj<0.1)

sig_para_12 <- future_purf_52%>%
  filter(n_para_padj<0.1)



clean_data %>%
  filter(targetName %in% unique(sig_para_12$targetName),
         timepoint=="52 weeks",
         !is.na(total_n_para_12_18),
         log_qpcr<0.001
  )%>%
  mutate(pardens_dich = if_else(log_qpcr>0.001, "some", "none"))%>%
  mutate(targetNamef=factor(targetName, levels=sig_mala_12$targetName))%>%
  ggplot(aes(x=factor(total_n_para_12_18), y=conc, fill=factor(total_n_para_12_18)))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_violin(outliers = FALSE)+
  facet_wrap(~targetNamef, scales = "free")+
  viridis::scale_fill_viridis(discrete = T)+
  xlab("number of parasitemic episodes in months 12-18")+
  ylab("concentration at 12 months")+
  facet_wrap(~targetName)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_blank())


## symp prob####
symp_prob_purf_52 <- clean_data%>%
  filter(mstatus==0, timepoint=="52 weeks")%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~glm(total_n_malaria_12_18/total_n_para_12_18~conc+gender_categorical+log_qpcr, data=., family = "binomial",  weights = total_n_para_12_18))) %>%
  # mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(total_n_malaria_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  # mutate(n_para_model=map(data, ~glmmTMB::glmmTMB(total_n_para_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[14]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"))


symp_mala_12_18 <- symp_prob_purf_52%>%
  filter(n_malaria_p<0.01)



clean_data %>%
  filter(targetName %in% symp_mala_12_24$targetName,
         timepoint=="52 weeks",
         !is.na(total_n_para_12_24),
         log_qpcr<0.001)%>%
  mutate(pardens_dich = if_else(log_qpcr>0.001, "some", "none"))%>%
  ggplot(aes(y=total_n_malaria_12_24/total_n_para_12_24, x=conc))+
  geom_point(aes(color=factor(total_n_malaria_12_24)))+
  facet_wrap(~targetName, scales = "free")+
  geom_smooth(method="gam")+
  ggpubr::stat_cor(method = "spearman", aes(group=1))+
  viridis::scale_color_viridis(discrete = T)+
  xlab("conc")+
  ylab("probability of symptoms, fiven infection")+
  facet_wrap(~targetName, scales="free")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_blank())





symp_prob_purf_52 <- clean_data%>%
  filter(mstatus==0, timepoint=="52 weeks")%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~glm(total_n_malaria_12_/total_n_para_12_24~conc+gender_categorical+log_qpcr, data=., family = "binomial",  weights = total_n_para_12_24))) %>%
  # mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(total_n_malaria_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  # mutate(n_para_model=map(data, ~glmmTMB::glmmTMB(total_n_para_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[14]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"))


symp_mala_12_24 <- symp_prob_purf_52%>%
  filter(n_malaria_padj<0.1)



clean_data %>%
  filter(targetName %in% symp_mala_12_24$targetName,
         timepoint=="52 weeks",
         !is.na(total_n_para_12_24),
         log_qpcr<0.001)%>%
  mutate(pardens_dich = if_else(log_qpcr>0.001, "some", "none"))%>%
  ggplot(aes(y=total_n_malaria_12_24/total_n_para_12_24, x=conc))+
  geom_point(aes(color=factor(total_n_malaria_12_24)))+
  facet_wrap(~targetName, scales = "free")+
  geom_smooth(method="gam")+
  ggpubr::stat_cor(method = "spearman", aes(group=1))+
  viridis::scale_color_viridis(discrete = T)+
  xlab("conc")+
  ylab("probability of symptoms, fiven infection")+
  facet_wrap(~targetName, scales="free")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_blank())




# sandbox LOOK! ####

future_purf_52 <- clean_data%>%
  filter(mstatus==0, timepoint=="52 weeks")%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~MASS::glm.nb(total_n_malaria_12_24~conc+gender_categorical+log_qpcr+total_n_para_12, data=.))) %>%
  mutate(n_para_model=map(data, ~MASS::glm.nb(total_n_para_12_24~conc+gender_categorical+log_qpcr+total_n_para_12, data=.))) %>%
  # mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(total_n_malaria_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  # mutate(n_para_model=map(data, ~glmmTMB::glmmTMB(total_n_para_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  mutate(n_para_model_summary=map(n_para_model, ~summary(.)))%>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[17]))%>%
  mutate(n_para_p=map_dbl(n_para_model_summary, ~coef(.)[17]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"),
         n_para_padj=p.adjust(n_para_p, method="BH"))

# negative binomial means coefficient is on natural log scale.
#so coef of -0.3 -> each 1-unit increase in conc is associated with ~26.6% decrease in the expected number of infections
# higher conc at month 12 is associated with fewer cumulative infections over the first year.


sig_mala_12 <- future_purf_52%>%
  filter(n_malaria_padj<0.1)

sig_para_12 <- future_purf_52%>%
  filter(n_para_padj<0.1)


