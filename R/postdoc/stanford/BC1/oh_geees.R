# testing interactions when excluding people with 0 malaria episodes as we can't guarantee they weren't exposed
#nothing when using glm.nb, poisson, flipping formula and using lm



six_purf <- six_months %>%
  group_by(antigen) %>%
  nest() %>%
  mutate(model_6_12=map(data, ~glm.nb(inf_6_12 ~ conc, data=.)))%>%
  mutate(model_6_12_symp=map(data, ~glm.nb(symp_6_12 ~ conc, data=.)))%>%
  mutate(model_6_12_symp_prob=map(data, ~glm(symp_6_12/inf_6_12~conc, data=., family = "binomial", weights = inf_6_12)))%>%
  mutate(model_0_6=map(data, ~glm.nb(inf_0_6 ~ conc, data=.)))%>%
  mutate(model_any_symp_6_12=map(data, ~glm(any_symp_6_12~conc, data=., family = "binomial")))%>%
  mutate(model_any_6_12=map(data, ~glm(any_inf_6_12~conc, data=., family = "binomial")))%>%
  
  
  mutate(summary_6_12=map(model_6_12, ~summary(.))) %>%
  mutate(summary_6_12_symp=map(model_6_12_symp, ~summary(.))) %>%
  mutate(summary_6_12_symp_prob=map(model_6_12_symp_prob, ~summary(.))) %>%
  mutate(summary_0_6=map(model_0_6, ~summary(.))) %>%
  mutate(summary_any_symp_6_12=map(model_any_symp_6_12, ~summary(.))) %>%
  mutate(summary_any_6_12=map(model_any_6_12, ~summary(.))) %>%
  
  mutate(summary_6_12_p=map_dbl(summary_6_12, ~unlist(.$coefficients[8])))%>%
  mutate(summary_6_12_symp_p=map_dbl(summary_6_12_symp, ~unlist(.$coefficients[8])))%>%
  mutate(summary_6_12_symp_prob_p=map_dbl(summary_6_12_symp_prob, ~unlist(.$coefficients[8])))%>%
  mutate(summary_0_6_p=map_dbl(summary_0_6, ~unlist(.$coefficients[8])))%>%
  mutate(summary_any_symp_6_12_p=map_dbl(summary_any_symp_6_12, ~unlist(.$coefficients[8])))%>%
  mutate(summary_any_6_12_p=map_dbl(summary_any_6_12, ~unlist(.$coefficients[8])))%>%
  
  ungroup()%>%
  mutate(summary_6_12_padj=p.adjust(summary_6_12_p))%>%
  mutate(summary_6_12_symp_padj=p.adjust(summary_6_12_symp_p))%>%
  mutate(summary_6_12_symp_prob_padj=p.adjust(summary_6_12_symp_prob_p))%>%   
  mutate(summary_0_6_padj=p.adjust(summary_0_6_p))%>%
  mutate(summary_any_symp_6_12_padj=p.adjust(summary_any_symp_6_12_p))%>%
  mutate(summary_any_6_12_padj=p.adjust(summary_any_6_12_p))





#nothing with glm.nb, poissonm, flipped formula lm

twelve_purf <- twelve_months %>%
  filter(symp_12_18>0)%>%
  group_by(antigen) %>%
  nest() %>%
  mutate(model_12_18=map(data, ~lm( conc ~inf_12_18, data=.)))%>%
  mutate(model_12_18_symp=map(data, ~lm(conc ~symp_12_18, data=.)))%>%
  mutate(model_12_18_symp_prob=map(data, ~glm(symp_12_18/inf_12_18~conc, data=., family = "binomial", weights = inf_12_18)))%>%
  mutate(model_6_12=map(data, ~lm( conc ~inf_6_12, data=.)))%>%
  mutate(model_0_12=map(data, ~lm( conc ~inf_0_12, data=.)))%>%
  mutate(model_any_symp_12_18=map(data, ~glm(any_symp_12_18~conc, data=., family = "binomial")))%>%
  mutate(model_any_12_18=map(data, ~glm(any_inf_12_18~conc, data=., family = "binomial")))%>%
  
  mutate(summary_12_18=map(model_12_18, ~summary(.))) %>%
  mutate(summary_12_18_symp=map(model_12_18_symp, ~summary(.))) %>%
  mutate(summary_12_18_symp_prob=map(model_12_18_symp_prob, ~summary(.))) %>%
  mutate(summary_6_12=map(model_6_12, ~summary(.))) %>%
  mutate(summary_0_12=map(model_0_12, ~summary(.))) %>%
  mutate(summary_any_symp_12_18=map(model_any_symp_12_18, ~summary(.))) %>%
  mutate(summary_any_12_18=map(model_any_12_18, ~summary(.))) %>%
  
  mutate(summary_12_18_p=map_dbl(summary_12_18, ~unlist(.$coefficients[8])))%>%
  mutate(summary_12_18_symp_p=map_dbl(summary_12_18_symp, ~unlist(.$coefficients[8])))%>%
  mutate(summary_12_18_symp_prob_p=map_dbl(summary_12_18_symp_prob, ~unlist(.$coefficients[8])))%>%
  mutate(summary_6_12_p=map_dbl(summary_6_12, ~unlist(.$coefficients[8])))%>%
  mutate(summary_0_12_p=map_dbl(summary_0_12, ~unlist(.$coefficients[8])))%>%
  mutate(summary_any_symp_12_18_p=map_dbl(summary_any_symp_12_18, ~unlist(.$coefficients[8])))%>%
  mutate(summary_any_12_18_p=map_dbl(summary_any_12_18, ~unlist(.$coefficients[8])))%>%
  
  mutate(model_12_24=map(data, ~lm( conc ~inf_12_24, data=.)))%>%
  mutate(model_12_24_symp=map(data, ~lm( conc ~symp_12_24, data=.)))%>%
  mutate(model_12_24_symp_prob=map(data, ~glm(symp_12_24/inf_12_24~conc, data=., family = "binomial", weights = inf_12_24)))%>%
  
  mutate(summary_12_24=map(model_12_24, ~summary(.))) %>%
  mutate(summary_12_24_symp=map(model_12_24_symp, ~summary(.))) %>%
  mutate(summary_12_24_symp_prob=map(model_12_24_symp_prob, ~summary(.))) %>%
  
  mutate(summary_12_24_p=map_dbl(summary_12_24, ~unlist(.$coefficients[8])))%>%
  mutate(summary_12_24_symp_p=map_dbl(summary_12_24_symp, ~unlist(.$coefficients[8])))%>%
  mutate(summary_12_24_symp_prob_p=map_dbl(summary_12_24_symp_prob, ~unlist(.$coefficients[8])))%>%
  ungroup()%>%
  mutate(summary_12_18_padj= p.adjust(summary_12_18_p))%>%
  mutate(summary_12_18_symp_padj=p.adjust(summary_12_18_symp_p))%>%
  mutate(summary_12_18_symp_prob_padj=p.adjust(summary_12_18_symp_prob_p))%>%
  mutate(summary_6_12_padj=p.adjust(summary_6_12_p))%>%
  mutate(summary_12_24_padj=p.adjust(summary_12_24_p))%>%
  mutate(summary_12_24_symp_padj=p.adjust(summary_12_24_symp_p))%>%
  mutate(summary_12_24_symp_prob_padj=p.adjust(summary_12_24_symp_prob_p))%>%
  mutate(summary_0_12_padj=p.adjust(summary_0_12_p))%>%
  mutate(summary_any_symp_12_18_padj=p.adjust(summary_any_symp_12_18_p))%>%
  mutate(summary_any_12_18_padj=p.adjust(summary_any_12_18_p))


twelve_18_sig <- twelve_purf %>%
  filter(summary_12_18_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_12_18_padj)

twelve_18_symp_sig <- twelve_purf %>%
  filter(summary_12_18_symp_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_12_18_symp_padj)

twelve_18_symp_prob_sig <- twelve_purf %>%
  filter(summary_12_18_symp_prob_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_12_18_symp_prob_padj)

twelve_6_12_sig <- twelve_purf %>%
  filter(summary_6_12_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_6_12_padj)


twelve_24_sig <- twelve_purf %>%
  filter(summary_12_24_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_12_24_padj)

twelve_24_symp_sig <- twelve_purf %>%
  filter(summary_12_24_symp_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_12_24_symp_padj)

twelve_24_symp_prob_sig <- twelve_purf %>%
  filter(summary_12_24_symp_prob_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_12_24_symp_prob_padj)


twelve_0_12_sig <- twelve_purf %>%
  filter(summary_0_12_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_0_12_padj)

twelve_18_any_sig <- twelve_purf %>%
  filter(summary_any_12_18_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_any_12_18_padj)

twelve_18_any_symp_sig <- twelve_purf %>%
  filter(summary_any_symp_12_18_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_any_symp_12_18_padj)


inf_12_18_data <- filter(twelve_months, inf_12_18>0)

gee_purf <- twelve_months %>%
  filter(inf_12_18>0)%>%
  group_by(antigen) %>%
  nest() %>%
  mutate(inf_gee=map(data, ~geepack::geeglm(inf_12_18 ~ conc,
                               data = ., 
                               id = id, 
                               family = "poisson",
                               corstr = "independence")))%>%
  mutate(inf_gee_summary = map(inf_gee, ~summary(.)))%>%
  mutate(inf_gee_p = map(inf_gee_summary, ~.$coefficients[2,4]))%>%
  mutate(inf_gee_padj = p.adjust(inf_gee_p))

gee_twelve_18_sig <- gee_purf %>%
  filter(inf_gee_padj<fdr_cutoff)%>%
  dplyr::select(antigen, inf_gee_padj)




gee_symp_purf <- twelve_months %>%
  filter(symp_12_18>0)%>%
  group_by(antigen) %>%
  nest() %>%
  mutate(inf_gee_symp=map(data, ~geepack::geeglm(symp_12_18 ~ conc,
                                            data = ., 
                                            id = id, 
                                            family = "binomial",
                                            corstr = "independence")))%>%
  mutate(inf_gee_symp_summary = map(inf_gee_symp, ~summary(.)))%>%
  mutate(inf_gee_symp_p = map(inf_gee_symp_summary, ~.$coefficients[2,4]))%>%
  mutate(inf_gee_symp_padj = p.adjust(inf_gee_symp_p))

gee_symp_twelve_18_sig <- gee_symp_purf %>%
  filter(inf_gee_symp_padj<fdr_cutoff)%>%
  dplyr::select(antigen, inf_gee_symp_padj)



gee_symp_purf <- twelve_months %>%
  filter(symp_12_18>0)%>%
  group_by(antigen) %>%
  nest() %>%
  mutate(inf_gee_symp=map(data, ~geepack::geeglm(symp_12_18/inf_12_18 ~ conc,
                                                 data = ., 
                                                 weights = inf_12_18,
                                                 id = id, 
                                                 family = "binomial",
                                                 corstr = "independence")))%>%
  mutate(inf_gee_symp_summary = map(inf_gee_symp, ~summary(.)))%>%
  mutate(inf_gee_symp_p = map(inf_gee_symp_summary, ~.$coefficients[2,4]))%>%
  mutate(inf_gee_symp_padj = p.adjust(inf_gee_symp_p))

gee_symp_twelve_18_sig <- gee_symp_purf %>%
  filter(inf_gee_symp_padj<0.2)%>%
  dplyr::select(antigen, inf_gee_symp_padj)


twelve_months %>%
  filter(symp_12_18>0)%>%
  filter(antigen %in% c("logHSP40", "logH103", "logGEXP"))%>%
  ggplot(aes(x=inf_12_18, y=conc, fill=factor(inf_12_18)))+
  geom_boxplot()+
  # geom_boxplot()+
  # geom_point(alpha=0.2)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
  xlab("Number of Infections in Months 12-18")+
  ylab("Concentration at 12 Months")+
  ggpubr::stat_cor(method = "pearson", label.x = 1.5, label.y = 1)+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)



twelve_months %>%
  # filter(antigen %in% twelve_18_sig$antigen)%>%
  filter(antigen %in% c("logHSP40", "logH103", "logGEXP"))%>%
  ggplot(aes(x=symp_12_18, y=conc, fill=factor(symp_12_18)))+
  geom_boxplot()+
  # geom_boxplot()+
  # geom_point(alpha=0.2)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
  xlab("Number of Infections in Months 12-18")+
  ylab("Concentration at 12 Months")+
  ggpubr::stat_cor(method = "pearson", label.x = 1.5, label.y = 1)+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)
