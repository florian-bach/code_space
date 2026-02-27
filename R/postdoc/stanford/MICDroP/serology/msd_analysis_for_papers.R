library(tidyr)
library(dplyr)
library(xlsx)
library(ggplot2)
library(purrr)
library(data.table)
library(glmmTMB)
library(emmeans)

mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")
maternal_treatment_arms <- haven::read_dta("~/Library/CloudStorage/Box-Box/DP+SP study/Databases and preliminary findings/Final database used for analyses/DPSP treatment allocation_FINAL.dta")
epi_data <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/micdrop_metadata_for_immunology.csv")
nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")
msd_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/msd_results_v2.csv")
vaccines = c("Diptheria",     "Measles" ,         "Pertussis",     "Polio",
             "Rotavirus" ,    "Rubella",       "Tetanus", "Pneumo.1.4.14")
vaccines_iga <- c("Diptheria_IgA", "Measles_IgA",  "Pertussis_IgA", "Polio_IgA",  "Rotavirus_IgA", "Rubella_IgA",  "Tetanus_IgA",  "Varicella_IgA")
viruses <- c("Mumps",    
             "Varicella",
             "CoV.2",    
             "EV.71",    
             "EV.D68",   
             "Flu.Mix",  
             'hMPV',     
             "PIV.Mix",  
             "RSV",      
             "RV.C",     
             "Flu.A.H1", 
             "Flu.A.H3", 
             "Flu.B",    
             "PIV.1",    
             "PIV.2",    
             'PIV.3',    
             "PIV.4")
cmv_data <-  read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/CMV_ELISA_Uganda.csv", skip = 1)


long_msd <- msd_data%>%
  mutate(sample=paste(SubjectID, "_", "tp", TimePt, sep=""))%>%
  mutate(id=SubjectID, timepoint=paste(TimePt, "weeks"))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))%>%
  select(-SubjectID, -TimePt)%>%
  pivot_longer(cols=-c(sample, id, timepoint), names_to = "antigen", values_to = "titer")%>%
  mutate(Ig_class=ifelse(grepl("*IgA$", antigen), "IgA", "IgG"))%>%
  filter(!is.na(titer))%>%
  mutate(CMV_at_52_weeks=cmv_data$Diagnosis[match(id, cmv_data$id)],)


# additional data from the moment of sampling
additional_clinical_data <- nulisa_data%>%
  distinct(id, sample, date, ageinwks, gender_categorical, mstatus, qPCRparsdens)

antibodies_and_epi <- epi_data%>%
  inner_join(., long_msd, by="id")%>%
  left_join(additional_clinical_data, by="sample")

# infection history ####
# no interaction between total_n_malaria_14w or total_n_para_14w with antibody levels at 8, 24, 52 weeks of age
# with or without adjustment for parasitemia

n_para14w_purf <- antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(id_cat=factor(id.x))%>%
  filter(timepoint=="8 weeks", total_n_para_14w!=-5)%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~lm(log_titer~total_n_malaria_14w, data=.))) %>%
  mutate(n_para_model=map(data, ~lm(log_titer~total_n_para_14w, data=.))) %>%
  # mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(total_n_malaria_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  # mutate(n_para_model=map(data, ~glmmTMB::glmmTMB(total_n_para_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  mutate(n_para_model_summary=map(n_para_model, ~summary(.)))%>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[8]))%>%
  mutate(n_para_p=map_dbl(n_para_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"),
         n_para_padj=p.adjust(n_para_p, method="BH"))


inf_sigs_14w <- n_para14w_purf%>%
  filter(n_malaria_padj<0.1| n_para_padj<0.1)



# Diptheria and/or Rotavirus titers are increased at 6 or 12 months, correlated with number of parasitemic / malaria episodes in the first 6 months of life 
# leaving out parasitemia coefficient lowers p value 
n_para6_purf <- antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(id_cat=factor(id.x))%>%
  filter(timepoint=="52 weeks")%>%
  group_by(antigen)%>%
  nest()%>%
  # mutate(n_malaria_model=map(data, ~lm(log_titer~total_n_malaria_6, data=.))) %>%
  # mutate(n_para_model=map(data, ~lm(log_titer~total_n_para_6, data=.))) %>%
  mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(total_n_malaria_12~log_titer+log_qpcr, data=., family=nbinom2))) %>%
  mutate(n_para_model=map(data, ~glmmTMB::glmmTMB(total_n_para_12~log_titer+log_qpcr, data=., family=nbinom2))) %>%
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  mutate(n_para_model_summary=map(n_para_model, ~summary(.)))%>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)$cond[11]))%>%
  mutate(n_para_p=map_dbl(n_para_model_summary, ~coef(.)$cond[11]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"),
         n_para_padj=p.adjust(n_para_p, method="BH"))

inf_sigs_6 <- n_para6_purf%>%
  filter(n_malaria_padj<0.1| n_para_padj<0.1)

antibodies_and_epi%>%
  filter(timepoint=="52 weeks", antigen %in% c("Rotavirus", "Diptheria","Rotavirus_IgA", "Diptheria_IgA"))%>%
  ggplot(., aes(x=factor(total_n_para_6), y=titer, fill=factor(total_n_para_6)))+
  geom_boxplot(outliers = F)+
  geom_point()+
  scale_y_log10()+
  #ggpubr::stat_compare_means(vjust = 0.5, comparisons = list(c("0", "1")))+
  facet_wrap(~antigen, scales="free")+
  scale_fill_manual(values=viridis::viridis(n=7))+
  theme_minimal()

antibodies_and_epi%>%
  filter(timepoint=="52 weeks", antigen %in% c("Rotavirus", "Diptheria","Rotavirus_IgA", "Diptheria_IgA"))%>%
  mutate(any_para_6=ifelse(total_n_para_6==0, "none", "some"))%>%
  ggplot(., aes(x=factor(any_para_6), y=titer, fill=factor(any_para_6)))+
  geom_boxplot(outliers = F)+
  geom_point()+
  scale_y_log10()+
  ggpubr::stat_compare_means(vjust = 0.5)+
  facet_wrap(~antigen, scales="free")+
  scale_fill_manual(values=viridis::viridis(n=2))+
  theme_minimal()



# kids with more malaria and/or more parasites in the first year of life have significantly higher varicella antibodies
n_para12_purf <- antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(id_cat=factor(id.x))%>%
  filter(timepoint=="52 weeks")%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~lm(log_titer~total_n_malaria_12, data=.))) %>%
  mutate(n_para_model=map(data, ~lm(log_titer~total_n_para_12, data=.))) %>%
  # mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(total_n_malaria_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  # mutate(n_para_model=map(data, ~glmmTMB::glmmTMB(total_n_para_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  mutate(n_para_model_summary=map(n_para_model, ~summary(.)))%>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[8]))%>%
  mutate(n_para_p=map_dbl(n_para_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"),
         n_para_padj=p.adjust(n_para_p, method="BH"))

inf_sigs_12 <- n_para12_purf%>%
  filter(n_malaria_padj<0.1| n_para_padj<0.1)


antibodies_and_epi%>%
  filter(timepoint=="52 weeks", antigen %in% c("Varicella", "Diptheria", "Rotavirus"))%>%
  ggplot(., aes(x=factor(total_n_malaria_12), y=titer))+
  geom_boxplot(outliers = F)+
  geom_point()+
  scale_y_log10()+
  #ggpubr::stat_compare_means(vjust = 0.5, comparisons = list(c("0", "1")))+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()


# kids treatment group ####
# increased Measles IgA and Rubella_IgA in kids that receive DP at 1 year
treatment_purf <- antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(id_cat=factor(id.x))%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(treatment_model=map(data, ~lme4::lmer(log_titer~treatmentarm*timepoint+log_qpcr+(1|id_cat), data=.))) %>%
  mutate(treatment_model_summary=map(treatment_model, ~summary(.)))%>%
  mutate(emm=map(treatment_model, ~emmeans(., specs = pairwise ~ treatmentarm | timepoint)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate("8 weeks"=map_dbl(emm_contrast, ~summary(.)$p.value[1])) %>%
  mutate("24 weeks"=map_dbl(emm_contrast, ~summary(.)$p.value[2])) %>%
  mutate("52 weeks"=map_dbl(emm_contrast, ~summary(.)$p.value[3])) %>%
  ungroup()%>%
  pivot_longer(cols=ends_with("weeks"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

treatment_sigs <- treatment_purf%>%
  filter(padj<0.1)%>%
  select(antigen, contrast, p, padj)

#
antibodies_and_epi%>%
  filter(timepoint=="52 weeks", antigen %in% c("Measles", "Rubella", "Measles_IgA", "Rubella_IgA"))%>%
  ggplot(., aes(x=factor(treatmentarm), y=titer))+
  geom_boxplot(outliers = F)+
  geom_point()+
  scale_y_log10()+
  ggpubr::stat_compare_means(vjust = 0.5)+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()

# correlate 8 weeks with 52 weeks
long_msd %>%
  filter(!is.na(timepoint), antigen %in% c(vaccines, vaccines_iga[c(4,5)]))%>%
  pivot_wider(names_from = "timepoint", values_from = "titer", id_cols = c("id", "antigen", "Ig_class"))%>%
  ggplot(., aes(x=`8 weeks`, y=`52 weeks`))+
  geom_point()+
  facet_wrap(~antigen)+
  scale_x_log10()+
  scale_y_log10()+
  ggpubr::stat_cor(method="spearman")+
  geom_smooth(method="lm")+
  theme_minimal()

long_msd %>%
  filter(!is.na(timepoint), antigen %in% c(vaccines, vaccines_iga[c(4,5)]))%>%
  pivot_wider(names_from = "timepoint", values_from = "titer", id_cols = c("id", "antigen", "Ig_class"))%>%
  ggplot(., aes(x=`8 weeks`, y=`24 weeks`))+
  geom_point()+
  facet_wrap(~antigen)+
  scale_x_log10()+
  scale_y_log10()+
  ggpubr::stat_cor(method="spearman")+
  geom_smooth(method="lm")+
  theme_minimal()
# parasitemia ####
# varicella, PIV1, PIV4 antibodies are positively associated with parasitemia
parasitemia_purf <- antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(id_cat=factor(id.x))%>%
  filter(timepoint=="24 weeks")%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(parasitemia_model=map(data, ~lm(log_titer~log_qpcr, data=.))) %>%
  mutate(parasitemia_model_summary=map(parasitemia_model, ~summary(.)))%>%
  #11 when additional covariate is included
  mutate(parasitemia_p=map_dbl(parasitemia_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  mutate(parasitemia_padj=p.adjust(parasitemia_p, method="BH"))

parasitemia_sigs <- parasitemia_purf%>%
  filter(parasitemia_padj<0.1)


antibodies_and_epi%>%
  filter(timepoint=="24 weeks", antigen %in% c("Varicella", "PIV.1", "PIV.4"))%>%
  ggplot(., aes(x=qPCRparsdens, y=titer))+
  geom_point(outliers = F)+
  geom_point()+
  geom_smooth(method="lm")+
  scale_x_log10()+
  scale_y_log10()+
  ggpubr::stat_compare_means(vjust = 0.5)+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()



# parasitemia alternative ####
long_antibodies_and_epi <- antibodies_and_epi%>%
  filter(antigen %in% c(vaccines, vaccines_iga[c(4, 5)]))%>%
  mutate(log_titer=log10(titer))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(id_cat=factor(id.x))%>%
  mutate(total_n_para_6_12=total_n_para_12-total_n_para_6,
         total_n_malaria_6_12=total_n_malaria_12-total_n_malaria_6)%>%
  pivot_longer(cols = c(total_n_para_14w, total_n_para_6, total_n_para_6_12, total_n_para_12,
                        total_n_malaria_14w, total_n_malaria_6, total_n_malaria_6_12, total_n_malaria_12), names_to = "incidence_type", values_to = "incidence_value")%>%
  filter(incidence_value>=0)

# zero inflation is not indicated; comparing AIC of identical models with and without ZI are not significantly different
# nbinom1 usually better
para_purf <- long_antibodies_and_epi%>%
  # filter(timepoint!="8 weeks")%>%
  filter(timepoint=="52 weeks"&incidence_type%in%c("total_n_para_14w", "total_n_para_6", "total_n_para_6_12", "total_n_para_12", "total_n_malaria_14w", "total_n_malaria_6", "total_n_malaria_6_12", "total_n_malaria_12") |
           timepoint=="24 weeks"&incidence_type%in%c("total_n_para_14w", "total_n_para_6", "total_n_malaria_14w", 'total_n_malaria_6'))%>%
  group_by(timepoint, antigen, incidence_type)%>%
  nest()%>%
  # mutate(n_malaria_model=map(data, ~lm(log_titer~total_n_malaria_14w, data=.))) %>%
  # mutate(n_para_model=map(data, ~lm(log_titer~total_n_para_14w, data=.))) %>%
  mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(incidence_value~log_titer+log_qpcr, data=., family=nbinom1))) %>%
  # mutate(n_malaria_model_bino=map(data, ~glmmTMB::glmmTMB(incidence_value/12~log_titer+log_qpcr, data=., family=binomial)))%>%
  # mutate(n_malaria_model_aic=map_dbl(n_malaria_model, ~AIC(.))) %>%
  # mutate(n_malaria_model_2_aic=map_dbl(n_malaria_model_2, ~AIC(.)))%>%
  # mutate(AIC_diff=n_malaria_model_aic-n_malaria_model_2_aic)
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  # mutate(n_malaria_model_bino_summary=map(n_malaria_model_bino, ~summary(.))) %>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~.$coefficients$cond[[11]]))%>%
  # mutate(n_malaria_bino_p=map_dbl(n_malaria_model_bino_summary, ~.$coefficients$cond[[11]]))%>%
  ungroup()%>%
  group_by(timepoint, incidence_type)%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"))
# mutate(n_malaria_bino_padj=p.adjust(n_malaria_bino_p, method="BH"))


# 12 with qPCR; 3 without
inf_sigs_14w <-para_purf%>%
  filter(n_malaria_padj<0.1)%>%
  select(antigen, timepoint, incidence_type,n_malaria_padj, n_malaria_padj)%>%
  arrange(antigen)

plot_list <- vector("list", length = nrow(inf_sigs_14w))

for(i in 1:nrow(inf_sigs_14w)){
  
  ag = inf_sigs_14w$antigen[i]
  tp = inf_sigs_14w$timepoint[i]
  inci_type = inf_sigs_14w$incidence_type[i]
  
  p = long_antibodies_and_epi%>%
    filter(antigen==ag, timepoint==tp, incidence_type==inci_type)%>%
    ggplot(.)+
    geom_point(aes(x=factor(incidence_value), y=titer))+
    # geom_smooth(aes(x=incidence_value, y=titer), method="lm")+
    geom_boxplot(aes(x=factor(incidence_value), y=titer, fill=factor(incidence_value)), outliers = F, inherit.aes = F)+
    scale_y_log10()+
    ggtitle(ag)+
    xlab(inci_type)+
    ylab(paste("titer at", tp))+
    #ggpubr::stat_compare_means(vjust = 0.5, comparisons = list(c("0", "1")))+
    theme_minimal()+
    scale_fill_manual(values = viridis::viridis(n=13))+
    theme(legend.position = "none")
  
  plot_list[[i]] = p
  
}

cowplot::plot_grid(plotlist = plot_list, ncol = 3)


antibodies_and_epi%>%
  filter(total_n_para_12>=0)%>%
  filter(timepoint=="52 weeks", antigen %in% c("Diptheria", "Rotavirus_IgA", "Pertussis"))%>%
  ggplot(., aes(x=factor(total_n_para_6), y=titer, fill=factor(total_n_para_6)))+
  geom_point()+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  theme_minimal()+
  facet_wrap(~antigen, ncol=5)+
  scale_fill_manual(values = viridis::viridis(n=12))+
  theme(legend.position = "none")


# small signal for rotavirus & pneumo
antibodies_and_epi%>%
  filter(timepoint=="52 weeks", antigen %in% c(vaccines, vaccines_iga[c(4, 5)]), CMV_at_52_weeks!="UNKNOWN")%>%
  ggplot(., aes(x=factor(CMV_at_52_weeks), y=titer, fill=factor(CMV_at_52_weeks)))+
  geom_point()+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ggpubr::stat_compare_means()+
  theme_minimal()+
  facet_wrap(~antigen, ncol=5)+
  scale_fill_manual(values = viridis::viridis(n=2))+
  theme(legend.position = "none")
# seropositivity ####

final_cutoff_frame <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/final_cutoff_frame.csv")
colnames(final_cutoff_frame)[2]="cut_off"

wide_long_msd <-  long_msd %>%
  filter(!is.na(timepoint))%>%
  pivot_wider(names_from = "timepoint", values_from = "titer", id_cols = c("id", "antigen", "Ig_class"))%>%
  mutate(log2fc_8_24=log2(`24 weeks`/`8 weeks`),
         log2fc_24_52=log2(`52 weeks`/`24 weeks`))%>%
  group_by(antigen)%>%
  mutate("z_fc_8_24"=scale(log2fc_8_24, center = T))%>%
  pivot_longer(cols = c("8 weeks", "24 weeks","52 weeks"), names_to = "timepoint", values_to = "titer")%>%
  pivot_longer(cols = c("log2fc_8_24", "log2fc_24_52"), names_to = "fc_flavour", values_to = "log2fc")%>%
  mutate(conversion=case_when(timepoint=="52 weeks" & fc_flavour=="log2fc_24_52" & log2fc > 2 ~ "converts 24 to 52",
                              timepoint=="24 weeks" & fc_flavour=="log2fc_8_24" & z_fc_8_24 > 2 | 
                                timepoint=="24 weeks" & fc_flavour=="log2fc_8_24" & log2fc > 1 ~ "converts 8 to 24",
                              .default="nothing"))%>%
  left_join(., final_cutoff_frame, by="antigen")%>%
  mutate(sero_positive=if_else(titer>=cut_off, 1, 0))%>%
  mutate(conversion=factor(conversion, levels=c("nothing", "converts 8 to 24", "converts 24 to 52")))

positivity_df <- wide_long_msd%>%
  distinct(id, timepoint, antigen, sero_positive, Ig_class)%>%
  mutate(mom_rx=maternal_treatment_arms$treatmentarm[match(id-10000, maternal_treatment_arms$id)],
         mom_rx=case_match(mom_rx,
                           1~"SP",
                           2~"DP",
                           3~"DPSP"))


antigens <- unique(positivity_df$antigen)
timepoints <- c("8 weeks", "24 weeks", "52 weeks")



# mad chatGPT model, comparing everything with everything between 6 and 12 months... ####
# should refine to only be about all viruses against all vaccines ####
df_pred <- antibodies_and_epi%>%
  distinct(id.x, antigen, timepoint, titer, qPCRparsdens) %>%
  filter(timepoint == "52 weeks", antigen %in% c(vaccines, vaccines_iga[c(4,5)])) %>%
  rename(
    vaccine_antigen = antigen,
    vaccine_titer = titer,
    id=id.x
  )

df_target <- wide_long_msd%>%
  filter(timepoint == "24 weeks" & fc_flavour=="log2fc_8_24" |
           timepoint == "52 weeks" & fc_flavour=="log2fc_24_52",
         antigen %in% viruses) %>%
  distinct(id, timepoint, antigen, conversion) %>%
  pivot_wider(values_from = conversion, names_from = timepoint)%>%
  mutate(sero_conversion=case_when(`24 weeks`=="nothing"&`52 weeks`=="nothing"~"nothing",
                   `24 weeks`=="nothing"&`52 weeks`=="converts 24 to 52"~"converts 24 to 52",
                   `24 weeks`=="converts 8 to 24"&`52 weeks`=="nothing"~"converts 8 to 24", .default = ))%>%
  distinct(id, antigen, sero_conversion)%>%
  rename(
    virus_antigen = antigen,
    virus_seroconversion = sero_conversion
  )%>%
  filter(!is.na(virus_seroconversion))
  


wide_long_msd_cross <- df_pred %>%
  inner_join(df_target, by = "id") %>%
  mutate(log2_titer = log2(vaccine_titer))%>%
  mutate(log_qPCRparsdens=log10(qPCRparsdens+0.001))

system.time(
            mod <- lme4::lmer(
              log2_titer ~ 
                virus_seroconversion *
                virus_antigen *
                vaccine_antigen +
                log_qPCRparsdens +
                (1 | id),
              data = wide_long_msd_cross
            )
)

system.time(emm <- emmeans(
  mod,
  ~ virus_seroconversion |
    virus_antigen * vaccine_antigen
  # pbkrtest.limit = 15070
))

emm_summary <- summary(
  contrast(emm, "revpairwise"),
  infer = TRUE,
  adjust = "fdr"
)

emm_summary%>%
  filter(p.value<0.05)

#polio impacted by PV1; Rotavirus impacted by RV.c; Rubella imapcted by Flu.A.H3, PIV.4, RV.C
virus_conversion_vaccine_plot <- conversion wide_long_msd_cross%>%
  filter(vaccine_antigen %in% c("Polio", "Rotavirus", "Rotavirus_IgA", "Rubella"),
         virus_antigen %in% c("RV.C", "PIV.4"))%>%
  mutate(virus_seroconversion=factor(virus_seroconversion, levels=c("converts 8 to 24",
                                                                    "converts 24 to 52",
                                                                    "nothing")))%>%
  ggplot(., aes(x=virus_seroconversion, y=log2_titer, fill=virus_seroconversion))+
  geom_violin()+
  geom_boxplot(outliers = F, width = 0.3)+
  theme_minimal()+
  ggpubr::stat_compare_means(vjust=0.2, comparisons = list(
                                                c("converts 8 to 24", "converts 24 to 52"),
                                                
                                                c("nothing", "converts 24 to 52"),
                                                c("nothing", "converts 8 to 24")))+
  facet_wrap(~vaccine_antigen+virus_antigen, nrow=4, scales="free_y")+
  ylab("antibody concentration at 52 weeks")+
  scale_fill_manual(values = c("red", "orange", "black"))+
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")


# seroconversion ####
line_data <- wide_long_msd %>%
  mutate(seroconver)
  distinct(id, timepoint, antigen, titer, conversion, Ig_class)%>%
  group_by(id, antigen, timepoint)%>%
  filter(if (n()>1) conversion!="nothing" else TRUE)%>%
  ungroup()%>%
  mutate(conversion=factor(conversion, levels=c("nothing", "converts 8 to 24", "converts 24 to 52")))%>%
  ungroup()%>%
  arrange(id, factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks"))) %>%
  group_by(id, antigen) %>%
  mutate(next_x = lead(timepoint),
         next_y = lead(titer),
         next_color = lead(conversion))%>%
  filter(!is.na(next_x))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))

## seroconversion statistics ####

conversion_purf <- wide_long_msd_cross%>%
  mutate(contrast=paste(vaccine_antigen, virus_antigen))%>%
  group_by(vaccine_antigen, virus_antigen)%>%
  nest()%>%
  mutate(time_model=purrr::map(data, ~lm(log2_titer~virus_seroconversion+log_qPCRparsdens, data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(contrast_p1=map_dbl(summary, ~coef(.)[14])) %>%
    mutate(contrast_p2=map_dbl(summary, ~coef(.)[15])) %>%
  ungroup()%>%
  group_by(virus_antigen)%>%
  mutate(contrast_padj1 = p.adjust(contrast_p1, method="fdr"))%>%
    mutate(contrast_padj2 = p.adjust(contrast_p2, method="fdr"))

  conversion_purf%>%
    select(vaccine_antigen, virus_antigen, contrast_padj1, contrast_padj2)%>%
    arrange(contrast_padj2)

big_plot <- wide_long_msd%>%
  filter(Ig_class=="IgG")%>%
  arrange(desc(conversion))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=titer))+
  geom_segment(data = filter(line_data, next_color=="nothing", Ig_class=="IgG"), aes(x = timepoint, y = titer, xend = next_x, yend = next_y, color = next_color)) +
  geom_segment(data = filter(line_data, next_color!="nothing", Ig_class=="IgG"), aes(x = timepoint, y = titer, xend = next_x, yend = next_y, color = next_color)) +
  geom_point(aes(color=conversion))+
  scale_y_continuous(trans='log10')+
  scale_color_manual(values = c("black", "red", "orange"))+
  theme_minimal()+
  facet_wrap(~antigen, scales="free")+
  theme(axis.title = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/big_seroconversion_plot.png", big_plot, width=15, height=9, dpi=444, bg="white")


converters_only824 <- wide_long_msd%>%
  group_by(id, antigen)%>%
  filter(any(conversion=="converts 8 to 24"))%>%
  arrange(desc(conversion))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=titer))+
  geom_line(aes(group=id), color="red")+
  scale_y_continuous(trans='log10')+
  theme_minimal()+
  facet_wrap(~antigen, scales="free")+
  theme(axis.title = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/big_seroconverters824_only_plot.png", converters_only824, width=15, height=9, dpi=444, bg="white")


converters_only2452 <- wide_long_msd%>%
  group_by(id, antigen)%>%
  filter(any(conversion=="converts 24 to 52"))%>%
  arrange(desc(conversion))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=titer))+
  geom_line(aes(group=id), color="orange")+
  scale_y_continuous(trans='log10')+
  theme_minimal()+
  facet_wrap(~antigen, scales="free")+
  theme(axis.title = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/big_seroconverters2452_only.png", converters_only2452, width=15, height=9, dpi=444, bg="white")


non_converters <- wide_long_msd%>%
  group_by(id, antigen)%>%
  filter(all(conversion=="nothing"))%>%
  arrange(desc(conversion))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=titer))+
  geom_line(aes(group=id), color="black")+
  scale_y_continuous(trans='log10')+
  theme_minimal()+
  facet_wrap(~antigen, scales="free")+
  theme(axis.title = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/big_serononconverters_only.png", non_converters, width=15, height=9, dpi=444, bg="white")


any_converters <- wide_long_msd%>%
  group_by(id, antigen)%>%
  filter(any(conversion!="nothing"))%>%
  arrange(desc(conversion))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=titer))+
  geom_line(aes(group=id, color=conversion))+
  scale_color_manual(values = c("black", "orange", "red"))+
  scale_y_continuous(trans='log10')+
  theme_minimal()+
  facet_wrap(~antigen, scales="free")+
  theme(axis.title = element_blank(),
        legend.position = "none")

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/big_any_serononconverters_only.png", any_converters, width=15, height=9, dpi=444, bg="white")

## correlate seroconversion with NULISA ####

# need to introduce fold change cutoff !!####

# kids who serconvert from 8 to 24 weeks to Teatnus_IgA
# have less CALCA, CTLA4, IL10RB, IL2RA, IL2RB, LILRB2, PDCD1, NFRSF1, TNFRSF8, TNFRSF1A, VCAM1

# kids who serconvert from 8 to 24 weeks to Measles, Measles IGa, Rotavirus and Rubella IgA
# have less IL2

# kids who serconvert from 24 to 52 weeks to PIV.Mix, PIV.3
# have less EPO
purf <- nulisa_and_seroconversion%>%
  mutate(bino_convert=if_else(conversion=="no", as.numeric(0), as.numeric(1)))%>%
  group_by(targetName, antigen, timepoint)%>%
  nest()%>%
  mutate(model = purrr::map(data, ~glm(bino_convert~conc,  family = "binomial",  data=.)))%>%
  mutate(model_summary=purrr::map(model, ~summary(.)))%>%
  mutate(p=purrr::map(model_summary, ~coef(.)[8]))%>%
  group_by(timepoint, targetName)%>%
  mutate(padj=p.adjust(p, method="fdr"))

sig_purf_conversion_nulisa <- purf%>%
  filter(padj<0.1)

Tetanus_IgA_analytes <- purf%>%
  filter(padj<0.1, antigen=="Tetanus_IgA")

tetanus_iga_hits <- nulisa_and_seroconversion%>%
  semi_join(., sig_purf_conversion_nulisa, by=c("targetName", "antigen"))%>%
  filter(targetName %in% Tetanus_IgA_analytes$targetName)%>%
  filter(antigen=="Tetanus_IgA")%>%
  ggplot(., aes(x=timepoint, y=conc, fill=conversion))+
  geom_boxplot(outliers = F)+
  ggpubr::stat_compare_means(label="p.format", size=2, vjust=.3)+
  facet_wrap(~targetName+antigen, scales="free", nrow=2)+
  scale_fill_manual(values = c("orange", "red", "black"))+
  theme_minimal()+
  theme(legend.position = "bottom")

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/seroconverter_nulisa_sigs.png", seroconverter_hits, width=8, height=4, bg="white", dpi=444)


il2_antigens <- purf%>%
  filter(padj<0.1, targetName=="IL2")

il2_hits <- nulisa_and_seroconversion%>%
  semi_join(., sig_purf_conversion_nulisa, by=c("targetName", "antigen"))%>%
  filter(antigen %in% il2_antigens$antigen, targetName=="IL2")%>%
  ggplot(., aes(x=timepoint, y=conc, fill=conversion))+
  geom_boxplot(outliers = F)+
  ggpubr::stat_compare_means(label="p.format", size=2, vjust=.3)+
  facet_wrap(~targetName+antigen, scales="free", nrow=2)+
  scale_fill_manual(values = c("orange", "red", "black"))+
  theme_minimal()+
  theme(legend.position = "bottom")


EPO_antigens <- purf%>%
  filter(padj<0.1, targetName=="EPO")

EPO_hits <- nulisa_and_seroconversion%>%
  semi_join(., sig_purf_conversion_nulisa, by=c("targetName", "antigen"))%>%
  filter(antigen %in% EPO_antigens$antigen, targetName=="EPO")%>%
  ggplot(., aes(x=timepoint, y=conc, fill=conversion))+
  geom_boxplot(outliers = F)+
  ggpubr::stat_compare_means(label="p.format", size=2, vjust=.3)+
  facet_wrap(~targetName+antigen, scales="free", nrow=2)+
  scale_fill_manual(values = c("orange", "red", "black"))+
  theme_minimal()+
  theme(legend.position = "bottom")



# WGCNA ####

traitRows <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/wgcna_MEs.csv")
long_traitRows <- traitRows%>%
  pivot_longer(cols = starts_with("ME"), names_to = "network", values_to = "network_value")

wgcna_purf <- long_msd%>%
  left_join(., long_traitRows, by="sample")%>%
  mutate(log_titer=log(titer+0.001))%>%
  filter(!is.na(network))%>%
  group_by(antigen, network, timepoint)%>%
  nest()%>%
  mutate(model = purrr::map(data, ~lm(log_titer~network_value, data=.)))%>%
  mutate(model_summary=purrr::map(model, ~summary(.)))%>%
  mutate(p=purrr::map(model_summary, ~coef(.)[8]))%>%
  group_by(antigen, network, timepoint)%>%
  mutate(padj=p.adjust(p, method="fdr"))

wgcna_purf_sigs <- wgcna_purf%>%
  filter(padj<0.05)
