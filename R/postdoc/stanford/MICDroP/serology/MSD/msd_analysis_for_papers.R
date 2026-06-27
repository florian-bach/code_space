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
maternal_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/DP+SP study/Databases and preliminary findings/Final database used for analyses/DPSP enrollment analysis database_FINAL.dta")

epi_data <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/micdrop_metadata_for_immunology.csv")
nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")
msd_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/msd_results_v2.csv")
vaccine_coverage <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/long_vaccine_coverage_with_recall.csv")
IU_cutoffs <- readxl::read_excel("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/iu_conversion.xlsx")

diagnoses_of_kids_in_msd <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/diagnoses_of_kids_in_msd.csv")


vaccines = c("Diphtheria",     "Measles" ,         "Pertussis",     "Polio",
             "Rotavirus" ,    "Rubella",       "Tetanus", "Pneumo.1.4.14")
vaccines_iga <- c("Diphtheria_IgA", "Measles_IgA",  "Pertussis_IgA", "Polio_IgA",  "Rotavirus_IgA", "Rubella_IgA",  "Tetanus_IgA",  "Varicella_IgA")
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
  mutate(mom_rx=maternal_treatment_arms$treatmentarm[match(id-10000, maternal_treatment_arms$id)],
         mom_age=maternal_data$AGE[match(id-10000, maternal_data$id)],
         CMV_at_52_weeks=cmv_data$Diagnosis[match(id, cmv_data$id)],
         antigen=ifelse(antigen=="Diptheria", "Diphtheria", antigen),
         vaccine_type=case_match(antigen,
                                 "Pneumo.1.4.14"~"PCV10",
                                 "Diphtheria"~"DTP_Hib",
                                 "Diphtheria_IgA"~"DTP_Hib",
                                 "Pertussis"~"DTP_Hib",
                                 "Pertussis_IgA"~"DTP_Hib",
                                 "Tetanus"~"DTP_Hib",
                                 "Measles"~"Measles Rubella",
                                 "Rubella"~"Measles Rubella",
                                 "Polio"~"Sabin Polio",
                                 "Polio_IgA"~"Oral Polio",
                                 "Rotavirus"~"Rotarix",
                                 "Rotavirus_IgA"~"Rotarix", .default = "no vaccine"
                                 ))%>%
  mutate(mom_rx=case_match(mom_rx, 1~"SP", 2~"DP", 3~"DPSP"))%>%
  left_join(., vaccine_coverage, by=c("id", "vaccine_type"))

  
# additional data from the moment of sampling
additional_clinical_data <- nulisa_data%>%
  distinct(id, sample, date, ageinwks, gender_categorical, mstatus, qPCRparsdens, pardens)

antibodies_and_epi <- epi_data%>%
  inner_join(., long_msd, by="id")%>%
  left_join(additional_clinical_data, by=c("sample", "id"))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))

write.csv(antibodies_and_epi, "~/postdoc/stanford/plasma_analytes/MICDROP/MSD/long_msd_and_epi.csv", row.names = F)
# model selection ####
## gender differences (not significant, no AIC gains) ####
# there are no significant gender differences; gender does not improve prediction
# auxilliary functions
get_lrt_p <- function(anova_tbl) {
  anova_tbl$`Pr(>Chisq)`[2]
}

get_gender_coef <- function(model) {
  coefs <- summary(model)$coefficients
  # assumes reference level is the first level of gender
  coefs[grep("gender_categorical", rownames(coefs)), "Estimate"]
}

get_gender_coef_p <- function(model) {
  coefs <- summary(model)$coefficients
  coefs[grep("gender_categorical", rownames(coefs)), "Pr(>|t|)"]
}


gender_purf <- antibodies_and_epi %>%
  mutate(log_titer=log10(titer))%>%
  filter(!is.na(gender_categorical))%>%
  group_by(antigen) %>%
  nest() %>%
  mutate(
    base_model = map(data, ~lmerTest::lmer(log_titer ~ timepoint + (1|id), data=.)),
    gender_add = map(data, ~lmerTest::lmer(log_titer ~ timepoint + gender_categorical + (1|id), data=.)),
    gender_interact = map(data, ~lmerTest::lmer(log_titer ~ timepoint * gender_categorical + (1|id), data=.))
  )%>%
  mutate(
    AIC_base = map_dbl(base_model, AIC),
    AIC_gender = map_dbl(gender_add, AIC),
    delta_AIC_gender = AIC_base - AIC_gender
  )%>%
  mutate(
    gender_main_test = purrr::map2(base_model, gender_add, anova),
    gender_interaction_test = purrr::map2(gender_add, gender_interact, anova)
  )%>%
  mutate(
    p_gender_main = map_dbl(gender_main_test, get_lrt_p),
    p_gender_interaction = map_dbl(gender_interaction_test, get_lrt_p),
    
    gender_coef = map_dbl(gender_add, get_gender_coef),
    gender_coef_p = map_dbl(gender_add, get_gender_coef_p)
  ) %>%
  # FDR correction across antigens
  ungroup() %>%
  mutate(
    p_gender_main_fdr = p.adjust(p_gender_main, method = "fdr"),
    p_gender_interaction_fdr = p.adjust(p_gender_interaction, method = "fdr"),
    gender_coef_p_fdr = p.adjust(gender_coef_p, method = "fdr")
  )


## wealth (not significant, no AIC gains) ####
# there are no significant wealth differences; wealth does not improve prediction

# auxilliary functions
get_lrt_p <- function(anova_tbl) {
  anova_tbl$`Pr(>Chisq)`[2]
}

get_wealth_coef <- function(model) {
  coefs <- summary(model)$coefficients
  # assumes reference level is the first level of wealth
  coefs[grep("wealthcat", rownames(coefs)), "Estimate"]
}

get_wealth_coef_p <- function(model) {
  coefs <- summary(model)$coefficients
  coefs[grep("wealthcat", rownames(coefs)), "Pr(>|t|)"]
}


wealth_purf <- antibodies_and_epi %>%
  mutate(log_titer=log10(titer))%>%
  filter(!is.na(wealthcat))%>%
  group_by(antigen) %>%
  nest() %>%
  mutate(
    base_model = map(data, ~lmerTest::lmer(log_titer ~ timepoint + (1|id), data=.)),
    wealth_add = map(data, ~lmerTest::lmer(log_titer ~ timepoint + wealthcat + (1|id), data=.)),
    wealth_interact = map(data, ~lmerTest::lmer(log_titer ~ timepoint * wealthcat + (1|id), data=.))
  )%>%
  mutate(
    AIC_base = map_dbl(base_model, AIC),
    AIC_wealth = map_dbl(wealth_add, AIC),
    delta_AIC_wealth = AIC_base - AIC_wealth
  )%>%
  mutate(
    wealth_main_test = purrr::map2(base_model, wealth_add, anova),
    wealth_interaction_test = purrr::map2(wealth_add, wealth_interact, anova)
  )%>%
  mutate(
    p_wealth_main = map_dbl(wealth_main_test, get_lrt_p),
    p_wealth_interaction = map_dbl(wealth_interaction_test, get_lrt_p),
    
    wealth_coef = map_dbl(wealth_add, get_wealth_coef),
    wealth_coef_p = map_dbl(wealth_add, get_wealth_coef_p)
  ) %>%
  # FDR correction across antigens
  ungroup() %>%
  mutate(
    p_wealth_main_fdr = p.adjust(p_wealth_main, method = "fdr"),
    p_wealth_interaction_fdr = p.adjust(p_wealth_interaction, method = "fdr"),
    wealth_coef_p_fdr = p.adjust(wealth_coef_p, method = "fdr")
  )

## educ (not significant, no AIC gains)  ####
# no significant interaction, does not improve prediction 
# auxilliary functions
get_lrt_p <- function(anova_tbl) {
  anova_tbl$`Pr(>Chisq)`[2]
}

get_educ_coef <- function(model) {
  coefs <- summary(model)$coefficients
  # assumes reference level is the first level of educ
  coefs[grep("educ", rownames(coefs)), "Estimate"]
}

get_educ_coef_p <- function(model) {
  coefs <- summary(model)$coefficients
  coefs[grep("educ", rownames(coefs)), "Pr(>|t|)"]
}


educ_purf <- antibodies_and_epi %>%
  mutate(log_titer=log10(titer))%>%
  filter(!is.na(educ))%>%
  group_by(antigen) %>%
  nest() %>%
  mutate(
    base_model = map(data, ~lmerTest::lmer(log_titer ~ timepoint + (1|id), data=.)),
    educ_add = map(data, ~lmerTest::lmer(log_titer ~ timepoint + educ + (1|id), data=.)),
    educ_interact = map(data, ~lmerTest::lmer(log_titer ~ timepoint * educ + (1|id), data=.))
  )%>%
  mutate(
    AIC_base = map_dbl(base_model, AIC),
    AIC_educ = map_dbl(educ_add, AIC),
    delta_AIC_educ = AIC_base - AIC_educ
  )%>%
  mutate(
    educ_main_test = purrr::map2(base_model, educ_add, anova),
    educ_interaction_test = purrr::map2(educ_add, educ_interact, anova)
  )%>%
  mutate(
    p_educ_main = map_dbl(educ_main_test, get_lrt_p),
    p_educ_interaction = map_dbl(educ_interaction_test, get_lrt_p),
    
    educ_coef = map_dbl(educ_add, get_educ_coef),
    educ_coef_p = map_dbl(educ_add, get_educ_coef_p)
  ) %>%
  # FDR correction across antigens
  ungroup() %>%
  mutate(
    p_educ_main_fdr = p.adjust(p_educ_main, method = "fdr"),
    p_educ_interaction_fdr = p.adjust(p_educ_interaction, method = "fdr"),
    educ_coef_p_fdr = p.adjust(educ_coef_p, method = "fdr")
  )

## log_qpcr (only matters for Varicella) ####

# auxilliary functions
get_lrt_p <- function(anova_tbl) {
  anova_tbl$`Pr(>Chisq)`[2]
}

get_log_qpcr_coef <- function(model) {
  coefs <- summary(model)$coefficients
  # assumes reference level is the first level of log_qpcr
  coefs[grep("log_qpcr", rownames(coefs)), "Estimate"]
}

get_log_qpcr_coef_p <- function(model) {
  coefs <- summary(model)$coefficients
  coefs[grep("log_qpcr", rownames(coefs)), "Pr(>|t|)"]
}


log_qpcr_purf <- antibodies_and_epi %>%
  mutate(log_titer=log10(titer))%>%
  filter(timepoint!="8 weeks")%>%
  filter(!is.na(log_qpcr))%>%
  group_by(antigen) %>%
  nest() %>%
  mutate(
    base_model = map(data, ~lmerTest::lmer(log_titer ~ timepoint + (1|id), data=.)),
    log_qpcr_add = map(data, ~lmerTest::lmer(log_titer ~ timepoint + log_qpcr + (1|id), data=.)),
    log_qpcr_interact = map(data, ~lmerTest::lmer(log_titer ~ timepoint * log_qpcr + (1|id), data=.))
  )%>%
  mutate(
    AIC_base = map_dbl(base_model, AIC),
    AIC_log_qpcr = map_dbl(log_qpcr_add, AIC),
    delta_AIC_log_qpcr = AIC_base - AIC_log_qpcr
  )%>%
  mutate(
    log_qpcr_main_test = purrr::map2(base_model, log_qpcr_add, anova),
    log_qpcr_interaction_test = purrr::map2(log_qpcr_add, log_qpcr_interact, anova)
  )%>%
  mutate(
    p_log_qpcr_main = map_dbl(log_qpcr_main_test, get_lrt_p),
    p_log_qpcr_interaction = map_dbl(log_qpcr_interaction_test, get_lrt_p),
    
    log_qpcr_coef = map_dbl(log_qpcr_add, get_log_qpcr_coef),
    log_qpcr_coef_p = map_dbl(log_qpcr_add, get_log_qpcr_coef_p)
  ) %>%
  # FDR correction across antigens
  ungroup() %>%
  mutate(
    p_log_qpcr_main_fdr = p.adjust(p_log_qpcr_main, method = "fdr"),
    p_log_qpcr_interaction_fdr = p.adjust(p_log_qpcr_interaction, method = "fdr"),
    log_qpcr_coef_p_fdr = p.adjust(log_qpcr_coef_p, method = "fdr")
  )

## treatmentarm (only significant for Varicella, borderline significant for Measles IgA & Rubella ) ####

# auxilliary functions
get_lrt_p <- function(anova_tbl) {
  anova_tbl$`Pr(>Chisq)`[2]
}

get_treatmentarm_coef <- function(model) {
  coefs <- summary(model)$coefficients
  # assumes reference level is the first level of treatmentarm
  coefs[grep("treatmentarm", rownames(coefs)), "Estimate"]
}

get_treatmentarm_coef_p <- function(model) {
  coefs <- summary(model)$coefficients
  coefs[grep("treatmentarm", rownames(coefs)), "Pr(>|t|)"]
}


treatmentarm_purf <- antibodies_and_epi %>%
  mutate(log_titer=log10(titer))%>%
  filter(timepoint!="8 weeks")%>%
  filter(!is.na(treatmentarm))%>%
  group_by(antigen) %>%
  nest() %>%
  mutate(
    base_model = map(data, ~lmerTest::lmer(log_titer ~ timepoint + (1|id), data=.)),
    treatmentarm_add = map(data, ~lmerTest::lmer(log_titer ~ timepoint + treatmentarm + (1|id), data=.)),
    treatmentarm_interact = map(data, ~lmerTest::lmer(log_titer ~ timepoint * treatmentarm + (1|id), data=.))
  )%>%
  mutate(
    AIC_base = map_dbl(base_model, AIC),
    AIC_treatmentarm = map_dbl(treatmentarm_add, AIC),
    delta_AIC_treatmentarm = AIC_base - AIC_treatmentarm
  )%>%
  mutate(
    treatmentarm_main_test = purrr::map2(base_model, treatmentarm_add, anova),
    treatmentarm_interaction_test = purrr::map2(treatmentarm_add, treatmentarm_interact, anova)
  )%>%
  mutate(
    p_treatmentarm_main = map_dbl(treatmentarm_main_test, get_lrt_p),
    p_treatmentarm_interaction = map_dbl(treatmentarm_interaction_test, get_lrt_p),
    
    treatmentarm_coef = map_dbl(treatmentarm_add, get_treatmentarm_coef),
    treatmentarm_coef_p = map_dbl(treatmentarm_add, get_treatmentarm_coef_p)
  ) %>%
  # FDR correction across antigens
  ungroup() %>%
  mutate(
    p_treatmentarm_main_fdr = p.adjust(p_treatmentarm_main, method = "fdr"),
    p_treatmentarm_interaction_fdr = p.adjust(p_treatmentarm_interaction, method = "fdr"),
    treatmentarm_coef_p_fdr = p.adjust(treatmentarm_coef_p, method = "fdr")
  )

# infection history ####
## number of parasitemic / malaria episodes ####
# no interaction between total_n_malaria_14w or total_n_para_14w with antibody levels at 8, 24, 52 weeks of age
# with or without adjustment for parasitemia

n_para14w_purf <- antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(id_cat=factor(id))%>%
  filter(timepoint=="24 weeks", total_n_para_14w!=-5)%>%
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
  filter(n_malaria_padj<0.05| n_para_padj<0.05)



# Diphtheria and/or Rotavirus titers are increased at 6 or 12 months, correlated with number of parasitemic / malaria episodes in the first 6 months of life 
# leaving out parasitemia coefficient lowers p value; makes sense because it's a mediator
# n_malaria6 has no influence
n_para6_purf <- antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(id_cat=factor(id))%>%
  filter(timepoint=="24 weeks")%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~lm(log_titer~para_prev_6, data=.))) %>%
  mutate(n_para_model=map(data, ~lm(log_titer~total_n_para_6, data=.))) %>%
  # mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(total_n_malaria_12~log_titer+log_qpcr, data=., family=nbinom2))) %>%
  # mutate(n_para_model=map(data, ~glmmTMB::glmmTMB(total_n_para_12~log_titer+log_qpcr, data=., family=nbinom2))) %>%
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  mutate(n_para_model_summary=map(n_para_model, ~summary(.)))%>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[8]))%>%
  mutate(n_para_p=map_dbl(n_para_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"),
         n_para_padj=p.adjust(n_para_p, method="BH"))

inf_sigs_6 <- n_para6_purf%>%
  filter(n_malaria_padj<0.05| n_para_padj<0.05)

antibodies_and_epi%>%
  filter(timepoint=="52 weeks", antigen %in% c("Rotavirus", "Diphtheria","Rotavirus_IgA", "Diphtheria_IgA"))%>%
  ggplot(., aes(x=factor(total_n_para_6), y=titer, fill=factor(total_n_para_6)))+
  geom_boxplot(outliers = F)+
  geom_point()+
  scale_y_log10()+
  #ggpubr::stat_compare_means(vjust = 0.5, comparisons = list(c("0", "1")))+
  facet_wrap(~antigen, scales="free")+
  scale_fill_manual(values=viridis::viridis(n=7))+
  theme_minimal()

antibodies_and_epi%>%
  filter(timepoint=="52 weeks", antigen %in% c("Rotavirus", "Diphtheria","Rotavirus_IgA", "Diphtheria_IgA"))%>%
  mutate(any_para_6=ifelse(total_n_para_6==0, "none", "some"))%>%
  ggplot(., aes(x=factor(any_para_6), y=titer, fill=factor(any_para_6)))+
  geom_boxplot(outliers = F)+
  geom_point()+
  scale_y_log10()+
  ggpubr::stat_compare_means(vjust = 0.5)+
  facet_wrap(~antigen, scales="free")+
  scale_fill_manual(values=viridis::viridis(n=2))+
  theme_minimal()



# kids with more malaria and more parasitemic months in the first year of life have significantly higher varicella antibodies;
# kids with more parasitemic months in the first year of life have significantly higher Diphtheria antibodies;

n_para12_purf <- antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(id_cat=factor(id))%>%
  filter(timepoint=="52 weeks")%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~lm(log_titer~total_n_malaria_12+ log_qpcr, data=.))) %>%
  mutate(n_para_model=map(data, ~lm(log_titer~total_n_para_12+log_qpcr, data=.))) %>%
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

inf_sigs_12 <- n_para12_purf%>%
  filter(n_malaria_padj<0.1| n_para_padj<0.1)


antibodies_and_epi%>%
  filter(timepoint=="52 weeks", antigen %in% c("Varicella", "Diphtheria", "Rotavirus"))%>%
  ggplot(., aes(x=factor(total_n_malaria_12), y=titer))+
  geom_boxplot(outliers = F)+
  geom_point()+
  scale_y_log10()+
  #ggpubr::stat_compare_means(vjust = 0.5, comparisons = list(c("0", "1")))+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()

## prevalent infections ####
### six months ####
n_preva6_purf <- antibodies_and_epi%>%
  mutate(id_cat=factor(id))%>%
  mutate(log_titer=log10(titer))%>%
  filter(timepoint=="24 weeks")%>%
  group_by(antigen)%>%
  nest()%>%
  # mutate(para_prev_6_bino_model=map(data, ~glm(any_para_6~log_titer+log_qpcr, data=., family = "binomial"))) %>%
  mutate(para_prev_6_model=map(data, ~lm(log_titer~total_n_para_6+log_qpcr+treatmentarm, data=.))) %>%
  # mutate(para_prev_6_bino_model_summary=map(para_prev_6_bino_model, ~summary(.))) %>%
  mutate(para_prev_6_model_summary=map(para_prev_6_model, ~summary(.)))%>%
  # mutate(para_prev_6_bino_p=map_dbl(para_prev_6_bino_model_summary, ~coef(.)11))%>%
  mutate(para_prev_6_p=map_dbl(para_prev_6_model_summary, ~coef(.)[14]))%>%
  ungroup()%>%
  mutate(#para_prev_6_bino_padj=p.adjust(para_prev_6_bino_p, method="BH"),
         para_prev_6_padj=p.adjust(para_prev_6_p, method="BH"))

preva_sigs_6 <- n_preva6_purf%>%
  filter(para_prev_6_padj<0.05)

### 12 months (significant for varicella; also for Diptheria, but only when logqpcr included ####
# when adjusting for qpcr: varicella, diphtheria and diptheria iga
# when adjusting for treatmentarm: only diphtheria IgG
# when adjusting for treatmentarm and qpcr: only diphtheria and diptheria iga

n_preva12_purf <- antibodies_and_epi%>%
  mutate(id_cat=factor(id))%>%
  mutate(log_titer=log10(titer))%>%
  filter(timepoint=="52 weeks")%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(para_prev_12_bino_model=map(data, ~glmmTMB::glmmTMB(total_n_para_12~log_titer+log_qpcr, data=., family=nbinom2))) %>%
  mutate(para_prev_12_model=map(data, ~lm(log_titer~total_n_para_12+log_qpcr, data=.))) %>%
  mutate(para_prev_12_bino_model_summary=map(para_prev_12_bino_model, ~summary(.))) %>%
  mutate(para_prev_12_model_summary=map(para_prev_12_model, ~summary(.)))%>%
  mutate(para_prev_12_bino_p=map_dbl(para_prev_12_bino_model_summary, ~coef(.)$cond[11]))%>%
  mutate(para_prev_12_p=map_dbl(para_prev_12_model_summary, ~coef(.)[11]))%>%
  ungroup()%>%
  mutate(para_prev_12_bino_padj=p.adjust(para_prev_12_bino_p, method="BH"),
         para_prev_12_padj=p.adjust(para_prev_12_p, method="BH"))

preva_sigs_12 <- n_preva12_purf%>%
  filter(para_prev_12_padj<0.1 | para_prev_12_bino_padj<0.1)

## geometric mean of parasitemia ####

# kids with more malaria and/or more parasites in the first year of life have significantly higher varicella antibodies
geom_para12_purf <- antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(id_cat=factor(id))%>%
  filter(timepoint=="52 weeks")%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(geom_para_model=map(data, ~lm(log_titer~log10(geometric_pardens12), data=.))) %>%
  mutate(geom_qpcr_model=map(data, ~lm(log_titer~log10(geometric_qPCRparsdens12), data=.))) %>%
  # mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(total_n_malaria_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  # mutate(geom_para_model=map(data, ~glmmTMB::glmmTMB(total_geom_para_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  mutate(geom_para_model_summary=map(geom_para_model, ~summary(.))) %>%
  mutate(geom_qpcr_model_summary=map(geom_qpcr_model, ~summary(.)))%>%
  #11 when additional covariate is included
  mutate(geom_para_model_p=map_dbl(geom_para_model_summary, ~coef(.)[8]))%>%
  mutate(geom_qpcr_model_p=map_dbl(geom_qpcr_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  mutate(geom_para_model_padj=p.adjust(geom_para_model_p, method="BH"),
         geom_qpcr_model_padj=p.adjust(geom_qpcr_model_p, method="BH"))


geom_sigs_12 <- geom_para12_purf%>%
  filter(geom_para_model_padj<0.1 | geom_qpcr_model_padj<0.1)



geom_para12_spear <- antibodies_and_epi %>%
  mutate(
    log_titer = log10(titer),
    log_para12 = log10(geometric_pardens12 + 0.001),
    log_qpcr12 = log10(geometric_qPCRparsdens12 + 0.001)
  ) %>%
  filter(timepoint == "52 weeks") %>%
  group_by(antigen) %>%
  nest() %>%
  mutate(
    spearman_para = map(data, ~cor.test(
      ~ log_titer + log_para12,
      data = .x,
      method = "spearman",
      exact = FALSE
    )),
    spearman_qpcr = map(data, ~cor.test(
      ~ log_titer + log_qpcr12,
      data = .x,
      method = "spearman",
      exact = FALSE
    )),
    
    rho_para = map_dbl(spearman_para, "estimate"),
    p_para = map_dbl(spearman_para, "p.value"),
    
    rho_qpcr = map_dbl(spearman_qpcr, "estimate"),
    p_qpcr = map_dbl(spearman_qpcr, "p.value")
  ) %>%
  ungroup() %>%
  mutate(
    q_para = p.adjust(p_para, method = "fdr"),
    q_qpcr = p.adjust(p_qpcr, method = "fdr")
  )

geom_spear_sigs_12 <- geom_para12_spear%>%
  filter(q_para<0.1 | q_qpcr<0.1)

antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(id_cat=factor(id))%>%
  # filter(timepoint=="52 weeks", antigen %in% c(vaccines, vaccines_iga[c(3:5)]))%>%
  filter(timepoint=="52 weeks", antigen %in% geom_spear_sigs_12$antigen)%>%
  ggplot(., aes(x=geometric_qPCRparsdens12, y=titer))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_x_log10()+
  scale_y_log10()+
  ggpubr::stat_cor(method = "spearman")+
  facet_wrap(~antigen, scales="free_y")+
  theme_minimal()

antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(id_cat=factor(id))%>%
  filter(timepoint=="52 weeks", antigen%in% viruses)%>%
  ggplot(., aes(x=geometric_pardens12, y=titer))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_x_log10()+
  scale_y_log10()+
  ggpubr::stat_cor(method = "spearman")+
  facet_wrap(~antigen)+
  theme_minimal()

# kids treatment group ####
# increased Measles IgA and Rubella_IgA in kids that receive DP at 1 year
treatment_purf <- antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(id_cat=factor(id))%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(treatment_model=map(data, ~lme4::lmer(log_titer~treatmentarm*timepoint+(1|id_cat), data=.))) %>%
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
  filter(padj<0.05)%>%
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


# maternal gravidity group ####
mom_gravidity_purf <- antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(id_cat=factor(id))%>%
  mutate(gravid_cat = case_match(
    gravidcat,
    1 ~ "Primigravid",
    2 ~ "Multigravid",
    3 ~ "Multigravid"
  ))%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(mom_rx_model=map(data, ~lme4::lmer(log_titer~timepoint*gravid_cat+(1|id_cat), data=.))) %>%
  mutate(mom_rx_model_summary=map(mom_rx_model, ~summary(.)))%>%
  mutate(emm=map(mom_rx_model, ~emmeans(., specs = pairwise ~ gravid_cat | timepoint)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate("8 weeks primi multi"=map_dbl(emm_contrast, ~summary(.)$p.value[1])) %>%
  mutate("24 weeks primi multi"=map_dbl(emm_contrast, ~summary(.)$p.value[2])) %>%
  mutate("52 weeks primi multi"=map_dbl(emm_contrast, ~summary(.)$p.value[3])) %>%
  # mutate("8 weeks DPSP vs SP"=map_dbl(emm_contrast, ~summary(.)$p.value[3])) %>%
  # mutate("24 weeks DPSP vs SP"=map_dbl(emm_contrast, ~summary(.)$p.value[6])) %>%
  # mutate("52 weeks DPSP vs SP"=map_dbl(emm_contrast, ~summary(.)$p.value[9])) %>%
  ungroup()%>%
  pivot_longer(cols=ends_with("multi"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

mom_gravidity_purf_sigs <- mom_gravidity_purf%>%
  filter(padj<0.1)%>%
  select(antigen, contrast, p, padj)


plot_df <- antibodies_and_epi %>%
  mutate(gravid_cat = case_match(
    gravidcat,
    1 ~ "Primigravid",
    2 ~ "Secundigravid",
    3 ~ "Multigravid"
  )) %>%
  filter(antigen %in% mom_gravidity_purf_sigs$antigen)




plot_df%>%
  filter(timepoint=="8 weeks")%>%
  ggplot(., aes(x = timepoint, y = titer, fill=factor(gravid_cat, levels=c("Primigravid", "Secundigravid", "Multigravid")))) +
  geom_boxplot(outliers = FALSE, position = position_dodge(0.75)) +
  scale_y_log10() +
  ylab("IgG concentration (MSD units)")+
  facet_wrap(~antigen, ncol=6, scales="free_y") +
  theme_minimal() +
  scale_fill_manual(values=c("#1E88E5", "#A798EC", "#582EE0"))+
  theme(axis.title.x = element_blank(),
        legend.title =  element_blank())

plot_df%>%
  # filter(timepoint=="8 weeks")%>%
  ggplot(., aes(x = timepoint, fill=factor(gravid_cat, levels=c("Primigravid", "Secundigravid", "Multigravid")), y = titer)) +
  geom_boxplot(outliers = FALSE, position = position_dodge(0.75)) +
  scale_y_log10() +
  ggpubr::stat_compare_means(method="wilcox.test", comparisons = list(c("Multigravid", "Primigravid")))+
  facet_wrap(~antigen, nrow=3) +
  theme_minimal() +
  # scale_fill_manual(values=c("#1E88E5", "#A798EC", "#582EE0"))+
  theme(legend.position = "none")



plot_df%>%
  filter(timepoint=="8 weeks")%>%
  ggplot(., aes(x = GAcomputed, y = titer)) +
  geom_point() +
  scale_y_log10() +
  ggpubr::stat_cor(method="spearman")+
  facet_wrap(~antigen+timepoint, nrow=1) +
  theme_minimal()

plot_df%>%
  filter(timepoint=="8 weeks")%>%
  ggplot(., aes(x = mom_age, y = titer)) +
  geom_point() +
  scale_y_log10() +
  ggpubr::stat_cor(method="spearman")+
  facet_wrap(~antigen+timepoint, nrow=1) +
  theme_minimal()

# maternal histopathology group ####
mom_hp_purf <- antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(id_cat=factor(id))%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(mom_rx_model=map(data, ~lme4::lmer(log_titer~timepoint*factor(anyHP)+(1|id_cat), data=.))) %>%
  mutate(mom_rx_model_summary=map(mom_rx_model, ~summary(.)))%>%
  mutate(emm=map(mom_rx_model, ~emmeans(., specs = pairwise ~ anyHP | timepoint)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate("8 weeks anyhp"=map_dbl(emm_contrast, ~summary(.)$p.value[1])) %>%
  mutate("24 weeks anyhp"=map_dbl(emm_contrast, ~summary(.)$p.value[2])) %>%
  mutate("52 weeks anyhp"=map_dbl(emm_contrast, ~summary(.)$p.value[3])) %>%
  # mutate("8 weeks DPSP vs SP"=map_dbl(emm_contrast, ~summary(.)$p.value[3])) %>%
  # mutate("24 weeks DPSP vs SP"=map_dbl(emm_contrast, ~summary(.)$p.value[6])) %>%
  # mutate("52 weeks DPSP vs SP"=map_dbl(emm_contrast, ~summary(.)$p.value[9])) %>%
  ungroup()%>%
  pivot_longer(cols=ends_with("anyhp"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

mom_hp_sigs <- mom_hp_purf%>%
  filter(padj<0.15)%>%
  select(antigen, contrast, p, padj)


antibodies_and_epi%>%
  filter(timepoint=="8 weeks", antigen %in% mom_gravidity_sigs$antigen)%>%
  ggplot(., aes(x=factor(gravidcat), y=titer))+
  geom_boxplot(outliers = F)+
  geom_point()+
  scale_y_log10()+
  ggpubr::stat_compare_means(vjust = 0.5)+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()


#maternal histopath + gravidity together ####
mom_gravidity_purf2 <- antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(id_cat=factor(id))%>%
  mutate(gravidcat=ifelse(gravid==1, 1, 2))%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(mom_rx_model=map(data, ~lme4::lmer(log_titer~timepoint*factor(gravidcat)+anyHP+(1|id_cat), data=.))) %>%
  mutate(mom_rx_model_summary=map(mom_rx_model, ~summary(.)))%>%
  mutate(emm=map(mom_rx_model, ~emmeans(., specs = pairwise ~ gravidcat | timepoint)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate("8 weeks primi multi"=map_dbl(emm_contrast, ~summary(.)$p.value[1])) %>%
  mutate("24 weeks primi multi"=map_dbl(emm_contrast, ~summary(.)$p.value[2])) %>%
  mutate("52 weeks primi multi"=map_dbl(emm_contrast, ~summary(.)$p.value[3])) %>%
  # mutate("8 weeks DPSP vs SP"=map_dbl(emm_contrast, ~summary(.)$p.value[3])) %>%
  # mutate("24 weeks DPSP vs SP"=map_dbl(emm_contrast, ~summary(.)$p.value[6])) %>%
  # mutate("52 weeks DPSP vs SP"=map_dbl(emm_contrast, ~summary(.)$p.value[9])) %>%
  ungroup()%>%
  pivot_longer(cols=ends_with("multi"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

mom_gravidity_sigs <- mom_gravidity_purf%>%
  filter(padj<0.05)%>%
  select(antigen, contrast, p, padj)


# maternal chemoprevention group ####
mom_rx_purf <- antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(id_cat=factor(id))%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(mom_rx_model=map(data, ~lme4::lmer(log_titer~mom_rx*timepoint+log_qpcr+wealthcat+educ+(1|id_cat), data=.))) %>%
  mutate(mom_rx_model_summary=map(mom_rx_model, ~summary(.)))%>%
  mutate(emm=map(mom_rx_model, ~emmeans(., specs = pairwise ~ mom_rx | timepoint)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate("8 weeks DP vs SP"=map_dbl(emm_contrast, ~summary(.)$p.value[2])) %>%
  mutate("24 weeks DP vs SP"=map_dbl(emm_contrast, ~summary(.)$p.value[5])) %>%
  mutate("52 weeks DP vs SP"=map_dbl(emm_contrast, ~summary(.)$p.value[8])) %>%
  mutate("8 weeks DPSP vs SP"=map_dbl(emm_contrast, ~summary(.)$p.value[3])) %>%
  mutate("24 weeks DPSP vs SP"=map_dbl(emm_contrast, ~summary(.)$p.value[6])) %>%
  mutate("52 weeks DPSP vs SP"=map_dbl(emm_contrast, ~summary(.)$p.value[9])) %>%
  ungroup()%>%
  pivot_longer(cols=ends_with("vs SP"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

mom_rx_sigs <- mom_rx_purf%>%
  filter(padj<0.05)%>%
  select(antigen, contrast, p, padj)

# correlate 8 weeks with 52 weeks ####
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
  mutate(log_pardens=log10(pardens+0.001))%>%
  # mutate(id_cat=factor(id))%>%
  filter(timepoint!="8 weeks")%>%
  group_by(antigen, timepoint)%>%
  nest()%>%
  mutate(parasitemia_corr=map(data, ~cor.test(.$log_titer, .$qPCRparsdens, method="spearman"))) %>%
  mutate(parasitemia_corr_p=map_dbl(parasitemia_corr, ~.$p.value))%>%
  mutate(parasitemia_corr_rho=map_dbl(parasitemia_corr, ~.$estimate))%>%
  #11 when additional covariate is included
  # ungroup()%>%
  group_by(timepoint)%>%
  mutate(parasitemia_padj=p.adjust(parasitemia_corr_p, method="BH"))

parasitemia_sigs <- parasitemia_purf%>%
  filter(parasitemia_padj<0.05)%>%
  arrange(desc(parasitemia_corr_rho))


antibodies_and_epi%>%
  filter(timepoint=="24 weeks", antigen %in% c("Diphtheria"))%>%
  ggplot(., aes(x=qPCRparsdens, y=titer))+
  geom_point()+
  geom_smooth(method="lm")+
  scale_x_log10()+
  scale_y_log10()+
  ggpubr::stat_cor(method="spearman")+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()



# infection history ####
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
  mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(incidence_value~log_titer, data=., family="nbinom1"))) %>%
  # mutate(n_malaria_model_bino=map(data, ~glmmTMB::glmmTMB(incidence_value/12~log_titer+log_qpcr, data=., family=binomial)))%>%
  # mutate(n_malaria_model_aic=map_dbl(n_malaria_model, ~AIC(.))) %>%
  # mutate(n_malaria_model_2_aic=map_dbl(n_malaria_model_2, ~AIC(.)))%>%
  # mutate(AIC_diff=n_malaria_model_aic-n_malaria_model_2_aic)
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  # mutate(n_malaria_model_bino_summary=map(n_malaria_model_bino, ~summary(.))) %>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~.$coefficients$cond[[8]]))%>%
  # mutate(n_malaria_bino_p=map_dbl(n_malaria_model_bino_summary, ~.$coefficients$cond[[11]]))%>%
  ungroup()%>%
  group_by(timepoint, incidence_type)%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"))
# mutate(n_malaria_bino_padj=p.adjust(n_malaria_bino_p, method="BH"))


# 12 with qPCR; 3 without
inf_sigs_14w <-para_purf%>%
  filter(n_malaria_padj<0.05)%>%
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
  filter(timepoint=="52 weeks", antigen %in% c("Diphtheria", "Rotavirus_IgA", "Pertussis"))%>%
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

## infection history model selection ####


para_purf <- long_antibodies_and_epi%>%
  # filter(timepoint!="8 weeks")%>%
  filter(timepoint=="52 weeks"&incidence_type%in%c("total_n_para_14w", "total_n_para_6", "total_n_para_6_12", "total_n_para_12", "total_n_malaria_14w", "total_n_malaria_6", "total_n_malaria_6_12", "total_n_malaria_12"))%>%
  # timepoint=="24 weeks"&incidence_type%in%c("total_n_para_14w", "total_n_para_6", "total_n_malaria_14w", 'total_n_malaria_6'))%>%
  group_by(timepoint, antigen, incidence_type)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~lm(log_titer~incidence_value+log_qpcr+gender_categorical, data=.))) %>%
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[14]))%>%
  mutate(AIC=map_dbl(n_malaria_model, ~AIC(.)))%>%
  ungroup()%>%
  group_by(timepoint, incidence_type)%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"))



# 7 when incidence_value~log_titer
# 4 when incidence_value~log_titer+gender_categorical
# 5 when incidence_value~log_titer+log_qpcr+gender_categorical


para_purf <- long_antibodies_and_epi%>%
  # filter(timepoint!="8 weeks")%>%
  filter(timepoint=="52 weeks"&incidence_type%in%c("total_n_para_14w", "total_n_para_6", "total_n_para_6_12", "total_n_para_12", "total_n_malaria_14w", "total_n_malaria_6", "total_n_malaria_6_12", "total_n_malaria_12"))%>%
  # timepoint=="24 weeks"&incidence_type%in%c("total_n_para_14w", "total_n_para_6", "total_n_malaria_14w", 'total_n_malaria_6'))%>%
  group_by(timepoint, antigen, incidence_type)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~lm(log_titer~treatmentarm, data=.))) %>%
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[8]))%>%
  mutate(AIC=map_dbl(n_malaria_model, ~AIC(.)))%>%
  ungroup()%>%
  group_by(timepoint, incidence_type)%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"))

# 3 when incidence_value~log_titer
# 4 when log_titer~incidence_value+gender_categorical
# 4 when log_titer~incidence_value+log_qpcr+gender_categorical


para_purf <- long_antibodies_and_epi%>%
  # filter(timepoint!="8 weeks")%>%
  distinct(id.x, log_titer, treatmentarm, timepoint, antigen)%>%
  group_by(timepoint, antigen)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~lm(log_titer~treatmentarm, data=.))) %>%
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[8]))%>%
  mutate(AIC=map_dbl(n_malaria_model, ~AIC(.)))%>%
  ungroup()%>%
  group_by(timepoint)%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"))

# 0 when log_titer~treatmentarm

# seropositivity ####

final_cutoff_frame <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/final_cutoff_frame.csv")
colnames(final_cutoff_frame)[2]="cut_off"

wide_long_msd <-  long_msd %>%
  filter(!is.na(timepoint))%>%
  pivot_wider(names_from = "timepoint", values_from = "titer", id_cols = c("id", "antigen", "Ig_class", "treatmentarm"))%>%
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
df_pred <- antibodies_and_epi%>%
  distinct(id, antigen, timepoint, titer, qPCRparsdens) %>%
  filter(timepoint %in% c("24 weeks", "52 weeks"), antigen %in% c(vaccines, vaccines_iga)) %>%
  rename(
    vaccine_antigen = antigen,
    vaccine_titer = titer
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
  mutate(id_cat=factor(as.character(id)))


system.time(
            mod <- lme4::lmer(
              log2_titer ~ virus_seroconversion * virus_antigen * vaccine_antigen + (1 | id),
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
  filter(p.value<0.1)%>%
  nrow()


#polio impacted by PV1; Rotavirus impacted by RV.c; Rubella imapcted by Flu.A.H3, PIV.4, RV.C
virus_conversion_vaccine_plot <- wide_long_msd_cross%>%
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
  filter(padj<0.05)

Tetanus_IgA_analytes <- purf%>%
  filter(padj<0.05, antigen=="Tetanus_IgA")

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
  filter(padj<0.05, targetName=="IL2")

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
  filter(padj<0.05, targetName=="EPO")

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

# vaccine coverage ####

schedule <- data.frame(vaccine_type=c("BCG at birth",
                                      "DTP_Hib",
                                      "Oral Polio at birth",
                                      "Measles Rubella",
                                      "PCV10",
                                      "Rotarix",
                                      "Oral Polio",
                                      "Sabin Polio"),
                       target_doses=c(1, 3, 1, 1, 3, 2, 3, 1))

df <- read.csv("~/Downloads/df_coverage_with_recall.csv")

long_vaccine_coverage2 <- df%>%
  mutate(rota2=ifelse(rota2=="NULL", 0, rota2))%>%
  mutate(rota2=as.numeric(rota2))%>%
  pivot_longer(cols = colnames(df)[3:19], names_to = "vaccine_dose", values_to="dose_received")%>%
  mutate(vaccine_type = case_when(grepl("polio0", vaccine_dose)~"Oral Polio at birth",
                                  grepl("polio*", vaccine_dose)~"Oral Polio",
                                  grepl("Bcg*", vaccine_dose)~"BCG at birth",
                                  grepl("pcv*", vaccine_dose, ignore.case = T)~"PCV10",
                                  grepl("penta*", vaccine_dose)~"DTP_Hib",
                                  grepl("rota*", vaccine_dose)~"Rotarix",
                                  grepl("mr", vaccine_dose)~"Measles Rubella",
                                  grepl("ipb*", vaccine_dose)~"Sabin Polio"))%>%
  group_by(vaccine_type, id)%>%
  summarise("number_of_doses_received"=sum(dose_received, na.rm = T))%>%
  left_join(., schedule)%>%
  mutate(course_outcome=case_when(number_of_doses_received>=target_doses~"complete",
                                  number_of_doses_received<target_doses&number_of_doses_received>0~"partial",
                                  number_of_doses_received==0~"none",
                                  is.na(number_of_doses_received)~"no record"
  ))

# write.csv(long_vaccine_coverage2, "~/postdoc/stanford/plasma_analytes/MICDROP/MSD/long_vaccine_coverage_with_recall.csv", row.names = F)



df2 <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/df_coverage_card_only.csv")

long_vaccine_coverage3 <- df2%>%
  mutate(rota2=ifelse(rota2=="NULL", 0, rota2))%>%
  mutate(rota2=as.numeric(rota2))%>%
  pivot_longer(cols = colnames(df)[3:19], names_to = "vaccine_dose", values_to="dose_received")%>%
  mutate(vaccine_type = case_when(grepl("polio0", vaccine_dose)~"Oral Polio at birth",
                                  grepl("polio*", vaccine_dose)~"Oral Polio",
                                  grepl("Bcg*", vaccine_dose)~"BCG at birth",
                                  grepl("pcv*", vaccine_dose, ignore.case = T)~"PCV10",
                                  grepl("penta*", vaccine_dose)~"DTP_Hib",
                                  grepl("rota*", vaccine_dose)~"Rotarix",
                                  grepl("mr", vaccine_dose)~"Measles Rubella",
                                  grepl("ipb*", vaccine_dose)~"Sabin Polio"))%>%
  group_by(vaccine_type, id)%>%
  summarise("number_of_doses_received"=sum(dose_received, na.rm = T))%>%
  left_join(., schedule)%>%
  mutate(course_outcome=case_when(number_of_doses_received>=target_doses~"complete",
                                  number_of_doses_received<target_doses&number_of_doses_received>0~"partial",
                                  number_of_doses_received==0~"none",
                                  is.na(number_of_doses_received)~"no record"
  ))

write.csv(long_vaccine_coverage3, "~/postdoc/stanford/plasma_analytes/MICDROP/MSD/long_vaccine_coverage_card_only.csv", row.names = F)

# causal mediation ####
run_med_long <- function(df) {
  
  # mediator model (same for every analyte, but that's fine)
  med.fit <- lm(log_qpcr ~ para_prev_12 + treatmentarm, data=df)
  
  # outcome model
  out.fit <- lm(log_titer~para_prev_12+log_qpcr+treatmentarm+wealthcat+educ, data=df)
  
  med.out <- mediate(
    med.fit,
    out.fit,
    treat = "para_prev_12",
    mediator = "log_qpcr",
    sims = 1000,
    boot = FALSE
  )
  
  tibble(
    ACME  = med.out$d.avg,
    ADE   = med.out$z.avg,
    TE    = med.out$tau.coef,
    PM    = med.out$n.avg,
    p_ACME = med.out$d.avg.p,
    p_ADE  = med.out$z.avg.p,
    p_TE   = med.out$tau.p
  )
}

system.time(
  results <- antibodies_and_epi %>%
    mutate(log_titer=log10(titer))%>%
    mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
    filter(timepoint=="52 weeks", mstatus==0, !is.na(wealthcat))%>%
    group_by(antigen) %>%
    nest() %>%
    mutate(med = map(data, run_med_long)) %>%
    unnest(med) %>%
    ungroup()
)


results <- results %>%
  mutate(
    q_ACME = p.adjust(p_ACME, method = "fdr"),
    q_ADE  = p.adjust(p_ADE,  method = "fdr"),
    q_TE   = p.adjust(p_TE,   method = "fdr"),
    prop_mediated = ACME / TE
  )

results <- results %>%
  mutate(
    mediation_class = case_when(
      q_ACME < 0.05 & q_ADE >= 0.05 ~ "Mediated",
      q_ADE  < 0.05 & q_ACME >= 0.05 ~ "Direct",
      q_ACME < 0.05 & q_ADE  < 0.05 ~ "Both",
      q_TE   < 0.05 ~ "Total only",
      TRUE ~ "None"
    )
  )
write.csv(results[,-2], "~/postdoc/stanford/plasma_analytes/MICDROP/MSD/msd_causal_mediation_results_total_n_para_12_with_confounders.csv", row.names = F)
