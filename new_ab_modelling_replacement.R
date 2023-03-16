# preamble & read in data ####

#palettes
time_palette <- colorspace::sequential_hcl(n=4, "RdPu")[1:3]
pc1_cols <- colorspace::sequential_hcl(23, palette = "Purple Yellow")
incidence_cols <- colorspace::sequential_hcl(11, palette = "Purple Yellow")

#libraries 
library(patchwork)
library(tidyr)
library(ggplot2)
library(purrr)
library(MASS)
library(knitr)
library(dplyr)

fdr_cutoff <- 0.1

`%notin%` <- Negate(`%in%`)

# from looking at the frequencies of observations only these antibodies can be modelled longitudinally; NB some drop off at timepoint 3 so only 1 and 2 can be compared properly
modelable_antigens <- c("Tet Tox", "SBP1", "Rh5", "PfSEA", "PfAMA1", "Hyp2", "HSP40 Ag1", "GST", "GEXP", "CSP GENOVA")


# read in antibody data
bc1 <- haven::read_dta("~/postdoc/stanford/clinical_data/BC1/MergedAntibodyData_ChildClinical.dta")

ab_columns <- grep("log", colnames(bc1), value = TRUE)



long_raw_dfff <- bc1 %>%
  dplyr::select(all_of(c("id", "timepoint", ab_columns)))%>%
  pivot_longer(cols=all_of(ab_columns), names_to = "antigen", values_to = "conc")%>%
  filter(antigen %notin% c("logpd", "logGST"))


# read in new visits database to look at correlations with malaria incidence
#clin_data <- haven::read_dta("~/postdoc/stanford/clinical_data/BC1/BC-1 childs routine visit database FINAL_ALL.dta")
clin_data <- haven::read_dta("~/postdoc/stanford/clinical_data/BC1/BC-1 childs routine visit database FINAL_REV.dta")
# it's a big table, so let's only include kids for whom we have any antibody measurements
clin_data <- clin_data %>%
  filter(id %in% long_raw_dfff$id)

# make a data frame where we add a bunch of columns that contain how many (symptomatic) infections were experienced by the child in the indicated time window

infs <- clin_data %>%
  group_by(id) %>%
  dplyr::select(id, age, anyinfection, sxifinfected) %>%
  mutate(inf_0_6   = sum(if_else(age<6, anyinfection, 0), na.rm = TRUE),
         inf_6_12  = sum(if_else(age>6  & age<12, anyinfection, 0), na.rm = TRUE),
         inf_12_18 = sum(if_else(age>12 & age<18, anyinfection, 0), na.rm = TRUE),
         inf_12_24 = sum(if_else(age>12 & age<24, anyinfection, 0), na.rm = TRUE),
         inf_0_12  = sum(if_else(age<12, anyinfection, 0), na.rm = TRUE),
         inf_0_24  = sum(if_else(age<24, anyinfection, 0), na.rm = TRUE),
         symp_0_6   = sum(if_else(age<6, sxifinfected, 0), na.rm = TRUE),
         symp_6_12  = sum(if_else(age>6  & age<12, sxifinfected, 0), na.rm = TRUE),
         symp_12_18 = sum(if_else(age>12 & age<18, sxifinfected, 0), na.rm = TRUE),
         symp_12_24 = sum(if_else(age>12 & age<24, sxifinfected, 0), na.rm = TRUE),
         symp_0_12  = sum(if_else(age<12, sxifinfected, 0), na.rm = TRUE),
         symp_0_24  = sum(if_else(age<24, sxifinfected, 0), na.rm = TRUE),
         any_inf_0_6 = ifelse(inf_0_6==0, 0, 1),
         any_inf_6_12 = ifelse(inf_6_12==0, 0, 1),
         any_inf_12_18 = ifelse(inf_12_18==0, 0, 1),
         any_symp_0_6 = ifelse(symp_0_6==0, 0, 1),
         any_symp_6_12 = ifelse(symp_6_12==0, 0, 1),
         any_symp_12_18 = ifelse(symp_12_18==0, 0, 1)
  ) %>%
  dplyr::select(-anyinfection, -sxifinfected, -age) %>%
  distinct()

long_raw_dfff <- bc1 %>%
  dplyr::select(all_of(c("id", "timepoint", ab_columns, "MomFinalRx", "anyHPfinal")))%>%
  pivot_longer(cols=all_of(ab_columns), names_to = "antigen", values_to = "conc")%>%
  filter(antigen %notin% c("logpd", "logGST"))


# read in new visits database to look at correlations with malaria incidence
clin_data <- haven::read_dta("~/postdoc/stanford/clinical_data/BC1/BC-1 childs routine visit database FINAL_ALL.dta")

# it's a big table, so let's only include kids for whom we have any antibody measurements
clin_data <- clin_data %>%
  filter(id %in% long_raw_dfff$id)

# make a data frame where we add a bunch of columns that contain how many (symptomatic) infections were experienced by the child in the indicated time window

infs <- clin_data %>%
  group_by(id) %>%
  dplyr::select(id, age, anyinfection, sxifinfected) %>%
  mutate(inf_0_6   = sum(if_else(age<6, anyinfection, 0), na.rm = TRUE),
         inf_6_12  = sum(if_else(age>6  & age<12, anyinfection, 0), na.rm = TRUE),
         inf_12_18 = sum(if_else(age>12 & age<18, anyinfection, 0), na.rm = TRUE),
         inf_12_24 = sum(if_else(age>12 & age<24, anyinfection, 0), na.rm = TRUE),
         inf_0_12  = sum(if_else(age<12, anyinfection, 0), na.rm = TRUE),
         inf_0_24  = sum(if_else(age<24, anyinfection, 0), na.rm = TRUE),
         symp_0_6   = sum(if_else(age<6, sxifinfected, 0), na.rm = TRUE),
         symp_6_12  = sum(if_else(age>6  & age<12, sxifinfected, 0), na.rm = TRUE),
         symp_12_18 = sum(if_else(age>12 & age<18, sxifinfected, 0), na.rm = TRUE),
         symp_12_24 = sum(if_else(age>12 & age<24, sxifinfected, 0), na.rm = TRUE),
         symp_0_12  = sum(if_else(age<12, sxifinfected, 0), na.rm = TRUE),
         symp_0_24  = sum(if_else(age<24, sxifinfected, 0), na.rm = TRUE),
         any_inf_0_6 = ifelse(inf_0_6==0, 0, 1),
         any_inf_6_12 = ifelse(inf_6_12==0, 0, 1),
         any_inf_12_18 = ifelse(inf_12_18==0, 0, 1),
         any_symp_0_6 = ifelse(symp_0_6==0, 0, 1),
         any_symp_6_12 = ifelse(symp_6_12==0, 0, 1),
         any_symp_12_18 = ifelse(symp_12_18==0, 0, 1)
  ) %>%
  dplyr::select(-anyinfection, -sxifinfected, -age) %>%
  distinct()


# combine antibody data with malaria incidence data
# the -1 removes the id column from the infs df, otherwise it's duplicated         
combo_data <- cbind(long_raw_dfff, infs[match(long_raw_dfff$id, infs$id),-1])

# for these kids we only have antibody measurements at birth and they're not in the clinical database so we'll cut them her
combo_data <- filter(combo_data, id %notin% c(11130, 11084, 11037))

raw_combo_data <- cbind(long_raw_dfff, infs[match(long_raw_dfff$id, infs$id),-1])

kids_with_complete_timecourses <- bc1 %>%
  group_by(id)%>%
  summarise("n_time"=n()) %>%
  filter(n_time==3)%>%
  dplyr::select(id)

combo_data <- filter(combo_data, id %in% kids_with_complete_timecourses$id)

birth <- combo_data %>%
  filter(timepoint==1, !is.na(conc))

six_months <- combo_data %>%
  filter(timepoint==2, !is.na(conc))

twelve_months <- combo_data %>%
  filter(timepoint==3, !is.na(conc))


fdr_cutoff <- 0.1

# 
# #inspect data
# combo_data %>%
#   filter(antigen %in% modelable_antigens)%>%
#   ggplot(aes(x=inf_0_12, y=conc, fill=factor(inf_0_12)))+
#   geom_boxplot()+
#   facet_grid(timepoint~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)))+
#   scale_y_log10(labels=scales::label_log())+
#   xlab("Number of Infections in First Year of Life")+
#   ylab("Concentration")+
#   theme_minimal()+
#   theme(panel.grid = element_blank(),
#         legend.position = "none")+
#   viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)
# 
# 
# summary_combo_data %>%
#   ggplot(aes(x=inf_0_12, y=sum, fill=factor(inf_0_12)))+
#   geom_boxplot()+
#   facet_grid(~timepoint)+
#   scale_y_log10(labels=scales::label_log())+
#   geom_smooth()+
#   xlab("Number of Infections in First Year of Life")+
#   ylab("Sum of all Measured IgG")+
#   theme_minimal()+
#   theme(panel.grid = element_blank(),
#         legend.position = "none")+
#   viridis::scale_fill_viridis(option="A", direction = -1, discrete = TRUE)
# 
# 

# linear regression on general temporal dynamics ####



sec_contrast <- t(matrix(c(0,1,0)))
ter_contrast <- t(matrix(c(0,0,1)))
sec_ter_contrast <- t(matrix(c(0,-1,1)))

purrrf <- combo_data %>%
  group_by(antigen) %>%
  filter(!is.na(factor(MomFinalRx)), !is.na(factor(anyHPfinal)))%>%
  nest() %>%
  mutate(model=map(data, ~lme4::lmer(conc~factor(timepoint)+(1|id), data=.))) %>%
  mutate(model_add_tx=map(data, ~lme4::lmer(conc~factor(timepoint)+factor(MomFinalRx)+(1|id), data=.))) %>%
  mutate(model_x_tx=map(data, ~lme4::lmer(conc~factor(timepoint)*factor(MomFinalRx)+(1|id), data=.))) %>%
  mutate(model_add_hp=map(data, ~lme4::lmer(conc~factor(timepoint)+factor(anyHPfinal)+(1|id), data=.))) %>%
  mutate(model_x_hp=map(data, ~lme4::lmer(conc~factor(timepoint)*factor(anyHPfinal)+(1|id), data=.))) %>%
  mutate(hp_inf_model=map(data, ~glm(conc~factor(timepoint)*factor(anyHPfinal)+(1|id), data=.))) %>%
  
 
  # mutate(summary_model_add_tx=map(model_add_tx, ~summary(.))) %>%
  # mutate(summary_x_tx=map(model_x_tx, ~summary(.)))%>%
  # mutate(summary_model_add_hp=map(model_add_hp, ~summary(.))) %>%
  # mutate(summary_model_x_hp=map(model_x_hp, ~summary(.)))%>%
  
  mutate(AIC_ref=map_dbl(model, ~AIC(.))) %>%
  mutate(AIC_model_add_tx=map_dbl(model_add_tx, ~AIC(.))) %>%
  mutate(AIC_x_tx=map_dbl(model_x_tx, ~AIC(.)))%>%
  mutate(AIC_model_add_hp=map_dbl(model_add_hp, ~AIC(.))) %>%
  mutate(AIC_model_x_hp=map_dbl(model_x_hp, ~AIC(.)))%>%
  # 
  mutate(AIC_delta_model_add_tx=AIC_model_add_tx-AIC_ref) %>%
  mutate(AIC_delta_x_tx=AIC_x_tx-AIC_ref)%>%
  mutate(AIC_delta_model_add_hp=AIC_model_add_hp-AIC_ref) %>%
  mutate(AIC_delta_model_x_hp=AIC_model_x_hp-AIC_ref)%>%
  
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(t2_t1=map(model, ~multcomp::glht(., sec_contrast)),
         t2_t1_p=map_dbl(t2_t1, ~summary(.)$test$pvalues)) %>%
  mutate(t3_t1=map(model, ~multcomp::glht(., ter_contrast)),
         t3_t1_p=map_dbl(t3_t1, ~summary(.)$test$pvalues)) %>%
  mutate(t3_t2=map(model, ~multcomp::glht(., sec_ter_contrast)),
         t3_t2_p=map_dbl(t3_t2, ~summary(.)$test$pvalues))%>%
  ungroup()%>%
  mutate(t2_t1_padj=p.adjust(t2_t1_p),
         t3_t1_padj=p.adjust(t3_t1_p),
         t3_t2_padj=p.adjust(t3_t2_p))
  

#make results table that includes the coef & p_adj values for the sec_ter_contrast
results_table <- purrrf %>%
  mutate(coef=map_dbl(model, ~-coef(.)$id[1,2]+coef(.)$id[1,3]))%>%
  dplyr::select(antigen, coef, t2_t1_padj, t3_t1_padj, t3_t2_padj) %>%
  ungroup()


sig_2_3_abs <- results_table %>%
  filter(t3_t2_padj<fdr_cutoff) %>%
  dplyr::select(antigen)

sig_1_2_abs <- results_table %>%
  filter(t2_t1_padj<fdr_cutoff) %>%
  dplyr::select(antigen)

sig_1_3_abs <- results_table %>%
  filter(t3_t1_padj<fdr_cutoff) %>%
  dplyr::select(antigen)


# vis modelling results
combo_data %>%
  filter(antigen %in% sig_2_3_abs$antigen) %>% 
  ggplot(., aes(x=factor(timepoint), y=conc))+
  facet_wrap(~antigen, scales = "free")+
  geom_violin(aes(fill=anyHPfinal), draw_quantiles = c(0.25, 0.5, 0.75))+
  theme_minimal()+
  scale_color_manual(values=pc1_cols)+
  scale_fill_manual(values=pc1_cols)+
  xlab("")+
  ylab("Concentration")+
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")

# linear regression to test associations between antibody concentrations and malaria incidence ####

# For correlations with protection, we should do the following:  
# 1 look at correlations between antibody levels at birth and anyinfection from birth to 12 months of age
# 2 look at correlations between antibody levels at 6 months of age and anyinfection from 6 to 12 months of age
# 3 look at correlations between antibody levels at 12 months of age and anyinfection from 12-24 months of age
# 4 look at correlations between antibody levels at 12 months of age and sxifinfected from 12-24 months of age
# 
# We could also look in the past to see if antibodies are associated with exposure


# check suitability for incidence data to be modelled by poisson.. doesn't look ideal so we go with negative binomial
# mean_var_summary <- long_combo %>%
#   dplyr::select(id, incidence_type, incidence_value)%>%
#   distinct()%>%
#   group_by(incidence_type)%>%
#   nest()%>%
#   mutate("pois_test"=map(data, ~energy::poisson.tests(as.integer(.$incidence_value), R=10, test = "all")))
# 
# summarise("variance"=var(incidence_value, na.rm = T), "mean"=mean(incidence_value, na.rm=T), "ratio"=variance/mean)
# 



# 1 look at correlations between antibody levels at birth and anyinfection from birth to 12 months of age
birth_purf <- birth %>%
  group_by(antigen)%>%
  nest() %>%
  mutate(model_0_12=map(data, ~glm.nb(inf_0_12 ~ conc, data=.)))%>%
  mutate(model_0_6=map(data, ~glm.nb(inf_0_6 ~ conc, data=.)))%>%
  mutate(model_0_12_symp=map(data, ~glm.nb(symp_0_12 ~ conc, data=.)))%>%
  mutate(model_0_6_symp=map(data, ~glm.nb(symp_0_12 ~ conc, data=.)))%>%
  mutate(model_0_12_symp_prob=map(data, ~glm(symp_0_12/inf_0_12~conc, data=., family = "binomial", weights = inf_0_12)))%>%
  mutate(model_0_6_symp_prob=map(data, ~glm(symp_0_6/inf_0_6~conc, data=., family = "binomial", weights = inf_0_6)))%>%
  mutate(model_any_symp_06=map(data, ~glm(any_symp_0_6~conc, data=., family = "binomial")))%>%
  mutate(model_any_0_6=map(data, ~glm(any_inf_0_6~conc, data=., family = "binomial")))%>%
  
  mutate(hp_model=map(data, ~lm(conc~factor(anyHPfinal), data=.)))%>%
  mutate(tx_model=map(data, ~lm(conc~factor(MomFinalRx), data=.)))%>%
  mutae(tx_inf=map(data, ~glm.nb(inf_0_6~factor(MomFinalRx), data=.)))%>%
  mutate(hp_inf=map(data, ~glm.nb(inf_0_6~factor(anyHPfinal), data=.)))%>%
  
  
  mutate(summary_0_12=map(model_0_12, ~summary(.))) %>%
  mutate(summary_0_6=map(model_0_6, ~summary(.))) %>%
  mutate(summary_0_12_symp=map(model_0_12_symp, ~summary(.))) %>%
  mutate(summary_0_6_symp=map(model_0_6_symp, ~summary(.))) %>%
  mutate(summary_0_12_symp_prob=map(model_0_12_symp_prob, ~summary(.))) %>%
  mutate(summary_0_6_symp_prob=map(model_0_6_symp_prob, ~summary(.))) %>%
  mutate(summary_any_symp_06=map(model_any_symp_06, ~summary(.))) %>%
  mutate(summary_any_0_6=map(model_any_0_6, ~summary(.))) %>%
  
  mutate(hp_model_summary=map(hp_model, ~summary(.))) %>%
  mutate(tx_model_summary=map(tx_model, ~summary(.))) %>%
  mutate(tx_inf_summary=map(tx_inf, ~summary(.))) %>%
  mutate(hp_inf_summary=map(hp_inf, ~summary(.))) %>%
  
  
  mutate(summary_0_12_p=map_dbl(summary_0_12, ~unlist(.$coefficients[8])))%>%
  mutate(summary_0_6_p=map_dbl(summary_0_6, ~unlist(.$coefficients[8])))%>%
  mutate(summary_0_12_symp_p=map_dbl(summary_0_12_symp, ~unlist(.$coefficients[8])))%>%
  mutate(summary_0_6_symp_p=map_dbl(summary_0_6_symp, ~unlist(.$coefficients[8])))%>%
  mutate(summary_0_12_symp_prob_p=map_dbl(summary_0_12_symp_prob, ~unlist(.$coefficients[8])))%>%
  mutate(summary_0_6_symp_prob_p=map_dbl(summary_0_6_symp_prob, ~unlist(.$coefficients[8])))%>%
  mutate(summary_any_symp_06_p=map_dbl(summary_any_symp_06, ~unlist(.$coefficients[8])))%>%
  mutate(summary_any_0_6_p=map_dbl(summary_any_0_6, ~unlist(.$coefficients[8])))%>%
  
  
  ungroup()%>%
  mutate(summary_0_12_padj=p.adjust(summary_0_12_p))%>%
  mutate(summary_0_6_padj=p.adjust(summary_0_6_p))%>%
  mutate(summary_0_12_symp_padj=p.adjust(summary_0_12_symp_p))%>%
  mutate(summary_0_6_symp_padj=p.adjust(summary_0_6_symp_p))%>%
  mutate(summary_0_12_symp_padj=p.adjust(summary_0_12_symp_prob_p))%>%
  mutate(summary_0_6_symp_padj=p.adjust(summary_0_6_symp_prob_p))%>%
  mutate(summary_any_symp_06_padj=p.adjust(summary_any_symp_06_p))%>%
  mutate(summary_any_0_6_padj=p.adjust(summary_any_0_6_p))
  

zero_six_sig <- birth_purf %>%
  filter(summary_0_6_p<fdr_cutoff)%>%
  dplyr::select(antigen, summary_0_6_p)

zero_twelve_sig <- birth_purf %>%
  filter(summary_0_12_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_0_12_padj)


zero_six_symp_sig <- birth_purf %>%
  filter(summary_0_6_symp_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_0_6_symp_padj)

zero_twelve_symp_sig <- birth_purf %>%
  filter(summary_0_12_symp_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_0_12_symp_padj)


zero_six_symp_prob_sig <- birth_purf %>%
  filter(summary_0_6_symp_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_0_6_symp_padj)

zero_twelve_symp_prob_sig <- birth_purf %>%
  filter(summary_0_12_symp_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_0_6_symp_padj)


zero_six_any_sig <- birth_purf %>%
  filter(summary_any_0_6_padj <fdr_cutoff)%>%
  dplyr::select(antigen, summary_any_0_6_padj)

zero_six_any_symp_sig <- birth_purf %>%
  filter(summary_any_symp_06_padj < fdr_cutoff)%>%
  dplyr::select(antigen, summary_any_symp_06_padj)


birth %>%
  filter(antigen %in% zero_six_sig$antigen)%>%
  ggplot(aes(x=inf_0_6, y=conc, fill=factor(inf_0_6)))+
  geom_boxplot()+
  geom_point(alpha=0.2)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
  xlab("Number of Infections in First Six Months of Life")+
  ylab("Concentration in Cord Blood")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)






# 2 look at correlations between antibody levels at 6 months of age and anyinfection from 6 to 12 months of age
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


six_12_sig <- six_purf %>%
  filter(summary_6_12_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_6_12_padj)

six_12_symp_sig <- six_purf %>%
  filter(summary_6_12_symp_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_6_12_symp_padj)

six_12_symp_prob_sig <- six_purf %>%
  filter(summary_6_12_symp_prob_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_6_12_symp_prob_padj)

six_0_6_sig <- six_purf %>%
  filter(summary_0_6_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_0_6_padj)

six_any_12_symp_sig <- six_purf %>%
  filter(summary_any_symp_6_12_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_0_6_padj)

six_any_12_sig <- six_purf %>%
  filter(summary_any_6_12_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_0_6_padj)


six_months %>%
  filter(antigen %in% six_0_6_sig$antigen)%>%
  ggplot(aes(x=inf_0_6, y=conc, fill=factor(inf_0_6)))+
  geom_boxplot()+
  geom_point(alpha=0.2)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)))+
  xlab("Number of Infections in Months 0-6")+
  ylab("Concentration at 6 Months")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)



# 3 look at correlations between antibody levels at 12 months of age and anyinfection from 12-24 months of age
# 4 look at correlations between antibody levels at 12 months of age and sxifinfected from 12-24 months of age


twelve_purf <- twelve_months %>%
  filter(inf_12_18!=0)%>%
  group_by(antigen) %>%
  nest() %>%
  mutate(model_12_18=map(data, ~glm.nb(inf_12_18 ~ conc, data=.)))%>%
  mutate(model_12_18_symp=map(data, ~glm.nb(symp_12_18 ~ conc, data=.)))%>%
  mutate(model_12_18_symp_prob=map(data, ~glm(symp_12_18/inf_12_18~conc, data=., family = "binomial", weights = inf_12_18)))%>%
  mutate(model_6_12=map(data, ~glm.nb(inf_6_12 ~ conc, data=.)))%>%
  mutate(model_0_12=map(data, ~glm.nb(inf_0_12 ~ conc, data=.)))%>%
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
  
  mutate(model_12_24=map(data, ~glm.nb(inf_12_24 ~ conc, data=.)))%>%
  mutate(model_12_24_symp=map(data, ~glm.nb(symp_12_24 ~ conc, data=.)))%>%
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



twelve_months %>%
  # filter(antigen %in% twelve_18_sig$antigen)%>%
  filter(antigen %in% c("logSEA", "logCSP"))%>%
  ggplot(aes(x=any_inf_12_18, y=conc, fill=factor(any_inf_12_18)))+
  geom_violin()+
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
  filter(antigen %in% twelve_18_symp_sig$antigen)%>%
  ggplot(aes(x=symp_12_18, y=conc, fill=factor(symp_12_18)))+
  geom_boxplot()+
  geom_point(alpha=0.2)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
  scale_y_log10(labels=scales::label_log())+
  xlab("Number of Symptomatic Infections in Months 12-18")+
  ylab("Concentration at 12 Months")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)

twelve_months %>%
  filter(antigen %in% twelve_24_sig$antigen)%>%
  ggplot(aes(x=inf_12_24, y=conc, fill=factor(inf_12_24)))+
  geom_boxplot()+
  geom_point(alpha=0.2)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
  scale_y_log10(labels=scales::label_log())+
  xlab("Number of Infections in Months 12-24")+
  ylab("Concentration at 12 Months")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)



twelve_months %>%
  filter(antigen %in% twelve_24_symp_sig$antigen)%>%
  ggplot(aes(x=symp_12_24, y=conc, fill=factor(symp_12_24)))+
  geom_boxplot()+
  geom_point(alpha=0.2)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
  scale_y_log10(labels=scales::label_log())+
  xlab("Number of Symptomatic Infections in Months 12-24")+
  ylab("Concentration at 12 Months")+
  ggpubr::stat_cor(method = "pearson", label.x = 1.5, label.y = 1)+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)



twelve_months %>%
  filter(antigen %in% c("PfSEA", "GST"))%>%
  ggplot(aes(x=inf_0_12, y=conc, fill=factor(inf_0_12)))+
  geom_boxplot()+
  geom_point(alpha=0.2)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
  scale_y_log10(labels=scales::label_log())+
  xlab("Number of Infections in Months 0-12")+
  ylab("Concentration at 12 Months")+
  ggpubr::stat_cor(method = "pearson", label.x = 1.5, label.y = 1)+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  scale_fill_manual(values=time_palette)


ggs_through_time <- combo_data %>%
  filter(antigen %in% c("GEXP", "PfSEA", "GST"))%>%
  ggplot(aes(x=timepoint, y=conc, fill=timepoint))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(alpha=0.2)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
  scale_y_log10(labels=scales::label_log())+
  ggpubr::stat_cor(method = "pearson", label.x = 1.5, label.y = 1)+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = 1, discrete = TRUE)

ggs_vs_exsure6_12 <- combo_data %>%
  filter(antigen %in% c("GEXP", "PfSEA", "GST"))%>%
  ggplot(aes(x=inf_6_12, y=conc, fill=factor(inf_6_12)))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(alpha=0.2)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
  scale_y_log10(labels=scales::label_log())+
  ggpubr::stat_cor(method = "pearson", label.x = 1.5, label.y = 1)+
  theme_minimal()+
  xlab("number of infections in months 6-12")+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = 1, discrete = TRUE)


ggs_vs_exsure_symp_6_12 <- combo_data %>%
  filter(antigen %in% c("GEXP", "PfSEA", "GST"))%>%
  ggplot(aes(x=symp_6_12, y=conc, fill=factor(symp_6_12)))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(alpha=0.2)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
  scale_y_log10(labels=scales::label_log())+
  ggpubr::stat_cor(method = "pearson", label.x = 1.5, label.y = 1)+
  theme_minimal()+
  xlab("number of symptomatic infections in months 6-12")+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = 1, discrete = TRUE)

ggs_combo_plot <- cowplot::plot_grid(ggs_through_time, ggs_vs_exsure6_12, ggs_vs_exsure_symp_6_12, ncol = 1)
# 
# birth %>%
#   filter(antigen %in% six_12_sig$antigen)%>%
#   ggplot(aes(x=inf_6_12, y=conc, fill=factor(inf_6_12)))+
#   geom_boxplot()+
#   geom_point(alpha=0.2)+
#   facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)))+
#   scale_y_log10(labels=scales::label_log())+
#   xlab("Number of Infections in Months 6-12")+
#   ylab("Concentration")+
#   theme_minimal()+
#   theme(panel.grid = element_blank(),
#         legend.position = "none")+
#   viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)
# 
# 
# six_months %>%
#   filter(antigen %in% six_12_symp_sig$antigen)%>%
#   ggplot(aes(x=symp_6_12, y=conc, fill=factor(symp_6_12)))+
#   geom_boxplot()+
#   geom_point(alpha=0.2)+
#   facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)))+
#   scale_y_log10(labels=scales::label_log())+
#   xlab("Number of Symptomatic Infections in Months 6-12")+
#   ylab("Concentration")+
#   theme_minimal()+
#   theme(panel.grid = element_blank(),
#         legend.position = "none")+
#   viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)

birth %>%
  ggplot(aes(x=inf_0_12, y=inf_12_24))+
  geom_point()+
  xlab("Number of Infections in Months 0-12")+
  ylab("Number of Infections in Months 12-24")+
  geom_smooth(formula = y~x, method = "lm")+
  ggpubr::stat_cor(method = "pearson")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_color_viridis(option="B", direction = -1, discrete = TRUE)



# We could also look in the past to see if antibodies are associated with exposure
# also make a column for the fraction of symp episodes to model those



try <- haven::read_dta("~/postdoc/stanford/clinical_data/BC1/MergedAntibodyData_ChildClinical.dta")


ab_columns <- grep("log", colnames(try), value = TRUE)

demo_columns <- c("id",
                  "wealthcat",
                  "age",
                  "totalmalariapreg",
                  "anymalariapreg",
                  "gender",
                  "gestage",
                  "anyparasitemia",
                  "timepoint",
                  "MomFinalRx",
                  "nonmalariafebrile",
                  "ChildFinalRx",
                  "malaria",
                  "mal0to6",
                  "febrile0to6",
                  "incidentmalaria",
                  "febrile6to12",
                  "mal6to12")

long_slim_bc <- try %>%
  dplyr::select(all_of(c("id", "timepoint", ab_columns)))%>%
  pivot_longer(cols=all_of(ab_columns), names_to = "antigen", values_to = "conc")%>%
  filter(antigen %notin% c("logpd", "logGST"))


# geeeeee ####


gee_inf_purf <- twelve_months %>%
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

#logGEXP, logHSP40
gee_twelve_18_sig <- gee_inf_purf %>%
  filter(inf_gee_padj<fdr_cutoff)%>%
  dplyr::select(antigen, inf_gee_padj)




gee_symp_purf <- twelve_months %>%
  filter(symp_12_18>0)%>%
  group_by(antigen) %>%
  nest() %>%
  mutate(inf_gee_symp=map(data, ~geepack::geeglm(symp_12_18 ~ conc,
                                                 data = ., 
                                                 id = id, 
                                                 family = "poisson",
                                                 corstr = "independence")))%>%
  mutate(inf_gee_symp_summary = map(inf_gee_symp, ~summary(.)))%>%
  mutate(inf_gee_symp_p = map(inf_gee_symp_summary, ~.$coefficients[2,4]))%>%
  mutate(inf_gee_symp_padj = p.adjust(inf_gee_symp_p))


#logH103, logHSP40
gee_symp_twelve_18_sig <- gee_symp_purf %>%
  filter(inf_gee_symp_padj<fdr_cutoff)%>%
  dplyr::select(antigen, inf_gee_symp_padj)



gee_any_symp_purf <- twelve_months %>%
  #filter(symp_12_18>0)%>%
  group_by(antigen) %>%
  nest() %>%
  mutate(inf_gee_any_symp=map(data, ~geepack::geeglm(any_symp_12_18 ~ conc,
                                                 data = ., 
                                                 id = id, 
                                                 family = "binomial",
                                                 corstr = "independence")))%>%
  mutate(inf_gee_any_symp_summary = map(inf_gee_any_symp, ~summary(.)))%>%
  mutate(inf_gee_any_symp_p = map(inf_gee_any_symp_summary, ~.$coefficients[2,4]))%>%
  mutate(inf_gee_any_symp_padj = p.adjust(inf_gee_any_symp_p))

#logMSP2_CH150
gee_any_symp_twelve_18_sig <- gee_any_symp_purf %>%
  filter(inf_gee_any_symp_padj<fdr_cutoff)%>%
  dplyr::select(antigen, inf_gee_any_symp_padj)



gee_any_inf_purf <- twelve_months %>%
  #filter(inf_12_18>0)%>%
  group_by(antigen) %>%
  nest() %>%
  mutate(inf_gee_any_inf=map(data, ~geepack::geeglm(any_inf_12_18 ~ conc,
                                                     data = ., 
                                                     id = id, 
                                                     family = "binomial",
                                                     corstr = "independence")))%>%
  mutate(inf_gee_any_inf_summary = map(inf_gee_any_inf, ~summary(.)))%>%
  mutate(inf_gee_any_inf_p = map(inf_gee_any_inf_summary, ~.$coefficients[2,4]))%>%
  mutate(inf_gee_any_inf_padj = p.adjust(inf_gee_any_inf_p))

#logMSP2_CH150
gee_any_inf_twelve_18_sig <- gee_any_inf_purf %>%
  filter(inf_gee_any_inf_padj<fdr_cutoff)%>%
  dplyr::select(antigen, inf_gee_any_inf_padj)



twelve_months %>%
  filter(inf_12_18>0)%>%
  filter(antigen %in% c("logHSP40", "logH103", "logMSP2_CH150", "logSEA", "logGEXP"))%>%
  ggplot(aes(x=inf_12_18, y=conc, fill=factor(inf_12_18)))+
  geom_boxplot()+
  geom_point(color="lightgrey")+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
  xlab("Number of Infections in Months 12-18")+
  ylab("Concentration at 12 Months")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)



twelve_months %>%
  filter(symp_12_18>0)%>%
  filter(antigen %in% c("logHSP40", "logH103", "logMSP2_CH150", "logSEA", "logGEXP"))%>%
  ggplot(aes(x=symp_12_18, y=conc, fill=factor(symp_12_18)))+
  geom_boxplot()+
  geom_point(color="lightgrey")+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
  xlab("Number of Symptomatic Infections in Months 12-18")+
  ylab("Concentration at 12 Months")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)



twelve_months %>%
  # filter(antigen %in% twelve_18_sig$antigen)%>%
  filter(antigen %in% c("logHSP40", "logH103", "logMSP2_CH150", "logSEA", "logGEXP"))%>%
  ggplot(aes(x=any_symp_12_18, y=conc, fill=factor(any_symp_12_18)))+
  geom_boxplot()+
  geom_point(color="lightgrey")+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
  xlab("Any Symptomatic Infections in Months 12-18")+
  ylab("Concentration at 12 Months")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)




twelve_months %>%
  filter(antigen %in% c("logHSP40", "logH103", "logMSP2_CH150", "logSEA", "logGEXP"))%>%
  ggplot(aes(x=any_inf_12_18, y=conc, fill=factor(any_inf_12_18)))+
  geom_boxplot()+
  geom_point(color="lightgrey")+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
  xlab("Any Infections in Months 12-18")+
  ylab("Concentration at 12 Months")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)




# sandbox ####




# 
# purrrf <- combo_data %>%
#   filter(!is.na(anyHPfinal))%>%
#   group_by(antigen) %>%
#   nest() %>%
#   mutate(model=map(data, ~lme4::lmer(conc~factor(timepoint)+(1|id), data=.))) %>%
#   mutate(model_add_tx=map(data, ~lme4::lmer(conc~factor(timepoint)+MomFinalRx+(1|id), data=.))) %>%
#   mutate(model_x_tx=map(data, ~lme4::lmer(conc~factor(timepoint)*MomFinalRx+(1|id), data=.))) %>%
#   mutate(model_add_hp=map(data, ~lme4::lmer(conc~factor(timepoint)+anyHPfinal+(1|id), data=.))) %>%
#   mutate(model_x_hp=map(data, ~lme4::lmer(conc~factor(timepoint)*anyHPfinal+(1|id), data=.))) %>%
#   
#   mutate(AIC_model=map(model, ~AIC(.))) %>%
#   mutate(AIC_model_add_tx=map(model_add_tx, ~AIC(.))) %>%
#   mutate(AIC_x_tx=map(model_x_tx, ~AIC(.)))%>%
#   mutate(AIC_model_add_hp=map(model_add_hp, ~AIC(.))) %>%
#   mutate(AIC_model_x_hp=map(model_x_hp, ~AIC(.)))%>%
#   
#   # mutate(AIC_diff_model_add_tx=map(AIC_model_add_tx, ~.-AIC(model))) %>%
#   # mutate(AIC_diff_x_tx=map(AIC_x_tx, ~.-AIC(model))) %>%
#   # mutate(AIC_diff_model_add_hp=map(AIC_model_add_hp, ~.-AIC(model))) %>%
#   # mutate(AIC_diff_model_x_hp=map(AIC_model_x_hp, ~.-AIC(model))) %>%
#   mutate(summary=map(model, ~summary(.))) %>%
#   mutate(t2_t1=map(model, ~multcomp::glht(., sec_contrast)),
#          t2_t1_p=map_dbl(t2_t1, ~summary(.)$test$pvalues)) %>%
#   mutate(t3_t1=map(model, ~multcomp::glht(., ter_contrast)),
#          t3_t1_p=map_dbl(t3_t1, ~summary(.)$test$pvalues)) %>%
#   mutate(t3_t2=map(model, ~multcomp::glht(., sec_ter_contrast)),
#          t3_t2_p=map_dbl(t3_t2, ~summary(.)$test$pvalues))%>%
#   ungroup()%>%
#   mutate(t2_t1_padj=p.adjust(t2_t1_p),
#          t3_t1_padj=p.adjust(t3_t1_p),
#          t3_t2_padj=p.adjust(t3_t2_p)
#   )%>%
#   mutate(delta_AIC_model=unlist(AIC_model)-unlist(AIC_model)) %>%
#   mutate(delta_AIC_model_add_tx=unlist(AIC_model_add_tx)-unlist(AIC_model)) %>%
#   mutate(delta_AIC_x_tx=unlist(AIC_x_tx)-unlist(AIC_model))%>%
#   mutate(delta_AIC_model_add_hp=unlist(AIC_model_add_hp)-unlist(AIC_model)) %>%
#   mutate(delta_AIC_model_x_hp=unlist(AIC_model_x_hp)-unlist(AIC_model))
# 
# 
# combo_data %>%
#   filter(!is.na(anyHPfinal))%>%
#   filter(antigen %in% c("logAMA1", "logTT", "logRh2"))%>%
#   ggplot(., aes(x=factor(timepoint), y=conc, fill=factor(anyHPfinal)))+
#   geom_boxplot(position = "dodge2", outlier.alpha = 0, na.rm = TRUE)+
#   geom_point(alpha=0.2, position=position_dodge(width = 0.75), shape=21, na.rm = TRUE)+
#   facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), nrow = 3, scales = "free")+
#   xlab("Timepoint")+
#   ylab("Concentration")+
#   theme_minimal()+
#   theme(panel.grid = element_blank(),
#         strip.text = element_text(size=13.5))+
#   scale_fill_manual(values=n_infection_cols)
# 
# combo_data %>%
#   filter(!is.na(MomFinalRx))%>%
#   #filter(antigen %in% c("logAMA1", "logTT", "logRh2"))%>%
#   ggplot(., aes(x=factor(timepoint), y=conc, fill=factor(MomFinalRx)))+
#   geom_boxplot(position = "dodge2", outlier.alpha = 0, na.rm = TRUE)+
#   #geom_point(alpha=0.2, position=position_dodge(width = 0.75), shape=21, na.rm = TRUE)+
#   facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), nrow = 3, scales = "free")+
#   xlab("Timepoint")+
#   ylab("Concentration")+
#   theme_minimal()+
#   theme(panel.grid = element_blank(),
#         strip.text = element_text(size=13.5))+
#   scale_fill_manual(values=n_infection_cols)





purrrf2 <- combo_data %>%
  group_by(antigen) %>%
  filter(!is.na(factor(MomFinalRx)), !is.na(factor(anyHPfinal)))%>%
  #filter(timepoint==1)%>%
  nest() %>%
  mutate(model=map(data, ~lm(conc~1, data=.))) %>%
  mutate(model_add_tx=map(data, ~lm(conc~1+factor(MomFinalRx), data=.))) %>%
  mutate(model_x_tx=map(data, ~lm(conc~1+factor(MomFinalRx):factor(MomFinalRx), data=.))) %>%
  mutate(model_add_hp=map(data, ~lm(conc~1+factor(anyHPfinal), data=.))) %>%
  mutate(model_x_hp=map(data, ~lm(conc~1+factor(anyHPfinal):factor(anyHPfinal), data=.))) %>%
  
  mutate(summary_model_add_tx=map(model_add_tx, ~summary(.))) %>%
  mutate(summary_x_tx=map(model_x_tx, ~summary(.)))%>%
  mutate(summary_model_add_hp=map(model_add_hp, ~summary(.))) %>%
  mutate(summary_model_x_hp=map(model_x_hp, ~summary(.)))%>%
  
  mutate(AIC_ref=map_dbl(model, ~AIC(.))) %>%
  mutate(AIC_model_add_tx=map_dbl(model_add_tx, ~AIC(.))) %>%
  mutate(AIC_x_tx=map_dbl(model_x_tx, ~AIC(.)))%>%
  mutate(AIC_model_add_hp=map_dbl(model_add_hp, ~AIC(.))) %>%
  mutate(AIC_model_x_hp=map_dbl(model_x_hp, ~AIC(.)))%>%
  # 
  mutate(AIC_delta_model_add_tx=AIC_model_add_tx-AIC_ref) %>%
  mutate(AIC_delta_x_tx=AIC_x_tx-AIC_ref)%>%
  mutate(AIC_delta_model_add_hp=AIC_model_add_hp-AIC_ref) %>%
  mutate(AIC_delta_model_x_hp=AIC_model_x_hp-AIC_ref)#%>%



View(purrrf2 %>%
       dplyr::select(AIC_delta_model_add_tx, AIC_delta_model_add_hp)
     )


# ew p values; same results as AIC, thankfully
# leaving out the filtering step returns four "significant" hits, logGEXP    logH103    logAMA1     logRh2, but nothing by AIC(which i think is more believable)
hp_pvals <- sapply(purrrf2$summary_model_add_hp, function(x) x$coefficients[8])
names(hp_pvals) <- purrrf2$antigen
hp_pvals_adj <- p.adjust(hp_pvals)
# 
# combo_data %>%
#   filter(!is.na(anyHPfinal))%>%
#   filter(antigen %in% c("logAMA1", "logTT", "logRh2", "logGEXP", "logGLURP", "logSEA", "logHyp2", "logSBP1"))%>%
#   ggplot(., aes(x=factor(timepoint), y=conc, fill=factor(anyHPfinal)))+
#   geom_boxplot(position = "dodge2", outlier.alpha = 0, na.rm = TRUE)+
#   geom_point(alpha=0.2, position=position_dodge(width = 0.75), shape=21, na.rm = TRUE)+
#   facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), nrow = 4, scales = "free")+
#   xlab("Timepoint")+
#   ylab("Concentration")+
#   theme_minimal()+
#   theme(panel.grid = element_blank(),
#         strip.text = element_text(size=13.5))+
#   scale_fill_manual(values=n_infection_cols)
# 
# combo_data %>%
#   filter(!is.na(MomFinalRx))%>%
#   filter(antigen %in% c("logAMA1", "logTT", "logRh2", "logGEXP", "logGLURP", "logSEA", "logHyp2", "logSBP1"))%>%
#   ggplot(., aes(x=factor(timepoint), y=conc, fill=factor(MomFinalRx)))+
#   geom_boxplot(position = "dodge2", outlier.alpha = 0, na.rm = TRUE)+
#   geom_point(alpha=0.2, position=position_dodge(width = 0.75), shape=21, na.rm = TRUE)+
#   facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), nrow = 4, scales = "free")+
#   xlab("Timepoint")+
#   ylab("Concentration")+
#   theme_minimal()+
#   theme(panel.grid = element_blank(),
#         strip.text = element_text(size=13.5))+
#   scale_fill_manual(values=n_infection_cols)


# covering my poisson ####
# 
# wide_infs <- combo_data%>%
#   select(-antigen, -timepoint)%>%
#   filter(!is.na(anyHPfinalx), !is.na(MomFinalRxx))%>%
#   distinct(id, inf_0_6, inf_6_12, symp_0_6, symp_6_12, anyHPfinalx, MomFinalRxx)
# 
# 
# 
# inf_model <- MASS::glm.nb(inf_0_6~anyHPfinalx, data =wide_infs)
# symp_model <- glm(symp_0_6~anyHPfinalx, family = "poisson", data =wide_infs)
# 
# inf_model2 <- MASS::glm.nb(inf_0_6~MomFinalRxx, data =wide_infs)
# symp_model2 <- glm(symp_0_6~MomFinalRxx, family = "poisson", data =wide_infs)
# 
# inf_model3 <- glm(inf_6_12~anyHPfinalx, family = "poisson", data =wide_infs)
# symp_model3 <- glm(symp_6_12~anyHPfinalx, family = "poisson", data =wide_infs)
# 
# inf_model4 <- glm(inf_6_12~MomFinalRxx, family = "poisson", data =wide_infs)
# symp_model4 <- glm(symp_6_12~MomFinalRxx, family = "poisson", data =wide_infs)
# 
# 
# ggplot(wide_infs, aes(x=MomFinalRxx, y=inf_0_6, fill=MomFinalRxx, group=MomFinalRxx))+
#   geom_violin()+
#   #scale_y_continuous(limits = c(-1,3))+
#   theme_minimal()+
#   theme(legend.position = "none")
# 
# 
# ggplot(wide_infs, aes(x=anyHPfinalx, y=inf_0_6, fill=anyHPfinalx))+
#   geom_violin()+
#   geom_dotplot(binaxis = "y", stackdir = "center")+
#   #scale_y_continuous(limits = c(-1,3))+
#   theme_minimal()+
#   theme(legend.position = "none")
# 
# 
# 
# 
# 
# 
# 
# suppressWarnings(birth_purf <- birth %>%
#                    group_by(antigen)%>%
#                    nest() %>%
#                    mutate(model_0_12=map(data, ~glm(family="poisson", inf_0_12 ~ conc, data=.)))%>%
#                    mutate(model_0_6=map(data, ~glm(family="poisson", inf_0_6 ~ conc, data=.)))%>%
#                    mutate(model_0_12_symp=map(data, ~glm(family="poisson", symp_0_12 ~ conc, data=.)))%>%
#                    mutate(model_0_6_symp=map(data, ~glm(family="poisson", symp_0_12 ~ conc, data=.)))%>%
#                    mutate(model_0_12_symp_prob=map(data, ~glm(symp_0_12/inf_0_12~conc, data=., family = "binomial", weights = inf_0_12)))%>%
#                    mutate(model_0_6_symp_prob=map(data, ~glm(symp_0_6/inf_0_6~conc, data=., family = "binomial", weights = inf_0_6)))%>%
#                    mutate(model_any_symp_06=map(data, ~glm(any_symp_0_6~conc, data=., family = "binomial")))%>%
#                    mutate(model_any_0_6=map(data, ~glm(any_inf_0_6~conc, data=., family = "binomial")))%>%
#                    
#                    mutate(summary_0_12=map(model_0_12, ~summary(.))) %>%
#                    mutate(summary_0_6=map(model_0_6, ~summary(.))) %>%
#                    mutate(summary_0_12_symp=map(model_0_12_symp, ~summary(.))) %>%
#                    mutate(summary_0_6_symp=map(model_0_6_symp, ~summary(.))) %>%
#                    mutate(summary_0_12_symp_prob=map(model_0_12_symp_prob, ~summary(.))) %>%
#                    mutate(summary_0_6_symp_prob=map(model_0_6_symp_prob, ~summary(.))) %>%
#                    mutate(summary_any_symp_06=map(model_any_symp_06, ~summary(.))) %>%
#                    mutate(summary_any_0_6=map(model_any_0_6, ~summary(.))) %>%
#                    
#                    
#                    mutate(summary_0_12_p=map_dbl(summary_0_12, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_0_6_p=map_dbl(summary_0_6, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_0_12_symp_p=map_dbl(summary_0_12_symp, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_0_6_symp_p=map_dbl(summary_0_6_symp, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_0_12_symp_prob_p=map_dbl(summary_0_12_symp_prob, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_0_6_symp_prob_p=map_dbl(summary_0_6_symp_prob, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_any_symp_06_p=map_dbl(summary_any_symp_06, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_any_0_6_p=map_dbl(summary_any_0_6, ~unlist(.$coefficients[8])))%>%
#                    
#                    
#                    ungroup()%>%
#                    mutate(summary_0_12_padj=p.adjust(summary_0_12_p))%>%
#                    mutate(summary_0_6_padj=p.adjust(summary_0_6_p))%>%
#                    mutate(summary_0_12_symp_padj=p.adjust(summary_0_12_symp_p))%>%
#                    mutate(summary_0_6_symp_padj=p.adjust(summary_0_6_symp_p))%>%
#                    mutate(summary_0_12_symp_padj=p.adjust(summary_0_12_symp_prob_p))%>%
#                    mutate(summary_0_6_symp_padj=p.adjust(summary_0_6_symp_prob_p))%>%
#                    mutate(summary_any_symp_06_padj=p.adjust(summary_any_symp_06_p))%>%
#                    mutate(summary_any_0_6_padj=p.adjust(summary_any_0_6_p))
# )
# 
# zero_six_sig <- birth_purf %>%
#   filter(summary_0_6_p<fdr_cutoff)%>%
#   dplyr::select(antigen, summary_0_6_p)
# 
# birth %>%
#   filter(antigen %in% zero_six_sig$antigen)%>%
#   ggplot(aes(x=inf_0_6, y=conc, fill=factor(inf_0_6)))+
#   geom_boxplot()+
#   geom_point(alpha=0.2)+
#   facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
#   xlab("Number of Infections in First Six Months of Life")+
#   ylab("Concentration in Cord Blood")+
#   theme_minimal()+
#   theme(panel.grid = element_blank(),
#         legend.position = "none",
#         strip.text = element_text(size = 14),
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size=16))+
#   scale_fill_manual(values=n_infection_cols)
# 
# 
# 
# 
# suppressWarnings(six_purf <- six_months %>%
#                    group_by(antigen) %>%
#                    nest() %>%
#                    mutate(model_6_12=map(data, ~glm(family="poisson", inf_6_12 ~ conc, data=.)))%>%
#                    mutate(model_6_12_symp=map(data, ~glm(family="poisson", symp_6_12 ~ conc, data=.)))%>%
#                    mutate(model_6_12_symp_prob=map(data, ~glm(symp_6_12/inf_6_12~conc, data=., family = "binomial", weights = inf_6_12)))%>%
#                    mutate(model_0_6=map(data, ~glm(family="poisson", inf_0_6 ~ conc, data=.)))%>%
#                    mutate(model_any_symp_6_12=map(data, ~glm(any_symp_6_12~conc, data=., family = "binomial")))%>%
#                    mutate(model_any_6_12=map(data, ~glm(any_inf_6_12~conc, data=., family = "binomial")))%>%
#                    
#                    
#                    mutate(summary_6_12=map(model_6_12, ~summary(.))) %>%
#                    mutate(summary_6_12_symp=map(model_6_12_symp, ~summary(.))) %>%
#                    mutate(summary_6_12_symp_prob=map(model_6_12_symp_prob, ~summary(.))) %>%
#                    mutate(summary_0_6=map(model_0_6, ~summary(.))) %>%
#                    mutate(summary_any_symp_6_12=map(model_any_symp_6_12, ~summary(.))) %>%
#                    mutate(summary_any_6_12=map(model_any_6_12, ~summary(.))) %>%
#                    
#                    mutate(summary_6_12_p=map_dbl(summary_6_12, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_6_12_symp_p=map_dbl(summary_6_12_symp, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_6_12_symp_prob_p=map_dbl(summary_6_12_symp_prob, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_0_6_p=map_dbl(summary_0_6, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_any_symp_6_12_p=map_dbl(summary_any_symp_6_12, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_any_6_12_p=map_dbl(summary_any_6_12, ~unlist(.$coefficients[8])))%>%
#                    
#                    ungroup()%>%
#                    mutate(summary_6_12_padj=p.adjust(summary_6_12_p))%>%
#                    mutate(summary_6_12_symp_padj=p.adjust(summary_6_12_symp_p))%>%
#                    mutate(summary_6_12_symp_prob_padj=p.adjust(summary_6_12_symp_prob_p))%>%   
#                    mutate(summary_0_6_padj=p.adjust(summary_0_6_p))%>%
#                    mutate(summary_any_symp_6_12_padj=p.adjust(summary_any_symp_6_12_p))%>%
#                    mutate(summary_any_6_12_padj=p.adjust(summary_any_6_12_p))
# )
# 
# 
# six_0_6_sig <- six_purf %>%
#   filter(summary_0_6_padj<fdr_cutoff)%>%
#   dplyr::select(antigen, summary_0_6_padj)
# 
# six_months %>%
#   filter(antigen %in% six_0_6_sig$antigen)%>%
#   ggplot(aes(x=inf_0_6, y=conc, fill=factor(inf_0_6)))+
#   geom_boxplot()+
#   geom_point(alpha=0.2)+
#   facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
#   xlab("Number of Infections in Months 0-6")+
#   ylab("Concentration at 6 Months")+
#   theme_minimal()+
#   theme(panel.grid = element_blank(),
#         legend.position = "none",
#         strip.text = element_text(size = 14),
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size=16))+
#   scale_fill_manual(values=n_infection_cols)
# 
# 
# 
# 
# suppressWarnings(twelve_purf <- twelve_months %>%
#                    group_by(antigen) %>%
#                    nest() %>%
#                    mutate(model_12_18=map(data, ~glm(family="poisson", inf_12_18 ~ conc, data=.)))%>%
#                    mutate(model_12_18_symp=map(data, ~glm(family="poisson", symp_12_18 ~ conc, data=.)))%>%
#                    mutate(model_12_18_symp_prob=map(data, ~glm(symp_12_18/inf_12_18~conc, data=., family = "binomial", weights = inf_12_18)))%>%
#                    mutate(model_6_12=map(data, ~glm(family="poisson", inf_6_12 ~ conc, data=.)))%>%
#                    mutate(model_0_12=map(data, ~glm(family="poisson", inf_0_12 ~ conc, data=.)))%>%
#                    mutate(model_any_symp_12_18=map(data, ~glm(any_symp_12_18~conc, data=., family = "binomial")))%>%
#                    mutate(model_any_12_18=map(data, ~glm(any_inf_12_18~conc, data=., family = "binomial")))%>%
#                    
#                    mutate(summary_12_18=map(model_12_18, ~summary(.))) %>%
#                    mutate(summary_12_18_symp=map(model_12_18_symp, ~summary(.))) %>%
#                    mutate(summary_12_18_symp_prob=map(model_12_18_symp_prob, ~summary(.))) %>%
#                    mutate(summary_6_12=map(model_6_12, ~summary(.))) %>%
#                    mutate(summary_0_12=map(model_0_12, ~summary(.))) %>%
#                    mutate(summary_any_symp_12_18=map(model_any_symp_12_18, ~summary(.))) %>%
#                    mutate(summary_any_12_18=map(model_any_12_18, ~summary(.))) %>%
#                    
#                    mutate(summary_12_18_p=map_dbl(summary_12_18, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_12_18_symp_p=map_dbl(summary_12_18_symp, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_12_18_symp_prob_p=map_dbl(summary_12_18_symp_prob, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_6_12_p=map_dbl(summary_6_12, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_0_12_p=map_dbl(summary_0_12, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_any_symp_12_18_p=map_dbl(summary_any_symp_12_18, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_any_12_18_p=map_dbl(summary_any_12_18, ~unlist(.$coefficients[8])))%>%
#                    
#                    mutate(model_12_24=map(data, ~glm(family="poisson", inf_12_24 ~ conc, data=.)))%>%
#                    mutate(model_12_24_symp=map(data, ~glm(family="poisson", symp_12_24 ~ conc, data=.)))%>%
#                    mutate(model_12_24_symp_prob=map(data, ~glm(symp_12_24/inf_12_24~conc, data=., family = "binomial", weights = inf_12_24)))%>%
#                    
#                    mutate(summary_12_24=map(model_12_24, ~summary(.))) %>%
#                    mutate(summary_12_24_symp=map(model_12_24_symp, ~summary(.))) %>%
#                    mutate(summary_12_24_symp_prob=map(model_12_24_symp_prob, ~summary(.))) %>%
#                    
#                    mutate(summary_12_24_p=map_dbl(summary_12_24, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_12_24_symp_p=map_dbl(summary_12_24_symp, ~unlist(.$coefficients[8])))%>%
#                    mutate(summary_12_24_symp_prob_p=map_dbl(summary_12_24_symp_prob, ~unlist(.$coefficients[8])))%>%
#                    ungroup()%>%
#                    mutate(summary_12_18_padj= p.adjust(summary_12_18_p))%>%
#                    mutate(summary_12_18_symp_padj=p.adjust(summary_12_18_symp_p))%>%
#                    mutate(summary_12_18_symp_prob_padj=p.adjust(summary_12_18_symp_prob_p))%>%
#                    mutate(summary_6_12_padj=p.adjust(summary_6_12_p))%>%
#                    mutate(summary_12_24_padj=p.adjust(summary_12_24_p))%>%
#                    mutate(summary_12_24_symp_padj=p.adjust(summary_12_24_symp_p))%>%
#                    mutate(summary_12_24_symp_prob_padj=p.adjust(summary_12_24_symp_prob_p))%>%
#                    mutate(summary_0_12_padj=p.adjust(summary_0_12_p))%>%
#                    mutate(summary_any_symp_12_18_padj=p.adjust(summary_any_symp_12_18_p))%>%
#                    mutate(summary_any_12_18_padj=p.adjust(summary_any_12_18_p))
# )
# 
# twelve_18_symp_sig <- twelve_purf %>%
#   filter(summary_12_18_symp_padj<fdr_cutoff)%>%
#   dplyr::select(antigen, summary_12_18_symp_padj)
# 
# 
# twelve_18_symp_sig <- twelve_purf %>%
#   filter(summary_12_18_padj<fdr_cutoff)%>%
#   dplyr::select(antigen, summary_12_18_symp_padj)
# 
# 
# 
# 
# twelve_months %>%
#   filter(antigen %in% twelve_6_12_sig$antigen)%>%
#   ggplot(aes(x=inf_6_12, y=conc, fill=factor(inf_6_12)))+
#   geom_boxplot()+
#   geom_point(alpha=0.2)+
#   facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
#   xlab("Number of Infections in Months 0-6")+
#   ylab("Concentration at 6 Months")+
#   theme_minimal()+
#   theme(panel.grid = element_blank(),
#         legend.position = "none",
#         strip.text = element_text(size = 14),
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size=16))+
#   scale_fill_manual(values=n_infection_cols)
# 
# 
# 
# long_infs <- infs %>%
#   pivot_longer(cols=colnames(infs)[2:19], names_to = "incidence_type", values_to = "incidence_value")%>%
#   group_by(incidence_type)%>%
#   summarise(mean_root=sqrt(mean(incidence_value)), sd=sd(incidence_value))
