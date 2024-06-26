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

bc1$age <- bc1$date -bc1$dob


#subset data to include ID, timepoint, all antigens, flags and concentrations
raw_data <- bc1[,c(1,67:132, 155)]

#shove all antigens in their own column, same for flags and concentrations
raw_conc_data <- raw_data %>%
  pivot_longer(cols = colnames(raw_data)[seq(3, 66, by=3)], values_to = "conc")%>%
  dplyr::select(conc)

raw_flag_data <- raw_data%>%
  pivot_longer(cols = colnames(raw_data)[seq(4, 67, by=3)], values_to = "flag")%>%
  dplyr::select(flag)

# make long data frame combining both conc and flag data
long_raw_df <- raw_data %>%
  pivot_longer(cols = colnames(raw_data)[seq(2, 65, by=3)], values_to = "antigen") %>%
  dplyr::select(id, timepoint, antigen)%>%
  mutate(timepoint=factor(timepoint))%>%
  mutate(conc=raw_conc_data$conc, flag=raw_flag_data$flag)%>%
  mutate(antigen=gsub(".", " ", antigen, fixed = TRUE))%>%
  mutate(antigen=gsub("  ", " ", antigen, fixed = TRUE))%>%
  mutate(id=factor(id))%>%
  mutate(flag_type = case_when(flag==1 ~ "AboveMaxStd",
                               flag==2 ~ "AboveUpperbound",
                               flag==3 ~ "Above_fitted_asymptote",
                               flag==4 ~ "BelowMinStd",
                               is.na(flag) ~ "No Flag"))%>%
  filter(antigen != "")
  
# remove all rows that have any flag
flagless_df <- long_raw_df %>%
  filter(is.na(flag))

# make a df that has the sum of all antibodies
long_raw_summary <- flagless_df %>%
  group_by(id, timepoint)%>%
  summarise("sum"=sum(conc, na.rm=TRUE))


# read in new visits database to look at correlations with malaria incidence
clin_data <- haven::read_dta("~/postdoc/stanford/clinical_data/BC1/BC-1 childs routine visit database FINAL_ALL.dta")

# it's a big table, so let's only include kids for whom we have any antibody measurements
clin_data <- clin_data %>%
  filter(id %in% raw_data$id)

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
         symp_0_24  = sum(if_else(age<24, sxifinfected, 0), na.rm = TRUE)
         ) %>%
  dplyr::select(-anyinfection, -sxifinfected, -age) %>%
  distinct()
  
# combine antibody data with malaria incidence data
# the -1 removes the id column from the infs df, otherwise it's duplicated         
combo_data <- cbind(flagless_df, infs[match(flagless_df$id, infs$id),-1])
summary_combo_data <- cbind(long_raw_summary, infs[match(long_raw_summary$id, infs$id),-1])

# for these kids we only have antibody measurements at birth and they're not in the clinical database so we'll cut them her
combo_data <- filter(combo_data, id %notin% c(11130, 11084, 11037))

raw_combo_data <- cbind(long_raw_df, infs[match(long_raw_df$id, infs$id),-1])

kids_with_complete_timecourses <- bc1 %>%
  group_by(id)%>%
  summarise("n_time"=n()) %>%
  filter(n_time==3)%>%
  dplyr::select(id)

combo_data <- filter(combo_data, id %in% kids_with_complete_timecourses$id)

#inspect data
combo_data %>%
  filter(antigen %in% modelable_antigens)%>%
  ggplot(aes(x=inf_0_12, y=conc, fill=factor(inf_0_12)))+
  geom_boxplot()+
  facet_grid(timepoint~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)))+
  scale_y_log10(labels=scales::label_log())+
  xlab("Number of Infections in First Year of Life")+
  ylab("Concentration")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)
  

summary_combo_data %>%
  ggplot(aes(x=inf_0_12, y=sum, fill=factor(inf_0_12)))+
  geom_boxplot()+
  facet_grid(~timepoint)+
  scale_y_log10(labels=scales::label_log())+
  geom_smooth()+
  xlab("Number of Infections in First Year of Life")+
  ylab("Sum of all Measured IgG")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="A", direction = -1, discrete = TRUE)



# linear regression on general temporal dynamics ####



sec_contrast <- t(matrix(c(0,1,0)))
ter_contrast <- t(matrix(c(0,0,1)))
sec_ter_contrast <- t(matrix(c(0,-1,1)))

purrrf <- combo_data %>%
  #filter(antigen %in% modelable_antigens)%>%
  filter(antigen %notin% c("Rh4 2", "EBA181 VIII V"))%>%
  group_by(antigen) %>%
  mutate(log_conc=log10(conc))%>%
  nest() %>%
  mutate(model=map(data, ~lme4::lmer(log_conc~timepoint+(1|id), data=.))) %>%
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
         t3_t2_padj=p.adjust(t3_t2_p)
         )

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


# vis modelling results
combo_data %>%
  filter(antigen %in% sig_2_3_abs$antigen & antigen %in% modelable_antigens) %>% 
  ggplot(., aes(x=timepoint, y=conc))+
  facet_wrap(~antigen, scales = "free")+
  geom_violin(aes(fill=antigen), draw_quantiles = c(0.25, 0.5, 0.75))+
  scale_y_log10(labels=scales::label_log())+
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


birth <- combo_data %>%
  filter(timepoint==1)
             
six_months <- combo_data %>%
  filter(timepoint==2)
             
twelve_months <- combo_data %>%
  filter(timepoint==3)
       

fdr_cutoff <- 0.1
      

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
  group_by(antigen) %>%
  mutate(log_conc=log10(conc))%>%
  nest() %>%
  mutate(model_0_12=map(data, ~glm.nb(inf_0_12 ~ log_conc, data=.)))%>%
  mutate(model_0_6=map(data, ~glm.nb(inf_0_6 ~ log_conc, data=.)))%>%
  mutate(model_0_12_symp=map(data, ~glm.nb(symp_0_12 ~ log_conc, data=.)))%>%
  mutate(model_0_6_symp=map(data, ~glm.nb(symp_0_12 ~ log_conc, data=.)))%>%
  mutate(model_0_12_symp_prob=map(data, ~glm(symp_0_12/inf_0_12~log_conc, data=., family = "binomial", weights = inf_0_12)))%>%
  mutate(model_0_6_symp_prob=map(data, ~glm(symp_0_6/inf_0_6~log_conc, data=., family = "binomial", weights = inf_0_6)))%>%
  
  mutate(summary_0_12=map(model_0_12, ~summary(.))) %>%
  mutate(summary_0_6=map(model_0_6, ~summary(.))) %>%
  mutate(summary_0_12_symp=map(model_0_12_symp, ~summary(.))) %>%
  mutate(summary_0_6_symp=map(model_0_6_symp, ~summary(.))) %>%
  mutate(summary_0_12_symp_prob=map(model_0_12_symp_prob, ~summary(.))) %>%
  mutate(summary_0_6_symp_prob=map(model_0_6_symp_prob, ~summary(.))) %>%
  
  mutate(summary_0_12_p=map_dbl(summary_0_12, ~unlist(.$coefficients[8])))%>%
  mutate(summary_0_6_p=map_dbl(summary_0_6, ~unlist(.$coefficients[8])))%>%
  mutate(summary_0_12_symp_p=map_dbl(summary_0_12_symp, ~unlist(.$coefficients[8])))%>%
  mutate(summary_0_6_symp_p=map_dbl(summary_0_6_symp, ~unlist(.$coefficients[8])))%>%
  mutate(summary_0_12_symp_prob_p=map_dbl(summary_0_12_symp_prob, ~unlist(.$coefficients[8])))%>%
  mutate(summary_0_6_symp_prob_p=map_dbl(summary_0_6_symp_prob, ~unlist(.$coefficients[8])))%>%
  ungroup()%>%
  mutate(summary_0_12_padj=p.adjust(summary_0_12_p))%>%
  mutate(summary_0_6_padj=p.adjust(summary_0_6_p))%>%
  mutate(summary_0_12_symp_padj=p.adjust(summary_0_12_symp_p))%>%
  mutate(summary_0_6_symp_padj=p.adjust(summary_0_6_symp_p))%>%
  mutate(summary_0_12_symp_padj=p.adjust(summary_0_12_symp_prob_p))%>%
  mutate(summary_0_6_symp_padj=p.adjust(summary_0_6_symp_prob_p))

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




birth %>%
  filter(antigen %in% zero_twelve_sig$antigen)%>%
  ggplot(aes(x=inf_0_12, y=conc, fill=factor(inf_0_12)))+
  geom_boxplot()+
  geom_point(alpha=0.2)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
  scale_y_log10(labels=scales::label_log())+
  xlab("Number of Infections in First Year of Life")+
  ylab("Concentration in Cord Blood")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)


birth %>%
  filter(antigen %in% zero_six_sig$antigen)%>%
  ggplot(aes(x=inf_0_6, y=conc, fill=factor(inf_0_6)))+
  geom_boxplot()+
  geom_point(alpha=0.2)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
  scale_y_log10(labels=scales::label_log())+
  xlab("Number of Infections in First Six Months of Life")+
  ylab("Concentration in Cord Blood")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)



birth %>%
  filter(antigen %in% zero_six_symp_prob_sig$antigen)%>%
  mutate(symp_prob_0_6=symp_0_6/inf_0_6)%>%
  #mutate(symp_prob_0_6=if_else(is.na(symp_prob_0_6), -0.25, symp_prob_0_6))%>%
  ggplot(aes(x=symp_prob_0_6, y=conc, color=factor(inf_0_6), shape=factor(symp_0_6)))+
  #geom_point()+
  geom_text(aes(y=conc, label= paste0("frac(",symp_0_6, ",", inf_0_6,")")),parse = TRUE, size=2.5)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
  scale_y_log10(labels=scales::label_log())+
  xlab("Fraction of Infections With Symptoms in First Six Months of Life")+
  ylab("Concentration in Cord Blood")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)







# 2 look at correlations between antibody levels at 6 months of age and anyinfection from 6 to 12 months of age
six_purf <- six_months %>%
  filter(antigen %in% modelable_antigens)%>%
  group_by(antigen) %>%
  mutate(log_conc=log10(conc))%>%
  nest() %>%
  mutate(model_6_12=map(data, ~glm.nb(inf_6_12 ~ log_conc, data=.)))%>%
  mutate(model_6_12_symp=map(data, ~glm.nb(symp_6_12 ~ log_conc, data=.)))%>%
  mutate(model_6_12_symp_prob=map(data, ~glm(symp_6_12/inf_6_12~log_conc, data=., family = "binomial", weights = inf_6_12)))%>%
  mutate(model_0_6=map(data, ~glm.nb(inf_0_6 ~ log_conc, data=.)))%>%
  mutate(summary_6_12=map(model_6_12, ~summary(.))) %>%
  mutate(summary_6_12_symp=map(model_6_12_symp, ~summary(.))) %>%
  mutate(summary_6_12_symp_prob=map(model_6_12_symp_prob, ~summary(.))) %>%
  mutate(summary_0_6=map(model_0_6, ~summary(.))) %>%
  mutate(summary_6_12_p=map_dbl(summary_6_12, ~unlist(.$coefficients[8])))%>%
  mutate(summary_6_12_symp_p=map_dbl(summary_6_12_symp, ~unlist(.$coefficients[8])))%>%
  mutate(summary_6_12_symp_prob_p=map_dbl(summary_6_12_symp_prob, ~unlist(.$coefficients[8])))%>%
  mutate(summary_0_6_p=map_dbl(summary_0_6, ~unlist(.$coefficients[8])))%>%
  ungroup()%>%
  mutate(summary_6_12_padj=p.adjust(summary_6_12_p))%>%
  mutate(summary_6_12_symp_padj=p.adjust(summary_6_12_symp_p))%>%
  mutate(summary_6_12_symp_prob_padj=p.adjust(summary_6_12_symp_prob_p))%>%   
  mutate(summary_0_6_padj=p.adjust(summary_0_6_p))
  

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


six_months %>%
  filter(antigen %in% six_12_sig$antigen)%>%
  ggplot(aes(x=inf_6_12, y=conc, fill=factor(inf_6_12)))+
  geom_boxplot()+
  geom_point(alpha=0.2)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)))+
  scale_y_log10(labels=scales::label_log())+
  xlab("Number of Infections in Months 6-12")+
  ylab("Concentration at 6 Months")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)


six_months %>%
  filter(antigen %in% six_12_symp_sig$antigen)%>%
  ggplot(aes(x=symp_6_12, y=conc, fill=factor(symp_6_12)))+
  geom_boxplot()+
  geom_point(alpha=0.2)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)))+
  scale_y_log10(labels=scales::label_log())+
  xlab("Number of Symptomatic Infections in Months 6-12")+
  ylab("Concentration at 6 Months")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)



six_months %>%
  filter(antigen %in% six_0_6_sig$antigen)%>%
  ggplot(aes(x=inf_0_6, y=conc, fill=factor(inf_0_6)))+
  geom_boxplot()+
  geom_point(alpha=0.2)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)))+
  scale_y_log10(labels=scales::label_log())+
  xlab("Number of Infections in Months 0-6")+
  ylab("Concentration at 6 Months")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)



# 3 look at correlations between antibody levels at 12 months of age and anyinfection from 12-24 months of age
# 4 look at correlations between antibody levels at 12 months of age and sxifinfected from 12-24 months of age


twelve_purf <- twelve_months %>%
  group_by(antigen) %>%
  mutate(log_conc=log10(conc))%>%
  nest() %>%
  mutate(model_12_18=map(data, ~glm.nb(inf_12_18 ~ log_conc, data=.)))%>%
  mutate(model_12_18_symp=map(data, ~glm.nb(symp_12_18 ~ log_conc, data=.)))%>%
  mutate(model_12_18_symp_prob=map(data, ~glm(symp_12_18/inf_12_18~log_conc, data=., family = "binomial", weights = inf_12_18)))%>%
  mutate(model_6_12=map(data, ~glm.nb(inf_6_12 ~ log_conc, data=.)))%>%
  mutate(model_0_12=map(data, ~glm.nb(inf_0_12 ~ log_conc, data=.)))%>%
  
  mutate(summary_12_18=map(model_12_18, ~summary(.))) %>%
  mutate(summary_12_18_symp=map(model_12_18_symp, ~summary(.))) %>%
  mutate(summary_12_18_symp_prob=map(model_12_18_symp_prob, ~summary(.))) %>%
  mutate(summary_6_12=map(model_6_12, ~summary(.))) %>%
  mutate(summary_0_12=map(model_0_12, ~summary(.))) %>%
  
  mutate(summary_12_18_p=map_dbl(summary_12_18, ~unlist(.$coefficients[8])))%>%
  mutate(summary_12_18_symp_p=map_dbl(summary_12_18_symp, ~unlist(.$coefficients[8])))%>%
  mutate(summary_12_18_symp_prob_p=map_dbl(summary_12_18_symp_prob, ~unlist(.$coefficients[8])))%>%
  mutate(summary_6_12_p=map_dbl(summary_6_12, ~unlist(.$coefficients[8])))%>%
  mutate(summary_0_12_p=map_dbl(summary_0_12, ~unlist(.$coefficients[8])))%>%
  
  mutate(model_12_24=map(data, ~glm.nb(inf_12_24 ~ log_conc, data=.)))%>%
  mutate(model_12_24_symp=map(data, ~glm.nb(symp_12_24 ~ log_conc, data=.)))%>%
  mutate(model_12_24_symp_prob=map(data, ~glm(symp_12_24/inf_12_24~log_conc, data=., family = "binomial", weights = inf_12_24)))%>%
  
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
  mutate(summary_0_12_padj=p.adjust(summary_0_12_p))


   


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



twelve_months %>%
  filter(antigen %in% twelve_18_sig$antigen)%>%
  ggplot(aes(x=inf_12_18, y=conc, fill=factor(inf_12_18)))+
  geom_boxplot()+
  geom_point(alpha=0.2)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales="free")+
  scale_y_log10(labels=scales::label_log())+
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


