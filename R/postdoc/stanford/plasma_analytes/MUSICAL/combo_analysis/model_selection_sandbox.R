# model selection ####
## lm() covariate inclusion####
model_selection <- clean_data %>%
  filter(timepoint!="day28", timepoint!="bad_baseline")%>%
  mutate(log_qpcr=log10(qpcr+0.1))%>%
  group_by(targetName)%>%
  nest() %>%
  mutate(simple_model=map(data,  ~lm(concentration~timepoint+infectiontype+id, data=.)), simple_AIC=map_dbl(simple_model, ~AIC(.)))%>%
  mutate(qpcr_model=map(data,    ~lm(concentration~timepoint+infectiontype+id+log_qpcr, data=.)), qpcr_AIC=map_dbl(qpcr_model, ~AIC(.)))%>%
  mutate(age_model=map(data,     ~lm(concentration~timepoint+infectiontype+id+ageyrs, data=.)), age_AIC=map_dbl(age_model, ~AIC(.)))%>%
  mutate(sex_model=map(data,     ~lm(concentration~timepoint+infectiontype+id+gender_categorical, data=.)), sex_AIC=map_dbl(sex_model, ~AIC(.)))%>%
  mutate(age_sex_model=map(data, ~lm(concentration~timepoint+infectiontype+id+gender_categorical+ageyrs, data=.)), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))%>%
  mutate(simple_model_summary=map(simple_model, ~summary(.))) %>%
  mutate(sex_model_summary=map(sex_model, ~summary(.))) %>%
  rowwise()%>%
  mutate("lowest"=min(simple_AIC, age_AIC, sex_AIC, age_sex_AIC, qpcr_AIC))%>%
  mutate("best_model"=list("simple_AIC", "age_AIC", "sex_AIC", "age_sex_AIC", "qpcr_AIC")[which(c(simple_AIC, age_AIC, sex_AIC, age_sex_AIC, qpcr_AIC)==lowest)])%>%
  mutate("better_than_base_model"=(lowest+4)<simple_AIC)%>%
  muate("aic_diff"=)



table(mixed_model_selection$best_model, mixed_model_selection$better_than_base_model==TRUE) #101
age

mixed_model_selection$targetName[mixed_model_selection$best_model=="qpcr_AIC" & mixed_model_selection$better_than_base_model==TRUE]

# --> inlcude age as covariate

## lm() model structure####

model_selection <- clean_data %>%
  filter(timepoint!="day28", timepoint!="bad_baseline")%>%
  group_by(targetName)%>%
  nest() %>%
  mutate(additive_model=map(data,  ~lm(concentration~timepoint+infectiontype+ageyrs+id, data=.)), additive_model_AIC=map_dbl(additive_model, ~AIC(.)))%>%
  mutate(age_interact=map(data,    ~lm(concentration~timepoint+infectiontype*ageyrs+id, data=.)), age_interact_AIC=map_dbl(age_interact, ~AIC(.)))%>%
  mutate(infectiontype_interact=map(data,     ~lm(concentration~timepoint*infectiontype+ageyrs+id, data=.)), infectiontype_interact_AIC=map_dbl(infectiontype_interact, ~AIC(.)))%>%
  # mutate(age_sex_model=map(data, ~lm(concentration~timepoint+infectiontype+ageyrs+gender_categorical+id, data=.)), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))%>%
  rowwise()%>%
  mutate(lowest=min(additive_model_AIC, age_interact_AIC, infectiontype_interact_AIC))%>%
  ungroup()%>%
  mutate(across(c(additive_model_AIC, age_interact_AIC, infectiontype_interact_AIC), ~ .x - lowest))

table(model_selection$additive_model_AIC==0)#129
table(model_selection$age_interact_AIC==0 & model_selection$additive_model_AIC>4)#4
table(model_selection$infectiontype_interact_AIC==0 & model_selection$additive_model_AIC>4)#66
table(model_selection$age_sex_AIC==0 & model_selection$age_AIC>4)#47

# --> simple additive model is best


## model selection plots ####
mixed_model_selection <- clean_data %>%
  filter(timepoint!="day28", timepoint!="bad_baseline", infectiontype%in% c("S", "A"))%>%
  group_by(targetName)%>%
  mutate(log_qpcr=log10(qpcr+0.1))%>%
  nest() %>%
  mutate(simple_model=map(data,  ~lme4::lmer(concentration~timepoint*infectiontype+(1|id), data=.)), simple_AIC=map_dbl(simple_model, ~AIC(.)))%>%
  mutate(qpcr_model=map(data,     ~lme4::lmer(concentration~timepoint*infectiontype+log_qpcr+(1|id), data=.)), qpcr_AIC=map_dbl(qpcr_model, ~AIC(.)))%>%
  # mutate(age_model=map(data,     ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+(1|id), data=.)), age_AIC=map_dbl(age_model, ~AIC(.)))%>%
  # mutate(sex_model=map(data,     ~lme4::lmer(concentration~timepoint*infectiontype+gender_categorical+(1|id), data=.)), sex_AIC=map_dbl(sex_model, ~AIC(.)))%>%
  # mutate(age_sex_model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id), data=.)), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))%>%
  mutate(full_model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+log_qpcr+(1|id), data=.)), full_model_AIC=map_dbl(full_model, ~AIC(.)))%>%
  mutate(no_timepoint_model=map(data, ~lme4::lmer(concentration~log_qpcr*infectiontype+ageyrs+gender_categorical+(1|id), data=.)), no_timepoint_model_AIC=map_dbl(no_timepoint_model, ~AIC(.)))
# rowwise()%>%
# mutate("lowest"=min(simple_AIC, age_AIC, sex_AIC, age_sex_AIC, qpcr_AIC))%>%
# mutate("best_model"=c("simple_AIC", "age_AIC", "sex_AIC", "age_sex_AIC", "qpcr_AIC")[which(c(simple_AIC, age_AIC, sex_AIC, age_sex_AIC, qpcr_AIC)==lowest)])%>%
# mutate("better_than_base_model"=(lowest+4)<simple_AIC)
# ungroup()%>%
# mutate(across(c(simple_AIC, age_AIC, sex_AIC, age_sex_AIC, qpcr_AIC), ~ .x - lowest))

table(mixed_model_selection$best_model, mixed_model_selection$better_than_base_model==TRUE) #101
#               FALSE TRUE
# age_AIC        15   19
# age_sex_AIC     0    1
# qpcr_AIC       22   79
# sex_AIC        14    3
# simple_AIC     97    0

for(i in 1:nrow(mixed_model_selection)){
  
  models <- list(
    "simple model" = mixed_model_selection$simple_model[[i]],
    "qpcr model" = mixed_model_selection$qpcr_model[[i]],
    # "age model" = mixed_model_selection$age_model[[i]],
    # "sex model" = mixed_model_selection$sex_model[[i]],
    # "age sex model" = mixed_model_selection$age_sex_model[[i]],
    "full model" = mixed_model_selection$full_model[[i]],
    "no timepoint model" = mixed_model_selection$no_timepoint_model[[i]]
  )
  
  plt <- ggstats::ggcoef_compare(models, type = "faceted")
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/model_plots/mixed_model_selection_", mixed_model_selection$targetName[[i]], ".png", sep=""), plt, height=6, width=18, dpi=444, bg="white")
  
}

# --> use simple mixed model
# > table(mixed_model_selection$sex_AIC==0 & mixed_model_selection$simple_AIC>4)
# FALSE  TRUE 
# 237    13 
# [1] "ANGPT1"  "CD40LG"  "CRP"     "CTF1"    "GZMB"    "IL12p70" "IL17F"   "IL1RN"   "IL27"    "IRAK4"   "MIF"     "PTX3"    "TNFSF14"

# > table(mixed_model_selection$age_AIC==0 & mixed_model_selection$simple_AIC>4)
# FALSE  TRUE
# 213    37 
# [1] "AGRP"          "BST2"          "CCL17"         "CCL28"         "CD70"          "CD80"          "CSF2"          "CTLA4"        
# [9] "CXCL2"         "EGF"           "FGF2"          "GDF15"         "GFAP"          "IFNA1; IFNA13" "IFNG"          "IKBKG"        
# [17] "IL17RB"        "IL22"          "IL2RA"         "IL33"          "IL7"           "KLRK1"         "LGALS9"        "MUC16"        
# [25] "NCR1"          "NGF"           "S100A9"        "SCG2"          "SELE"          "SELP"          "TEK"           "TGFB1"        
# [33] "THBS2"         "TIMP1"         "TNFSF15"       "VEGFD"         "VSNL1" 



## model selection ####

mixed_model_selection <- clean_data %>%
  filter(timepoint!="day28", timepoint!="bad_baseline", infectiontype%in% c("S"))%>%
  group_by(targetName)%>%
  mutate(log_qpcr=log10(qpcr+0.1))%>%
  nest() %>%
  mutate(simple_model=map(data,  ~lme4::lmer(concentration~timepoint+(1|id), data=.)), simple_AIC=map_dbl(simple_model, ~AIC(.)))%>%
  # mutate(qpcr_model=map(data,     ~lme4::lmer(concentration~timepoint+infectiontype+log_qpcr+(1|id), data=.)), qpcr_AIC=map_dbl(qpcr_model, ~AIC(.)))%>%
  mutate(age_model=map(data,     ~lme4::lmer(concentration~timepoint+ageyrs+(1|id), data=.)), age_AIC=map_dbl(age_model, ~AIC(.)))%>%
  mutate(sex_model=map(data,     ~lme4::lmer(concentration~timepoint+gender_categorical+(1|id), data=.)), sex_AIC=map_dbl(sex_model, ~AIC(.)))%>%
  mutate(age_sex_model=map(data, ~lme4::lmer(concentration~timepoint+ageyrs+gender_categorical+(1|id), data=.)), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))%>%
  mutate(summary=map(simple_model, ~summary(.))) %>%
  mutate(simple_day0_coef=map(simple_model,  ~coef(summary(.))[2]), simple_AIC=map_dbl(simple_model, ~AIC(.)))%>%
  # mutate(qpcr_day0_coef=map(data,     ~lme4::lmer(concentration~timepoint+infectiontype+log_qpcr+(1|id), data=.)), qpcr_AIC=map_dbl(qpcr_model, ~AIC(.)))%>%
  mutate(age_day0_coef=map(age_model,     ~coef(summary(.))[2]), age_AIC=map_dbl(age_model, ~AIC(.)))%>%
  mutate(sex_day0_coef=map(sex_model,     ~coef(summary(.))[2]), sex_AIC=map_dbl(sex_model, ~AIC(.)))%>%
  mutate(age_sex_day0_coef=map(age_sex_model, ~coef(summary(.))[2]), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))%>%
  rowwise()%>%
  mutate("lowest"=min(simple_AIC, age_AIC, sex_AIC, age_sex_AIC))%>%
  mutate("best_model"=c("simple_AIC", "age_AIC", "sex_AIC", "age_sex_AIC")[which(c(simple_AIC, age_AIC, sex_AIC, age_sex_AIC)==lowest)])%>%
  mutate("better_than_base_model"=(lowest+4)<simple_AIC)



slim_mixed_model_selection <- mixed_model_selection %>%
  select(ends_with("_AIC"), ends_with("coef"))%>%
  rowwise()%>%
  mutate("lowest_coef"=min(simple_day0_coef, age_day0_coef, sex_day0_coef, age_sex_day0_coef))%>%
  mutate("highest_coef"=max(simple_day0_coef, age_day0_coef, sex_day0_coef, age_sex_day0_coef))%>%
  mutate("range"=lowest_coef/highest_coef)%>%
  filter()

# for S0, there are only ~10 that vary by 20% or more and the coefs are < | 0.1 |
# for S7 there are only ~10 that vary by 20% or more and the coefs are < | 0.1 |
# for A0 there are ~30 that vary by 20% or more, ~5 that vary by more than 2 fold  
# for A14 there are ~60 that vary by 20% or more, ~10 2 fold, couple of sign changes

# ungroup()%>%
# mutate(across(c(simple_AIC, age_AIC, sex_AIC, age_sex_AIC, qpcr_AIC), ~ .x - lowest))
# 
# table(mixed_model_selection$best_model, mixed_model_selection$better_than_base_model==TRUE) #101



mixed_model_selection <- clean_data %>%
  filter(timepoint!="day28", timepoint!="bad_baseline", infectiontype%in% c("S"))%>%
  group_by(targetName)%>%
  mutate(log_qpcr=log10(qpcr+0.1))%>%
  nest() %>%
  mutate(simple_model=map(data,  ~lme4::lmer(concentration~timepoint+(1|id), data=.)), simple_AIC=map_dbl(simple_model, ~AIC(.)))%>%
  mutate(qpcr_model=map(data,     ~lme4::lmer(concentration~timepoint+day0_qpcr+(1|id), data=.)), qpcr_AIC=map_dbl(qpcr_model, ~AIC(.)))%>%
  mutate(qpcr_i_model=map(data,     ~lme4::lmer(concentration~timepoint*day0_qpcr+(1|id), data=.)), qpcr_AIC=map_dbl(qpcr_model, ~AIC(.)))%>%
  # mutate(age_model=map(data,     ~lme4::lmer(concentration~timepoint+ageyrs+(1|id), data=.)), age_AIC=map_dbl(age_model, ~AIC(.)))%>%
  # mutate(sex_model=map(data,     ~lme4::lmer(concentration~timepoint+gender_categorical+(1|id), data=.)), sex_AIC=map_dbl(sex_model, ~AIC(.)))%>%
  # mutate(age_sex_model=map(data, ~lme4::lmer(concentration~timepoint+ageyrs+gender_categorical+(1|id), data=.)), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))%>%
  mutate(summary=map(qpcr_model, ~summary(.))) %>%
  mutate(simple_day0_coef=map(simple_model,  ~coef(summary(.))[2]), simple_AIC=map_dbl(simple_model, ~AIC(.)))%>%
  mutate(qpcr_day0_coef=map(qpcr_model,     ~coef(summary(.))[2]), qpcr_AIC=map_dbl(qpcr_model, ~AIC(.)))%>%
  mutate(qpcr_day0_i_coef=map(qpcr_i_model,     ~coef(summary(.))[2]), qpcr_i_AIC=map_dbl(qpcr_i_model, ~AIC(.)))%>%
  rowwise()%>%
  mutate("lowest"=min(simple_AIC, qpcr_AIC, qpcr_i_AIC))%>%
  mutate("best_model"=c("simple_AIC", "qpcr_AIC", "qpcr_i_AIC")[which(c(simple_AIC, qpcr_AIC, qpcr_i_AIC)==lowest)])%>%
  mutate("better_than_base_model"=(lowest+4)<simple_AIC)



slim_mixed_model_selection <- mixed_model_selection %>%
  select(ends_with("_AIC"), ends_with("coef"))%>%
  rowwise()%>%
  mutate("lowest_coef"=min(simple_day0_coef, qpcr_day0_coef, qpcr_day0_i_coef))%>%
  mutate("highest_coef"=max(simple_day0_coef, qpcr_day0_coef, qpcr_day0_i_coef))%>%
  mutate("log2FC"=log2(abs(lowest_coef/highest_coef)))%>%
  mutate("lowest"=min(simple_AIC, qpcr_AIC, qpcr_i_AIC))%>%
  mutate("best_model"=list("simple_AIC", "qpcr_AIC", "qpcr_i_AIC")[which(c(simple_AIC, qpcr_AIC, qpcr_i_AIC)==lowest)])%>%
  mutate("better_than_base_model"=(lowest+4)<simple_AIC)


table(abs(slim_mixed_model_selection$log2FC)>log2(1.2)) 
table(mixed_model_selection$best_model, mixed_model_selection$better_than_base_model==TRUE) #101
# FALSE TRUE
# qpcr_AIC      40  160
# qpcr_i_AIC     4   20
# simple_AIC    26    0

# mixed_model_selection$targetName[mixed_model_selection$best_model=="qpcr_i_AIC"]
# [1] "CCL2"      "CCL20"     "CCL3"      "CCL4"      "CD80"      "CSF1R"     "CXCL11"    "CXCL16"   
# [9] "FLT1"      "ICAM1"     "IL10"      "IL1B"      "IL1RN"     "IL34"      "IL6"       "PDGFB"    
# [17] "PTX3"      "TNF"       "TNFRSF11A" "TNFRSF1A"  "TNFRSF1B"  "TREM1"     "TSLP"      "VSTM1" 

#substantial change in day0 coef: log qpcr  #240
#substantial change in day0 coef: day0 qpcr #220

# examining effect of covariates on timepoint ####

# sandbox ####
indie_plot <- clean_data%>%
  filter(infectiontype%in%c("A","S"),
         targetName %in% c("IL10", "CXCL10", "TNF"),
         !is.na(timepoint),
         timepoint %in% c("baseline", "day0"))%>%
  ggplot(aes(x=timepoint, y=concentration, fill=timepoint, color=infectiontype))+
  geom_point()+#
  geom_line(aes(group=id))+
  # geom_boxplot(outliers = FALSE)+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  # ggtitle("regulated during asymptomatic parasitemia")+
  facet_wrap(infectiontype~targetName, scales = "free")+
  scale_color_manual(values=viridis::magma(3))+
  # scale_color_gradientn(colors=viridis::magma(8)[-8])+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/indie_plot.png", indie_plot, width=6, height=6, bg="white", dpi=444)

clean_data%>%
  filter(infectiontype%in%c("A","S"),
         targetName %in% c("CXCL10"),
         !is.na(timepoint))%>%
  ggplot(aes(x=qpcr+00.1, y=concentration, color=infectiontype))+
  geom_point()+#
  geom_smooth(method=lm)+
  ggpubr::stat_cor()+
  scale_x_log10()+
  # geom_line(aes(group=id))+
  # geom_boxplot(outliers = FALSE)+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  # ggtitle("regulated during asymptomatic parasitemia")+
  facet_wrap(~targetName, scales = "free")+
  scale_color_manual(values=viridis::magma(5)[-2])+
  # scale_color_gradientn(colors=viridis::magma(8)[-8])+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/indie_plot.png", indie_plot, width=6, height=6, bg="white", dpi=444)


# misc ####



target_vars <- clean_data %>%
  filter(targetName%notin%c("CTSS"))%>%
  group_by(targetName)%>%
  summarise("mean"=mean(concentration), "var"=var(concentration), "ratio"=var/mean)

most_variable <- slice_max(target_vars, ratio, prop = 0.5)

clean_data%>%
  filter(targetName %in% most_variable_features$targetName,
         infectiontype %notin% c("PS", "A2", "NM"))%>%
  ggplot(aes(x=timepoint, y=concentration, fill=infectiontype))+
  geom_violin()+
  facet_wrap(~targetName)+
  theme_minimal()




model_selection <- clean_data %>%
  filter(timepoint!="day28", timepoint!="bad_baseline")%>%
  group_by(targetName)%>%
  nest() %>%
  mutate(simple_model=map(data,  ~lm(concentration~timepoint*infectiontype+id, data=.)), simple_AIC=map_dbl(simple_model, ~AIC(.)))%>%
  mutate(age_model=map(data,     ~lm(concentration~timepoint*infectiontype+ageyrs+id, data=.)), age_AIC=map_dbl(age_model, ~AIC(.)))%>%
  mutate(sex_model=map(data,     ~lm(concentration~timepoint*infectiontype+gender_categorical+id, data=.)), sex_AIC=map_dbl(sex_model, ~AIC(.)))%>%
  mutate(age_sex_model=map(data, ~lm(concentration~timepoint*infectiontype+ageyrs+gender_categorical+id, data=.)), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))%>%
  rowwise()%>%
  mutate(lowest=min(simple_AIC, age_AIC, sex_AIC, age_sex_AIC))%>%
  ungroup()%>%
  mutate(simplest_is_best=lowest==simple_AIC,
         age_is_best=lowest==age_AIC,
         sex_is_best=lowest==sex_AIC,
         age_sex_is_best=lowest==age_sex_AIC)

model_selection$targetName[model_selection$simplest_is_best]
model_selection$targetName[model_selection$age_is_best]
model_selection$targetName[model_selection$sex_is_best]
model_selection$targetName[model_selection$age_sex_is_best]






#most variable features

topVarProts <- head(order(rowVars(cpm(d_norm)), decreasing = TRUE))



# clean_data <- big_df %>%
#   mutate(id=substr(sample_id, 2,4),
#          infectiontype=substr(sample_id, 6,6),
#          barcode= substr(sample_id, nchar(sample_id)-5, nchar(sample_id)))%>%
#   rowwise()%>%
#   mutate("sample_type"=substr(sample_id,
#                               stringi::stri_locate_all(pattern = "_", sample_id, fixed=TRUE)[[1]][1]+1,
#                               stringi::stri_locate_all(pattern = "_", sample_id, fixed=TRUE)[[1]][2]-1))%>%
#   ungroup()%>%
#   mutate("timepoint"=case_when(sample_type=="A.1"~"baseline",
#                                sample_type=="A.2"~"drop",
#                                sample_type=="A0"~"day0",
#                                sample_type=="A14"~"day14",
#                                sample_type=="A2.1"~"drop",
#                                sample_type=="A20"~"day0",
#                                sample_type=="A214"~"day14",
#                                sample_type=="A28"~"day28",
#                                sample_type=="NM.1"~"baseline",
#                                sample_type=="NM0"~"day0",
#                                sample_type=="NM14"~"day14",
#                                sample_type=="NM7"~"day7",
#                                sample_type=="PS0"~"PS0",
#                                sample_type=="S.1"~"baseline",
#                                sample_type=="S.2"~"drop",
#                                sample_type=="S0"~"day0",
#                                sample_type=="S14"~"day14",
#                                sample_type=="S28"~"day28",
#                                sample_type=="S7"~"day7")
#   )%>%
#   mutate(timepoint=factor(timepoint, levels=c("baseline", "day0", "day7", "day14", "day28")))%>%
#   group_by(targetName) %>%
#   mutate(z_conc=scale(concentration, center = TRUE, scale = TRUE))%>%
#   ungroup()%>%
#   group_by(id)%>%
#   mutate("mean_z_conc"=median(z_conc))%>%
#   ungroup()%>%
#   filter(mean_z_conc > -0.35)
# parasitemia inclussion?


base_zero_contrast <- t(matrix(c(0,1,0,0)))
base_14_contrast <- t(matrix(c(0,0,1,0)))
# base_28_contrast <- t(matrix(c(0,0,0)))
zero_14_contrast <- t(matrix(c(0,-1,1,0)))

base_zero_para_contrast <- t(matrix(c(0,1,0,1)))
base_14_para_contrast <- t(matrix(c(0,0,1,1)))

asymp_only_purff <- clean_data %>%
  filter(infectiontype=="A", timepoint!="day28", timepoint!="bad_baseline")%>%
  group_by(targetName)%>%
  mutate(log_pcr=log10(qpcr+0.1),
         log_pcr_dich=if_else(log_pcr>=2.69, "over 500k/ml", "under 500k/ml"))%>%
  nest() %>%
  # mutate(model=map(data, ~lm(concentration~timepoint+ageyrs+id, data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~timepoint+log_pcr+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_zero=map(model, ~multcomp::glht(., base_zero_contrast)),
         base_zero_p=map_dbl(base_zero, ~summary(.)$test$pvalues)) %>%
  mutate(base_zero_para=map(model, ~multcomp::glht(., base_zero_para_contrast)),
         base_zero_para_p=map_dbl(base_zero_para, ~summary(.)$test$pvalues),
         base_zero_para_coef=map_dbl(summary, ~coef(.)[5]),
         base_zero_para_coef_p=map_dbl(summary, ~coef(.)[5]))%>%
  mutate(base_14_para=map(model, ~multcomp::glht(., base_14_para_contrast)),
         base_14_para_p=map_dbl(base_14_para, ~summary(.)$test$pvalues),
         base_14_para_coef=map_dbl(summary, ~coef(.)[5]),
         base_14_para_coef_p=map_dbl(summary, ~coef(.)[5]))%>%
  ungroup()%>%
  mutate(base_zero_padj=p.adjust(base_zero_p, method="BH"),
         base_14_para_padj=p.adjust(base_14_para_p, method="BH"),
         base_zero_para_padj=p.adjust(base_zero_para_p, method="BH")
  )



asymp_results_table <- asymp_only_purff %>%
  dplyr::select(targetName,
                base_zero_padj,
                base_zero_para_padj,
                base_zero_para_coef,
                base_14_para_padj)%>%
  ungroup()

#79 with lm(); 159 with lme4;
asymp_sig_base_zero <- asymp_results_table %>%
  filter(base_zero_padj<fdr_cutoff)


#0 with lm(); 36 with lme4
asymp_sig_base_0_para <- asymp_only_purff %>%
  filter(base_zero_para_padj<fdr_cutoff)%>%
  arrange(desc(base_zero_para_coef))

asymp_sig_base_14_para <- asymp_only_purff %>%
  filter(base_14_para_padj<fdr_cutoff)%>%
  arrange(desc(base_14_para_padj))


#2 with lme4
symp_sig_base_14 <- symp_results_table %>%
  filter(base_14_padj<fdr_cutoff)%>%
  select(targetName)




clean_data %>%
  filter(infectiontype=="A", timepoint!="day28", timepoint!="bad_baseline")%>%
  filter(targetName %in% head(asymp_sig_base_14_para$targetName, n=29))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  mutate(log_pcr=log10(qpcr+0.1),
         log_pcr_dich=if_else(log_pcr>=2.69, "over 500k/ml", "under 500k/ml"))%>%
  ggplot(aes(x=log_pcr, y=concentration))+
  geom_point(aes(color=timepoint))+#
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman", color="red")+
  # geom_line(aes(group=id, color=qpcr_cat))+
  # geom_boxplot(aes(fill=log_pcr_dich), outliers = FALSE)+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  # ggtitle("regulated during asymptomatic parasitemia")+
  facet_wrap(~targetName, scales = "free")+
  scale_color_manual(values=viridis::magma(6))+
  theme_minimal()





