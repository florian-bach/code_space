# preamble ####
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_minimal())
library(patchwork)

quadratic_glm_with_covar <- function(model, x_name, x_range, coef_name=NULL, coef_value=NULL){
  
  model_fun = function(x){
    exp(model$coefficients[1])*
      exp(model$coefficients[2])^x*
      exp(model$coefficients[3])^(x^2)}
  
  model_prd = data.frame(seq(from = x_range[1], to = x_range[2], length.out = 100),
                         coef_value)
  colnames(model_prd)=c(x_name, coef_name)
  model_err = predict(model, newdata = model_prd, se.fit = TRUE)
  model_prd$lci = model_err$fit - 1.96 * model_err$se.fit
  model_prd$fit = model_err$fit
  model_prd$uci = model_err$fit + 1.96 * model_err$se.fit
  
  output=list("fun"=model_fun,
              "se"=model_prd)
  
  return(output)
  
}

simple_glmer_visualiser <- function(model, x_name, x_range){
  
  model_fun = function(x){
    exp(coef(model)$id[1,1])*
      exp(coef(model)$id[1,2])^x*
      exp(coef(model)$id[1,3])^(x^2)}
  
  model_prd = data.frame(seq(from = x_range[1], to = x_range[2], length.out = 100))
  colnames(model_prd)=c(x_name)
  model_err = predict(model, newdata = model_prd, se.fit = TRUE, re.form=NA)
  model_prd$lci = model_err$fit - 1.96 * model_err$se.fit
  model_prd$fit = model_err$fit
  model_prd$uci = model_err$fit + 1.96 * model_err$se.fit
  
  output=list("fun"=model_fun,
              "se"=model_prd)
  
  return(output)
}

quadratic_glm_visualiser <- function(model, x_name, x_range){
  
  model_fun = function(x){
    exp(model$coefficients[1])*
      exp(model$coefficients[2])^x*
      exp(model$coefficients[3])^(x^2)}
  
  
  model_prd = data.frame(seq(from = x_range[1], to = x_range[2], length.out = 100))
  colnames(model_prd)=c(x_name)
  model_err = predict(model, newdata = model_prd, se.fit = TRUE, re.form=NA)
  model_prd$lci = model_err$fit - 1.96 * model_err$se.fit
  model_prd$fit = model_err$fit
  model_prd$uci = model_err$fit + 1.96 * model_err$se.fit
  
  output=list("fun"=model_fun,
              "se"=model_prd)
  
  return(output)
  
}

# load data ####
meta_cols <- c("id", "date", "flo_age_in_wks", "flo_age_in_months", "mstatus", "pardens", "study", "gender", "treatmentarm")

## promote data  ####
promote1 <- haven::read_dta("~/postdoc/stanford/clinical_data/PROMOTE/BC-3 childs all visit database FINAL.dta")
promote2 <- haven::read_dta("~/Library/CloudStorage/Box-Box/BC3NaturalHistoryPaper/BC-3 children individual level database FINAL.dta")
#1 male, 2 female
promote_data <- promote1 %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7,
         "flo_age_in_months"=as.numeric(date-dob)%/%30.5)%>%
  mutate(pardens=parsdens)%>%
  mutate("study"="promote")%>%
  mutate(treatmentarm="No DP")%>%
  mutate(gender=promote2$gender[match(id, promote2$id)])%>%
  select(all_of(meta_cols))
  
## mic drop ####
mic_drop <-  haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/MICDroP Data/MICDROP all visit database through July 31st 2025.dta")
mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")

mic_drop_data <- mic_drop %>%
  mutate("study"="micdrop")%>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7,
         "flo_age_in_months"=as.numeric(date-dob)%/%30.5)%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         treatmentarm=case_match(treatmentarm,
                                          1~"No DP",
                                          2~"DP 1 year",
                                          3~"DP 2 years"))%>%
  select(all_of(meta_cols))


## impact ####
impact <- haven::read_dta("~/Library/CloudStorage/Box-Box/IMPACT Study (Busia)/Data/IMPACT all visits database FINAL.dta")

impact_data <- impact %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7,
         "flo_age_in_months"=as.numeric(date-dob)%/%30.5)%>%
  mutate("study"="impact")%>%
  mutate(treatmentarm="No DP")%>%
  select(all_of(meta_cols))

# categorize each month based on infection outcome ####
# plan: categorize every 28 week period as 1 of 4 outcomes: nothing, asymp, symp, comp;

combo_data <- bind_rows(mic_drop_data, promote_data, impact_data)%>%
  mutate(treatmentarm=factor(treatmentarm, levels=c("No DP", "DP 1 year", "DP 2 years")))%>%
  mutate(case_when(gender==1~"male",
                   gender==2~"female"))%>%
  mutate(mom_rx=case_when(study=="impact"~maternal_treatment_arms$treatmentarm[match(id-20000, maternal_treatment_arms$id)],
                          study=="micdrop"~maternal_treatment_arms$treatmentarm[match(id-10000, maternal_treatment_arms$id)],
                          study=="promote"~maternal_treatment_arms$treatmentarm[match(id-30000, maternal_treatment_arms$id)]))
  
duration_df <- combo_data%>%
  group_by(id)%>%
  summarise(duration=max(date)-min(date))

person_year_df <- duration_df%>%
  ungroup()%>%
  summarise(as.numeric(sum(duration))/365)


three_month_labels <- paste0(seq(0, 57, by=3), " to ", seq(3, 60, by=3), " months")
six_month_labels <- paste0(seq(0, 54, by=6), " to ", seq(6, 60, by=6), " months")
twelve_month_labels <- paste0(seq(0, 48, by=12), " to ", seq(12, 60, by=12), " months")

monthly_categories <- combo_data %>%
  filter(!is.na(mstatus), !is.na(flo_age_in_wks), mstatus!=3, mstatus!=4)%>%
  mutate("four_week_increment"=floor(flo_age_in_wks/4))%>%
  mutate(age_quarter=cut(flo_age_in_months, breaks = seq(0,60,by=3), labels=three_month_labels))%>%
  mutate(age_semi=cut(flo_age_in_months, breaks = seq(0,60,by=6), labels=six_month_labels))%>%
  arrange(id, four_week_increment)%>%
  group_by(id, four_week_increment)%>%
  mutate(outcome=ifelse(any(mstatus==2), "complicated",
                        ifelse(any(mstatus==1), "uncomplicated",
                               #ifelse(any(mstatus%in%c(3,4)), "treatment_failure",
                               ifelse(all(mstatus==0)&any(pardens>0), "asymptomatic",
                                      ifelse(all(mstatus==0)&all(pardens==0 | is.na(pardens)), "nothing", "something_else")))))%>%
  mutate(has_comp=any(mstatus==2))%>%
  mutate(has_uncomp=any(mstatus==1))%>%
  mutate(has_hyper = ifelse(any(pardens>100000), "hyperparasitemia", "normal"))%>%
  filter((n() == 1) |
    (has_comp & mstatus == 2) |
      (has_comp==FALSE & has_uncomp & mstatus == 1) |
      (has_hyper=="hyperparasitemia" & pardens==max(pardens))
  )%>%
  group_by(id, outcome)%>%
  mutate(order_of_outcome = seq(1, n()))

  

## count instances of each outcome for each individual ####
long_monthly_categories <- monthly_categories%>%
  tidyr::pivot_wider(names_from = outcome, values_fill = 0,values_from = order_of_outcome,
                     id_cols = c("id", "date", "flo_age_in_wks", "flo_age_in_months", "age_quarter", "age_semi", "four_week_increment", "mstatus", "pardens", "has_hyper", "study", "treatmentarm", "gender", "mom_rx"))%>%
  
  ungroup()%>%
  arrange(id, date)%>%
  group_by(id)%>%
  #these n are not the order at the time, but the number of preceding episodes
  mutate(n_nothing=cummax(nothing))%>%
  mutate(n_complicated=cummax(complicated))%>%
  mutate(n_uncomplicated=cummax(uncomplicated))%>%
  mutate(n_asymptomatic=cummax(asymptomatic))%>%
  mutate(n_any_para=n_uncomplicated+n_asymptomatic+n_complicated)%>%
  mutate(n_any_malaria=n_uncomplicated+n_complicated)%>%
  mutate(age_at_first_malaria=min(flo_age_in_months[n_any_malaria==1]))%>%
  mutate(age_at_first_para=min(flo_age_in_months[n_any_para==1]))%>%
  # mutate(age_at_first_malariaf=ifelse(any(flo_age_in_months==age_at_first_malaria), paste(age_semi), "none"))%>%
  mutate(anyDP=if_else(treatmentarm=="No DP", "no", "yes"),
         bino_complicated=ifelse(mstatus==2, 1, 0),
         bino_symp=ifelse(mstatus%in%c(1,2), 1, 0))%>%
  mutate(increment_since_cessation=case_when(treatmentarm=="No DP"~four_week_increment,
                                             treatmentarm=="DP 1 year"~four_week_increment-13,
                                             treatmentarm=="DP 2 years"~four_week_increment-26))




# summarise complicated malaria incidence####
## grouped by malaria episodes ####
## grouped by parasitemic months ####
all_malaria <- long_monthly_categories%>%
  group_by(id)%>%
  filter(complicated!=0|uncomplicated!=0)%>%
  arrange(id, date)

all_para <- long_monthly_categories%>%
  group_by(id)%>%
  filter(complicated!=0 | uncomplicated!=0 | asymptomatic != 0)%>%
  arrange(id, date)


# complicated malaria ####
## individual level models , all data####

malaria_null_model <- lme4::glmer(bino_complicated ~ 1 + mom_rx + (1|id), family = "binomial", data=all_malaria)
malaria_linear_model <- lme4::glmer(bino_complicated ~ n_any_malaria + mom_rx +  (1|id), family = "binomial", data=all_malaria)
malaria_model <- lme4::glmer(bino_complicated ~ n_any_malaria + I(n_any_malaria^2) + mom_rx +  (1|id), family = "binomial", data=all_malaria)
malaria_model_para <- lme4::glmer(bino_complicated ~ n_any_malaria + I(n_any_malaria^2)+ log10(pardens) + mom_rx +  (1|id), family = "binomial", data=all_malaria)
malaria_model_age <- lme4::glmer(bino_complicated ~ n_any_malaria + I(n_any_malaria^2) + flo_age_in_months + mom_rx +  (1|id), family = "binomial", data=all_malaria)
malaria_model_age_gender <- lme4::glmer(bino_complicated ~ n_any_malaria + I(n_any_malaria^2) + flo_age_in_months+ gender+ mom_rx +  (1|id), family = "binomial", data=all_malaria)
malaria_model_age_gender_para <- lme4::glmer(bino_complicated ~ n_any_malaria + I(n_any_malaria^2) + log10(pardens) + flo_age_in_months+ gender+  mom_rx + (1|id), family = "binomial", data=all_malaria)
AIC(malaria_null_model, malaria_linear_model, malaria_model, malaria_model_para, malaria_model_age, malaria_model_age_gender, malaria_model_age_gender_para)

para_null_model <- lme4::glmer(bino_complicated ~ 1 +  mom_rx + (1|id),  family = "binomial", data=all_para)
para_linear_model <- lme4::glmer(bino_complicated ~ n_any_para + mom_rx +  (1|id),  family = "binomial", data=all_para)
para_model <- lme4::glmer(bino_complicated ~ n_any_para + I(n_any_para^2) + mom_rx +  (1|id),  family = "binomial", data=all_para)
para_model_para <- lme4::glmer(bino_complicated ~ n_any_para + I(n_any_para^2) + log10(pardens)+ mom_rx + (1|id),  family = "binomial", data=all_para)
para_model_age <- lme4::glmer(bino_complicated ~ n_any_para + I(n_any_para^2) + flo_age_in_months +  mom_rx + (1|id),  family = "binomial", data=all_para)
para_model_age_gender <- lme4::glmer(bino_complicated ~ n_any_para + I(n_any_para^2) + flo_age_in_months + gender + mom_rx +  (1|id),  family = "binomial", data=all_para)
para_model_age_gender_para <- lme4::glmer(bino_complicated ~ n_any_para + I(n_any_para^2) + log10(pardens)+flo_age_in_months + gender +  mom_rx + (1|id),  family = "binomial", data=all_para)
AIC(para_null_model, para_linear_model, para_model, para_model_para, para_model_age, para_model_age_gender, para_model_age_gender_para)


## summary level models, all data ####
all_malaria_summary <- all_malaria%>%
  group_by(bino_complicated, n_any_malaria)%>%
  summarise("n_cases"=n())%>%
  tidyr::pivot_wider(names_from = bino_complicated, values_from = n_cases, names_prefix = "comp_", values_fill = 0)%>%
  mutate(total_infections=comp_1+comp_0,
         risk=comp_1/total_infections)%>%
  arrange(n_any_malaria)

all_para_summary <- all_para%>%
  group_by(bino_complicated, n_any_para)%>%
  summarise("n_cases"=n())%>%
  tidyr::pivot_wider(names_from = bino_complicated, values_from = n_cases, names_prefix = "comp_", values_fill = 0)%>%
  mutate(total_infections=comp_1+comp_0,
         risk=comp_1/total_infections)%>%
  arrange(n_any_para)

all_malaria_age_summary <- all_malaria%>%
  group_by(bino_complicated, age_quarter, treatmentarm)%>%
  summarise("n_cases"=n())%>%
  tidyr::pivot_wider(names_from = bino_complicated, values_from = n_cases, names_prefix = "comp_",  values_fill = NA)%>%
  mutate(total_infections=comp_1+comp_0,
         risk=comp_1/total_infections)%>%
  arrange(age_quarter)%>%
  mutate(age_quarter_num=as.numeric(factor(age_quarter)))

all_malaria_maternal_rx_age_summary <- all_malaria%>%
  filter(!is.na(mom_rx))%>%
  group_by(bino_complicated, age_quarter, mom_rx)%>%
  summarise("n_cases"=n())%>%
  tidyr::pivot_wider(names_from = bino_complicated, values_from = n_cases, names_prefix = "comp_",  values_fill = NA)%>%
  mutate(total_infections=comp_1+comp_0,
         risk=comp_1/total_infections)%>%
  arrange(age_quarter)%>%
  mutate(age_quarter_num=as.numeric(factor(age_quarter)))

# plot summary level models ####
malaria_summary_model <- glm(risk ~ n_any_malaria + I(n_any_malaria^2), family = "binomial", weights = total_infections, data=all_malaria_summary)
para_summary_model <- glm(risk ~ n_any_para + I(n_any_para^2) ,  family = "binomial", weights = total_infections, data=all_para_summary)
age_summary_model <- glm(risk ~ age_quarter_num,  family = "binomial", weights = total_infections, data=all_malaria_age_summary)

malaria_summary_model_prediction <- quadratic_glm_visualiser(model=malaria_summary_model,
                                                   x_name =  "n_any_malaria",
                                                   x_range = c(1,20))

(uncomplicated_comp_plot <- all_malaria_summary%>%
    ggplot(., aes(x=n_any_malaria, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_ribbon(data=malaria_summary_model_prediction$se, aes(x=n_any_malaria, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = malaria_summary_model_prediction$fun, colour="black")+
  geom_text(aes(y=0.17, label= paste0("frac(",comp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
  scale_x_continuous(breaks = 1:50, limits=c(1,20))+
  # scale_y_continuous(limits = c(0,0.2), labels = scales::label_percent())+
  # ggtitle("Placebo")+
  ggtitle(paste(sum(all_malaria_summary$comp_0+all_malaria_summary$comp_1), "malaria episodes (", sum(all_malaria_summary$comp_1), "complicated)"))+
  xlab("Order of Malaria Episode")+
  ylab("Risk of Complicated Malaria, Given a Malaria Episode"))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/all_uncomplicated_comp_plot.png", uncomplicated_comp_plot, width=7, height=6, dpi=444, bg="white")

para_summary_model_prediction <- quadratic_glm_visualiser(model=para_summary_model,
                                                   x_name =  "n_any_para",
                                                   x_range = c(1, 25))


(para_comp_plot <- ggplot(all_para_summary, aes(x=n_any_para, y=risk))+
    geom_point(color="darkred")+
    theme_minimal()+
    geom_ribbon(data=para_summary_model_prediction$se, aes(x=n_any_para, ymin = exp(lci), ymax = exp(uci)),
                alpha = 0.2, inherit.aes = FALSE)+
    geom_function(fun = para_summary_model_prediction$fun, colour="black")+
    geom_text(aes(y=0.1, label= paste0("frac(",comp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
    scale_x_continuous(breaks = 1:50, limits=c(1,25))+
    scale_y_continuous(limits = c(0,0.12), labels = scales::label_percent())+
    ggtitle(paste(sum(all_para_summary$comp_0+all_para_summary$comp_1), "parasitemic months (", sum(all_para_summary$comp_1), "complicated)"))+
    xlab("Number of Months with Parasitemia")+
    ylab("Risk of Complicated Malaria, Given Parasitemia"))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/para_comp_plot.png", para_comp_plot, width=7, height=6, dpi=444, bg="white")



(age_comp_plot <- all_malaria_age_summary%>%
    filter(!is.na(age_quarter))%>%
    ggplot(., aes(x=age_quarter, y=risk))+
    geom_point(color="darkred")+
    geom_smooth(method="lm", color="black", aes(x=as.numeric(factor(age_quarter))))+
    theme_minimal()+
    # geom_ribbon(data=para_summary_model_prediction$se, aes(x=n_any_para, ymin = exp(lci), ymax = exp(uci)),
    #             alpha = 0.2, inherit.aes = FALSE)+
    # geom_function(fun = para_summary_model_prediction$fun, colour="black")+
    geom_text(aes(y=0.17, label= paste0("frac(",comp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
    # scale_x_continuous(breaks = 1:50, limits=c(1,12))+
    # scale_y_continuous(limits = c(0,0.2), labels = scales::label_percent())+
    ggtitle("")+
    facet_wrap(~treatmentarm)+
    xlab("Age")+
    ylab("Risk of Complicated Malaria")+
    theme(axis.text.x = element_text(angle=90, hjust=0.5, vjust=1)))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/age_comp_plot.png", age_comp_plot, width=12, height=4, dpi=444)



(mo_rx_age_comp_plot <- all_malaria_maternal_rx_age_summary%>%
    filter(!is.na(age_quarter))%>%
    ggplot(., aes(x=age_quarter, y=risk))+
    geom_point(color="darkred")+
    geom_smooth(method="lm", color="black", aes(x=as.numeric(factor(age_quarter))))+
    theme_minimal()+
    # geom_ribbon(data=para_summary_model_prediction$se, aes(x=n_any_para, ymin = exp(lci), ymax = exp(uci)),
    #             alpha = 0.2, inherit.aes = FALSE)+
    # geom_function(fun = para_summary_model_prediction$fun, colour="black")+
    geom_text(aes(y=0.17, label= paste0("frac(",comp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
    # scale_x_continuous(breaks = 1:50, limits=c(1,12))+
    # scale_y_continuous(limits = c(0,0.2), labels = scales::label_percent())+
    ggtitle("")+
    facet_wrap(~mom_rx)+
    xlab("Age")+
    ylab("Risk of Complicated Malaria")+
    theme(axis.text.x = element_text(angle=90, hjust=0.5, vjust=1)))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/age_comp_plot.png", age_comp_plot, width=12, height=4, dpi=444)




# individual level models, placebo only

## individual level models , all data####
placebo_malaria <- all_malaria%>%
  filter(treatmentarm=="No DP")

placebo_para <- all_para%>%
  filter(treatmentarm=="No DP")

malaria_null_model <- lme4::glmer(bino_complicated ~ 1 + mom_rx + (1|id), family = "binomial", data=placebo_malaria)
malaria_linear_model <- lme4::glmer(bino_complicated ~ n_any_malaria + mom_rx +  (1|id), family = "binomial", data=placebo_malaria)
malaria_model <- lme4::glmer(bino_complicated ~ n_any_malaria + I(n_any_malaria^2) + mom_rx +  (1|id), family = "binomial", data=placebo_malaria)
malaria_model_para <- lme4::glmer(bino_complicated ~ n_any_malaria + I(n_any_malaria^2)+ log10(pardens) + mom_rx +  (1|id), family = "binomial", data=placebo_malaria)
malaria_model_age <- lme4::glmer(bino_complicated ~ n_any_malaria + I(n_any_malaria^2) + flo_age_in_months + mom_rx +  (1|id), family = "binomial", data=placebo_malaria)
malaria_model_age_gender <- lme4::glmer(bino_complicated ~ n_any_malaria + I(n_any_malaria^2) + flo_age_in_months+ gender+ mom_rx +  (1|id), family = "binomial", data=placebo_malaria)
malaria_model_age_gender_para <- lme4::glmer(bino_complicated ~ n_any_malaria + I(n_any_malaria^2) + log10(pardens) + flo_age_in_months+ gender+  mom_rx + (1|id), family = "binomial", data=placebo_malaria)
AIC(malaria_null_model, malaria_linear_model, malaria_model, malaria_model_para, malaria_model_age, malaria_model_age_gender, malaria_model_age_gender_para)

para_null_model <- lme4::glmer(bino_complicated ~ 1 +  mom_rx + (1|id),  family = "binomial", data=placebo_para)
para_linear_model <- lme4::glmer(bino_complicated ~ n_any_para + mom_rx +  (1|id),  family = "binomial", data=placebo_para)
para_model <- lme4::glmer(bino_complicated ~ n_any_para + I(n_any_para^2) + mom_rx +  (1|id),  family = "binomial", data=placebo_para)
para_model_para <- lme4::glmer(bino_complicated ~ n_any_para + I(n_any_para^2) + log10(pardens)+ mom_rx + (1|id),  family = "binomial", data=placebo_para)
para_model_age <- lme4::glmer(bino_complicated ~ n_any_para + I(n_any_para^2) + flo_age_in_months +  mom_rx + (1|id),  family = "binomial", data=placebo_para)
para_model_age_gender <- lme4::glmer(bino_complicated ~ n_any_para + I(n_any_para^2) + flo_age_in_months + gender + mom_rx +  (1|id),  family = "binomial", data=placebo_para)
para_model_age_gender_para <- lme4::glmer(bino_complicated ~ n_any_para + I(n_any_para^2) + log10(pardens)+flo_age_in_months + gender +  mom_rx + (1|id),  family = "binomial", data=placebo_para)
AIC(para_null_model, para_linear_model, para_model, para_model_para, para_model_age, para_model_age_gender, para_model_age_gender_para)


# probability of symptoms ####
para_model <- lme4::glmer(bino_symp ~ n_any_para*treatmentarm + flo_age_in_months+ mom_rx + (1|id),  family = "binomial", data=all_para)
car::Anova(para_model, type = "III")


all_para_symp_summary <- all_para%>%
  group_by(bino_symp, n_any_para, treatmentarm)%>%
  summarise("n_cases"=n())%>%
  tidyr::pivot_wider(names_from = bino_symp, values_from = n_cases, names_prefix = "symp_", values_fill = 0)%>%
  mutate(total_infections=symp_1+symp_0,
         risk=symp_1/total_infections)%>%
  arrange(n_any_para)

(all_para_symp_summary_plot <- all_para_symp_summary%>%
    filter(n_any_para<13)%>%
    ggplot(., aes(x=n_any_para, y=risk))+
    theme_minimal()+
    geom_smooth(method = "lm", color="black")+
    geom_point(color="darkred")+
    geom_text(aes(y=0.3, label= paste0("frac(",symp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
    scale_x_continuous(breaks = 1:50)+
    scale_y_continuous(labels = scales::label_percent())+
    facet_wrap(~treatmentarm, scales="free_x")+
    coord_cartesian(ylim = c(0.25, 0.85))+
    xlab("order of parasitemic month")+
    ylab("P(symptoms | parsitemia)"))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/all_para_symp_summary_plot.png", all_para_symp_summary_plot, width=10, height=5, bg="white")




age_symp_summary <- all_para%>%
  group_by(bino_symp, age_quarter, treatmentarm)%>%
  summarise("n_cases"=n())%>%
  tidyr::pivot_wider(names_from = bino_symp, values_from = n_cases, names_prefix = "symp_", values_fill = 0)%>%
  mutate(total_infections=symp_1+symp_0,
         risk=symp_1/total_infections)%>%
  arrange(age_quarter)

(age_symp_summary_plot <- age_symp_summary%>%
    filter(!is.na(age_quarter))%>%
    mutate(quartn=as.numeric(age_quarter))%>%
    filter(treatmentarm=="No DP" | (treatmentarm=="DP 1 year" & quartn>4) | (treatmentarm=="DP 2 years" & quartn>8))%>%
    # filter(treatmentarm=="No DP" | (treatmentarm=="DP 1 year" & quartn>4) | (treatmentarm=="DP 2 years" & quartn>8))%>%
    ggplot(., aes(x=quartn, y=risk))+
    theme_minimal()+
    geom_smooth(aes(x=quartn, y=risk),inherit.aes = FALSE, method = "lm", color="black")+
    geom_point(color="darkred")+
    geom_text(aes(y=0.3, label= paste0("frac(",symp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
    scale_x_continuous(breaks = unique(as.numeric(age_symp_summary$age_quarter[!is.na(age_symp_summary$age_quarter)])),
                       labels = levels(age_symp_summary$age_quarter)[1:14])+
    scale_y_continuous(labels = scales::label_percent())+
    facet_wrap(~treatmentarm, scales="free_x")+
    coord_cartesian(ylim = c(0.2, 0.9))+
    xlab("order of parasitemic month")+
    ylab("P(symptoms | parsitemia)")+
    theme(axis.text.x = element_text(angle=90)))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/age_symp_summary_plot.png", age_symp_summary_plot, width=10, height=5.5, bg="white")



mom_rx_age_symp_summary <- all_para%>%
  group_by(bino_symp, age_quarter, mom_rx)%>%
  summarise("n_cases"=n())%>%
  tidyr::pivot_wider(names_from = bino_symp, values_from = n_cases, names_prefix = "symp_", values_fill = 0)%>%
  mutate(total_infections=symp_1+symp_0,
         risk=symp_1/total_infections)%>%
  arrange(age_quarter)

(mom_rx_age_symp_summary_plot <- mom_rx_age_symp_summary%>%
    filter(!is.na(age_quarter), !is.na(mom_rx))%>%
    mutate(quartn=as.numeric(age_quarter))%>%
    # filter(treatmentarm=="No DP" | (treatmentarm=="DP 1 year" & quartn>4) | (treatmentarm=="DP 2 years" & quartn>8))%>%
    # filter(treatmentarm=="No DP" | (treatmentarm=="DP 1 year" & quartn>4) | (treatmentarm=="DP 2 years" & quartn>8))%>%
    ggplot(., aes(x=quartn, y=risk))+
    theme_minimal()+
    geom_smooth(aes(x=quartn, y=risk),inherit.aes = FALSE, method = "lm", color="black")+
    geom_point(color="darkred")+
    geom_text(aes(y=0.3, label= paste0("frac(",symp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
    scale_x_continuous(breaks = unique(as.numeric(age_symp_summary$age_quarter[!is.na(age_symp_summary$age_quarter)])),
                       labels = levels(age_symp_summary$age_quarter)[1:14])+
    scale_y_continuous(labels = scales::label_percent())+
    facet_wrap(~mom_rx, scales="free_x")+
    coord_cartesian(ylim = c(0.2, 0.9))+
    xlab("order of parasitemic month")+
    ylab("P(symptoms | parsitemia)")+
    theme(axis.text.x = element_text(angle=90)))
ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/mom_rx_age_symp_summary_plot.png", mom_rx_age_symp_summary_plot, width=10, height=5.5, bg="white")



mom_rx_all_para_symp_summary <- all_para%>%
  filter(!is.na(mom_rx))%>%
  group_by(bino_symp, n_any_para, mom_rx)%>%
  summarise("n_cases"=n())%>%
  tidyr::pivot_wider(names_from = bino_symp, values_from = n_cases, names_prefix = "symp_", values_fill = 0)%>%
  mutate(total_infections=symp_1+symp_0,
         risk=symp_1/total_infections)%>%
  arrange(n_any_para)

(mom_rx_all_para_symp_summary_plot <- mom_rx_all_para_symp_summary%>%
    filter(n_any_para<13)%>%
    ggplot(., aes(x=n_any_para, y=risk))+
    theme_minimal()+
    geom_smooth(method = "lm", color="black")+
    geom_point(color="darkred")+
    geom_text(aes(y=0.3, label= paste0("frac(",symp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
    scale_x_continuous(breaks = 1:50)+
    scale_y_continuous(labels = scales::label_percent())+
    facet_wrap(~mom_rx, scales="free_x")+
    coord_cartesian(ylim = c(0.25, 0.85))+
    xlab("order of parasitemic month")+
    ylab("P(symptoms | parsitemia)"))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/mom_rx_all_para_symp_summary_plot.png", mom_rx_all_para_symp_summary_plot, width=10, height=5, bg="white")

# parasite density ####
malaria_pardens_plot <- all_malaria%>%
  ggplot(., aes(x=factor(n_any_malaria), y=pardens, fill=factor(n_any_malaria)))+
  scale_y_continuous(trans = "log10", labels = scales::label_log(), breaks = c(10^seq(1,6)))+
  geom_boxplot(outliers = F)+
  stat_summary(geom="crossbar", colour="white", fun = "median", width=0.65, linewidth=0.19 )+
  theme_minimal()+
  ylab("parasites / μL")+
  xlab("")+
  scale_fill_manual(values=viridis::rocket(n=20))+
  facet_wrap(~treatmentarm, scales="free")+
  xlab("order of malaria episode")+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/n_malaria_parasite_load.png", malaria_pardens_plot, width=12, height=5, bg="white", dpi=444)


para_pardens_plot <- all_para%>%
  ggplot(., aes(x=factor(n_any_para), y=pardens, fill=factor(n_any_para)))+
  scale_y_continuous(trans = "log10", labels = scales::label_log(), breaks = c(10^seq(1,6)))+
  geom_boxplot(outliers = F)+
  stat_summary(geom="crossbar", colour="white", fun = "median", width=0.65, linewidth=0.19 )+
  theme_minimal()+
  ylab("parasites / μL")+
  xlab("")+
  scale_fill_manual(values=viridis::rocket(n=25))+
  facet_wrap(~treatmentarm, scales="free")+
  xlab("order of parasitemic month")+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/n_para_parasite_load.png", para_pardens_plot, width=12, height=5, bg="white", dpi=444)

all_asymp <- long_monthly_categories%>%
  group_by(id)%>%
  filter(bino_symp==0 & pardens!=0, !is.na(pardens))%>%
  arrange(id, date)

asymp_para_pardens_plot <- all_asymp%>%
  ggplot(., aes(x=factor(n_asymptomatic), y=pardens, fill=factor(n_asymptomatic)))+
  scale_y_continuous(trans = "log10", labels = scales::label_log(), breaks = c(10^seq(1,6)))+
  geom_boxplot(outliers = F)+
  stat_summary(geom="crossbar", colour="white", fun = "median", width=0.65, linewidth=0.19 )+
  theme_minimal()+
  ggtitle("asymptomatic infections")+
  ylab("parasites / μL")+
  xlab("")+
  scale_fill_manual(values=viridis::rocket(n=25))+
  facet_wrap(~treatmentarm, scales="free")+
  xlab("order of parasitemic month")+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/n_asymp_parasite_load.png", asymp_para_pardens_plot, width=12, height=5, bg="white", dpi=444)



age_para_pardens_plot <- all_para%>%
  filter(treatmentarm=="No DP" | (treatmentarm=="DP 1 year" & flo_age_in_wks >= 54) | (treatmentarm=="DP 2 years" & flo_age_in_wks >= 106),
         !is.na(age_quarter))%>%
  ggplot(., aes(x=factor(age_quarter), y=pardens, fill=factor(age_quarter)))+
  scale_y_continuous(trans = "log10", labels = scales::label_log(), breaks = c(10^seq(1,6)))+
  geom_boxplot(outliers = F)+
  stat_summary(geom="crossbar", colour="white", fun = "median", width=0.65, linewidth=0.19 )+
  theme_minimal()+
  ggtitle("age")+
  ylab("parasites / μL")+
  xlab("")+
  scale_fill_manual(values=viridis::rocket(n=25))+
  facet_wrap(~treatmentarm, scales="free")+
  xlab("order of parasitemic month")+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle=90))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/age_parasite_load.png", age_para_pardens_plot, width=12, height=5, bg="white", dpi=444)

# plot idea: timeline of all parasitemias of all kids who ever get complicated malaria
kids_with_comp <- combo_data%>%
  filter(mstatus==2)%>%
  distinct(id)

kids_with_two_comp <- combo_data%>%
  filter(mstatus==2)%>%
  filter(duplicated(id))

kids_with_three_comp <- combo_data%>%
  filter(mstatus==2)%>%
  filter(duplicated(id))%>%
  filter(duplicated(id))

comp_pardens_timeline_plot <- long_monthly_categories%>%
  filter(mstatus %in% c(0:2))%>%
  filter(id %in%kids_with_comp$id)%>%
  filter(!is.na(pardens))%>%
  ggplot(., aes(x=n_any_para, y=pardens))+
  geom_line(aes(group=id), alpha=0.2)+
  geom_point(aes(color=factor(mstatus)))+
  # geom_smooth(method="lm", aes(group=factor(mstatus), color=factor(mstatus)))+
  scale_y_continuous(trans="log10", labels=scales::label_log())+
  facet_wrap(~id, scales="free")+
  scale_color_manual(values=c("grey", "black", "red"))+
  theme_minimal()
ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/comp_pardens_age_plot.png", comp_pardens_timeline_plot, width=24, height=24, bg="white", dpi=444)


## make goncalves like arrow plot ####
parsdens_arrow_data <- all_malaria %>%
  filter(id %in% kids_with_comp$id)%>%
  group_by(id)%>%
  filter(any(pardens>100000))%>%
  mutate(age_at_first_comp=flo_age_in_wks[which(mstatus==2)[1]])%>%
  group_by(id)%>%
  arrange(desc(pardens)) %>%
  mutate(rank = row_number())%>%
  filter(
    mstatus == 2 |
      rank <= 2
  )%>%
  group_by(id, mstatus)%>%
  filter(rank==min(rank) | mstatus==2)%>%
  group_by(id)%>%
  #there are 12 singletons
  filter(n()>=2)%>%
  # mild high density infection after severe malaria;
  # mild high density infection before severe malaria;
  mutate(up_down=case_when(any(pardens == max(pardens) & mstatus==1 & flo_age_in_wks>age_at_first_comp)~"after",
                           any(pardens == max(pardens) & mstatus==1 & flo_age_in_wks<age_at_first_comp)~"before",
                           .default = "complicated malaria is highest parasitemia"))%>%
  select(id, flo_age_in_wks, age_at_first_comp, pardens, mstatus, up_down)%>%
  arrange(id)


parsdens_arrow_counts <- parsdens_arrow_data %>%
  ungroup()%>%
  distinct(id, up_down)%>%
  count(up_down)

labels <- setNames(
  paste0(stringr::str_wrap(parsdens_arrow_counts$up_down, width = 25), " (n=", parsdens_arrow_counts$n, ")"),
  parsdens_arrow_counts$up_down
)

parsdens_arrow_plot <- parsdens_arrow_data %>%
  ggplot(aes(x=flo_age_in_wks, y=pardens, group=factor(id)))+
  geom_point(aes(color=factor(mstatus)))+
  geom_line(arrow = arrow(length=unit(0.10,"cm"), ends="last", type = "closed"))+
  # geom_segment(aes(x=age_at_comp_episode, y=pars_dens_comp, xend=age, yend=parsdens),
  #              arrow = arrow(length=unit(0.10,"cm"), ends="last", type = "open"),
  #              lineend = "butt", # See available arrow types in example above
  #              linejoin = "mitre")+
  scale_y_log10(limits=c(10^2, 10^6))+
  annotation_logticks(sides = "l", color = "grey")+
  # ggtitle("n = 19")+
  ylab("parasites /  μL")+
  xlab("age (weeks)")+
  facet_wrap(~up_down,
             labeller = labeller(up_down = labels))+
  scale_color_manual(values = c("1"="black",
                                "2"="red"))+
  theme_minimal()+
  theme(legend.position = "none")
ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/100k_parasites_only_up_arrow_plot.png", parsdens_arrow_plot, height = 4, width=8, dpi=444, bg="white")
# ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/next_highest_infs_up_arrow_plot.png", parsdens_arrow_plot, height = 4, width=8, dpi=444, bg="white")


parsdens_down_arrow <- parsdens_arrow_data %>%
  filter(up_down=="down")%>%
  ggplot(aes(x=age*52, y=parsdens, group=factor(id)))+
  geom_point(aes(color=complicatedmalaria))+
  geom_line(arrow = arrow(length=unit(0.10,"cm"), ends="last", type = "closed"))+
  # geom_segment(aes(x=age_at_comp_episode, y=pars_dens_comp, xend=age, yend=parsdens),
  #              arrow = arrow(length=unit(0.10,"cm"), ends="last", type = "open"),
  #              lineend = "butt", # See available arrow types in example above
  #              linejoin = "mitre")+
  # facet_wrap(~up_down)+
  scale_y_log10(limits=c(10^2, 10^6))+
  annotation_logticks(sides = "l", color = "grey")+
  side_by_side_theme+
  ggtitle("n = 11")+
  ylab("parasites /  μL")+
  xlab("age (weeks)")+
  scale_color_manual(values = c("uncomplicated"="#FC6A03",
                                "complicated"="darkred"))+
  theme_minimal()+
  theme(legend.position = "none")


ggsave("~/postdoc/stanford/clinical_data/PROMOTE/figures/parsdens_down_arrow_plot.png", parsdens_down_arrow, height = 3, width=4, dpi=444, bg="white")


# no chemoprevention ####
no_dp_malaria_summary <- all_malaria%>%
  filter(treatmentarm=="No DP")%>%
  group_by(bino_complicated, n_any_malaria, treatmentarm)%>%
  summarise("n_cases"=n())%>%
  tidyr::pivot_wider(names_from = bino_complicated, values_from = n_cases, names_prefix = "comp_", values_fill = 0)%>%
  mutate(total_infections=comp_1+comp_0,
         risk=comp_1/total_infections)%>%
  arrange(n_any_malaria)

no_dp_para_summary <- all_para%>%
  filter(treatmentarm=="No DP")%>%
  group_by(bino_complicated, n_any_para, treatmentarm)%>%
  summarise("n_cases"=n())%>%
  tidyr::pivot_wider(names_from = bino_complicated, values_from = n_cases, names_prefix = "comp_", values_fill = 0)%>%
  mutate(total_infections=comp_1+comp_0,
         risk=comp_1/total_infections)%>%
  arrange(n_any_para)

no_dp_malaria_summary_model <- glm(risk ~ n_any_malaria + I(n_any_malaria^2), family = "binomial", weights = total_infections, data=no_dp_malaria_summary)
no_dp_malaria_summary_model_prd <- quadratic_glm_visualiser(no_dp_malaria_summary_model, x_range=c(1,20), x_name="n_any_malaria")

no_dp_para_summary_model <- glm(risk ~ n_any_para + I(n_any_para^2), family = "binomial", weights = total_infections, data=no_dp_para_summary)
no_dp_para_summary_model_prd <- quadratic_glm_visualiser(no_dp_para_summary_model, x_range=c(1,20), x_name="n_any_para")

(no_dp_uncomplicated_comp_plot <- no_dp_malaria_summary%>%
    ggplot(., aes(x=n_any_malaria, y=risk))+
    geom_point(color="darkred")+
    theme_minimal()+
    geom_ribbon(data=no_dp_malaria_summary_model_prd$se, aes(x=n_any_malaria, ymin = exp(lci), ymax = exp(uci)),
                alpha = 0.2, inherit.aes = FALSE)+
    geom_function(fun = malaria_summary_model_prediction$fun, colour="black")+
    geom_text(aes(y=0.17, label= paste0("frac(",comp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
    scale_x_continuous(breaks = 1:50, limits=c(1,20))+
    ggtitle(paste(sum(no_dp_malaria_summary$comp_0+no_dp_malaria_summary$comp_1), "malaria episodes (", sum(no_dp_malaria_summary$comp_1), "complicated)"))+
    xlab("Order of Malaria Episode")+
    ylab("Risk of Complicated Malaria, Given a Malaria Episode"))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/no_drug_uncomplicated_comp_plot.png", no_dp_uncomplicated_comp_plot, width=7, height=6, dpi=444, bg="white")


(no_dp_para_comp_plot <- no_dp_para_summary%>%
    ggplot(., aes(x=n_any_para, y=risk))+
    geom_point(color="darkred")+
    theme_minimal()+
    geom_ribbon(data=no_dp_para_summary_model_prd$se, aes(x=n_any_para, ymin = exp(lci), ymax = exp(uci)),
                alpha = 0.2, inherit.aes = FALSE)+
    geom_function(fun = no_dp_para_summary_model_prd$fun, colour="black")+
    geom_text(aes(y=0.1, label= paste0("frac(",comp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
    scale_x_continuous(breaks = 1:50, limits=c(1,20))+
    scale_y_continuous(limits = c(0,0.15), labels = scales::label_percent())+
    # ggtitle("Placebo")+
    ggtitle(paste(sum(no_dp_para_summary$comp_0+no_dp_para_summary$comp_1), "parasitemic months (", sum(no_dp_para_summary$comp_1), "complicated)"))+
    xlab("order of parasitemic month")+
    ylab("risk of complicated malaria, given parasitemia"))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/no_dp_para_comp_plot.png", no_dp_para_comp_plot, width=7, height=6, dpi=444, bg="white")


# age at first ####
# plot idea: remake with placebo only
all_malaria_age_at_first_summary <- all_malaria%>%
  mutate(age_at_first_bins=cut(age_at_first_malaria, breaks = seq(0,60,by=6), labels=six_month_labels))%>%
  group_by(bino_complicated, n_any_malaria, age_at_first_bins)%>%
  summarise("n_cases"=n())%>%
  tidyr::pivot_wider(names_from = bino_complicated, values_from = n_cases, names_prefix = "comp_",  values_fill = NA)%>%
  mutate(total_infections=comp_1+comp_0,
         risk=comp_1/total_infections)


all_malaria_age_at_first_plot <- all_malaria_age_at_first_summary%>%
  filter(!is.na(age_at_first_bins))%>%
  # filter(age_at_first_bins %in% three_month_labels[1:3])%>%
  ggplot(., aes(x=n_any_malaria, y=risk))+
  geom_text(aes(y=0.2, label= paste0("frac(",comp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
  geom_point()+
  scale_y_continuous(labels = scales::label_percent())+
  scale_x_continuous(breaks = seq(1,18))+
  coord_cartesian(ylim=c(0,0.33))+
  ggtitle("age at first malaria episode in months")+
  ylab("risk of complicated malaria")+
  xlab("order of malaria episode")+
  facet_wrap(~age_at_first_bins, ncol=4, scales="free_x")+
  theme_minimal()

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/all_malaria_age_at_first_plot.png", all_malaria_age_at_first_plot, width=16, height=8, dpi=444, bg="white")

all_para_age_at_first_summary <- all_para%>%
  mutate(age_at_first_bins=cut(age_at_first_para, breaks = seq(0,60,by=6), labels=six_month_labels))%>%
  group_by(bino_complicated, n_any_para, age_at_first_bins)%>%
  summarise("n_cases"=n())%>%
  tidyr::pivot_wider(names_from = bino_complicated, values_from = n_cases, names_prefix = "comp_",  values_fill = NA)%>%
  mutate(total_infections=comp_1+comp_0,
         risk=comp_1/total_infections)


all_para_age_at_first_plot <- ggplot(all_para_age_at_first_summary, aes(x=n_any_para, y=risk))+
  geom_text(aes(y=0.2, label= paste0("frac(",comp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
  geom_point()+
  scale_y_continuous(labels = scales::label_percent())+
  scale_x_continuous(breaks = seq(1,18), limits=c(1,18))+
  coord_cartesian(ylim=c(0,0.33))+
  ggtitle("age at first para episode in months")+
  ylab("risk of complicated malaria")+
  xlab("order of parasitemic month")+
  facet_wrap(~age_at_first_bins, ncol=4, scales="free_x")+
  theme_minimal()

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/all_para_age_at_first_plot.png", all_para_age_at_first_plot, width=16, height=8, dpi=444, bg="white")



all_para_age_at_first_summary2 <- all_para%>%
  mutate(age_at_first_bins=cut(age_at_first_malaria, breaks = seq(0,60,by=6), labels=six_month_labels))%>%
  group_by(bino_symp, n_any_para, age_at_first_bins)%>%
  summarise("n_cases"=n())%>%
  tidyr::pivot_wider(names_from = bino_symp, values_from = n_cases, names_prefix = "symp_",  values_fill = NA)%>%
  mutate(total_infections=symp_1+symp_0,
         risk=symp_1/total_infections)


all_para_age_at_first_plot2 <- ggplot(all_para_age_at_first_summary2, aes(x=n_any_para, y=risk))+
  geom_text(aes(y=0.2, label= paste0("frac(",symp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
  geom_point()+
  scale_y_continuous(labels = scales::label_percent())+
  # scale_x_continuous(breaks = seq(1,18), limits=c(1,18))+
  # coord_cartesian(ylim=c(0,0.33))+
  ggtitle("age at first malaria episode in months")+
  ylab("risk of malaria, given parasitemia")+
  xlab("order of parasitemic month")+
  facet_wrap(~age_at_first_bins, ncol=4, scales="free_x")+
  theme_minimal()

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/clinical_immunity_all_para_age_at_first_plot2.png", all_para_age_at_first_plot2, width=16, height=8, dpi=444, bg="white")


placebo_all_malaria_age_at_first_summary <- all_malaria%>%
  filter(treatmentarm=="No DP")%>%
  mutate(age_at_first_bins=cut(age_at_first_malaria, breaks = seq(0,60,by=3), labels=three_month_labels))%>%
  group_by(bino_complicated, n_any_malaria, age_at_first_bins)%>%
  summarise("n_cases"=n())%>%
  tidyr::pivot_wider(names_from = bino_complicated, values_from = n_cases, names_prefix = "comp_",  values_fill = NA)%>%
  mutate(total_infections=comp_1+comp_0,
         risk=comp_1/total_infections)


placebo_all_malaria_age_at_first_plot <- placebo_all_malaria_age_at_first_summary%>%
  filter(!is.na(age_at_first_bins))%>%
  filter(age_at_first_bins %in% three_month_labels[1:3])%>%
  ggplot(., aes(x=n_any_malaria, y=risk))+
  geom_text(aes(y=0.2, label= paste0("frac(",comp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
  geom_point()+
  scale_y_continuous(labels = scales::label_percent())+
  scale_x_continuous(breaks = seq(1,18))+
  coord_cartesian(ylim=c(0,0.33))+
  ggtitle("age at first malaria episode in months")+
  ylab("risk of complicated malaria")+
  xlab("order of malaria episode")+
  facet_wrap(~age_at_first_bins, ncol=4, scales="free_x")+
  theme_minimal()

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/placebo_all_malaria_age_at_first_plot.png", placebo_all_malaria_age_at_first_plot, width=16, height=4, dpi=444, bg="white")



placebo_all_para_age_at_first_summary <- all_para%>%
  filter(treatmentarm=="No DP")%>%
  mutate(age_at_first_bins=cut(age_at_first_para, breaks = seq(0,60,by=3), labels=three_month_labels))%>%
  group_by(bino_complicated, n_any_para, age_at_first_bins)%>%
  summarise("n_cases"=n())%>%
  tidyr::pivot_wider(names_from = bino_complicated, values_from = n_cases, names_prefix = "comp_",  values_fill = NA)%>%
  mutate(total_infections=comp_1+comp_0,
         risk=comp_1/total_infections)


placebo_all_para_age_at_first_plot <- placebo_all_para_age_at_first_summary%>%
  filter(!is.na(age_at_first_bins))%>%
  filter(age_at_first_bins %in% three_month_labels[1:3])%>%
  ggplot(., aes(x=n_any_para, y=risk))+
  geom_text(aes(y=0.2, label= paste0("frac(",comp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
  geom_point()+
  scale_y_continuous(labels = scales::label_percent())+
  scale_x_continuous(breaks = seq(1,18))+
  coord_cartesian(ylim=c(0,0.33))+
  ggtitle("age at first para episode in months")+
  ylab("risk of complicated malaria")+
  xlab("order of parasitemic month")+
  facet_wrap(~age_at_first_bins, ncol=4, scales="free_x")+
  theme_minimal()

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/placebo_all_para_age_at_first_plot.png", placebo_all_para_age_at_first_plot, width=16, height=4, dpi=444, bg="white")


### no DP, individual ####
no_dp_malaria <- filter(all_malaria, treatmentarm=="No DP")
no_dp_para <- filter(all_para, treatmentarm=="No DP")

no_dp_malaria_model <- lme4::glmer(bino_complicated ~ n_any_malaria + I(n_any_malaria^2) + (1|id), family = "binomial", data=no_dp_malaria)
no_dp_para_model <- lme4::glmer(bino_complicated ~ n_any_para + I(n_any_para^2) + (1|id),  family = "binomial", data=no_dp_para)

no_dp_age_malaria_model <- lme4::glmer(bino_complicated ~ flo_age_in_months + (1|id),  family = "binomial", data=no_dp_malaria)
no_dp_age_para_model <- lme4::glmer(bino_complicated ~ flo_age_in_months + (1|id),  family = "binomial", data=no_dp_para)

pardens_malaria_model <- lme4::lmer(log10(pardens) ~ n_any_malaria * treatmentarm+ (1|id),  data=all_malaria)

para_summary_model_prediction <- quadratic_glm_visualiser(model=para_summary_model,
                                                          x_name =  "n_any_para",
                                                          x_range = c(1, 25))


(para_comp_plot <- ggplot(all_para_summary, aes(x=n_any_para, y=risk))+
    geom_point(color="darkred")+
    theme_minimal()+
    geom_ribbon(data=para_summary_model_prediction$se, aes(x=n_any_para, ymin = exp(lci), ymax = exp(uci)),
                alpha = 0.2, inherit.aes = FALSE)+
    geom_function(fun = para_summary_model_prediction$fun, colour="black")+
    geom_text(aes(y=0.1, label= paste0("frac(",comp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
    scale_x_continuous(breaks = 1:50, limits=c(1,25))+
    scale_y_continuous(limits = c(0,0.12), labels = scales::label_percent())+
    ggtitle(paste(sum(all_para_summary$comp_0+all_para_summary$comp_1), "parasitemic months (", sum(all_para_summary$comp_1), "complicated)"))+
    xlab("Number of Months with Parasitemia")+
    ylab("Risk of Complicated Malaria, Given Parasitemia"))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/para_comp_plot.png", para_comp_plot, width=7, height=6, dpi=444, bg="white")

### DP 1 only ####
one_dp_malaria_summary <- all_malaria%>%
  filter(treatmentarm=="DP 1 year")%>%
  group_by(bino_complicated, n_any_malaria, treatmentarm)%>%
  summarise("n_cases"=n())%>%
  pivot_wider(names_from = bino_complicated, values_from = n_cases, names_prefix = "comp_", values_fill = 0)%>%
  mutate(total_infections=comp_1+comp_0,
         risk=comp_1/total_infections)%>%
  arrange(n_any_malaria)

one_dp_malaria_summary_model <- glm(risk ~ n_any_malaria + I(n_any_malaria^2), family = "binomial", weights = total_infections, data=one_dp_malaria_summary)

one_dp_malaria_summary_model_prd <- quadratic_glm_visualiser(one_dp_malaria_summary_model, x_range=c(1,12), x_name="n_any_malaria")


(one_dp_malaria_comp_plot <- one_dp_malaria_summary%>%
    ggplot(., aes(x=n_any_malaria, y=risk))+
    geom_point(color="darkred")+
    theme_minimal()+
    # geom_ribbon(data=one_dp_malaria_summary_model_prd$se, aes(x=n_any_malaria, ymin = exp(lci), ymax = exp(uci)),
    # alpha = 0.2, inherit.aes = FALSE)+
    # geom_function(fun = one_dp_malaria_summary_model_prd$fun, colour="black")+
    geom_text(aes(y=0.17, label= paste0("frac(",comp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
    scale_x_continuous(breaks = 1:50, limits=c(1,10))+
    scale_y_continuous(limits = c(0,0.2), labels = scales::label_percent())+
    # ggtitle("Placebo")+
    xlab("Order of Infection")+
    ylab("Risk of Complicated"))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/one_dp_uncomplicated_comp_plot.png", no_dp_uncomplicated_comp_plot, width=4, height=4, dpi=444)

### DP 1 only ####
any_dp_malaria_summary <- all_malaria%>%
  filter(anyDP=="yes")%>%
  group_by(bino_complicated, n_any_malaria)%>%
  summarise("n_cases"=n())%>%
  pivot_wider(names_from = bino_complicated, values_from = n_cases, names_prefix = "comp_", values_fill = 0)%>%
  mutate(total_infections=comp_1+comp_0,
         risk=comp_1/total_infections)%>%
  arrange(n_any_malaria)

any_dp_malaria_summary_model <- glm(risk ~ n_any_malaria + I(n_any_malaria^2), family = "binomial", weights = total_infections, data=any_dp_malaria_summary)

any_dp_malaria_summary_model_prd <- quadratic_glm_visualiser(any_dp_malaria_summary_model, x_range=c(1,12), x_name="n_any_malaria")


(any_dp_malaria_comp_plot <- any_dp_malaria_summary%>%
    ggplot(., aes(x=n_any_malaria, y=risk))+
    theme_minimal()+
    geom_ribbon(data=any_dp_malaria_summary_model_prd$se, aes(x=n_any_malaria, ymin = exp(lci), ymax = exp(uci)),
                alpha = 0.2, inherit.aes = FALSE)+
    geom_function(fun = any_dp_malaria_summary_model_prd$fun, colour="black")+
    geom_text(aes(y=0.17, label= paste0("frac(",comp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
    geom_point(color="darkred")+
    scale_x_continuous(breaks = 1:50, limits=c(1,5))+
    scale_y_continuous(limits = c(0,0.2), labels = scales::label_percent())+
    # ggtitle("Placebo")+
    xlab("Order of Infection")+
    ylab("Risk of Complicated"))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/any_dp_uncomplicated_comp_plot.png", any_dp_malaria_comp_plot, width=4, height=4, dpi=444)

# sandbox ####

all_mala_symp_summary <- all_para%>%
  # filter(treatmentarm!= "DP 2 years")%>%
  group_by(bino_symp, n_any_malaria, treatmentarm)%>%
  summarise("n_cases"=n())%>%
  tidyr::pivot_wider(names_from = bino_symp, values_from = n_cases, names_prefix = "symp_", values_fill = 0)%>%
  mutate(total_infections=symp_1+symp_0,
         risk=symp_1/total_infections)%>%
  arrange(n_any_malaria)


(all_mala_symp_summary_plot <- ggplot(all_mala_symp_summary, aes(x=n_any_malaria, y=risk))+
    theme_minimal()+
    geom_smooth(method = "lm", color="black")+
    geom_point(color="darkred")+
    # geom_ribbon(data=para_summary_model_prediction$se, aes(x=n_any_para, ymin = exp(lci), ymax = exp(uci)),
    #             alpha = 0.2, inherit.aes = FALSE)+
    # geom_function(fun = para_summary_model_prediction$fun, colour="black")+
    geom_text(aes(y=0.8, label= paste0("frac(",symp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
    scale_x_continuous(breaks = 1:50, limits=c(1,12))+
    scale_y_continuous(limits=c(0.3, 0.85), labels = scales::label_percent())+
    ggtitle("")+
    facet_wrap(~treatmentarm)+
    xlab("Months with any parasitemia")+
    ylab("Probability of Symptoms"))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/all_mala_symp_summary.png", all_mala_symp_summary, width=5, height=15, bg="white", dpi=444)



categories_matrix <- monthly_categories%>%
  filter(study=="micdrop", treatmentarm=="No DP")%>%
  mutate(mstatus_num=as.numeric(mstatus))%>%
  group_by(id, four_week_increment, mstatus_num)%>%
  slice_head(n = 1)%>%
  ungroup()%>%
  tidyr::pivot_wider(names_from = four_week_increment,values_fill = 0, values_from = mstatus_num, id_cols = c("id"))%>%
  arrange(id)

no_zero <- as.matrix(categories_matrix[,2:43])
baseline_dist <- dist(no_zero, method = "euclidean")

baseline_hclust <- hclust(baseline_dist)

inferno <- colorspace::sequential_hcl("inferno", n=9)
col_inferno <- circlize::colorRamp2(seq(min(no_zero), max(no_zero), length.out=length(inferno)), inferno)

internal_split_plot <- Heatmap(matrix = no_zero,
                               #right_annotation=right_anno,
                               #cluster_rows = baseline_hclust,
                               cluster_rows=TRUE,
                               cluster_columns=FALSE,
                               show_row_dend = TRUE,
                               cluster_row_slices = TRUE,
                               show_heatmap_legend = TRUE,
                               # name = "geometric mean log10 parasites in time window",
                               #cluster_columns = FALSE,
                               row_split=8,
                               column_names_gp = gpar(fontsize = 6),
                               column_names_centered=TRUE,
                               # heatmap_legend_param = list(title = "geometric mean log10\nparasitaemia in time window", title_position = "topleft", grid_width=unit(7, "mm")),
                               row_names_gp = gpar(fontsize = 6),
                               row_names_side = "left",
                               col = col_inferno,
                               column_names_rot = 0
)