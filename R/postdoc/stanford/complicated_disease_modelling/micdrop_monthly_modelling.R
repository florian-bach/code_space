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

## promote data  ####
promote <- haven::read_dta("~/postdoc/stanford/clinical_data/PROMOTE/BC-3 childs all visit database FINAL.dta")

promote_data <- promote %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7,
         "flo_age_in_months"=as.numeric(date-dob)%/%30.5)%>%
  mutate(pardens=parsdens)%>%
  mutate("study"="promote")%>%
  mutate(treatmentarm="No DP")%>%
  select(id, date, flo_age_in_wks, flo_age_in_months, mstatus, pardens, study, treatmentarm)
  
## mic drop ####
mic_drop <-  haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/May25/MICDSpecimenBoxMay25_withclinical.dta")
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
  select(id, date, flo_age_in_wks, flo_age_in_months, mstatus, pardens, study, treatmentarm)
  

## impact ####

impact <- haven::read_dta("~/Library/CloudStorage/Box-Box/IMPACT Study (Busia)/Data/Specimens/Mar25/IMPASpecimenBoxMar25_withclinical.dta")

impact_data <- impact %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7,
         "flo_age_in_months"=as.numeric(date-dob)%/%30.5)%>%
  mutate("study"="impact")%>%
  mutate(treatmentarm="No DP")%>%
  select(id, date, flo_age_in_wks, flo_age_in_months, mstatus, pardens, study, treatmentarm)

## categorize each month based on infection outcome ####
combo_data <- bind_rows(mic_drop_data, promote_data, impact_data)%>%
  mutate(treatmentarm=factor(treatmentarm, levels=c("No DP", "DP 1 year", "DP 2 years")))
  
three_month_labels <- paste0(seq(0, 57, by=3), " to ", seq(3, 60, by=3), " months")
# plan: categorize every 28 week period as 1 of 4 outcomes: nothing, asymp, symp, comp;
monthly_categories <- combo_data %>%
  filter(!is.na(pardens), !is.na(mstatus), !is.na(flo_age_in_wks), mstatus!=3, mstatus!=4)%>%
  mutate("four_week_increment"=round(flo_age_in_wks/4))%>%
  mutate(age_quarter=cut(flo_age_in_months, breaks = seq(0,60,by=3), labels=three_month_labels))%>%
  group_by(id, four_week_increment)%>%
  mutate(outcome=ifelse(any(mstatus==2), "complicated",
                        ifelse(any(mstatus==1), "uncomplicated",
                               #ifelse(any(mstatus%in%c(3,4)), "treatment_failure",
                               ifelse(all(mstatus==0)&any(pardens>0), "asymptomatic",
                                      ifelse(all(mstatus==0)&all(pardens==0), "nothing", "something_else")))))%>%
  mutate(has_comp=any(mstatus==2))%>%
  mutate(has_uncomp=any(mstatus==1))%>%
  filter((n() == 1) |
    (has_comp & mstatus == 2) |
      (!has_comp & has_uncomp & mstatus == 1)
  )%>%
  group_by(id, outcome)%>%
  mutate(order_of_outcome = seq(1, n()))
  

## count instances of each outcome for each individual ####
long_monthly_categories <- monthly_categories%>%
  pivot_wider(names_from = outcome, values_fill = 0,values_from = order_of_outcome, id_cols = c("id", "date", "flo_age_in_wks", "flo_age_in_months", "age_quarter", "four_week_increment", "mstatus", "pardens", "study", "treatmentarm"))%>%
    ungroup()%>%
  group_by(id)%>%
  arrange(date)%>%
  #these n are not the order at the time, but the number of preceding episodes
  mutate(n_nothing=cummax(ifelse(is.na(nothing), 0, nothing)))%>%
  mutate(n_complicated=cummax(ifelse(is.na(complicated), 0, complicated)))%>%
  mutate(n_uncomplicated=cummax(ifelse(is.na(uncomplicated), 0, uncomplicated)))%>%
  mutate(n_asymptomatic=cummax(ifelse(is.na(asymptomatic), 0, asymptomatic)))%>%
  mutate(n_any_para=n_uncomplicated+n_asymptomatic+n_complicated)%>%
  mutate(n_any_malaria=n_uncomplicated+n_complicated)%>%
  mutate(anyDP=if_else(treatmentarm=="No DP", "no", "yes"),
         bino_complicated=ifelse(mstatus==2, 1, 0),
         bino_symp=ifelse(mstatus%in%c(1,2), 1, 0))%>%
  mutate(increment_since_cessation=case_when(treatmentarm=="No DP"~four_week_increment,
                                             treatmentarm=="DP 1 year"~four_week_increment-13,
                                             treatmentarm=="DP 2 years"~four_week_increment-26))




## summarise complicated malaria incidence, grouped by number of malaria episodes, parasitemia and age ####
all_malaria <- long_monthly_categories%>%
  group_by(id)%>%
  filter(!is.na(complicated)|!is.na(uncomplicated))%>%
  arrange(id, date)

all_para <- long_monthly_categories%>%
  group_by(id)%>%
  filter(!is.na(complicated) | !is.na(uncomplicated) | !is.na(asymptomatic))%>%
  arrange(id, date)


# complicated malaria ####
## individual level models ####
### all data ####
malaria_model <- lme4::glmer(bino_complicated ~ n_any_malaria + I(n_any_malaria^2) + (1|id), family = "binomial", data=all_malaria)
para_model <- lme4::glmer(bino_complicated ~ n_any_para + I(n_any_para^2) + (1|id),  family = "binomial", data=all_para)

age_malaria_model <- lme4::glmer(bino_complicated ~ flo_age_in_months + (1|id),  family = "binomial", data=all_malaria)
age_para_model <- lme4::glmer(bino_complicated ~ flo_age_in_months + (1|id),  family = "binomial", data=all_para)

### no DP only ####
no_dp_malaria <- filter(all_malaria, treatmentarm=="No DP")
no_dp_para <- filter(all_para, treatmentarm=="No DP")

no_dp_malaria_model <- lme4::glmer(bino_complicated ~ n_any_malaria + I(n_any_malaria^2) + (1|id), family = "binomial", data=no_dp_malaria)
no_dp_para_model <- lme4::glmer(bino_complicated ~ n_any_para + I(n_any_para^2) + (1|id),  family = "binomial", data=no_dp_para)

no_dp_age_malaria_model <- lme4::glmer(bino_complicated ~ flo_age_in_months + (1|id),  family = "binomial", data=no_dp_malaria)
no_dp_age_para_model <- lme4::glmer(bino_complicated ~ flo_age_in_months + (1|id),  family = "binomial", data=no_dp_para)

pardens_malaria_model <- lme4::lmer(log10(pardens) ~ n_any_malaria * treatmentarm+ (1|id),  data=all_malaria)
# coefs <- data.frame(coef(pardens_malaria_model)$id[1,])
# coefs$p.z <- 2 * pnorm(abs(coefs$n_any_malaria),lower.tail=FALSE)


## summary level models ####
all_malaria_summary <- all_malaria%>%
  group_by(bino_complicated, n_any_malaria)%>%
  summarise("n_cases"=n())%>%
  pivot_wider(names_from = bino_complicated, values_from = n_cases, names_prefix = "comp_")%>%
  mutate(total_infections=comp_1+comp_0,
         risk=comp_1/total_infections)%>%
  arrange(n_any_malaria)

all_para_summary <- all_para%>%
  group_by(bino_complicated, n_any_para)%>%
  summarise("n_cases"=n())%>%
  pivot_wider(names_from = bino_complicated, values_from = n_cases, names_prefix = "comp_")%>%
  mutate(total_infections=comp_1+comp_0,
         risk=comp_1/total_infections)%>%
  arrange(n_any_para)

all_malaria_age_summary <- all_para%>%
  # filter(treatmentarm!="DP 2 years")%>%
  group_by(bino_complicated, as.numeric(age_quarter))%>%
  summarise("n_cases"=n())%>%
  pivot_wider(names_from = bino_complicated, values_from = n_cases, names_prefix = "comp_")%>%
  mutate(total_infections=comp_1+comp_0,
         risk=comp_1/total_infections)%>%
  arrange(age_quarter)

### all data ####
malaria_summary_model <- glm(risk ~ n_any_malaria + I(n_any_malaria^2), family = "binomial", weights = total_infections, data=all_malaria_summary)
para_summary_model <- glm(risk ~ age_quarter ,  family = "binomial", weights = total_infections, data=all_malaria_age_summary)


#simple
# uncomplicated_model <- glm(risk ~ n_any_malaria + I(n_any_malaria^2), weights=total_infections, family = "binomial", data=all_malaria_summary)
all_malaria_summary_model <- glm(risk ~ n_any_malaria + I(n_any_malaria^2), family = "binomial", weights = total_infections, data=all_malaria_summary)
para_summary_model <- glm(risk ~ n_any_para + I(n_any_para^2) ,  family = "binomial", weights = total_infections, data=all_para_summary)

malaria_summary_model_prediction <- quadratic_glm_visualiser(model=malaria_summary_model,
                                                   x_name =  "n_any_malaria",
                                                   x_range = c(1,15))

(uncomplicated_comp_plot <- all_malaria_summary%>%
    ggplot(., aes(x=n_any_malaria, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_ribbon(data=malaria_summary_model_prediction$se, aes(x=n_any_malaria, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = malaria_summary_model_prediction$fun, colour="black")+
  geom_text(aes(y=0.17, label= paste0("frac(",comp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
  # scale_x_continuous(breaks = 1:50, limits=c(1,10))+
  scale_y_continuous(limits = c(0,0.2), labels = scales::label_percent())+
  # ggtitle("Placebo")+
  xlab("Order of Infection")+
  ylab("Risk of Complicated"))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/all_uncomplicated_comp_plot.png", uncomplicated_comp_plot, width=4, height=4, dpi=444)

para_summary_model_prediction <- quadratic_glm_visualiser(model=para_summary_model,
                                                   x_name =  "n_any_para",
                                                   x_range = c(1, 12))

(para_comp_plot <- ggplot(all_para_summary, aes(x=n_any_para, y=risk))+
    geom_point(color="darkred")+
    theme_minimal()+
    # geom_ribbon(data=para_summary_model_prediction$se, aes(x=n_any_para, ymin = exp(lci), ymax = exp(uci)),
    #             alpha = 0.2, inherit.aes = FALSE)+
    # geom_function(fun = para_summary_model_prediction$fun, colour="black")+
    geom_text(aes(y=0.17, label= paste0("frac(",comp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
    scale_x_continuous(breaks = 1:50, limits=c(1,12))+
    scale_y_continuous(limits = c(0,0.2), labels = scales::label_percent())+
    ggtitle("")+
    xlab("Order of Infection")+
    ylab("Risk of Complicated Malaria"))




(age_comp_plot <- ggplot(all_malaria_age_summary, aes(x=age_quarter, y=risk))+
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
    # facet_wrap(~treatmentarm)+
    xlab("Age")+
    ylab("Risk of Complicated Malaria")+
    theme(axis.text.x = element_text(angle=90, hjust=0.5, vjust=1)))


### no DP only ####
no_dp_malaria_summary <- all_malaria%>%
  filter(treatmentarm=="No DP")%>%
  group_by(bino_complicated, n_any_malaria, treatmentarm)%>%
  summarise("n_cases"=n())%>%
  pivot_wider(names_from = bino_complicated, values_from = n_cases, names_prefix = "comp_")%>%
  mutate(total_infections=comp_1+comp_0,
         risk=comp_1/total_infections)%>%
  arrange(n_any_malaria)

no_dp_malaria_summary_model <- glm(risk ~ n_any_malaria + I(n_any_malaria^2), family = "binomial", weights = total_infections, data=no_dp_malaria_summary)

no_dp_malaria_summary_model_prd <- quadratic_glm_visualiser(no_dp_malaria_summary_model, x_range=c(1,12), x_name="n_any_malaria")


(no_dp_uncomplicated_comp_plot <- no_dp_malaria_summary%>%
    ggplot(., aes(x=n_any_malaria, y=risk))+
    geom_point(color="darkred")+
    theme_minimal()+
    geom_ribbon(data=no_dp_malaria_summary_model_prd$se, aes(x=n_any_malaria, ymin = exp(lci), ymax = exp(uci)),
                alpha = 0.2, inherit.aes = FALSE)+
    geom_function(fun = no_dp_malaria_summary_model_prd$fun, colour="black")+
    geom_text(aes(y=0.17, label= paste0("frac(",comp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
    scale_x_continuous(breaks = 1:50, limits=c(1,10))+
    scale_y_continuous(limits = c(0,0.2), labels = scales::label_percent())+
    # ggtitle("Placebo")+
    xlab("Order of Infection")+
    ylab("Risk of Complicated"))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/no_dp_uncomplicated_comp_plot.png", no_dp_uncomplicated_comp_plot, width=4, height=4, dpi=444)

### DP 1 only ####
one_dp_malaria_summary <- all_malaria%>%
  filter(treatmentarm=="DP 1 year")%>%
  group_by(bino_complicated, n_any_malaria, treatmentarm)%>%
  summarise("n_cases"=n())%>%
  pivot_wider(names_from = bino_complicated, values_from = n_cases, names_prefix = "comp_")%>%
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
  pivot_wider(names_from = bino_complicated, values_from = n_cases, names_prefix = "comp_")%>%
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

# probability of symptoms ####
para_model <- lme4::glmer(bino_symp ~ n_any_para*treatmentarm + (1|id),  family = "binomial", data=all_para)

all_para_symp_summary <- all_para%>%
  filter(treatmentarm!= "DP 2 years")%>%
  group_by(bino_symp, n_any_para, treatmentarm)%>%
  summarise("n_cases"=n())%>%
  pivot_wider(names_from = bino_symp, values_from = n_cases, names_prefix = "symp_")%>%
  mutate(total_infections=symp_1+symp_0,
         risk=symp_1/total_infections)%>%
  arrange(n_any_para)

all_mala_symp_summary <- all_para%>%
  filter(treatmentarm!= "DP 2 years")%>%
  group_by(bino_symp, n_any_malaria, treatmentarm)%>%
  summarise("n_cases"=n())%>%
  pivot_wider(names_from = bino_symp, values_from = n_cases, names_prefix = "symp_")%>%
  mutate(total_infections=symp_1+symp_0,
         risk=symp_1/total_infections)%>%
  arrange(n_any_malaria)


(all_para_symp_summary_plot <- all_para_symp_summary%>%
    filter(treatmentarm=="No DP")%>%
    ggplot(., aes(x=n_any_para, y=risk))+
    theme_minimal()+
    geom_smooth(method = "lm", color="black")+
    geom_point(color="darkred")+
    geom_text(aes(y=0.8, label= paste0("frac(",symp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
    scale_x_continuous(breaks = 1:50, limits=c(1,18))+
    scale_y_continuous(labels = scales::label_percent())+
    ggtitle("")+
    xlab("")+
    ylab("P(symptoms | parsitemia)"))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/all_para_symp_summary_plot.png", all_para_symp_summary_plot, width=5, height=5, bg="white")


(all_mala_symp_summary <- ggplot(all_mala_symp_summary, aes(x=n_any_malaria, y=risk))+
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

# parasite density ####
pardens_plot <- all_malaria%>%
  filter(anyDP=="no")%>%
  ggplot(., aes(x=factor(n_any_malaria), y=pardens, fill=factor(n_any_malaria)))+
  scale_y_log10()+
  geom_boxplot(outliers = F)+
  stat_summary(geom="crossbar", colour="white", fun = "median", width=0.65, linewidth=0.19 )+
  theme_minimal()+
  ylab("parasites / μL")+
  xlab("")+
  scale_fill_manual(values=viridis::rocket(n=18))+
  # facet_wrap(~anyDP)+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/n_infection_age_load.png", pardens_plot, width=5, height=5, bg="white", dpi=444)



# heatmap library(ComplexHeatmap)

categories_matrix <- monthly_categories%>%
  filter(study=="micdrop", treatmentarm=="No DP")%>%
  mutate(mstatus_num=as.numeric(mstatus))%>%
  group_by(id, four_week_increment, mstatus_num)%>%
  slice_head(n = 1)%>%
  ungroup()%>%
  pivot_wider(names_from = four_week_increment,values_fill = 0, values_from = mstatus_num, id_cols = c("id"))%>%
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
