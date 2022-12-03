library(dplyr)
library(ggplot2)
library(tidyr)
library(ggside)
library(visreg)
library(mediation)


model_visualiser <- function(model=NULL, xvar=NULL){
  model_data = model$data
  xvar_range = range(model_data[[xvar]])
  pre_df = data.frame(seq(from = xvar_range[1], to = xvar_range[2], length.out = 100))
  colnames(pre_df)=xvar
  error = predict(model, newdata = pre_df, se.fit = TRUE)
  pre_df$lci = error$fit - 1.96 * error$se.fit
  pre_df$fit = error$fit
  pre_df$uci = error$fit + 1.96 * error$se.fit
  return(pre_df)
}


# data generation ####
promote_data <- haven::read_dta("~/postdoc/stanford/clinical_data/PROMOTE/BC-3 childs all visit database FINAL.dta")

# restrict data to only include incident malaria (i.e. diagnosis more than 14 days after last treatment)
malaria_only <- filter(promote_data, incidentmalaria==1)

# make a dataframe with only the kids that did get complicated malaria
kids_with_complicated_malaria <- unique(subset(promote_data, promote_data$complicatedmalaria==1, select = id))

# select relevant columns, recode disease columns to something understandable
smaller_data <- promote_data %>%
  # filter(id %in% kids_with_complicated_malaria$id ) %>%
  select(id, date, dob, age, temp, hb, uniqueid, incidentmalaria, complicatedmalaria, severe, parsdens) %>%
  mutate(complicatedmalaria=ifelse(is.na(complicatedmalaria), "mild", "complicated")) %>%
  mutate(severe=ifelse(is.na(severe), "non_severe", "severe")) %>%
  mutate(incidentmalaria=ifelse(is.na(incidentmalaria), "not_incident", "incident")) 
  
# modelling aggreagated risk data ####

# age_at_first_infection <- smaller_data %>%
#   filter(incidentmalaria=="incident")%>%
#   group_by(id)%>%
#   slice_min(n = 1, order_by = age)

# look at only incident episodes; split data by individual, make column indicating the total number of infections for each kid;
# add column for order of infection, add single column for disease status
complicated_data <- smaller_data %>%
  filter(incidentmalaria=="incident")%>%
  mutate(id=factor(id))%>%
  group_by(id) %>%
  add_count(name="total_n_infection") %>%
  arrange(age) %>%
  mutate(n_infection = seq(1, max(total_n_infection))) %>%
  mutate(disease=ifelse(severe=="severe", "severe",
                        ifelse(complicatedmalaria=="complicated", "complicated", "mild")
    ))
  

total_infections <- complicated_data %>%
  group_by(n_infection) %>%
  count(name = "total_infections")

complicated_df <- as.data.frame(table(complicated_data$n_infection, complicated_data$complicatedmalaria))
colnames(complicated_df) <- c("n_infection", "disease", "complicated_episodes")
complicated_df <- subset(complicated_df, complicated_df$disease=="complicated")
complicated_df$total_infections <- total_infections$total_infections
complicated_df$risk <- complicated_df$complicated_episodes/complicated_df$total_infections
complicated_df$n_infection <- as.numeric(complicated_df$n_infection)

complicated_df$n_infection_squared <- complicated_df$n_infection^2
# binomial regression
comp_model <- glm(risk~n_infection, family = "binomial", weights = total_infections, data = complicated_df)
# weighted linear regression with quadratic term
comp_model2 <- lm(risk~n_infection+I(n_infection^2), data = complicated_df, weights = total_infections)
summary(comp_model2)

visreg::visreg(comp_model2, ylab="risk")

#save model as function for plotting
comp_model_fun <- function(x){exp(comp_model$coefficients[1])*exp(comp_model$coefficients[2])^x}


#calculate SE for plotting
prd <- data.frame(n_infection = seq(from = 1, to = 10, length.out = 100))
err <- predict(comp_model, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit



ggplot(complicated_df, aes(x=n_infection, y=risk))+
  geom_point(color="red")+
  theme_minimal()+
  geom_ribbon(data=prd, aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = comp_model_fun, colour="black")+
  ylab("Risk of Complicated Episodes")+
  scale_x_continuous(breaks = 1:10)+
  scale_y_continuous(limits = c(0,0.6))+
  xlab("Order of Infection")



observation_model <- lm(total_infections~n_infection+I(n_infection^2), data=complicated_df)
visreg::visreg(observation_model)

severe_df <- as.data.frame(table(complicated_data$n_infection, complicated_data$severe))
colnames(severe_df) <- c("n_infection", "disease", "severe_episodes")
severe_df <- subset(severe_df, severe_df$disease=="severe")
severe_df$total_infections <- total_infections$total_infections
severe_df$risk <- severe_df$severe_episodes/severe_df$total_infections
severe_df$n_infection <- as.numeric(severe_df$n_infection)

severe_model <- glm(risk~n_infection, family = "binomial", weights = total_infections, data = severe_df)

# individual level data ####

ggplot(complicated_data, aes(x=factor(n_infection), y=parsdens, fill=factor(n_infection)))+
  geom_violin()+
  geom_point(color="red", aes(alpha=complicatedmalaria))+
  theme_minimal()+
  ylab("Parasites /  μL")+
  scale_fill_manual(values = colorspace::sequential_hcl(10, palette = "Purple Yellow"))+
  scale_alpha_manual(values = c("mild"=0, "complicated"=1))+
  scale_y_log10()+
  theme(legend.position = "none")
  xlab("Order of Infection")
  
  
ggplot(complicated_data, aes(x=factor(n_infection), y=parsdens, fill=complicatedmalaria))+
    geom_boxplot()+
    theme_minimal()+
    ylab("Parasites /  μL")+
    scale_fill_manual(values = colorspace::sequential_hcl(10, palette = "Purple Yellow"))+
    scale_y_log10()+
    #theme(legend.position = "none")
    xlab("Order of Infection")
  
parasitaemia_model <- lm(log10(parsdens)~n_infection, data=complicated_data)

complicated_data$comp_num <- ifelse(complicated_data$complicatedmalaria=="complicated", 1, 0)
disease_model <-glm(comp_num~age, data=complicated_data, family = "binomial")

disease_model2 <-glm(comp_num~age+log10(parsdens), data=complicated_data, family = "binomial")



disease_model <-glm(comp_num~age, data=complicated_data, family = "binomial")


#  figure ####


disease_model <-glm(comp_num~age, data=complicated_data, family = "binomial")
disease_model_fun <- function(x){exp(disease_model$coefficients[1])*exp(disease_model$coefficients[2])^x}

ggplot(complicated_data, aes(x=age, y=comp_num))+
  geom_point(color="red")+
  theme_minimal()+
  geom_ribbon(data=model_visualiser(disease_model, "age"), aes(x=age, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = comp_model_fun, colour="black")+
  ylab("Risk of Complicated Episodes")+
  xlab("Age in years")




# streamlined code ####


promote_data <- haven::read_dta("~/postdoc/stanford/clinical_data/PROMOTE/BC-3 childs all visit database FINAL.dta")

# restrict data to only include incident malaria 
malaria_only <- filter(promote_data, incidentmalaria==1)

# make a dataframe with only the kids that did get complicated malaria
kids_with_complicated_malaria <- unique(subset(promote_data, promote_data$complicatedmalaria==1, select = id))

# select relevant columns, recode disease columns to something human-readable
smaller_data <- promote_data %>%
  dplyr::select(id, date, dob, age, temp, hb, uniqueid, incidentmalaria, complicatedmalaria, severe, parsdens) %>%
  mutate(complicatedmalaria=ifelse(is.na(complicatedmalaria), "uncomplicated", "complicated")) %>%
  mutate(severe=ifelse(is.na(severe), "non_severe", "severe")) %>%
  mutate(incidentmalaria=ifelse(is.na(incidentmalaria), "not_incident", "incident")) %>%
  mutate(age_cat=ifelse(age>0.5, "under_six_months", "over_six_months"))%>%
  mutate(age_cat_num=ifelse(age>0.5, 0, 1))%>%
  mutate(age_months= round(age*12, digits = 0))
  

# look at only incident episodes; split data by individual, make column indicating the total number of infections for each kid;
# add column for order of infection, add single column for disease status
complicated_data <- smaller_data %>%
  filter(incidentmalaria=="incident")%>%
  group_by(id) %>%
  add_count(name="total_n_infection") %>%
  arrange(age) %>%
  mutate(n_infection = seq(1, max(total_n_infection))) %>%
  mutate(disease=ifelse(severe=="severe", "severe",
                        ifelse(complicatedmalaria=="complicated", "complicated", "uncomplicated")),
         comp_num=ifelse(complicatedmalaria=="complicated", 1, 0)
  )

# calculate total number of observations for 1st, 2nd infection and so on
total_infections <- complicated_data %>%
  group_by(n_infection) %>%
  count(name = "total_infections")

# put together summary data frame 
complicated_df <- as.data.frame(table(complicated_data$n_infection, complicated_data$complicatedmalaria))
colnames(complicated_df) <- c("n_infection", "disease", "complicated_episodes")

# only include complicated (not uncomplicated) cases, inlcude more summary stats
complicated_df <- subset(complicated_df, complicated_df$disease=="complicated")
complicated_df$total_infections <- total_infections$total_infections
complicated_df$risk <- complicated_df$complicated_episodes/complicated_df$total_infections
complicated_df$n_infection <- as.numeric(complicated_df$n_infection)


complicated_data$n_infection2 <- complicated_data$n_infection^2

outcome_model <- lme4::glmer(comp_num~n_infection+n_infection2+age_cat_num+(1|id), data=complicated_data, family = "binomial")
age_mediator <- lme4::glmer(age_cat_num~n_infection+n_infection2+(1|id), data=complicated_data, family = "binomial")

med_out <- mediate(age_mediator, outcome_model, treat="age_cat_num", mediator = "n_infection", robustSE = TRUE, sims = 1000, group.out = "age_cat_num")


outcome_model <- lme4::glmer(comp_num~n_infection*+age_cat_num+(1|id), data=complicated_data, family = "binomial")
age_mediator <- lme4::glmer(age_cat_num~n_infection+n_infection2+(1|id), data=complicated_data, family = "binomial")

med_out <- mediate(age_mediator, outcome_model, treat="age_cat_num", mediator = "n_infection", sims = 1000, group.out = NULL)


med_out2 <- mediate(age_mediator, outcome_model, treat="n_infection", mediator = "age", robustSE = TRUE, sims = 1000, group.out = NULL, control.value = 1, treat.value = 2)
med_out3 <- mediate(age_mediator, outcome_model, treat="n_infection", mediator = "age", robustSE = TRUE, sims = 1000, group.out = NULL, control.value = 1, treat.value = 3)
med_out4 <- mediate(age_mediator, outcome_model, treat="n_infection", mediator = "age", robustSE = TRUE, sims = 1000, group.out = NULL, control.value = 1, treat.value = 4)
med_out5 <- mediate(age_mediator, outcome_model, treat="n_infection", mediator = "age", robustSE = TRUE, sims = 1000, group.out = NULL, control.value = 1, treat.value = 5)
med_out6 <- mediate(age_mediator, outcome_model, treat="n_infection", mediator = "age", robustSE = TRUE, sims = 1000, group.out = NULL, control.value = 1, treat.value = 6)



complicated_data %>%
  filter(n_infection<7)%>%
  ggplot(., aes(y=age, x=complicatedmalaria, group=factor(complicatedmalaria, levels = c("uncomplicated", "complicated")), fill=complicatedmalaria))+
  facet_wrap(~n_infection)+
  geom_boxplot()+
  scale_fill_manual(values = c("uncomplicated"="#FC6A03",
                               "complicated"="darkred"))+
  theme_minimal()+
  theme(legend.title = element_blank(),
        axis.text.x = element_blank())

# calculate total number of observations for 1st, 2nd infection and so on
total_infections_by_age <- complicated_data %>%
  group_by(age_months) %>%
  count(name = "total_infections")


complicated_df_age <- as.data.frame(table(complicated_data$age_months, complicated_data$complicatedmalaria))
colnames(complicated_df_age) <- c("age_months", "disease", "complicated_episodes")

# only include complicated (not uncomplicated) cases, inlcude more summary stats
complicated_df_age <- subset(complicated_df_age, complicated_df_age$disease=="complicated")
complicated_df_age$total_infections <- total_infections_by_age$total_infections
complicated_df_age$risk <- complicated_df_age$complicated_episodes/complicated_df_age$total_infections



ggplot(complicated_df_age, aes(x=as.numeric(age_months), y=risk))+
  geom_point(color="red")+
  theme_minimal()+
  # geom_ribbon(data=prd, aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
  #             alpha = 0.2, inherit.aes = FALSE)+
  # geom_function(fun = comp_model_fun, colour="black")+
  ylab("Risk of Complicated Episodes")+
  geom_smooth(method="glm", color="black")+
  # scale_x_continuous(breaks = 1:10)+
  #scale_y_continuous(limits = c(0,0.2))+
  xlab("Age (months)")




complicated_data$is_first <- ifelse(complicated_data$n_infection==1, 1, 0)

age_at_fist_glm <- complicated_data%>%
  lme4::glmer(comp_num ~ n_infection+I(n_infection^2)+age_cat+is_first+log(parsdens)+(1|id), family = "binomial", data=.)


#Intercept) n_infection I(n_infection^2) age_catover_six_months log(parsdens)
age_cat_contrast <- t(matrix(c(0,0,0,1,1,0)))
age_cat_test <- multcomp::glht(age_at_fist_glm, linfct = age_cat_contrast)
summary(age_cat_test)


age_at_fist_glm2 <- complicated_data%>%
  lme4::glmer(comp_num ~ age_cat+is_first+log(parsdens)+(1|id), family = "binomial", data=.)

#(Intercept) age_catover_six_months   is_first log(parsdens)
age_cat_contrast2 <- t(matrix(c(0,1,0,0)))
age_cat_test2 <- multcomp::glht(age_at_fist_glm2, linfct = age_cat_contrast2)
summary(age_cat_test2)


visreg(age_at_fist_glm, xvar="n_infection", by="age_cat", trans=exp)
