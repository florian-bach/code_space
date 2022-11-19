library(dplyr)
library(ggplot2)
library(tidyr)
library(ggside)
library(visreg)

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

ggplot(complicated_df, aes(x=n_infection, y=total_infections))+
  geom_point(color="red")+
  theme_minimal()+
  geom_smooth(method = glm, formula = y~x+I(x^2))+
  ylab("Risk of Complicated Episodes")+
  scale_x_continuous(breaks = 1:10)+
  #scale_y_continuous(limits = c(0,0.6))+
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





comp_model_fun2 <- function(x){
  exp(comp_model2$coefficients[1])*
    exp(comp_model2$coefficients[2])^x*
    exp(comp_model2$coefficients[3])^x^2}


comp_model_fun2(5)
