library(dplyr)
library(ggplot2)
library(tidyr)
library(ggside)
library(visreg)

promote_data <- haven::read_dta("~/postdoc/stanford/clinical_data/PROMOTE/BC-3 childs all visit database FINAL.dta")

# restrict data to only include incident malaria (i.e. diagnosis more than 14 days after last treatment)
malaria_only <- filter(promote_data, incidentmalaria==1)

# make a dataframe with only the kids that did get complicated malaria
kids_with_complicated_malaria <- unique(subset(promote_data, promote_data$complicatedmalaria==1, select = id))

# select relevant columns, recode disease columns to something understandable
smaller_data <- promote_data %>%
  # filter(id %in% kids_with_complicated_malaria$id ) %>%
  select(id, date, dob, age, uniqueid, incidentmalaria, complicatedmalaria, severe) %>%
  mutate(complicatedmalaria=ifelse(is.na(complicatedmalaria), "mild", "complicated")) %>%
  mutate(severe=ifelse(is.na(severe), "non_severe", "severe")) %>%
  mutate(incidentmalaria=ifelse(is.na(incidentmalaria), "not_incident", "incident")) 
  
  
# age_at_first_infection <- smaller_data %>%
#   filter(incidentmalaria=="incident")%>%
#   group_by(id)%>%
#   slice_min(n = 1, order_by = age)

# look at only incident episiodes; split data by individual, make column indicating the total number of infections for each kid;
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

# binomial regression
comp_model <- glm(risk~n_infection, family = "binomial", weights = total_infections, data = complicated_df)

# negative binomial is an alternative to poisson, not binomial.
# comp_model2 <- MASS::glm.nb(complicated_episodes~n_infection, weights = total_infections, data = complicated_df)
visreg::visreg(comp_model)
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
  scale_y_continuous(limits = c(0,0.2))+
  xlab("Order of Infection")
