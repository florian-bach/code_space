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


side_by_side_theme <- theme(axis.text = element_text(size = 20),
                            axis.title = element_text(size = 22))


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
comp_model <- glm(risk~n_infection+I(n_infection^2), family = "binomial", weights = total_infections, data = complicated_df)
# weighted linear regression with quadratic term

#save model as function for plotting
comp_model_fun2 <- function(x){
  exp(comp_model2$coefficients[1])*
    exp(comp_model2$coefficients[2])^x*
    exp(comp_model2$coefficients[3])^x^2}



#calculate SE for plotting
prd <- data.frame(n_infection = seq(from = 1, to = 10, length.out = 100))
err <- predict(comp_model, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit



ggplot(complicated_df, aes(x=n_infection, y=risk))+
  geom_point(color="darkred")+
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
  geom_point(color="darkred", aes(alpha=complicatedmalaria))+
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
  geom_point(color="darkred")+
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
  mutate(age_at_first=age[n_infection==1])%>%
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


comp_model <- glm(risk~n_infection+I(n_infection^2), family = "binomial", weights = total_infections, data = complicated_df)
# weighted linear regression with quadratic term

#save model as function for plotting
comp_model_fun <- function(x){
  exp(comp_model$coefficients[1])*
    exp(comp_model$coefficients[2])^x*
    exp(comp_model$coefficients[3])^x^2}



#calculate SE for plotting

comp_risk_n_infection_plot <- ggplot(complicated_df, aes(x=n_infection, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_ribbon(data=model_visualiser(comp_model, "n_infection"), aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = comp_model_fun, colour="black")+
  ylab("Risk of Complicated Malaria")+
  geom_text(aes(y=0.19, label= paste0("frac(",complicated_episodes, ",", total_infections,")")),parse = TRUE, size=2.5)+
  scale_x_continuous(breaks = 1:10)+
  scale_y_continuous(limits = c(0,0.2), labels = scales::label_percent())+
  xlab("Order of Infection")

ggsave("~/postdoc/stanford/clinical_data/PROMOTE/figures/comp_risk_n_infection.png", comp_risk_n_infection_plot, height = 3, width=5, dpi=444, bg="white")




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
  geom_point(color="darkred")+
  theme_minimal()+
  # geom_ribbon(data=prd, aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
  #             alpha = 0.2, inherit.aes = FALSE)+
  # geom_function(fun = comp_model_fun, colour="black")+
  ylab("Risk of Complicated Episodes")+
  geom_smooth(method="glm", color="black")+
  # scale_x_continuous(breaks = 1:10)+
  #scale_y_continuous(limits = c(0,0.2))+
  xlab("Age (months)")



complicated_data %>%
  filter(id %in% kids_with_complicated_malaria$id)%>%
  ggplot(., aes(x=n_infection, y=parsdens, group=n_infection))+
  geom_boxplot(aes(fill=factor(complicatedmalaria)))+
  scale_fill_manual(values = c("uncomplicated"="#FC6A03",
                               "complicated"="darkred"))+
  theme_minimal()+
  theme(legend.title = element_blank())



# comp_model2 <- glm(risk~n_infection+I(n_infection^2), data = complicated_df, weights = total_infections, family="binomial")
# 
# comp_model_fun2 <- function(x){
#   exp(comp_model2$coefficients[1])*
#     exp(comp_model2$coefficients[2])^x*
#     exp(comp_model2$coefficients[3])^x^2}
# 

risk_by_n_infection_plot <- ggplot(complicated_df, aes(x=n_infection, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  side_by_side_theme+
  geom_ribbon(data=model_visualiser(comp_model2, "n_infection"), aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_text(aes(y=0.15, label= paste0("frac(",complicated_episodes, ",", total_infections,")")),parse = TRUE, vjust= -0.2, size=3.5)+
  geom_function(fun = comp_model_fun2, colour="black")+
  scale_x_continuous(breaks = seq(1, 10))+
  ylab("Risk of Complicated Episodes")+
  xlab("Order of Infection")



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


# visualise risk~n_infection bucketed by an age cutoff ####
# age_cutoff <- 2/3

for(i in seq(2/12, 1, by=1/12)){

age_cutoff <- i
 
#looks the same 
# old_complicated_data <- filter(complicated_data, age_at_first>age_cutoff)
old_complicated_data <- filter(complicated_data, age>age_cutoff)
old_complicated_df <- as.data.frame(table(old_complicated_data$n_infection, old_complicated_data$complicatedmalaria))
colnames(old_complicated_df) <- c("n_infection", "disease", "complicated_episodes")

total_old_infections <- old_complicated_data %>%
  group_by(n_infection) %>%
  count(name = "total_infections")

# only include complicated (not uncomplicated) cases, inlcude more summary stats
old_complicated_df <- subset(old_complicated_df, old_complicated_df$disease=="complicated")
old_complicated_df$total_infections <- total_old_infections$total_infections
old_complicated_df$risk <- old_complicated_df$complicated_episodes/old_complicated_df$total_infections
old_complicated_df$n_infection <- as.numeric(old_complicated_df$n_infection)


comp_model_old <- glm(risk~n_infection+I(n_infection^2), data = old_complicated_df, weights = total_old_infections$total_infections, family="binomial")

comp_model_fun_old <- function(x){
  exp(comp_model_old$coefficients[1])*
    exp(comp_model_old$coefficients[2])^x*
    exp(comp_model_old$coefficients[3])^x^2}


old_complicated_n_infection <- ggplot(old_complicated_df, aes(x=n_infection, y=risk))+
  geom_point(color="red")+
  theme_minimal()+
  side_by_side_theme+
  geom_ribbon(data=model_visualiser(comp_model_old, "n_infection"), aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = comp_model_fun_old, colour="black")+
  scale_x_continuous(breaks = seq(1, 10))+
  scale_y_continuous(limits=c(0, 0.225))+
  geom_text(aes(y=0.2, label= paste0("frac(",complicated_episodes, ",", total_infections,")")),parse = TRUE, vjust=-0.2, size=3.5)+
  ggtitle(paste("over", 12*age_cutoff, "months"))+
  ylab("Risk of Complicated Episodes")+
  xlab("Order of Infection")

# looks the same
# young_complicated_data <- filter(complicated_data, age_at_first<age_cutoff)
young_complicated_data <- filter(complicated_data, age<age_cutoff)
young_complicated_data <- filter(complicated_data, age<age_cutoff)
young_complicated_df <- as.data.frame(table(young_complicated_data$n_infection, young_complicated_data$complicatedmalaria))
colnames(young_complicated_df) <- c("n_infection", "disease", "complicated_episodes")

total_young_infections <- young_complicated_data %>%
  group_by(n_infection) %>%
  count(name = "total_infections")

# only include complicated (not uncomplicated) cases, inlcude more summary stats
young_complicated_df <- subset(young_complicated_df, young_complicated_df$disease=="complicated")
young_complicated_df <- young_complicated_df[1:nrow(total_young_infections),]
young_complicated_df$total_infections <- total_young_infections$total_infections
young_complicated_df$risk <- young_complicated_df$complicated_episodes/young_complicated_df$total_infections
young_complicated_df$n_infection <- as.numeric(young_complicated_df$n_infection)


comp_model_young <- glm(risk~n_infection+I(n_infection^2), data = young_complicated_df, weights = total_young_infections$total_infections, family="binomial")

comp_model_fun_young <- function(x){
  exp(comp_model_young$coefficients[1])*
    exp(comp_model_young$coefficients[2])^x*
    exp(comp_model_young$coefficients[3])^x^2}


young_complicated_n_infection <- ggplot(young_complicated_df, aes(x=n_infection, y=risk))+
  geom_point(color="red")+
  theme_minimal()+
  side_by_side_theme+
  geom_ribbon(data=model_visualiser(comp_model_young, "n_infection"), aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = comp_model_fun_young, colour="black")+
  geom_text(aes(y=0.2, label= paste0("frac(",complicated_episodes, ",", total_infections,")")),parse = TRUE, vjust= -0.2, size=3.5)+
  scale_x_continuous(breaks = seq(1, 10))+
  scale_y_continuous(limits=c(0, 0.225))+
  ggtitle(paste("under", 12*age_cutoff, "months"))+
  # geom_smooth(formula = y~x+I(x^2))+
  ylab("Risk of Complicated Episodes")+
  xlab("Order of Infection")

combo_plot <- young_complicated_n_infection+old_complicated_n_infection

file_name <- paste("~/postdoc/stanford/clinical_data/PROMOTE/figures/risk_n_infection_split_", i*12, "_months.png", sep="")


ggsave(filename = file_name, height=8, width=12, bg="white", dpi=500)

}

# sandbox / stuff that's cut out ####

# restricting data to only kids that will eventually get complicated disease ##

# restricted_complicated_data <- complicated_data %>%
#   filter(id %in% kids_with_complicated_malaria$id)
# 
# restricted_total_infections <- restricted_complicated_data %>%
#   group_by(n_infection) %>%
#   count(name = "total_infections")
# 
# 
# restricted_complicated_df <- as.data.frame(table(restricted_complicated_data$n_infection, restricted_complicated_data$complicatedmalaria))
# colnames(restricted_complicated_df) <- c("n_infection", "disease", "complicated_episodes")
# restricted_complicated_df <- subset(restricted_complicated_df, restricted_complicated_df$disease=="complicated")
# restricted_complicated_df$total_infections <- restricted_total_infections$total_infections
# restricted_complicated_df$risk <- restricted_complicated_df$complicated_episodes/restricted_complicated_df$total_infections
# restricted_complicated_df$n_infection <- as.numeric(restricted_complicated_df$n_infection)
# 
# 
# 
# 
# comp_model2a <- glm(risk~n_infection+I(n_infection^2), data = restricted_complicated_df, weights = restricted_complicated_df$total_infections, family="binomial", )
# 
# comp_model_fun2a <- function(x){
#   exp(comp_model2a$coefficients[1])*
#     exp(comp_model2a$coefficients[2])^x*
#     exp(comp_model2a$coefficients[3])^x^2}
# 
# 
# risk_by_n_infection_plot <- ggplot(restricted_complicated_df, aes(x=n_infection, y=risk))+
#   geom_point(color="darkred")+
#   geom_text(aes(y=0.15, label= paste0("frac(",complicated_episodes, ",", total_infections,")")),parse = TRUE, vjust= -0.2, size=3.5)+
#   theme_minimal()+
#   side_by_side_theme+
#   geom_ribbon(data=model_visualiser(comp_model2a, "n_infection"), aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
#               alpha = 0.4, inherit.aes = FALSE)+
#   geom_function(fun = comp_model_fun2a, colour="black")+
#   scale_x_continuous(breaks = seq(1, 8))+
#   ylab("Risk of Complicated Episodes")+
#   xlab("Order of Infection")}


# Goncalves=esque arrow plot of parasitaemia ####

kids_with_exactly_one_complicated_episode <- complicated_data %>%
  group_by(id) %>%
  filter(complicatedmalaria=="complicated")%>%
  count(name="number_of_complicated_episodes")%>%
  filter(number_of_complicated_episodes==1)





parsdens_arrow_data <- complicated_data %>%
  filter(id %in% kids_with_exactly_one_complicated_episode$id)%>%
  mutate(age_at_comp_episode=age[comp_num==1])%>%
  filter(age_at_comp_episode<=age)%>%
  group_by(id, complicatedmalaria) %>%
  slice_max(parsdens, n=1)%>%
  group_by(id)%>%
  filter(n()==2)%>%
  mutate(pars_dens_comp=parsdens[comp_num==1])%>%
  group_by(id) %>%
  mutate(up_down=if_else(max(parsdens)>pars_dens_comp, "up", "down"))





parsdens_up_arrow <- parsdens_arrow_data %>%
  filter(up_down=="up")%>%
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
  ggtitle("n = 19")+
  ylab("parasites /  μL")+
  xlab("age (weeks)")+
  scale_color_manual(values = c("uncomplicated"="#FC6A03",
                               "complicated"="darkred"))+
  theme_minimal()+
  theme(legend.position = "none")


ggsave("~/postdoc/stanford/clinical_data/PROMOTE/figures/parsdens_up_arrow_plot.png", parsdens_up_arrow, height = 3, width=4, dpi=444, bg="white")


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
