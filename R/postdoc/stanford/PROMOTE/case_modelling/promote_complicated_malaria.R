library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(ggside)
library(visreg)
library(mediation)
library(patchwork)


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


age_cat_palette <- c(rgb(5,50,80, maxColorValue = 255),
                     rgb(250, 100, 0, maxColorValue = 255))



# data generation ####
promote_data <- haven::read_dta("~/postdoc/stanford/clinical_data/PROMOTE/BC-3 childs all visit database FINAL.dta")

# restrict data to only include incident malaria (i.e. diagnosis more than 14 days after last treatment)
malaria_only <- filter(promote_data, mstatus %in% c(1,2))

# make a dataframe with only the kids that did get complicated malaria
kids_with_complicated_malaria <- unique(subset(promote_data, promote_data$complicatedmalaria==1, select = id))

# select relevant columns, recode disease columns to something understandable
smaller_data <- promote_data %>%
  # filter(id %in% kids_with_complicated_malaria$id ) %>%
  dplyr::select(id, date, dob, age, temp, hb, uniqueid, mstatus, incidentmalaria, complicatedmalaria, severe, parsdens) %>%
  mutate(complicatedmalaria=ifelse(is.na(complicatedmalaria), "mild", "complicated")) %>%
  mutate(severe=ifelse(is.na(severe), "non_severe", "severe")) %>%
  mutate(incidentmalaria=ifelse(is.na(incidentmalaria), "not_incident", "incident")) 
  
# modelling aggreagated risk data ####

# age_at_first_infection <- smaller_data %>%
#   filter(incidentmalaria=="incident")%>%
#   group_by(id)%>%
#   slice_min(n = 1, order_by = age)

# look at only incident episodes; split data by individual, make column indicating the total number of infections for each kid;
# add column for order of infection, add single column for disease status (comp and severe into one)
complicated_data <- smaller_data %>%
  filter(mstatus %in% c(1,2))%>%
  mutate(id=factor(id))%>%
  group_by(id) %>%
  add_count(name="total_n_infection") %>%
  arrange(age) %>%
  mutate(n_infection = seq(1, max(total_n_infection)))%>%
  mutate(disease=ifelse(severe=="severe", "complicated",
                        ifelse(complicatedmalaria=="complicated", "complicated", "uncomplicated")
    ))%>%
  # mutate(dummy=as.numeric(if_else(any(incidentmalaria=="not_incident"), which(incidentmalaria=="not_incident"), 100)))%>%
  group_by(id) %>%
  #there are 2 complicated episodes that are non-incident, and 11 total. we keep the later one for all, allowing resolution of complicated episodes.
  mutate(when_not_incident=ifelse(any(incidentmalaria=="not_incident"), which(incidentmalaria=="not_incident"), 200))%>%
  mutate(n_infection=ifelse(n_infection>=when_not_incident, n_infection-1, n_infection))%>%
  group_by(id, n_infection)%>%
  slice_max(order_by = date, n=1, with_ties = FALSE)



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
  geom_boxplot()+
  geom_point(color="darkred", aes(alpha=factor(disease)))+
  theme_minimal()+
  ylab("Parasites /  μL")+
  scale_fill_manual(values = colorspace::sequential_hcl(10, palette = "Purple Yellow"))+
  scale_alpha_manual(values = c("uncomplicated"=0.1, "complicated"=1))+
  scale_y_log10()+
  theme(legend.position = "none")
  xlab("Order of Infection")
  
  
  
comp_counts <- complicated_data%>%
  group_by(disease, n_infection)%>%
  count()

ggplot(complicated_data, aes(x=factor(n_infection), y=parsdens, fill=complicatedmalaria))+
    geom_boxplot()+
    theme_minimal()+
    geom_text(data = comp_counts, aes(x=factor(n_infection), y=5*10^5, label=n), inherit.aes = FALSE,
              #position = position_dodge(width = 0.75)
              )+
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
malaria_only <- filter(promote_data, mstatus %in% c(1,2))

# make a dataframe with only the kids that did get complicated malaria
kids_with_complicated_malaria <- unique(subset(promote_data, promote_data$complicatedmalaria==1, select = id))

# select relevant columns, recode disease columns to something human-readable
smaller_data <- promote_data %>%
  dplyr::select(id, date, dob, age, temp, hb, uniqueid, mstatus, incidentmalaria, complicatedmalaria, severe, parsdens, tempgrade) %>%
  mutate(complicatedmalaria=ifelse(is.na(complicatedmalaria), "uncomplicated", "complicated")) %>%
  mutate(severe=ifelse(is.na(severe), "non_severe", "severe")) %>%
  mutate(incidentmalaria=ifelse(is.na(incidentmalaria), "not_incident", "incident")) %>%
  mutate(age_cat=ifelse(age>0.5, "under_six_months", "over_six_months"))%>%
  mutate(age_cat_num=ifelse(age>0.5, 0, 1))%>%
  mutate(age_months= round(age*12, digits = 0))
  

# look at only incident episodes; split data by individual, make column indicating the total number of infections for each kid;
# add column for order of infection, add single column for disease status
complicated_data <- smaller_data %>%
  filter(mstatus %in% c(1,2))%>%
  group_by(id) %>%
  add_count(name="total_n_infection") %>%
  arrange(age) %>%
  mutate(n_infection = seq(1, max(total_n_infection))) %>%
  mutate(age_at_first=age[n_infection==1])%>%
  mutate(disease=ifelse(severe=="severe", "severe",
                        ifelse(complicatedmalaria=="complicated", "complicated", "uncomplicated")),
         comp_num=ifelse(complicatedmalaria=="complicated", 1, 0)
  )%>%
  mutate(first_inf_before_two = if_else(age_at_first<2/12, "before 4 months", "after 4 months"))%>%
  mutate(first_inf_before_three = if_else(age_at_first<3/12, "before 4 months", "after 4 months"))%>%
  mutate(first_inf_before_four = if_else(age_at_first<4/12, "before 4 months", "after 4 months"))%>%
  mutate(first_inf_before_six = if_else(age_at_first<0.5, "before 6 months", "after 6 months"))%>%
  mutate(when_not_incident=ifelse(any(incidentmalaria=="not_incident"), which(incidentmalaria=="not_incident"), 200))%>%
  mutate(n_infection=ifelse(n_infection>=when_not_incident, n_infection-1, n_infection))%>%
  group_by(id, n_infection)%>%
  slice_max(order_by = date, n=1, with_ties = FALSE)
  

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
  ylab("Risk of Complicated Malaria,\n Given an Infection\n")+
  #geom_text(aes(y=0.19, label= paste0("frac(",complicated_episodes, ",", total_infections,")")),parse = TRUE, size=2.5)+
  scale_x_continuous(breaks = 1:10)+
  scale_y_continuous(limits = c(0,0.2), labels = scales::label_percent())+
  xlab("Order of Infection")

ggsave("~/postdoc/stanford/clinical_data/PROMOTE/figures/comp_risk_n_infection.png", comp_risk_n_infection_plot, height = 3, width=5, dpi=444, bg="white")

 
promote_comp_cases_plot <- complicated_data %>%
  arrange(desc(factor(complicatedmalaria)))%>%
  ggplot(., aes(x=n_infection, y=parsdens))+
  geom_point(aes(color=factor(complicatedmalaria, levels = c("complicated", "uncomplicated"))))+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, color="darkred")+
  # geom_line(alpha=0.3, aes(group=id))+
  # geom_boxplot(aes(group=factor(AGE)))+
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", width = 0.8, color="white")+
  ggtitle("")+
  ylab("qPCR Parasite Density")+
  xlab("n_infection")+
  scale_y_log10(limits=c(0.9, 10^6), breaks=c(10^0, 10^2, 10^4, 10^6), labels=scales::label_number())+
  scale_x_continuous(limits = c(0,10), breaks = seq(1,10))+
  annotation_logticks(sides = "l", )+
  theme_minimal()+
  # facet_wrap(~mstatus)+
  theme(legend.title = element_blank())+
  scale_color_manual(values = c(comp_pal, "quinine"="purple"))




promote_comp <-  promote_comp_cases_plot + comp_risk_n_infection_plot + plot_annotation(
  title = 'Complicated Disease in PROMOTE',
  subtitle = '')+
  theme(plot.tag = element_text(size = 10, hjust = 0, vjust = 0))

ggsave("~/postdoc/stanford/clinical_data/PROMOTE/figures/promote_comp.png", promote_comp, width=8, height=4, bg="white", dpi=444)





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
old_complicated_data <- filter(complicated_data, age_at_first>age_cutoff)

old_complicated_df <- complicated_data %>%
  group_by(n_infection, complicatedmalaria)%>%
  summarise("count"=n())
  

total_old_infections <- old_complicated_data %>%
  group_by(n_infection) %>%
  count(name = "total_infections")

# only include complicated (not uncomplicated) cases, inlcude more summary stats
old_complicated_df <- subset(old_complicated_df, old_complicated_df$complicatedmalaria=="complicated")
old_complicated_df$total_infections <- total_old_infections$total_infections
old_complicated_df$risk <- old_complicated_df$complicated_episodes/old_complicated_df$total_infections
old_complicated_df$n_infection <- as.numeric(old_complicated_df$n_infection)


comp_model_old <- glm(risk~n_infection+I(n_infection^2), data = old_complicated_df, weights = total_old_infections$total_infections, family="binomial")

comp_model_fun_old <- function(x){
  exp(comp_model_old$coefficients[1])*
    exp(comp_model_old$coefficients[2])^x*
    exp(comp_model_old$coefficients[3])^x^2}


old_complicated_n_infection <- ggplot(old_complicated_df, aes(x=n_infection, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  side_by_side_theme+
  geom_ribbon(data=model_visualiser(comp_model_old, "n_infection"), aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = comp_model_fun_old, colour="black")+
  scale_x_continuous(breaks = seq(1, 10))+
  scale_y_continuous(limits=c(0, 0.225), labels=scales::label_percent())+
  geom_text(aes(y=0.2, label= paste0("frac(",complicated_episodes, ",", total_infections,")")),parse = TRUE, vjust=-0.2, size=3.5)+
  ggtitle(paste("first infection over", 12*age_cutoff, "months"))+
  ylab("Risk of Complicated Episodes")+
  xlab("Order of Infection")

# looks the same
# young_complicated_data <- filter(complicated_data, age_at_first<age_cutoff)
young_complicated_data <- filter(complicated_data, age<age_cutoff)

young_complicated_df <- complicated_data %>%
  group_by(n_infection, complicatedmalaria)%>%
  summarise("count"=n())

total_young_infections <- young_complicated_data %>%
  group_by(n_infection) %>%
  count(name = "total_infections")

# only include complicated (not uncomplicated) cases, inlcude more summary stats
young_complicated_df <- subset(young_complicated_df, young_complicated_df$complicatedmalaria=="complicated")
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
  geom_point(color="darkred")+
  theme_minimal()+
  side_by_side_theme+
  geom_ribbon(data=model_visualiser(comp_model_young, "n_infection"), aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = comp_model_fun_young, colour="black")+
  geom_text(aes(y=0.2, label= paste0("frac(",complicated_episodes, ",", total_infections,")")),parse = TRUE, vjust= -0.2, size=3.5)+
  scale_x_continuous(breaks = seq(1, 10))+
  scale_y_continuous(limits=c(0, 0.225), labels=scales::label_percent())+
  ggtitle(paste("first infection under", 12*age_cutoff, "months"))+
  # geom_smooth(formula = y~x+I(x^2))+
  ylab("Risk of Complicated Episodes")+
  xlab("Order of Infection")

combo_plot <- young_complicated_n_infection+old_complicated_n_infection

# file_name <- paste("~/postdoc/stanford/clinical_data/PROMOTE/figures/risk_n_infection_split_", i*12, "_months.png", sep="")
file_name <- paste("~/postdoc/stanford/clinical_data/PROMOTE/figures/risk_n_infection_first_infection_split_", i*12, "_months.png", sep="")


ggsave(filename = file_name, height=8, width=12, bg="white", dpi=500)

}

# tempgrade distribution ####

tempgrade_plot <- complicated_data %>%
  group_by(n_infection, tempgrade)%>%
  summarise("n_grade"=n())%>%
  ungroup()%>%
  group_by(n_infection)%>%
  add_count(name="n_count", wt=n_grade)%>%
  mutate("temp_grade_perc"=n_grade/n_count)%>%
  ggplot(aes(x=factor(n_infection), y=temp_grade_perc))+
  geom_bar(stat = "identity", aes(fill=factor(tempgrade)))+
  scale_y_continuous(labels = scales::label_percent())+
  scale_fill_manual(values=viridis::inferno(n=5))+
  ggtitle("fever distribution stable\nin the first year of life")+
  theme(axis.title = element_blank(), 
        legend.title = element_blank())
  

ggsave("~/postdoc/stanford/clinical_data/PROMOTE/figures/tempgrade_plot.png", width=4, height=3, bg='white')

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




# calculate confidence interval for complicated malaria risk based on whether first infections happen
# before or after 2 months

age_first_nest2 <- complicated_data %>%
  group_by(first_inf_before_four) %>%
  nest() %>%
  mutate(confy_model=purrr::map(data, ~glm(comp_num~1, data=., family="binomial")),
         coef=purrr::map(confy_model, ~exp(coef(.))),
         upper=purrr::map(confy_model, ~exp(confint(.))[2]),
         lower=purrr::map(confy_model, ~exp(confint(.))[1]))%>%
  unnest(data)%>%
  dplyr::select(n_infection, age_at_first, coef, upper, lower, first_inf_before_four)%>%
  unnest()


# first_inf_model <- lme4::glmer(comp_num~first_inf_before_four+(1|id), data=complicated_data, family = "binomial")
# first_inf_model <- glm(comp_num~first_inf_before_four, data=complicated_data, family = "binomial")


overall_risk_age_at_first_confint4<- ggplot(age_first_nest2, aes(x = first_inf_before_four, y = coef, color=first_inf_before_four))+
  geom_point(aes(x = first_inf_before_four, y = coef)) + 
  geom_linerange(aes(x = first_inf_before_four,  ymin = lower, ymax = upper), linewidth = 1)+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 0.25), labels = scales::label_percent())+
  xlab("Order of Infection")+
  ylab("Risk of Complicated Malaria")+
  scale_color_manual(values=age_cat_palette)+
  theme(axis.text.x = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/clinical_data/PROMOTE/figures/overall_risk_age_2_at_first_confint.png", overall_risk_age_at_first_confint4, height = 3, width=6, dpi=444, bg="white")




# calculate confidence interval for complicated malaria risk based on whether first infections happen before or after 4
# months

age_first_nest4 <- complicated_data %>%
  group_by(first_inf_before_four) %>%
  nest() %>%
  mutate(confy_model=purrr::map(data, ~glm(comp_num~1, data=., family="binomial")),
         coef=purrr::map(confy_model, ~exp(coef(.))),
         upper=purrr::map(confy_model, ~exp(confint(.))[2]),
         lower=purrr::map(confy_model, ~exp(confint(.))[1]))%>%
  unnest(data)%>%
  dplyr::select(n_infection, age_at_first, coef, upper, lower, first_inf_before_four)%>%
  unnest()


# first_inf_model <- lme4::glmer(comp_num~first_inf_before_four+(1|id), data=complicated_data, family = "binomial")
# first_inf_model <- glm(comp_num~first_inf_before_four, data=complicated_data, family = "binomial")


overall_risk_age_at_first_confint4<- ggplot(age_first_nest4, aes(x = first_inf_before_four, y = coef, color=first_inf_before_four))+
  geom_point(aes(x = first_inf_before_four, y = coef)) + 
  geom_linerange(aes(x = first_inf_before_four,  ymin = lower, ymax = upper), linewidth = 1)+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 0.25), labels = scales::label_percent())+
  xlab("Order of Infection")+
  ylab("Risk of Complicated Malaria")+
  scale_color_manual(values=age_cat_palette)+
  theme(axis.text.x = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/clinical_data/PROMOTE/figures/overall_risk_age_4_at_first_confint.png", overall_risk_age_at_first_confint4, height = 3, width=6, dpi=444, bg="white")







age_n_inf_nest4 <- complicated_data %>%
  filter(n_infection<7)%>%
  group_by(n_infection, first_inf_before_four) %>%
  nest() %>%
  mutate(confy_model=purrr::map(data, ~glm(comp_num~1, data=., family="binomial")),
         coef=purrr::map(confy_model, ~exp(coef(.))),
         upper=purrr::map(confy_model, ~exp(confint(.))[2]),
         lower=purrr::map(confy_model, ~exp(confint(.))[1]))%>%
  unnest(data)%>%
  dplyr::select(n_infection, age_at_first, coef, upper, lower, first_inf_before_four)%>%
  unnest()


risk_age_at_first_confint4 <- ggplot(age_n_inf_nest4, aes(x = n_infection, y = coef, color=first_inf_before_four))+
  geom_point(position = position_dodge(width=1)) + 
  geom_linerange(aes(x = n_infection,  ymin = lower, ymax = upper), linewidth = 1, position = position_dodge(width=1))+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 0.42), labels = scales::label_percent())+
  scale_x_continuous(breaks = seq(1, 6))+
  xlab("Order of Infection")+
  ylab("Risk of Complicated Malaria")+
  scale_color_manual(values=age_cat_palette)+
  theme(legend.title = element_blank())


ggsave("~/postdoc/stanford/clinical_data/PROMOTE/figures/risk_age_4_at_first_confint.png", risk_age_at_first_confint4, height = 3, width=6, dpi=444, bg="white")



# do it again but split data at first infection by 6 month mark



age_first_nest6 <- complicated_data %>%
  group_by(first_inf_before_six) %>%
  nest() %>%
  mutate(confy_model=purrr::map(data, ~glm(comp_num~1, data=., family="binomial")),
         coef=purrr::map(confy_model, ~exp(coef(.))),
         upper=purrr::map(confy_model, ~exp(confint(.))[2]),
         lower=purrr::map(confy_model, ~exp(confint(.))[1]))%>%
  unnest(data)%>%
  dplyr::select(n_infection, age_at_first, coef, upper, lower, first_inf_before_six)%>%
  unnest()


# first_inf_model <- lme4::glmer(comp_num~first_inf_before_six+(1|id), data=complicated_data, family = "binomial")
# first_inf_model <- glm(comp_num~first_inf_before_six, data=complicated_data, family = "binomial")


overall_risk_age_at_first_confint6 <- ggplot(age_first_nest6, aes(x = first_inf_before_six, y = coef, color=first_inf_before_six))+
  geom_point(aes(x = first_inf_before_six, y = coef)) + 
  geom_linerange(aes(x = first_inf_before_six,  ymin = lower, ymax = upper), linewidth = 1)+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 0.25), labels = scales::label_percent())+
  xlab("Order of Infection")+
  ylab("Risk of Complicated Malaria")+
  scale_color_manual(values=age_cat_palette)+
  theme(axis.text.x = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/clinical_data/PROMOTE/figures/overall_risk_age_6_at_first_confint.png", overall_risk_age_at_first_confint6, height = 3, width=4, dpi=444, bg="white")



age_n_inf_nest6 <- complicated_data %>%
  filter(n_infection<7)%>%
  mutate(first_inf_before_six = if_else(age_at_first<0.5, "before 6 months", "after 6 months"))%>%
  group_by(n_infection, first_inf_before_six) %>%
  nest() %>%
  mutate(confy_model=purrr::map(data, ~glm(comp_num~1, data=., family="binomial")),
         coef=purrr::map(confy_model, ~exp(coef(.))),
         upper=purrr::map(confy_model, ~exp(confint(.))[2]),
         lower=purrr::map(confy_model, ~exp(confint(.))[1]))%>%
  unnest(data)%>%
  dplyr::select(n_infection, age_at_first, coef, upper, lower, first_inf_before_six)%>%
  unnest()


risk_age_at_first_confint6 <- ggplot(age_n_inf_nest6, aes(x = n_infection, y = coef, color=first_inf_before_six))+
  geom_point(position = position_dodge(width=1)) + 
  geom_linerange(aes(x = n_infection,  ymin = lower, ymax = upper), linewidth = 1, position = position_dodge(width=1))+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 0.42), labels = scales::label_percent())+
  scale_x_continuous(breaks = seq(1, 6))+
  xlab("Order of Infection")+
  ylab("Risk of Complicated Malaria")+
  scale_color_manual(values=age_cat_palette)+
  theme(legend.title = element_blank())


ggsave("~/postdoc/stanford/clinical_data/PROMOTE/figures/risk_age_6_at_first_confint.png", risk_age_at_first_confint6, height = 3, width=6, dpi=444, bg="white")



p1 <- promote_data%>%
  filter(complicatedmalaria==1)%>%
  ggplot(., aes(x=temp, y=parsdens))+
  geom_point(color="red")+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = "none")

p1m <- ggExtra::ggMarginal(p1)


p2 <- promote_data%>%
  filter(is.na(complicatedmalaria))%>%
  ggplot(., aes(x=temp, y=parsdens))+
  geom_point(color="darkgrey")+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = "none")

p2m <- ggExtra::ggMarginal(p2, type = "densigram", groupColour = )

cowplot::plot_grid(p1m, p2m)



p <- promote_data%>%
  #filter(complicatedmalaria==1)%>%
  ggplot(., aes(x=temp, y=parsdens, color=factor(complicatedmalaria)))+
  geom_point()+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = "none")

ggExtra::ggMarginal(p, groupColour = TRUE, type = "boxplot", )




#sandbox ####

# complicated_data%>%
#   arrange(desc(disease))%>%
#   ggplot(aes(x=factor(n_infection), y=parsdens, fill=factor(n_infection)))+
#   geom_boxplot()+
#   # geom_point(color="darkred", aes(alpha=factor(disease)))+
#   geom_point(aes(color=factor(disease)))+
#   theme_minimal()+
#   ylab("Parasites /  μL")+
#   scale_fill_manual(values = colorspace::sequential_hcl(10, palette = "Purple Yellow"))+
#   scale_color_manual(values = c("darkred", "orange"))+
#   scale_alpha_manual(values = c("uncomplicated"=0.1, "complicated"=1))+
#   scale_y_log10()+
#   theme(legend.position = "none")
# xlab("Order of Infection")
# 
# complicated_data %>%
#   group_by(id)%>%
#   mutate("high_before_six"=if_else(any(parsdens>10000) & age>=0.5, "yes", "no"))%>%
#   arrange(desc(disease))%>%
#   ggplot(aes(x=factor(n_infection), y=parsdens, fill=factor(n_infection)))+
#   geom_boxplot()+
#   # geom_point(color="darkred", aes(alpha=factor(disease)))+
#   geom_point(aes(color=factor(disease)))+
#   facet_wrap(~high_before_six)+
#   theme_minimal()+
#   ylab("Parasites /  μL")+
#   scale_fill_manual(values = colorspace::sequential_hcl(10, palette = "Purple Yellow"))+
#   scale_color_manual(values = c("darkred", "orange"))+
#   scale_alpha_manual(values = c("uncomplicated"=0.1, "complicated"=1))+
#   scale_y_log10()+
#   theme(legend.position = "none")
# xlab("Order of Infection")
# 
# 
# complicated_data <- complicated_data %>%
#   group_by(id)%>%
#   mutate("high_before_six"=if_else(any(parsdens>10000) & age <= 0.5, "yes", "no"))
# 
# df <- complicated_data %>%
#   group_by(high_before_six, n_infection)%>%
#   count(disease)%>%
#   pivot_wider(values_from = n, names_from = disease)%>%
#   mutate(complicated=replace_na(complicated, 0))%>%
#   mutate(risk=complicated/(uncomplicated+complicated))
# 
# 
# 
# ggplot(df, aes(x=n_infection, y=risk, fill=high_before_six))+
#   geom_bar(stat="identity")+
#   geom_text(aes(label= paste0("frac(",complicated, ",", uncomplicated+complicated,")")),parse = TRUE, vjust= -0.2, size=3.5)+
#   facet_wrap(~high_before_six)+
#   scale_x_continuous(breaks = seq(1, 7), limits=c(0, 7))+
#   scale_y_continuous(labels=scales::label_percent(), limits = c(0,0.15))+
#   scale_fill_manual(values=c("magenta4", "darkred"))+
#   theme_minimal()
# 


# jason models ####


firsts <- complicated_data %>%
  # filter(n_infection <=3)%>%
  mutate(binomial=if_else(complicatedmalaria=="complicated", 1, 0))
# 
# arrange(desc(disease))%>%
# ggplot(., aes(x=age, y=disease, color=disease))+
# geom_point()+
# theme_minimal()+
# theme()

model <- glm(binomial~age, data=firsts, family="binomial")
summary(model)


complicated_data %>%
  filter(n_infection ==2)%>%
  # mutate(binomial=if_else(complicatedmalaria=="complicated", 1, 0))%>%
  arrange(desc(disease))%>%
  ggplot(., aes(x=age, y=disease, color=disease))+
  geom_point()+
  theme_minimal()+
  theme()



complicated_df <- as.data.frame(table(firsts$age_months, firsts$complicatedmalaria))
colnames(complicated_df) <- c("age_months", "disease", "complicated_episodes")

# only include complicated (not uncomplicated) cases, inlcude more summary stats
complicated_df <- subset(complicated_df, complicated_df$disease=="complicated")
complicated_df$total_infections <- total_infections_by_age$total_infections
complicated_df$risk <- complicated_df$complicated_episodes/complicated_df$total_infections



# comp_risk_n_infection_plot <-
ggplot(complicated_df, aes(x=age_months, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_smooth()+
  ylab("Risk of Complicated Malaria\n")+
  ggtitle("in the first 3 infections of life")+
  #geom_text(aes(y=0.19, label= paste0("frac(",complicated_episodes, ",", total_infections,")")),parse = TRUE, size=2.5)+
  # scale_x_continuous(breaks = 1:12)+
  # scale_y_continuous(limits = c(0,0.2), labels = scales::label_percent())+
  xlab("Age in Months")

