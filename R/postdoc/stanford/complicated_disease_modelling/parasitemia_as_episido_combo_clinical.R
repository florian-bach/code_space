# preamble ####
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(ggside)
library(visreg)
# library(mediation)
# library(patchwork)


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

# promote ####

promote <- haven::read_dta("~/postdoc/stanford/clinical_data/PROMOTE/BC-3 childs all visit database FINAL.dta")
bc3 <- haven::read_dta("~/postdoc/stanford/clinical_data/PROMOTE/BC-3 childs all visit database FINAL_withIncidentNMF.dta")


promote_data <- promote %>%
  filter(!is.na(parsdens), !is.na(mstatus))%>%
  group_by(id) %>%
  add_count(name="total_n_visits") %>%
  mutate(n_visit = seq(1, max(total_n_visits)))%>%
  mutate("total_n_para"=sum(parsdens!=0, na.rm=TRUE),
         "total_n_malaria"=sum(mstatus!=0, , na.rm=TRUE),
         "n_para"=if_else(parsdens!=0, cumsum(parsdens!=0), NA),
         "n_malaria"=if_else(mstatus!=0, cumsum(mstatus!=0), NA))%>%
  mutate(mstatus = case_match(mstatus,
                              0~"no malaria",
                              1~"uncomplicated",
                              2~"complicated",
                              3~"quinine for AL failure",
                              4~"Q/AS failure"))%>%
  mutate("hbs"=bc3$hbs[match(id, bc3$id)],
         hbs=case_match(hbs,
                        1~"HbAA",
                        2~"HbAS",
                        3~"HbSS"))

# mic drop ####

mic_drop <-  haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2024_04/MICDROP expanded database through April 30th 2024.dta")

mic_drop_hbs <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/MICDROP SickleTr final.dta")

mic_drop_data <- mic_drop %>%
  filter(!is.na(pardens), !is.na(mstatus))%>%
  group_by(id) %>%
  add_count(name="total_n_visits") %>%
  mutate(n_visit = seq(1, max(total_n_visits)))%>%
  mutate("total_n_para"=sum(pardens!=0),
         "total_n_malaria"=sum(mstatus!=0),
         "n_para"=if_else(pardens!=0, cumsum(pardens!=0), NA),
         "n_malaria"=if_else(mstatus!=0, cumsum(mstatus!=0), NA))%>%
  mutate(mstatus = case_match(mstatus,
                              0~"no malaria",
                              1~"uncomplicated",
                              2~"complicated",
                              3~"quinine for AL failure",
                              4~"Q/AS failure"))%>%
  ungroup()%>%
  mutate(hbs=mic_drop_hbs$HbS[match(as.numeric(id), mic_drop_hbs$id)],
         hbs=case_match(hbs,
                        1~"HbAA",
                        2~"HbAS",
                        3~"HbSS"))

# impact ####

impact <- haven::read_dta("~/postdoc/stanford/clinical_data/IMPACT/IMPACT expanded database through April 30th 2024.dta")

impact_data <- impact %>%
  filter(!is.na(pardens), !is.na(mstatus))%>%
  group_by(id) %>%
  add_count(name="total_n_visits") %>%
  mutate(n_visit = seq(1, max(total_n_visits)))%>%
  mutate("total_n_para"=sum(pardens!=0),
         "total_n_malaria"=sum(mstatus!=0),
         "n_para"=if_else(pardens!=0, cumsum(pardens!=0), NA),
         "n_malaria"=if_else(mstatus!=0, cumsum(mstatus!=0), NA))%>%
  mutate(mstatus = case_match(mstatus,
                              0~"no malaria",
                              1~"uncomplicated",
                              2~"complicated",
                              3~"quinine for AL failure",
                              4~"Q/AS failure"),
         hbs=NA)


# merge it all ####


mic_drop_for_merge <- mic_drop_data %>%
  dplyr::select(id, date, dob, hbs, mstatus, pardens, total_n_visits, total_n_para, total_n_malaria, n_para, n_malaria, temp)%>%
  mutate("study"="micdrop")%>%
  mutate("age"=date-dob)

promote_for_merge <- promote_data %>%
  dplyr::select(id, date, dob, hbs, mstatus, parsdens, total_n_visits, total_n_para, total_n_malaria, n_para, n_malaria, temp)%>%
  mutate("study"="promote", "pardens"=parsdens)%>%
  dplyr::select(-parsdens)%>%
  mutate("age"=date-dob)

impact_for_merge <- impact_data %>%
  dplyr::select(id, date, dob, hbs, mstatus, pardens, total_n_visits, total_n_para, total_n_malaria, n_para, n_malaria, temp)%>%
  mutate("study"="impact")%>%
  mutate("age"=date-dob)

combo_data <- rbind(mic_drop_for_merge , promote_for_merge, impact_for_merge)

three_month_labels <- paste0(seq(0, 21, by=3), " to ", seq(3, 24, by=3), " months")

all_malaria <- combo_data %>%
  ungroup()%>%
  filter(total_n_malaria>=2)%>%
  group_by(id)%>%
  mutate("age_at_first"=na.omit(age[n_malaria==2]))%>%
  ungroup()%>%
  mutate("treatment_failure"=if_else(mstatus%in% c("quinine for AL failure", "Q/AS failure"), 1, 0))%>%
  mutate(agebins=cut(as.numeric(age), breaks = seq(0, max(as.numeric(age)), by=30)))%>%
  mutate(age_at_first_bins=cut(as.numeric(age_at_first), breaks = seq(0, 730, by=90), labels = three_month_labels))%>%
  mutate(disease=if_else(mstatus=="complicated", "complicated", if_else(mstatus=="no malaria", "asymptomatic", "uncomplicated")),
         complicated=if_else(mstatus=="complicated", 1, 0), 
         symptoms=if_else(disease != "asymptomatic", 1, 0))%>%
  mutate(id=factor(id),
         age=as.numeric(age),
         age_months=as.numeric(factor(agebins)),
         log_pardens=log10(pardens+0.1))


  ## add test-positivity rate ####
all_malaria <- all_malaria %>%
  group_by(id)%>%
  mutate("para_positivity_rate"=total_n_para/total_n_visits,
         "malaria_positivity_rate"=total_n_malaria/total_n_visits,
         "para_risk"=total_n_malaria/para_positivity_rate,
         "malaria_risk"=total_n_malaria/malaria_positivity_rate)

# complicated disease modelling ####
n_para_comp <- all_malaria %>%
  group_by(disease, n_para)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated),
         total_infections=complicated+uncomplicated+asymptomatic,
         risk=complicated/total_infections,
         asymp_prob=asymptomatic/total_infections)


comp_model <- glm(risk~n_para+I(n_para^2), family = "binomial", weights = total_infections, data = n_para_comp)

comp_model2 <- glm(complicated~n_para+age, family = "binomial", data = all_malaria)

comp_model3 <- lme4::glmer(complicated~n_para+I(n_para^2)+age_at_first+(1|id), family = "binomial", data = all_malaria)

  ## risk~n_infection ####

#save model as function for plotting
comp_model_fun <- function(x){
  exp(comp_model$coefficients[1])*
    exp(comp_model$coefficients[2])^x*
    exp(comp_model$coefficients[3])^x^2}

#calculate SE for plotting
prd <- data.frame(n_para = seq(from = 1, to = 19, length.out = 100))
err <- predict(comp_model, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit


n_para_comp_plot <- ggplot(n_para_comp[1:19,], aes(x=n_para, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_ribbon(data=prd, aes(x=n_para, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = comp_model_fun, colour="black")+
  geom_text(aes(y=0.10, label= paste0("frac(",complicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  scale_x_continuous(breaks = 1:50, limits=c(1,19))+
  scale_y_continuous(limits = c(0,0.12), labels = scales::label_percent())+
  xlab("Order of Infection")+
  ylab("Risk of Complicated")+
  ggtitle("827 children, aged 8 weeks - 120 weeks
3541 parasitemic episodes; at least two malaria episodes")

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/all_impact_promote_and_micdrop_comp_plot.png", n_para_comp_plot, height = 4.5, width=7.5, dpi=444, bg="white")



n_para_asymp_plot <- ggplot(n_para_comp[1:15,], aes(x=n_para, y=1-asymp_prob))+
  geom_point(color="darkred")+
  theme_minimal()+
  # geom_ribbon(data=prd, aes(x=n_para, ymin = exp(lci), ymax = exp(uci)),
  #             alpha = 0.2, inherit.aes = FALSE)+
  # geom_function(fun = comp_model_fun, colour="black")+
  # facet_wrap(~study)+
  geom_text(aes(y=0.8, label= paste0("frac(",asymptomatic, ",", total_infections,")")),parse = TRUE, size=2.5)+
  geom_smooth(method="lm", color="black")+
  # scale_x_continuous(breaks = 1:50, limits=c(1,15))+
  # scale_y_continuous(limits = c(0,0.12), labels = scales::label_percent())+
  xlab("Order of Infection")+
  ylab("Risk of Symtoms Given Parasitemia")+
  ggtitle("827 children, aged 8 weeks - 2 years
3181 parasitemic episodes; at least one malaria episode")

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/all_para_impact_promote_and_micdrop_asymp_plot.png", n_para_asymp_plot, height = 4.5, width=7.5, dpi=444, bg="white")



age_para <- all_malaria %>%
  filter(pardens!=0)%>%
  mutate(agebins=cut(as.numeric(age), breaks = seq(min(0), max(as.numeric(age)), by=30)))%>%
  ggplot(aes(x=as.numeric(agebins), y=pardens+0.1))+
  geom_point()+
  # see::geom_violindot(aes(fill=factor(agebins)))+
  geom_boxplot(aes(fill=factor(agebins)), outlier.shape = NA)+
  scale_y_log10()+
  scale_x_continuous(breaks=seq(0, 24, by=3), minor_breaks = seq(0, 12))+
  scale_fill_manual(values=colorspace::sequential_hcl(n=26, palette = "Lajolla"))+
  xlab("age in months")+
  ylab("parasites / μL")+
  theme(#axis.text.x = element_text(angle=90),
    legend.position = "none")


n_infection_para <- all_malaria %>%
  filter(n_para<15)%>%
  # filter(pardens>=500)%>%
  # mutate(quarter=cut(as.numeric(factor(agebins)),breaks = seq(0, 24, by=3)))%>%
  ggplot(aes(x=n_para, y=pardens+0.1))+
  # geom_point(aes(color=age_at_first_bins))+
  geom_boxplot(aes(fill=factor(n_para)), outlier.shape = NA)+
  scale_y_log10()+
  # geom_smooth(method="lm", formula = y~x+I(x^2))+
  # facet_wrap(~age_at_first_bins)+
  scale_fill_manual(values=colorspace::sequential_hcl(n=25, palette = "Lajolla"))+
  xlab("order of parasitemic event")+
  ylab("parasites / μL")+
  theme(#axis.text.x = element_text(angle=90),
    legend.position = "none")

n_infection_age_para <- cowplot::plot_grid(age_para, n_infection_para, nrow = 1)
ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/n_para_infection_age_load.png", n_infection_age_para, width=8, height=4)
  ## risk~age####
age_para_comp <- all_malaria %>%
  group_by(disease, agebins)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated),
         total_infections=complicated+uncomplicated+asymptomatic,
         risk=complicated/total_infections,
         asymp_prob=asymptomatic/total_infections,
         age_months=as.numeric(factor(agebins)))

comp_model_age <- glm(risk~age_months+I(age_months^2), family = "binomial", weights = total_infections, data = age_para_comp)

#save model as function for plotting
comp_model_age_fun <- function(x){
  exp(comp_model$coefficients[1])*
    exp(comp_model$coefficients[2])^x*
    exp(comp_model$coefficients[3])^x^2}

#calculate SE for plotting
prd <- data.frame(age_months = seq(from = 1, to = 24, length.out = 100))
err <- predict(comp_model, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit


comp_age <- ggplot(age_para_comp, aes(x=age_months, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_ribbon(data=prd, aes(x=age_months, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = comp_model_age_fun, colour="black")+
  geom_text(aes(y=0.10, label= paste0("frac(",complicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  scale_x_continuous(breaks = seq(0,24, by=3), limits=c(1,24))+
  # scale_y_continuous(limits = c(0,0.12), labels = scales::label_percent())+
  xlab("Age (Months)")+
  ylab("Risk of Complicated")+
  ggtitle("827 children, aged 8 weeks - 120 weeks
3541 parasitemic episodes; at least two malaria episodes")

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/all_impact_promote_and_micdrop_comp_age_plot.png", comp_age, height = 4.5, width=7.5, dpi=444, bg="white")

# complicated statistics ####

comp_model2 <- glm(complicated~n_infection+I(n_infection^2)+age, family = "binomial", data = all_malaria)
comp_model2a <- lme4::glmer(complicated~n_infection+I(n_infection^2)+age+log_pardens+(1|id), family = "binomial", data = all_malaria)

comp_model2b <- MASS::glmmPQL(complicated~n_infection+I(n_infection^2)+age_at_first+log_pardens, random=~1 | id, family = "binomial", data = all_malaria)


combo_comp_hbs <- all_malaria %>%
  filter(!is.na(hbs))%>%
  group_by(disease, n_para, hbs, study)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated),
         uncomplicated=if_else(is.na(uncomplicated), 0, uncomplicated),
         asymptomatic=if_else(is.na(asymptomatic), 0, asymptomatic))%>%
  mutate(total_infections=complicated+uncomplicated+asymptomatic,
         risk=complicated/total_infections)


hbs_comp <- combo_comp_hbs%>%
  filter(hbs!="HbSS")%>%
  ggplot(aes(x=n_para, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  # geom_ribbon(data=prd, aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
  #             alpha = 0.2, inherit.aes = FALSE)+
  # geom_function(fun = comp_model_fun, colour="black")+
  geom_text(aes(y=0.15, label= paste0("frac(",complicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  scale_x_continuous(breaks = 1:50, limits=c(1,15))+
  scale_y_continuous(limits = c(0,0.17), labels = scales::label_percent())+
  geom_smooth(method="glm", formula = y~x+I(x^2))+
  xlab("Order of Infection")+
  ylab("Risk of Complicated")+
  facet_wrap(study~hbs)+
  ggtitle("610 HbAA, 177 HbAS children, aged 8 weeks - 2 years\n3051 parasitemic episodes")

# probability of symptoms, given any parasitemia ####


combo_clin <- all_malaria %>%
  group_by(disease, n_para)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  # group_by(disease)%>%
  group_by(n_para)%>%
  mutate(#complicated=if_else(is.na(complicated), 0, complicated),
    total_infections=sum(c(complicated,uncomplicated,asymptomatic),na.rm = TRUE),
    risk=sum(c(complicated,uncomplicated, na.rm = TRUE))/total_infections
  )


combo_clin_age <- all_malaria %>%
  group_by(disease, agebins)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  # group_by(disease)%>%
  group_by(agebins)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated),
         total_infections=sum(c(complicated,uncomplicated,asymptomatic),na.rm = TRUE),
         risk=sum(c(complicated,uncomplicated, na.rm = TRUE))/total_infections
  )


combo_clin_hbs <- all_malaria %>%
  filter(!is.na(hbs))%>%
  group_by(disease, n_para, hbs)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated))%>%
  mutate(uncomplicated=if_else(is.na(uncomplicated), 0, uncomplicated))%>%
  mutate(asymptomatic=if_else(is.na(asymptomatic), 0, asymptomatic))%>%
  ungroup()%>%
  group_by(hbs, n_para)%>%
  mutate(#complicated=if_else(is.na(complicated), 0, complicated),
    total_infections=sum(complicated,uncomplicated,asymptomatic,na.rm = TRUE),
    risk=sum(complicated,uncomplicated, na.rm = TRUE)/total_infections
  )

# 
# clin_model_lin <- glm(risk~n_infection, family = "binomial", weights = total_infections, data = combo_clin[1:10,])
# 
# clin_model <- glm(risk~agebins, family = "binomial", weights = total_infections, data = combo_clin[1:8,])
# clin_model2 <- glm(risk~as.numeric(agebins), family = "binomial", weights = total_infections, data = combo_clin_age[1:8,])
# 
# clin_model2a <- glm(risk~n_infection+age, family = "binomial", weights = total_infections, data = combo_clin[1:10,])

micdrop_only <- filter(all_malaria, study=="micdrop")
clin_model2a <- lme4::glmer(symptoms~age+hbs+(1|id), family = "binomial", data = micdrop_only)
summary(clin_model2a)


study_contrast <- t(matrix(c(-1,0,0,0,1)))
study_test <- multcomp::glht(clin_model2a, linfct = study_contrast)
summary(study_test)

# clin_model <- lme4::glmer(complicated~n_infection+I(n_infection^2)+(1|id), family = "binomial", data = all_malaria)

# weighted linear regression with quadratic term

#save model as function for plotting
clin_model_fun <- function(x){
  (clin_model$coefficients[1])*
    (clin_model$coefficients[2])^x*
    (clin_model$coefficients[3])^x^2}

#calculate SE for plotting
prd <- data.frame(n_infection = seq(from = 1, to = 10, length.out = 100))
err <- predict(clin_model, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit


clin_n_infection <- ggplot(combo_clin, aes(x=n_para, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  # geom_ribbon(data=prd, aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
  #             alpha = 0.2, inherit.aes = FALSE)+
  # geom_function(fun = clin_model_fun, colour="black")+
  geom_smooth(method="lm", color="black")+
  # geom_text(aes(y=0.19, label= paste0("frac(",complicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  scale_x_continuous(breaks = 1:50, limits=c(1,10))+
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, NA))+
  # facet_wrap(~age)+
  xlab("Order of Parasitemic Episode")+
  ylab("Risk of Malaria When Parasitemic")+
  ggtitle("903 children, aged 8 weeks - 2 years\n3438 parsitemic episodes")


clin_n_age <- ggplot(combo_clin_age, aes(x=age_months, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  # geom_ribbon(data=prd, aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
  #             alpha = 0.2, inherit.aes = FALSE)+
  # geom_function(fun = clin_model_fun, colour="black")+
  geom_smooth(method="lm", color="black")+
  # geom_text(aes(y=0.19, label= paste0("frac(",complicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  scale_x_continuous(breaks = 1:50, limits=c(1,18))+
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, NA))+
  # facet_wrap(~hbs)+
  xlab("Age in Months")+
  ylab("Risk of Malaria When Parasitemic")+
  ggtitle("903 children, aged 8 weeks - 2 years\n3438 parsitemic episodes")


clin_n_infection_age <- cowplot::plot_grid(clin_n_age, clin_n_infection, nrow = 1)
ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/symptoms_per_n_infection_age.png", clin_n_infection_age, width=8, height=4, bg="white")



combo_clin_hbs %>%
  filter(hbs != "HbSS")%>%
  ggplot(aes(x=n_infection, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_smooth(method="glm", method.args=list(family="binomial", weights=no_hbss$total_infections))+
  facet_wrap(~hbs)+
  scale_x_continuous(breaks = 1:50, limits=c(1,10))+
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, NA))+
  geom_text(aes(y=1.1, label= paste0("frac(",complicated+uncomplicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  # ggtitle("903 children, aged 8 weeks - 2 years\n3438 parsitemic episodes")+
  xlab("Order of Infection")+
  ylab("Risk of Malaria When Parasitemic")


  ## probability of symptoms by order of infection given any parasitemia ####

combo_clin_hbs_study <- all_malaria %>%
  filter(!is.na(hbs), pardens>50)%>%
  group_by(disease, n_infection, hbs, study)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated))%>%
  mutate(uncomplicated=if_else(is.na(uncomplicated), 0, uncomplicated))%>%
  mutate(asymptomatic=if_else(is.na(asymptomatic), 0, asymptomatic))%>%
  ungroup()%>%
  group_by(hbs, n_infection, study)%>%
  mutate(#complicated=if_else(is.na(complicated), 0, complicated),
    total_infections=sum(complicated,uncomplicated,asymptomatic,na.rm = TRUE),
    risk=sum(complicated,uncomplicated, na.rm = TRUE)/total_infections
  )

combo_clin_hbs_study %>%
  filter(hbs != "HbSS")%>%
  ggplot(aes(x=n_infection, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  # geom_smooth(method="lm")+
  facet_grid(study~hbs)+
  scale_x_continuous(breaks = 1:50, limits=c(1,10))+
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, NA))+
  geom_text(aes(y=0.9, label= paste0("frac(",complicated+uncomplicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  # ggtitle("903 children, aged 8 weeks - 2 years\n3438 parsitemic episodes")+
  xlab("Order of Infection")+
  ylab("Risk of Malaria When Parasitemic")

all_malaria %>%
  filter(!is.na(hbs), hbs!="HbSS", pardens>50)%>%
  filter(study=="promote")%>%
  ggplot(aes(x=factor(n_infection), y=log_pardens, fill=hbs))+
  geom_boxplot(position = position_dodge())+
  theme_minimal()+
  # geom_smooth(method="lm")+
  facet_grid(~study)+
  # scale_x_continuous(breaks = 1:50, limits=c(1,10))+
  scale_y_log10()+
  scale_fill_manual(values=colorspace::sequential_hcl(n=3, palette = "Lajolla"))+
  # geom_text(aes(y=0.9, label= paste0("frac(",complicated+uncomplicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  # ggtitle("903 children, aged 8 weeks - 2 years\n3438 parsitemic episodes")+
  xlab("Order of Infection")+
  ylab("Parasitemia")
  

## probability of symptoms by age  given parasitemia >500 ####

combo_clin_hbs_age_study <- all_malaria %>%
  filter(!is.na(hbs))%>%
  group_by(disease, agebins, hbs, study)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated))%>%
  mutate(uncomplicated=if_else(is.na(uncomplicated), 0, uncomplicated))%>%
  mutate(asymptomatic=if_else(is.na(asymptomatic), 0, asymptomatic))%>%
  ungroup()%>%
  group_by(hbs, agebins, study)%>%
  mutate(#complicated=if_else(is.na(complicated), 0, complicated),
    total_infections=sum(complicated,uncomplicated,asymptomatic,na.rm = TRUE),
    risk=sum(complicated,uncomplicated, na.rm = TRUE)/total_infections
  )

combo_clin_hbs_age_study %>%
  filter(hbs != "HbSS")%>%
  ggplot(aes(x=as.numeric(agebins), y=risk))+
  geom_point(color="darkred")+
  # geom_line(aes(group=hbs))+
  theme_minimal()+
  geom_smooth(method="lm")+
  facet_grid(study~hbs)+
  scale_x_continuous(breaks = 1:50)+
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, NA))+
  geom_text(aes(y=0.9, label= paste0("frac(",complicated+uncomplicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  # ggtitle("903 children, aged 8 weeks - 2 years\n3438 parsitemic episodes")+
  xlab("Age in Months")+
  ylab("Risk of Malaria When Parasitemic")







all_malaria %>%
  filter(!is.na(hbs), hbs!="HbSS", pardens>50)%>%
  filter(study=="micdrop")%>%
  ggplot(aes(x=factor(agebins), y=log_pardens, fill=hbs))+
  geom_boxplot(position = position_dodge())+
  theme_minimal()+
  # geom_smooth(method="lm")+
  # facet_grid(~study)+
  # scale_x_continuous(breaks = 1:50, limits=c(1,10))+
  scale_y_log10()+
  scale_fill_manual(values=colorspace::sequential_hcl(n=3, palette = "Lajolla"))+
  # geom_text(aes(y=0.9, label= paste0("frac(",complicated+uncomplicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  # ggtitle("903 children, aged 8 weeks - 2 years\n3438 parsitemic episodes")+
  xlab("Order of Infection")+
  ylab("Parasitemia")




para_temp_plot <- all_malaria %>%
  mutate(arm=if_else(age_at_first>350, "IPT", "unknown"))%>%
  filter(total_n_infection>3, hbs !="HbSS", !is.na(hbs))%>%
  ggplot(., aes(x=as.numeric(agebins), y=pardens, group=id))+
  geom_vline(xintercept = 12, linetype="dashed")+
  geom_point(aes(shape=arm, color=factor(round(temp))))+
  geom_line(alpha=0.2)+
  facet_wrap(hbs~id, nrow=3)+
  # ggtitle("")+
  ylab("Parasite Density")+
  xlab("Age in Months")+
  scale_y_log10(limits=c(0.9, 10^6), breaks=c(10^0, 10^2, 10^4, 10^6), labels=scales::label_number())+
  # scale_x_continuous(breaks=seq(c(0,24, by=3)))+
  theme_minimal()+
  # facet_wrap(~mstatus)+
  theme(legend.title = element_blank())+
  scale_color_manual(values=c("darkgreen", "green", "yellow", "orange", "red", "darkred", "#380000"))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/complicated/para_temp_plot2.png", para_temp_plot, width=120, height=4, bg="white", limitsize = FALSE)

