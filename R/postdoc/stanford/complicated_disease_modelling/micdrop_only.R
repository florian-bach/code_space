# preamble ####
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_minimal())
library(patchwork)
mic_drop <-  haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/Mar25/MICDSpecimenBoxMar25_withclinical.dta")

mic_drop_hbs <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/MICDROP SickleTr final.dta")

mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")


# parasitemia as episode ####
pardens_mic_drop_data <- mic_drop %>%
  filter(!is.na(pardens), !is.na(mstatus))%>%
  group_by(id) %>%
  add_count(name="total_n_visits") %>%
  mutate(n_visit = seq(1, max(total_n_visits)))%>%
  mutate("total_n_para"=sum(pardens!=0, na.rm = T),
         "total_n_malaria"=sum(mstatus!=0, na.rm = T),
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
         treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         hbs=case_match(hbs,
                        1~"HbAA",
                        2~"HbAS",
                        3~"HbSS"),
         treatmentarm=case_match(treatmentarm,
                                1~"Placebo",
                                2~"DP 1 year",
                                3~"DP 2 years"))


slim_pardens_mic_drop_data <- pardens_mic_drop_data %>%
  dplyr::select(id, treatmentarm, date, dob, hbs, mstatus, pardens, total_n_visits, total_n_para, total_n_malaria, n_para, n_malaria, temp, anyDP, gender)%>%
  mutate("study"="micdrop")%>%
  mutate("age"=date-dob)

three_month_labels <- paste0(seq(0, 57, by=3), " to ", seq(3, 60, by=3), " months")

all_para <- slim_pardens_mic_drop_data %>%
  ungroup()%>%
  filter(total_n_malaria>=1)%>%
  group_by(id)%>%
  mutate("age_at_first"=ifelse(any(n_malaria==1), age[n_malaria==1], NULL))%>%
  ungroup()%>%
  mutate("treatment_failure"=if_else(mstatus%in% c("quinine for AL failure", "Q/AS failure"), 1, 0))%>%
  mutate(agebins=cut(as.numeric(age), breaks = seq(0, max(as.numeric(age)), by=30)))%>%
  mutate(age_at_first=cut(as.numeric(age_at_first), breaks = seq(0, 1800, by=90), labels = three_month_labels))%>%
  mutate(disease=if_else(mstatus=="complicated", "complicated", if_else(mstatus=="no malaria", "asymptomatic", "uncomplicated")),
         complicated=if_else(mstatus=="complicated", 1, 0), 
         symptoms=if_else(disease != "asymptomatic", 1, 0))%>%
  mutate(id=factor(id),
         age=as.numeric(age),
         age_months=as.numeric(factor(agebins)),
         log_pardens=log10(pardens+0.1))

#all para did not have the filtering on total_n_malaria
# write.csv(all_para, "~/postdoc/stanford/clinical_data/all_para_mar25.csv", row.names = F)

n_para_comp <- all_para %>%
  group_by(disease, n_para, anyDP)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated),
         total_infections=complicated+uncomplicated+asymptomatic,
         risk=complicated/total_infections,
         asymp_prob=asymptomatic/total_infections)


# individual level models
comp_model3 <- lme4::glmer(complicated~n_para+I(n_para^2)+(1|id), family = "binomial", data = all_para)

#do this to add two comp infections in 2 year DP group
para_placebo_data <- all_para %>%
  filter(anyDP=="no")%>%
  filter(!is.na(n_para))

para_dp_data <- all_para%>%
  filter(anyDP=="yes")%>%
  filter(!is.na(n_para))

# individual level models
para_dp_model <- lme4::glmer(complicated~n_para+I(n_para^2)+(1|id), family = "binomial", data = para_dp_data)
para_placebo_model <- lme4::glmer(complicated~n_para+I(n_para^2)+(1|id), family = "binomial", data = para_placebo_data)

# summary model
para_placebo_summary <- n_para_comp%>%
  filter(anyDP=="no")

para_dp_summary <- n_para_comp%>%
  filter(anyDP=="yes")

placebo_comp_model <- glm(risk~n_para+I(n_para^2), family = "binomial", weights=total_infections, data = para_placebo_summary)
dp_comp_model <- glm(risk~n_para+I(n_para^2), family = "binomial", weights=total_infections, data = para_dp_summary)

para_dp_model <-  glm(risk~n_para+I(n_para^2), family = "binomial", weights=total_infections, data = para_dp_summary)
para_placebo_model <- glm(risk~n_para+I(n_para^2), family = "binomial", weights=total_infections, data = para_placebo_summary)

## risk ~ n_para plots ####

#save model as function for plotting
para_placebo_model_fun <- function(x){
  exp(para_placebo_model$coefficients[1])*
    exp(para_placebo_model$coefficients[2])^x*
    exp(para_placebo_model$coefficients[3])^x^2}

para_dp_model_fun <- function(x){
  exp(para_dp_model$coefficients[1])*
    exp(para_dp_model$coefficients[2])^x*
    exp(para_dp_model$coefficients[3])^x^2}

#calculate SE for plotting
para_placebo_model_prd <- data.frame(n_para = seq(from = 1, to = 20, length.out = 100))
para_placebo_model_err <- predict(para_placebo_model, newdata = para_placebo_model_prd, se.fit = TRUE)

para_placebo_model_prd$lci <- para_placebo_model_err$fit - 1.96 * para_placebo_model_err$se.fit
para_placebo_model_prd$fit <- para_placebo_model_err$fit
para_placebo_model_prd$uci <- para_placebo_model_err$fit + 1.96 * para_placebo_model_err$se.fit


para_dp_model_prd <- data.frame(n_para = seq(from = 1, to = 20, length.out = 100))
para_dp_model_err <- predict(para_dp_model, newdata = para_dp_model_prd, se.fit = TRUE, re.form=NA)

para_dp_model_prd$lci <- para_dp_model_err$fit - 1.96 * para_dp_model_err$se.fit
para_dp_model_prd$fit <- para_dp_model_err$fit
para_dp_model_prd$uci <- para_dp_model_err$fit + 1.96 * para_dp_model_err$se.fit


placebo_n_para_comp_plot <- ggplot(para_placebo_summary, aes(x=n_para, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_ribbon(data=para_placebo_model_prd, aes(x=n_para, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = placebo_comp_model_fun, colour="black")+
  geom_text(aes(y=0.10, label= paste0("frac(",complicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  scale_x_continuous(breaks = 1:50, limits=c(1,10))+
  scale_y_continuous(limits = c(0,0.22), labels = scales::label_percent())+
  xlab("Order of Infection")+
  ylab("Risk of Complicated")
#   ggtitle("827 children, aged 8 weeks - 120 weeks
# 3541 parasitemic episodes")


dp_n_para_comp_plot <- ggplot(para_dp_summary, aes(x=n_para, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_ribbon(data=para_dp_model_prd, aes(x=n_para, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = dp_comp_model_fun, colour="black")+
  geom_text(y=0.15, aes(label= paste0("frac(",complicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  # scale_x_continuous(breaks = 1:50, limits=c(1,10))+
  scale_y_continuous(limits = c(0,0.22), labels = scales::label_percent())+
  xlab("Order of Infection")+
  ylab("Risk of Complicated")

placebo_n_para_comp_plot+dp_n_para_comp_plot
ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/micdrop_only_qpcr_n_para_comp_plot.png", n_para_comp_plot, height = 4.5, width=7.5, dpi=444, bg="white")

# malaria episode as n_infection ####
n_malaria_mic_drop_data <- mic_drop %>%
  filter(mstatus != 0, !is.na(mstatus))%>%
  group_by(id) %>%
  add_count(name="total_n_malaria") %>%
  mutate("n_malaria"=if_else(mstatus!=0, cumsum(mstatus!=0), NA))%>%
  mutate(mstatus = case_match(mstatus,
                              0~"no malaria",
                              1~"uncomplicated",
                              2~"complicated",
                              3~"quinine for AL failure",
                              4~"Q/AS failure"))%>%
  ungroup()%>%
  mutate(hbs=mic_drop_hbs$HbS[match(as.numeric(id), mic_drop_hbs$id)],
         treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         hbs=case_match(hbs,
                        1~"HbAA",
                        2~"HbAS",
                        3~"HbSS"),
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))


slim_n_malaria_mic_drop_data <- n_malaria_mic_drop_data %>%
  dplyr::select(id, treatmentarm, date, dob, hbs, mstatus, pardens, total_n_malaria, n_malaria, temp, anyDP)%>%
  mutate("study"="micdrop")%>%
  mutate("age"=date-dob)

three_month_labels <- paste0(seq(0, 57, by=3), " to ", seq(3, 60, by=3), " months")

all_malaria <- slim_n_malaria_mic_drop_data %>%
  ungroup()%>%
  filter(total_n_malaria>=1)%>%
  group_by(id)%>%
  mutate("age_at_first"=ifelse(any(n_malaria==1), age[n_malaria==1], NULL))%>%
  ungroup()%>%
  mutate("treatment_failure"=if_else(mstatus%in% c("quinine for AL failure", "Q/AS failure"), 1, 0))%>%
  mutate(agebins=cut(as.numeric(age), breaks = seq(0, max(as.numeric(age)), by=30)))%>%
  mutate(age_at_first=cut(as.numeric(age_at_first), breaks = seq(0, 1800, by=90), labels = three_month_labels))%>%
  mutate(disease=if_else(mstatus=="complicated", "complicated", if_else(mstatus=="no malaria", "asymptomatic", "uncomplicated")),
         complicated=if_else(mstatus=="complicated", 1, 0), 
         symptoms=if_else(disease != "asymptomatic", 1, 0))%>%
  mutate(id=factor(id),
         age=as.numeric(age),
         age_months=as.numeric(factor(agebins)),
         log_pardens=log10(pardens+0.1))


comp_model3 <- lme4::glmer(complicated~n_malaria+I(n_malaria^2)+(1|id), family = "binomial", data = all_malaria)
comp_model <- glm(risk~n_malaria+I(n_malaria^2), family = "binomial", weights=total_n_malaria, data = n_malaria_comp)

n_malaria_comp <- all_malaria %>%
  group_by(disease, n_malaria, anyDP)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated),
         total_n_malaria=complicated+uncomplicated,
         risk=complicated/total_n_malaria)

placebo_data <- all_malaria%>%
  filter(anyDP=="no")

dp_data <- all_malaria%>%
  filter(anyDP=="yes")

# individual level models
dp_model <- lme4::glmer(complicated~n_malaria+I(n_malaria^2)+(1|id), family = "binomial", data = dp_data)
placebo_model <- lme4::glmer(complicated~n_malaria+I(n_malaria^2)+(1|id), family = "binomial", data = placebo_data)

# summary model
placebo_summary <- n_malaria_comp%>%
  filter(anyDP=="no")

dp_summary <- n_malaria_comp%>%
  filter(anyDP=="yes")

placebo_comp_model <- glm(risk~n_malaria+I(n_malaria^2), family = "binomial", weights=total_n_malaria, data = placebo_summary)
dp_comp_model <- glm(risk~n_malaria+I(n_malaria^2), family = "binomial", weights=total_n_malaria, data = dp_summary)

## risk ~ n_malaria plots ####

#save model as function for plotting
placebo_comp_model_fun <- function(x){
  exp(placebo_comp_model$coefficients[1])*
    exp(placebo_comp_model$coefficients[2])^x*
    exp(placebo_comp_model$coefficients[3])^x^2}

dp_comp_model_fun <- function(x){
  exp(dp_comp_model$coefficients[1])*
    exp(dp_comp_model$coefficients[2])^x*
    exp(dp_comp_model$coefficients[3])^x^2}

#calculate SE for plotting
placebo_comp_model_prd <- data.frame(n_malaria = seq(from = 1, to = 20, length.out = 100))
placebo_comp_model_err <- predict(placebo_comp_model, newdata = placebo_comp_model_prd, se.fit = TRUE)

placebo_comp_model_prd$lci <- placebo_comp_model_err$fit - 1.96 * placebo_comp_model_err$se.fit
placebo_comp_model_prd$fit <- placebo_comp_model_err$fit
placebo_comp_model_prd$uci <- placebo_comp_model_err$fit + 1.96 * placebo_comp_model_err$se.fit


dp_comp_model_prd <- data.frame(n_malaria = seq(from = 1, to = 20, length.out = 100))
dp_comp_model_err <- predict(dp_comp_model, newdata = dp_comp_model_prd, se.fit = TRUE)

dp_comp_model_prd$lci <- dp_comp_model_err$fit - 1.96 * dp_comp_model_err$se.fit
dp_comp_model_prd$fit <- dp_comp_model_err$fit
dp_comp_model_prd$uci <- dp_comp_model_err$fit + 1.96 * dp_comp_model_err$se.fit


placebo_n_malaria_comp_plot <- ggplot(placebo_summary, aes(x=n_malaria, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_ribbon(data=placebo_comp_model_prd, aes(x=n_malaria, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = placebo_comp_model_fun, colour="black")+
  geom_text(aes(y=0.10, label= paste0("frac(",complicated, ",", total_n_malaria,")")),parse = TRUE, size=2.5)+
  scale_x_continuous(breaks = 1:50, limits=c(1,10))+
  scale_y_continuous(limits = c(0,0.22), labels = scales::label_percent())+
  ggtitle("Placebo")+
  xlab("Order of Infection")+
  ylab("Risk of Complicated")
#   ggtitle("827 children, aged 8 weeks - 120 weeks
# 3541 parasitemic episodes")


dp_n_malaria_comp_plot <- ggplot(dp_summary, aes(x=n_malaria, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_ribbon(data=dp_comp_model_prd, aes(x=n_malaria, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = dp_comp_model_fun, colour="black")+
  geom_text(aes(y=0.10, label= paste0("frac(",complicated, ",", total_n_malaria,")")),parse = TRUE, size=2.5)+
  scale_x_continuous(breaks = 1:50, limits=c(1,10))+
  scale_y_continuous(limits = c(0,0.22), labels = scales::label_percent())+
  ggtitle("any DP")+
  xlab("Order of Infection")+
  ylab("Risk of Complicated")

n_para_comp_plot <- placebo_n_malaria_comp_plot+dp_n_malaria_comp_plot
ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/micdrop_only_n_malaria_comp_plot.png", n_para_comp_plot, height = 4.5, width=7.5, dpi=444, bg="white")
