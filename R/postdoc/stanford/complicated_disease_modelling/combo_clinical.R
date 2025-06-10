# preamble ####
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(ggside)
library(visreg)
# library(mediation)
library(patchwork)

`%notin%`=Negate(`%in%`)

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

promote_data <- promote %>%
  filter(mstatus != 0, !is.na(mstatus), incidentmalaria==1 | is.na(incidentmalaria) & mstatus==2 | incidentmalaria==0 & mstatus==2)%>%
  group_by(id) %>%
  add_count(name="total_n_infection") %>%
  mutate(n_infection = seq(1, max(total_n_infection)),
         mstatus = case_match(mstatus,
                              0~"no malaria",
                              1~"uncomplicated",
                              2~"complicated",
                              3~"quinine for AL failure",
                              4~"Q/AS failure"))%>%
  mutate(n_infection=ifelse(incidentmalaria==0 & mstatus=="complicated" | is.na(incidentmalaria) & mstatus=="complicated", n_infection-1, n_infection))

# mic drop ####

mic_drop <-  haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/May25/MICDSpecimenBoxMay25_withclinical.dta")
mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")

mic_drop_data <- mic_drop %>%
  filter(mstatus != 0, !is.na(mstatus), incidentmalaria==1 | is.na(incidentmalaria) & mstatus==2 | incidentmalaria==0 & mstatus==2)%>%
  group_by(id) %>%
  add_count(name="total_n_infection") %>%
  mutate(n_infection = seq(1, max(total_n_infection)),
         mstatus = case_match(mstatus,
                              0~"no malaria",
                              1~"uncomplicated",
                              2~"complicated",
                              3~"quinine for AL failure",
                              4~"Q/AS failure"))%>%
  mutate(n_infection=ifelse(incidentmalaria==0 & mstatus=="complicated" | is.na(incidentmalaria) & mstatus=="complicated", n_infection-1, n_infection))

# impact ####

impact <- haven::read_dta("~/Library/CloudStorage/Box-Box/IMPACT Study (Busia)/Data/Specimens/Mar25/IMPASpecimenBoxMar25_withclinical.dta")


impact_data <- impact %>%
  filter(mstatus != 0, !is.na(mstatus), incidentmalaria==1 | is.na(incidentmalaria) & mstatus==2 | incidentmalaria==0 & mstatus==2)%>%
  group_by(id) %>%
  add_count(name="total_n_infection") %>%
  mutate(n_infection = seq(1, max(total_n_infection)),
         mstatus = case_match(mstatus,
                                     0~"no malaria",
                                     1~"uncomplicated",
                                     2~"complicated",
                                     3~"quinine for AL failure",
                                     4~"Q/AS failure"))%>%
  mutate(n_infection=ifelse(incidentmalaria==0 & mstatus=="complicated" | is.na(incidentmalaria) & mstatus=="complicated", n_infection-1, n_infection))


# merge it all ####

clin_cols <- c("id", "date", "dob", "mstatus", "pardens",
               "total_n_infection", "n_infection", "incidentmalaria",
               "hbgrade", "hb")


mic_drop_for_merge <- mic_drop_data %>%
  select(all_of(clin_cols))%>%
  mutate(study="micdrop")%>%
  mutate("age"=date-dob)

promote_for_merge <- promote_data %>%
  mutate("pardens"=parsdens)%>%
  select(all_of(clin_cols))%>%
  mutate("age"=date-dob, "study"="promote")

impact_for_merge <- impact_data %>%
  select(all_of(clin_cols))%>%
  mutate("study"="impact")%>%
  mutate("age"=date-dob)

combo_data <- rbind(mic_drop_for_merge , promote_for_merge, impact_for_merge)

combo_data <- combo_data%>%
  group_by(id, n_infection)%>%
  # mutate(has_comp=any(mstatus==2))%>%
  # mutate(has_uncomp=any(mstatus==1))%>%
  filter((n() == 1) |
           (mstatus == "complicated")
  )%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(id, mic_drop_key$id)])%>%
  mutate(treatmentarm=case_when(
                                is.na(treatmentarm)&study!="micdrop"~"No DP",
                                treatmentarm==1~"No DP",
                                treatmentarm==2~"DP 1 year",
                                treatmentarm==3~"DP 2 years"))



three_month_labels <- paste0(seq(0, 45, by=3), " to ", seq(3, 48, by=3))
six_month_labels <- paste0(seq(0, 42, by=6), " to ", seq(6, 48, by=6))

dp_ids <- mic_drop_key$id[mic_drop_key$treatmentarm%in%c(2,3)]


all_malaria <- combo_data %>%
  filter(mstatus != "no malaria", !is.na(mstatus))%>%
  group_by(id)%>%
  mutate(age_at_first=age[n_infection==1])%>%
  ungroup()%>%
  mutate("treatment_failure"=if_else(mstatus%in% c("quinine for AL failure", "Q/AS failure"), 1, 0))%>%
  # mutate(agebins=cut(as.numeric(age), breaks = seq(0, max(as.numeric(age)), by=30)))%>%
  mutate(agebins=cut(as.numeric(age), breaks = seq(0, 1460, by=90), labels = three_month_labels))%>%
  mutate(disease=if_else(mstatus=="complicated", "complicated", "uncomplicated"),
         complicated=if_else(mstatus=="complicated", 1, 0))%>%
  mutate(age_at_first_bins=cut(as.numeric(age_at_first),  breaks = seq(0, 1440, by=180), labels = six_month_labels))%>%
  mutate(id=factor(id),
         age=as.numeric(age),
         ageyrs=age/365,
         log_pardens=log10(pardens+0.1))


  
  

combo_comp <- all_malaria %>%
  group_by(disease, n_infection)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated),
         risk=complicated/uncomplicated,
         total_infections=complicated+uncomplicated)
  


comp_model <- glm(risk~n_infection+I(n_infection^2), family = "binomial", weights = total_infections, data = combo_comp)

# # weighted linear regression with quadratic term
# comp_model2 <- lme4::glmer(complicated~n_infection+I(n_infection^2)+(1|id), family = "binomial", data = all_malaria)
# comp_model3 <- lme4::glmer(complicated~n_infection+I(n_infection^2)+ageyrs+(1|id), family = "binomial", data = all_malaria)

#save model as function for plotting
comp_model_fun <- function(x){
  exp(comp_model$coefficients[1])*
    exp(comp_model$coefficients[2])^x*
    exp(comp_model$coefficients[3])^x^2}

#calculate SE for plotting
prd <- data.frame(n_infection = seq(from = 1, to = 17, length.out = 100))
err <- predict(comp_model, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit





combo_comp_plot <- ggplot(combo_comp, aes(x=n_infection, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_ribbon(data=prd, aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = comp_model_fun, colour="black")+
  geom_text(aes(y=0.19, label= paste0("frac(",complicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  scale_x_continuous(breaks = 1:17, limits=c(1,17))+
  scale_y_continuous(labels = scales::label_percent())+
  coord_cartesian(ylim = c(0,0.2))+
  xlab("Order of Malaria Episode")+
  ylab("Risk of Complicated Malaria")+
  # ggtitle(paste(n_distinct(all_malaria$id), " children, aged 8 weeks - 3 years\n", nrow(all_malaria), " malaria episodes", sep=""))+
  theme(panel.grid.minor = element_blank(),
        #axis.title.x = element_blank()
        )

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/no_dp_impact_promote_and_micdrop_comp_plot.png", combo_comp_plot, height = 4.5, width=6, dpi=444, bg="white")

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/impact_promote_and_micdrop_comp_plot.png", combo_comp_plot, height = 4.5, width=6, dpi=444, bg="white")


combo_comp_age <- ggplot(combo_comp2, aes(x=factor(agebins), y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_smooth(aes(x=as.numeric(factor(agebins)), y=risk), method="lm", color="black", inherit.aes = F)+
  # geom_ribbon(data=prd, aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
  #             alpha = 0.2, inherit.aes = FALSE)+
  # geom_function(fun = comp_model_fun, colour="black")+
  geom_text(aes(y=0.19, label= paste0("frac(",complicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  # scale_x_continuous(breaks = 1:10, limits=c(1,10))+
  scale_y_continuous(labels = scales::label_percent())+
  coord_cartesian(ylim=c(0,0.2))+
  # xlab("Age in months")+
  ylab("Risk of Complicated Malaria")+
  # ggtitle(paste(n_distinct(all_malaria$id), " children, aged 8 weeks - 3 years\n", nrow(all_malaria), " malaria episodes", sep=""))+
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/impact_promote_and_micdrop_comp_age.png", combo_comp_age, height = 4.5, width=6, dpi=444, bg="white")

library(patchwork)
combo_combo <- combo_comp_plot+combo_comp_age
ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/impact_promote_and_micdrop_comp_age_combo.png", combo_combo, height = 4.5, width=12, dpi=444, bg="white")


age_para <- all_malaria %>%
  mutate(agebins=cut(as.numeric(age), breaks = seq(min(0), max(as.numeric(age)), by=30)))%>%
  ggplot(aes(x=as.numeric(agebins), y=pardens+0.1))+
  # geom_point()+
  geom_boxplot(aes(fill=factor(agebins)), outlier.shape = NA)+
  scale_y_log10()+
  scale_fill_manual(values=colorspace::sequential_hcl(n=37, palette = "Lajolla"))+
  xlab("age in months")+
  ylab("parasites / μL")+
  theme(#axis.text.x = element_text(angle=90),
        legend.position = "none")
  

n_infection_para <- all_malaria %>%
  # mutate(quarter=cut(as.numeric(factor(agebins)),breaks = seq(0, 24, by=3)))%>%
  ggplot(aes(x=factor(n_infection), y=pardens+0.1))+
  # geom_point(aes(color=age_at_first_bins))+
  geom_boxplot(aes(fill=factor(n_infection)), outlier.shape = NA)+
  scale_y_log10()+
  # facet_wrap(~n_infection)+
  scale_fill_manual(values=colorspace::sequential_hcl(n=17, palette = "Lajolla"))+
  xlab("order of infection")+
  ylab("parasites / μL")+
  theme(#axis.text.x = element_text(angle=90),
    legend.position = "none")

n_infection_age_para <- cowplot::plot_grid(age_para, n_infection_para, nrow = 1)
ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/n_infection_age_load.png", n_infection_age_para, width=8, height=4)

# DP only ####
one_year_dp_ids <- mic_drop_key$id[mic_drop_key$treatmentarm%in%c(2)]

one_dp_malaria <- combo_data %>%
  filter(mstatus != "no malaria", !is.na(mstatus))%>%
  group_by(id)%>%
  mutate(age_at_first=age[n_infection==1])%>%
  ungroup()%>%
  mutate("treatment_failure"=if_else(mstatus%in% c("quinine for AL failure", "Q/AS failure"), 1, 0))%>%
  mutate(agebins=cut(as.numeric(age), breaks = seq(0, max(as.numeric(age)), by=30)))%>%
  mutate(agebins=cut(as.numeric(age), breaks = seq(0, 1460, by=90), labels = three_month_labels))%>%
  mutate(disease=if_else(mstatus=="complicated", "complicated", "uncomplicated"),
         complicated=if_else(mstatus=="complicated", 1, 0))%>%
  mutate(age_at_first_bins=cut(as.numeric(age_at_first), breaks = seq(0, 1460, by=90), labels = three_month_labels))%>%
  mutate(id=factor(id),
         age=as.numeric(age),
         ageyrs=age/365,
         log_pardens=log10(pardens+0.1))%>%
  filter(id %in% dp_ids)



dp_combo_comp <- one_dp_malaria %>%
  group_by(disease, n_infection)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated),
         risk=complicated/uncomplicated,
         total_infections=complicated+uncomplicated)


dp_comp_model <- glm(risk~n_infection+I(n_infection^2), family = "binomial", weights = total_infections, data = dp_combo_comp)

# # weighted linear regression with quadratic term
# comp_model2 <- lme4::glmer(complicated~n_infection+I(n_infection^2)+(1|id), family = "binomial", data = all_malaria)
# comp_model3 <- lme4::glmer(complicated~n_infection+I(n_infection^2)+ageyrs+(1|id), family = "binomial", data = all_malaria)

#save model as function for plotting
dp_comp_model_fun <- function(x){
  exp(dp_comp_model$coefficients[1])*
    exp(dp_comp_model$coefficients[2])^x*
    exp(dp_comp_model$coefficients[3])^x^2}

#calculate SE for plotting
prd <- data.frame(n_infection = seq(from = 1, to = 17, length.out = 100))
err <- predict(dp_comp_model, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit





dp_combo_comp_plot <- ggplot(dp_combo_comp, aes(x=n_infection, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_ribbon(data=prd, aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = dp_comp_model_fun, colour="black")+
  geom_text(aes(y=0.19, label= paste0("frac(",complicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  scale_x_continuous(breaks = 1:17, limits=c(1,17))+
  scale_y_continuous(labels = scales::label_percent())+
  coord_cartesian(ylim = c(0,0.2))+
  # xlab("Order of Malaria Episode")+
  ylab("Risk of Complicated Malaria")+
  # ggtitle(paste(n_distinct(all_malaria$id), " children, aged 8 weeks - 3 years\n", nrow(all_malaria), " malaria episodes", sep=""))+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/dp_impact_promote_and_micdrop_comp_plot.png", dp_combo_comp_plot, height = 4.5, width=6, dpi=444, bg="white")




conf_int_data <- dp_combo_comp %>%
  group_by(n_infection) %>%
  nest() %>%
  mutate(confy_model=purrr::map(data, ~glm(risk~1, data=., weights=total_infections, family="binomial")),
         coef=purrr::map_dbl(confy_model, ~exp(coef(.))),
         upper=purrr::map_dbl(confy_model, ~exp(confint(.))[2]),
         lower=purrr::map_dbl(confy_model, ~exp(confint(.))[1]))


conf_int_plot <- ggplot(conf_int_data, aes(x = n_infection, y = coef))+
  geom_linerange(aes(x = n_infection,  ymin = lower, ymax = upper), linewidth = 1, color="grey", alpha=0.5)+
  geom_point(color="darkred")+
  geom_text(aes(y=0.25, label= paste0("frac(",complicated, ",", total_infections,")")),parse = TRUE, size=2.5, data=dp_combo_comp)+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 0.25), labels = scales::label_percent())+
  xlab("Order of Infection")+
  ylab("Risk of Complicated Malaria")+
  # scale_color_manual(values=age_cat_palette)+
  theme(legend.title = element_blank())

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/dp_impact_promote_and_micdrop_comp_plot.png", conf_int_plot, height = 4.5, width=6, dpi=444, bg="white")


# age ####


age_combo_comp <- all_malaria %>%
  group_by(disease, agebins, )%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated),
         risk=complicated/uncomplicated,
         total_infections=complicated+uncomplicated)


age_comp_model <- glm(risk~as.numeric(agebins), family = "binomial", weights = total_infections, data = age_combo_comp)

# # weighted linear regression with quadratic term
# comp_model2 <- lme4::glmer(complicated~n_infection+I(n_infection^2)+(1|id), family = "binomial", data = all_malaria)
# comp_model3 <- lme4::glmer(complicated~n_infection+I(n_infection^2)+ageyrs+(1|id), family = "binomial", data = all_malaria)

#save model as function for plotting
comp_model_fun <- function(x){
  exp(comp_model$coefficients[1])*
    exp(comp_model$coefficients[2])^x*
    exp(comp_model$coefficients[3])^x^2}

#calculate SE for plotting
prd <- data.frame(n_infection = seq(from = 1, to = 17, length.out = 100))
err <- predict(comp_model, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit





(age_combo_comp_plot <- ggplot(age_combo_comp, aes(x=factor(agebins), y=risk))+
    geom_point(color="darkred")+
    theme_minimal()+
    # geom_ribbon(data=prd, aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
    #             alpha = 0.2, inherit.aes = FALSE)+
    # geom_function(fun = comp_model_fun, colour="black")+
    geom_text(aes(y=0.19, label= paste0("frac(",complicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
    # scale_x_continuous(breaks = 1:17, limits=c(1,17))+
    # scale_y_continuous(labels = scales::label_percent())+
    geom_smooth(method="lm", aes(x=as.numeric(factor(agebins)), y=risk))+
    coord_cartesian(ylim = c(0,0.2))+
    xlab("Age in Months")+
    ylab("Risk of Complicated Malaria")+
    # ggtitle(paste(n_distinct(all_malaria$id), " children, aged 8 weeks - 3 years\n", nrow(all_malaria), " malaria episodes", sep=""))+
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)
    ))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/no_dp_impact_promote_and_micdrop_comp_plot.png", combo_comp_plot, height = 4.5, width=6, dpi=444, bg="white")

# statistics ####

comp_model2 <- glm(complicated~n_infection+I(n_infection^2)+age+age_at_first, family = "binomial", data = all_malaria)

comp_model2a <- lme4::glmer(complicated~n_infection+I(n_infection^2)+log_pardens+(1|id), family = "binomial", data = all_malaria)

comp_model2b <- MASS::glmmPQL(complicated~n_infection+I(n_infection^2)+age_at_first+log_pardens, random=~1 | id, family = "binomial", data = all_malaria)



# biomalpar thoughts ####

all_malaria %>%
  filter(age_at_first<540)%>%
  group_by(disease, n_infection, age_at_first_bins)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated),
         risk=complicated/uncomplicated,
         total_infections=complicated+uncomplicated)%>%
  ggplot(., aes(x=n_infection, y=risk))+
  geom_text(aes(y=0.2, label= paste0("frac(",complicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  geom_point()+
  scale_y_continuous(labels = scales::label_percent())+
  scale_x_continuous(breaks = seq(1,18), limits=c(1,18))+
  coord_cartesian(ylim=c(0,0.33))+
  ggtitle("age at first malaria episode in months")+
  ylab("risk of complicated malaria")+
  xlab("order of incident infection")+
  facet_wrap(~age_at_first_bins, ncol=1)+
  theme_minimal()

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/n_infection_age_at_first.png", width = 4.5, height=9, dpi=444, bg="white")


all_malaria %>%
  mutate(treatmentarm=factor(treatmentarm, levels=c("No DP", "DP 1 year", "DP 2 years")))%>%
  group_by(disease, agebins, treatmentarm)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  group_by(agebins)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated),
         risk=complicated/uncomplicated,
         total_infections=complicated+uncomplicated)%>%
  ggplot(., aes(x=factor(agebins), y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_text(aes(y=0.19, label= paste0("frac(",complicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  # scale_x_continuous(breaks = 1:10, limits=c(1,10))+
  scale_y_continuous(labels = scales::label_percent())+
  coord_cartesian(ylim=c(0,0.33))+
  xlab("age in months")+
  ylab("risk of complicated malaria")+
  facet_wrap(~treatmentarm, ncol=1)+
  # ggtitle(paste(n_distinct(all_malaria$id), " children, aged 8 weeks - 3 years\n", nrow(all_malaria), " malaria episodes", sep=""))+
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/age_age_at_first.png", width = 4.5, height=9, dpi=444, bg="white")


all_malaria %>%
  filter(age_at_first<540)%>%
  mutate(treatmentarm=factor(treatmentarm, levels=c("No DP", "DP 1 year", "DP 2 years")))%>%
  group_by(disease, treatmentarm, n_infection)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  group_by(n_infection, treatmentarm)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated),
         risk=complicated/uncomplicated,
         total_infections=complicated+uncomplicated)%>%
  ggplot(., aes(x=factor(n_infection), y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  # geom_smooth(aes(x=as.numeric(factor(agebins)), y=risk), method="lm", color="black", inherit.aes = F)+
  # geom_ribbon(data=prd, aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
  #             alpha = 0.2, inherit.aes = FALSE)+
  # geom_function(fun = comp_model_fun, colour="black")+
  geom_text(aes(y=0.19, label= paste0("frac(",complicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  # scale_x_continuous(breaks = 1:10, limits=c(1,10))+
  scale_y_continuous(labels = scales::label_percent())+
  coord_cartesian(ylim=c(0,0.33))+
  xlab("order of incident infection")+
  ylab("risk of complicated malaria")+
  facet_wrap(~treatmentarm, ncol=1)+
  # ggtitle(paste(n_distinct(all_malaria$id), " children, aged 8 weeks - 3 years\n", nrow(all_malaria), " malaria episodes", sep=""))+
  theme()
ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/n_infection_treatment.png", width = 4.5, height=9, dpi=444, bg="white")
