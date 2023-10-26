library(dplyr)
library(tidyr)
library(ggplot2)
library(ggside)
# library(visreg)
# library(mediation)
library(patchwork)

#promote
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
  )%>%
  mutate(first_inf_before_four = if_else(age_at_first<4/12, "before 4 months", "after 4 months"))%>%
  mutate(first_inf_before_six = if_else(age_at_first<0.5, "before 6 months", "after 6 months"))

promote <- complicated_data%>%
  group_by(id)%>%
  summarise("n_comp"=sum(complicatedmalaria=="complicated"))%>%
  group_by(n_comp)%>%
  count()

# Groups:   n_comp [4]
# n_comp     n
# <int> <int>
#   1      0   367
# 2      1    53
# 3      2     3
# 4      3     2

promote_modeled <- dpois(0:5, 58/425)

promote_kids_with_pf=425
promote_kids_with_any=58
promote_kids_with_one=53
promote_kids_with_two=3
promote_kids_with_three=2

#goncalves ####
#see supplementary table s6 

gonc_kids_with_pf=715
gonc_kids_with_any=102
n_comp=122

gonc_kids_with_one=87
gonc_kids_with_two=12
gonc_kids_with_three=1
gonc_kids_with_four=2

gonc_modeled <- dpois(0:5, kids_with_any/kids_with_pf)

poisson <- data.frame("n_infection"=0:5,
                      "gonc_modeled"=gonc_modeled,
                      "promote_modeled"=promote_modeled,
                      "gonc_actual"=c(
                        (gonc_kids_with_pf-kids_with_any)/gonc_kids_with_pf,
                        gonc_kids_with_one/gonc_kids_with_pf,
                        gonc_kids_with_two/gonc_kids_with_pf,
                      gonc_kids_with_three/gonc_kids_with_pf,
                      gonc_kids_with_four/gonc_kids_with_pf,
                      0),
                      "promote_actual"= c((promote_kids_with_pf-promote_kids_with_any)/promote_kids_with_pf,
                                          promote_kids_with_one/promote_kids_with_pf,
                                          promote_kids_with_two/promote_kids_with_pf,
                                          promote_kids_with_three/promote_kids_with_pf,
                                          0,
                                          0)
                      ) 



long_poisson <- poisson %>%
  pivot_longer(cols = c(gonc_modeled, promote_modeled, gonc_actual, promote_actual), names_to = "source", values_to = "case_perc")%>%
  mutate("study"=if_else(grepl("gonc*", source), "Gonçalves", "Promote"),
         "type"=if_else(grepl("*modeled", source), "modeled", "actual"))




log_poiss <- ggplot(long_poisson, aes(x=n_infection, y=case_perc, color=type))+
  geom_point()+
  geom_line()+
  # geom_ribbon()+
  facet_wrap(~study)+
  xlab("Number of Children with n Severe Episodes")+
  ylab("Fraction of all Children")+
  scale_color_manual(values = c("darkred", "darkblue"))+
  scale_y_log10(limits=c(0.00001,1), labels = scales::label_percent(accuracy = 0.01), breaks=c(0.0001, 0.001, 0.01, 0.1, 1))+
  annotation_logticks(sides = "l")+
  theme_minimal()+
  theme(legend.title = element_blank())

ggsave("~/postdoc/stanford/clinical_data/PROMOTE/figures/log_complicated_poisson.png", log_poiss, width=8, height=6, bg = "white")


non_log_poiss <- ggplot(long_poisson, aes(x=n_infection, y=case_perc, color=type))+
  geom_point()+
  geom_line()+
  facet_wrap(~study)+
  xlab("Number of Children with n Severe Episodes")+
  ylab("Fraction of all Children")+
  scale_color_manual(values = c("darkred", "darkblue"))+
  theme_minimal()+
  theme(legend.title = element_blank())

ggsave("~/postdoc/stanford/clinical_data/PROMOTE/figures/complicated_poisson.png", non_log_poiss, width=8, height=6, bg = "white")

# confint  ####

gonc_simul <- list()

for(i in 1:1000){
  single_simul=table(rpois(715, kids_with_severe/kids_with_pf))
  gonc_simul[[i]]=single_simul
}

gonc_simul_df <- do.call(bind_rows, gonc_simul)


gonc_simul_df <- data.frame(apply(gonc_simul_df, 2, FUN = function(x)ifelse(is.na(x), 0, x)))



gonc_simul_means <- apply(gonc_simul_df, 2, FUN = function(x)mean(x))
gonc_simul_sd <- apply(gonc_simul_df, 2, FUN = function(x)sd(x))
gonc_simul_se <- gonc_simul_sd/sqrt(nrow(gonc_simul_df))

gonc_stats <- data.frame("n_infection"=0:4,
                         "mean"=gonc_simul_means,
                         "upper"=gonc_simul_means+(gonc_simul_sd))



