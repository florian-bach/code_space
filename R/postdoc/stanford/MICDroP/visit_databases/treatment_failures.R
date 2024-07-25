library(tidyr)
library(dplyr)
library(ggplot2)

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

comp_pal <- c("asymptomatic"="lightgrey", "uncomplicated"="black", "complicated"="orange", "severe"="darkred")
mstatus_pal <- c("no malaria"="lightgrey", "uncomplicated"="black", "complicated"="orange", "quinine for AL failure"="violet", "Q/AS failure"="purple")

mic_drop <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2024_04/MICDROP expanded database through April 30th 2024.dta")

# merge parasitemia data so that qPCR takes precedent when both slide and qPCR are present
mic_drop <- mic_drop %>%
  mutate(mstatus = case_match(mstatus,
                              0~"no malaria",
                              1~"uncomplicated",
                              2~"complicated",
                              3~"quinine for AL failure",
                              4~"Q/AS failure"))%>%
  mutate(visit_id = paste(id, date, sep=""))%>%
  mutate("parasitaemia_method" = if_else(qPCRdich==1, "qPCR", if_else(BSdich==1, "smear", "dunno")))%>%
  mutate(any_parsdens = if_else(is.na(qPCRparsdens) & !is.na(pardens), pardens, qPCRparsdens))%>%
  mutate(parasitaemia_method = if_else(is.na(qPCRparsdens) & !is.na(pardens), "smear", parasitaemia_method))


# count the number of days between malaria episodes
malaria_calendar <- mic_drop %>%
  select(id, date, AGE, dob, mstatus, any_parsdens, parasitaemia_method)%>%
  filter(mstatus!="no malaria")%>%
  group_by(id) %>%
  add_count(name="number_of_episodes")%>%
  filter(number_of_episodes>1)%>%
  group_by(id) %>%
  arrange(date)%>%
  mutate(
    "previous_symptomatic_date" = lag(date),
    "next_symptomatic_date" = lead(date),
    "time_to_next_symptomatic_date" = next_symptomatic_date-date,
    "time_to_previous_symptomatic_date" = lag(time_to_next_symptomatic_date),
  )



failure_plus_previous_episode <- malaria_calendar %>%
  # filter(date>"2023-06-30")%>%
  filter(time_to_previous_symptomatic_date <= 14 | time_to_next_symptomatic_date <= 14) %>%
  # filter(time_to_next_symptomatic_date > 1) %>%
  mutate(visit_id = paste(id, date, sep=""))
# 
# 
for_abel <- failure_plus_previous_episode %>%
  select(id, date, mstatus, any_parsdens, time_to_next_symptomatic_date)%>%
  filter(date>"2023-07-08")%>%
  arrange(id, date)
  

# write.csv(for_abel, "~/postdoc/stanford/clinical_data/MICDROP/20240424_repeat_infections_for_abel.csv")
# 
# mic_drop %>%
#   filter(visit_id %in% failure_plus_previous_episode$visit_id)%>%
#   select(id, gender, mstatus, smal, date, AGE, heartrate, respirate, temp)



two_week_symptoms_plot <- mic_drop %>%
  # filter(date>"2023-06-30")%>%
  filter(id %in% unique(failure_plus_previous_episode$id), !is.na(any_parsdens))%>%
  group_by(id, date)%>%
  mutate("same_day"=ifelse(n()>1, "same day", "no"))%>%
  ungroup()%>%
  # mutate(id_dob=paste(id, dob))
  # mutate(same_day=ifelse(time_to_next_symptomatic_date==0, "same day", "no"))%>%
  ggplot(., aes(x=date, y=as.numeric(any_parsdens)+0.001))+
  geom_point(aes(color=factor(mstatus, #levels=c("0",
                                         #       "1",
                                          #      "2",
                                           #     "3")
                              ),
                 shape=same_day))+
  geom_line(alpha=0.3, aes(group=id))+
  # ggrepel::geom_text_repel(data=label_df, aes(x=AGE, y=as.numeric(any_parsdens), label=age_in_days))+
  facet_wrap(~ id,nrow = 8)+
  ylab("qPCR parasites / μl\n")+
  xlab("Date")+
  scale_y_log10(breaks=c(1/10, 10, 10^3, 10^5))+
  scale_x_date(breaks="2 week",
               minor_breaks = NULL
               # labels = seq(0,72,by=8)
  )+
  theme_minimal()+
  # scale_shape_manual(values=c(16,15))+
  scale_color_manual(values=mstatus_pal)+
  guides(color=guide_legend(title=""))+
  theme(axis.text.x = element_text(size=5, angle=90, vjust=0.5),
        )


ggsave("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2024_01/figures/14_days_symptoms.png", two_week_symptoms_plot, width = 20, height=12, bg="white")

# since july 2023 ####

failure_plus_previous_episode <- malaria_calendar %>%
  filter(date>"2023-06-30")%>%
  filter(time_to_previous_symptomatic_date <= 14 | time_to_next_symptomatic_date <= 14) %>%
  # filter(time_to_next_symptomatic_date > 1) %>%
  mutate(visit_id = paste(id, date, sep=""))


two_week_symptoms_plot <- mic_drop %>%
  filter(date>"2023-06-30")%>%
  filter(id %in% unique(failure_plus_previous_episode$id), !is.na(any_parsdens))%>%
  group_by(id, date)%>%
  mutate("same_day"=ifelse(n()>1, "same day", "no"))%>%
  ungroup()%>%
  # mutate(id_dob=paste(id, dob))
  # mutate(same_day=ifelse(time_to_next_symptomatic_date==0, "same day", "no"))%>%
  ggplot(., aes(x=date, y=as.numeric(any_parsdens)+0.001))+
  geom_point(aes(color=factor(mstatus, #levels=c("0",
                              #       "1",
                              #      "2",
                              #     "3")
  ),
  shape=same_day))+
  geom_line(alpha=0.3, aes(group=id))+
  # ggrepel::geom_text_repel(data=label_df, aes(x=AGE, y=as.numeric(any_parsdens), label=age_in_days))+
  facet_wrap(~ id,nrow = 4)+
  ylab("qPCR parasites / μl\n")+
  xlab("Date")+
  scale_y_log10(breaks=c(1/10, 10, 10^3, 10^5))+
  scale_x_date(breaks="2 week",
               minor_breaks = NULL
               # labels = seq(0,72,by=8)
  )+
  theme_minimal()+
  # scale_shape_manual(values=c(16,15))+
  scale_color_manual(values=mstatus_pal)+
  guides(color=guide_legend(title=""))+
  theme(axis.text.x = element_text(size=5, angle=90, vjust=0.5),
  )

ggsave("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2024_01/figures/14_days_symptoms_after_Jul23.png", two_week_symptoms_plot, width = 8, height=4, bg="white")



# since january 2024

failure_plus_previous_episode <- malaria_calendar %>%
  filter(date>"2024-01-31")%>%
  filter(time_to_previous_symptomatic_date <= 14 | time_to_next_symptomatic_date <= 14) %>%
  # filter(time_to_next_symptomatic_date > 1) %>%
  mutate(visit_id = paste(id, date, sep=""))


two_week_symptoms_plot <- mic_drop %>%
  filter(date>"2024-01-31")%>%
  filter(id %in% unique(failure_plus_previous_episode$id), !is.na(any_parsdens))%>%
  group_by(id, date)%>%
  mutate("same_day"=ifelse(n()>1, "same day", "no"))%>%
  ungroup()%>%
  # mutate(id_dob=paste(id, dob))
  # mutate(same_day=ifelse(time_to_next_symptomatic_date==0, "same day", "no"))%>%
  ggplot(., aes(x=date, y=as.numeric(any_parsdens)+0.001))+
  geom_point(aes(color=factor(mstatus, #levels=c("0",
                              #       "1",
                              #      "2",
                              #     "3")
  ),
  shape=same_day))+
  geom_line(alpha=0.3, aes(group=id))+
  # ggrepel::geom_text_repel(data=label_df, aes(x=AGE, y=as.numeric(any_parsdens), label=age_in_days))+
  facet_wrap(~ id,nrow = 4)+
  ylab("qPCR parasites / μl\n")+
  xlab("Date")+
  scale_y_log10(breaks=c(1/10, 10, 10^3, 10^5))+
  scale_x_date(breaks="2 week",
               minor_breaks = NULL
               # labels = seq(0,72,by=8)
  )+
  theme_minimal()+
  # scale_shape_manual(values=c(16,15))+
  scale_color_manual(values=mstatus_pal)+
  guides(color=guide_legend(title=""))+
  theme(axis.text.x = element_text(size=5, angle=90, vjust=0.5),
  )

ggsave("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2024_03/figures/14_days_symptoms_after_Jan24.png", two_week_symptoms_plot, width = 8, height=4, bg="white")


# incidence of treatment failure ####


tf_calc <- mic_drop %>% 
  filter(!is.na(mstatus))%>%
  mutate(datebins=cut(date, breaks = seq(min(mic_drop$date), max(mic_drop$date), by=30)))%>%
  filter(!is.na(datebins))%>%
  group_by(datebins, mstatus)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = mstatus, values_from = n)%>%
  mutate(`quinine for AL failure`=if_else(is.na(`quinine for AL failure`), 0, `quinine for AL failure`))%>%
  mutate(all_cases=sum(uncomplicated, complicated, `quinine for AL failure`, na.rm = TRUE))%>%
  mutate(tf_frac=`quinine for AL failure`/all_cases)%>%
  mutate(month_number=as.numeric(datebins))
  

treatment_failure_per_month <- ggplot(tf_calc, aes(datebins, y=tf_frac))+
  geom_bar(stat="identity", fill="darkred")+
  geom_text(aes(label= paste0("frac(",`quinine for AL failure`, ",", all_cases,")")),parse = TRUE, vjust= -0.2, size=3.5)+
  scale_y_continuous(limits = c(0, 0.17), labels = scales::label_percent())+
  theme_minimal()+
  ggtitle("Treatment Failure in MICDRoP (overall 21/818)")+
  xlab("")+
  ylab("AL failure / all cases")+
  theme(axis.text.x = element_text(angle=90))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2024_03/figures/treatment_failure_per_month.png", treatment_failure_per_month, height=5, width=8, bg="white")


tf_model <- glm(tf_frac~month_number, weights = all_cases, data=tf_calc, family = "binomial")

tf_fun <- function(x){exp(tf_model$coefficients[1])*exp(tf_model$coefficients[2])^x}


smooth_treatment_failure_per_month <-ggplot(tf_calc, aes(x=datebins, y=tf_frac))+
  geom_point(color="darkred")+
  geom_ribbon(data= model_visualiser(tf_model, "month_number"), aes(x=month_number, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = tf_fun, colour="black")+
  scale_y_continuous(limits = c(0, 0.17), labels = scales::label_percent())+
  theme_minimal()+
  xlab("")+
  ylab("AL failure / all cases")+
  theme(axis.text.x = element_text(angle=90))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2024_01/figures/smooth_treatment_failure_per_month.png", smooth_treatment_failure_per_month, height=5, width=8, bg="white")



# treatment failures in BC3 ###

promote_data <- haven::read_dta("~/postdoc/stanford/clinical_data/PROMOTE/BC-3 childs all visit database FINAL.dta")

promote_data <- promote_data %>%
  mutate(mstatus = case_match(mstatus,
                            0~"no malaria",
                            1~"uncomplicated",
                            2~"complicated",
                            3~"quinine for AL failure",
                            4~"Q/AS failure"))%>%
  mutate(visit_id = paste(id, date, sep=""))%>%
  mutate("parasitaemia_method" = "smear")


bc3_tf_calc <- promote_data %>% 
  filter(!is.na(mstatus))%>%
  mutate(datebins=cut(date, breaks = seq(min(promote_data$date), max(promote_data$date), by=30)))%>%
  filter(!is.na(datebins))%>%
  group_by(datebins, mstatus)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = mstatus, values_from = n)%>%
  mutate(`quinine for AL failure`=if_else(is.na(`quinine for AL failure`), 0, `quinine for AL failure`))%>%
  mutate(all_cases=sum(uncomplicated, complicated, `quinine for AL failure`, na.rm = TRUE))%>%
  mutate(tf_frac=`quinine for AL failure`/all_cases)%>%
  mutate(month_number=as.numeric(datebins))


bc3_treatment_failure_per_month <- ggplot(bc3_tf_calc, aes(datebins, y=tf_frac))+
  geom_bar(stat="identity", fill="darkred")+
  geom_text(aes(label= paste0("frac(",`quinine for AL failure`, ",", all_cases,")")),parse = TRUE, vjust= -0.2, size=3.5)+
  scale_y_continuous(limits = c(0, 0.17), labels = scales::label_percent())+
  theme_minimal()+
  ggtitle("Treatment Failure in Promote (overall 16/1078)")+
  xlab("")+
  ylab("AL failure / all cases")+
  theme(axis.text.x = element_text(angle=90))

  
library(patchwork)
combo_plot <- treatment_failure_per_month + bc3_treatment_failure_per_month
ggsave("~/postdoc/stanford/clinical_data/treatment_failure/bc3_and_micdrop.png", combo_plot, height=5, width=12, bg="white", dpi=444)


# statistics ####

mic_drop_for_merge <- mic_drop %>%
  select(id, date, dob, mstatus)%>%
  mutate("study"="micdrop")


promote_for_merge <- promote_data %>%
  select(id, date, dob, mstatus)%>%
  mutate("study"="promote")

combo_data <- full_join(mic_drop_for_merge, promote_for_merge)

combo_data$age <- as.numeric(combo_data$date - combo_data$dob)

all_malaria <- combo_data %>%
  filter(mstatus != "no malaria", !is.na(mstatus))%>%
  mutate("treatment_failure"=if_else(mstatus%in% c("quinine for AL failure", "Q/AS failure"), 1, 0))

glm1 <- glm(treatment_failure~age*study, family="binomial", data=all_malaria)
glm1a <- glm(treatment_failure~age+date, family="binomial", data=all_malaria)

glm2 <- lme4::glmer(treatment_failure~study+(1|id), family="binomial", data=all_malaria)



tf_calc_age <- all_malaria %>% 
  filter(!is.na(mstatus))%>%
  mutate(agebins=cut(age, breaks = seq(min(0), max(age), by=30)))%>%
  filter(!is.na(agebins))%>%
  group_by(agebins, mstatus)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = mstatus, values_from = n)%>%
  mutate(`quinine for AL failure`=if_else(is.na(`quinine for AL failure`), 0, `quinine for AL failure`))%>%
  mutate(all_cases=sum(uncomplicated, complicated, `quinine for AL failure`, na.rm = TRUE))%>%
  mutate(tf_frac=`quinine for AL failure`/all_cases)


ggplot(tf_calc_age, aes(x=as.numeric(agebins), y=tf_frac))+
  geom_bar(stat="identity", fill="darkred")+
  geom_text(aes(label= paste0("frac(",`quinine for AL failure`, ",", all_cases,")")),parse = TRUE, vjust= -0.2, size=3.5)+
  scale_y_continuous(limits = c(0, 0.11), labels = scales::label_percent())+
  scale_x_continuous(breaks = seq(0,24,by=3))+
  theme_minimal()+
  ggtitle("Treatment Failure in Combined MC3 and MICDRoP")+
  xlab("Age in Months")+
  ylab("AL failure / all cases")+
  theme()


grouped_tf_calc_age <- all_malaria %>% 
  filter(!is.na(mstatus))%>%
  mutate(agebins=cut(age, breaks = seq(min(0), max(age), by=30)))%>%
  filter(!is.na(agebins))%>%
  group_by(agebins, mstatus, study)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = mstatus, values_from = n)%>%
  mutate(`quinine for AL failure`=if_else(is.na(`quinine for AL failure`), 0, `quinine for AL failure`))%>%
  mutate(all_cases=sum(uncomplicated, complicated, `quinine for AL failure`, na.rm = TRUE))%>%
  mutate(tf_frac=`quinine for AL failure`/all_cases)



tf_by_age_and_study <- ggplot(grouped_tf_calc_age, aes(x=as.numeric(agebins), y=tf_frac))+
  geom_bar(stat="identity", fill="darkred")+
  geom_text(aes(label= paste0("frac(",`quinine for AL failure`, ",", all_cases,")")),parse = TRUE, vjust= -0.2, size=3.5)+
  scale_y_continuous(limits = c(0, 0.11), labels = scales::label_percent())+
  scale_x_continuous(breaks = seq(0,24,by=3))+
  theme_minimal()+
  facet_wrap(~study)+
  ggtitle("Treatment Failure by Age")+
  xlab("Age in Months")+
  ylab("AL failure / all cases")+
  theme()
ggsave("~/postdoc/stanford/clinical_data/treatment_failure/tf_by_age_and_study.png", tf_by_age_and_study, height=5, width=12, bg="white", dpi=444)

