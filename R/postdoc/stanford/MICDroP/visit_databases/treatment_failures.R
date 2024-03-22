library(dplyr)
library(tidyr)
library(ggplot2)

comp_pal <- c("asymptomatic"="lightgrey", "uncomplicated"="black", "complicated"="orange", "severe"="darkred")
mstatus_pal <- c("no malaria"="lightgrey", "uncomplicated"="black", "complicated"="orange", "quinine for AL failure"="violet", "Q/AS failure"="purple")

mic_drop <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2023_07/MICDROP all visit database through July 31st 2023.dta")

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
  filter(time_to_previous_symptomatic_date < 21 | time_to_next_symptomatic_date < 21) %>%
  # filter(time_to_next_symptomatic_date > 1) %>%
  mutate(visit_id = paste(id, date, sep=""))


for_abel <- failure_plus_previous_episode %>%
  select(id, date, mstatus, any_parsdens, time_to_next_symptomatic_date)
write.csv(for_abel, "~/postdoc/stanford/clinical_data/MICDROP/repeat_infections_for_abel.csv")

mic_drop %>%
  filter(visit_id %in% failure_plus_previous_episode$visit_id)%>%
  select(id, gender, mstatus, smal, date, AGE, heartrate, respirate, temp)



two_week_symptoms_plot <- mic_drop %>%
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
  scale_x_date(breaks="2 week"
               # labels = seq(0,72,by=8)
  )+
  theme_minimal()+
  # scale_shape_manual(values=c(16,15))+
  scale_color_manual(values=mstatus_pal)+
  guides(color=guide_legend(title=""))+
  theme(axis.text.x = element_text(size=5, angle=90, vjust=0.5))


ggsave("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2023_07/figures/21_days_symptoms.png", two_week_symptoms_plot, width = 24, height=8, bg="white")
