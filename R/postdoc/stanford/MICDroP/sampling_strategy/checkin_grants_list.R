library(ggplot2)
library(tidyr)
library(ggplot2)


grant_list <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/Samples selected for Florian Oct 30 2024.dta")

comp_pal <- c("asymptomatic"="lightgrey", "uncomplicated"="black", "complicated"="orange", "severe"="darkred")
mstatus_pal <- c("no malaria"="lightgrey", "uncomplicated"="black", "complicated"="orange", "quinine for AL failure"="violet", "Q/AS failure"="purple")

mic_drop <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2024_07/MICDROP expanded database through July 31st 2024.dta")

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
  mutate(parasitaemia_method = if_else(is.na(qPCRparsdens) & !is.na(pardens), "smear", parasitaemia_method))%>%
  group_by(id)



mic_drop %>%
  filter(id %in% grant_list$id, !is.na(lymph), AGE<=52)%>%
  ggplot(aes(x=AGE, y=lymph))+
  geom_point(aes(color=mstatus))+
  geom_line(alpha=0.3, aes(group=id))+
  # ggrepel::geom_text_repel(data=label_df, aes(x=AGE, y=as.numeric(any_parsdens), label=age_in_days))+
  facet_wrap(~ id,nrow = 10)+
  ylab("qPCR parasites / μl\n")+
  xlab("Date")+
  # scale_y_log10(breaks=c(1/10, 10, 10^3, 10^5))+
  theme_minimal()+
  # scale_shape_manual(values=c(16,15))+
  scale_color_manual(values=mstatus_pal)+
  guides(color=guide_legend(title=""))+
  theme(axis.text.x = element_text(size=5, angle=90, vjust=0.5))


grant_kids_para <- mic_drop %>%
  filter(id %in% grant_list$id, !is.na(any_parsdens))%>%
  # select(id, date, AGE, dob, mstatus, any_parsdens, parasitaemia_method)%>%
  filter(AGE<=52)%>%
  group_by(id, mstatus) %>%
  add_count(name="number_of_episodes_in_year1")%>%
  distinct(id, mstatus, number_of_episodes_in_year1)

table(grant_kids_para$mstatus, grant_kids_para$number_of_episodes_in_year1)

comp_cases <- mic_drop %>%
  filter(mstatus == "complicated")%>%
  filter(AGE<=60)%>%
  reframe("id"=id)

kids_with_comp <- comp_cases %>%
  reframe("id"=id)

grant_kids <- mic_drop %>%
  filter(id %in% grant_list$id, !is.na(any_parsdens))%>%
  filter(!is.na(sickletr))
  distinct(id, gender, sickletr)


hello <- mic_drop %>%
  filter(AGE<=52)%>%
  mutate(include=if_else(id %in% grant_list$id, "include", "delay"))%>%
  pivot_longer(cols= grep("grade$", colnames(.), value=TRUE), names_to="AE", values_to="grade")%>%
  filter(grade>0)%>%
  group_by(AE, include)%>%
  summarise("n"=n())%>%
  group_by(AE)%>%
  mutate(sum=sum(n), frac=n/sum)

hi <- mic_drop %>%
  filter(AGE<=52)%>%
  mutate(include=if_else(id %in% grant_list$id, "include", "delay"))%>%
  pivot_longer(cols= grep("grade$", colnames(.), value=TRUE), names_to="AE", values_to="grade")%>%
  filter(grade>0, AE=="respgrade")%>%
  group_by(AE, id, include)%>%
  summarise(n=n())

hi %>%
  group_by(include)%>%
  summarise(n())

mic_drop %>%
  mutate(include=if_else(id %in% grant_list$id, "include", "delay"))%>%
  mutate(any_placenta=rogerson>0)%>%
  # select(id, date, AGE, dob, mstatus, any_parsdens, parasitaemia_method)%>%
  filter(AGE<=52)%>%
  group_by(include, sallergy)%>%
  distinct(id)%>%
  summarise(n=n())%>%
  group_by(sallergy)%>%
  mutate(sum=sum(n), frac=n/sum)
  

mic_drop %>%
  mutate(include=if_else(id %in% grant_list$id, "include", "delay"))%>%
  # select(id, date, AGE, dob, mstatus, any_parsdens, parasitaemia_method)%>%
  filter(AGE<=52)%>%
  group_by(include)%>%
  summarise(mean(birthweight))

all_kids_for_project <- data.frame("id"=unique(c(kids_with_comp$id, grant_kids$id)))
write.csv(all_kids_for_project, file = "~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/all_209_for_project.csv", row.names = F)
