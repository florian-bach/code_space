data <- read.delim("~/Downloads/dantu_chmi.txt")


data %>%
  group_by(SubjectID)%>%
  mutate("endpoint"=if_else(any(qPCR_Parasitaemia>500), "above threshold", "below threshold"))%>%
  ungroup()%>%
  #filter(Phenotypes %in% c("febrile", "susceptible"))%>%
  ggplot(., aes(x=Day, y=qPCR_Parasitaemia, color=factor(SubjectID), group=SubjectID))+
    facet_wrap(Dantu~Phenotypes+endpoint)+
    scale_y_log10()+
    ylab("parasites / µL")+
    geom_hline(yintercept = 10, linetype="dashed", color="grey")+
    geom_hline(yintercept = 500, linetype="dashed", color="black")+
    scale_x_continuous(limits=c(7, NA))+
    geom_point(alpha=0.7)+
    geom_line(alpha=0.7)+
    theme_minimal()+
    theme(legend.position = "none")

data %>%
  group_by(SubjectID)%>%
  mutate("endpoint"=if_else(any(qPCR_Parasitaemia>500), "above threshold", "below threshold"))%>%
  ungroup()%>%
  group_by(Dantu, Phenotypes, endpoint)%>%
  summarise("n"=n_distinct(SubjectID))%>%
  group_by(Dantu)%>%
  mutate("perc"=n/sum(n))

Dantu~Phenotypes+endpoint