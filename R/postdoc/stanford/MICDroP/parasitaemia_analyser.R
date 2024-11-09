library(dplyr)
library(tidyr)
library(ggplot2)

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
  mutate(parasitaemia_method = if_else(is.na(qPCRparsdens) & !is.na(pardens), "smear", parasitaemia_method))



counter_16 <- mic_drop %>%
  select(id, date, AGE, dob, mstatus, any_parsdens, parasitaemia_method)%>%
  filter(mstatus!="no malaria", AGE<=16)%>%
  group_by(id) %>%
  add_count(name="number_of_episodes_in_16")%>%
  filter(number_of_episodes_in_16>=1)

counter_16 <- mic_drop %>%
  select(id, date, AGE, dob, mstatus, any_parsdens, parasitaemia_method)%>%
  filter(mstatus!="no malaria", AGE<=16)%>%
  group_by(id) %>%
  add_count(name="number_of_episodes_in_16")%>%
  filter(number_of_episodes_in_16>=1)



mic_drop %>%
  filter(#id %in% counter_16$id, AGE<=16,
         !is.na(any_parsdens))%>%
  ggplot(., aes(x=AGE, y=as.numeric(any_parsdens)+0.001))+
  geom_point(aes(color=factor(mstatus, #levels=c("0",
                              #       "1",
                              #      "2",
                              #     "3")
  )))+
  geom_line(alpha=0.3, aes(group=id))+
  # ggrepel::geom_text_repel(data=label_df, aes(x=AGE, y=as.numeric(any_parsdens), label=age_in_days))+
  facet_wrap(~ id,nrow = 4)+
  ylab("qPCR parasites / μl\n")+
  xlab("Date")+
  scale_y_log10(breaks=c(1/10, 10, 10^3, 10^5))+
  theme_minimal()+
  # scale_shape_manual(values=c(16,15))+
  scale_color_manual(values=mstatus_pal)+
  guides(color=guide_legend(title=""))+
  theme(axis.text.x = element_text(size=5, angle=90, vjust=0.5))



mic_drop %>%
  filter(!is.na(hb))%>%
  group_by(id)%>%
  mutate(is_ever_anaemic=any(hb<8))%>%
  ungroup()%>%
  filter(is_ever_anaemic==TRUE)%>%
  ggplot(., aes(x=AGE, y=hb))+
  geom_point(aes(color=factor(mstatus, #levels=c("0",
                              #       "1",
                              #      "2",
                              #     "3")
  )))+
  geom_line(alpha=0.3, aes(group=id))+
  geom_hline(yintercept = 8, linetype="dashed")+
  # ggrepel::geom_text_repel(data=label_df, aes(x=AGE, y=as.numeric(any_parsdens), label=age_in_days))+
  facet_wrap(~ id)+
  theme_minimal()+
  # scale_shape_manual(values=c(16,15))+
  scale_color_manual(values=mstatus_pal)+
  guides(color=guide_legend(title=""))+
  theme(axis.text.x = element_text(size=5, angle=90, vjust=0.5))


mic_drop %>%
  filter(!is.na(hb), !is.na(any_parsdens))%>%
  group_by(id)%>%
  mutate(is_ever_anaemic=any(hb<8))%>%
  ungroup()%>%
  filter(is_ever_anaemic==TRUE)%>%
  ggplot(., aes(x=hb, y=any_parsdens))+
  geom_point(aes(color=factor(mstatus, #levels=c("0",
                              #       "1",
                              #      "2",
                              #     "3")
  )))+
  geom_line(alpha=0.3, aes(group=id))+
  geom_vline(xintercept = 8, linetype="dashed")+
  # ggrepel::geom_text_repel(data=label_df, aes(x=AGE, y=as.numeric(any_parsdens), label=age_in_days))+
  # facet_wrap(~ id)+
  
  theme_minimal()+
  scale_y_log10(breaks=c(1/10, 10, 10^3, 10^5))+
  
  # scale_shape_manual(values=c(16,15))+
  scale_color_manual(values=mstatus_pal)+
  guides(color=guide_legend(title=""))+
  theme(axis.text.x = element_text(size=5, angle=90, vjust=0.5))
