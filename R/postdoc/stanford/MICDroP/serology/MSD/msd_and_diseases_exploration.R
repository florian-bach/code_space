diagnoses_of_kids <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/diagnoses_of_kids_in_msd.csv")
raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/MICDroP Data/MICDROP screening database_FINAL.dta")

micdrop_dobs <- raw_data%>%
  distinct(id, dob)

diagnoses_of_kids_in_msd <- diagnoses_of_kids%>%
  left_join(., micdrop_dobs, by="id")%>%
  mutate(date=as.Date(date))%>%
  mutate(flo_age_in_weeks=date-dob)%>%
  filter(flo_age_in_weeks <=365)

personal_summary_of_illnesses <- diagnoses_of_kids_in_msd%>%
  group_by(id, Diagnosis)%>%
  count()%>%
  arrange(desc(n))

wide_personal_summary_of_illnesses <- personal_summary_of_illnesses%>%
  pivot_wider(names_from = Diagnosis, values_from = n)


antibodies_and_illnesses <- long_msd%>%
  filter(timepoint=="52 weeks")%>%
  left_join(., wide_personal_summary_of_illnesses, by="id")

antibodies_and_illnesses%>%
  filter(antigen %in% viruses)%>%
  ggplot(., aes(x=`Anaemia - general`, y=titer))+
  geom_smooth(method="lm")+
  geom_point()+
  scale_y_log10()+
  ggpubr::stat_cor(method = "spearman", color="red")+
  facet_wrap(~antigen)+
  theme_minimal()

antibodies_and_illnesses%>%
  filter(antigen %in% c(vaccines, vaccines_iga[3:5]))%>%
  ggplot(., aes(x=`Common cold / flu`, y=titer))+
  geom_smooth(method="lm")+
  geom_point()+
  scale_y_log10()+
  ylab("antibody concentration at 52 weeks")+
  ggpubr::stat_cor()+
  facet_wrap(~antigen)+
  theme_minimal()



general_summary_of_illnesses <- personal_summary_of_illnesses%>%
  group_by(Diagnosis)%>%
  summarise(n=sum(n))%>%
  arrange(desc(n))



number_of_seroconversions <- wide_long_msd%>%
  filter(antigen %in% viruses, conversion !="nothing")%>%
  group_by(antigen)%>%
  count()%>%
  arrange(desc(n))

number_of_repsiratory_conversions <- number_of_seroconversions%>%
  filter(antigen %in% c("RV.C", "CoV.2", "RSV", "PIV.3", "EV.D68", "hMPV", "PIV.4", "Flu.A.H1", "Flu.B", "Flu.A.H3", "PIV.1", "PIV.2"))%>%
  ungroup()%>%
  summarise(sum(n))

number_of_diarrheal_conversions <- number_of_seroconversions%>%
  filter(antigen %in% c("RV.C", "CoV.2", "RSV", "PIV.3", "EV.D68", "hMPV", "PIV.4", "Flu.A.H1", "Flu.B", "Flu.A.H3", "PIV.1", "PIV.2"))%>%
  ungroup()%>%
  summarise(sum(n))

number_of_repsiratory_illnesses <- personal_summary_of_illnesses%>%
  filter(Diagnosis %in% c())
  
