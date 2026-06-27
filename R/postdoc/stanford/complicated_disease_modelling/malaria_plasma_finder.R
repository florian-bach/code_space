micdrop_nulisa <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")


mic_drop_data <- micdrop_specimens %>%
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
  mutate(date=as.character(date))%>%
  select(id, date, n_infection)

infections <- micdrop_nulisa%>%
  filter(mstatus==1)%>%
  distinct(id, date, ageinwks)


n_infections <- left_join(infections, mic_drop_data, by=c("id", "date"))%>%
  mutate(date=as.Date(date))
  



micdrop_nulisa_samples <- micdrop_nulisa%>%
  filter(mstatus!=0)%>%
  distinct(id, date, mstatus)%>%
  mutate(date=as.Date(date))%>%
  mutate(sample_id=paste(id, date, sep="_"))


potential_malaria_samples <- long_specimen_data%>%
  # filter(Specimen_Type%in%c("PlasmaPK", "Plasma"), mstatus!=0)%>%
  filter(Specimen_Type%in%c("Plasma"), mstatus!=0)%>%
  filter(Specimen_ID!="")%>%
  left_join(., slim_inf_data, by=c("id", "date"))%>%
  # select(id, treatmentarm, date, flo_age_in_wks, mstatus.x, n_malaria, n_para, RandomNumber4, BoxNumber4, PositionColumn4, PositionRow4, visittype)%>%
  anti_join(., n_infections, by = c("id", "date"))



first_infections_of_in_pilots <- potential_malaria_samples%>%
  filter(n_malaria ==1, flo_age_in_wks>=66)%>%
  filter(id %in% micdrop_nulisa$id)

non_first_infections_under_two <- potential_malaria_samples%>%
  filter(n_malaria >1, flo_age_in_wks<=106)

first_infections_of_in_pilots %>%
  bind_rows(non_first_infections_under_two)%>%
  select(id, date, BoxNumber3, PositionColumn3, PositionRow3)%>%
  arrange(BoxNumber3)%>%
  write.csv("~/Downloads/mstatus_one.csv", row.names = F)

actual_list <- read.csv("~/Downloads/mstatus_one_25_for_nulisa.csv")

actual_list%>%
  mutate("date"=as.Date(lubridate::parse_date_time(.$date, orders="%m/%d/%y")))%>%
  left_join(., potential_malaria_samples, by=c("id", "date"))%>%
  select(id, date, RandomNumber3)%>%
  write.csv("~/Downloads/mstatus_one_25_for_nulisa_for_kylie.csv", row.names = F)
