#BoxNumber1=PBMC
#BoxNumber2=Paxgene
#BoxNumber3=Plasma
#BoxNumber4=PlasmaPK
#BoxNumber5=CellStabiliser
#BoxNumber7=CellStabiliser


## read in databases
specimen_database <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/Jun25/MICDSpecimenBoxJun25_withclinical.dta")
mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")

sample_ages <- c(8, 24, 52, 56, 60, 64, 68, 84, 104, 120)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1

sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))


long_specimen_data <- specimen_database %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7)%>%
  mutate(stool=ifelse(stool==1, "ordered", "not"))%>%
  select(id, dob, date, flo_age_in_wks, mstatus, qPCRparsdens, ageinwks, SampleDate, starts_with(c("BoxNumber", "PositionColumn", "PositionRow", "RandomNumber")), PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, stool, qPCR, visittype, withdrawaldate) %>%
  mutate("visit_id"=paste(id, date, sep="_"))%>%
  pivot_longer(cols = c(PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR, stool), names_to = "Specimen_Type", values_to = "Specimen_ID")%>%
  mutate(subject_id=id)%>%
  #Specimen_IDs are shared between specimen types, so let's create a unique code
  mutate(Specimen_ID_ID=paste(Specimen_Type, visit_id, sep="_"))%>%
  mutate("Timepoint_in_weeks"=if_else(
    flo_age_in_wks %in% sample_ages, flo_age_in_wks, ifelse(
      flo_age_in_wks %in% sample_ages_minus, flo_age_in_wks+1, if_else(
        flo_age_in_wks %in% sample_ages_plus, flo_age_in_wks-1, 999)
    )
  )
  )%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))


all_pk_samples <- long_specimen_data%>%
  filter(Specimen_Type=="PlasmaPK", Timepoint_in_weeks %in% c(24, 52, 60, 64))%>%
  filter(Specimen_ID!="")%>%
  select(id, treatmentarm, date, Timepoint_in_weeks, Specimen_ID_ID, RandomNumber4, BoxNumber4, PositionColumn4, PositionRow4)%>%
  ungroup()%>%
  group_by(id)%>%
  filter(n()>=4)

fifty_people <- all_pk_samples%>%
  distinct(id, treatmentarm)%>%
  filter(treatmentarm != "DP 2 years")%>%
  group_by(treatmentarm)%>%
  sample_n(size = 25)

(fifty_samples <- all_pk_samples%>%
  filter(id %in% fifty_people$id)%>%
  group_by(id, Timepoint_in_weeks)%>%
  slice_max(order_by = date))

de_identified <- fifty_samples %>%
  select(-treatmentarm)
write.csv(de_identified, "~/postdoc/stanford/plasma_pk/micdrop/plasma_pk_to_pick_aug25.csv", row.names = F)

# plasma pk during malaria ####


## add n_infection data ####
pardens_mic_drop_data <- specimen_database %>%
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


slim_inf_data <- pardens_mic_drop_data%>%
  select(id, date, n_malaria, n_para, mstatus)%>%
  mutate(most_recent_malaria=date)

malaria_pk <- long_specimen_data%>%
  filter(Specimen_Type=="PlasmaPK", mstatus!=0)%>%
  filter(Specimen_ID!="")%>%
  left_join(., slim_inf_data, by=c("id", "date"))%>%
  select(id, treatmentarm, date, flo_age_in_wks, mstatus.x, n_malaria, n_para, RandomNumber4, BoxNumber4, PositionColumn4, PositionRow4, visittype)

kids_with_comp_pk <- malaria_pk%>%
  filter(mstatus.x==2)

kids_with_uncomp_timecourses <- malaria_pk%>%
  filter(n_malaria%in%c(1,3,5) | n_malaria >=7)%>%
  group_by(id)%>%
  slice_min(order_by = n_malaria, n = 4)%>%
  filter(n()==4)


# make new list with clinical only ####

clinical_database <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/MICDroP Data/MICDROP all visit database through July 31st 2025.dta")

mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")

sample_ages <- c(8, 24, 52, 56, 60, 64, 68, 84, 104, 120)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1

sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))


long_malaria <- clinical_database %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7)%>%
  select(id, dob, date, visittype, flo_age_in_wks, mstatus, qPCRparsdens, ageinwks) %>%
  mutate("visit_id"=paste(id, date, sep="_"))%>%
  mutate(subject_id=id)%>%
  filter(visittype==1, !is.na(mstatus), mstatus!=0)%>%
  #Specimen_IDs are shared between specimen types, so let's create a unique code
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))


#