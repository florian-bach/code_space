clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")
mic_drop <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/May25/MICDSpecimenBoxMay25_withclinical.dta")

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

counter <- mic_drop %>%
  select(id, date, AGE, dob, mstatus, any_parsdens, parasitaemia_method)%>%
  filter(mstatus!="no malaria", AGE<=107)%>%
  filter(id %in% kids_without_full_one_year$id)%>%
  group_by(id) %>%
  add_count(name="number_of_episodes_in_year1")%>%
  filter(number_of_episodes_in_year1>=3)




mic_drop %>%
  filter(id %in% counter$id, !is.na(mstatus), !is.na(any_parsdens))%>%
  ggplot(., aes(x=date, y=any_parsdens+0.01, group=id))+
  geom_line()+
  geom_point(aes(color=mstatus))+
  scale_y_log10()+
  facet_wrap(~id)+
  scale_color_manual(values = mstatus_pal)+
  theme_minimal()

samples_for_scott <- long_specimen_data %>%
  filter(id %in% counter$id)%>%
  filter(id %notin% unique(clean_data$id))%>%
  filter(Timepoint_in_weeks %in% c(104),
         Specimen_ID!="", Specimen_Type %in% c("PBMC"))%>%
  # filter(!is.na(withdrawaldate))%>%
  arrange(id, Specimen_Type)%>%
  select(id, date, Timepoint_in_weeks, BoxNumber1, PositionColumn1, PositionRow1, BoxNumber2, PositionColumn2, PositionRow2)%>%
  arrange(BoxNumber1)

write.csv(samples_for_scott, "~/Downloads/samples_for_scott2.csv", row.names = F)


# samples for Alex B cell 10x pilot:
# 5 placebo & 5 DP1 kids, 52 and 104 weeks -> 20 samples
# ideally matched 12 to 24 month number of infections
msd_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/batch_one.csv")
mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")
specimen_database <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/Jun25/MICDSpecimenBoxJun25_withclinical.dta")

kids_with_comp <- specimen_database%>%
  filter(mstatus==2)%>%
  distinct(id)


dobs <- specimen_database%>%
  filter(!is.na(dob))%>%
  distinct(id, dob)

infs <- specimen_database%>%
  mutate(dob2=dobs$dob[match(id, dobs$id)])%>%
  mutate("flo_age_in_wks"=as.numeric(date-dob2)%/%7)%>%
  # filter(ageinwks<54)%>%
  group_by(id) %>%
  mutate("total_n_para_12"=sum(pardens>1&flo_age_in_wks<53, na.rm=TRUE),
         "total_n_malaria_12"=sum(mstatus!=0&flo_age_in_wks<53, na.rm=TRUE),
         "total_n_para_24"=sum(pardens>1&flo_age_in_wks<105, na.rm=TRUE),
         "total_n_malaria_24"=sum(mstatus!=0&flo_age_in_wks<105, na.rm=TRUE),
         "total_n_para_12_24"=sum(pardens>1&flo_age_in_wks<105&flo_age_in_wks>52, na.rm=TRUE),
         "total_n_malaria_12_24"=sum(mstatus!=0&flo_age_in_wks<105&flo_age_in_wks>52, na.rm=TRUE),
         "total_n_para_6"=sum(pardens>1&flo_age_in_wks<27, na.rm=TRUE),
         "total_n_malaria_6"=sum(mstatus!=0&flo_age_in_wks<27, na.rm=TRUE),
         "any_malar_6"=if_else(total_n_malaria_6==0, FALSE, TRUE),
         "any_malar_12"=if_else(total_n_malaria_12==0, FALSE, TRUE))%>%
  distinct(id, total_n_para_6, total_n_malaria_6, total_n_para_12, total_n_malaria_12, total_n_para_24, total_n_malaria_24, total_n_para_12_24, total_n_malaria_12_24)

msd_kids <- msd_data%>%
  mutate(id=SubjectID)%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))%>%
  dplyr::distinct(id, id)%>%
  left_join(., infs, by="id")


msd_kids_select <- msd_kids%>%
  filter(total_n_para_12_24%in%c(3,4,5))%>%
  select(id, treatmentarm, total_n_para_12, total_n_para_12_24, total_n_malaria_12, total_n_malaria_12_24)%>%
  filter(!id %in% kids_with_comp$id)

table(msd_kids_select$id %in% kids_with_comp$id)
# FALSE 
# 10 

long_luminex <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/long_luminex.csv")


luminex_kids <- long_luminex%>%
  dplyr::distinct(id)%>%
  left_join(., infs, by="id")%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))
  

luminex_kids_select <- luminex_kids %>%
  filter(id %in% msd_kids$id)%>%
  select(id, treatmentarm, total_n_para_12, total_n_para_12_24, total_n_malaria_12, total_n_malaria_12_24)%>%
  filter(total_n_para_12_24 %in% c(3,4,5))%>%
  filter(!id %in% kids_with_comp$id)


long_luminex%>%
  filter(antigen %in% c("schizont", "tetanus", "CIDRg3 IT4var08", "CIDRa1.7 2083-1"))%>%
  mutate(treatmentarm_x_pilot=if_else(id %in% c(luminex_kids_select$id, msd_kids_select$id),
                                      paste(treatmentarm, "in pilot"), 
                                      paste(treatmentarm, "not in pilot")))%>%
  ggplot(., aes(x = factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks", "104 weeks")), y=log_mfi, fill=treatmentarm_x_pilot))+
  geom_boxplot(outliers = F)+
  scale_fill_manual(values=c("turquoise", "#00555A", "red", "darkred"))+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()+
  theme(axis.title=element_blank())


potentials <- bind_rows(luminex_kids_select, msd_kids_select)%>%
  mutate("serology"=case_when(id %in% long_luminex$id & id %in% msd_kids$id~"both",
                              id %in% long_luminex$id~"luminex",
                              id %in% msd_kids$id~"msd"
  ))%>%
  arrange(treatmentarm, serology)%>%
  filter(!(total_n_para_12==0&treatmentarm=="Placebo"))


# ids_for_alex_pilot <- potentials[c(1:4, 6,
#                                    7:9, 11, 14),]

long_specimen_database <- specimen_database%>%
  mutate(dob2=dobs$dob[match(id, dobs$id)])%>%
  mutate("flo_age_in_wks"=as.numeric(date-dob2)%/%7)%>%
  select(id, dob, date, flo_age_in_wks, mstatus, qPCRparsdens, ageinwks, SampleDate, starts_with(c("BoxNumber", "PositionColumn", "PositionRow", "RandomNumber")), PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, , qPCR, visittype, withdrawaldate) %>%
  mutate("visit_id"=paste(id, date, sep="_"))%>%
  tidyr::pivot_longer(cols = c(PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR), names_to = "Specimen_Type", values_to = "Specimen_ID")%>%
  mutate(subject_id=id)%>%
  #Specimen_IDs are shared between specimen types, so let's create a unique code
  mutate(Specimen_ID_ID=paste(Specimen_Type, visit_id, sep="_"))%>%
  mutate("Timepoint_in_weeks"=if_else(
    flo_age_in_wks %in% sample_ages, flo_age_in_wks, ifelse(
      flo_age_in_wks %in% sample_ages_minus, flo_age_in_wks+1, if_else(
        flo_age_in_wks %in% sample_ages_plus, flo_age_in_wks-1, 999)
    )
  )
  )

sample_locations_for_pilot <- long_specimen_database%>%
  filter(id %in% potentials$id, Specimen_Type=="PBMC", flo_age_in_wks > 48 & flo_age_in_wks < 56 | flo_age_in_wks > 102 & flo_age_in_wks < 108)%>%
  filter(Specimen_ID!="")%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))%>%
  select(id, flo_age_in_wks, treatmentarm, date, RandomNumber1, BoxNumber1, PositionColumn1, PositionRow1)%>%
  arrange(treatmentarm, id, date)%>%
  group_by(id)%>%
  filter(n()==2)

write.csv(sample_locations_for_pilot, file = "~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/pbmc_locations_for_alex_b_cell_pilot.csv", row.names = F)



metadata_for_alex <- sample_locations_for_pilot%>%
  select(-treatmentarm)%>%
  filter(id %notin% c(11637, 10669))%>%
  left_join(potentials, by=c("id"))%>%
  distinct(id, flo_age_in_wks, RandomNumber1, total_n_para_12, total_n_para_12_24,
           total_n_malaria_12, total_n_malaria_12_24)

write.csv(metadata_for_alex, file = "~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/metadata_for_alex_b_cell_pilot.csv", row.names = F)
