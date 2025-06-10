clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")
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
  filter(mstatus!="no malaria", AGE<=104)%>%
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
