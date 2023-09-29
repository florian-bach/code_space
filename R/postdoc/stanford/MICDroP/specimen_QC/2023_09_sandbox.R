no_blank <- raw_data[!is.blank(raw_data$PBMC),]
hist(no_blank$ageinwks)

enrolment <- routine_ish_visits %>%
  filter(visittype %in% c(2))%>%
  group_by(Timepoint_in_weeks)%>%
  summarise("Number_of_Individuals"=n_distinct(Specimen_ID))


Number_of_Visits_per_Timepoint <- routine_ish_visits %>%
  filter(Specimen_ID != "") %>%
  group_by(Timepoint_in_weeks, subject_id) %>%
  summarise("Number_of_Visits_at_Timepoint"=n_distinct(visit_id))

duplications <- Number_of_Visits_per_Timepoint %>%
  filter(Number_of_Visits_at_Timepoint!=1)

duplicate_visits <- routine_ish_visits %>%
  filter(subject_id %in% duplications$subject_id & Timepoint_in_weeks %in% duplications$Timepoint_in_weeks)

duplicate_visits %>%
  # filter(Specimen_ID != "" & visittype != 1)
  filter(Specimen_ID != "")

visits_without_type <-  routine_ish_visits %>%
  filter(is.na(visittype), Specimen_ID != "")




enrolment_and_routine <- routine_ish_visits %>%
  # filter(Specimen_ID != "") %>%
  filter(visittype %in% c(1, 0))%>%
  group_by(Timepoint_in_weeks) %>%
  summarise("Number_of_Individuals"=n_distinct(subject_id))

no_sample_visits <- routine_ish_visits %>%
  filter(Specimen_ID == "" & visittype %notin% c(0,2,NA))

no_samples_2 <- no_sample_visits %>%
  group_by(visit_id) %>%
  summarise("samples_missed_per_visit"=n())%>%
  summarise(table(samples_missed_per_visit))


samples_missed_per_visit <- routine_ish_visits %>%
  filter(Specimen_ID == "" & visittype %in% c(1,0))%>%
  group_by(visit_id) %>%
  summarise("samples_missed_per_visit"=n())%>%
  summarise(table(samples_missed_per_visit))


n_distinct(no_sample_visits$visit_id)


# table(table(no_sample_visits$subject_id))
# 1   2   3   4   5   6   7   8   9 
# 243 138  60  35  15   7   3   1   1 
sum(table(table(no_sample_visits$subject_id)))
# [1] 503

no_visit_type_recorded <- routine_ish_visits %>%
  filter(Specimen_ID != "" & visittype %notin% c(1,0) )

no_visit_type_recorded <- routine_ish_visits %>%
  filter(Specimen_ID  "" & visittype %notin% c(1,0) )

no_samples_2 <- no_sample_visits %>%
  group_by(visit_id) %>%
  summarise("samples_missed_per_visit"=n())%>%
  ungroup()

none_samples <- no_samples_2 %>%
 filter(samples_missed_per_visit==4)



routine_ish_visits %>%
    filter(Timepoint_in_weeks==8)%>%
    group_by(subject_id)%>%
    summarise("n"=n()
              
n_rows <- routine_ish_visits %>%
  filter(Timepoint_in_weeks==8)%>%
  group_by(subject_id)%>%
  summarise("n"=n())

table(n_rows$n)


samples_missing_visits <- routine_ish_visits %>%
  mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
  filter(Specimen_ID == "" & visittype==1)

all_samples_missing_visits <- samples_missing_visits %>%
  group_by(visit_id)%>%
  summarise(id, "n_missing"=n(), Timepoint_in_weeks)%>%
  filter(n_missing==4, !duplicated(visit_id))

blank_visits <- table(all_samples_missing_visits$Timepoint_in_weeks)






#which individuals have more than one visit at the 8 week point
dupli_people <- routine_ish_visits %>%
  mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
  filter(visittype_edit %notin% c(2, NA))%>%
  group_by(subject_id, Timepoint_in_weeks)%>%
  summarise("number_of_rows"=n(), visit_id)

dupli_visits <- subset(dupli_people, dupli_people$number_of_rows==8)
# the 11 extra visits exist because the enrolment visit was like a day before the routine visit 



individuals_not_in_routine <- unique(raw_data$id)[unique(raw_data$id) %notin% unique(routine_ish_visits$id)]

those_visits <- raw_data %>%
  filter(id %in% individuals_not_in_routine)%>%
  dplyr::select(id, dob, ageinwks, SampleDate, PBMC, CellStabilizer, Plasma, Paxgene)

malaria_visits <- raw_data %>%
  filter(incidentmalaria ==1)
