# mysterious third sampls ####
people_with_three_samples <- long_specimen_data %>%
  filter(Specimen_Type %in% c("PBMC", "Paxgene", "CellStabiliser"))%>%
  group_by(Specimen_Type, id) %>%
  mutate("Number_of_Samples" = n())%>%
  filter(Number_of_Samples==3)%>%
  arrange(Specimen_Type, id)

third_sample_visits <- people_with_three_samples %>%
  filter(flo_age_in_wks %notin% sample_ranges)

all_incorrect_visits <- long_specimen_data %>%
  filter(flo_age_in_wks %notin% sample_ranges & Specimen_ID != "" & Specimen_Type %in% c("PBMC", "Paxgene", "CellStabiliser"))

# kids with broken birthdays ####
bbd_id <- raw_data[is.na(raw_data$dob),"id"]

table(c(raw_data[raw_data$id%in%bbd_id$id, "dob"], raw_data[raw_data$id%in%bbd_id$id, "id"]), useNA = "ifany")

# 14 kids have unknown birthday, 
# 10443 10447 10449 10458 10548 10667 10682 10755 10875 10883 10945 10987 11382
# 10883 has a birthday it's just missing in some entries

# long_specimen_data <- long_specimen_data %>%
#   mutate(dob=ifelse(id==10883, as.Date("2022-05-29"), dob))
         