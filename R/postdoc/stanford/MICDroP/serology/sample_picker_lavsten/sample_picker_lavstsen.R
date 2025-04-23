msd_map <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/msd_map_edit.csv")
non_msd_map <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/non_msd_edit.csv")

slim_msd_map <- msd_map%>%
  mutate(reorder.column=rep_len(1:9, nrow(.)),
         reorder.row=rep_len(rep(1:9,each=9), nrow(.)),
         reorder.box=rep_len(rep(1:9, each=81), nrow(.)))%>%
  select(id, date, Timepoint_in_weeks, reorder.column, reorder.row, reorder.box)
  
slim_non_msd_map <- non_msd_map%>%
  mutate("Timepoint_in_weeks"=flo_age_in_wks)%>%
  mutate(reorder.box=reorder.box+3)%>%
  select(id, date, Timepoint_in_weeks, reorder.column, reorder.row, reorder.box)


combo_map <- rbind(slim_non_msd_map, slim_msd_map)

write.csv(combo_map, "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/existing_aliquot_map.csv")

grant_list_100 <- readxl::read_excel("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/Samples selected for Florian April 17 2025.xls")

# samples to plate for sweden ####
swedes <- combo_map%>%
  filter(Timepoint_in_weeks %in% c(51, 52))%>%
  arrange(reorder.box)

# samples to plate for denmark
danes_plate <- combo_map%>%
  filter(id %in% grant_list_100$id)%>%
  arrange(reorder.box)

# samples to pick for denmark
danes_pick <- long_specimen_data %>%
  filter(id %in% grant_list_100$id, Specimen_Type=="Plasma")%>%
  filter(Specimen_ID!="", flo_age_in_wks<106, flo_age_in_wks>102, flo_age_in_wks %in% sample_ranges)%>%
  select(id, date, Timepoint_in_weeks, BoxNumber3, PositionColumn3, PositionRow3)%>%
  arrange(BoxNumber3)%>%
  mutate("new column"=rep_len(1:9, nrow(.)),
         "new row"=rep_len(rep(1:9,each=9), nrow(.)),
         "new box"=rep_len(rep(1:9, each=81), nrow(.)))

write.csv(danes_pick, "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/104_week_timepoints_to_pick.csv", row.names = F)
