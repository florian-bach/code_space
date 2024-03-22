
long_specimen_data <- raw_data %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7)%>%
  # select(id, dob, date, flo_age_in_wks, mstatus, qPCRparsdens, ageinwks, SampleDate, PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR, visittype, withdrawaldate) %>%
  select(id, dob, date, flo_age_in_wks, mstatus, qPCRparsdens, ageinwks, SampleDate, PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR, visittype, withdrawaldate, starts_with(c("BoxNumber", "PositionColumn", "PositionRow"))) %>%
  mutate("visit_id"=paste(id, date, sep="_"))%>%
  pivot_longer(cols = c(PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR), names_to = "Specimen_Type", values_to = "Specimen_ID")%>%
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



plasma_boxes_in_stanford_march24 <- paste("MICD", c(1001:1026, 1028, 1029, 1031:1038), sep="-")

#BoxNumber3 = plasmas
withdrawn_plasmas <- long_specimen_data %>%
  filter(Specimen_Type == "Plasma",
    BoxNumber3 %in% plasma_boxes_in_stanford_march24,
    !is.na(withdrawaldate))%>%
  group_by(id)%>%
  add_count()
  
samples_to_pick <- withdrawn_plasmas%>%
  filter(n==3)%>%
  select(id, ageinwks, BoxNumber3, PositionColumn3, PositionRow3)

# table(withdrawn_plasmas$ageinwks)                  
# 8  9 24 52 
# 31  2 21  1 

