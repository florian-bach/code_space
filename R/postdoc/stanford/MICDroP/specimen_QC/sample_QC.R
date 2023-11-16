library(dplyr)
library(tidyr)
library(xlsx)

`%notin%` <- Negate(`%in%`)
is.blank <- function(x){sapply(x, function(y) {ifelse(y=="", TRUE, FALSE)})}

raw_data <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_10/MICDSpecimenBoxOct23_withclinical.dta")

# for finding putative sampling visits we'll filter the database to only include visits around the proper sampling timepoint dates, with a week plus/minus
sample_ages <- c(8, 24, 52, 68, 84, 104, 120)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1

sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))

#turn the data into long format, include a couple of convenience variables, drop irrelevant columns
long_specimen_data <- raw_data %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7)%>%
  select(id, dob, date, flo_age_in_wks, mstatus, qPCRparsdens, ageinwks, SampleDate, PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR, visittype) %>%
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


# subset data to include all putative sampling visits, regardless of coded visit type
routine_ish_visits <- long_specimen_data %>%
  filter(flo_age_in_wks %in% sample_ranges)

# contrast with subset of data coded with the appropriate routine visit type
routine_visits <- long_specimen_data %>%
  filter(visittype==1)

# subset routine-ish visits to only include 8 week visits that
putative_eight <- routine_ish_visits %>%
  # mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
  filter(visittype %notin% c(2, NA))%>%
  filter(Timepoint_in_weeks==8) %>%
  filter(!duplicated(visit_id))%>%
  select(subject_id, date, visit_id)


definitive_eight <- routine_visits %>%
  filter(Timepoint_in_weeks==8) %>%
  filter(!duplicated(visit_id))%>%
  select(subject_id, date, visit_id, visittype)




uncoded_visits <- putative_eight %>%
  filter(visit_id %notin% definitive_eight$visit_id)
write.csv(uncoded_visits, "~/Box Sync/MIC_DroP IPTc Study/Florian_Specimen_QC/2023_10/uncoded_visits.csv")

unknown_ids <- raw_data %>%
  filter(id %notin% putative_eight$subject_id)%>%
  distinct(id)
write.csv(unknown_ids, "~/Box Sync/MIC_DroP IPTc Study/Florian_Specimen_QC/2023_10/unknown_ids.csv")

unknown_visits <- long_specimen_data %>%
  filter(id %notin% putative_eight$subject_id)%>%
  filter(!duplicated(visit_id))
write.csv(unknown_visits, "~/Box Sync/MIC_DroP IPTc Study/Florian_Specimen_QC/2023_10/unknown_id_visits.csv")

missing_birthdays <- long_specimen_data %>%
  filter(is.na(dob))%>%
  filter(!duplicated(visit_id))%>%
  select(id, date)
write.csv(missing_birthdays, "~/Box Sync/MIC_DroP IPTc Study/Florian_Specimen_QC/2023_10/missing_birthdays.csv")