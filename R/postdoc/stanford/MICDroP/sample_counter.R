library(dplyr)
library(tidyr)

raw_data <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/MICDropSpecimensJan2023_withClinical.dta")

sample_ages <- seq(8, 88, by=16)
# 
sample_ranges <- sort(c(sample_ages, sample_ages-1, sample_ages+1))

long_specimen_data <- raw_data |>
  select(id, date, ageinwks, SampleDate, SpecimenID1, SpecimenID2, SpecimenID3, SpecimenID4, SpecimenID5, SpecimenID6, visittype) |>
  pivot_longer(cols = c(SpecimenID1, SpecimenID2, SpecimenID3, SpecimenID4, SpecimenID5, SpecimenID6), names_to = "Specimen_Type", values_to = "Specimen_ID")|>
  filter(Specimen_ID!="")|>
  mutate(subject_id=id)|>
  mutate(Specimen_Type=recode(Specimen_Type,
                            "SpecimenID1"="PBMC",
                            "SpecimenID2"="Paxgene",
                            "SpecimenID3"="Plasma",
                            "SpecimenID4"="PlasmaPK",
                            "SpecimenID5"="CellStabiliser",
                            "SpecimenID6"="qPCR"))

participant_summary <- long_specimen_data |>
  filter(ageinwks %in% sample_ranges)|>
  group_by(ageinwks)|>
  summarise("Number_of_Individuals"=n_distinct(id))


sample_summary <- long_specimen_data |>
  filter(ageinwks %in% sample_ranges)|>
  group_by(Specimen_Type, ageinwks) |>
  summarise("Number_of_Samples"=n())


sample_summary$Collection_Rate <- sample_summary$Number_of_Samples/participant_summary$Number_of_Individuals[match(sample_summary$ageinwks, participant_summary$ageinwks)]


sample_summary <- sample_summary |>
  mutate("Collection_Failure_Rate"=1-Collection_Rate,
         "Missed_Samples"=participant_summary$Number_of_Individuals[match(ageinwks, participant_summary$ageinwks)]-Number_of_Samples)


n_sample_summary <- long_specimen_data |>
  group_by(Specimen_Type, id) |>
  summarise("Number_of_Samples" = n()) |>
  group_by(Specimen_Type,  Number_of_Samples) |>
  summarise("Number_of_Individuals_with_n_Samples"=n())