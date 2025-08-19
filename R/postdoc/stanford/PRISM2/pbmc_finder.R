library(dplyr)

prism_visits <- haven::read_dta("~/postdoc/stanford/clinical_data/PRISM_2/PRISM 2 all visits database through July 31st 2018.dta")

prism_kids <- subset(prism_visits, prism_visits$ageyrs<5)

prism_specimens <- haven::read_dta("~/postdoc/stanford/clinical_data/PRISM_2/PRISM 2 sample database FINAL.dta")

prism_pbmcs <- subset(prism_specimens, prism_specimens$sampletype %in% c("PBMC - 2"))
colnames(prism_pbmcs)[9] <- "Barcode"

prism_locations <- readxl::read_excel("~/postdoc/stanford/clinical_data/PRISM_2/PRISM2PBMC_Aliquot2.xlsx")
prism_locations <- prism_locations%>%
  filter(cohortid %in% prism_kids$cohortid)%>%
  mutate("combo_location", paste(BoxNumber, BoxRow, BoxColumn, sep="_"))

merged_locations <- prism_locations%>%
  left_join(., prism_pbmcs, by = c("Barcode", "date"))
  

samples_for_experiment <- filter(prism_pbmcs, Barcode %in% c("P22BE5", "P2EGWG", "P2DCQX", "P2MEKZ"))

birthdays <- prism_visits %>%
  filter(date %in% samples_for_experiment$date, cohortid %in% samples_for_experiment$cohortid)








prism_adults <- subset(prism_visits, prism_visits$ageyrs>20)

prism_pbmcs <- subset(prism_specimens, prism_specimens$sampletype %in% c("PBMC - 2"))
colnames(prism_pbmcs)[9] <- "Barcode"

prism_locations <- readxl::read_excel("~/postdoc/stanford/clinical_data/PRISM_2/PRISM2PBMC_Aliquot2.xlsx")
prism_locations <- prism_locations%>%
  filter(cohortid %in% prism_adults$cohortid)%>%
  mutate("combo_location", paste(BoxNumber, BoxRow, BoxColumn, sep="_"))

merged_locations <- prism_locations%>%
  left_join(., prism_pbmcs, by = c("Barcode", "date"))%>%
  filter(cellcount2>5, cellcount2<8)%>%
  filter(!duplicated(cohortid.x))

merged_locations%>%
  select(cohortid.x, Barcode, BoxNumber, BoxRow, BoxColumn)

# samples for ICS trial = P2WH6L, P28HKW, P2L43A, P2Q88G


