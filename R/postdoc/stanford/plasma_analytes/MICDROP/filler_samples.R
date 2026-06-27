library(tidyr)
library(dplyr)
library(xlsx)
library(ggplot2)

`%notin%` <- Negate(`%in%`)

is.blank <- function(x){sapply(x, function(y) {ifelse(y=="", TRUE, FALSE)})}


extra_samples <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/extra_24/extra_24_sample_manifest_for_tran.csv")

raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/Nov24/MICDSpecimenBoxNov24_withclinical.dta")

grant_and_comp <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/all_209_for_project.csv")

nulisa_pilot <- readxl::read_excel("/Users/fbach/postdoc/stanford/plasma_analytes/MICDROP/mic_drop_samples_for_nulisa.xlsx")

sample_ages <- c(8, 24, 52, 68, 84, 104, 120)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1

sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))

#turn the data into long format, include a couple of convenience variables, drop irrelevant columns
long_specimen_data <- raw_data %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7)%>%
  select(id, dob, date, flo_age_in_wks, mstatus, qPCRparsdens, pardens, ageinwks, SampleDate, starts_with(c("BoxNumber", "PositionColumn", "PositionRow", "RandomNumber")), PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR, visittype, withdrawaldate) %>%
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

day0_plasmas <- long_specimen_data %>%
  filter(!is.na(mstatus))%>%
  group_by(id)%>%
  # make variable "days_since_malaria" that calculates the days a malaria episode and resets to 0 when mstatus!=0
  mutate(malaria_episode_num_date = if_else(mstatus != 0, as.numeric(date), -Inf),
         # malaria_episode_dich = if_else(mstatus != 0, 1, 0))%>%
         most_recent_malaria = as.Date(cummax(malaria_episode_num_date)),
         malaria_episode_id = paste(id, most_recent_malaria,sep="_"),
         days_since_malaria = date-most_recent_malaria)%>%
  # filter(Specimen_Type=="PBMC", Specimen_ID!="", days_since_malaria > 3 & days_since_malaria<14)%>%
  filter(Specimen_Type%in%c("Plasma"), Specimen_ID!="", days_since_malaria==0)%>%
  select(id, date, flo_age_in_wks, malaria_episode_id, mstatus, Specimen_ID, days_since_malaria, qPCRparsdens, pardens)

day0_plasmas_not_run <- day0_plasmas%>%
  filter(id %notin% grant_and_comp$id)%>%
  filter(id %notin% extra_samples$id)%>%
  #nulisa pilot
  filter(id %notin% c(10766, 10794, 10842))%>%
  #already padded
  filter(id %notin% c(11807, 10950))

nulisa_fill_ids <- day0_plasmas_not_run%>%
  filter(!is.na(qPCRparsdens))%>%
  filter(id %in% c(11807, 10950, 10766, 10651))

nulisa_fill_locations <- long_specimen_data %>%
  filter(id %in% day0_plasmas_not_run$id)%>%
  dplyr::filter(id %in% c(11651, 10947, 11700), Timepoint_in_weeks<70)%>%
  dplyr::filter(Specimen_ID!="", Specimen_Type %in% c("Plasma"))%>%
  arrange(id, Specimen_Type)%>%
  select(RandomNumber3, id, Timepoint_in_weeks, BoxNumber3, PositionColumn3, PositionRow3)%>%
  arrange(BoxNumber3)

write.csv(nulisa_fill_locations, "~/postdoc/stanford/plasma_analytes/MICDROP/dpsp_filler_samples.csv", row.names = F)



