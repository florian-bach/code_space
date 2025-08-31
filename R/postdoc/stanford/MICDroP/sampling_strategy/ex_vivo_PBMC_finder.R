library(tidyr)
library(dplyr)
library(xlsx)
library(ggplot2)

`%notin%` <- Negate(`%in%`)
is.blank <- function(x){sapply(x, function(y) {ifelse(y=="", TRUE, FALSE)})}

# all except lavstsen plasmas were done with this
# raw_data <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2024_06/MICDSpecimenBoxJun24_withclinical.dta")
mic_drop_hbs <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/MICDROP SickleTr final.dta")

mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")

raw_data <-  haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/Jun25/MICDSpecimenBoxJun25_withclinical.dta")


# for finding putative sampling visits we'll filter the database to only include visits around the proper sampling timepoint dates, with a week plus/minus
sample_ages <- c(8, 24, 52, 68, 84, 104, 120)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1

sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))

#turn the data into long format, include a couple of convenience variables, drop irrelevant columns
long_specimen_data <- raw_data %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7)%>%
  mutate(stool=ifelse(stool==1, "ordered", "not"))%>%
  select(id, dob, date, flo_age_in_wks, mstatus, qPCRparsdens, ageinwks, SampleDate, starts_with(c("BoxNumber", "PositionColumn", "PositionRow", "RandomNumber")), PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, stool, qPCR, visittype, withdrawaldate) %>%
  mutate("visit_id"=paste(id, date, sep="_"))%>%
  pivot_longer(cols = c(PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR, stool), names_to = "Specimen_Type", values_to = "Specimen_ID")%>%
  mutate(subject_id=id)%>%
  #Specimen_IDs are shared between specimen types, so let's create a unique code
  mutate(Specimen_ID_ID=paste(Specimen_Type, visit_id, sep="_"))%>%
  mutate("Timepoint_in_weeks"=if_else(
    flo_age_in_wks %in% sample_ages, flo_age_in_wks, ifelse(
      flo_age_in_wks %in% sample_ages_minus, flo_age_in_wks+1, if_else(
        flo_age_in_wks %in% sample_ages_plus, flo_age_in_wks-1, 999)
    )
  )
  )%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(id, mic_drop_key$id)],
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))
  
# finding "ex vivo" samples ####
ex_vivo <- long_specimen_data %>%
  filter(!is.na(mstatus))%>%
  group_by(id)%>%
  # make variable "days_since_malaria" that calculates the days a malaria episode and resets to 0 when mstatus!=0
  mutate(malaria_episode_num_date = if_else(mstatus != 0, as.numeric(date), -Inf),
         # malaria_episode_dich = if_else(mstatus != 0, 1, 0))%>%
         most_recent_malaria = as.Date(cummax(malaria_episode_num_date)),
         malaria_episode_id = paste(id, most_recent_malaria,sep="_"),
         days_since_malaria = date-most_recent_malaria)







# parasitemia as episode ####

pardens_mic_drop_data <- raw_data %>%
  filter(!is.na(pardens), !is.na(mstatus))%>%
  group_by(id) %>%
  add_count(name="total_n_visits") %>%
  mutate(n_visit = seq(1, max(total_n_visits)))%>%
  mutate("total_n_para"=sum(pardens!=0, na.rm = T),
         "total_n_malaria"=sum(mstatus!=0, na.rm = T),
         "n_para"=if_else(pardens!=0, cumsum(pardens!=0), NA),
         "n_malaria"=if_else(mstatus!=0, cumsum(mstatus!=0), NA))%>%
  mutate(mstatus = case_match(mstatus,
                              0~"no malaria",
                              1~"uncomplicated",
                              2~"complicated",
                              3~"quinine for AL failure",
                              4~"Q/AS failure"))%>%
  ungroup()%>%
  mutate(hbs=mic_drop_hbs$HbS[match(as.numeric(id), mic_drop_hbs$id)],
         treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         hbs=case_match(hbs,
                        1~"HbAA",
                        2~"HbAS",
                        3~"HbSS"),
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))


slim_inf_data <- pardens_mic_drop_data%>%
  select(id, date, n_malaria, n_para, mstatus)%>%
  mutate(most_recent_malaria=date)
  

ids_with_two <- ex_vivo %>%
  filter(days_since_malaria>=3, days_since_malaria<=8, Specimen_Type=="PBMC", BoxNumber1!="")%>%
  select(id, flo_age_in_wks, days_since_malaria, most_recent_malaria, BoxNumber1, PositionColumn1, PositionRow1)%>%
  left_join(., slim_inf_data, by=c("id", "most_recent_malaria"))%>%
  group_by(id)%>%
  add_count(name = "n_ex_vivo")%>%
  filter(n_ex_vivo>1)




