
stool_samples1 <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/microbiome/micdrop_stool_samples_dec_2025.csv")

stool_samples2 <- readxl::read_excel("~/postdoc/stanford/clinical_data/MICDROP/microbiome/LIST OF STOOL SAMPLES ALIQUOT 1.xlsx")

bdays <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2024_06/MICDSpecimenBoxJun24_withclinical.dta", col_select = c(id, dob))
bdays <- bdays[!duplicated(bdays),]


full_stools <- stool_samples%>%
  mutate(Date = as.Date(lubridate::parse_date_time(.$Date, orders="%m/%d/%y")))%>%
  bind_rows(., stool_samples2)%>%
  select(SubjectID, Date, RandomNumber, BoxNumber, PositionRow, PositionColumn)%>%
  mutate(dob= bdays$dob[match(.$SubjectID, bdays$id)])

mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")


sample_ages <- c(8, 24, 52, 68, 84, 104, 120, 136, 156, 172, 188, 208)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1

sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))


full_stool_details <- full_stools %>%
  mutate("age_in_weeks"=as.numeric(as.Date(Date)-dob)%/%7)%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(SubjectID, mic_drop_key$id)],
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))


kids_with_8_24_52 <- full_stool_details %>%
  filter(age_in_weeks<60)%>%
  group_by(SubjectID)%>%
  add_count(name = "n_samples")%>%
  filter(n_samples>=3)

kids_with_24_and_52 <- full_stool_details %>%
  filter(age_in_weeks>20 & age_in_weeks<60)%>%
  group_by(SubjectID)%>%
  add_count(name = "n_samples")%>%
  filter(n_samples>=2)

kids_with_6_12_24 <- full_stool_details %>%
  filter(age_in_weeks>22 & age_in_weeks<26 | age_in_weeks>50 & age_in_weeks<54 |age_in_weeks>102 & age_in_weeks<106)%>%
  group_by(SubjectID)%>%
  add_count(name = "n_samples")%>%
  filter(n_samples>=3)








# old and bad code ####


mic_drop <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2024_09/MICDROP expanded database through September 30th 2024.dta")

dobs <- mic_drop%>%
  distinct(id, dob)

stool_samples1 <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/microbiome/micdrop_stool_samples_dec_2025.csv")

grant_and_comp <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/all_209_for_project.csv")

sample_ages <- c(8, 24, 52, 68, 84, 104, 120)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1
sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))

stool_samples <- stool_samples %>%
  mutate("dob"=dobs$dob[match(.$SubjectID, dobs$id)])%>%
  mutate("Date"=as.Date(lubridate::parse_date_time(.$Date, orders="%m/%d/%y")))%>%
  mutate("age"=as.numeric(as.Date(Date)-dob)%/%7)%>%
  # this selects only aliquot 1s
  filter(grepl("^MICD-ST-1", BoxNumber))

stool_sample_match <- stool_samples%>%
  filter(SubjectID %in% grant_and_comp$id)


# all_stool_samples <- long_specimen_data %>%
#   filter(id %in% grant_and_comp$id, Specimen_Type=="stool")%>%
#   select(id, flo_age_in_wks, BoxNumber6)%>%
#   filter(BoxNumber6!="")

kids_with_6_and_12 <- stool_samples %>%
  filter(age>22 & age<54)%>%
  group_by(SubjectID)%>%
  add_count(name = "n_samples")%>%
  filter(n_samples>=2)
  
kids_with_8 <- stool_samples %>%
  filter(age>6 & age<10)

kids_with_84<- stool_samples %>%
  filter(age>82 & age<86)




# bakars new data ####

library(tidyr)
library(dplyr)
library(xlsx)

`%notin%` <- Negate(`%in%`)
is.blank <- function(x){sapply(x, function(y) {ifelse(y=="", TRUE, FALSE)})}

raw_data <- readxl::read_excel("~/postdoc/stanford/clinical_data/MICDROP/microbiome/LIST OF STOOL SAMPLES ALIQUOT 1.xlsx")

stool_samples%>%
  mutate(Date = as.Date(lubridate::parse_date_time(.$Date, orders="%m/%d/%y")))%>%
  bind_rows(., stool_samples2)%>%
  select(SubjectID, Date, RandomNumber, BoxNumber, PositionRow, PositionColumn)

bdays <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2024_06/MICDSpecimenBoxJun24_withclinical.dta", col_select = c(id, dob))
bdays <- bdays[!duplicated(bdays),]
# 
raw_data$dob <- bdays$dob[match(raw_data$SubjectID, bdays$id)]

mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")


sample_ages <- c(8, 24, 52, 68, 84, 104, 120, 136, 156, 172, 188, 208)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1

sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))


stool_samples2 <- raw_data %>%
  mutate(dob= bdays$dob[match(SubjectID, bdays$id)])%>%
  mutate("age_in_weeks"=as.numeric(as.Date(Date)-dob)%/%7)%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(SubjectID, mic_drop_key$id)],
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))



kids_with_8_24_52 <- stool_samples %>%
  filter(age_in_weeks<60)%>%
  group_by(SubjectID)%>%
  add_count(name = "n_samples")%>%
  filter(n_samples>=3)

kids_with_24_and_52 <- stool_samples %>%
  filter(age_in_weeks>20 & age_in_weeks<60)%>%
  group_by(SubjectID)%>%
  add_count(name = "n_samples")%>%
  filter(n_samples>=2)

kids_with_6_12_24 <- stool_samples %>%
  filter(age_in_weeks>22 & age_in_weeks<106)%>%
  group_by(SubjectID)%>%
  add_count(name = "n_samples")%>%
  filter(n_samples>=3)


kids_with_two_or_more <- stool_samples %>%
  group_by(SubjectID)%>%
  add_count(name = "n_samples")%>%
  filter(n_samples>=2)

table(stool_samples1$RandomNumber %in% stool_samples2)

stool_samples%>%
  mutate(Date = as.Date(lubridate::parse_date_time(.$Date, orders="%m/%d/%y")))%>%
  bind_rows(., stool_samples2)

  