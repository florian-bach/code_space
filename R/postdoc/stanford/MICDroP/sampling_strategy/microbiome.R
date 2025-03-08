mic_drop <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2024_09/MICDROP expanded database through September 30th 2024.dta")

dobs <- mic_drop%>%
  distinct(id, dob)

stool_samples <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/microbiome/micdrop_stool_samples_dec_2025.csv")

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
