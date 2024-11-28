basset_list <- readxl::read_xlsx("~/Downloads/TLRc_unique_cohort_ids.xlsx")

comp_pal <- c("asymptomatic"="lightgrey", "uncomplicated"="black", "complicated"="orange", "severe"="darkred")
mstatus_pal <- c("no malaria"="lightgrey", "uncomplicated"="black", "complicated"="orange", "quinine for AL failure"="violet", "Q/AS failure"="purple")

mic_drop <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2024_07/MICDROP expanded database through July 31st 2024.dta")

# merge parasitemia data so that qPCR takes precedent when both slide and qPCR are present
mic_drop <- mic_drop %>%
  mutate(mstatus = case_match(mstatus,
                              0~"no malaria",
                              1~"uncomplicated",
                              2~"complicated",
                              3~"quinine for AL failure",
                              4~"Q/AS failure"))%>%
  mutate(visit_id = paste(id, date, sep=""))%>%
  mutate("parasitaemia_method" = if_else(qPCRdich==1, "qPCR", if_else(BSdich==1, "smear", "dunno")))%>%
  mutate(any_parsdens = if_else(is.na(qPCRparsdens) & !is.na(pardens), pardens, qPCRparsdens))%>%
  mutate(parasitaemia_method = if_else(is.na(qPCRparsdens) & !is.na(pardens), "smear", parasitaemia_method))
# 
# basset_subset <- impact%>%
#   filter(id %in% (basset_list$Cohort_IDs+10000))
# 





impact <- haven::read_dta("~/postdoc/stanford/clinical_data/IMPACT/IMPACT expanded database through July 31st 2024.dta")



tlr_metadata <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/MICDroP Data/TLRC_SPECIMENS/MICDTLRC.DTA")

basset_samples <- c("VPRGI", "Q9A5I", "HAG6I", "YJXRI", "HAGGI")

sample_subset <- tlr_metadata %>%
  filter(tblBoxDetails_RandomNumber %in% basset_samples)%>%
  distinct(SubjectID, Date)

accidental_micdrop <- mic_drop %>%
  filter(id %in% (basset_list$Cohort_IDs-1000))%>%
  distinct(id, dob)

accidental_impact <- impact %>%
  filter(id %in% (basset_list$Cohort_IDs+10000))%>%
  distinct(id, dob)

