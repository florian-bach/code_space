library(tidyr)
library(dplyr)
library(xlsx)

`%notin%` <- Negate(`%in%`)
is.blank <- function(x){sapply(x, function(y) {ifelse(y=="", TRUE, FALSE)})}

raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/May25/MICDSpecimenBoxMay25_withclinical.dta")
# bdays <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_10/MICDSpecimenBoxDec23_withclinical.dta", col_select = c(id, dob))
# bdays <- bdays[!duplicated(bdays),]
# 
# raw_data$dob <- bdays$id[match(raw_data$id, bdays$id)]

sample_ages <- c(8, 24, 52, 68, 84, 104, 120, 136, 156, 172, 188, 208)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1

sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))

long_specimen_data <- raw_data %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7)%>%
  select(id, dob, date, flo_age_in_wks, SampleDate, PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR, visittype) %>%
  mutate("visit_id"=paste(id, date, sep="_"))%>%
  # mutate(stool=case_when(stool==1~paste(id, "poop", date, sep="_"),
  #                        stool==0~"",
  #                        stool==NA~""))%>%
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


routine_ish_visits <- long_specimen_data %>%
  filter(flo_age_in_wks %in% sample_ranges)


# how many individuals have attended each timepoint where samples were taken;
# sometimes samples are taken but visit-type is not recorded
# the 12 extra visits exist because the enrolment visit was like a day before the routine visit; why?
participant_summary <- long_specimen_data %>%
  # mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
  filter(visittype %in% c(0, 1))%>%
  group_by(Timepoint_in_weeks) %>%
  summarise("Number_of_Individuals"=n_distinct(subject_id), "Number_of_Visits"=n_distinct(visit_id))


#correct for withdrawal
withdrawn <- raw_data %>%
  filter(!is.na(raw_data$withdrawaldate))%>%
  mutate("withdrawal_age" = withdrawaldate-dob,
         "withdrawal_age_in_wks" = round(as.integer(withdrawal_age)/7))%>%
  dplyr::select(id, withdrawal_age, withdrawal_age_in_wks)%>%
  distinct()

withdrawn_before_8 <- withdrawn%>%
  filter(withdrawal_age_in_wks<10)
# ten people, exactly how many 8 week visits are missing

visits_with_samples_taken <- routine_ish_visits %>%
  filter(Specimen_ID != "") %>%
  group_by(Timepoint_in_weeks) %>%
  summarise("Number_of_Individuals"=n_distinct(subject_id),
            "Number_of_Visits"=n_distinct(visit_id))


# routine_ish_visits %>%
#   filter(Timepoint_in_weeks %in% sample_ages)%>%
#   group_by(Specimen_Type, Timepoint_in_weeks) %>%
#   summarise("Number_of_Samples"=n_distinct(Specimen_ID_ID))

# number of samples of each type at each timepoint, including collection rate
sample_counts <- long_specimen_data %>%
  filter(Timepoint_in_weeks %in% sample_ages, Specimen_ID != "", visittype %in% c(0,1))%>%
  group_by(Specimen_Type, Timepoint_in_weeks) %>%
  summarise("Number_of_Samples"=n_distinct(Specimen_ID))%>%
  ungroup()%>%
  mutate("Number_of_Participants"=participant_summary$Number_of_Individuals[match(.$Timepoint_in_weeks, participant_summary$Timepoint_in_weeks)],
         "Collection_Rate"=.$Number_of_Samples/participant_summary$Number_of_Individuals[match(.$Timepoint_in_weeks, participant_summary$Timepoint_in_weeks)],
         "Collection_Failure_Rate"=1-Collection_Rate,
         "Missed_Samples"=participant_summary$Number_of_Individuals[match(Timepoint_in_weeks, participant_summary$Timepoint_in_weeks)]-Number_of_Samples)


collection_rate_plot <-   sample_counts%>%
  filter(Specimen_Type %in% c("CellStabilizer", "Paxgene", "PBMC", "Plasma", "Stool"))%>%
  ggplot(aes(x=factor(Timepoint_in_weeks), y=Collection_Rate, fill=Specimen_Type))+
  geom_bar(stat="identity")+
  facet_grid(~Specimen_Type)+
  geom_text(aes(y=1.07, label= paste0("frac(",Number_of_Samples, ",", Number_of_Participants,")")),parse = TRUE, size=3)+
  scale_y_continuous(labels = scales::label_percent(), breaks = seq(0,1,0.2), limits=c(0,1.11))+
  scale_fill_manual(values = colorspace::sequential_hcl(n=7, "RdPu")[1:6])+
  xlab("Timepoint (weeks)")+
  ylab("Collection Rate")+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2025_05/collection_rate_plot.png", collection_rate_plot, bg="white", dpi=444, width=10, height=4)



n_sample_summary <- routine_ish_visits %>%
  filter(Specimen_ID != "")%>%
  group_by(Specimen_Type, subject_id) %>%
  summarise("Number_of_Samples_of_Individual" = n()) %>%
  ungroup()%>%
  group_by(Specimen_Type,  Number_of_Samples_of_Individual) %>%
  summarise("Number_of_Individuals_with_n_Samples"=n())


tight_n_sample_summary_plot  <- n_sample_summary %>%
  filter(Specimen_Type %in% c("CellStabilizer", "Paxgene", "PBMC", "Plasma"))%>%
  ggplot(aes(x=factor(Number_of_Samples_of_Individual), y=Number_of_Individuals_with_n_Samples, fill=Specimen_Type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=Number_of_Individuals_with_n_Samples),vjust= -0.2, size=3)+
  facet_grid(~Specimen_Type, scales = "free_x")+
  scale_fill_manual(values = colorspace::sequential_hcl(n=7, "RdPu")[1:6])+
  xlab("Number of Samples of Individual")+
  ylab("Number of Individuals with n Samples")+
  theme_minimal()+
  theme(legend.position = "none",
        # axis.text.x = element_text(angle=90)
  )

ggsave("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2025_05/tight_n_sample_summary_plot.png", tight_n_sample_summary_plot, bg="white", dpi=444, width=10, height=4)
