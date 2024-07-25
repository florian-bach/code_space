library(tidyr)
library(dplyr)
library(xlsx)

`%notin%` <- Negate(`%in%`)
is.blank <- function(x){sapply(x, function(y) {ifelse(y=="", TRUE, FALSE)})}

raw_data <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_12/MICDSpecimenBoxDec23_withclinical.dta")
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


no_visit_type_recorded <- routine_ish_visits %>%
  filter(Specimen_ID != "" & visittype %notin% c(0, 1, 2))

n_distinct(no_visit_type_recorded$visit_id)
n_distinct(no_visit_type_recorded$id)
#776; 746


# subset routine visits with one or more samples missing
visits_missing_samples <- routine_ish_visits %>%
  filter(Specimen_ID == "" & visittype==1)

n_distinct(visits_missing_samples$visit_id)
n_distinct(visits_missing_samples$id)
#1632 visits; 854 children


# count the number of missing samples per visit, subset to only include visits where all are missing
visits_missing_all_samples <- visits_missing_samples %>%
  group_by(visit_id)%>%
  mutate("n_missing"=n())%>%
  filter(n_missing==4, !duplicated(visit_id))
# 45 visits miss all samples 


# how many individuals have attended each timepoint where samples were taken;
# sometimes samples are taken but visit-type is not recorded
# the 12 extra visits exist because the enrolment visit was like a day before the routine visit; why?
participant_summary <- routine_ish_visits %>%
  # mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
  filter(visittype %in% c(0, 1))%>%
  group_by(Timepoint_in_weeks) %>%
  summarise("Number_of_Individuals"=n_distinct(subject_id), "Number_of_Visits"=n_distinct(visit_id))



# n_distinct(multiple_8$id)
# [1] 12

long_specimen_data %>%
  group_by(Timepoint_in_weeks)%>%
  summarise("n"=n_distinct(id), "n_vis"=n_distinct(visit_id))

# Timepoint_in_weeks     n
# <dbl> <int>
#   1                  8   919
# 2                 24   901
# 3                 52   612
# 4                 68   436
# 5                 84   218
# 6                104     4
# 7                999   952 

# participants_with_no_8 <- raw_data$id[raw_data$id %notin% participants_with_8$id]

# n_distinct(participants_with_8$id)
# 919

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
sample_counts <- routine_ish_visits %>%
  filter(Timepoint_in_weeks %in% sample_ages, Specimen_ID != "")%>%
  group_by(Specimen_Type, Timepoint_in_weeks) %>%
  summarise("Number_of_Samples"=n_distinct(Specimen_ID))%>%
  ungroup()%>%
  mutate("Collection_Rate"=.$Number_of_Samples/participant_summary$Number_of_Individuals[match(.$Timepoint_in_weeks, participant_summary$Timepoint_in_weeks)],
         "Collection_Failure_Rate"=1-Collection_Rate,
         "Missed_Samples"=participant_summary$Number_of_Individuals[match(Timepoint_in_weeks, participant_summary$Timepoint_in_weeks)]-Number_of_Samples)


collection_rate_plot <-   sample_counts%>%
  filter(Specimen_Type %in% c("CellStabilizer", "Paxgene", "PBMC", "Plasma"))%>%
  ggplot(aes(x=factor(Timepoint_in_weeks), y=Collection_Rate, fill=Specimen_Type))+
  geom_bar(stat="identity")+
  facet_grid(~Specimen_Type)+
  scale_y_continuous(labels = scales::label_percent())+
  scale_fill_manual(values = colorspace::sequential_hcl(n=7, "RdPu")[1:6])+
  xlab("Timepoint (weeks)")+
  ylab("Collection Rate")+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_12/collection_rate_plot.png", collection_rate_plot, bg="white", dpi=444, width=8, height=4)


indie_samples <- routine_ish_visits %>%
  filter(Timepoint_in_weeks %in% sample_ages, Specimen_ID!="")%>%
  group_by(id, Specimen_Type, Timepoint_in_weeks) %>%
  summarise("Number_of_Samples"=n_distinct(Specimen_ID))

table(indie_samples$Number_of_Samples, indie_samples$Specimen_Type)

routine_ish_visits %>%
  filter(id==11101, Timepoint_in_weeks==24)%>%
  select(id, date, Specimen_ID, Specimen_Type)%>%
  arrange(desc(Specimen_ID), date)

routine_ish_visits %>%
  filter(id==10680, Timepoint_in_weeks==52)%>%
  select(id, date, Specimen_ID, Specimen_Type)%>%
  arrange(desc(Specimen_ID), date)





# subset visit database to be only "duplicate" visits and send IDRC

# 
# # how many individuals do we have n samples available
# n_sample_summary <- routine_ish_visits %>%
#   filter(Specimen_ID != "")%>%
#   group_by(Specimen_Type, subject_id) %>%
#   summarise("Number_of_Samples_of_Individual" = n()) %>%
#   ungroup()%>%
#   group_by(Specimen_Type,  Number_of_Samples_of_Individual) %>%
#   summarise("Number_of_Individuals_with_n_Samples"=n())
# 
# 
# 
# n_sample_summary_plot  <- n_sample_summary %>%
#   filter(Specimen_Type %in% c("CellStabilizer", "Paxgene", "PBMC", "Plasma"))%>%
#   ggplot(aes(x=factor(Number_of_Samples_of_Individual), y=Number_of_Individuals_with_n_Samples, fill=Specimen_Type))+
#   geom_bar(stat="identity")+
#   geom_text(aes(label=Number_of_Individuals_with_n_Samples),vjust= -0.2, size=3)+
#   facet_grid(~Specimen_Type, scales = "free_x")+
#   scale_fill_manual(values = colorspace::sequential_hcl(n=7, "RdPu")[1:6])+
#   xlab("Number of Samples of Individual")+
#   ylab("Number of Individuals with n Samples")+
#   theme_minimal()+
#   theme(legend.position = "none",
#         # axis.text.x = element_text(angle=90)
#   )
# 
# ggsave("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_12/n_sample_summary_plot.png", n_sample_summary_plot, bg="white", dpi=444, width=8, height=4)
# 

# visits without samples


samples_missing_visits <- routine_ish_visits %>%
  mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
  filter(Specimen_ID == "" & visittype==1, !duplicated(visit_id))


all_samples_missing_visits <- samples_missing_visits %>%
  group_by(visit_id)%>%
  summarise(id, "n_missing"=n(), Timepoint_in_weeks)%>%
  filter(n_missing==4, !duplicated(visit_id))

blank_visits <- table(all_samples_missing_visits$Timepoint_in_weeks)


# loose vs tight visit definition ####
## loose####
n_sample_summary <- routine_ish_visits %>%
  filter(Specimen_ID != "")%>%
  group_by(Specimen_Type, subject_id) %>%
  summarise("Number_of_Samples_of_Individual" = n()) %>%
  ungroup()%>%
  group_by(Specimen_Type,  Number_of_Samples_of_Individual) %>%
  summarise("Number_of_Individuals_with_n_Samples"=n())



loose_n_sample_summary_plot  <- n_sample_summary %>%
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


## tight ####
n_sample_summary <- routine_ish_visits %>%
  filter(Specimen_ID != "", visittype%in%c(0,1))%>%
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

loose_vs_tight_plot = loose_n_sample_summary_plot + tight_n_sample_summary_plot


ggsave("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_12/loose_vs_tight_plot_n_sample_summary_plot.png", loose_vs_tight_plot, bg="white", dpi=444, width=12, height=4)



# sand box ####


View(
  routine_ish_visits %>%
    filter(Specimen_ID%in%c("11101-AB", "11101-AA"))
)

View(
  routine_ish_visits %>%
    filter(Specimen_ID%in%c("10680-AG", "10680-AA"))
)


View(
  routine_ish_visits %>%
    filter(Specimen_ID%in%c("10561-AA"))
)




try <- routine_ish_visits %>%
  filter(Specimen_ID!="")%>%
  group_by(id, Specimen_ID, Specimen_Type)%>%
  summarise(n=n(), Timepoint_in_weeks)

multiples <- try %>%
  filter(n>1)


hi <- long_specimen_data %>%
  filter(Specimen_ID != "")
table(hi$visittype, useNA = "ifany")

routine_visits <- long_specimen_data %>%
  filter(visittype==1)%>%
  group_by(Timepoint_in_weeks)%>%
  filter(!duplicated(id))%>%
  summarise(n=n())

table(routine_visits$Timepoint_in_weeks)

routine_ish_visits %>%
  filter(visittype %in% c(0,1))%>%
  group_by(Timepoint_in_weeks)%>%
  filter(!duplicated(id))%>%
  summarise(n=n())



withdrawn_before_24 <- withdrawn%>%
  filter(withdrawal_age_in_wks<52)

withdrawn_at_24 <- withdrawn%>%
  filter(withdrawal_age_in_wks<25, withdrawal_age_in_wks>22)


try <- routine_ish_visits %>%
  filter(id %in% withdrawn_at_24$id, , Specimen_ID!="")

try <- routine_ish_visits %>%
  # mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
  # filter(visittype %in% c(0, 1))%>%
  group_by(Timepoint_in_weeks, id)%>%
  mutate("n_visits"=n_distinct(visit_id))%>%
  filter(n_visits>1, Timepoint_in_weeks==104)

#11 at week 104??

long_specimen_data %>%
  group_by(Timepoint_in_weeks)%>%
  summarise("n"=n_distinct(id))





tight <- routine_ish_visits %>%
  filter(visittype%in%c(0,1), Specimen_ID != "")%>%
  group_by(visittype, Timepoint_in_weeks, Specimen_Type)%>%
  filter(!duplicated(visit_id))%>%
  summarise(n=n())


loose <- routine_ish_visits %>%
  filter(Timepoint_in_weeks %in% sample_ages, Specimen_ID != "")%>%
  group_by(visittype, Timepoint_in_weeks, Specimen_Type)%>%
  filter(!duplicated(visit_id))%>%
  summarise(n=n())




# files for IDRC team ####
  ## doublet 8 week visits####

# this has actually been resolved, this duplicates are enrolment visits, labeled correctly, but counted incorrectly by me.

# multiple_8 <- routine_ish_visits %>%
#   # mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
#   filter(visittype %in% c(0, 1), Timepoint_in_weeks==8, Specimen_Type=="PlasmaPK")%>%
#   group_by(Timepoint_in_weeks, id)%>%
#   filter(n_distinct(visit_id)>1)%>%
#   select(id, dob, date, Timepoint_in_weeks, visit_id)
# 
# write.csv(multiple_8, file = "~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_12/doublet_7_8_week_visits.csv", row.names = FALSE)

  ## stray 104 week visits ####
weird_104 <- routine_ish_visits %>%
  # mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
  select(id, dob, date, Timepoint_in_weeks, visit_id)%>%
  filter(Timepoint_in_weeks==104, !duplicated(visit_id))
  
write.csv(weird_104, file = "~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_12/empty_104_week_visits.csv", row.names = FALSE)

## duplicate specimen ids####

duplicate_specimens <- routine_ish_visits %>%
  filter(Specimen_ID!="")%>%
  group_by(id, Specimen_ID, Specimen_Type)%>%
  select(id)

duplicate_specimens <- routine_visits %>%
  filter(Specimen_ID!="")%>%
  group_by(id, Specimen_ID, Specimen_Type)%>%
  add_count(name = "n")%>%
  filter(n>1)%>%
  select(id, dob, date, visittype, Timepoint_in_weeks, visit_id)

write.csv(duplicate_specimens, file = "~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_12/duplicate_specimen_ids.csv", row.names = FALSE)


