library(tidyr)
library(dplyr)
library(xlsx)

`%notin%` <- Negate(`%in%`)
is.blank <- function(x){sapply(x, function(y) {ifelse(y=="", TRUE, FALSE)})}

raw_data <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2024_06/MICDSpecimenBoxJun24_withclinical.dta")
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

routine_visits <- long_specimen_data %>%
  filter(visittype %in% c(0,1))

kids_with_full_one_year <- routine_visits %>%
  filter(Timepoint_in_weeks %in% c(8, 24, 52), Specimen_ID!="", Specimen_Type %in% c("Plasma", "PBMC", "Paxgene"))%>%
  group_by(subject_id, Specimen_Type)%>%
  filter(n()==3)%>%
  ungroup()%>%
  group_by(subject_id)%>%
  filter(n()==9)

write.csv(kids_with_full_one_year, "~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/kids_with_full_first_year_plasma_pbmc.csv", row.names = FALSE)

kids_with_full_year_and_neurocog <- kids_with_full_one_year %>%
  filter(subject_id %in% neuro_cog$ID)

write.csv(kids_with_full_year_and_neurocog, "~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/kids_with_full_first_year_plasma_pbmc_neurocog.csv", row.names = FALSE)




kids_with_full_set %>%
  filter(Specimen_Type == "Plasma")%>%
  group_by(subject_id)

kids_with_full_set %>%
  filter(Specimen_Type == "PBMC"))



# subset visit database to be only "duplicate" visits and send IDRC

# 
# # how many individuals do we have n samples available
# n_sample_summary <- routine_visits %>%
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


samples_missing_visits <- routine_visits %>%
  mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
  filter(Specimen_ID == "" & visittype==1, !duplicated(visit_id))


all_samples_missing_visits <- samples_missing_visits %>%
  group_by(visit_id)%>%
  summarise(id, "n_missing"=n(), Timepoint_in_weeks)%>%
  filter(n_missing==4, !duplicated(visit_id))

blank_visits <- table(all_samples_missing_visits$Timepoint_in_weeks)


# loose vs tight visit definition ####
## loose####
n_sample_summary <- routine_visits %>%
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
n_sample_summary <- routine_visits %>%
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
  routine_visits %>%
    filter(Specimen_ID%in%c("11101-AB", "11101-AA"))
)

View(
  routine_visits %>%
    filter(Specimen_ID%in%c("10680-AG", "10680-AA"))
)


View(
  routine_visits %>%
    filter(Specimen_ID%in%c("10561-AA"))
)




try <- routine_visits %>%
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

routine_visits %>%
  filter(visittype %in% c(0,1))%>%
  group_by(Timepoint_in_weeks)%>%
  filter(!duplicated(id))%>%
  summarise(n=n())



withdrawn_before_24 <- withdrawn%>%
  filter(withdrawal_age_in_wks<52)

withdrawn_at_24 <- withdrawn%>%
  filter(withdrawal_age_in_wks<25, withdrawal_age_in_wks>22)

withdrawn_at_52 <- withdrawn %>%
  filter(withdrawal_age_in_wks<25, withdrawal_age_in_wks>22)



try <- routine_visits %>%
  filter(id %in% withdrawn_at_24$id, , Specimen_ID!="")

try <- routine_visits %>%
  # mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
  # filter(visittype %in% c(0, 1))%>%
  group_by(Timepoint_in_weeks, id)%>%
  mutate("n_visits"=n_distinct(visit_id))%>%
  filter(n_visits>1, Timepoint_in_weeks==104)

#11 at week 104??

long_specimen_data %>%
  group_by(Timepoint_in_weeks)%>%
  summarise("n"=n_distinct(id))





tight <- routine_visits %>%
  filter(visittype%in%c(0,1), Specimen_ID != "")%>%
  group_by(visittype, Timepoint_in_weeks, Specimen_Type)%>%
  filter(!duplicated(visit_id))%>%
  summarise(n=n())


loose <- routine_visits %>%
  filter(Timepoint_in_weeks %in% sample_ages, Specimen_ID != "")%>%
  group_by(visittype, Timepoint_in_weeks, Specimen_Type)%>%
  filter(!duplicated(visit_id))%>%
  summarise(n=n())




# files for IDRC team ####
## doublet 8 week visits####

# this has actually been resolved, this duplicates are enrolment visits, labeled correctly, but counted incorrectly by me.

# multiple_8 <- routine_visits %>%
#   # mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
#   filter(visittype %in% c(0, 1), Timepoint_in_weeks==8, Specimen_Type=="PlasmaPK")%>%
#   group_by(Timepoint_in_weeks, id)%>%
#   filter(n_distinct(visit_id)>1)%>%
#   select(id, dob, date, Timepoint_in_weeks, visit_id)
# 
# write.csv(multiple_8, file = "~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_12/doublet_7_8_week_visits.csv", row.names = FALSE)

## stray 104 week visits ####
weird_104 <- routine_visits %>%
  # mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
  select(id, dob, date, Timepoint_in_weeks, visit_id)%>%
  filter(Timepoint_in_weeks==104, !duplicated(visit_id))

write.csv(weird_104, file = "~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_12/empty_104_week_visits.csv", row.names = FALSE)

## duplicate specimen ids####

duplicate_specimens <- routine_visits %>%
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


