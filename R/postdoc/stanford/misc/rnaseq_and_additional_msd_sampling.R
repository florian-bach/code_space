library(tidyr)
library(dplyr)
`%notin%`=Negate(`%in%`)

#BoxNumber1=PBMC
#BoxNumber2=Paxgene
#BoxNumber3=Plasma
#BoxNumber4=PlasmaPK
#BoxNumber5=CellStabiliser
#BoxNumber7=CellStabiliser

raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/May25/MICDSpecimenBoxMay25_withclinical.dta")

mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")

msd_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/batch_one.csv")

long_luminex <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/long_luminex.csv")%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks", "104 weeks")))

msd_samples <- msd_data%>%
  ungroup()%>%
  mutate(timepoint=paste(TimePt, "weeks"),
         id=SubjectID)%>%
  distinct(id, timepoint)%>%
  mutate(sample_id=paste(id, timepoint, sep="_"))

var_samples <- long_luminex%>%
  ungroup()%>%
  distinct(id, timepoint)%>%
  mutate(sample_id=paste(id, timepoint, sep="_"))

# actually locate paxgene samples ####

## read in databases
specimen_database <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/Jun25/MICDSpecimenBoxJun25_withclinical.dta")

sample_ages <- c(8, 24, 52, 68, 84, 104, 120)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1

sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))


#turn the data into long format, include a couple of convenience variables, drop irrelevant columns
long_specimen_data <- specimen_database %>%
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
  )



dobs <- raw_data%>%
  filter(!is.na(dob))%>%
  distinct(id, dob)


metadata_columns <- c("id", "dob", "date", "ageinwks", "gender", "mstatus", "qPCRparsdens", "visittype", "fever", "febrile", "rogerson", "anyHP", "GAcomputed", "gi", "SGA", "qPCRdich", "mqPCRparsdens")

#merge based on id and date
metadata <- raw_data%>%
  select(all_of(metadata_columns))%>%
  mutate(timepoint_num=as.numeric(ageinwks), id=as.numeric(id))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.01))

infs_and_meta <- specimen_database%>%
  mutate(dob2=dobs$dob[match(id, dobs$id)])%>%
  mutate("flo_age_in_wks"=as.numeric(date-dob2)%/%7)%>%
  # filter(ageinwks<54)%>%
  group_by(id) %>%
  mutate("total_n_para_12"=sum(pardens>1&flo_age_in_wks<53, na.rm=TRUE),
         "total_n_malaria_12"=sum(mstatus!=0&flo_age_in_wks<53, na.rm=TRUE),
         "total_n_para_24"=sum(pardens>1&flo_age_in_wks<105, na.rm=TRUE),
         "total_n_malaria_24"=sum(mstatus!=0&flo_age_in_wks<105, na.rm=TRUE),
         "total_n_para_12_24"=sum(pardens>1&flo_age_in_wks<105&flo_age_in_wks>52, na.rm=TRUE),
         "total_n_malaria_12_24"=sum(mstatus!=0&flo_age_in_wks<105&flo_age_in_wks>52, na.rm=TRUE),
         "total_n_para_6"=sum(pardens>1&flo_age_in_wks<27, na.rm=TRUE),
         "total_n_malaria_6"=sum(mstatus!=0&flo_age_in_wks<27, na.rm=TRUE),
         "any_malar_6"=if_else(total_n_malaria_6==0, FALSE, TRUE),
         "any_malar_12"=if_else(total_n_malaria_12==0, FALSE, TRUE))%>%
  select(id, date, total_n_para_6, total_n_malaria_6, total_n_para_12, total_n_malaria_12, total_n_para_24, total_n_malaria_24, total_n_para_12_24, total_n_malaria_12_24, -gender)%>%
  left_join(., metadata, by=c("id", "date"))%>%
  group_by(id)%>%
  mutate(total_n_para=max(total_n_para_12, na.rm = T),
         total_n_malaria=max(total_n_malaria_12, na.rm = T))%>%
  mutate(timepoint_num=case_when(timepoint_num==9~8,
                                 timepoint_num==25~24,
                                 timepoint_num==53~52,
                                 .default=timepoint_num))%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))

kids_with_comp <- raw_data%>%
  filter(mstatus==2)%>%
  pull(unique(id))
## subset database to include correct cohort of 100 individuals, 50:50 treatment:placebo, 8 & 52 weeks

# 70
n_distinct(na.omit(msd_samples$id))

#19
potential_msd <- long_luminex%>%
  distinct(id, total_n_para_12, total_n_malaria_12)%>%
  filter(id %notin% unique(msd_samples$id))%>%
  filter(total_n_para_12>=2 | total_n_malaria_12>1)%>%
  distinct(id)

msd_with_additional <- c(na.omit(unique(msd_samples$id)), potential_msd$id)

# 54 = placebo; 35 DP
mic_drop_key %>%
  filter(id %in% msd_with_additional)%>%
  group_by(treatmentarm)%>%
  summarise(n())

# add 11 DP kids from NULISA & var
potential_additional <- long_luminex%>%
  filter(treatmentarm=="DP 1 year")%>%
  filter(id %notin% msd_with_additional)%>%
  distinct(id, total_n_malaria_12, total_n_para_12_24, total_n_malaria_12_24)%>%
  arrange(desc(total_n_malaria_12_24))%>%
  slice_head(n = 11)

full_hundred <- c(msd_with_additional, potential_additional$id)

infs_and_meta%>%
  # filter(treatmentarm=="DP 1 year")%>%
  filter(id %in% full_hundred)%>%
  distinct(id, treatmentarm, total_n_para_12_24, total_n_malaria_12_24)%>%
  group_by(total_n_malaria_12_24, treatmentarm)%>%
  summarise(n=n())%>%
  pivot_wider(values_from = n, names_from = treatmentarm)


infs_and_meta%>%
  # filter(treatmentarm=="DP 1 year")%>%
  distinct(id, treatmentarm)%>%
  filter(id %in% c(msd_with_additional))%>%
  group_by(treatmentarm)%>%
  summarise(n())


long_luminex%>%
  # filter(treatmentarm=="DP 1 year")%>%
  filter(id %in% full_hundred)%>%
  distinct(id, treatmentarm)%>%
  group_by(treatmentarm)%>%
  summarise(n=n())%>%
  pivot_wider(values_from = n, names_from = treatmentarm)


paxgenes_to_pick <- long_specimen_data %>%
  filter(id %in% full_hundred)%>%
  filter(Timepoint_in_weeks %in% c(8, 52))%>%
  filter(Specimen_ID!="")%>%
  filter(Specimen_Type=="Paxgene")%>%
  select(id, date, RandomNumber2, BoxNumber2, PositionColumn1, PositionRow2)%>%
  arrange(BoxNumber2)
  
# write.csv(paxgenes_to_pick, "~/postdoc/stanford/rna_seq/micdrop/paxgenes_to_pick_aug25.csv", row.names = F)


# venn diagrams ####
library(ggvenn)


micdrop_sample_venn <- list("MSD"=unique(msd_samples$id),
                            "NULISA"=unique(clean_data$id),
                            "PfEMP1 serology"=unique(var_samples$id)
)

## what it is right now ####
# make data a list (of assyas) of lists (of individuals) 
(current_venn_diagram <- ggvenn(
  micdrop_sample_venn, 
  show_stats = "c",
  text_color = "black",
  text_size = 6,
  fill_color = c("darkred", "orange", "darkblue"),
  # stroke_size = 0.5, set_name_size = 4
)+
    ggtitle("current sampling")+
    theme(plot.title = element_text(hjust = 0.5, size=24)))




# what it could be ####

potential_micdrop_sample_venn <- list("MSD"=msd_with_additional,
                                      "RNAseq"=full_hundred,
                                      "PfEMP1 serology"=unique(var_samples$id)[unique(var_samples$id) %in% full_hundred]
                                      # "RNAseq"=unique(c(msd_samples$id, potential_msd$id))
)


# make data a list (of assyas) of lists (of individuals) 
(potential_venn_diagram <- ggvenn(
  potential_micdrop_sample_venn, 
  show_stats = "c",
  text_color = "black",
  text_size = 6,
  fill_color = c("darkred", "orange", "darkblue"),
  # stroke_size = 0.5, set_name_size = 4
)+
    ggtitle("potential sampling")+
    theme(plot.title = element_text(hjust = 0.5, size=24)))

treatment_group_distribution <- mic_drop_key%>%
  filter(id %in% c(msd_samples$id, potential_msd$id))%>%
  group_by(treatmentarm)%>%
  summarise(n())

# locate additional MSD plasma samples ####

long_specimen_data %>%
  filter(id %in% potential_msd$id)%>%
  filter(Timepoint_in_weeks %in% c(8, 24, 52))%>%
  filter(Specimen_ID!="")%>%
  filter(Specimen_Type=="Plasma")%>%
  select(id, date, RandomNumber3, BoxNumber3, PositionColumn3, PositionRow3)%>%
  arrange(BoxNumber3) -> additional_MSD_to_pick

write.csv(additional_MSD_to_pick, "~/postdoc/stanford/plasma_analytes/MICDROP/MSD/additional_samples_for_msd.csv", row.names = F)




potential_micdrop_sample_venn <- list("vaccination_and_viral_serology"=msd_with_additional,
                                      "whole_blood_bulk_RNAseq"=full_hundred,
                                      "Pf_EMP1_serology"=unique(var_samples$id)[unique(var_samples$id) %in% full_hundred],
                                      "plasma_proteomics"=unique(micdrop_nulisa_with_meta$id)
)

immuno_df_for_anna <- data.frame("id"=unlist(potential_micdrop_sample_venn))
immuno_df_for_anna$assay <- rownames(immuno_df_for_anna) 
immuno_df_for_anna$assay <- gsub("[0-9]", "", immuno_df_for_anna$assay)

existing_immunology_samples_for_anna.csv <- immuno_df_for_anna%>%
  mutate(value=TRUE)%>%
  pivot_wider(names_from = assay, values_from = value, values_fill=F)

write.csv(existing_immunology_samples_for_anna.csv,
          "~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/existing_immunology_samples_for_anna.csv",
          row.names = F)



sweden_map <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/swedes_for_plating.csv")

one_year_samples <- data.frame("sample_id"=rep(paste(sweden_map$id, sweden_map$date), each=2))

one_year_samples <- one_year_samples%>%
  mutate("plate"=rep_len(rep(1:9, each=80), nrow(.)),
         "plate row"=rep_len(rep(LETTERS[1:8],each=10), nrow(.)),
         "plate column"=rep_len(3:12, nrow(.)))

write.csv(one_year_samples, "~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/big_csp_wlisa_plate_map.csv", row.names = F)


## 2 year samples ####

two_year_samples <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/104_week_timepoints_to_pick.csv")
two_year_platemap <- data.frame("id"=rep(two_year_samples$id, each=2))%>%
  mutate("plate"=rep_len(rep(1:9, each=80), nrow(.)),
         "plate_row"=rep_len(rep(LETTERS[1:8],each=10), nrow(.)),
         "plate_column"=rep_len(3:12, nrow(.)))

write.csv(two_year_platemap, "~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/two_year_platemap.csv", row.names = F)


additional_MSD_to_pick=read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/additional_samples_for_msd.csv")

plate_map_for_george_box_4 <- additional_MSD_to_pick%>%
  arrange(id)%>%
  mutate("plate"="TS05021786",
         "plate row"=rep_len(rep(LETTERS[1:8],each=12), nrow(.)),
         "plate column"=rep_len(1:12, nrow(.)),
         "age_in_weeks"=rep_len(c("8 weeks", "24 weeks", "52 weeks"),nrow(.)))%>%
  select(id, age_in_weeks, plate, `plate row`, `plate column`)

write.csv(plate_map_for_george_box_4, "~/postdoc/stanford/plasma_analytes/MICDROP/MSD/additional_samples_platemap.csv", row.names = F)

