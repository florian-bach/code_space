library(tidyr)
library(dplyr)
library(xlsx)
library(ggplot2)

`%notin%` <- Negate(`%in%`)
is.blank <- function(x){sapply(x, function(y) {ifelse(y=="", TRUE, FALSE)})}

# all except lavstsen plasmas were done with this
# raw_data <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2024_06/MICDSpecimenBoxJun24_withclinical.dta")


raw_data <-  haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/Mar25/MICDSpecimenBoxMar25_withclinical.dta")


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
  )
# finding PBMC dropout samples ####
# 461 individuals, 646 samples
pbmc_samples <- long_specimen_data %>%
  filter(BoxNumber1 %in% paste("MICD-", c(3001:3013, 3015:3040), sep=""))
  # filter(Specimen_Type =="PBMC") #%>%
  # select(id, date, Specimen_Type, ageinwks, BoxNumber1)

# 49 dropouts
dropout <- long_specimen_data %>%
  filter(!is.na(withdrawaldate))%>%
  mutate(age_at_withdrawal=withdrawaldate-dob)%>%
  select(id, withdrawaldate, age_at_withdrawal)%>%
  filter(!duplicated(id))

# 31 individuals; 43 tubes; 27 8 weeks; 15 24 weeks
pbmc_dropouts <- pbmc_samples%>%
  filter(id %in% dropout$id)

old_pbmc_dropouts <- pbmc_samples%>%
  filter(id %in% dropout$id, flo_age_in_wks>10)

kids_with_clinical <- long_specimen_data %>%
  filter(id %in% pbmc_samples$id, mstatus!=0 & ageinwks<24) %>%
  select(id, date, Specimen_Type, ageinwks, BoxNumber1, mstatus, withdrawaldate)

# finding "ex vivo" samples ####
ex_vivo <- long_specimen_data %>%
  filter(!is.na(mstatus))%>%
  group_by(id)%>%
  # make variable "days_since_malaria" that calculates the days a malaria episode and resets to 0 when mstatus!=0
  mutate(malaria_episode_num_date = if_else(mstatus != 0, as.numeric(date), -Inf),
         # malaria_episode_dich = if_else(mstatus != 0, 1, 0))%>%
         most_recent_malaria = as.Date(cummax(malaria_episode_num_date)),
         malaria_episode_id = paste(id, most_recent_malaria,sep="_"),
         days_since_malaria = date-most_recent_malaria)%>%
  # filter(Specimen_Type=="PBMC", Specimen_ID!="", days_since_malaria > 3 & days_since_malaria<14)%>%
  # filter(Specimen_Type=="Plasma", Specimen_ID!="", days_since_malaria==0)%>%
  select(id, date, flo_age_in_wks, malaria_episode_id, mstatus, Specimen_ID, days_since_malaria)
  

# 
# n_inf_summary <- ex_vivo_pbmc %>%
#   group_by(id)%>%
#   summarise("n_inf"=max(n_infection))


mic_drop %>%
  filter(id %in% unique(pbmc_dropouts$id),
         !is.na(any_parsdens),
         # ageinwks < 25
         )%>%
  ggplot(., aes(x=date, y=as.numeric(any_parsdens)+0.001))+
  geom_point(aes(color=factor(mstatus, #levels=c("0",
                              #       "1",
                              #      "2",
                              #     "3")
  )))+
  geom_line(alpha=0.3, aes(group=id))+
  # ggrepel::geom_text_repel(data=label_df, aes(x=AGE, y=as.numeric(any_parsdens), label=age_in_days))+
  facet_wrap(~ id)+
  ylab("qPCR parasites / μl\n")+
  xlab("Date")+
  scale_y_log10(breaks=c(1/100, 1, 10^2, 10^4, 10^6))+
  scale_x_date(breaks="2 week"
               # labels = seq(0,24,by=8)
  )+
  theme_minimal()+
  geom_vline(xintercept = 24)+
  # scale_shape_manual(values=c(16,15))+
  # scale_color_manual(values=mstatus_pal)+
  guides(color=guide_legend(title=""))+
  theme(axis.text.x = element_text(size=5, angle=90, vjust=0.5))




mic_drop %>%
  filter(id %in% unique(old_pbmc_dropouts$id),
         !is.na(any_parsdens),
         # ageinwks < 25
  )%>%
  ggplot(., aes(x=ageinwks, y=as.numeric(any_parsdens)+0.001))+
  geom_point(aes(color=factor(mstatus, #levels=c("0",
                              #       "1",
                              #      "2",
                              #     "3")
  )))+
  geom_line(alpha=0.3, aes(group=id))+
  # ggrepel::geom_text_repel(data=label_df, aes(x=AGE, y=as.numeric(any_parsdens), label=age_in_days))+
  facet_wrap(~ id)+
  ylab("qPCR parasites / μl\n")+
  xlab("Date")+
  scale_y_log10(breaks=c(1/100, 1, 10^2, 10^4, 10^6))+
  theme_minimal()+
  geom_vline(xintercept = 24)+
  geom_vline(xintercept = 52, linetype="dashed")+
  # scale_shape_manual(values=c(16,15))+
  # scale_color_manual(values=mstatus_pal)+
  guides(color=guide_legend(title=""))+
  theme(axis.text.x = element_text(size=5, angle=90, vjust=0.5))


ggsave("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2023_07/figures/21_days_symptoms.png", two_week_symptoms_plot, width = 24, height=8, bg="white")



# fei gao hla typing ####
# need 100 24 week samples where 8 week samples exist
mic_drop_clin <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2024_09/MICDROP expanded database through September 30th 2024.dta")

kids_with_comp <- mic_drop_clin %>%
  filter(mstatus==2)%>%
  group_by(id)%>%
  filter(!duplicated(date))

eight_week_pbmcs <- pbmc_samples %>%
  filter(Timepoint_in_weeks==8, Specimen_Type=="PBMC")%>%
  filter(!duplicated(id))

set.seed(123)
samples_for_fei <- pbmc_samples %>%
  filter(Timepoint_in_weeks==24, qPCRparsdens==0, Specimen_Type=="PBMC", id %in% eight_week_pbmcs$id, id %notin% kids_with_comp)%>%
  slice_sample(n=100)%>%
  select(id, BoxNumber1, PositionColumn1, PositionRow1)

write.csv(samples_for_fei,"~/postdoc/stanford/clinical_data/MICDROP/samples_for_hla.csv")  

  
forty_more <- pbmc_samples %>%
  filter(Timepoint_in_weeks==24, qPCRparsdens==0, Specimen_Type=="PBMC", id %in% eight_week_pbmcs$id, id %notin% kids_with_comp, id %notin%samples_for_fei$id)%>%
  slice_sample(n=40)%>%
  select(id, BoxNumber1, PositionColumn1, PositionRow1)

write.csv(forty_more,"~/postdoc/stanford/clinical_data/MICDROP/forty_more.csv")  

# validation of DRB11 hits against neurocognitive testing ####
library(tidyr)
library(dplyr)


neuro_cog <- readxl::read_excel("~/postdoc/stanford/clinical_data/MICDROP/neurocog_19092024.xlsx")
neuro_cog_ids <- neuro_cog$ID

hla_data <- readxl::read_excel("~/postdoc/stanford/clinical_data/MICDROP/HLA typing results_Pras&Davis lab_BCGinfant.xlsx")

drb11 <- hla_data%>%
  mutate(`Sample ID2`=lag(`Sample ID`))%>%
  mutate(`Sample ID3`=ifelse(is.na(`Sample ID`), `Sample ID2`, `Sample ID`))%>%
  filter(hla_data$`IMGT/DRB1`=="11:01:02:--")%>%
  select(`Sample ID3`)

drb11_ids <- substr(drb11$`Sample ID3`, 10,15)
drb11_ids <- gsub(" ", "", drb11_ids)
drb11_ids <- unique(c(drb11_ids, 11451, 11455))


sample_ages <- c(8, 24, 52, 68, 84, 104, 120)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1

sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))

mic_drop <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2024_06/MICDSpecimenBoxJun24_withclinical.dta")%>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7)%>%
  mutate("Timepoint_in_weeks"=if_else(
    flo_age_in_wks %in% sample_ages, flo_age_in_wks, ifelse(
      flo_age_in_wks %in% sample_ages_minus, flo_age_in_wks+1, if_else(
        flo_age_in_wks %in% sample_ages_plus, flo_age_in_wks-1, 999)))
  )
  

table(drb11_ids %in% neuro_cog_ids)
drb11_ids[!drb11_ids%in%neuro_cog_ids]


more_for_fei <- mic_drop %>%
  filter(id %in% drb11_ids[!drb11_ids%in%neuro_cog_ids], Timepoint_in_weeks==8)%>%
  select(id, Timepoint_in_weeks, RandomNumber1, BoxNumber1, PositionRow1, PositionColumn1)%>%
  filter(RandomNumber1!="")
  
write.csv(more_for_fei, "~/postdoc/stanford/clinical_data/MICDROP/fei_gao_mark_davis/8weeks_for_fei.csv")



# finding 10 paxgene and PBMC for Scott and Oliver
# let's find 12 month old samples of kids with multiple malaria episodes, but without full sample sets, as to not collide with the big MICDROP experiment

mstatus_pal <- c("no malaria"="lightgrey", "uncomplicated"="black", "complicated"="orange", "quinine for AL failure"="violet", "Q/AS failure"="purple")

routine_visits <- long_specimen_data %>%
  filter(visittype %in% c(0,1))

kids_with_full_one_year <- routine_visits %>%
  filter(Timepoint_in_weeks %in% c(8, 24, 52), Specimen_ID!="", Specimen_Type %in% c("Plasma", "PBMC", "Paxgene"))%>%
  group_by(subject_id, Specimen_Type)%>%
  filter(n()==3)%>%
  ungroup()%>%
  group_by(subject_id)%>%
  filter(n()==9)

kids_without_full_one_year <- routine_visits %>%
  filter(id %notin% kids_with_full_one_year$id)%>%
  select(id)%>%
  distinct()


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

counter <- mic_drop %>%
  select(id, date, AGE, dob, mstatus, any_parsdens, parasitaemia_method)%>%
  filter(mstatus!="no malaria", AGE<=52)%>%
  filter(id %in% kids_without_full_one_year$id)%>%
  group_by(id) %>%
  add_count(name="number_of_episodes_in_year1")%>%
  filter(number_of_episodes_in_year1>=3)

mic_drop %>%
  filter(id %in% counter$id, !is.na(mstatus), !is.na(any_parsdens))%>%
  ggplot(., aes(x=date, y=any_parsdens+0.01, group=id))+
  geom_line()+
  geom_point(aes(color=mstatus))+
  scale_y_log10()+
  facet_wrap(~id)+
  scale_color_manual(values = mstatus_pal)+
  theme_minimal()

samples_for_scott <- long_specimen_data %>%
  filter(id %in% counter$id)%>%
  filter(Timepoint_in_weeks %in% c(52),
         Specimen_ID!="", Specimen_Type %in% c("PBMC", "Paxgene"))%>%
  arrange(id, Specimen_Type)%>%
  select(id, Timepoint_in_weeks, BoxNumber1, PositionColumn1, PositionRow1, BoxNumber2, PositionColumn2, PositionRow2)

write.csv(samples_for_scott, "~/postdoc/stanford/clinical_data/MICDROP/oliver_wirz_scott_boyd/eleven_kids_with_3_or_more_malaria_episodes_but_excluded_from_big_study.csv", row.names = FALSE)

#BoxNumber1=PBMC
#BoxNumber2=Paxgene
#BoxNumber3=Plasma
#BoxNumber4=PlasmaPK
#BoxNumber5=CellStabiliser
#BoxNumber7=CellStabiliser
table(long_specimen_data$Specimen_Type[long_specimen_data$Specimen_ID!=""], long_specimen_data$BoxNumber7[long_specimen_data$Specimen_ID!=""]=="")

grant_and_comp <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/all_209_for_project.csv")

list_for_petter <- long_specimen_data %>%
  filter(id %in% grant_list$id | id %in% kids_with_comp$id)%>%
  filter(Specimen_Type=="CellStabilizer")%>%
  select(id, flo_age_in_wks, Specimen_ID, BoxNumber5, PositionColumn5, PositionRow5)%>%
  filter(Specimen_ID!="")
  
write.csv(list_for_petter, "~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/list_for_petter.csv", row.names = F)

list_for_kenneth <- long_specimen_data %>%
  filter(id %in% grants_list$id, Specimen_Type=="Plasma")%>%
  filter(Specimen_ID!="", flo_age_in_wks<54)%>%
  select(id, flo_age_in_wks, date)
  

write.csv(list_for_kenneth, "~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/list_for_kenneth.csv", row.names = F)


# 70 individuals MSD mesoscale ####
grant_list <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/msd/Samples selected for Florian Nov 7 2024.dta")
grant_and_comp <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/all_209_for_project.csv")

grant_list_subset <- grant_list%>%
  filter(subset==1)

samples_to_pick_for_mesoscale <- long_specimen_data %>%
  filter(id %in% grant_list_subset$id, Specimen_Type=="Plasma")%>%
  filter(Specimen_ID!="", flo_age_in_wks<54, flo_age_in_wks %in% sample_ranges)%>%
  select(id, date, Timepoint_in_weeks, BoxNumber3, PositionColumn3, PositionRow3)%>%
  arrange(BoxNumber3)%>%
  mutate("new column"=rep_len(1:9, 210), "new row"=rep_len(rep(1:9,each=9), 210))
# write.csv(samples_to_pick_for_mesoscale, "~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/msd/samples_to_pick_for_mesoscale.csv", row.names = F)

reorganise_samples <- samples_to_pick_for_mesoscale%>%
  mutate("new.box"=rep_len(rep(1:9, each=81), nrow(.)))%>%
  select(id, date, Timepoint_in_weeks, "new column", "new row", new.box)%>%
  arrange(id)

#write.csv(reorganise_samples, "~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/msd/reorganise_for_mesoscale.csv", row.names = F)
write.csv(reorganise_samples, "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/msd_map_edit.csv", row.names = F)


plate_map <- reorganise_samples%>%
  mutate("plate_number"=rep_len(rep(c("TS05021694", "TS05021774", "TS05021687"), each=84), nrow(.)))%>%
  mutate("plate_column"=rep_len(rep(1:12, 7), nrow(.)))%>%
  mutate("plate_row"=rep_len(rep(LETTERS[1:7], each=12), nrow(.)))%>%
  select(id, date, Timepoint_in_weeks, plate_number, plate_column, plate_row)

write.csv(plate_map, "~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/msd/plate_map_for_mesoscale2.csv", row.names = F)


label_maker <- long_specimen_data %>%
  filter(id %in% grant_list_subset$id, Specimen_Type=="Plasma")%>%
  filter(Specimen_ID!="", flo_age_in_wks<54, flo_age_in_wks %in% sample_ranges)%>%
  select(id, date, RandomNumber3, BoxNumber3, PositionColumn3, PositionRow3)%>%
  arrange(BoxNumber3)%>%
  mutate("new column"=rep_len(1:9, 210), "new row"=rep_len(rep(1:9,each=9), 210))%>%
  mutate("new.box"=rep_len(rep(1:9, each=81), nrow(.)))%>%
  select(id, date, RandomNumber3)%>%
  arrange(id)

label_maker_with_space <- label_maker[rep(1:nrow(label_maker), each = 2), ]
label_maker_with_space$id <- as.character(label_maker_with_space$id)
label_maker_with_space$date <- as.character(label_maker_with_space$date)
label_maker_with_space[1:nrow(label_maker_with_space) %% 2 == 0, ] <- " "

write.csv(label_maker_with_space, "~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/msd/label_maker_with_space.csv", row.names = F)

grant_list_rest <- grant_list %>%
  filter(subset==0)

# grant_and_comp <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/all_209_for_project.csv")



write.csv(samples_to_pick_for_mesoscale, "~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/msd/samples_to_pick_for_mesoscale.csv", row.names = F)

# 200 individuals for Nulisa ####
grant_and_comp <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/all_209_for_project.csv")
samples_to_pick_for_mesoscale <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/msd/samples_to_pick_for_mesoscale.csv")

# micdrop_samples_to_pick_for_nulisa <- long_specimen_data %>%
#   # filter(id %in% grant_list$id | id %in% kids_with_comp$id, Specimen_Type=="Plasma")%>%
#   filter(id %in% grant_and_comp$id,  Specimen_Type=="Plasma")%>%
#   filter(Specimen_ID!="", flo_age_in_wks<54, flo_age_in_wks %in% sample_ranges)%>%
#   select(id, date, Timepoint_in_weeks)

rest_of_samples_to_pick <- long_specimen_data %>%
  filter(id %notin% samples_to_pick_for_mesoscale$id & id %in% grant_and_comp$id, Specimen_Type=="Plasma")%>%
  filter(Specimen_ID!="", flo_age_in_wks<54, flo_age_in_wks %in% sample_ranges)%>%
  select(id, flo_age_in_wks, date, BoxNumber3, PositionColumn3, PositionRow3)%>%
  arrange(BoxNumber3)%>%
  mutate("new column"=rep_len(1:9, nrow(.)),
         "new row"=rep_len(rep(1:9,each=9), nrow(.)),
         "new box"=rep_len(rep(1:9, each=81), nrow(.)))
  
write.csv(rest_of_samples_to_pick, "~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/rest_of_samples_to_pick_for_NULISA.csv", row.names = F)

reordered_rest_of_samples_to_pick <- rest_of_samples_to_pick%>%
  select(id, date, flo_age_in_wks, `new column`, `new row`, `new box`)%>%
  arrange(id)%>%
  mutate("reorder column"=rep_len(1:9, nrow(.)),
          "reorder row"=rep_len(rep(1:9,each=9), nrow(.)),
          "reorder box"=rep_len(rep(1:9, each=81), nrow(.)),
         "nulisa plate"=rep_len(rep(4:9, each=84), nrow(.)))

# write.csv(reordered_rest_of_samples_to_pick, "~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/reordered_rest_of_samples_to_pick_for_NULISA.csv", row.names = F)
write.csv(reordered_rest_of_samples_to_pick, "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/non_msd_edit.csv", row.names = F)



label_maker <- long_specimen_data %>%
  filter(id %notin% samples_to_pick_for_mesoscale$id & id %in% grant_and_comp$id, Specimen_Type=="Plasma")%>%
  filter(Specimen_ID!="", flo_age_in_wks<54, flo_age_in_wks %in% sample_ranges)%>%
  select(id, date, RandomNumber3)%>%
  arrange(id)

label_maker_with_space <- label_maker[rep(1:nrow(label_maker), each = 2), ]
label_maker_with_space$id <- as.character(label_maker_with_space$id)
label_maker_with_space$date <- as.character(label_maker_with_space$date)
label_maker_with_space[1:nrow(label_maker_with_space) %% 2 == 0, ] <- " "

write.csv(label_maker_with_space, "~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/rest_of_samples_label_maker_with_space.csv", row.names = F)



non_msd_plate_map <- rest_of_samples_to_pick%>%
  select(id, date, flo_age_in_wks)%>%
  add_row(id=10852, date=as.Date("2023-09-22"), flo_age_in_wks=68)%>%
  arrange(id)%>%
  mutate("plate_number"=rep_len(rep(4:12, each=84), nrow(.)))%>%
  mutate("plate_column"=rep_len(rep(1:12, 7), nrow(.)))%>%
  mutate("plate_row"=rep_len(rep(LETTERS[1:7], each=12), nrow(.)))

write.csv(non_msd_plate_map, "~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/plate_map_for_non_msd.csv", row.names = F)



# DPSP mums for Kattria####

neuro_cog <- readxl::read_excel("~/postdoc/stanford/clinical_data/MICDROP/neurocog_19092024.xlsx")
dpsp_clin <- haven::read_dta("~/postdoc/stanford/clinical_data/DPSP/DPSPSpecimenBoxOct23_withclinical.dta")
# dpsp_specimen <- readxl::read_excel("~/postdoc/stanford/clinical_data/DPSP/tblSpecimenDetails_immunology.xlsx")
dpsp_specimen <- haven::read_dta("~/Downloads/DPSPSpecimenBoxJul24_withclinical_plasmapbmc.dta")
#dpsp samples in stanfod

dpsp_samples <- readxl::read_excel("~/postdoc/stanford/clinical_data/DPSP/tblBoxDetails_immunology.xlsx")

new_db <- haven::read_dta("~/Downloads/DPSPSpecimenBoxJul24_withclinical_plasmapbmc.dta")

#144 kids
kids_for_kat <- neuro_cog$ID[neuro_cog$ID %in% grant_and_comp$id]
moms <- kids_for_kat-10000

small_dpsp_clin <- dpsp_clin%>%
  filter(id %in% moms, visittype %in% c(0, 3))

# samples Y8AQI, FYCGI have been thawed twice already.
samples_to_pick <- dpsp_specimen %>%
  filter(id %in% moms)%>%
  filter(RandomNumber9 %in% dpsp_samples$RandomNumber)%>%
  # filter(Aliquot%in% c(1,2))%>%
  select(id, date, RandomNumber9, BoxNumber9, PositionRow9, PositionColumn9)
  
kattria_sample_locations <- samples_to_pick%>%
  arrange(BoxNumber9)%>%
  mutate("new column"=rep_len(1:9, nrow(.)),
         "new row"=rep_len(rep(1:9,each=9), nrow(.)),
         "new box"=rep_len(rep(1:9, each=81), nrow(.)))

write.csv(kattria_sample_locations, "~/postdoc/stanford/clinical_data/DPSP/kattria_sample_locations.csv", row.names = F)


kattria_reorder_sample_locations <- kattria_sample_locations%>%
  mutate("column"=new.column, "row"=new.row, "box"=new.box)%>%
  select(id, date, column, row, box)%>%
  arrange(id)%>%
  mutate("new column"=rep_len(1:9, nrow(.)),
         "new row"=rep_len(rep(1:9,each=9), nrow(.)),
         "new box"=rep_len(rep(1:9, each=81), nrow(.)))

write.csv(kattria_reorder_sample_locations, "~/postdoc/stanford/clinical_data/DPSP/reorder_kattria_sample_locations.csv", row.names = F)

kattria_sample_locations <- read.csv("~/Downloads/kattria_sample_locations.csv")

kattria_label_maker <- kattria_sample_locations%>%
  filter()
  select(id, date, RandomNumber9)%>%
  arrange(id)

  
kattria_label_maker_with_space <- kattria_label_maker[rep(1:nrow(label_maker), each = 2), ]
kattria_label_maker_with_space$id <- as.character(kattria_label_maker_with_space$id)
kattria_label_maker_with_space$date <- as.character(kattria_label_maker_with_space$date)
kattria_label_maker_with_space[1:nrow(kattria_label_maker_with_space) %% 2 == 0, ] <- " "

write.csv(kattria_label_maker_with_space, "~/postdoc/stanford/clinical_data/DPSP/kattria_label_maker_with_space.csv", row.names = F)


write.csv(kattria_sample_locations, "~/Downloads/kattria_sample_locations.csv", row.names = F)


kattria_reorder_sample_locations <- read.csv("~/postdoc/stanford/clinical_data/DPSP/reorder_kattria_sample_locations.csv")

dpsp_platemap <- kattria_reorder_sample_locations %>%
  # filter(!(id==508 & date =="2022-01-27"))%>%
  select(id, date)%>%
  mutate("plate_number"=rep_len(rep(9:11, each=84), nrow(.)))%>%
  mutate("plate_column"=rep_len(rep(1:12, 7), nrow(.)))%>%
  mutate("plate_row"=rep_len(rep(LETTERS[1:7], each=12), nrow(.)))

write.csv(dpsp_platemap, "~/postdoc/stanford/clinical_data/DPSP/dpsp_plate_map.csv", row.names = F)


# DROPPED BOX 19

# box19 <- dpsp_specimen %>%
#   filter(BoxNumber9=="DPSP-1019")%>%
#   select(id, date, RandomNumber9, BoxNumber9, PositionRow9, PositionColumn9)


# extra 24 NULISA ####
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
  select(id, date, flo_age_in_wks, malaria_episode_id, mstatus, Specimen_ID, days_since_malaria, qPCRparsdens)

day0_plasmas_not_grant <- day0_plasmas%>%
  filter(id %notin% grant_and_comp$id)%>%
  filter(flo_age_in_wks < 53)

day0_plasmas_grant <- day0_plasmas%>%
  filter(id %in% grant_and_comp$id)%>%
  filter(flo_age_in_wks < 53)
#11343 has two, ~3 months apart
mic_drop <-  haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2024_09/MICDROP expanded database through September 30th 2024.dta")

mic_drop %>%
  filter(id %in% c(day0_plasmas_not_grant$id))%>%
  ggplot(aes(x=ageinwks, y=qPCRparsdens+0.01))+
  geom_line(alpha=0.3, aes(group=id))+
  facet_wrap(~id, scales = "free_x")+
  geom_vline(xintercept = c(8, 24, 52), linetype="dashed")+
  geom_point(aes(color=factor(mstatus)))+
  scale_y_log10()+
  theme_minimal()

mic_drop %>%
  filter(id %in% c(day0_plasmas_grant$id))%>%
  ggplot(aes(x=ageinwks, y=qPCRparsdens+0.01))+
  geom_line(alpha=0.3, aes(group=id))+
  facet_wrap(~id, scales = "free_x")+
  geom_vline(xintercept = c(8, 24, 52), linetype="dashed")+
  geom_point(aes(color=factor(mstatus)))+
  scale_y_log10()+
  theme_minimal()
  
mic_drop %>%
  filter(id %in% c(day0_plasmas_not_grant$id))%>%
  ggplot(aes(x=ageinwks, y=qPCRparsdens+0.01))+
  geom_line(alpha=0.3, aes(group=id))+
  geom_vline(xintercept = c(8, 24, 52, 68), linetype="dashed")+
  geom_point(aes(color=factor(mstatus)))+
  facet_wrap(~id, scales = "free_x")+
  scale_y_log10()+
  theme_minimal()

View(long_specimen_data %>%
  filter(id==11343, Specimen_ID!=""))

# individuals to include

# 11685; malaria at 8 weeks; parasitemic at 24
# 10501; malaria at 8 weeks; parasitemic at 24
# 10857; malaria at 8 weeks; mildly parasitemic at 52
# 11831; malaria at 8 weeks; no parasites after
# 11721; malaria at 24 weeks; parasitemic at 52
# 11462; malaria at 24 weeks; parasitemic at 52
# 10766; malaria at 24 weeks; parasitemic at 8 and 52 weeks
# 11343; malaria at 52 weeks; paraistemic at 24

#11651; malaria at 52 weeks; no parasites other visits
#11622; malaria at 24 weeks, no parasites ever again
#10950; malaria at 52 weeks, no paarasites other visits

extra_individuals <- c(11685, 10501, 10857, 11831, 11721, 11462, 11622, 11343)

longer_timecourses <- c(11343, 10857)

extra_samples <- long_specimen_data %>%
  dplyr::filter(id %in% extra_individuals) %>%
  dplyr::filter(flo_age_in_wks < 54 | id %in% longer_timecourses,
         Specimen_ID!="", Specimen_Type %in% c("Plasma"))%>%
  arrange(id, Specimen_Type)%>%
  select(RandomNumber3, id, Timepoint_in_weeks, BoxNumber3, PositionColumn3, PositionRow3)%>%
  arrange(BoxNumber3)

write.csv(extra_samples, "~/postdoc/stanford/plasma_analytes/MICDROP/extra_24/extra_24_sample_manifest_for_tran.csv", row.names = F)

#11343
#11685, 11831 only two samples each

# 10460 only one sample
# 10950 complete
# 11622 complete
# 11651 complete

# 100 for lavstsen

