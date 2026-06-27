library(tidyr)
library(dplyr)

## infection data ####
raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/MICDroP Data/MICDROP all visit database through March 31st 2026.dta")

mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")
mic_drop_hbs <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/MICDROP SickleTr final.dta")


three_month_labels <- paste0(seq(0, 57, by=3), " to ", seq(3, 60, by=3), " months")
six_month_labels <- paste0(seq(0, 54, by=6), " to ", seq(6, 60, by=6), " months")
twelve_month_labels <- paste0(seq(0, 48, by=12), " to ", seq(12, 60, by=12), " months")


monthly_categories <- raw_data %>%
  filter(visittype!=0)%>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7,
         "flo_age_in_months"=as.numeric(date-dob)%/%30.5)%>%
  filter(!is.na(mstatus), !is.na(flo_age_in_wks), mstatus!=3, mstatus!=4)%>%
  mutate("four_week_increment"=floor(flo_age_in_months))%>%
  mutate(age_quarter=cut(flo_age_in_wks, breaks = seq(0,60,by=3), labels=three_month_labels))%>%
  mutate(age_semi=cut(flo_age_in_wks, breaks = seq(0,60,by=6), labels=six_month_labels))%>%
  arrange(id, four_week_increment)%>%
  group_by(id, four_week_increment)%>%
  mutate(outcome=case_when(any(mstatus==2)~"complicated",
                           any(mstatus==1)~"uncomplicated",
                           #any(mstatus==3)~"treatment_failure",
                           all(mstatus==0)&any(qPCRparsdens>0, pardens>0)~"asymptomatic",
                           all(mstatus==0)&all(qPCRparsdens==0|is.na(qPCRparsdens), pardens==0)~"uninfected",
                           all(mstatus==0)&all(qPCRparsdens==0|is.na(qPCRparsdens), pardens==0|is.na(qPCRparsdens))~"undetermined",
                           .default="something_else"))%>%
  group_by(id, four_week_increment)%>%
  mutate(has_comp=any(mstatus==2))%>%
  mutate(has_uncomp=any(mstatus==1))%>%
  mutate(has_hyper = ifelse(any(qPCRparsdens>100000), "hyperparasitemia", "normal"))%>%
  filter(
    (outcome=="complicated" & mstatus == 2) |
      (outcome=="uncomplicated" & mstatus == 1) |
      outcome%in%c("uninfected", "asymptomatic", "undetermined")
    #(has_hyper=="hyperparasitemia" & pardens==max(pardens))
  )%>%
  group_by(id, four_week_increment)%>%
  slice_max(n = 1, order_by = pardens, with_ties = F)%>%
  group_by(id, outcome)%>%
  arrange(date)%>%
  mutate(order_of_outcome = seq(1, n()))

## count instances of each outcome for each individual ####
long_monthly_categories <- monthly_categories%>%
  tidyr::pivot_wider(names_from = outcome, values_fill = 0,values_from = order_of_outcome,
                     id_cols = c("id", "date", "flo_age_in_wks", "flo_age_in_months", "age_quarter", "age_semi", "four_week_increment", "mstatus", "pardens", "qPCRparsdens", "has_hyper"))%>%
  
  ungroup()%>%
  arrange(id, date)%>%
  group_by(id)%>%
  #these n are not the order at the time, but the number of preceding episodes
  mutate(n_nothing=cummax(uninfected))%>%
  mutate(n_complicated=cummax(complicated))%>%
  mutate(n_uncomplicated=cummax(uncomplicated))%>%
  mutate(n_asymptomatic=cummax(asymptomatic))%>%
  mutate(n_any_para=n_uncomplicated+n_asymptomatic+n_complicated)%>%
  mutate(n_any_malaria=n_uncomplicated+n_complicated)%>%
  mutate(age_at_first_malaria=min(flo_age_in_months[n_any_malaria==1]))%>%
  mutate(age_at_first_para=min(flo_age_in_months[n_any_para==1]))



# extract outcome incidence at pre-defined 
outcomes <- long_monthly_categories%>%
  group_by(id)%>%
  # filter(flo_age_in_wks < 106)%>%
  # filter(four_wee k_increment %in% c(4, 6, 13, 19, 26))%>%
  reframe(
          "total_n_para_6"=max(n_any_para[flo_age_in_wks>=20&flo_age_in_wks<=28]),
          "total_n_malaria_6"=max(n_any_malaria[flo_age_in_wks>=20&flo_age_in_wks<=28]),
          "total_n_para_12"=max(n_any_para[flo_age_in_wks>=48&flo_age_in_wks<=56]),
          "total_n_malaria_12"=max(n_any_malaria[flo_age_in_wks>=48&flo_age_in_wks<=56]),
          "total_n_para_18"=max(n_any_para[flo_age_in_wks>=74&flo_age_in_wks<=82]),
          "total_n_malaria_18"=max(n_any_malaria[flo_age_in_wks>=74&flo_age_in_wks<=82]),
          "total_n_para_24"=max(n_any_para[flo_age_in_wks>=100&flo_age_in_wks<=108]),
          "total_n_malaria_24"=max(n_any_malaria[flo_age_in_wks>=100&flo_age_in_wks<=108]),
          
          "total_n_para_36"=max(n_any_para[flo_age_in_wks>=152&flo_age_in_wks<=160]),
          "total_n_malaria_36"=max(n_any_malaria[flo_age_in_wks>=152&flo_age_in_wks<=160]),
          
          
          "total_n_para_12_24"=total_n_para_24-total_n_para_12,
          "total_n_malaria_12_24"=total_n_malaria_24-total_n_malaria_12,
          "total_n_para_24_36"=total_n_para_36-total_n_para_24,
          "total_n_malaria_24_36"=total_n_malaria_36-total_n_malaria_24,
          
          "total_n_para_12_18"=total_n_para_18-total_n_para_12,
          "total_n_malaria_12_18"=total_n_malaria_18-total_n_malaria_12,
          "total_n_para_12_24"=total_n_para_24-total_n_para_12,
          "total_n_malaria_12_24"=total_n_malaria_24-total_n_malaria_12,
          "symp_prop_6"=total_n_malaria_6/total_n_para_6,
          "symp_prop_12"=total_n_malaria_12/total_n_para_12,
          "symp_prop_24"=total_n_malaria_24/total_n_para_24,
          "symp_prop_12_24"=total_n_malaria_12_24/total_n_para_12_24,
          "any_malar_6"=if_else(total_n_malaria_6==0, 0, 1),
          "any_malar_12"=if_else(total_n_malaria_12==0, 0, 1),
          "any_para_6"=if_else(total_n_para_6==0, 0, 1),
          "any_para_12"=if_else(total_n_para_12==0, 0, 1))%>%
  mutate(across(everything(), ~ ifelse(is.infinite(.)|is.nan(.), NA, .)))





# static metadata: demographics, genotypes, birth outcomes etc ####
meta_raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/Oct25/MICDSpecimenBoxOct25_withclinical.dta")
maternal_treatment_arms <- haven::read_dta("~/Library/CloudStorage/Box-Box/DP+SP study/Databases and preliminary findings/Final database used for analyses/DPSP treatment allocation_FINAL.dta")

meta <- meta_raw_data%>%
  distinct(id, dob, withdrawaldate, gender, rogerson, anyHP, gravid, gravidcat,
           preterm, birthweight, LBWdich, GAcomputed, SGA)%>%
  filter(!is.na(LBWdich))%>%
  mutate(withdrawalage=(withdrawaldate-dob))%>%
  select(-withdrawaldate, -dob)%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"),
         mom_rx=maternal_treatment_arms$treatmentarm[match(id-10000, maternal_treatment_arms$id)],
         mom_rx=case_match(mom_rx,
                           1~"SP",
                           2~"DP",
                           3~"DPSP"))%>%
  mutate(hbs=mic_drop_hbs$HbS[match(as.numeric(id), mic_drop_hbs$id)],
         hbs=case_match(hbs,
                        1~"HbAA",
                        2~"HbAS",
                        3~"HbSS"))

infs_and_meta <- meta%>%
  left_join(., outcomes, by="id")

write.csv(infs_and_meta, "~/postdoc/stanford/clinical_data/MICDROP/micdrop_metadata_for_immunology.csv", row.names = F)

