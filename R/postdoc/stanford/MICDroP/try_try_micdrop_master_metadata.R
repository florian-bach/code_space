library(tidyr)
library(dplyr)

# infection data ####

## read in data ####
raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/MICDroP Data/MICDROP all visit database through October 31st 2025.dta")

mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")
mic_drop_hbs <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/MICDROP SickleTr final.dta")

three_month_labels <- paste0(seq(0, 57, by=3), " to ", seq(3, 60, by=3), " months")
six_month_labels <- paste0(seq(0, 54, by=6), " to ", seq(6, 60, by=6), " months")
twelve_month_labels <- paste0(seq(0, 48, by=12), " to ", seq(12, 60, by=12), " months")

## prevalent infections ####
# we count parasitemic months first and then count incidence malaria episodes in a different way;
# the reason for that is that multiple measurements of asymptomatic parasitemia within a four week window will be counted as a single
# parasitemic epiodes; in contrast, two malaria episodes, so long as they're more than 14 days apart post-treatment, will be counted
# as seperate episodes, in conflict with the parasitemia counting algorithm 

### calculate monthly categories for parasite growth ####
monthly_categories <- raw_data %>%
  filter(visittype!=0)%>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)/7,
         "flo_age_in_months"=as.numeric(date-dob)%/%30.5)%>%
  filter(!is.na(mstatus), !is.na(flo_age_in_wks), mstatus!=3, mstatus!=4)%>%
  mutate("four_week_increment"=floor(flo_age_in_wks/4))%>%
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

### count instances of each outcome for each individual ####
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



# extract outcome incidence at pre-defined four week increments
prevalent_outcomes <- long_monthly_categories%>%
  group_by(id)%>%
  # filter(flo_age_in_wks < 106)%>%
  # filter(four_wee k_increment %in% c(4, 6, 13, 19, 26))%>%
  reframe("total_n_para_6"=max(n_any_para[four_week_increment==5], na.rm = T),
          "total_n_para_12"=max(n_any_para[four_week_increment==12], na.rm = T),
          "total_n_para_24"=max(c(n_any_para[four_week_increment==25], n_any_para[four_week_increment==26]), na.rm = T),
          "total_n_para_36"=max(c(n_any_para[four_week_increment==37], n_any_para[four_week_increment==38]), na.rm = T),
          "total_n_para_12_24"=total_n_para_24-total_n_para_12,
          "total_n_para_24_36"=total_n_para_36-total_n_para_24)%>%
  # -5 = truncated data
  # -10 = impossible match (denominator is 0 / NA)
  mutate(across(starts_with("total_n"), ~ case_when(is.infinite(.)~ -5,
                                          is.nan(.)~ -10,
                                          .default = .)))
  


# incident infections 

all_malaria <- raw_data %>%
  filter(visittype!=0)%>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)/7,
         "flo_age_in_months"=as.numeric(date-dob)%/%30.5)%>%
  group_by(id)%>%
  mutate(max_age=max(flo_age_in_wks, na.rm = T))%>%
  filter(!is.na(mstatus), !is.na(flo_age_in_wks), mstatus!=3, mstatus!=4)%>%
  # sometimes randomly a lot of time passes so we "skip" a window..
  mutate("four_week_increment"=floor(flo_age_in_wks/4))%>%
  filter(mstatus%in%c(1,2))%>%
  select(id, date, flo_age_in_wks, four_week_increment, mstatus, pardens, qPCRparsdens, max_age)

# this misses all kids that never have any malaria, though we should count these somehow
incident_outcomes <- all_malaria%>%
  group_by(id)%>%
  reframe(
    "total_n_malaria_6"=ifelse(max_age>=26, sum(mstatus%in%c(1,2)&flo_age_in_wks<=26), NA),
    "total_n_malaria_12"=ifelse(max_age>=52, sum(mstatus%in%c(1,2)&flo_age_in_wks<=52), NA),
    "total_n_malaria_24"=ifelse(max_age>=104, sum(mstatus%in%c(1,2)&flo_age_in_wks<=104), NA),
    "total_n_malaria_36"=ifelse(max_age>=156, sum(mstatus%in%c(1,2)&flo_age_in_wks<=156), NA),
    "total_n_malaria_12_24"=ifelse(total_n_malaria_24==0&total_n_malaria_12==0, 0, total_n_malaria_24-total_n_malaria_12),
    "total_n_malaria_24_36"=ifelse(total_n_malaria_36==0&total_n_malaria_24==0, 0, total_n_malaria_36-total_n_malaria_24))%>%
  filter(!duplicated(id))

`%notin%`=Negate(`%in%`)

incident_and_prevalent_outcomes <- prevalent_outcomes%>%
  left_join(., incident_outcomes, by="id")%>%
  mutate("total_n_malaria_6"=ifelse(id %notin% incident_outcomes$id, 0, total_n_malaria_6),
         "total_n_malaria_12"=ifelse(id %notin% incident_outcomes$id & total_n_para_12%notin%c(-10, -5), 0, total_n_malaria_12),
         "total_n_malaria_24"=ifelse(id %notin% incident_outcomes$id & total_n_para_24%notin%c(-10, -5), 0, total_n_malaria_24),
         "total_n_malaria_36"=ifelse(id %notin% incident_outcomes$id & total_n_para_36%notin%c(-10, -5), 0, total_n_malaria_36),
         "total_n_malaria_12_24"=ifelse(id %notin% incident_outcomes$id & total_n_para_12_24%notin%c(-10, -5), 0, total_n_malaria_12_24),
         "total_n_malaria_24_36"=ifelse(id %notin% incident_outcomes$id & total_n_para_24_36%notin%c(-10, -5), 0, total_n_malaria_24_36))%>%
  mutate("symp_prop_6"=total_n_malaria_6/total_n_para_6,
         "symp_prop_12"=total_n_malaria_12/total_n_para_12,
         "symp_prop_24"=total_n_malaria_24/total_n_para_24,
         "symp_prob_36"=total_n_malaria_36/total_n_para_36,
         "symp_prop_12_24"=total_n_malaria_12_24/total_n_para_12_24,
         "symp_prop_24_36"=total_n_malaria_24_36/total_n_para_24_36,
         "any_malar_6"=if_else(total_n_malaria_6==0, 0, 1),
         "any_malar_12"=if_else(total_n_malaria_12==0, 0, 1),
         "any_para_6"=if_else(total_n_para_6==0, 0, 1),
         "any_para_12"=if_else(total_n_para_12==0, 0, 1))




long_monthly_categories%>%
  filter(four_week_increment==26)%>%
  ungroup()%>%
  summarise("mean"=mean(flo_age_in_wks), "min"=min(flo_age_in_wks), "max"=max(flo_age_in_wks))
# raw_data%>%
#        mutate(age=date-dob, withdrawal_age=withdrawaldate-dob)%>%
#        distinct(id, withdrawal_age)%>%
#        filter(withdrawal_age<1095)%>%
#        summarise(n())

#64 withdrawals under 12 months; 126 withdrawals under 24 months; 165 withdrawals under 36 months
# sapply(outcomes, function(x)sum(is.na(x)))
sapply(incident_and_prevalent_outcomes, function(x)sum(x%in%c(-10, -5)))
sapply(incident_and_prevalent_outcomes, function(x)sum(x%in%c(-5)))




sample_kids <- sample(unique(outcomes$id), 15)

year_three_raw_data <- raw_data%>%
  mutate(flo_age_in_wks = as.numeric(floor((date-dob)/7)))%>%
  filter(id %in% sample_kids, !is.na(mstatus))%>%
  select(id, flo_age_in_wks, mstatus, pardens, qPCRparsdens)


outcomes_df <- incident_and_prevalent_outcomes%>%
  filter(id %in% sample_kids)



year_three_raw_data <- year_three_raw_data%>%
  arrange(desc(mstatus))


ggplot()+
  geom_vline(xintercept = c(52, 104), linetype="dashed", color="grey")+
  geom_point(data=year_three_raw_data, aes(x=flo_age_in_wks, y=qPCRparsdens+0.001, color=factor(mstatus)), shape=2)+
  geom_point(data=year_three_raw_data, aes(x=flo_age_in_wks, y=pardens+0.001, color=factor(mstatus)), shape=1)+
  # geom_point(data=year_three_raw_data, aes(color=factor(mstatus)))+
  # geom_line(data=year_three_raw_data, aes(group=id))+
  geom_label(size=3, data=outcomes_df, color="black", aes(x=52, y=3*0.00001, label=paste("para_12_24=", total_n_para_12_24, sep="")))+
  geom_label(size=3,data=outcomes_df, color="red", aes(x=52, y=3*0.0001, label=paste("mala_12_24=", total_n_malaria_12_24, sep="")))+
  facet_wrap(~id, scales="free_y", nrow=3)+
  xlab("parasite density")+
  scale_color_manual(values=c("darkgrey", "black", "red", "purple"))+
  scale_y_log10()+
  scale_x_continuous(limits = c(48, 108), breaks = c(seq(48, 108, by=8)))+
  theme_minimal()+
  theme(legend.position = "bottom")



long_monthly_categories%>%
  filter(id==11750)%>%
  select(id, flo_age_in_wks, four_week_increment, pardens, qPCRparsdens)%>%
  print(n=45)



# export ####

# static metadata: demographics, genotypes, birth outcomes etc ####
meta_raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/Oct25/MICDSpecimenBoxOct25_withclinical.dta")
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
                                 3~"DP 2 years"))%>%
  mutate(hbs=mic_drop_hbs$HbS[match(as.numeric(id), mic_drop_hbs$id)],
         hbs=case_match(hbs,
                        1~"HbAA",
                        2~"HbAS",
                        3~"HbSS"))

infs_and_meta <- meta%>%
  left_join(., incident_and_prevalent_outcomes, by="id")

write.csv(infs_and_meta, "~/postdoc/stanford/clinical_data/MICDROP/micdrop_metadata_for_immunology.csv", row.names = F)



# sand box ####



year_three_raw_data%>%
  arrange(desc(mstatus))%>%
  ggplot(., aes(x=flo_age_in_wks, y=pardens+0.001))+
  geom_vline(xintercept = c(104, 156), linetype="dashed", color="grey")+
  geom_point(aes(color=factor(mstatus)))+
  geom_line(aes(group=id))+
  geom_label(size=3, data=outcomes_df, color="black", aes(x=10, y=0.0001, label=paste("para_24_36=", total_n_para_24_36, sep="")))+
  geom_label(size=3,data=outcomes_df, color="red", aes(x=10, y=0.001, label=paste("mala_24_36=", total_n_malaria_24_36, sep="")))+
  facet_wrap(~id, scales="free_y", nrow=3)+
  scale_color_manual(values=c("darkgrey", "black", "red"))+
  scale_y_log10()+
  scale_x_continuous(limits = c(-30, NA))+
  theme_minimal()+
  theme(legend.position = "bottom")


