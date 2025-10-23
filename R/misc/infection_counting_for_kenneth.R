# read in clinical data
raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/May25/MICDSpecimenBoxMay25_withclinical.dta")

# define metadata columns that we want to keep
metadata_columns <- c("id", "dob", "date", "ageinwks", "gender", "mstatus", "qPCRparsdens", "visittype", "fever", "febrile", "rogerson", "anyHP", "GAcomputed", "gi", "SGA", "qPCRdich", "mqPCRparsdens")

# count infection prevalences ####
infs_and_meta <- raw_data%>%
  # filter(ageinwks<54)%>%
  group_by(id) %>%
  # in R, you can sum logic vectors (vectors that only contain true or false); TRUE =1, FALSE=8;
  # so by summing the logic vectors associated our filtering conditions, we are essentially counting the isntances where all conditions are TRUE
  # you should adapt these to fit with the intervals that we have discussed.
  mutate("total_n_para_12"=sum(pardens>1&ageinwks<53, na.rm=TRUE),
         "total_n_malaria_12"=sum(mstatus!=0&ageinwks<53, na.rm=TRUE),
         "total_n_para_24"=sum(pardens>1&ageinwks<105, na.rm=TRUE),
         "total_n_malaria_24"=sum(mstatus!=0&ageinwks<105, na.rm=TRUE),
         "total_n_para_12_24"=sum(pardens>1&ageinwks<105&ageinwks>52, na.rm=TRUE),
         "total_n_malaria_12_24"=sum(mstatus!=0&ageinwks<105&ageinwks>52, na.rm=TRUE),
         "symp_prob_12_24"=total_n_malaria_12_24/total_n_para_12_24,
         "total_n_para_6"=sum(pardens>1&ageinwks<27, na.rm=TRUE),
         "total_n_malaria_6"=sum(mstatus!=0&ageinwks<27, na.rm=TRUE),
         "any_malar_6"=if_else(total_n_malaria_6==0, FALSE, TRUE),
         "any_malar_12"=if_else(total_n_malaria_12==0, FALSE, TRUE),
         "any_malar_12_24"=if_else(total_n_malaria_12_24==0, FALSE, TRUE))%>%
  select(all_of(metadata_columns), any_malar_6, any_malar_12, any_malar_12_24, total_n_para_6, total_n_malaria_6, total_n_para_12, total_n_malaria_12, total_n_para_24, total_n_malaria_24, total_n_para_12_24,total_n_malaria_12_24)%>%
  right_join(., metadata, by=c("id", "date"))


# categorize each month based on infection outcome ####
# plan: categorize every 28 week period as 1 of 4 outcomes: nothing, asymptomatic, symptomatic, complicated malaria;
# if multiple events happen in the same 4 week window, it gets assigne the more severe outcome

monthly_categories <- raw_data %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7)%>%
  # remove NAs and treatment failures
  filter(!is.na(mstatus), !is.na(flo_age_in_wks), mstatus!=3, mstatus!=4)%>%
  # assign four week increment by rounding down age in weeks divided by 4
  mutate("four_week_increment"=floor(flo_age_in_wks/4))%>%
  arrange(id, four_week_increment)%>%
  group_by(id, four_week_increment)%>%
  #define hierarchy of outcomes;
  mutate(outcome=ifelse(any(mstatus==2), "complicated",
                        ifelse(any(mstatus==1), "uncomplicated",
                               ifelse(all(mstatus==0)&any(pardens>0), "asymptomatic",
                                      ifelse(all(mstatus==0)&all(pardens==0 | is.na(pardens)), "nothing", "something_else")))))%>%
  mutate(has_comp=any(mstatus==2))%>%
  mutate(has_uncomp=any(mstatus==1))%>%
  mutate(has_hyper = ifelse(any(pardens>100000), "hyperparasitemia", "normal"))%>%
  # in cases where there are multiple events, keep only most severe outcome
  filter((n() == 1) |
           (has_comp & mstatus == 2) |
           (has_comp==FALSE & has_uncomp & mstatus == 1) |
           (has_hyper=="hyperparasitemia" & pardens==max(pardens))
  )%>%
  group_by(id, outcome)%>%
  mutate(order_of_outcome = seq(1, n()))



## count instances of each outcome for each individual ####
# this makes the data wider, giving a column for each outcome; it also allows us retain information about all
# outcomes e.g. we will know easily how many asymptomatic infections a person had before their 5th symptomatic infection etc.
# that's helpful for modelling e.g. the probabiltiy of complciated malaria given an 5 previous infections etc.
long_monthly_categories <- monthly_categories%>%
  tidyr::pivot_wider(names_from = outcome, values_fill = 0,values_from = order_of_outcome,
                     id_cols = c("id", "date", "flo_age_in_wks", "flo_age_in_months", "age_quarter", "age_semi", "four_week_increment", "mstatus", "pardens", "has_hyper", "study", "treatmentarm", "gender", "mom_rx"))%>%
  
  ungroup()%>%
  arrange(id, date)%>%
  group_by(id)%>%
  #these n are not the order at the time, but the number of preceding episodes
  mutate(n_nothing=cummax(nothing))%>%
  mutate(n_complicated=cummax(complicated))%>%
  mutate(n_uncomplicated=cummax(uncomplicated))%>%
  mutate(n_asymptomatic=cummax(asymptomatic))%>%
  mutate(n_any_para=n_uncomplicated+n_asymptomatic+n_complicated)%>%
  mutate(n_any_malaria=n_uncomplicated+n_complicated)%>%
  mutate(age_at_first_malaria=min(flo_age_in_months[n_any_malaria==1]))%>%
  mutate(age_at_first_para=min(flo_age_in_months[n_any_para==1]))%>%
  # mutate(age_at_first_malariaf=ifelse(any(flo_age_in_months==age_at_first_malaria), paste(age_semi), "none"))%>%
  mutate(anyDP=if_else(treatmentarm=="No DP", "no", "yes"),
         bino_complicated=ifelse(mstatus==2, 1, 0),
         bino_symp=ifelse(mstatus%in%c(1,2), 1, 0))%>%
  mutate(increment_since_cessation=case_when(treatmentarm=="No DP"~four_week_increment,
                                             treatmentarm=="DP 1 year"~four_week_increment-13,
                                             treatmentarm=="DP 2 years"~four_week_increment-26))

