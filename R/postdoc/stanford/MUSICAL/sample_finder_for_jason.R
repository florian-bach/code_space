prism_data_base <- haven::read_dta("~/Downloads/PRISM Border Cohort Study all visits database_FINAL.dta")
prism_pbmc <-  haven::read_dta("~/Downloads/PRISMBC_pbmc_final.dta")

not_individuals <- c(178, 410, 354, 149, 657, 706, 612, 658, 587, 559, 636)

prism_malaria_episodes <- prism_data_base%>%
  filter(malaria%in%c(2,3), ageyrs>3, ageyrs<15)%>%
  distinct(cohortid, date_numeric, ageyrs, malaria, temperature)%>%
  mutate(episode_date=date_numeric)%>%
  select(-date_numeric)%>%
  group_by(cohortid)%>%
  mutate(n_infection=paste("infection_", seq(1,n()), sep=""))%>%
  pivot_wider(names_from = n_infection, values_from = episode_date)

slim_prism_pbmc <- prism_pbmc%>%
  mutate(pbmc_date=date_numeric, cohortid=id)%>%
  select(-date_numeric,-id)%>%
  select(cohortid,pbmc_date, barcode, BoxNumber, boxrow, boxcolumn, Freezer, Shelf)

combo_df <- slim_prism_pbmc%>%
  inner_join(., prism_malaria_episodes, by="cohortid")%>%
  mutate(days_since_malaria1=pbmc_date-infection_1,
         days_since_malaria2=pbmc_date-infection_2,
         days_since_malaria3=pbmc_date-infection_3)
  
select_combo_df <- combo_df%>%
  filter(
    (days_since_malaria1 >3 &days_since_malaria1<10) | (days_since_malaria2>3& days_since_malaria1<10 )| (days_since_malaria3 >3& days_since_malaria1<10 ))%>%
  filter(cohortid %notin% not_individuals)%>%
  filter(grepl("^P1", barcode))%>%
  select(cohortid, barcode, ageyrs, temperature, days_since_malaria1, Freezer, Shelf, BoxNumber, boxrow, boxcolumn)%>%
  arrange(BoxNumber)



write.csv(select_combo_df, "~/Downloads/prism_day7s.csv", row.names = F)


