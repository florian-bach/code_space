vax <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/vaccines/infant_vaccination_data.csv")

long_vax <- vax %>%
  pivot_longer(cols=c("Bcg", "pcv1", "pCV2", "pCV3"), names_to = "vaccine", values_to = "received")

complete <- long_vax %>%
  filter(received==1)%>%
  select(-received, -sched)%>%
  pivot_wider(names_from = "vaccine", values_from = "date", values_fill = NA, id_cols = id)
  # summarise("full_pcv"=ifelse(any(pcv1==1,
  #                                 pCV2[vaccine=="pCV2"]==1,
  #                                 pCV3[vaccine=="pCV3"]==1), TRUE, FALSE))%>%
  ungroup()
         
  
 duplicates <- long_vax |>
    filter(received==1)|>
    summarise(n = dplyr::n(), .by = c(id, vaccine)) |>
    filter(n > 1L) 
  
 odd_ones <- long_vax %>%
   filter(id %in% duplicates$id, received==1)%>%
   arrange(id, vaccine)
 
 write.csv(odd_ones, "~/postdoc/stanford/clinical_data/MICDROP/vaccines/duplicate_vaccinations.csv")
 