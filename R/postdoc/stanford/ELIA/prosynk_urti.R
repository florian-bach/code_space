prosynk_ae <- readRDS("~/Library/CloudStorage/Box-Box/Florian Bach's Externally Shareable Files/ELIA/PROSYNK_data_for_Stanford_ELIA/Data/Main/rds/AE_Data.rds")

prosynk_enrolment_data <- readRDS("~/Library/CloudStorage/Box-Box/Florian Bach's Externally Shareable Files/ELIA/PROSYNK_data_for_Stanford_ELIA/Data/Main/rds/Enrolment_Data.rds")

prosynk_dobs <- prosynk_enrolment_data%>%
  mutate(id=IDNumber, dob=enr_DtDel)%>%
  mutate(treatmentarm=case_match(enr_EnlAlloc,
                                 1~"Labinic",
                                 2~"Lab4b",
                                 3~"Probiotic",
                                 4~"Placebo"))%>%
  select(id, dob, treatmentarm)
  
  

# filter out kids with incomplete followup
prosynk_urtis <- prosynk_ae%>%
  mutate(id=IDNumber, URTI=aes_URTI, date=aes_DTAESRT)%>%
  select(id, date, aes_NAEDx, URTI)%>%
  left_join(., prosynk_dobs, by="id")%>%
  mutate(flo_age_in_wks=(date-dob)/7)%>%
  filter(flo_age_in_wks<=104 & flo_age_in_wks>0)%>%
  group_by(id)%>%
  filter(max(flo_age_in_wks)<100)

prosynk_individual_cases <- prosynk_urtis%>%
  group_by(id, treatmentarm)%>%
  summarise(total_n_URTI_24=sum(URTI[flo_age_in_wks<=104]),
            total_n_URTI_12=sum(URTI[flo_age_in_wks<=52]),
            total_n_URTI_6=sum(URTI[flo_age_in_wks<=24]),
            total_n_malaria_24=sum(grepl("MALARIA", aes_NAEDx)[flo_age_in_wks<=104]),
            total_n_malaria_12=sum(grepl("MALARIA", aes_NAEDx)[flo_age_in_wks<=52]),
            total_n_malaria_6=sum(grepl("MALARIA", aes_NAEDx)[flo_age_in_wks<=24]))%>%
  mutate(study="PROSYNK")


prosynk_individual_cases%>%
  ggplot(., aes( x=factor(treatmentarm), y=total_n_malaria_24, fill = treatmentarm))+
  geom_point(position = position_jitter(width=0.2, height=0.1))+
  geom_violin(alpha=0.6)+
  geom_smooth(method="lm")+
  ggpubr::stat_compare_means(method = "wilcox.test", ref.group = "Placebo")+
  theme_minimal()

prosynk_individual_cases%>%
  ggplot(., aes(y=total_n_URTI_24, x=total_n_malaria_24, colour = treatmentarm))+
  geom_point(position = position_jitter(width=0.2, height=0.1))+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  scale_x_continuous(breaks=seq(1,10))+
  theme_minimal()

prosynk_individual_cases%>%
  ggplot(., aes(y=total_n_URTI_6, x=total_n_malaria_6))+
  geom_point(position = position_jitter(width=0.2, height=0.1))+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  scale_x_continuous(breaks=seq(1,10))+
  theme_minimal()


prosynk_micdrop_urti <- para_counts%>%
  left_join(., urti_count2, by="id")%>%
  mutate("study"="MICDROP", total_n_URTI_24=n_urti, id=as.character(id))%>%
  filter(treatmentarm=="Placebo")%>%
  select(id, treatmentarm, total_n_malaria_24, total_n_URTI_24, study)%>%
  bind_rows(prosynk_individual_cases%>%
              select(id, treatmentarm, total_n_malaria_24, total_n_URTI_24, study))
  
  
ggplot(prosynk_micdrop_urti, aes(y=total_n_URTI_24, x=total_n_malaria_24, colour = study))+
  geom_point(position = position_jitter(width=0.2, height=0.1))+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  scale_x_continuous(breaks=seq(1,10))+
  theme_minimal()
