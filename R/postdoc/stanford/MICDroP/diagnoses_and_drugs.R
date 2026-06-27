raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/MICDroP Data/MICDROP expanded database through April 30th 2026.dta")

diagnostic_codes <- readxl::read_excel("~/Library/CloudStorage/Box-Box/DP+SP study/Trial clinical documents/Supporting documents/Clinical Support Documents/3_Diagnostic codes_FINAL.xlsx")

slim_diag_data <- raw_data%>%
  select(id, dob, date, starts_with("diag"))%>%
  pivot_longer(cols = starts_with("diag"), names_to = "diag_num", values_to="Code")%>%
  filter(!is.na(Code))%>%
  distinct(id, dob, date, diag_num, Code)%>%
  mutate(flo_age_in_weeks=as.numeric(date-dob)/7)%>%
  left_join(., diagnostic_codes, by="Code")

write.csv(slim_diag_data, "~/postdoc/stanford/clinical_data/MICDROP/all_diagnoses.csv", row.names = F)

# medication_codes_dpsp <- readxl::read_excel("~/Library/CloudStorage/Box-Box/DP+SP study/Trial clinical documents/Supporting documents/Clinical Support Documents/4_Medication codes_FINAL.xlsx")
medication_codes_micdrop <- readxl::read_excel("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study//Trial clinical documents/Supporting documents/Clinical Support Documents/4_Medication codes_FINAL.xlsx")


slim_drug_data <- raw_data%>%
  select(id, date, starts_with(c("med","nMed")))%>%
  pivot_longer(cols = starts_with(c("med","nMed")), names_to = "drug_num", values_to="Code")%>%
  filter(!is.na(Code))%>%
  distinct(id, date, drug_num, Code)%>%
  left_join(., medication_codes_micdrop, by="Code")

write.csv(slim_drug_data, "~/postdoc/stanford/clinical_data/MICDROP/all_drugs.csv", row.names = F)

slim_drug_data%>%
  group_by(Medicaton)%>%
  count()%>%
  arrange(desc(n))%>%
  print(n=25)



slim_diag_data <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/all_diagnoses.csv")
epi_data <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/micdrop_metadata_for_immunology.csv")

msd_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/msd_results_v2.csv")




mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")
maternal_treatment_arms <- haven::read_dta("~/Library/CloudStorage/Box-Box/DP+SP study/Databases and preliminary findings/Final database used for analyses/DPSP treatment allocation_FINAL.dta")
maternal_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/DP+SP study/Databases and preliminary findings/Final database used for analyses/DPSP enrollment analysis database_FINAL.dta")

epi_data <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/micdrop_metadata_for_immunology.csv")
nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")
msd_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/msd_results_v2.csv")
vaccine_coverage <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/long_vaccine_coverage_with_recall.csv")
IU_cutoffs <- readxl::read_excel("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/iu_conversion.xlsx")






diagnoses_of_kids_in_msd <- slim_diag_data%>%
  filter(id %in% msd_data$SubjectID)
write.csv(diagnoses_of_kids_in_msd, "~/postdoc/stanford/clinical_data/MICDROP/diagnoses_of_kids_in_msd.csv", row.names = F)





df <- slim_diag_data%>%
  left_join(epi_data, by="id")%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(id, mic_drop_key$id)])%>%
  group_by(treatmentarm, Code)%>%
  count()%>%
  pivot_wider(names_from = treatmentarm, values_from = n)

colnames(df) <-  c("Code", "placebo", "DP1", "DP2")      

df_clean <- df%>%
  # select(placebo, DP1, DP1)%>%
  filter(!is.na(placebo) & !is.na(DP1) & !is.na(DP2))


mat <- as.matrix(df_clean[,  c("placebo", "DP1", "DP2")])

results <- df_clean %>%
  rowwise() %>%
  mutate(
    test = list(
      chisq.test(c_across(c("placebo", "DP1", "DP2")))
    ),
    p_value = test$p.value
  ) %>%
  ungroup()%>%
  mutate(padj=p.adjust(p_value, method="BH"))%>%
  mutate(disease=slim_diag_data$Diagnosis[match(Code, slim_diag_data$Code)])%>%
  arrange(padj)


slim_diag_data%>%
  left_join(epi_data, by="id")%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(id, mic_drop_key$id)])%>%
  group_by(treatmentarm, Code, total_n_para_36)%>%
  count()%>%
  filter(Code==183)%>%
  filter(!is.na(total_n_para_36))%>%
  ggplot(., aes(x=factor(total_n_para_36), y= n))+
  geom_boxplot()+
  # facet_wrap(~treatmentarm)+
  theme_minimal()


urti_count3 <- slim_diag_data%>%
  filter(flo_age_in_weeks<=156)%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(id, mic_drop_key$id)])%>%
  filter(Code==183)%>%
  group_by(id)%>%
  count(name = "n_urti")
  
urti_count2<- slim_diag_data%>%
  filter(flo_age_in_weeks<=104)%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(id, mic_drop_key$id)])%>%
  filter(Code==183)%>%
  group_by(id)%>%
  count(name = "n_urti")

urti_count1 <- slim_diag_data%>%
  filter(flo_age_in_weeks<=52)%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(id, mic_drop_key$id)])%>%
  filter(Code==183)%>%
  group_by(id)%>%
  count(name = "n_urti")






para_counts <- epi_data%>%
  filter(is.na(withdrawalage))%>%
  select(id, treatmentarm, starts_with("total_n"))%>%
  distinct()
  
para_counts%>%
  left_join(., urti_count3, by="id")%>%
  ggplot(., aes(x=total_n_para_36, y= n_urti, color=treatmentarm))+
  geom_point(position = position_jitter(height=0.3, width=0.3))+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  theme_minimal()

para_counts%>%
  left_join(., urti_count3, by="id")%>%
  ggplot(., aes(x=total_n_malaria_36, y= n_urti, color=treatmentarm))+
  geom_point(position = position_jitter(height=0.3, width=0.3))+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  theme_minimal()

para_counts%>%
  left_join(., urti_count3, by="id")%>%
  ggplot(., aes(x=total_n_para_36-total_n_malaria_36, y= n_urti, color=treatmentarm))+
  geom_point(position = position_jitter(height=0.3, width=0.3))+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  theme_minimal()




 para_counts%>%
  left_join(., urti_count2, by="id")%>%
  ggplot(., aes(x=total_n_malaria_24, y= n_urti, color=treatmentarm))+
  geom_point(position = position_jitter(height=0.3, width=0.3))+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  theme_minimal()

para_counts%>%
  left_join(., urti_count2, by="id")%>%
  ggplot(., aes(x=total_n_para_24, y= n_urti, color=treatmentarm))+
  geom_point(position = position_jitter(height=0.3, width=0.3))+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  theme_minimal()


para_counts%>%
  left_join(., urti_count2, by="id")%>%
  ggplot(., aes(x=total_n_para_24-total_n_malaria_24, y= n_urti, color=treatmentarm))+
  geom_point(position = position_jitter(height=0.3, width=0.3))+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  theme_minimal()





para_counts%>%
  left_join(., urti_count1, by="id")%>%
  ggplot(., aes(x=total_n_malaria_12, y= n_urti, color=treatmentarm))+
  geom_point(position = position_jitter(height=0.3, width=0.3))+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  theme_minimal()

para_counts%>%
  left_join(., urti_count1, by="id")%>%
  ggplot(., aes(x=total_n_para_12, y= n_urti, color=treatmentarm))+
  geom_point(position = position_jitter(height=0.3, width=0.3))+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  theme_minimal()

para_counts%>%
  left_join(., urti_count1, by="id")%>%
  ggplot(., aes(x=total_n_para_12-total_n_malaria_12, y= n_urti, color=treatmentarm))+
  geom_point(position = position_jitter(height=0.3, width=0.3))+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  theme_minimal()







data_filtered <- slim_diag_data %>%
  filter(flo_age_in_weeks <= 156, Code == 183) %>%
  mutate(treatmentarm = mic_drop_key$treatmentarm[match(id, mic_drop_key$id)],
         anyDP=ifelse(treatmentarm==1, "no DP", "DP"))

bootstrap_density <- function(x, grid, R = 500){
  dens_mat <- replicate(R, {
    samp <- sample(x, replace = TRUE)
    density(samp, from=min(grid), to=max(grid), n=length(grid))$y
  })
  
  tibble(
    x = grid,
    y = apply(dens_mat, 1, median),
    ymin = apply(dens_mat, 1, quantile, .025),
    ymax = apply(dens_mat, 1, quantile, .975)
  )
}


grid <- seq(0, 156, length.out = 200)

density_ci <- data_filtered %>%
  group_by(treatmentarm) %>%
  group_modify(~ bootstrap_density(.x$flo_age_in_weeks, grid)) %>%
  ungroup()


ggplot(density_ci,
       aes(x = x, y = y, color = factor(treatmentarm), fill = factor(treatmentarm))) +
  
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.25, color = NA) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = c(52,104), linetype="dashed") +
  theme_minimal() +
  labs(y = "Density", x = "Age (weeks)", color = "Arm", fill = "Arm")


