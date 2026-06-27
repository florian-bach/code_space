micdrop_files <- list.files(path="~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/Data Files/", pattern = "*Report.xlsx", full.names = T)

micdrop_nulisa_controls <- data.frame()

for(file in micdrop_files){
  
  sheet <- readxl::read_excel(file, sheet=1)
  
  control_well_names <- tail(colnames(sheet), n = 10)
  print(control_well_names)
  actual_sample_names <- colnames(sheet)[-1][!colnames(sheet)[-1]%in%control_well_names]

  long_sheet <- sheet %>%
    select(targetName, control_well_names)%>%
    mutate("file_name"=file)%>%
    mutate(plate=unique(substr(file_name, 126, 132)))%>%
    mutate(study="micdrop")

  micdrop_nulisa_controls <- rbind(micdrop_nulisa_controls, long_sheet)
  
}



prosynk_data <- readxl::read_excel("/Users/fbach/Library/CloudStorage/Box-Box/Florian Bach's Externally Shareable Files/ELIA/PROSYNK_data_for_Stanford_ELIA/Data/main/xlsx/NULISAseq_Data.xlsx")

prosynk_nulisa_controls <- prosynk_data[,c(4, 6:18)]%>%
  mutate(study="prosynk")
  

combined_controls <- bind_rows(micdrop_nulisa_controls, prosynk_nulisa_controls)

long_combined_controls <- combined_controls%>%
  pivot_longer(control_well_names, names_to = "control_well", values_to = "conc")

long_cvs <- long_combined_controls%>%
  group_by(targetName, control_well, study)%>%
  summarise("mean"=mean(conc), "sd"=sd(conc), "cv"=mean/sd)

wide_cvs <- long_cvs%>%
  pivot_wider(names_from = "study", values_from = c(mean, sd, cv))%>%
  mutate(mean_ratio=mean_micdrop/mean_prosynk,
         sd_ratio=sd_micdrop/sd_prosynk,
         cv_ratio=cv_micdrop/cv_prosynk)

wide_cvs%>%
  filter(control_well %in% c("IPC_1", "IPC_2", "IPC_3"))%>%
  arrange(desc(mean_ratio))
  
  
