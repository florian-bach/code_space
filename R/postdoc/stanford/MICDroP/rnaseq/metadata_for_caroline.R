raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/May25/MICDSpecimenBoxMay25_withclinical.dta")

rnaseq_samples <- readxl::read_xlsx("/Users/fbach/postdoc/stanford/rna_seq/micdrop/MIC-DroP RNA-seq Aug 2025.xlsx")

static_metadata <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/micdrop_metadata_for_immunology.csv")

# write.csv(all_samples, "~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/all_samples.csv", row.names = F)

metadata_columns <- c("id","date", "ageinwks", "mstatus", "pardens", "qPCRparsdens", "fever", "febrile", "temp")

#merge based on id and date
visit_and_static_metadata <- raw_data%>%
  mutate(date=as.Date(lubridate::parse_date_time(date, "%y/%m/%d")))%>%
  select(all_of(metadata_columns))%>%
  right_join(., rnaseq_samples, by=c("id", "date"))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(timepoint_num=as.numeric(ageinwks))%>%
  mutate(timepoint_num=case_when(timepoint_num==9~8,
                                 timepoint_num==25~24,
                                 timepoint_num==53~52,
                                 .default=timepoint_num))%>%
  left_join(., static_metadata, by=c("id"))%>%
  mutate(gender_categorical=if_else(gender==1, "male", "female"))

write.csv(visit_and_static_metadata, "/Users/fbach/postdoc/stanford/rna_seq/micdrop/proper_metadata_feb26.csv", row.names = F)
