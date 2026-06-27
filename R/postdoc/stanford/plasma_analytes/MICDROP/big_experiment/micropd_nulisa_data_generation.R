library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

`%notin%` <- Negate(`%in%`)

# read in pilots NULISA data ####
pilot_nulisa <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/extra_24/combined_pilot_data.csv")

# get data in same shape as big experiemnt
pilot_nulisa_edit <- pilot_nulisa%>%
  mutate("sample"=gsub("w", "", sample_id))%>%
  mutate("sample"=gsub("_", "_tp", sample))%>%
  mutate("file_name"="doesnt_matter")%>%
  mutate(timepoint2=stringr::str_split(sample, pattern = "_", simplify = T)[,2])%>%
  mutate(plate=study)%>%
  mutate(conc=concentration)%>%
  mutate(study="MICDROP")%>%
  dplyr::select(targetName, sample, conc, file_name, id, timepoint2, study, plate)
  
# read in big NULISA data ####
excel_files <- list.files(path="~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/Data Files/", pattern = "*Report.xlsx", full.names = T)

big_data <- data.frame()

for(file in excel_files){
  
  sheet <- readxl::read_excel(file, sheet=1)
  
  long_sheet <- sheet %>%
    pivot_longer(cols=colnames(sheet)[-1], names_to = "sample", values_to = "conc")%>%
    mutate("file_name"=file)%>%
    mutate(id=stringr::str_split(sample, pattern = "_", simplify = T)[,1],
           timepoint2=stringr::str_split(sample, pattern = "_", simplify = T)[,2])%>%
    mutate(id=as.numeric(id))%>%
    filter(!is.na(id))%>%
    mutate(study=if_else(timepoint2 %in% c("pregnant", "post", "delivery"), "DPSP", "MICDROP"))%>%
    mutate(plate=unique(substr(file_name, 126, 132)))
  
  big_data <- rbind(big_data, long_sheet)
  
}

very_big_data <- rbind(big_data, pilot_nulisa_edit)

micdrop_nulisa <- very_big_data%>%
  filter(study=="MICDROP")%>%
  filter(sample %notin% c("10720_tp24_1", "10640_tp24_1"))%>%
  mutate(timepoint=case_when(timepoint2=="tp7"~"8 weeks",
                              timepoint2=="tp8"~"8 weeks",
                              timepoint2=="tp23"~"24 weeks",
                              timepoint2=="tp24"~"24 weeks",
                              timepoint2=="tp51"~"52 weeks",
                              timepoint2=="tp52"~"52 weeks",
                              timepoint2=="tp68"~"68 weeks"))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks", "68 weeks")))%>%
  mutate(timepoint_num=as.numeric(substr(timepoint, 1, 2)))


# read in metadata####
## clinical data at sampling visit####
raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/May25/MICDSpecimenBoxMay25_withclinical.dta")

master_plate_map <- readxl::read_excel("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/master_plate_map.xlsx", col_types = c("numeric", "date", "numeric", "text", "guess", "guess", "guess"))
mic_drop_hbs <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/MICDROP SickleTr final.dta")

big_experiment_samples <- master_plate_map%>%
  dplyr::select(id, date)%>%
  filter(!is.na(id))%>%
  mutate(date=as.Date(date))

pilot_samples <- raw_data%>%
  filter(id %in% c(unique(pilot_nulisa_edit$id)))%>%
  filter(Plasma!="")%>%
  dplyr::select(id, date, ageinwks, Plasma)%>%
  mutate(sample=paste(id, "_tp", ageinwks, sep=""))%>%
  filter(sample %in% unique(pilot_nulisa_edit$sample))%>%
  add_row(id=10766, date=as.Date("2022-11-04"), ageinwks=25, Plasma="10766-AB", sample="10766_tp24")

all_samples <- bind_rows(pilot_samples, big_experiment_samples)%>%
  dplyr::select(id, date)%>%
  mutate(hbs=mic_drop_hbs$HbS[match(as.numeric(id), mic_drop_hbs$id)],
         hbs=case_match(hbs,
                        1~"HbAA",
                        2~"HbAS",
                        3~"HbSS"))

# write.csv(all_samples, "~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/all_samples.csv", row.names = F)

metadata_columns <- c("id","date", "ageinwks", "mstatus", "pardens", "qPCRparsdens", "fever", "febrile", "gender", "temp")

#merge based on id and date
visit_metadata <- raw_data%>%
  mutate(date=as.Date(lubridate::parse_date_time(date, "%y/%m/%d")))%>%
  dplyr::select(all_of(metadata_columns))%>%
  right_join(., all_samples, by=c("id", "date"))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(timepoint_num=as.numeric(ageinwks))%>%
  mutate(timepoint_num=case_when(timepoint_num==9~8,
                                 timepoint_num==25~24,
                                 timepoint_num==53~52,
                                 .default=timepoint_num))%>%
  mutate(gender_categorical=if_else(gender==1, "male", "female"))
    
## infection data ####
static_metadata <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/micdrop_metadata_for_immunology.csv")

## merge all metadata, apply qc ####
qc_results <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/NULISA_QC_summary.csv")

micdrop_nulisa_with_meta <- micdrop_nulisa%>%
  inner_join(., static_metadata, by=c("id"))%>%
  left_join(., visit_metadata, by=c("id", "timepoint_num"), multiple = "first")%>%
  mutate(qc_failed=if_else(sample %in% qc_results$sample_id, TRUE, FALSE))


#13% dropout; 642 - 86 samples; mostly plate 8 (37)

write.csv(micdrop_nulisa_with_meta, "~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv", row.names = F)


dpsp_nulisa <- very_big_data%>%
  filter(study=="DPSP")

write.csv(dpsp_nulisa, "~/postdoc/stanford/plasma_analytes/DPSP/dpsp_nulisa_data.csv", row.names = F)

# 
# micdrop_nulisa%>%
#   select(-targetName, -conc)%>%
#   distinct() -> micdrop_nulisa_metadata
# 
# write.csv(micdrop_nulisa_metadata, "~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/NULISA/micdrop_nulisa_metadata.csv")
# 

# 
# clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")
# data_for_petter <- clean_data%>%
#   select(-dob, -date, -gi, -file_name, -timepoint2, -gender, -qPCRdich, -mqPCRparsdens, -total_n_para, -total_n_malaria)
# 
# write.csv(data_for_petter, "~/Downloads/data_for_petter.csv", row.names = F)
# 
# col_names <- data_for_petter%>%
#   filter(id=="h")
# 
# write.csv(col_names, "~/Downloads/col_names.csv", row.names = F)
# 

