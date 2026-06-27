
# proteomics ####
nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/tr1_paper/revised_baseline_clean_musical_combo_with_metadata.csv")

nulisa_samples <- nulisa_data%>%
  distinct(id, date, infectiontype, timepoint_imm, plasma1_barcode)%>%
  mutate(barcode=plasma1_barcode)%>%
  select(-plasma1_barcode)%>%
  mutate(id=as.numeric(id))%>%
  mutate(date=as.Date(lubridate::parse_date_time(date, "%y/%m/%d")))%>%
  mutate(timepoint_imm=as.numeric(timepoint_imm))%>%
  mutate(data_type="NULISA")
  
# nk cells ####
nk_sample_files <- list.files(path = "~/postdoc/stanford/cytometry/musical_nk/drive-download-20260210T012003Z-1-001/", pattern = "xlsx", full.names = T)

list_of_dfs <- list(vector(length = length(nk_sample_files)))

for(i in 1:length(nk_sample_files)){
  temp_file <- readxl::read_excel(nk_sample_files[i])
  temp_file$file_name <- basename(nk_sample_files[i])
  list_of_dfs[[i]] <- temp_file
}

big_nk_df <- do.call(rbind, list_of_dfs)

slim_big_nk_df <- big_nk_df%>%
  mutate(id=as.character(ID))%>%
  distinct(id, timepoint)%>%
  mutate(infectiontype=substr(timepoint, 1, 1))%>%
  mutate(timepoint_imm=as.numeric(substr(timepoint, 2, 10)))%>%
  filter(!is.na(id))

nk_metadata <- readxl::read_excel("~/postdoc/stanford/cytometry/musical_nk/samples_to_analyze_KPedits_Final_2ndaliquot.xlsx")

slim_nk_metadata <- nk_metadata%>%
  select(id, date,  infectiontype, timepoint_imm, pbmc1_barcode)%>%
  filter(!is.na(id))%>%
  mutate(id=as.character(id))

nk_samples <- slim_big_nk_df%>%
  left_join(., slim_nk_metadata, by=c("id", "timepoint_imm", "infectiontype"))%>%
  select(id, date, infectiontype, timepoint_imm, pbmc1_barcode)%>%
  mutate(barcode=pbmc1_barcode)%>%
  select(-pbmc1_barcode)%>%
  mutate(id=as.numeric(id))%>%
  mutate(date=as.Date(date))%>%
  mutate(timepoint_imm=as.numeric(timepoint_imm))%>%
  mutate(data_type="ADCC")

# t cell stim ####

slim_cell_count_data <- read.csv("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/final_data_tables/t_cell_stim.csv")

pre_t_cell_samples <- slim_cell_count_data%>%
  mutate(id=cohortid, timepoint_imm=timepoint)%>%
  mutate(infectiontype=substr(infectiontype, 1, 1))%>%
  distinct(id, date, infectiontype, timepoint_imm, pbmc1_barcode)%>%
  mutate(date=as.Date(lubridate::parse_date_time(date, "%m/%d/%y")))%>%
  mutate(barcode=pbmc1_barcode)%>%
  mutate("data_type"="iRBC T cell stim")
  
pre_t_cell_samples2 <- pre_t_cell_samples%>%
  filter(!is.na(infectiontype), !is.na(date))
  
t_cell_sample_nas <- pre_t_cell_samples%>%
  filter(is.na(infectiontype) | is.na(date))%>%
  left_join(., pras_metadata, by="pbmc1_barcode")%>%
  select(id.y, date.y, infectiontype.y, timepoint_imm.y, barcode)%>%
  filter(!is.na(id.y))%>%
  mutate("id"=id.y,
         "date"=date.y,
         "infectiontype"=infectiontype.y,
         "timepoint_imm"=timepoint_imm.y,
         "data_type"="iRBC T cell stim")%>%
  select(id, date, infectiontype, timepoint_imm, barcode, data_type)
  

t_cell_samples <- bind_rows(t_cell_sample_nas, pre_t_cell_samples2)%>%
  select(-pbmc1_barcode)%>%
  mutate(date=as.Date(date))


# rnaseq ####
caro_data <- readRDS("~/Downloads/NULISA_genes_expression.rds")

rnaseq_samples <- caro_data$pheno%>%
  mutate(id=cohortid, infectiontype=toupper(substr(infection_clean, 1, 1)), timepoint_imm=timepoint, date=as.Date(date_clean))%>%
  distinct(id, date, infectiontype, timepoint_imm, barcode)%>%
  mutate(data_type="whole blood RNAseq")

# spectral flow ####
spectral_flow_metadata <- read.csv("~/postdoc/stanford/cytometry/spectral/MUSICAL/big_experiment/edited_metadata.csv")

pre_spectral_samples <- spectral_flow_metadata%>%
  filter(control==FALSE)%>%
  mutate(timepoint_imm=substr(timepoint, 2, 10))%>%
  mutate(timepoint_imm=case_when(timepoint_imm=="m-1"~"-1",
                                 timepoint_imm=="m0"~"0",
                                 timepoint_imm=="m14"~"14",
                                 timepoint_imm=="m7"~"7",
                                 timepoint_imm=="s0"~"PS",
                                 .default=timepoint_imm))%>%
  filter(id!="<NA>", timepoint_imm!="PS")%>%
  mutate(timepoint_imm=as.numeric(timepoint_imm),
         id=as.numeric(id))%>%
  select(id, infectiontype, timepoint_imm)

pras_metadata <- readxl::read_excel("~/postdoc/stanford/plasma_analytes/MUSICAL/big_data/Immunology List _MUSICAL.xlsx")

pre_spectral_samples2 <- pre_spectral_samples%>%
  left_join(., pras_metadata, by=c("id", "timepoint_imm", "infectiontype"))%>%
  distinct(id, date, infectiontype, timepoint_imm, pbmc1_barcode)

pre_spectral_samples2_nas <- pre_spectral_samples2%>%
  filter(is.na(pbmc1_barcode))

pre_spectral_samples3 <- pre_spectral_samples2%>%
  filter(!is.na(pbmc1_barcode))

pre_spectral_samples4 <- pre_spectral_samples2_nas%>%
  left_join(., t_cell_samples, by=c("id", "timepoint_imm", "infectiontype"))%>%
  mutate(pbmc1_barcode=pbmc1_barcode.y)%>%
  select(pbmc1_barcode)%>%
  left_join(., pras_metadata, by="pbmc1_barcode")%>%
  distinct(id, date, infectiontype, timepoint_imm, pbmc1_barcode)

pre_spectral_samples2 <- rbind(pre_spectral_samples3, pre_spectral_samples4)

spectral_samples <- pre_spectral_samples2%>%
  mutate(data_type="spectral flow")%>%
  mutate(date=as.Date(date))%>%
  mutate(barcode=pbmc1_barcode)%>%
  select(-pbmc1_barcode)

# cytof ####

cytof_metadata <- read.csv("/Users/fbach/Library/CloudStorage/Box-Box/Florian Bach's Externally Shareable Files/stanford/cytometry/CyTOF/MUSICAL/pilot75/single_cell_metadata.csv")

cytof_samples <- cytof_metadata%>%
  filter(class!="control", class!="immune")%>%
  mutate(id=subject_id)%>%
  mutate(infectiontype=toupper(substr(class, 1, 1)))%>%
  mutate(crap_timepoint_imm=substr(timepoint, nchar(timepoint)-2, nchar(timepoint)))%>%
  mutate(timepoint_imm=as.numeric(case_when(crap_timepoint_imm=="ine"~"-1",
                                            crap_timepoint_imm=="y 0"~"0",
                                            crap_timepoint_imm=="y 7"~"7",
                                            crap_timepoint_imm==" 14"~"14")))%>%
  left_join(., pras_metadata, by=c("id", "timepoint_imm", "infectiontype"))%>%
  select(id, date, infectiontype, timepoint_imm, pbmc1_barcode)%>%
  mutate(barcode=pbmc1_barcode)%>%
  select(-pbmc1_barcode)%>%
  mutate(data_type="CyTOF")

  
# aggregate ####

aggregate_sample_df <- rbind(nulisa_samples, nk_samples, t_cell_samples, rnaseq_samples, spectral_samples, cytof_samples)

wide_aggregate_sample_df <- aggregate_sample_df%>%
  filter(!is.na(id))%>%
  pivot_wider(names_from = data_type, values_from = barcode)

write.csv(wide_aggregate_sample_df, "~/postdoc/stanford/clinical_data/MUSICAL/all_samples.csv", row.names = F)
