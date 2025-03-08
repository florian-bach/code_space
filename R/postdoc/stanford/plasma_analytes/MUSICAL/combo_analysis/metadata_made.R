prism_metadata <- haven::read_dta("/Users/fbach/Library/CloudStorage/Box-Box/Jagannathan_Lab_Folder/PRISM_Specimens/BorderCohort/PBMC_Samples/PRISMBC_samples_final/PRISMBC_pbmc_plasma_paxgene_final.dta")

three_six_three_specimen <- prism_metadata %>%
  filter(id==363, date=="2022-03-02 00:00:00")

big_musical_metadata <- haven::read_dta("/Users/fbach/Library/CloudStorage/Box-Box/Jagannathan_Lab_Folder/PRISM_Specimens/BorderCohort/PBMC_Samples/PRISMBC_samples_final/PRISM Border Cohort Study all visits database_FINAL.dta")

pras_metadata <- readxl::read_excel("~/postdoc/stanford/plasma_analytes/MUSICAL/big_data/Immunology List _MUSICAL.xlsx")

three_six_three <- big_musical_metadata%>%
  filter(cohortid==363, date=="02mar2022")%>%
  # mutate(masterbarcode=barcode, fever=as.character(fever), qpcr=as.character(qpcr))%>%
  mutate(id=cohortid, 
         date=as.Date("2022-03-02"),
         infectiontype="S",
         timepoint_imm=14,
         masterbarcode=barcode,
         pbmc1_barcode = three_six_three_specimen$barcode[three_six_three_specimen$sampletype=="PBMC" & three_six_three_specimen$aliquot==1],
         plasma1_barcode = three_six_three_specimen$barcode[three_six_three_specimen$sampletype=="Plasma" & three_six_three_specimen$aliquot==1],
         pbmc_date=date,
         ageyrs=ageyrs,
         gender_categorical=as.character(gender),
         fever=as.character(fever),
         qpcr=as.character(qpcr),
         parasitedensity=as.character(parasitedensity),
         temperature=as.character(temperature),
         presymptomatic_sample=NA,
         baseline_date=pras_metadata$date[pras_metadata$id==363 & pras_metadata$infectiontype=="S" & pras_metadata$timepoint_imm==-1],
         cleaner_baseline_date=NA)%>%
  select(colnames(pras_metadata))
# fix gender thing


pras_metadata2 <- bind_rows(pras_metadata, three_six_three)


#update parsitemia
prism_metadata <- haven::read_dta("/Users/fbach/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/Cohort data/PRISM Border Cohort Study all visits database_FINAL.dta")
new_parasitemias <- prism_metadata%>%
  mutate(id=cohortid, new_qpcr=qpcr, new_parasitedensity=parasitedensity)%>%
  select(id, date, new_qpcr, new_parasitedensity)

pras_metadata3 <- left_join(pras_metadata2, new_parasitemias, by=c("id", "date"))

write.csv(pras_metadata3, "~/postdoc/stanford/plasma_analytes/MUSICAL/big_data/musical_metadata_new_qpcr.csv", row.names = F)
write.csv(pras_metadata3, "~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/musical_metadata_new_qpcr.csv", row.names = F)

# work in progress ####
# filter to only include nulisa samples

nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")

# pras_metadata2 <- pras_metadata2 %>%
#   filter(id %in%)

# write.csv(pras_metadata2, "~/postdoc/stanford/plasma_analytes/MUSICAL/big_data/fixed_musical_metadata_final.csv", row.names = F)

# sandbox
