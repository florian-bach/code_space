library(tidyr)
library(dplyr)

# get UCSF shipped samples into R ####
files <- list.files("~/postdoc/stanford/clinical_data/MUSICAL/UCSF_Plasma_Locations", pattern = "csv", full.names = TRUE)

box_list <- vector("list", length = length(files))

for(i in 1:length(files)){
  boxv <- read.csv(files[i], header=FALSE)
  boxv$box <- paste("Box", i) 
  boxv$row <- 1:nrow(boxv)
  box_list[[i]] <- boxv[1:9, 1:11]
}

box_df <- do.call(rbind, box_list)
box_df <- na.omit(box_df)

ucsf_boxes_df <- box_df %>%
  pivot_longer(starts_with("V"), names_to = "column", values_to = "barcode")%>%
  mutate(column=gsub("V", "", .$column), row=as.character(row))

# load in other databases ####
# arefin_data <- read.csv("/Users/fbach/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Arefin for Isabel and Saki/Analysis using data until 2023April23/Primary samples/primary_samples_from_paired_AS/primary_paired_AS_flo_edit.csv")
arefin_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/primary_paired_AS_flo_edit2.csv")
arefin_data$inf_type_dich <- ifelse(grepl("^S", arefin_data$inf_type), "S", ifelse(grepl("^A", arefin_data$inf_type), "A", "?"))

big_metadata <- read.csv("~/postdoc/stanford/clinical_data/MUSICAL/PRISMBC_lab_samples.csv")
plasma_metadata <- filter(big_metadata, sampletype=="Plasma - 1")
# jason_metadata <- read.csv("~/postdoc/stanford/clinical_data/MUSICAL/plasma_metadata.csv")

samples_at_stanford <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/musical_samples_already_at_stanford.csv")

samples_at_stanford <- samples_at_stanford %>%
  mutate(barcode=plasma_metadata$barcode[match(.$masterbarcode, plasma_metadata$Masterbarcode)])%>%
  select(combined_id, combined_date, masterbarcode, barcode)

# Kylie DM'd this file on slack on April 29th; it links plasma barcodes with Olink plate numbers and locations
olink_locations <- read.csv("~/postdoc/stanford/clinical_data/MUSICAL/olink_locations.csv")

samples_in_plates <- inner_join(olink_locations, samples_at_stanford, by="barcode")
samples_in_plates <- samples_in_plates %>%
  mutate("box"=Plate.Name, "row"=as.character(Row), "column"=as.character(Column))%>%
  select(box, row, column, barcode)


# put it together ####
# these are all the samples at stanford
extant_tubes <- ucsf_boxes_df%>%
  bind_rows(samples_in_plates)%>%
  inner_join(., plasma_metadata, by="barcode")%>%
  mutate("masterbarcode"=Masterbarcode,
         "box_row_column"=paste(box, row, column, sep="_"))

# all 405 samples at stanford are contained in arefin's sheet
overlap <- inner_join(extant_tubes, arefin_data, by="masterbarcode") %>%
  mutate("timepoint"=paste(inf_type_dich, day_annotation, sep="_"),
         "timepoint2"=paste(inf_type, day_annotation, sep="_"),
         "box_row_column"=paste(box, row, column, sep="_"),
         "sample_id"=paste(cohortid, timepoint, sep="_"))
  
# overlap <- left_join(overlap, arefin_data, by="masterbarcode")


wide_locations <- overlap %>%
  filter(inf_type_dich %in% c("S", "A"), day_annotation %in% c(-1, 0, 7, 14))%>%
  group_by(timepoint, cohortid) %>%
  mutate(instance = row_number())%>%
  # select(cohortid, inf_type, timepoint, box_row_column)%>%
  pivot_wider(names_from = timepoint, values_from = c(box_row_column), id_cols = c(cohortid, instance))%>%
  select(-instance)%>%
  select(cohortid, `S_-1`, S_0, S_7, S_14,`A_-1`,A_0,A_14)%>%
  arrange(cohortid)

wide_barcodes <- overlap %>%
  filter(day_annotation %in% c(-1, 0, 7, 14))%>%
  group_by(timepoint, cohortid) %>%
  mutate(instance = row_number())%>%
  # select(cohortid, inf_type, timepoint, box_row_column)%>%
  pivot_wider(names_from = timepoint, values_from = c(barcode), id_cols = c(cohortid, instance))%>%
  select(-instance)%>%
  select(cohortid, `S_-1`, S_0, S_7, S_14,`A_-1`,A_0,A_14)%>%
  arrange(cohortid)

write.csv(wide_locations, "~/postdoc/stanford/clinical_data/MUSICAL/wide_locations.csv", row.names = FALSE)
write.csv(wide_barcodes, "~/postdoc/stanford/clinical_data/MUSICAL/wide_barcodes.csv", row.names = FALSE)

# 
# duplicates <- overlap %>%
#   group_by(sample_id)%>%
#   add_count(name="n")%>%
#   filter(n>1)%>%
#   select(cohortid, timepoint, combined_date, barcode)
# 
# write.csv(duplicates, "~/postdoc/stanford/clinical_data/MUSICAL/duplciate_plasmas_in_arefin.csv")
# 
# # duplicates <- overlap|>
# #   dplyr::summarise(n = dplyr::n(), .by = c(cohortid, inf_type, timepoint)) |>
# #   dplyr::filter(n > 1L) 
# 
# duplicate_context <- overlap%>%
#   filter(cohortid %in% duplicates$cohortid & timepoint %in% duplicates$timepoint)%>%
#   select(cohortid, inf_type, timepoint, combined_date, barcode)%>%
#   arrange(cohortid, timepoint)
# 
# write.csv(duplicate_context, "~/postdoc/stanford/clinical_data/MUSICAL/duplciate_plasmas_context_in_arefin.csv")
# 


# non malarial fevers ####

slim_metadata <- big_metadata %>%
  filter(sampletype=="Plasma - 1")
select(masterbarcode, boxnumber, boxrow, boxcolumn)

nmf <- read.csv("/Users/fbach/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Arefin for Isabel and Saki/Analysis using data until 2023April23/Primary samples/nmf_samples_from_paired_AS/nmf_samples_for_12_individuals_outof_51AS_paired.csv")

nmf_plasmas <- nmf %>%
  left_join(., slim_metadata, by="masterbarcode")%>%
  filter(day_annotation %in% c(-1, 0, 7))%>%
  select(combined_id, combined_date, masterbarcode, day_annotation, barcode, boxnumber, boxrow, boxcolumn)

write.csv(nmf_plasmas, "~/postdoc/stanford/plasma_analytes/MUSICAL/nmf_plasmas.csv")




# confusion ####

big_metadata <- read.csv("~/postdoc/stanford/clinical_data/MUSICAL/PRISMBC_lab_samples.csv")

slim_metadata <- big_metadata %>%
  filter(sampletype=="Plasma - 1")%>%
  mutate(masterbarcode=Masterbarcode)
  
consistent_meta <- inner_join(arefin_data, slim_metadata, by = "masterbarcode")

realism <- consistent_meta %>%
  # mutate("timepoint"=paste(inf_type_dich, day_annotation, sep="_"))%>%
  # group_by(cohortid, inf_type_dich, day_annotation)%>%
  # mutate(instance = row_number())%>%
  # select(cohortid, inf_type, timepoint, box_row_column)%>%
  pivot_wider(names_from = day_annotation, values_from = c(barcode), id_cols = c(cohortid))%>%
  # select(-instance)%>%
  select(cohortid, `-1`, `0`, `7`, `14`)%>%
  arrange(cohortid)


here <- overlap %>%
  # mutate("timepoint"=paste(inf_type_dich, day_annotation, sep="_"))%>%
  # group_by(cohortid, inf_type_dich, day_annotation)%>%
  # mutate(instance = row_number())%>%
  # select(cohortid, inf_type, timepoint, box_row_column)%>%
  pivot_wider(names_from = day_annotation, values_from = c(barcode), id_cols = c(cohortid))%>%
  # select(-instance)%>%
  select(cohortid, `-1`, `0`, `7`, `14`)%>%
  arrange(cohortid)


# try again without symptom info

wide_locations2 <- overlap %>%
  filter(day_annotation %in% c(-1, 0, 7, 14))%>%
  group_by(day_annotation, cohortid) %>%
  mutate(instance = row_number())%>%
  # select(cohortid, inf_type, timepoint, box_row_column)%>%
  pivot_wider(names_from = day_annotation, values_from = c(box_row_column), id_cols = c(cohortid, instance))%>%
  select(-instance)%>%
  select(cohortid, `-1`, `0`, `7`, `14`)%>%
  arrange(cohortid)

wide_barcodes2 <- overlap %>%
  filter(day_annotation %in% c(-1, 0, 7, 14))%>%
  group_by(day_annotation, cohortid) %>%
  mutate(instance = row_number())%>%
  # select(cohortid, inf_type, timepoint, box_row_column)%>%
  pivot_wider(names_from = day_annotation, values_from = c(barcode), id_cols = c(cohortid, instance))%>%
  select(-instance)%>%
  select(cohortid, `-1`, `0`, `7`, `14`)%>%
  arrange(cohortid)

write.csv(wide_locations2, "~/postdoc/stanford/clinical_data/MUSICAL/wide_locations2.csv", row.names = FALSE)
write.csv(wide_barcodes2, "~/postdoc/stanford/clinical_data/MUSICAL/wide_barcodes2.csv", row.names = FALSE)



#pras coding of timepoints and infections 


pras_data <- read.csv("~/Downloads/Immunology List _MUSICAL.csv")
pilot_metadata <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/MASTER_METADATA.csv")

pras_data_with_location <- left_join(pras_data, extant_tubes, by = "masterbarcode")

wide_pras_barcodes <- pras_data_with_location %>%
  mutate("inf_type_dich"=ifelse(grepl("^S", infectiontype), "S", ifelse(grepl("^A", infectiontype), "A", "?")),
         "timepoint"=paste(inf_type_dich, timepoint_imm, sep=" "))%>%
  # group_by(cohortid, inf_type_dich, day_annotation)%>%
  # mutate(instance = row_number())%>%
  # select(cohortid, inf_type, timepoint, box_row_column)%>%
  pivot_wider(names_from = timepoint, values_from = c(plasma1_barcode), id_cols = c(id))%>%
  # select(-instance)%>%
  # select(cohortid, `-1`, `0`, `7`, `14`)%>%
  arrange(id)

wide_pras_locations <- pras_data_with_location %>%
  mutate("inf_type_dich"=ifelse(grepl("^S", infectiontype), "S", ifelse(grepl("^A", infectiontype), "A", "?")),
         "timepoint"=paste(inf_type_dich, timepoint_imm, sep=" "))%>%
  # group_by(cohortid, inf_type_dich, day_annotation)%>%
  # mutate(instance = row_number())%>%
  # select(cohortid, inf_type, timepoint, box_row_column)%>%
  pivot_wider(names_from = timepoint, values_from = c(box_row_column), id_cols = c(id))%>%
  # select(-instance)%>%
  # select(cohortid, `-1`, `0`, `7`, `14`)%>%
  arrange(id)


`%notin`=Negate(`%in%`)

hi <- pras_data_with_location%>%
  filter("id" %notin% unique(pilot_metadata$icombined_id))%>%
  group_by(id)%>%
  summarise(n=n())


wide_pras_locations <- pras_data_with_location %>%
  mutate("inf_type_dich"=ifelse(grepl("^S", infectiontype), "S",
                                ifelse(grepl("^A", infectiontype), "A", infectiontype)),
         "timepoint"=paste(inf_type_dich, timepoint_imm, sep=" "),
         "barcode_location"=paste(barcode, box_row_column))%>%
  group_by(cohortid, inf_type_dich, timepoint)%>%
  mutate(instance = row_number())%>%
  pivot_wider(names_from = timepoint, values_from = barcode_location, id_cols = c(id, instance))%>%
  dplyr::select(-instance)%>%
  # select(cohortid, `-1`, `0`, `7`, `14`)%>%
  arrange(id)

write.csv(wide_pras_locations, "~/postdoc/stanford/clinical_data/MUSICAL/pras_barcode_locations.csv")



#sanity check

View(wide_pras_locations %>%
       pivot_longer(cols=colnames(wide_pras_locations)[2:ncol(wide_pras_locations)], names_to = "timepoint", values_to = "barcode_location")%>%
       filter(!is.na(barcode_location)))



long_pras_locations <- pras_data_with_location %>%
  mutate(#"inf_type_dich"=ifelse(grepl("^S", infectiontype), "S",
    # ifelse(grepl("^A", infectiontype), "A", infectiontype)),
    #"timepoint"=paste(inf_type_dich, timepoint_imm, sep=" "),
    "barcode_location"=paste(barcode, box_row_column),
    "in_pilot"=ifelse(id %in% musical_metadata$combined_id, TRUE, FALSE),
    timepoint_imm=factor(timepoint_imm, levels=c(-1, 0, 7, 14, -2)),
    infectiontype=factor(infectiontype, levels=c("S", "A", "S2", "A2", "PS", "CHECK SAMPLE")),
    cohortid=id)%>%
  filter(infectiontype!="NM", in_pilot==FALSE)%>%
  group_by(cohortid)%>%
  arrange(infectiontype, timepoint_imm)%>%
  mutate(new_column = row_number())%>%
  dplyr::select(cohortid, infectiontype, timepoint_imm, barcode_location, new_column)%>%
  # select(cohortid, `-1`, `0`, `7`, `14`)%>%
  arrange(cohortid, infectiontype, timepoint_imm)

write.csv(long_pras_locations,  "~/postdoc/stanford/clinical_data/MUSICAL/long_pras_barcode_locations_fixed.csv")


#check whether each individual has each timepoint

wanted_timepoints <- c("S-1", "S0", "S7", "S14", "A-1", "A0", "A14")
searching <- pras_data %>%
  filter(infectiontype!="NM")%>%
  mutate("inf_type_dich"=ifelse(grepl("^S", infectiontype), "S",
                                ifelse(grepl("^A", infectiontype), "A", infectiontype)),
         flo_timepoint=paste(inf_type_dich, timepoint_imm, sep=""))%>%
  group_by(id)%>%
  mutate(is_complete=ifelse(all(wanted_timepoints %in% flo_timepoint), TRUE, FALSE))


search_list <- split(searching, searching$id)
search_list <- sapply(search_list, function(x) wanted_timepoints[wanted_timepoints%notin% x$flo_timepoint])
cleaned_search_list <- search_list[lapply(search_list, length)>0]

ucsf_void <- read.table("~/Downloads/samples_missing_at_ucsf.csv", sep = "\t")

long_wide_pras_location <- pras_data_with_location %>%
  mutate(#"inf_type_dich"=ifelse(grepl("^S", infectiontype), "S",
    # ifelse(grepl("^A", infectiontype), "A", infectiontype)),
    #"timepoint"=paste(inf_type_dich, timepoint_imm, sep=" "),
    "barcode_location"=paste(barcode, box_row_column),
    "in_pilot"=ifelse(id %in% musical_metadata$combined_id, TRUE, FALSE),
    timepoint_imm=factor(timepoint_imm, levels=c(-1, 0, 7, 14, -2)),
    infectiontype=factor(infectiontype, levels=c("S", "A", "S2", "A2", "PS", "CHECK SAMPLE")),
    cohortid=id)%>%
  filter(infectiontype!="NM", in_pilot==FALSE)%>%
  group_by(cohortid)%>%
  arrange(infectiontype, timepoint_imm)%>%
  mutate(new_column = row_number())%>%
  dplyr::select(cohortid, masterbarcode, date, barcode, barcode_location, new_column)



no_location <- long_wide_pras_location%>%
  filter(barcode_location =="NA NA")

no_location$barcode <- slim_metadata$barcode[match(no_location$masterbarcode, slim_metadata$Masterbarcode)]

subset(no_location, no_location$masterbarcode %notin% ucsf_void$V3)