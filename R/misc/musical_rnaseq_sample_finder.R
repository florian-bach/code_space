metadata <- read.csv("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/fixed_musical_metadata_final.csv")

unclean_nulisa_meta <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/unclean_musical_combo_with_metadata.csv")

prism_metadata <- haven::read_dta("/Users/fbach/Library/CloudStorage/Box-Box/Jagannathan_Lab_Folder/PRISM_Specimens/BorderCohort/PBMC_Samples/PRISMBC_samples_final/PRISMBC_pbmc_plasma_paxgene_final.dta")

#missing day28
musical_rnaseq_samples <- unclean_nulisa_meta%>%
  filter(timepoint %in% c("baseline", "day0") | timepoint=="day14" & infectiontype=="A")%>%
  filter(infectiontype!="NM")%>%
  mutate(date=as.Date(date))%>%
  distinct(id, date, infectiontype, timepoint, timepoint_imm)

#add day28s
asymp_day28s <- musical_rnaseq_samples%>%
  filter(timepoint=="day0", infectiontype=="A")%>%
  select(id, date, infectiontype, timepoint)%>%
  mutate(date=date+28, timepoint="day28", timepoint_imm=28)

non_empty <- prism_metadata %>%
  filter(Freezer !="", aliquot=="1", sampletype=="Paxgene reagent")%>%
  mutate(date=as.Date(date))%>%
  select(id, date, barcode, BoxNumber, boxrow, boxcolumn, aliquot)

combo_musical_rnaseq_samples <- rbind(musical_rnaseq_samples, asymp_day28s)

paxgene_locations <- left_join(combo_musical_rnaseq_samples, non_empty, by=c("id", "date"))

# View(paxgene_locations %>%
#        arrange(id, infectiontype, timepoint))

missing_samples <- paxgene_locations%>%
  filter(is.na(boxrow))%>%
  mutate(missing="yes")

asymp_28_ranges <- musical_rnaseq_samples%>%
  filter(id %in% missing_samples$id)%>%
  filter(timepoint=="day0", infectiontype=="A")%>%
  select(id, date, infectiontype, timepoint)%>%
  mutate(start_range=date+20, end_range=date+35)


# write a for loop to look for missing samples:
found28s <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("id", "date", "barcode", "BoxNumber", "boxrow", "boxcolumn", "aliquot"))

for(i in 1:nrow(asymp_28_ranges)){
  
  d28_ranges <- asymp_28_ranges[i,]
  
  possible_sample <- non_empty %>%
    filter(id==d28_ranges$id, date >= d28_ranges$start_range & date <= d28_ranges$end_range)
  
  found28s[i,] <- possible_sample
}

found28s <- found28s%>%
  mutate(date=as.Date(date))%>%
  mutate(infectiontype="A", timepoint="day28", timepoint_imm=28)

# missing_samples2 <- missing_samples%>%
#   select(colnames(found28s))

all_samples <- rbind(found28s, paxgene_locations)
all_samples2 <- all_samples%>%
  filter(!is.na(BoxNumber))
  
all_day_28s <- all_samples2%>%
  filter(timepoint=="day28")%>%
  filter(!duplicated(id))%>%
  arrange(id)

write.csv(all_day_28s, "~/Downloads/all_day_28s_for_kylie.csv", row.names = F)
