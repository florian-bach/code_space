

# metadata for gated pops ####

omiq_metadata <- read.csv("~/postdoc/stanford/cytometry/spectral/MUSICAL/big_experiment/OMIQ_metadata-complete MUSICAL.csv", header=T)


kylie_metadata <- readxl::read_excel("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/Aurora Spectral Flow/For Arefin/MUSICALSpectralFlowMetadata.xlsx")
colnames(kylie_metadata)[1] <- "file_name" 

kylie_metadata$live_name <- paste(kylie_metadata$experiment, kylie_metadata$file_name, sep="_")
# kylie_metadata$live_name <- gsub("mus", "livecells_mus", kylie_metadata$live_name, fixed = T)
kylie_metadata$live_name <- gsub(".fcs", "_Live Cells.fcs", kylie_metadata$live_name, fixed = T)


big_omiq_metadata <- omiq_metadata%>%
  mutate(live_name=Filename)%>%
  select(OmiqID, Filename,live_name )%>%
  left_join(., kylie_metadata, by="live_name")%>%
  select(-file_name)%>%
  mutate(infectiontype=case_when(grepl("^a", timepoint)~"A",
                                 grepl("^s", timepoint)~"S",
                                 grepl("^n", timepoint)~"NM"),
         control=if_else(grepl("^lrs", id), TRUE, FALSE))
  

write.csv(big_omiq_metadata, "~/postdoc/stanford/cytometry/spectral/MUSICAL/big_experiment/edited_metadata.csv", row.names = F)

# metadata for ungated pops ####


omiq_metadata <- read.csv("~/postdoc/stanford/cytometry/spectral/MUSICAL/big_experiment/ungated/OMIQ_metadata-ungated MUSICAL.csv", header=T)

kylie_metadata <- readxl::read_excel("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/Aurora Spectral Flow/For Arefin/MUSICALSpectralFlowMetadata.xlsx")

ungated_omiq_metadata <- kylie_metadata%>%
  mutate(Filename=paste(experiment,`Sample:`, sep="_"),
         infectiontype=case_when(grepl("^a", timepoint)~"A",
                               grepl("^s", timepoint)~"S",
                               grepl("^n", timepoint)~"NM"),
         control=if_else(grepl("^lrs", id), TRUE, FALSE),
         timepoint=if_else(control==TRUE, "lrs", timepoint))%>%
  select(-`Sample:`)%>%
  left_join(omiq_metadata, ., by="Filename")

write.csv(ungated_omiq_metadata, "~/postdoc/stanford/cytometry/spectral/MUSICAL/big_experiment/ungated/ungated_omiq_metadata.csv", row.names = F)
