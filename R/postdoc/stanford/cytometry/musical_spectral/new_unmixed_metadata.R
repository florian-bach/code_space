library(dplyr)

omiq_metadata <- read.csv("~/Downloads/OMIQ_metadata-MUSICAL new unmixed 3.csv", header=T)
slim_omiq_metadata <- omiq_metadata%>%
  select(OmiqID, Filename, OriginalFileName, Upload.Folder)

kylie_metadata <- readxl::read_excel("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/Aurora Spectral Flow/For Arefin/MUSICALSpectralFlowMetadata.xlsx")

ungated_omiq_metadata <- kylie_metadata%>%
  mutate(OriginalFileName=if_else(experiment!="mus5", paste(toupper(experiment), " unmixed new_", `Sample:`, sep=""), `Sample:`),
         infectiontype=case_when(grepl("^a", timepoint)~"A",
                                 grepl("^s", timepoint)~"S",
                                 grepl("^n", timepoint)~"NM"),
         control=if_else(grepl("^lrs", id), TRUE, FALSE),
         timepoint=if_else(control==TRUE, "lrs", timepoint))%>%
  select(-`Sample:`)%>%
  left_join(slim_omiq_metadata, ., by="OriginalFileName")

write.csv(ungated_omiq_metadata, "~/postdoc/stanford/cytometry/spectral/MUSICAL/big_experiment/ungated/new_unmixed_omiq_metadata_no_duds_with5.csv", row.names = F)


musical_panel <- readxl::read_excel("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/Aurora Spectral Flow/For Arefin/MUSICALSpectralFlowPanel.xlsx")
musical_panel$class <- "type"
musical_panel$Fluorophore <- paste(musical_panel$Fluorophore, "A", sep="-")

omiq_primaries <- read.csv("~/Downloads/new_unmixed_primary_names.csv", header=F)
omiq_primaries$antigen <- musical_panel$Marker[match(omiq_primaries$V1, musical_panel$Fluorophore)]
write.csv(omiq_primaries, "~/Downloads/new_unmixed_secondary_names.csv")

