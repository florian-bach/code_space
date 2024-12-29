library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(ComplexHeatmap)

`%notin%` <- Negate(`%in%`)

# read in Nulisa and basic meta data ####

nulisa_part1 <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/pilot/nulisa_data.csv")
metadata_part1 <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/pilot_metadata.csv")

nulisa_part2 <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/extra_24/12-19-24_NULISAseq_Inflammation_FlorianBach_24 samples.csv", check.names = F)
nulisa_part2 <- nulisa_part2[,1:25]
metadata_part2 <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/extra_24/extra_24_sample_manifest_for_tran.csv")


long_nulisa_part1 <- nulisa_part1 %>%
  pivot_longer(cols = colnames(nulisa_part1)[2:ncol(nulisa_part1)], names_to = "plasma.barcode", values_to = "concentration")%>%
  filter(plasma.barcode %in% metadata_part1$plasma.barcode, plasma.barcode!="PD50")%>%
  left_join(metadata_part1, by="plasma.barcode")%>%
  select(-X, -plasma.barcode)%>%
  mutate("study"="pilot")
  

long_nulisa_part2 <- nulisa_part2 %>%
  pivot_longer(cols = colnames(nulisa_part2)[2:ncol(nulisa_part2)], names_to = "sample_id", values_to = "concentration")%>%
  mutate(id=as.integer(substr(sample_id, 1, 5)),
         Timepoint_in_weeks=substr(sample_id, nchar(sample_id)-1, nchar(sample_id)),
         Timepoint_in_weeks=as.integer(gsub("_", "", Timepoint_in_weeks)))%>%
  left_join(metadata_part2, by=c("id", "Timepoint_in_weeks"))%>%
  mutate("timepoint"=paste(Timepoint_in_weeks, "w", sep=""),
         "study"="extra 24")%>%
  select(colnames(long_nulisa_part1))
  

combo_nulisa <- rbind(long_nulisa_part1, long_nulisa_part2)%>%
  mutate(sample_id=paste(id, timepoint, sep="_"))%>%
  group_by(targetName) %>%
  mutate(z_conc=scale(concentration, center = TRUE, scale = TRUE))%>%
  group_by(sample_id)%>%
  mutate(mean_z_conc=mean(z_conc))
# add more metadata ####
raw_data <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2024_06/MICDSpecimenBoxJun24_withclinical.dta")

# for finding putative sampling visits we'll filter the database to only include visits around the proper sampling timepoint dates, with a week plus/minus
sample_ages <- c(8, 24, 52, 68, 84, 104, 120)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1

sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))

long_specimen_data <- raw_data %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7)%>%
  select(id, dob, date, flo_age_in_wks, mstatus, qPCRparsdens, ageinwks, SampleDate, starts_with(c("BoxNumber", "PositionColumn", "PositionRow", "RandomNumber")), PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR, visittype, withdrawaldate) %>%
  mutate("visit_id"=paste(id, date, sep="_"))%>%
  pivot_longer(cols = c(PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR), names_to = "Specimen_Type", values_to = "Specimen_ID")%>%
  mutate(subject_id=id)%>%
  #Specimen_IDs are shared between specimen types, so let's create a unique code
  mutate(Specimen_ID_ID=paste(Specimen_Type, visit_id, sep="_"))%>%
  mutate("Timepoint_in_weeks"=if_else(
    flo_age_in_wks %in% sample_ages, flo_age_in_wks, ifelse(
      flo_age_in_wks %in% sample_ages_minus, flo_age_in_wks+1, if_else(
        flo_age_in_wks %in% sample_ages_plus, flo_age_in_wks-1, 999)
    )
  )
  )%>%
  mutate("timepoint"=paste(Timepoint_in_weeks, "w", sep=""))

extra_individuals <- c(11685, 10501, 10857, 11831, 11721, 11462, 11622, 11343)

longer_timecourses <- c(11343, 10857)

slim_micdrop_metadata <- long_specimen_data %>%
  dplyr::filter(id %in% unique(combo_nulisa$id))%>%
  dplyr::filter(flo_age_in_wks < 54 | id %in% longer_timecourses,
                Specimen_ID!="", Specimen_Type %in% c("Plasma"))%>%
  mutate("timepoint"=paste(Timepoint_in_weeks, "w", sep=""),
         id=as.character(id))%>%
  select(id, timepoint, mstatus,qPCRparsdens)

full_data <- combo_nulisa%>%
  left_join(., slim_micdrop_metadata, by=c("id", "timepoint"))%>%
  mutate(timepoint=factor(timepoint, levels=c("8w", "24w", "52w", "68w", "84w")),
         infectiontype=case_when(mstatus==1 ~"S",
                                 mstatus==0 & qPCRparsdens>0 ~ "A"))


full_data%>%
  filter(targetName %in% c("IFNG", "CX3CL1", "CRP", "IL6", "IL12B", "IL10", "IL27", "LAG3", "GZMA"))%>%
  ggplot(., aes(x=timepoint, y=concentration, fill=factor(mstatus)))+
  geom_boxplot()+
  # geom_point(aes(color=factor(mstatus)), position = position_dodge(width = 0.5))+
  facet_wrap(~targetName, scales = "free")+
  # viridis::scale_fill_viridis(discrete = TRUE)+
  scale_fill_manual(values=c("black", "darkred"))+
  scale_color_manual(values=c("black", "darkred"))+
  theme_minimal()+
  theme(legend.position = "bottom")


full_data%>%
  filter(targetName %in% c("IFNG", "IL10", "LAG3", "GZMA"))%>%
  filter(timepoint %in% c("8w", "24w", "52w"))%>%
  ggplot(., aes(x=timepoint, y=concentration, fill=factor(mstatus)))+
  geom_boxplot()+
  geom_point(aes(color=factor(mstatus)), position = position_dodge(width = 0.75))+
  facet_wrap(~targetName, scales = "free")+
  # viridis::scale_fill_viridis(discrete = TRUE)+
  scale_fill_manual(values=c("black", "darkred"))+
  scale_color_manual(values=c("black", "darkred"))+
  theme_minimal()+
  theme(legend.position = "none")

full_data%>%
  filter(targetName %in% c("IL7R", "CCL16", "CCL24", "LILRB2"))%>%
  filter(timepoint %in% c("8w", "24w", "52w"))%>%
  ggplot(., aes(x=timepoint, y=concentration, fill=factor(mstatus)))+
  geom_boxplot()+
  geom_point(aes(color=factor(mstatus)), position = position_dodge(width = 0.75))+
  facet_wrap(~targetName, scales = "free")+
  # viridis::scale_fill_viridis(discrete = TRUE)+
  scale_fill_manual(values=c("black", "darkred"))+
  scale_color_manual(values=c("black", "darkred"))+
  theme_minimal()+
  theme(legend.position = "none")


full_data%>%
  filter(targetName %in% c("CX3CL1"))%>%
  ggplot(., aes(x=timepoint, y=concentration, fill=factor(mstatus)))+
  # geom_boxplot()+
  geom_point(aes(color=factor(mstatus)), position = position_dodge(width = 0.75))+
  facet_wrap(~id, scales = "free")+
  # viridis::scale_fill_viridis(discrete = TRUE)+
  scale_fill_manual(values=c("black", "darkred"))+
  scale_color_manual(values=c("black", "darkred"))+
  theme_minimal()+
  theme(legend.position = "none")
