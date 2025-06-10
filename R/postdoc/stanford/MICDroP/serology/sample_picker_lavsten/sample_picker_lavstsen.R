library(tidyr)
library(dplyr)
library(xlsx)
library(ggplot2)

`%notin%` <- Negate(`%in%`)
is.blank <- function(x){sapply(x, function(y) {ifelse(y=="", TRUE, FALSE)})}

# all except lavstsen plasmas were done with this
# raw_data <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2024_06/MICDSpecimenBoxJun24_withclinical.dta")


raw_data <-  haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/Mar25/MICDSpecimenBoxMar25_withclinical.dta")


# for finding putative sampling visits we'll filter the database to only include visits around the proper sampling timepoint dates, with a week plus/minus
sample_ages <- c(8, 24, 52, 68, 84, 104, 120)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1

sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))

#turn the data into long format, include a couple of convenience variables, drop irrelevant columns
long_specimen_data <- raw_data %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7)%>%
  mutate(stool=ifelse(stool==1, "ordered", "not"))%>%
  select(id, dob, date, flo_age_in_wks, mstatus, qPCRparsdens, ageinwks, SampleDate, starts_with(c("BoxNumber", "PositionColumn", "PositionRow", "RandomNumber")), PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, stool, qPCR, visittype, withdrawaldate) %>%
  mutate("visit_id"=paste(id, date, sep="_"))%>%
  pivot_longer(cols = c(PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR, stool), names_to = "Specimen_Type", values_to = "Specimen_ID")%>%
  mutate(subject_id=id)%>%
  #Specimen_IDs are shared between specimen types, so let's create a unique code
  mutate(Specimen_ID_ID=paste(Specimen_Type, visit_id, sep="_"))%>%
  mutate("Timepoint_in_weeks"=if_else(
    flo_age_in_wks %in% sample_ages, flo_age_in_wks, ifelse(
      flo_age_in_wks %in% sample_ages_minus, flo_age_in_wks+1, if_else(
        flo_age_in_wks %in% sample_ages_plus, flo_age_in_wks-1, 999)
    )
  )
  )

msd_map <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/msd_map_edit.csv")
non_msd_map <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/non_msd_edit.csv")

slim_msd_map <- msd_map%>%
  mutate(reorder.column=rep_len(1:9, nrow(.)),
         reorder.row=rep_len(rep(1:9,each=9), nrow(.)),
         reorder.box=rep_len(rep(1:9, each=81), nrow(.)))%>%
  select(id, date, Timepoint_in_weeks, reorder.column, reorder.row, reorder.box)
  
slim_non_msd_map <- non_msd_map%>%
  mutate("Timepoint_in_weeks"=flo_age_in_wks)%>%
  mutate(reorder.box=reorder.box+3)%>%
  select(id, date, Timepoint_in_weeks, reorder.column, reorder.row, reorder.box)


combo_map <- rbind(slim_non_msd_map, slim_msd_map)

write.csv(combo_map, "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/existing_aliquot_map.csv")

grant_list_100 <- readxl::read_excel("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/Samples selected for Florian April 17 2025.xls")

# samples to plate for sweden ####
swedes <- combo_map%>%
  filter(Timepoint_in_weeks %in% c(51, 52))%>%
  arrange(reorder.box)%>%
  mutate("plate"=rep_len(rep(1:9, each=96), nrow(.)),
         "plate row"=rep_len(rep(LETTERS[1:8],each=12), nrow(.)),
         "plate column"=rep_len(1:12, nrow(.)))


  

# samples to plate for denmark
danes_plate <- combo_map%>%
  filter(id %in% grant_list_100$id)%>%
  arrange(reorder.box)%>%
  mutate("plate"=rep_len(rep(1:9, each=96), nrow(.)),
         "plate row"=rep_len(rep(LETTERS[1:8],each=12), nrow(.)),
         "plate column"=rep_len(1:12, nrow(.)))



write.csv(swedes, "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/swedes_for_plating.csv")
write.csv(danes_plate, "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/danes_for_plating.csv")



# samples to pick for denmark
danes_pick <- long_specimen_data %>%
  filter(id %in% grant_list_100$id, Specimen_Type=="Plasma")%>%
  filter(Specimen_ID!="", flo_age_in_wks<106, flo_age_in_wks>102, flo_age_in_wks %in% sample_ranges)%>%
  select(id, date, RandomNumber3, Timepoint_in_weeks, BoxNumber3, PositionColumn3, PositionRow3)%>%
  arrange(BoxNumber3)%>%
  mutate("new column"=rep_len(1:9, nrow(.)),
         "new row"=rep_len(rep(1:9,each=9), nrow(.)),
         "new box"=rep_len(rep(c("104 Box 1", "104 Box 1"), each=81), nrow(.)))

write.csv(danes_pick, "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/104_week_timepoints_to_pick.csv", row.names = F)

label_maker <- danes_pick%>%
  select(id, date, RandomNumber3)

label_maker_with_space <- label_maker[rep(1:nrow(label_maker), each = 2), ]
label_maker_with_space$id <- as.character(label_maker_with_space$id)
label_maker_with_space$date <- as.character(label_maker_with_space$date)
label_maker_with_space[1:nrow(label_maker_with_space) %% 2 == 0, ] <- " "

write.csv(label_maker_with_space, "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/104_week_timepoints_label_maker_with_space.csv", row.names = F)


# plate map for aliquoting ####

slim_danes_pick <- danes_pick%>%
  select(id, date, `new box`, `new column`, `new row`)%>%
  mutate(date=as.character(date))

slim_danes_plate <- danes_plate%>%
  mutate(`new box`=reorder.box, `new column`=reorder.column, `new row`=reorder.row)%>%
  select(id, date, `new box`, `new column`, `new row`)%>%
  mutate(`new box`=as.character(`new box`))


day_one_danes <- slim_danes_plate%>%
  filter(`new box`%in%1:4)%>%
  distinct(id)


plate_for_danes_day_one <- bind_rows(slim_danes_pick, slim_danes_plate)%>%
  filter(id %in% day_one_danes$id)%>%
  arrange(id, date)%>%
  mutate("dane plate"=rep_len(rep(1:9, each=96), nrow(.)),
         "plate row"=rep_len(rep(LETTERS[1:8],each=12), nrow(.)),
         "plate column"=rep_len(1:12, nrow(.)))

colnames(plate_for_danes_day_one)[3:5] <- c("box", "box column", "box row") 
write.csv(plate_for_danes_day_one, "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/plate_plan_day_one.csv", row.names = F)


day_two_danes <- slim_danes_plate%>%
  filter(`new box`>4)%>%
  distinct(id)

plate_for_danes_day_two <- bind_rows(slim_danes_pick, slim_danes_plate)%>%
  filter(id %in% day_two_danes$id)%>%
  arrange(id, date)%>%
  mutate("dane plate"=rep_len(rep(3:9, each=96), nrow(.)),
         "plate row"=rep_len(rep(LETTERS[1:8],each=12), nrow(.)),
         "plate column"=rep_len(1:12, nrow(.)))

colnames(plate_for_danes_day_two)[3:5] <- c("box", "box column", "box row") 
write.csv(plate_for_danes_day_two, "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/plate_plan_day_two.csv", row.names = F)



