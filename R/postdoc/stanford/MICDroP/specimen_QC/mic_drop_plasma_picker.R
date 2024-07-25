library(tidyr)
library(dplyr)
library(xlsx)
library(ggplot2)

`%notin%` <- Negate(`%in%`)
is.blank <- function(x){sapply(x, function(y) {ifelse(y=="", TRUE, FALSE)})}

raw_data <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_12/MICDSpecimenBoxDec23_withclinical.dta")


sample_ages <- c(8, 24, 52, 68, 84, 104, 120)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1



long_specimen_data <- raw_data %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7)%>%
  # select(id, dob, date, flo_age_in_wks, mstatus, qPCRparsdens, ageinwks, SampleDate, PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR, visittype, withdrawaldate) %>%
  select(id, dob, date, flo_age_in_wks, mstatus,pardens, qPCRparsdens, ageinwks, SampleDate, PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR, visittype, withdrawaldate, starts_with(c("RandomNumber", "BoxNumber", "PositionColumn", "PositionRow"))) %>%
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
  )



plasma_boxes_in_stanford_april24 <- paste("MICD", c(1001:1026, 1028, 1029, 1031:1038), sep="-")

#BoxNumber3 = plasmas
withdrawn_plasmas <- long_specimen_data %>%
  filter(Specimen_Type == "Plasma",
    BoxNumber3 %in% plasma_boxes_in_stanford_april24,
    !is.na(withdrawaldate))%>%
  group_by(id)%>%
  add_count()
  
samples_to_pick <- withdrawn_plasmas%>%
  filter(n>=3)%>%
  select(id, ageinwks, SampleDate, RandomNumber3, BoxNumber3, PositionColumn3, PositionRow3)

samples_for_nulisa <- samples_to_pick %>%
  filter(id %in% c(10766, 10794, 10842))

kids_with_samples <- distinct(samples_to_pick, id)
# table(withdrawn_plasmas$ageinwks)                  
# 8  9 24 52 
# 31  2 21  1 

kids_with_nulisa <- c(10766, 10794, 10842)

nulisa_pilot_qpcr_parasitemias <- long_specimen_data %>%
  filter(id %in% kids_with_nulisa, !is.na(mstatus))%>%
  ggplot(aes(x=flo_age_in_wks, y=qPCRparsdens+0.1, color=factor(mstatus)))+
  geom_vline(xintercept = c(8, 24, 52))+
  geom_point()+
  scale_y_log10(limits=c(0.1, 10e5))+
  facet_grid(~id)+
  xlab("age in weeks")+
  ylab("parasites / μL")+
  annotation_logticks(sides = "l")+
  scale_color_manual(values=c("black", "darkred"))+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/nulisa_pilot_qpcr_parasitemias.png", nulisa_pilot_qpcr_parasitemias, width=8, height=4, bg="white", dpi=444)



nulisa_pilot_smear_parasitemias <- long_specimen_data %>%
  filter(id %in% kids_with_nulisa, !is.na(mstatus))%>%
  ggplot(aes(x=flo_age_in_wks, y=pardens+0.1, color=factor(mstatus)))+
  geom_vline(xintercept = c(8, 24, 52))+
  geom_point()+
  scale_y_log10(limits=c(0.1, 10e5))+
  facet_grid(~id)+
  xlab("age in weeks")+
  ylab("parasites / μL")+
  annotation_logticks(sides = "l")+
  scale_color_manual(values=c("black", "darkred"))+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/nulisa_pilot_smear_parasitemias.png", nulisa_pilot_smear_parasitemias, width=8, height=4, bg="white", dpi=444)
