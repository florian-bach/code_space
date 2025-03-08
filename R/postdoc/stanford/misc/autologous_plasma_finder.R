


`%notin%` <- Negate(`%in%`)
#mic_drop
mic_drop <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_10/MICDSpecimenBoxOct23_withclinical.dta")

list_of_micdrop_children <- unique(mic_drop$id)
list_of_micdrop_children_minus <- list_of_micdrop_children-10000

#impact
impact <- haven::read_dta("~/postdoc/stanford/clinical_data/IMPACT/IMPACT expanded database through October 31st 2023.dta")

list_of_impact_children <- unique(impact$id)
list_of_impact_children_minus <- list_of_impact_children-20000

#dpsp
dpsp <- haven::read_dta("~/postdoc/stanford/clinical_data/DPSP/DPSPSpecimenBoxOct23_withclinical.dta")

# 
smaller_dpsp <- dpsp %>%
  filter(dod>"2023-01-01")%>%
  select(id, dob, dod, enrolldate)

list_of_mothers <- unique(smaller_dpsp$id)
list_of_children <- unique(c(list_of_impact_children,list_of_micdrop_children))

# list of mothers without children in micdrop or impact who gave birth in 2023
list_of_mothers_without_children <- list_of_mothers[list_of_mothers %notin% list_of_children]


#dpsp specimens
dpsp_specimen <- readxl::read_excel("~/postdoc/stanford/clinical_data/DPSP/tblSpecimenDetails_immunology.xlsx")

small_dpsp_specimen <- dpsp_specimen %>%
  filter(SpecimenType=="PBMC", SubjectID %in% list_of_mothers_without_children)

#dpsp samples in stanfod
dpsp_samples <- readxl::read_excel("~/postdoc/stanford/clinical_data/DPSP/tblBoxDetails_immunology.xlsx")

dpsp_sample_locations <- dpsp_samples %>%
  filter(RandomNumber %in% small_dpsp_specimen$RandomNumber)%>%
  filter(Entrydate>"2023-01-01")

unique(dpsp_sample_locations$BoxNumber)

# picked_samples <- c("P6ZHT", "PLCAF", "PDYNR", "PLDZ7", "PFV48", "P5XME")

picked_samples <- dpsp_sample_locations$RandomNumber[1:16]
#find specimen IDs of PBMC samples picked
try <- dpsp_specimen %>%
  filter(RandomNumber %in% picked_samples)%>%
  select(SpecimenID)
 
#find random number of plasma samples
plasma_rn <- dpsp_specimen %>%
  filter(SpecimenID %in% try$SpecimenID, SpecimenType=="Plasma")%>%
  select(RandomNumber)

# find sample location of plasma samples
dpsp_samples %>%
  filter(RandomNumber %in% plasma_rn$RandomNumber)
