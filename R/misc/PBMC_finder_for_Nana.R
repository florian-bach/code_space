library(tidyr)
library(dplyr)

# DPSP PBMC finder
`%notin%` <- Negate(`%in%`)

# the lowest priority DPSP samples are samples from secundigravid mothers whose kids were enrolled neither in MICDROP nor impact;
# these clinical trials are linked such that for MICDROP child id = maternal id + 10000; IMPACT child id = maternal id + 20000

# so to ensure that the samples we are picking don't come from women that gave birth to MICDROP / IMPACT children, we read in those
# databases and extract all child IDs and then subset the DPSP database to exclude all IDs that match

# read in data ####
# mic_drop
mic_drop <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/Mar25/MICDSpecimenBoxMar25_withclinical.dta")
list_of_micdrop_children <- unique(mic_drop$id)
list_of_micdrop_children_minus <- list_of_micdrop_children-10000

# impact
impact <- haven::read_dta("~/Library/CloudStorage/Box-Box/IMPACT Study (Busia)/Data/Specimens/Dec24/IMPASpecimenBoxDec24.dta")
list_of_impact_children <- unique(impact$id)
list_of_impact_children_minus <- list_of_impact_children-20000

# dpsp
dpsp <- haven::read_dta("~/Library/CloudStorage/Box-Box/DP+SP study/DPSP_SPECIMENS/DPSP Specimens/Jul24/DPSPSpecimenBoxJul24_withclinical.dta")

#only keep secundigravidsl ditch unecessary columns
smaller_dpsp <- dpsp %>%
  filter(gravidcat==2)%>%
  select(id, dob, dod, enrolldate)

# select IDs of secundigravid mothers
list_of_mothers <- unique(smaller_dpsp$id)
# subset list of secundis to exclude mothers of other study kids
list_of_mothers_of_impact_or_micrdrop_children <- unique(c(list_of_micdrop_children_minus,list_of_impact_children_minus))

# list of mothers without children in micdrop or impact who gave birth in 2023
list_of_mothers_without_impact_or_micrdrop_children <- list_of_mothers[list_of_mothers %notin% list_of_mothers_of_impact_or_micrdrop_children]


# read in sample location databases
#dpsp specimens
dpsp_specimen <- readxl::read_excel("~/Library/CloudStorage/Box-Box/DP+SP study/DPSP_SPECIMENS/DPSP Vaginal and Stool Samples/Databases/tblSpecimenDetails_06.17.2024.xlsx")

small_dpsp_specimen <- dpsp_specimen %>%
  filter(SpecimenType=="PBMC", SubjectID %in% list_of_mothers_without_impact_or_micrdrop_children)

#dpsp samples in stanfod
dpsp_samples <- readxl::read_excel("~/Library/CloudStorage/Box-Box/DP+SP study/DPSP_SPECIMENS/DPSP Vaginal and Stool Samples/Databases/tblBoxDetails_06.17.2024.xlsx")

dpsp_sample_locations <- dpsp_samples %>%
  filter(RandomNumber %in% small_dpsp_specimen$RandomNumber)%>%
  filter(grepl("^DPSP-3", BoxNumber))%>%
  # filter(Entrydate>"2023-01-01")%>%
  arrange(BoxNumber)

# the box number can then be found on the freezer sheet in Scar or Rafiki LN
