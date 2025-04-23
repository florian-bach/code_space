dpsp_clinical <- haven::read_dta("~/postdoc/stanford/clinical_data/DPSP/DPSPSpecimenBoxOct23_withclinical.dta")

all_samples <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/all_samples.csv")
all_kids <- unique(all_samples$id[nchar(all_samples$id)==5])

maternal_columns <- c("gravidcat")
  maternal_malaria <- dpsp_clinical %>%
  filter(id %in% (all_kids-10000))%>%
  select(id, )
  