library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

blood_counts <- readxl::read_xls("~/postdoc/stanford/clinical_data/MUSICAL/PRISMBC LYMPHOCYTE TRUCOUNTS.xls")

slim_blood_counts <- blood_counts%>%
  select(Id, `Draw date`, `CD3+  cells/ul blood`, `lymphocytes/ul blood`)%>%
  filter(!is.na(Id))%>%
  mutate(date=as.Date(`Draw date`), id=Id)%>%
  select(-`Draw date`, -Id)

ciber_sort_data <- read.csv("~/Downloads/CIBERSORTx_AllKids_Results.csv")
ciber_sort_data <- ciber_sort_data %>%
  mutate(Novogene.sample.name=Mixture)

rnaseq_metadata <- read.csv("~/Downloads/new_rna-seq_metadata.csv")
rnaseq_metadata <- rnaseq_metadata %>%
  mutate(date = as.Date(date))

combo_data <- rnaseq_metadata %>%
  left_join(slim_blood_counts, by=c("id", "date"))%>%
  left_join(ciber_sort_data, by="Novogene.sample.name")

long_combo_data <- combo_data %>%
  pivot_longer(colnames(combo_data)[33:54], names_to = "cibersort_pop", values_to = "ciber_sort_freq")

long_combo_data%>%
  filter(cibersort_pop %in% c(unique(grep("T.cells", long_combo_data$cibersort_pop, value=T))
                              #unique(grep("B.cells", long_combo_data$cibersort_pop, value=T))
                              ))%>%
  filter(infectiontype %in% c("S", "A"), timepoint %notin% c("bad_baseline", "day28"))%>%
  group_by(Novogene.sample.name)%>%
  reframe(T_cell_ciber_sum=sum(ciber_sort_freq), `CD3+  cells/ul blood`)%>%
  filter(!duplicated(Novogene.sample.name))%>%
  ggplot(., aes(x=T_cell_ciber_sum, y=`CD3+  cells/ul blood`))+
  geom_point()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method = "spearman")+
  theme_minimal()


long_combo_data%>%
  filter(cibersort_pop %in% c(unique(grep("B.cells", long_combo_data$cibersort_pop, value=T))
                              #unique(grep("B.cells", long_combo_data$cibersort_pop, value=T))
  ))%>%
  filter(infectiontype %in% c("S", "A"), timepoint %notin% c("bad_baseline", "day28"))%>%
  group_by(Novogene.sample.name)%>%
  reframe(B_cell_ciber_sum=sum(ciber_sort_freq), `CD3+  cells/ul blood`, `lymphocytes/ul blood`)%>%
  filter(!duplicated(Novogene.sample.name))%>%
  mutate(lymph_minus_T_cbc=`lymphocytes/ul blood`-`CD3+  cells/ul blood`)%>%
  ggplot(., aes(x=B_cell_ciber_sum, y=lymph_minus_T_cbc))+
  geom_point()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method = "spearman")+
  theme_minimal()

long_combo_data%>%
  filter(cibersort_pop %in% c(unique(grep("B.cells", long_combo_data$cibersort_pop, value=T)),
                              unique(grep("T.cells", long_combo_data$cibersort_pop, value=T)),
                              unique(grep("Dendritic.cells", long_combo_data$cibersort_pop, value=T)),
                              unique(grep("Monocytes", long_combo_data$cibersort_pop, value=T)),
                              unique(grep("NK.cells", long_combo_data$cibersort_pop, value=T)),
                              unique(grep("Plasma.cells", long_combo_data$cibersort_pop, value=T))
                              
  ))%>%
  filter(infectiontype %in% c("S", "A"), timepoint %notin% c("bad_baseline", "day28"))%>%
  group_by(Novogene.sample.name)%>%
  reframe(lymph_ciber_sum=sum(ciber_sort_freq), `lymphocytes/ul blood`)%>%
  filter(!duplicated(Novogene.sample.name))%>%
  ggplot(., aes(x=lymph_ciber_sum, y=`lymphocytes/ul blood`))+
  geom_point()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method = "spearman")+
  theme_minimal()


long_combo_data%>%
  filter(cibersort_pop %in% c(unique(grep("T.cells", long_combo_data$cibersort_pop, value=T))
                              #unique(grep("B.cells", long_combo_data$cibersort_pop, value=T))
                              ))%>%
  filter(infectiontype %in% c("S", "A"), timepoint %notin% c("bad_baseline", "day28"))%>%
  group_by(Novogene.sample.name)%>%
  reframe(T_cell_ciber_sum=sum(ciber_sort_freq), `CD3+  cells/ul blood`)%>%
  filter(!duplicated(Novogene.sample.name))%>%
  ggplot(., aes(x=T_cell_ciber_sum, y=`CD3+  cells/ul blood`))+
  geom_point()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method = "spearman")+
  theme_minimal()
  
