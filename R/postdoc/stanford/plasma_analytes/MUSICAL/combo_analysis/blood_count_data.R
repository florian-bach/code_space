library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")

blood_counts <- readxl::read_xls("~/postdoc/stanford/clinical_data/MUSICAL/PRISMBC LYMPHOCYTE TRUCOUNTS.xls")

slim_blood_counts <- blood_counts%>%
  select(Id, `Draw date`, `CD3+  cells/ul blood`, `lymphocytes/ul blood`)%>%
  filter(!is.na(Id))%>%
  mutate(date=as.Date(`Draw date`), id=as.character(Id))%>%
  select(-`Draw date`, -Id)


clean_data_with_cell_counts <- clean_data %>%
  mutate(id=as.character(substr(id, nchar(id)-2,nchar(id))))%>%
  mutate(date=as.Date(date))%>%
  left_join(., slim_blood_counts, by=c("id", "date"))%>%
  pivot_longer(cols=c( `CD3+  cells/ul blood`, `lymphocytes/ul blood`), names_to = "blood_count", values_to = "cbc_per_ul")

na_counts <- clean_data_with_cell_counts%>%
  filter(is.na(cbc_per_ul))

n_distinct(na_counts$sample_id)

na_counts_for_kenneth <- na_counts%>%
  distinct(id, date, timepoint)

write.csv(na_counts_for_kenneth,"~/Downloads/visits_without_cbc.csv", row.names = F)

good_counts <- clean_data_with_cell_counts%>%
  filter(!is.na(cbc_per_ul))

n_distinct(good_counts$sample_id)


cbc_boxplot <- clean_data_with_cell_counts %>%
  filter(infectiontype %in% c("S", "A"), timepoint %in% c("baseline", "day0", "day7", "day14"))%>%
  ggplot(., aes(x=infectiontype, y=cbc_per_ul, fill=timepoint))+
  geom_boxplot(outliers = F)+
  facet_wrap(~blood_count)+
  theme_minimal()+
  scale_fill_manual(values=viridis::magma(n=5))+
  da_boxplot_theme

ggsave("~/postdoc/stanford/clinical_data/MUSICAL/cbc_boxplot.png", cbc_boxplot, height = 6, width = 8, dpi=444, bg="white")  


cbc_boxplot2 <- clean_data_with_cell_counts %>%
  filter(infectiontype %in% c("S", "A"), timepoint %in% c("baseline", "day0", "day7", "day14"))%>%
  distinct(blood_count, cbc_per_ul, id, timepoint, infectiontype)%>%
  pivot_wider(names_from = blood_count, values_from = cbc_per_ul, id_cols = c("id", "timepoint", "infectiontype"))%>%
  mutate("non T lymphocytes"=`lymphocytes/ul blood`-`CD3+  cells/ul blood`)%>%
  pivot_longer(cols=c(`lymphocytes/ul blood`, `CD3+  cells/ul blood`, `non T lymphocytes`), names_to = "blood_count", values_to = "cbc_per_ul" )%>%
  ggplot(., aes(x=timepoint, y=cbc_per_ul, fill=timepoint))+
  #geom_point()+
  geom_line(aes(group = id))+
  geom_boxplot()+
  theme_minimal()+
  facet_wrap(~blood_count+infectiontype, ncol=2, scales="free")+
  scale_fill_manual(values=viridis::magma(n=5))+
  da_boxplot_theme


clean_data_with_cell_counts %>%
  filter(infectiontype %in% c("S", "A"), timepoint %in% c("baseline", "day0", "day7", "day14"))%>%
  distinct(blood_count, cbc_per_ul, id, timepoint, infectiontype, parasitedensity)%>%
  pivot_wider(names_from = blood_count, values_from = cbc_per_ul)%>%
  mutate("non T lymphocytes"=`lymphocytes/ul blood`-`CD3+  cells/ul blood`)%>%
  pivot_longer(cols=c(`lymphocytes/ul blood`, `CD3+  cells/ul blood`, `non T lymphocytes`), names_to = "blood_count", values_to = "cbc_per_ul" )%>%
  ggplot(., aes(x=parasitedensity+0.1, y=cbc_per_ul+0.1))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  geom_smooth(method="lm")+
  theme_minimal()+
  scale_x_log10()+
  facet_wrap(~blood_count+infectiontype, ncol=2, scales="free")+
  scale_fill_manual(values=viridis::magma(n=5))+
  da_boxplot_theme


# add flow cytometry data ####

jason_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/jason_cytometry_241124.csv")

jason_data <- jason_data %>%
  mutate(id=as.character(cohortid), timepoint_imm=as.numeric(timepoint), infectiontype=case_when(inf_type=="asymp"~"A",
                                                                                               inf_type=="symp"~"S",
                                                                                               inf_type=="nmf"~"NM"))%>%
  select(-timepoint)%>%
  filter(!is.na(cohortid))

combo_data <- clean_data_with_cell_counts%>%
  filter(!is.na(cbc_per_ul))%>%
  inner_join(., jason_data, by = c("id", "timepoint_imm", "infectiontype"))


cell_count_data <- combo_data%>%
  filter(blood_count=="CD3+  cells/ul blood")%>%
  # distinct(id, infectiontype, timepoint,parasitedensity, qpcr, blood_count, cbc_per_ul, CD4_T_Cell_Frequency, Tr1_Frequency)%>%
  mutate(Tr1_count=cbc_per_ul*CD4_T_Cell_Frequency/100*Tr1_Frequency/100,
         Tr1_IFNg_count=Tr1_count*Tr1_IFNg_Frequency/100,
         Tr1_IL10_count=Tr1_count*Tr1_IL10_Frequency/100
         )%>%
  filter(!is.na(Tr1_count), !is.na(Tr1_IFNg_count), !is.na(Tr1_IL10_count))%>%
  pivot_longer(cols=c(Tr1_count, Tr1_IFNg_count, Tr1_IL10_count), names_to = "cell_pop", values_to = "cell_count")%>%
  distinct(id, infectiontype, timepoint,parasitedensity, qpcr, stim, blood_count, cbc_per_ul, cell_count, cell_pop)


tr1_freq_plot <- cell_count_data %>%
  filter(infectiontype %in% c("S", "A"), timepoint!="bad_baseline", stim=="unstim", cell_pop=="Tr1_count")%>%
  ggplot(., aes(x=timepoint, y=cell_count*1000, fill=timepoint))+
  geom_point()+
  geom_line(aes(group=id))+
  geom_boxplot(outliers = F)+
  theme_minimal()+
  # scale_y_log10()+
  ylab("Tr1 cells / mL")+
  facet_wrap(~infectiontype, scales = "free")+
  scale_fill_manual(values=viridis::magma(n=5))
  

tr1_function_plot <- cell_count_data%>%
  filter(infectiontype %in% c("S", "A"), timepoint!="bad_baseline", stim!="unstim", cell_pop!="Tr1_count")%>%
  ggplot(., aes(x=timepoint, y=cell_count*1000, fill=timepoint))+
  geom_point()+
  geom_line(aes(group=id))+
  geom_boxplot(outliers = F)+
  theme_minimal()+
  # scale_y_log10()+
  ylab("cyotkine+ Tr1 cells / mL")+
  facet_wrap(~infectiontype+cell_pop+stim, scales = "free", ncol=4)+
  scale_fill_manual(values=viridis::magma(n=5))

day_
