# data generation ####

cell_count_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/df_jason_analysis.csv")

gates_of_interest <- c("Tr1_Frequency", "Memory_CD4_T_Cell_Frequency", "Th1_Frequency", "Th2_Frequency", "Th17_Frequency", "Tregs_Frequency",
                       "Tr1_IFNg_Frequency", "Tr1_IL10_Frequency","Tr1_IL21_Frequency",
                       "IFNg_Frequency", "IL10_Frequency", "IL21_Frequency",
                       "cTfh_IFNg_Frequency", "cTfh_IL10_Frequency" , "cTfh_IL21_Frequency",
                       "Th1_IFNg_Frequency", "Th1_IL10_Frequency", "Th1_IL21_Frequency",
                       "Tregs_IFNg_Frequency", "Tregs_IL10_Frequency", "Tregs_IL21_Frequency")

cell_subset_gates <- c("Tr1_Frequency", "Tregs_Frequency", "Memory_CD4_T_Cell_Frequency", "Th1_Frequency", "Th2_Frequency", "Th17_Frequency", "Tregs_Frequency")

cytokine_gates <- c("Tr1_IFNg_Frequency", "Tr1_IL10_Frequency","IFNg_Frequency",
                    "IFNg_Frequency", "IL10_Frequency", "IL21_Frequency",
                    "cTfh_IFNg_Frequency", "cTfh_IL10_Frequency" , "cTfh_IL21_Frequency",
                    "Th1_IFNg_Frequency", "Th1_IL10_Frequency", "Th1_IL21_Frequency",
                    "Tregs_IFNg_Frequency", "Tregs_IL10_Frequency", "Tregs_IL21_Frequency")

treg_gates <- c("Tregs_IFNg_Frequency", "Tregs_IL10_Frequency", "Tregs_IL21_Frequency")
cd4_gates <- c("IFNg_Frequency", "IL10_Frequency", "IL21_Frequency")
tr1_gates <- c("Tr1_IFNg_Frequency", "Tr1_IL10_Frequency","Tr1_IL21_Frequency")
th1_gates <- c("Th1_IFNg_Frequency", "Th1_IL10_Frequency", "Th1_IL21_Frequency")

tfh_gates <- c("cTfh_IFNg_Frequency", "cTfh_IL10_Frequency", "cTfh_IL21_Frequency")

slim_cell_count_data <- cell_count_data%>%
  pivot_longer(cols = ends_with("Frequency"), names_to = "gate", values_to = "freq")%>%
  mutate(timepoint_imm=timepoint)%>%
  mutate("timepoint"=case_when(timepoint_imm==-2 & id %notin% c(176, 363, 577) ~"bad_baseline",
                               timepoint_imm==-2 & id %in% c(176, 363, 577) ~"baseline",
                               timepoint_imm==-1~"baseline",
                               timepoint_imm==0~"day0",
                               timepoint_imm==7~"day7",
                               timepoint_imm==14~"day14",
                               timepoint_imm==28~"day28"))%>%
  mutate(timepoint=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  select(id, timepoint, timepoint_imm, MusicalID, infectiontype, gate, stim, freq)%>%
  filter(gate %in% gates_of_interest, infectiontype %in% c("S", "A"))


# overview plots ####
## subset data ####

subset_freq_overview <- slim_cell_count_data %>%
  filter(gate %in% cell_subset_gates)%>%
  ggplot(., aes(x=timepoint, y=freq/100, fill=infectiontype))+
  geom_boxplot(outliers = F)+
  facet_wrap(~gate+stim, scales="free", ncol = 3)+
  scale_fill_manual(values=c("pink", "darkred"))+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  theme(axis.title = element_blank())

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/musical_t_cell_reanalysis/subset_freq_overview.png", subset_freq_overview, width=8, height=16, bg="white")

## cd4 cytokines ####

cd4_cytokines_overview <- slim_cell_count_data %>%
  filter(gate %in% cd4_gates)%>%
  ggplot(., aes(x=timepoint, y=freq/100, fill=infectiontype))+
  geom_boxplot(outliers = F)+
  facet_wrap(~gate+stim, scales="free", ncol = 3)+
  scale_fill_manual(values=c("pink", "darkred"))+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  theme(axis.title = element_blank())

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/musical_t_cell_reanalysis/cd4_cytokines_overview.png", cd4_cytokines_overview, width=8, height=8, bg="white")

## th1 cytokines ####

th1_cytokines_overview <- slim_cell_count_data %>%
  filter(gate %in% th1_gates)%>%
  ggplot(., aes(x=timepoint, y=freq/100, fill=infectiontype))+
  geom_boxplot(outliers = F)+
  facet_wrap(~gate+stim, scales="free", ncol = 3)+
  scale_fill_manual(values=c("pink", "darkred"))+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  theme(axis.title = element_blank())

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/musical_t_cell_reanalysis/th1_cytokines_overview.png", th1_cytokines_overview, width=8, height=8, bg="white")

## tr1 cytokines ####

tr1_cytokines_overview1 <- slim_cell_count_data %>%
  filter(gate %in% tr1_gates)%>%
  ggplot(., aes(x=infectiontype, y=freq/100, fill=timepoint))+
  geom_boxplot(outliers = F)+
  facet_wrap(~gate+stim, scales="free", ncol = 3)+
  viridis::scale_fill_viridis(discrete = T)+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  theme(axis.title = element_blank())

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/musical_t_cell_reanalysis/tr1_cytokines_overview.png", tr1_cytokines_overview1, width=8, height=8, bg="white")

tr1_cytokines_overview2 <- slim_cell_count_data %>%
  filter(gate %in% tr1_gates)%>%
  ggplot(., aes(x=timepoint, y=freq/100, fill=infectiontype))+
  geom_boxplot(outliers = F)+
  facet_wrap(~gate+stim, scales="free", ncol = 3)+
  scale_fill_manual(values=c("pink", "darkred"))+
  # scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  theme(axis.title = element_blank())

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/musical_t_cell_reanalysis/tr1_cytokines_overview2.png", tr1_cytokines_overview2, width=8, height=8, bg="white")


tr1_cytokines_overview3 <- slim_cell_count_data %>%
  filter(gate %in% tr1_gates, stim!="PMA")%>%
  ggplot(., aes(x=timepoint, y=freq/100, fill=factor(stim, levels=c("unstim", "iRBCs"))))+
  geom_boxplot(outliers = F)+
  facet_wrap(~gate+infectiontype, scales="free", ncol = 2)+
  scale_fill_manual(values=c("darkblue", "darkred"))+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  theme(axis.title = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/musical_t_cell_reanalysis/tr1_cytokines_overview3.png", tr1_cytokines_overview3, width=8, height=8, bg="white")

## treg cytokines ####

treg_cytokines_overview <- slim_cell_count_data %>%
  filter(gate %in% treg_gates)%>%
  ggplot(., aes(x=timepoint, y=freq/100, fill=infectiontype))+
  geom_boxplot(outliers = F)+
  facet_wrap(~gate+stim, scales="free", ncol = 3)+
  scale_fill_manual(values=c("pink", "darkred"))+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  theme(axis.title = element_blank())

treg_cytokines_overview2 <- slim_cell_count_data %>%
  filter(gate %in% tr1_gates, stim %in% c("unstim", "iRBCs"), timepoint=="baseline")%>%
  ggplot(., aes(x=stim, y=freq/100))+
  geom_boxplot(outliers = F)+
  facet_wrap(~gate+timepoint, scales="free", ncol = 3)+
  # scale_fill_manual(values=c("pink", "darkred"))+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  theme(axis.title = element_blank())

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/musical_t_cell_reanalysis/treg_cytokines_overview.png", treg_cytokines_overview, width=8, height=8, bg="white")



# absolute cell counts ####

blood_counts <- readxl::read_xls("~/postdoc/stanford/clinical_data/MUSICAL/PRISMBC LYMPHOCYTE TRUCOUNTS.xls")
metadata <- read.csv("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/musical_metadata_new_qpcr.csv")

slim_blood_counts <- blood_counts%>%
  select(Id, `Draw date`, `CD3+  cells/ul blood`, `lymphocytes/ul blood`)%>%
  filter(!is.na(Id))%>%
  mutate(date=as.Date(`Draw date`), id=as.character(Id))%>%
  select(-`Draw date`, -Id)

clean_data_with_cell_counts <- metadata %>%
  mutate(date=as.Date(date), id=as.character(id))%>%
  left_join(., slim_blood_counts, by=c("id", "date"))%>%
  pivot_longer(cols=c( `CD3+  cells/ul blood`, `lymphocytes/ul blood`), names_to = "blood_count", values_to = "cbc_per_ul")


jason_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/mus1_mus8_flow_data_edit.csv")

combo_data <- clean_data_with_cell_counts%>%
  mutate(id=as.integer(id))%>%
  filter(!is.na(cbc_per_ul))%>%
  right_join(., jason_data, by = c("id", "timepoint_imm", "infectiontype"))%>%
  mutate("timepoint"=case_when(timepoint_imm==-2 & id %notin% c("176", "363", "577") ~"bad_baseline",
                               timepoint_imm==-2 & id %in% c("176", "363", "577") ~"baseline",
                               timepoint_imm==-1~"baseline",
                               timepoint_imm==0~"day0",
                               timepoint_imm==7~"day7",
                               timepoint_imm==14~"day14",
                               timepoint_imm==28~"day28"))

subset_count_data <- combo_data%>%
  filter(blood_count=="CD3+  cells/ul blood" | is.na(blood_count))%>%
  mutate("absolute_Tr1"=cbc_per_ul*Memory_CD4_T_Cell_Frequency/100*Tr1_Frequency/100)%>%
  mutate("absolute_Tregs"=cbc_per_ul*Memory_CD4_T_Cell_Frequency/100*Tregs_Frequency/100)%>%
  mutate("absolute_CD4_Memory"=cbc_per_ul*Memory_CD4_T_Cell_Frequency/100)%>%
  mutate("absolute_Th1"=cbc_per_ul*Memory_CD4_T_Cell_Frequency/100*Th1_Frequency/100)%>%
  mutate("absolute_Th2"=cbc_per_ul*Memory_CD4_T_Cell_Frequency/100*Th2_Frequency/100)%>%
  mutate("absolute_Th17"=cbc_per_ul*Memory_CD4_T_Cell_Frequency/100*Th17_Frequency/100)%>%
  pivot_longer(cols = starts_with("absolute"), names_to = "gate", values_to = "count")%>%
  mutate(timepoint=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))


subset_count_data %>%
  filter(stim=="unstim", infectiontype!="NM")%>% 
  ggplot(., aes(x=timepoint, y=count, color=timepoint))+
  geom_point(alpha=0.5)+
  facet_wrap(~infectiontype+gate, scales="free", ncol = 3)+
  viridis::scale_fill_viridis(discrete = T)+
  scale_y_continuous()+
  theme_minimal()+
  theme(axis.title = element_blank())

subset_count_data %>%
  filter(stim=="unstim", infectiontype!="NM")%>% 
  ggplot(., aes(x=timepoint, y=count, color=timepoint))+
  geom_point(alpha=0.5)+
  facet_wrap(~infectiontype+gate, scales="free", ncol = 3)+
  viridis::scale_fill_viridis(discrete = T)+
  scale_y_continuous()+
  theme_minimal()+
  theme(axis.title = element_blank())
