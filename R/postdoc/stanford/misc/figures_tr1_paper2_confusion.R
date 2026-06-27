# correlation figure #####
nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/tr1_paper/revised_baseline_clean_musical_combo_with_metadata.csv")
nulisa_data <- nulisa_data%>%
  mutate(infectiontype=substr(infectiontype, 1, 1))


slim_nulisa_data <- nulisa_data %>%
  mutate(id=as.character(id))%>%
  select(id, date, timepoint, timepoint_imm, infectiontype, targetName, concentration, log_qpcr, ageyrs, gender_categorical)%>%
  filter(timepoint!="bad_baseline")%>%
  mutate(timepoint_imm=if_else(timepoint=="baseline", -1, timepoint_imm))

slim_cell_count_data <- read.csv("~/Downloads/df_jason_analysis (1).csv")

slim_cell_count_data <- slim_cell_count_data%>%
  pivot_longer(cols = ends_with("Frequency"), names_to = "gate", values_to = "freq")%>%
  mutate(id=cohortid, timepoint_imm=timepoint)%>%
  mutate(infectiontype=substr(infectiontype, 1, 1))%>%
  mutate("timepoint"=case_when(
    timepoint_imm==-1~"baseline",
    timepoint_imm==0~"day0",
    timepoint_imm==7~"day7",
    timepoint_imm==14~"day14",
    timepoint_imm==28~"day28"))%>%
  select(id, date, timepoint, timepoint_imm, infectiontype, gate, stim, freq)


long_combo <- slim_nulisa_data%>%
  mutate(id=as.integer(id))%>%
  full_join(., slim_cell_count_data, by = c("id", "infectiontype", "timepoint_imm"))%>%
  mutate(timepoint=timepoint.x)%>%
  select(-timepoint.x, -timepoint.y)%>%
  filter(!is.na(freq))

target_cor_plot <- long_combo%>%
  filter(gate=="Tr1_Frequency", stim=="unstim", infectiontype%in%c("A", "S"), targetName %in%  c("IL10", "LAG3", "GZMA", "IFNG"))%>%
  distinct(id, timepoint, infectiontype, targetName, gate, freq, stim, concentration)%>%
  ggplot(., aes(x=freq, y=concentration))+
  geom_point(aes(fill=timepoint, shape=infectiontype), stroke=0.1)+
  ggpubr::stat_cor(method="spearman")+
  geom_smooth(method="lm", color="black")+
  xlab("Tr1% of non-naive CD4")+
  scale_fill_manual(values=time_cols)+
  scale_shape_manual(values=c(21, 25))+
  facet_wrap(~targetName, scales = "free")+
  guides(fill=guide_none())+
  theme_minimal()
ggsave("~/Downloads/actual_target_cor_plot.png", target_cor_plot, width=5.33333, height=5.33333, bg="white")

# cell count figure ####
metadata <- read.csv("~/Downloads/musical_metadata_new_qpcr.csv")


metadata <- metadata%>%
  mutate("timepoint"=case_when(timepoint_imm==-2 & id %in% c(134, 164, 176, 317, 331, 363, 379, 398, 402, 442, 577) ~"baseline",
                               timepoint_imm==-1 & id %notin% c(134, 164, 176, 317, 331, 363, 379, 398, 402, 442, 577) ~"baseline",
                               timepoint_imm==-1 & infectiontype %in% c("A", "S") & id %in% c(134, 164, 176, 317, 331, 363, 379, 398, 402, 442, 577) ~"bad_baseline",
                               timepoint_imm==0~"day0",
                               timepoint_imm==7~"day7",
                               timepoint_imm==14~"day14",
                               timepoint_imm==28~"day28"))
blood_counts <- readxl::read_xls("~/postdoc/stanford/clinical_data/MUSICAL/PRISMBC LYMPHOCYTE TRUCOUNTS.xls")

slim_blood_counts <- blood_counts%>%
  select(Id, `Draw date`, `CD3+  cells/ul blood`, `lymphocytes/ul blood`)%>%
  filter(!is.na(Id))%>%
  mutate(date=as.Date(`Draw date`), id=as.character(Id))%>%
  select(-`Draw date`, -Id)


clean_data_with_cell_counts <- metadata %>%
  mutate(date=as.Date(date), id=as.character(id))%>%
  left_join(., slim_blood_counts, by=c("id", "date"))%>%
  pivot_longer(cols=c( `CD3+  cells/ul blood`, `lymphocytes/ul blood`), names_to = "blood_count", values_to = "cbc_per_ul")%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))

jason_data <- read.csv("~/Downloads/df_jason_analysis (1).csv")

jason_data <- jason_data %>%
  mutate(id=as.numeric(cohortid),
         "timepoint"=case_when(
           timepoint==-1~"baseline",
           timepoint==0~"day0",
           timepoint==7~"day7",
           timepoint==14~"day14",
           timepoint==28~"day28"))%>%
  filter(!is.na(cohortid))%>%
  select(-c(X.2, index, Unnamed..0, X))

# write.csv(jason_data, "~/postdoc/stanford/plasma_analytes/MUSICAL/combo/mus1_mus8_flow_data_edit.csv")

combo_data <- clean_data_with_cell_counts%>%
  mutate(id=as.numeric(id))%>%
  filter(!is.na(cbc_per_ul))%>%
  right_join(., jason_data, by = c("id", "timepoint", "infectiontype"))%>%
  mutate("timepoint"=case_when(timepoint_imm==-2 & id %notin% c("176", "363", "577") ~"bad_baseline",
                               timepoint_imm==-2 & id %in% c("176", "363", "577") ~"baseline",
                               timepoint_imm==-1~"baseline",
                               timepoint_imm==0~"day0",
                               timepoint_imm==7~"day7",
                               timepoint_imm==14~"day14",
                               timepoint_imm==28~"day28"))


cell_count_data <- combo_data%>%
  filter(blood_count=="CD3+  cells/ul blood" | is.na(blood_count))%>%
  mutate("absolute_Tr1"=cbc_per_ul*CD4_T_Cell_Frequency/100* Memory_CD4_T_Cell_Frequency/100 *Tr1_Frequency/100)%>%
  mutate("absolute_CD4_memory"=cbc_per_ul*CD4_T_Cell_Frequency/100* Memory_CD4_T_Cell_Frequency/100)%>%
  pivot_longer(cols = ends_with("Frequency"), names_to = "gate", values_to = "freq")%>%
  mutate(count = ifelse(grepl("Tr1", gate), absolute_Tr1*freq/100, 
                        ifelse(gate %in% c("IFNg_Frequency","IL10_Frequency", "IL21_Frequency"), absolute_CD4_memory*freq/100, NA)))

# write.csv(cell_count_data, "~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/t_cell_cytometry_count_data.csv", row.names = F)

tr1_counts_plot <- cell_count_data %>%
  filter(infectiontype %in% c("S"), timepoint!="bad_baseline", stim=="unstim")%>%
  distinct(id, infectiontype, timepoint, absolute_Tr1)%>%
  ggplot(., aes(x=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")), y=absolute_Tr1, fill=timepoint))+
  geom_line(color="grey",aes(group=id))+
  geom_boxplot(outliers = F)+
  # scale_y_log10()+
  ylab("Tr1 cells / μL")+
  facet_wrap(~infectiontype, scales = "free_x")+
  scale_fill_manual(values=time_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")

ggsave("~/Downloads/actual_tr1_count_plot.png", tr1_counts_plot, width =2.666667, height=3, bg="white")



target_cor_plot2 <- long_combo%>%
  filter(gate=="Th1_Frequency", stim=="unstim", infectiontype%in%c("A", "S"), targetName %in%  c("IL10", "LAG3", "GZMA", "IFNG"))%>%
  distinct(id, timepoint, infectiontype, targetName, gate, freq, stim, concentration)%>%
  ggplot(., aes(x=freq, y=concentration))+
  geom_point(aes(fill=timepoint, shape=infectiontype), stroke=0.1)+
  ggpubr::stat_cor(method="spearman")+
  geom_smooth(method="lm", color="black")+
  xlab("Th1% of non-naive CD4")+
  scale_fill_manual(values=time_cols)+
  scale_shape_manual(values=c(21, 25))+
  facet_wrap(~targetName, scales = "free")+
  # guides(fill=guide_none())+
  theme_minimal()

ggsave("~/Downloads/actual_correlation_th1_proteins.png", target_cor_plot2, width=5.33333, height=5.33333, bg="white")

