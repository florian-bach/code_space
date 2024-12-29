clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")
clean_data <- clean_data %>%
  mutate(timepoint = factor(timepoint, levels=c("bad_baseline", "baseline", "day0", "day7", "day14", "day28")))


# analytes that vary by parasitemia, but not class; 59/95
# sig_base_para$targetName[sig_base_para$targetName%notin%sig_base_zero_class$targetName]
clean_data %>%
  filter(class %in% c("A", "S"))%>%
  filter(targetName %in% c("CTLA4", "IL2RA", "TNFRSF1A", "FASLG", "CCL4", "GZMA"))%>%
  ggplot(., aes(x=log_qpcr, y=concentration, color=class))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~targetName)+
  scale_color_manual(values=viridis::magma(n=3))+
  theme_minimal()

# analytes that vary by timepoint, but not parasitemia;
try <- sig_base_zero$targetName[sig_base_zero$targetName%notin%sig_base_para$targetName]
clean_data %>%
  filter(class %in% c("A", "S"))%>%
  filter(targetName %in% try[1:6])%>%
  ggplot(., aes(x=log_qpcr, y=concentration, color=class))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~targetName)+
  scale_color_manual(values=viridis::magma(n=3))+
  theme_minimal()

clean_data %>%
  filter(class %in% c("A", "S"), timepoint %notin% c("bad_baseline", "day28", "day7"))%>%
  filter(targetName %in% try[1:6])%>%
  ggplot(., aes(x=class, y=concentration, fill=timepoint))+
  geom_boxplot()+
  # geom_point()+
  # geom_smooth(method="lm")+
  facet_wrap(~targetName)+
  scale_fill_manual(values=viridis::magma(n=3))+
  theme_minimal()

# analytes that vary through time and depend on parasitemia
sig_base_para$targetName[sig_base_para$targetName %in% sig_base_zero$targetName]

clean_data %>%
  filter(class %in% c("A", "S"))%>%
  filter(targetName %in% c("IL10", "CRP", "CXCL10", "LAG3", "IL27", "IFNG"))%>%
  ggplot(., aes(x=log_qpcr, y=concentration, color=class, shape=temperature_cat))+
  geom_point()+
  # geom_smooth(method="lm")+
  facet_wrap(~targetName, scales="free")+
  scale_color_manual(values=viridis::magma(n=3))+
  theme_minimal()

clean_data %>%
  filter(class %in% c("A", "S"), timepoint %notin% c("bad_baseline", "day28", "day7"))%>%
  filter(targetName %in% c("IL10", "CRP", "CXCL10", "LAG3", "IL27", "IFNG"))%>%
  ggplot(., aes(x=class, y=concentration, fill=timepoint))+
  geom_boxplot()+
  # geom_point()+
  # geom_smooth(method="lm")+
  facet_wrap(~targetName)+
  scale_fill_manual(values=viridis::magma(n=3))+
  theme_minimal()

#class vs parasitemia#class vs parasitemiatemperature_cat

clean_data %>%
  filter(class %in% c("A", "S"), timepoint %notin% c("bad_baseline", "day28"))%>%
  select(-targetName)%>%
  filter(!duplicated(sample_id))%>%
  ggplot(., aes(x=timepoint, y=log_qpcr, fill=class))+
  geom_violin(draw_quantiles = seq(0,1,0.25), color="grey")+
  scale_fill_manual(values=viridis::magma(n=3))+
  theme_minimal()

clean_data %>%
  filter(class %in% c("A", "S"), timepoint %notin% c("bad_baseline", "day28"))%>%
  select(-targetName)%>%
  filter(!duplicated(sample_id))%>%
  ggplot(., aes(x=timepoint, y=log_qpcr, fill=class, shape=class))+
  geom_dotplot(aes(shape=class), binaxis="y", stackdir = "center")+
  scale_fill_manual(values=viridis::magma(n=3))+
  theme_minimal()

# temperature

clean_data %>%
  filter(class %in% c("A", "S"), timepoint %notin% c("bad_baseline", "day28", "day7"))%>%
  # filter(temperature>38)%>%
  select(-targetName)%>%
  filter(!duplicated(sample_id))%>%
  # filter(targetName %in% c("IL10", "CRP", "CXCL10", "LAG3", "IL27", "IFNG"))%>%
  arrange(temperature_cat)%>%
  ggplot(., aes(x=log_qpcr, y=temperature, color=temperature_cat))+
  geom_point()+
  # geom_smooth(method="lm")+
  # facet_wrap(~targetName)+
  scale_color_manual(values=c("black", "pink", "orange", "red"))+
  theme_minimal()

clean_data %>%
  filter(class %in% c("A", "S"), timepoint %notin% c("bad_baseline", "day28", "day7"))%>%
  filter(targetName %in% c("IL10", "CRP", "CXCL10", "IL6", "TNF", "IFNG"))%>%
  mutate(temperature_cat=if_else(is.na(temperature_cat), "<38", temperature_cat))%>%
  arrange(temperature_cat)%>%
  ggplot(., aes(x=log_qpcr, y=concentration, color=temperature_cat))+
  geom_point()+
  # geom_smooth(method="lm")+
  facet_wrap(~targetName)+
  scale_color_manual(values=c("black", "pink", "orange", "red"))+
  theme_minimal()


clean_data %>%
  filter(class %in% c("A", "S"), timepoint %notin% c("bad_baseline", "day28", "day7"))%>%
  # filter(targetName %in% c("IL10", "CRP", "CXCL10", "IL6", "GZMA", "IFNG"))%>%
  arrange(temperature_cat)%>%
  ggplot(., aes(x=temperature, y=concentration, color=temperature_cat))+
  geom_point()+
  # geom_smooth(method="lm")+
  facet_wrap(~targetName)+
  scale_color_manual(values=viridis::magma(n=6)[-2])+
  theme_minimal()





temp_purf <- clean_data %>%
  filter(class=="S", timepoint!="day28", timepoint!="bad_baseline")%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  group_by(targetName)%>%
  nest() %>%
  # mutate(model=map(data, ~lm(concentration~timepoint+ageyrs+id, data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~temperature+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.)))




pcr_redo <- haven::read_dta("~/Downloads/Revised qPCR results July 25 2024.dta")


hi <- clean_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"), infectiontype %in% c("A", "S"))%>%
  distinct(date, sample_id, timepoint, parasitedensity, qpcr, infectiontype)%>%
  mutate(cohortid = as.numeric(substr(sample_id, 1,3)),
         date = as.Date(date))

hi2 <- dplyr::left_join(hi, pcr_redo, by=c("cohortid", "date"))



clean_data%>%
  filter(timepoint=="day0", infectiontype%in% c("A", "S"))%>%
  ggplot(., aes(x=qpcr+0.1, y=parasitedensity+0.1))+
  geom_point(aes(color=infectiontype))+
  ggpubr::stat_cor()+
  geom_smooth(method="lm")+
  scale_color_manual(values=viridis::magma(3))+
  scale_y_log10(limits=c(0.1, 10^6))+
  scale_x_log10(limits=c(0.1, 10^6))+
  # facet_wrap(~infectiontype)+
  theme_minimal()+
  theme(legend.position = "none")





# sanchita


metadata <- readxl::read_excel("~/postdoc/stanford/plasma_analytes/MUSICAL/big_data/Immunology List _MUSICAL.xlsx")

pilot_data <- metadata %>%
  filter(id %in% c(268, 324, 137, 176, 353, 161, 363, 571))%>%
  # mutate(timepoint = factor(timepoint, levels=c("bad_baseline", "baseline", "day0", "day7", "day14", "day28")))%>%
  # filter(grepl("pilot", id))%>%
  distinct(id, timepoint_imm, infectiontype, qpcr, parasitedensity, fever, temperature)%>%
  mutate(qpcr=as.numeric(qpcr),
         timepoint_imm=factor(timepoint_imm), 
         temperature=as.numeric(temperature),
         parasitedensity=as.numeric(parasitedensity))


pilot_data%>%
  filter(infectiontype%in%c("A", "S"), timepoint_imm %in% c(-1, 0, 7, 14))%>%
  ggplot(., aes(x=timepoint_imm, y=qpcr+0.1, group=factor(id)))+
  geom_point()+
  geom_line()+
  scale_y_log10()+
  ggtitle("qPCR")+
  facet_wrap(~infectiontype)+
  theme_minimal()


pilot_data%>%
  filter(infectiontype%in%c("A", "S"), timepoint_imm %in% c(-1, 0, 7, 14))%>%
  ggplot(., aes(x=timepoint_imm, y=parasitedensity+0.1, group=factor(id)))+
  geom_point()+
  geom_line()+
  scale_y_log10()+
  ggtitle("blood smear")+
  facet_wrap(~infectiontype)+
  theme_minimal()


summary_table <- pilot_data%>%
  filter(infectiontype%in%c("A", "S"), timepoint_imm %in% c(-1, 0, 7, 14))%>%
  group_by(infectiontype, timepoint_imm)%>%
  summarise("mean_qpcr"=mean(qpcr),
            "sd_qpcr"=sd(qpcr),
            "mean_pardens"=mean(parasitedensity),
            "sd_pardens"=sd(parasitedensity),
            "n_fever"=sum(as.numeric(fever)))

write.csv(summary_table, "~/Downloads/summary_table_pilot.csv")