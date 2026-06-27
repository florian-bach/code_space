
jason_data <- read.csv("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/final_data_tables/t_cell_stim.csv")

jason_data <- jason_data %>%
  mutate(id=as.numeric(cohortid),
         "timepoint"=case_when(
           timepoint==-1~"baseline",
           timepoint==0~"day0",
           timepoint==7~"day7",
           timepoint==14~"day14",
           timepoint==28~"day28"))%>%
  filter(!is.na(cohortid))%


tr1_freqs_with_networks2 <- jason_data%>%
  mutate(infectiontype=case_when(infectiontype=="A2"~"A",
                                 infectiontype=="S2"~"S", .default = infectiontype))%>%
  mutate(sample_id=paste(id, timepoint, infectiontype, sep = " "))%>%
  filter(stim=="unstim")%>%
  select(sample_id, Tr1_Frequency)%>%
  inner_join(., musical_with_networks, by="sample_id")


(tr1_freqs_with_networks_plot <- tr1_freqs_with_networks2 %>%
  filter(infectiontype %in% c("A", "S"),
         timepoint%in%c("baseline", "day0", "day14"),
         network %in% c("MEtomato3", "MEturquoise3"))%>%
  mutate(network=case_match(network, "MEturquoise3"~"module 2", "MEtomato3"~"module 1", "MEthistle3"~"module 3", "MEwhite"~"module 4"))%>%
  mutate(day14_para=if_else(timepoint=="day14" & new_parasitedensity > 10 & infectiontype=="A", "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id, infectiontype)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), ">10 parasites / μL on day 14", "no parasites day14"))%>%
  mutate(infectiontype=if_else(infectiontype=="A", "asymptomatic", "malaria"))%>%
  ggplot(., aes(x=Tr1_Frequency/100, y=network_value, color=infectiontype))+
  geom_smooth(method="lm", fill="grey")+
  geom_point()+
  ggpubr::stat_cor(size=7)+
  ylab("module eigengene")+
  xlab("Tr1% of non-naive CD4")+
  scale_x_continuous(labels = scales::label_percent())+
  scale_color_manual(values = c("darkgrey", "darkred"))+
  facet_wrap(~network, scales = "free_x")+
    guides(color=guide_legend(override.aes = list(label="")))+
  theme_minimal(base_size = 32)+
    theme(legend.position = "bottom",
          strip.text = element_text(size=32)))
ggsave("~/postdoc/stanford/abstracts/immunology_retreat_2025/tr1_freqs_with_networks_plot.png", tr1_freqs_with_networks_plot, height=10, width=12, dpi=444)
# unname(unlist(infectiontype_cols))


target_cor_plot <- long_combo%>%
  filter(gate=="Tr1_Frequency", stim=="unstim", infectiontype%in%c("A", "S"), targetName %in%  c("CTLA4", "LAG3", "IFNG", "IL10"))%>%
  distinct(id, timepoint, infectiontype, targetName, gate, freq, stim, concentration)%>%
  ggplot(., aes(x=freq, y=concentration))+
  geom_smooth(method="lm", color="black")+
  geom_point(aes(fill=timepoint, shape=infectiontype), stroke=0.1)+
  ggpubr::stat_cor(method="spearman", size=9)+
  xlab("Tr1% of non-naive CD4")+
  scale_fill_manual(values=time_cols)+
  scale_shape_manual(values=c(21, 25))+
  facet_wrap(~targetName, scales = "free", ncol=2)+
  theme_minimal(base_size = 32)+
  guides(fill=guide_none())+
  theme(legend.position = "bottom")
ggsave("~/postdoc/stanford/abstracts/immunology_retreat_2025/target_cor_plot.png", target_cor_plot, height=12, width=10, dpi=444)



(tr1_freqs_ctrl <- tr1_freqs_with_networks2 %>%
    filter(infectiontype %in% c("A"),
           timepoint%in%c("baseline", "day0", "day14"),
           network %in% c("MEtomato3", "MEturquoise3"))%>%
    mutate(network=case_match(network, "MEturquoise3"~"module 2", "MEtomato3"~"module 1", "MEthistle3"~"module 3", "MEwhite"~"module 4"))%>%
    mutate(day14_para=if_else(timepoint=="day14" & new_qpcr > 10 & infectiontype=="A", "parasitemic_day14", "no_parasites_day14"))%>%
    group_by(id, infectiontype)%>%
    mutate(class2= if_else(any(day14_para=="parasitemic_day14"), ">10 parasites / μL on day 14", "no parasites day14"))%>%
    mutate(infectiontype=if_else(infectiontype=="A", "asymptomatic", "malaria"))%>%
    filter(!is.na(class2))%>%
    ggplot(., aes(x=timepoint, y=Tr1_Frequency/100, fill=class2))+
    geom_boxplot(outliers = F)+
    geom_point(aes(color=timepoint), position=position_dodge(width=0.75))+
    ggpubr::stat_compare_means(label = "p.signif", vjust=1)+
    ylab("Tr1% of non-naive CD4")+
    scale_y_continuous(labels = scales::label_percent())+
    scale_fill_manual(values=c("#636363", "darkgrey"))+
    scale_color_manual(values=unname(unlist(time_cols)))+
    guides(color=guide_legend(override.aes = list(label="")))+
    # theme_minimal(base_size = 32)+
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.title.x = element_blank()))
  
ggsave("~/postdoc/stanford/abstracts/immunology_retreat_2025/tr1_freqs_ctrl.png", height = 7, width = 14)


clearance_detail <- clean_data%>%
  filter(infectiontype %in% c("A"),
         timepoint%in%c("baseline", "day0", "day14"),
         targetName %in% c("BDNF", "CCL5", "TNFRSF14", "TNFSF13"))%>%
  mutate(day14_para=if_else(timepoint=="day14" & new_qpcr > 10 & infectiontype=="A", "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id, infectiontype)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), ">10 parasites / μL on day 14", "no parasites day14"))%>%
  filter(!is.na(class2))%>%
  ggplot(., aes(x=timepoint, y=concentration, fill=class2))+
  geom_smooth(method="lm", color="black")+
  geom_point(aes(color=timepoint), position = position_dodge(width=0.75))+
  geom_boxplot(aes(color=timepoint))+
  ggpubr::stat_compare_means(size=12, label="p.signif", hide.ns = T, vjust=1)+
  scale_fill_manual(values=c("#636363", "darkgrey"))+
  scale_color_manual(values=unname(unlist(time_cols)))+
  scale_shape_manual(values=c(21, 25))+
  facet_wrap(~targetName, scales = "free", ncol=2)+
  theme_minimal(base_size = 32)+
  guides(fill=guide_none())+
  theme(legend.position = "none", axis.title.x = element_blank())

ggsave("~/postdoc/stanford/abstracts/immunology_retreat_2025/clearance_detail_plot.png", clearance_detail, height=10, width=10, dpi=444)




micdrop_nulisa <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")

musical_nulisa <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")

slim_musical_nulisa <- musical_nulisa %>%
  filter(infectiontype%in%c("S", "A"))%>%
  mutate(mstatus=if_else(infectiontype=="S"&timepoint=="day0", 1, 0),
         conc=concentration)%>%
  select(id, ageyrs, timepoint, targetName, conc, log_qpcr, mstatus)%>%
  mutate(study="musical")

slim_micdrop <- micdrop_nulisa%>%
  mutate(ageyrs=ageinwks/52)%>%
  select(id, ageyrs, timepoint, targetName, conc, log_qpcr, mstatus)%>%
  mutate(study="micdrop")

combo_data <- bind_rows(slim_micdrop, slim_musical_nulisa)

combo_data%>%
  mutate(age_bracket=case_when(ageyrs<2~"< 2 years",
                               ageyrs>2&ageyrs<4~"2-4 years",
                               ageyrs>4&ageyrs<6~"4-6 years",
                               ageyrs>6~">6 years"))%>%
  filter(targetName%in%c("EPO", "TREM1"), !is.na(mstatus))%>%
  filter(!is.na(age_bracket))%>%
  ggplot(., aes(x=factor(age_bracket, levels = c("< 2 years", "2-4 years", "4-6 years", ">6 years")), y=conc, fill=factor(mstatus)))+
  geom_boxplot(outliers = F)+
  scale_fill_manual(values=c("black", "red"))+
  
  facet_wrap(~targetName, scale="free")+
  ggpubr::stat_compare_means(aes(group=mstatus), label = "p.signif")+
  xlab("")+
  theme_minimal(base_size = 20)+
  theme(legend.position = "none")


micdrop_nulisa_and_counts <- micdrop_nulisa%>%
  mutate(date=as.Date(date))%>%
  left_join(., blood_counts, by=c('id', 'date'))


epo_hb_plot <- micdrop_nulisa_and_counts%>%
  filter(targetName=="EPO", cell_type=="hb")%>%
  arrange(mstatus.x)%>%
  ggplot(., aes(x = conc, y=cell_freq, color=factor(mstatus.x)))+
  geom_smooth(method="lm")+
  geom_point()+
  ggpubr::stat_cor(method="spearman", hjust=-1.7, size=5)+
  scale_color_manual(values=c("black", "red"))+
  theme_minimal(base_size = 20)+
  xlab("EPO concentration")+
  ylab("hb concentration")+
  theme(legend.position = "none")

epo_qpcr_plot <- micdrop_nulisa_and_counts%>%
  filter(targetName=="EPO", cell_type=="hb")%>%
  arrange(mstatus.x)%>%
  ggplot(., aes(x = conc, y=10^log_qpcr, color=factor(mstatus.x)))+
  # geom_smooth(method="lm")+
  geom_point(aes(color=factor(mstatus.x)))+
  ggpubr::stat_cor(method="spearman", hjust=-1.7, size=5)+
  scale_color_manual(values=c("black", "red"))+
  theme_minimal(base_size = 20)+
  scale_y_log10()+
  xlab("EPO concentration")+
  ylab("malaria")+
  theme(legend.position = "none")

epo_hb_plot+epo_qpcr_plot

