library(patchwork)

p1 <- nulisa_and_seroconversion%>%
  filter(antigen %in% c("Rotavirus_IgA"), targetName %in% c("TLR3", "IL23"), timepoint=="8 weeks")%>%
  mutate(converts=ifelse(sero_conversion=="nothing", "no", "yes"))%>%
  ggplot(., aes(x=converts, y=conc, fill=factor(converts)))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.2, color="lightgrey")+
  xlab("seroconversion to Rotavirus IgA")+
  ylab("log2 concentration at 8 weeks [NPQ]")+
  ggpubr::stat_compare_means(method="wilcox.test", label = "p.format")+
  facet_wrap(~targetName, scales="free_y")+
  scale_y_continuous(expand=expansion(mult=0.1))+
  scale_fill_manual(values=c("#233875", "#FF5A00"))+
  theme_minimal(base_size = 14)+
  theme(legend.position = "none", 
        plot.title = element_text(hjust=0.5))


p2 <- nulisa_and_seroconversion%>%
  mutate(converts=ifelse(sero_conversion=="nothing", "no", "yes"))%>%
  filter(antigen %in% c("Pertussis"), targetName %in% c("AREG", "CSF1"), timepoint=="8 weeks")%>%
  ggplot(., aes(x=converts, y=conc, fill=factor(converts)))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.2, color="lightgrey")+
  ggpubr::stat_compare_means(method="wilcox.test", label = "p.format")+
  xlab("seroconversion to Pertussis IgG")+
  ylab("log2 concentration at 8 weeks [NPQ]")+
  facet_wrap(~targetName, scales="free_y")+
  scale_y_continuous(expand=expansion(mult=0.1))+
  scale_fill_manual(values=c("#233875", "#FF5A00"))+
  theme_minimal(base_size = 14)+
  theme(legend.position = "none", 
        plot.title = element_text(hjust=0.5))


# add seroconversion figure!! ####

line_data <- line_data%>%
  mutate(antigen=gsub("_", " ", antigen))

p3 <- wide_long_msd%>%
  mutate(antigen=gsub("_", " ", antigen))%>%
  filter(antigen%in%c("Pertussis", "Rotavirus IgA"))%>%
  arrange(desc(conversion))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=titer))+
  geom_segment(data = filter(line_data, next_color=="nothing", antigen%in%c("Pertussis", "Rotavirus IgA")), aes(x = timepoint, y = titer, xend = next_x, yend = next_y, color = next_color), alpha=0.7) +
  geom_segment(data = filter(line_data, next_color!="nothing", antigen%in%c("Pertussis", "Rotavirus IgA")), aes(x = timepoint, y = titer, xend = next_x, yend = next_y, color = next_color), alpha=0.7) +
  geom_point(aes(color=conversion))+
  scale_y_continuous(trans='log10')+
  scale_color_manual(values = c("#233875", "#B80F0A", "orange"))+
  theme_minimal(base_size = 14)+
  xlab("")+
  ylab("antibody concentration [MSD units]")+
  facet_wrap(~antigen, scales="free_y", nrow=2)+
  theme(
        legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(size=13)
        #axis.text.x = element_text(angle=30, vjust=1, hjust=1)
  )

nulisa <- (p2 / p1)+
  plot_layout(axis_titles = "collect")

p4 <- (p3 | nulisa)

ggsave("~/Library/CloudStorage/Box-Box/Jagannathan_Lab_Folder/PROJECTS/Florian_Grant_Applications/2026/Rosenkrantz/figures/seroconversion.png", p4, width=7, height=6, dpi=444, bg="white")
  
#sandbox ####

seroconversion_df3 <- wide_long_msd%>%
  distinct(id, conversion, antigen)%>%
  group_by(id, antigen)%>%
  mutate(sero_conversion=case_when(any(conversion=="converts 8 to 24")~"converts 8 to 24",
                                   any(conversion=="converts 24 to 52")~"converts 24 to 52",
                                   all(conversion=="nothing")~"nothing",
                                   .default = "somthing_else"))%>%
  distinct(id, antigen, sero_conversion)


nulisa_and_seroconversion <-  nulisa_data%>%
  inner_join(., seroconversion_df3, by=c("id"))


purf_test%>%
  filter(antigen %in% c(vaccines[4:5], vaccines_iga[4:5]))%>%
  filter(padj_8_24<0.1 | padj_24_52<0.1)%>%
  arrange(padj_24_52)


purf_test <- nulisa_and_seroconversion%>%
  mutate(sero_conversion=factor(sero_conversion, levels=c("nothing", "converts 8 to 24", "converts 24 to 52")))%>%
  group_by(targetName, antigen, timepoint)%>%
  nest()%>%
  mutate(model = purrr::map(data, ~lm(conc~sero_conversion, data=.)))%>%
  mutate(model_summary=purrr::map(model, ~summary(.)))%>%
  mutate(p_8_24=purrr::map(model_summary, ~coef(.)[11]))%>%
  mutate(p_24_52=purrr::map(model_summary, ~coef(.)[12]))%>%
  group_by(timepoint, targetName)%>%
  mutate(padj_8_24=p.adjust(p_8_24, method="fdr"))%>%
  mutate(padj_24_52=p.adjust(p_24_52, method="fdr"))

purf_test2 <- nulisa_and_seroconversion%>%
  filter(antigen %in% c(vaccines, vaccines_iga[4:5]))%>%
  mutate(sero_conversion=factor(sero_conversion, levels=c("nothing", "converts 8 to 24", "converts 24 to 52")))%>%
  mutate(converts=ifelse(sero_conversion=="nothing", "no", "yes"))%>%
  group_by(targetName, antigen, timepoint)%>%
  nest()%>%
  mutate(model = purrr::map(data, ~lm(conc~converts, data=.)))%>%
  mutate(model_summary=purrr::map(model, ~summary(.)))%>%
  mutate(p=purrr::map(model_summary, ~coef(.)[8]))%>%
  group_by(timepoint, targetName)%>%
  mutate(padj=p.adjust(p, method="fdr"))


nulisa_and_seroconversion%>%
  filter(antigen %in% c("Diptheria"), targetName %in% c("IFNA2", "IFNW1", "IFNA1; IFNA13", "IFNW1", "IFNG"))%>%
  ggplot(., aes(x=timepoint, y=conc, fill=factor(sero_conversion)))+
  geom_boxplot()+
  facet_wrap(~targetName+antigen)+
  theme_minimal()


nulisa_and_seroconversion%>%
  filter(antigen %in% c("Rotavirus_IgA"), targetName %in% c("CNTF", "CCL16"))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=conc, fill=factor(sero_conversion)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName+antigen, scales="free_y")+
  theme_minimal()



nulisa_and_seroconversion%>%
  filter(antigen %in% c("Rotavirus_IgA"), targetName %in% c("TLR3", "IL23"), timepoint=="8 weeks")%>%
  mutate(converts=ifelse(sero_conversion=="nothing", "no", "yes"))%>%
  ggplot(., aes(x=timepoint, y=conc, fill=factor(converts)))+
  geom_boxplot(outliers = F)+
  ggpubr::stat_compare_means(method="wilcox.test")+
  facet_wrap(~targetName+antigen, scales="free_y")+
  theme_minimal()+
  theme(legend.position = "none")

nulisa_and_seroconversion%>%
  mutate(converts=ifelse(sero_conversion=="nothing", "no", "yes"))%>%
  filter(antigen %in% c("Pertussis"), targetName %in% c("AREG", "CSF1"), timepoint=="8 weeks")%>%
  ggplot(., aes(x=timepoint, y=conc, fill=factor(converts)))+
  geom_boxplot(outliers = F)+
  ggpubr::stat_compare_means(method="wilcox.test")+
  facet_wrap(~antigen+targetName, scales="free_y")+
  theme_minimal()+
  theme(legend.position = "none", axis.title.x = element_blank())



nulisa_and_seroconversion%>%
  filter(antigen %in% c("Polio"), targetName %in% c("S100A9"))%>%
  mutate(converts=ifelse(sero_conversion=="nothing", "no", "yes"))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=conc, fill=factor(converts)))+
  geom_boxplot(outliers = F)+
  ggpubr::stat_compare_means(method="wilcox.test")+
  facet_wrap(~targetName+antigen, scales="free_y")+
  theme_minimal()


nulisa_and_seroconversion%>%
  filter(antigen %in% c(vaccines[4:5], vaccines_iga[4:5]), targetName %in% c("CCL16", "CNTF", "IL3RA"))%>%
  ggplot(., aes(x=timepoint, y=conc, fill=factor(sero_conversion)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName+antigen, scales="free_y")+
  theme_minimal()


nulisa_and_seroconversion%>%
  filter(antigen %in% c("Rotavirus", "Rotavirus_IgA"), targetName %in% c("CCL16", "CNTF", "CHI3L1", "IFNG"))%>%
  ggplot(., aes(x=timepoint, y=conc, fill=factor(sero_conversion)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName+antigen, scales="free_y")+
  theme_minimal()

nulisa_and_seroconversion%>%
  filter(antigen %in% c("Rotavirus", "Rotavirus_IgA"), targetName %in% c("IL2", "CNTF", "IL23", "IL19"))%>%
  ggplot(., aes(x=timepoint, y=conc, fill=factor(sero_conversion)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName+antigen, scales="free_y")+
  theme_minimal()

nulisa_and_seroconversion%>%
  filter(antigen %in% c("Rotavirus_IgA"), targetName %in% c("CNTF", "IL23", "CCL16"))%>%
  ggplot(., aes(x=timepoint, y=conc, fill=factor(sero_conversion)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName+antigen, scales="free_y")+
  theme_minimal()

nulisa_and_seroconversion%>%
  filter(antigen %in% c("Polio", "Polio_IgA"), targetName %in% c("TLR3", "S100A9", "LCN2"))%>%
  ggplot(., aes(x=timepoint, y=conc, fill=factor(sero_conversion)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName+antigen, scales="free_y")+
  theme_minimal()






# seropositivity
colnames(final_cutoff_frame) <- c("antigen", "cutoff_titer", "method" )

seropos_frame <- wide_long_msd%>%
  left_join(., final_cutoff_frame, by=c("antigen"))%>%
  mutate(sero_positive=ifelse(titer>cutoff_titer, 0, 1))%>%
  distinct(id, timepoint, antigen, sero_positive)



seropos_frame_plus_nulisa <- nulisa_data %>%
  left_join(., seropos_frame, by=c("id", "timepoint"))

try_purf_2 <- seropos_frame_plus_nulisa %>%
  filter(antigen %in% c(vaccines, vaccines_iga[4:5]))%>%
  group_by(targetName, antigen, timepoint)%>%
  nest()%>%
  mutate(model = purrr::map(data, ~glm(sero_positive~conc+log_qpcr+gender_categorical, data=., family = "binomial")))%>%
  mutate(model_summary=purrr::map(model, ~summary(.)))%>%
  mutate(pos_p=purrr::map(model_summary, ~coef(.)[14]))%>%
  group_by(timepoint, targetName)%>%
  mutate(pos_padj=p.adjust(pos_p, method="fdr"))

try_purf_2_sigs <- try_purf_2%>%
  filter(pos_padj<0.15)%>%
  arrange(pos_padj)

try_purf_2_sigs%>%
  filter(antigen %in% c(vaccines[4:5], vaccines_iga[4:5]))

seropos_frame_plus_nulisa%>%
  filter(targetName %in% c("CD200",  "IFNB1", "TGFB3"),
         antigen%in%c("Rotavirus_IgA"))%>%
  ggplot(., aes(x=timepoint, y=conc, fill=factor(sero_positive)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName+antigen, scales="free_y")+
  theme_minimal()

seropos_frame_plus_nulisa%>%
  filter(targetName %in% c("TNFSF9", "BMP7", "IL17A"),
         antigen%in%c("Tetanus", "Pertussis"))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=conc, fill=factor(sero_positive)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName+antigen, scales="free_y")+
  theme_minimal()

p1 <- seropos_frame_plus_nulisa%>%
  filter(targetName %in% c("CSF2RB", "VSNL1", "IL17A"),
         antigen%in%c("Tetanus"))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=conc, fill=factor(sero_positive)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName+antigen, scales="free_y")+
  theme_minimal()

seropos_frame_plus_nulisa%>%
  filter(targetName %in% c("TNFSF9"),
         antigen%in%c("Pertussis"))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=conc, fill=factor(sero_positive)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName+antigen, scales="free_y")+
  theme_minimal()

seropos_frame_plus_nulisa%>%
  filter(targetName %in% c("NAMPT"),
         antigen%in%c("Polio_IgA"))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=conc, fill=factor(sero_positive)))+
  geom_boxplot(outliers = F)+
  ggpubr::stat_compare_means(method = "wilcox.test")+
  facet_wrap(~targetName+antigen, scales="free_y")+
  theme_minimal()




seropos_frame_plus_nulisa%>%
  filter(targetName %in% c("TNFSF9"),
         antigen%in%c("Pertussis", "Tetanus", "Diptheria"))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=conc, fill=factor(sero_positive)))+
  geom_boxplot(outliers = F)+
  ggpubr::stat_compare_means(method = "wilcox.test")+
  facet_wrap(~targetName+antigen, scales="free_y")+
  theme_minimal()





`# p1 <- nulisa_and_seroconversion%>%
#   filter(antigen %in% c("Tetanus"), targetName %in% c("CSF2RB", "VSNL1", "IL17A"))%>%
#   ggplot(., aes(x=timepoint, y=conc, fill=factor(sero_conversion)))+
#   geom_boxplot(outliers = F)+
#   facet_wrap(~targetName+antigen, scales="free_y")+
#   theme_minimal()




