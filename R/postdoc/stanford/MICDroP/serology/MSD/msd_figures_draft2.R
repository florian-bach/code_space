library(tidyr)
library(dplyr)
library(xlsx)
library(ggplot2)
library(purrr)
library(data.table)
library(patchwork)
library(emmeans)



mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")
maternal_treatment_arms <- haven::read_dta("~/Library/CloudStorage/Box-Box/DP+SP study/Databases and preliminary findings/Final database used for analyses/DPSP treatment allocation_FINAL.dta")
epi_data <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/micdrop_metadata_for_immunology.csv")
nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")
msd_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/msd_results_v2.csv")
vaccine_coverage <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/long_vaccine_coverage_with_recall.csv")
IU_cutoffs <- readxl::read_excel("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/iu_conversion.xlsx")

diagnoses_of_kids_in_msd <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/diagnoses_of_kids_in_msd.csv")


vaccines = c("Diphtheria",     "Measles" ,         "Pertussis",     "Polio",
             "Rotavirus" ,    "Rubella",       "Tetanus", "Pneumo.1.4.14")
vaccines_iga <- c("Diphtheria_IgA", "Measles_IgA",  "Pertussis_IgA", "Polio_IgA",  "Rotavirus_IgA", "Rubella_IgA",  "Tetanus_IgA",  "Varicella_IgA")
viruses <- c("Mumps",    
             "Varicella",
             "CoV.2",    
             "EV.71",    
             "EV.D68",   
             "Flu.Mix",  
             'hMPV',     
             "PIV.Mix",  
             "RSV",      
             "RV.C",     
             "Flu.A.H1", 
             "Flu.A.H3", 
             "Flu.B",    
             "PIV.1",    
             "PIV.2",    
             'PIV.3',    
             "PIV.4")
cmv_data <-  read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/CMV_ELISA_Uganda.csv", skip = 1)


gravidity_cols <- c("#1E88E5", "#A798EC", "#582EE0")
names(gravidity_cols) <- c("Secundigravid", "Primigravid", "Multigravid")

time_cols <- viridis::magma(n=5)[2:4]
# seroconversion_cols <- c("red", "orange", "black")
seroconversion_cols <- c("#0061c3", "#acff24", "#36454f")
names(seroconversion_cols) <- c("converts 24 to 52", "converts 8 to 24", "nothing")

long_msd <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/long_msd_and_epi.csv")

vaccine_coverage <- readxl::read_excel("~/postdoc/stanford/clinical_data/MICDROP/vaccines/Vaccination data_20251215.xlsx")

schedule <- data.frame(vaccine_type=c("BCG at birth",
                                      "DTP_Hib",
                                      "Oral Polio at birth",
                                      "Measles Rubella",
                                      "PCV10",
                                      "Rotarix",
                                      "Oral Polio",
                                      "Sabin Polio"),
                       target_doses=c(1, 3, 1, 2, 3, 2, 3, 1))

long_vaccine_coverage <- vaccine_coverage%>%
  mutate(rota2=ifelse(rota2=="NULL", 0, rota2))%>%
  mutate(rota2=as.numeric(rota2))%>%
  select(-inputdate, -lastmod, -clean)%>%
  pivot_longer(cols = colnames(vaccine_coverage)[4:19], names_to = "vaccine_dose", values_to="dose_received")%>%
  mutate(vaccine_type = case_when(grepl("polio*", vaccine_dose)&sched==0~"Oral Polio at birth",
                                  grepl("polio*", vaccine_dose)~"Oral Polio",
                                  grepl("Bcg*", vaccine_dose)~"BCG at birth",
                                  grepl("pcv1*", vaccine_dose, ignore.case = T)~"PCV10",
                                  grepl("dPT*", vaccine_dose)~"DTP_Hib",
                                  grepl("rota*", vaccine_dose)~"Rotarix",
                                  grepl("measles*", vaccine_dose)~"Measles Rubella",
                                  grepl("iPV*", vaccine_dose)~"Sabin Polio"))%>%
  group_by(vaccine_type, id)%>%
  summarise("number_of_doses_received"=sum(dose_received, na.rm = T))%>%
  left_join(., schedule)%>%
  mutate(course_outcome=case_when(number_of_doses_received>=target_doses~"complete",
                                  number_of_doses_received<target_doses&number_of_doses_received>0~"partial",
                                  number_of_doses_received==0~"none",
                                  is.na(number_of_doses_received)~"no record"
  ))
cmv_data <-  read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/CMV_ELISA_Uganda.csv", skip = 1)

# write.csv(long_vaccine_coverage, "~/postdoc/stanford/plasma_analytes/MICDROP/MSD/cleaned_vaccine_coverage.csv", row.names = F)  

long_msd <- msd_data%>%
  mutate(sample=paste(SubjectID, "_", "tp", TimePt, sep=""))%>%
  mutate(id=SubjectID, timepoint=paste(TimePt, "weeks"))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))%>%
  select(-SubjectID, -TimePt)%>%
  pivot_longer(cols=-c(sample, id, timepoint), names_to = "antigen", values_to = "titer")%>%
  mutate(antigen=ifelse(antigen=="Diptheria", "Diphtheria", antigen))%>%
  mutate(Ig_class=ifelse(grepl("*IgA$", antigen), "IgA", "IgG"))%>%
  filter(!is.na(titer))%>%
  mutate(CMV_at_52_weeks=cmv_data$Diagnosis[match(id, cmv_data$id)],)%>%
  mutate(vaccine_type=case_match(antigen,
                                 "Diphtheria"~"DTP_Hib",
                                 "Tetanus"~"DTP_Hib",
                                 "Pertussis"~"DTP_Hib",
                                 "Polio_IgA"~"Oral Polio",
                                 "Polio"~"Sabin Polio",
                                 "Measles"~"Measles Rubella",
                                 "Rubella"~"Measles Rubella",
                                 "Rotavirus"~"Rotarix",
                                 "Pneumo.1.4.14"~"PCV10",
                                 "Rotavirus_IgA"~"Rotarix",.default = "none")
  )%>%
  left_join(., long_vaccine_coverage, by=c("id", "vaccine_type"))%>%
  mutate(course_outcome=ifelse(!id%in%vaccine_coverage$id, "no record", course_outcome))


# additional data from the moment of sampling
additional_clinical_data <- nulisa_data%>%
  distinct(id, sample, date, ageinwks, gender_categorical, mstatus, qPCRparsdens, pardens)

antibodies_and_epi <- epi_data%>%
  inner_join(., long_msd, by="id")%>%
  left_join(additional_clinical_data, by="sample")






# Figure 1 longitudinal data overview vaccines ####

vax_overview <- long_msd %>%
  filter(!is.na(timepoint), antigen %in% c(vaccines, vaccines_iga[c(3:5)]))%>%
  ggplot(., aes(timepoint, y=titer, fill=timepoint))+
  geom_violin()+
  geom_boxplot(outliers = F, width=0.2, color="white")+
  # geom_point(shape=21, color="lightgrey")+
  scale_y_log10()+
  ylab("antibody concentration (MSD units)")+
  #ggpubr::stat_compare_means(vjust = 0.5, comparisons = list(c("0", "1")))+
  facet_wrap(~antigen, scales="free_y", ncol=4)+
  scale_fill_manual(values=time_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=8, angle = 90, vjust=0.5),
        legend.position = "none")
        #axis.text.x = element_text(angle=45, vjust=0.5))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/vaccine_time_overview.png", vax_overview, width=8, height=6, dpi=444, bg="white")


(vax_protection <-  long_msd %>%
  filter(!is.na(timepoint), antigen %in% c(IU_cutoffs$antigen[-c(3,7)]))%>%
  left_join(IU_cutoffs, by="antigen")%>%
  mutate(IU=titer*IU_conversion)%>%
  ggplot(., aes(timepoint, y=IU, fill=timepoint))+
  geom_hline(data=IU_cutoffs[-c(3,7),], aes(yintercept = as.numeric(sigal_cutoff)))+
  geom_boxplot(outliers = F, width=0.2, color="black")+
  geom_point(shape=21, color="black", position=position_jitter(width=0.1), alpha=0.35)+
  # scale_y_log10()+
  ylab("antibody concentration (IU / mL)")+
  #ggpubr::stat_compare_means(vjust = 0.5, comparisons = list(c("0", "1")))+
  facet_wrap(~antigen, scales="free_y", nrow=2)+
  scale_fill_manual(values=time_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=8, angle = 90, vjust=0.5),
        legend.position = "none"))
# ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/vax_protection.png", vax_protection, width=8, height=2.5, dpi=444, bg="white")


# FIgure 2 gravidity

gravid_plot <- antibodies_and_epi%>%
  mutate(gravid_cat = case_match(
    gravidcat,
    1 ~ "Primigravid",
    2 ~ "Secundigravid",
    3 ~ "Multigravid"
  )) %>%
  filter(antigen %in% mom_gravidity_purf_sigs$antigen)%>%
  filter(timepoint=="8 weeks")%>%
  ggplot(., aes(x = timepoint, y = titer, fill=factor(gravid_cat, levels=c("Primigravid", "Secundigravid", "Multigravid")))) +
  geom_boxplot(outliers = FALSE, position = position_dodge(0.75)) +
  scale_y_log10() +
  ylab("IgG concentration (MSD units)")+
  facet_wrap(~antigen, ncol=6, scales="free_y") +
  theme_minimal() +
  scale_fill_manual(values=gravidity_cols)+
  theme(axis.title.x = element_blank(),
        legend.title =  element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/gravidity_plot.png", gravid_plot, width=8, height=2.4)
# Figure 2 protection summary at 1 year ####
protection_barchart <- long_msd %>%
  filter(timepoint=="52 weeks", antigen %in% c(IU_cutoffs$antigen[-c(3,7)]))%>%
  left_join(IU_cutoffs, by="antigen")%>%
  mutate(IU=titer*IU_conversion)%>%
  mutate(above_sigal=ifelse(IU>=sigal_cutoff, "protected", "not protected"))%>%
  mutate(above_plotkin=ifelse(IU>=plotkin_cutoff, "protected", "not protected"))%>%
  distinct(id, timepoint, antigen, above_sigal)%>%
  group_by(antigen, timepoint, above_sigal)%>%
  count()%>%
  # pivot_wider(names_from = above_sigal, values_from = n)%>%
  # print(n=21)%>%
  ggplot(., aes(x=antigen, y=n, fill=above_sigal))+
  geom_bar(stat="identity", position="stack")+
  scale_fill_manual(values=c("darkred", "darkgrey"))+
  ylab("% Protected at 52 weeks")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=8, angle = 90, vjust=0.5),
        legend.title = element_blank())

# ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/protection_barchart.png", protection_barchart, width=4, height=4, dpi=444, bg="white")

protection_figure <- vax_protection + protection_barchart + plot_layout(widths = c(3,1)) + plot_annotation(tag_levels = "A")

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/protection_figure.png", protection_figure, width=8, height=4, dpi=444, bg="white")

## protection table ####
long_msd %>%
  filter(timepoint=="52 weeks", antigen %in% c(IU_cutoffs$antigen[-c(3,7)]))%>%
  left_join(IU_cutoffs, by="antigen")%>%
  mutate(IU=titer*IU_conversion)%>%
  mutate(above_sigal=ifelse(IU>=sigal_cutoff, "protected", "not protected"))%>%
  distinct(id, timepoint, antigen, above_sigal)%>%
  group_by(antigen, timepoint, above_sigal)%>%
  count()%>%
  pivot_wider(names_from = above_sigal, values_from = n)%>%
  print(n=21)%>%
  mutate(protected_perc=protected/n_distinct(long_msd$id))



# Figure 3 seroconversions to viral diseases ####

wide_long_msd <-  antibodies_and_epi %>%
  filter(!is.na(timepoint))%>%
  pivot_wider(names_from = "timepoint", values_from = "titer", id_cols = c("id", "antigen", "vaccine_type", "course_outcome", "number_of_doses_received", "treatmentarm"))%>%
  mutate(log2fc_8_24=log2(`24 weeks`/`8 weeks`),
         log2fc_24_52=log2(`52 weeks`/`24 weeks`))%>%
  group_by(antigen)%>%
  mutate("z_fc_8_24"=scale(log2fc_8_24, center = T))%>%
  pivot_longer(cols = c("8 weeks", "24 weeks","52 weeks"), names_to = "timepoint", values_to = "titer")%>%
  pivot_longer(cols = c("log2fc_8_24", "log2fc_24_52"), names_to = "fc_flavour", values_to = "log2fc")%>%
  mutate(conversion=case_when(timepoint=="52 weeks" & fc_flavour=="log2fc_24_52" & log2fc > 2 & titer > 10 ~ "converts 24 to 52",
                              timepoint=="24 weeks" & fc_flavour=="log2fc_8_24" & z_fc_8_24 > 2 & titer > 10 | 
                                timepoint=="24 weeks" & fc_flavour=="log2fc_8_24" & log2fc > 1 & titer > 10  ~ "converts 8 to 24",
                              .default="nothing"))%>%
  mutate(mom_rx=maternal_treatment_arms$treatmentarm[match(id.x-10000, maternal_treatment_arms$id)],
         mom_rx=case_match(mom_rx,
                           1~"SP",
                           2~"DP",
                           3~"DPSP"))


line_data <- wide_long_msd %>%
  distinct(id, timepoint, antigen, titer, conversion)%>%
  group_by(id, antigen, timepoint)%>%
  filter(if (n()>1) conversion!="nothing" else TRUE)%>%
  ungroup()%>%
  mutate(conversion=factor(conversion, levels=c("nothing", "converts 8 to 24", "converts 24 to 52")))%>%
  ungroup()%>%
  arrange(id, factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks"))) %>%
  group_by(id, antigen) %>%
  mutate(next_x = lead(timepoint),
         next_y = lead(titer),
         next_color = lead(conversion))%>%
  filter(!is.na(next_x))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))

select_viruses <- c("Mumps", "RSV", "CoV.2", "EV.71","Flu.Mix", "RV.C")

select_virus_seroconversion_plot <- wide_long_msd%>%
  filter(antigen%in%select_viruses)%>%
  arrange(desc(conversion))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=titer))+
  geom_segment(data = filter(line_data, next_color=="nothing", antigen%in%select_viruses), aes(x = timepoint, y = titer, xend = next_x, yend = next_y, color = next_color)) +
  geom_segment(data = filter(line_data, next_color!="nothing", antigen%in%select_viruses), aes(x = timepoint, y = titer, xend = next_x, yend = next_y, color = next_color)) +
  geom_point(aes(color=conversion))+
  scale_y_continuous(trans='log10')+
  scale_color_manual(values = seroconversion_cols)+
  theme_minimal()+
  facet_wrap(~antigen, scales="free_y", ncol=3)+
  theme(axis.title = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(size=8, angle = 90, vjust=0.5)
        #axis.text.x = element_text(angle=30, vjust=1, hjust=1)
        )

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/select_virus_seroconversion_plot.png", virus_seroconversion_plot, width=8, height=8, dpi=444, bg="white")



vaccine_seroconversion_plot <- wide_long_msd%>%
  filter(antigen%in%c(vaccines, vaccines_iga[4:5]))%>%
  arrange(desc(conversion))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=titer))+
  geom_segment(data = filter(line_data, next_color=="nothing", antigen%in%c(vaccines, vaccines_iga[4:5])), aes(x = timepoint, y = titer, xend = next_x, yend = next_y, color = next_color)) +
  geom_segment(data = filter(line_data, next_color!="nothing", antigen%in%c(vaccines, vaccines_iga[4:5])), aes(x = timepoint, y = titer, xend = next_x, yend = next_y, color = next_color)) +
  geom_point(aes(color=conversion))+
  scale_y_continuous(trans='log10')+
  scale_color_manual(values = seroconversion_cols)+
  theme_minimal()+
  facet_wrap(~antigen, scales="free_y", nrow=6)+
  theme(axis.title = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(size=8, angle = 90, vjust=0.5)
        #axis.text.x = element_text(angle=30, vjust=1, hjust=1)
  )

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/vaccine_seroconversion_plot.png", vaccine_seroconversion_plot, width=8, height=11, dpi=444, bg="white")

seroconversion_table <- wide_long_msd%>%
  distinct(id, antigen, conversion)%>%
  group_by(antigen, conversion)%>%
  filter(conversion!="nothing", antigen %in% viruses)%>%
  count()%>%
  arrange()%>%
  group_by(antigen)%>%
  mutate(conversion_perc=n/n_distinct(wide_long_msd$id))%>%
  arrange(desc(conversion_perc))


extra_rows <- seroconversion_table %>%
  group_by(antigen)%>%
  summarise(conversion_perc=sum(conversion_perc))%>%
  mutate(
    conversion = "nothing",
    conversion_perc = 1 - conversion_perc,
    n = 89*conversion_perc,
    
  )%>%
  arrange(conversion_perc)

seroconversion_table2 <- seroconversion_table%>%
  bind_rows(., extra_rows)%>%
  arrange(antigen)

seroconvversion_barchart <- seroconversion_table2%>%
  ggplot(., aes(x=factor(antigen, levels=extra_rows$antigen), y=conversion_perc, fill=conversion))+
  geom_bar(stat="identity", position="stack")+
  scale_fill_manual(values = seroconversion_cols)+
  scale_y_continuous(labels = scales::label_percent())+
  ylab("Seroconversion")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=8, angle = 90, vjust=0.5),
        legend.position = "none")

viral_seroconversion_figure <- select_virus_seroconversion_plot+seroconvversion_barchart + plot_layout(widths = c(2.7,1)) + plot_annotation(tag_levels = "A")
ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/viral_seroconversion_figure.png", viral_seroconversion_figure, width=8, height=4, dpi=444, bg="white")


# Figure 4 antibody concentrations to both vaccines and viruses correlate with parasitemia ####

antigen_parasitemia_corr_plot <- antibodies_and_epi %>%
  filter(antigen %in% c("PIV.1" ,    "PIV.2",      "Varicella", "Measles", "PIV.4", "Mumps",  "Rotavirus", "Rotavirus_IgA", "RV.C", "Diphtheria"),
         timepoint!="8 weeks")%>%
  ggplot(., aes(x=qPCRparsdens+0.001 , y=titer+0.001, color=timepoint))+
  geom_smooth(method="lm", )+
  geom_point(alpha=0.68)+
  facet_wrap(~antigen, nrow=2)+
  scale_x_log10()+
  scale_y_log10()+
  xlab("qPCR parasitemia")+
  ylab("antibody concentration (MSD units)")+
  ggpubr::stat_cor(aes(label = ..r.label..),method="spearman",
                   size=3.5, label.y=c(4.9, 5.25),label.x.npc=0.1, r.accuracy = 0.01, p.accuracy = 0.01)+
  scale_color_manual(values = time_cols[2:3])+
  theme_minimal(base_size = 15)+
  theme(axis.text.x = element_text(angle=90, vjust=0.5),
        legend.position = "bottom")


ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/fig4_combo_plot.png", antigen_parasitemia_corr_plot, width=12, height = 6, dpi=444, bg="white")


# Figure 5 vaccine titers correlated with plasmodium transmission ####

incidence_plot2 <- antibodies_and_epi%>%
  filter(antigen%in%c("Diphtheria", "Varicella", "Diptheria_IgA"), timepoint=="52 weeks")%>%
  ggplot(., aes(x=factor(total_n_para_12), y=titer))+
  geom_point(aes(color=factor(total_n_para_12)))+
  # geom_smooth(aes(x=incidence_value, y=titer), method="lm")+
  geom_boxplot(aes(fill=factor(total_n_para_12)), outliers = F)+
  stat_summary(geom="crossbar", colour="white", fun = "median", width=0.65, linewidth=0.19 )+
  scale_y_log10()+
  # xlab(inci_type)+
  xlab("number of parasitemic months in the first year of life")+
  facet_wrap(~antigen, nrow=1, scales="free")+
  # ggpubr::stat_compare_means()+
  theme_minimal(base_size = 16)+
  scale_fill_manual(values = viridis::viridis(n=13))+
  scale_color_manual(values = viridis::viridis(n=13))+
  theme(legend.position = "none",
        axis.title.y = element_blank())

chemoprevention_vaccine <- antibodies_and_epi%>%
    filter(timepoint=="52 weeks", antigen %in% c("Diphtheria", "Diptheria_IgA", "Varicella", "Measles_IgA", "Rubella_IgA"))%>%
    ggplot(., aes(x=factor(treatmentarm), y=titer, fill=treatmentarm))+
    geom_boxplot(outliers = F)+
    geom_point(shape=21)+
    scale_y_log10()+
    xlab("IPTc treatmentarm")+
  ggpubr::stat_compare_means(label = "p.signif", vjust=1, label.x = 1.44, size=4)+
    facet_wrap(~factor(antigen, levels=c("Diphtheria", "Diptheria_IgA", "Varicella", "Measles_IgA", "Rubella_IgA")), scales="free_y")+
    scale_fill_manual(values=c("darkred", "#283593"))+
    theme_minimal(base_size = 15)+
  theme(legend.position = "none",
        axis.title.y = element_blank())

incidence_plots <- incidence_plot2 / chemoprevention_vaccine +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(
    tag_levels = "A"
    # left = "antibody concentration at 52 weeks [MSD units]"
  )


ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/incidence_plots.png", incidence_plots, width=9, height=8, dpi=444, bg="white")



antibodies_and_epi%>%
  filter(antigen%in%c("Rotavirus"), timepoint=="52 weeks")%>%
  ggplot(., aes(x=factor(total_n_malaria_6), y=titer))+
  geom_point()+
  # geom_smooth(aes(x=incidence_value, y=titer), method="lm")+
  geom_boxplot(aes(fill=factor(total_n_malaria_6)), outliers = F)+
  stat_summary(geom="crossbar", colour="white", fun = "median", width=0.65, linewidth=0.19 )+
  scale_y_log10()+
  ylab("titer at 52 weeks")+
  xlab("total_n_malaria_6")+
  facet_wrap(~antigen)+
  ggpubr::stat_compare_means(vjust = 0.5, comparisons = list(c("0", "1")))+
  theme_minimal()+
  scale_fill_manual(values = viridis::viridis(n=13))+
  theme(legend.position = "none")


antibodies_and_epi%>%
  filter(antigen%in%c("Varicella"))%>%
  ggplot(., aes(x=factor(treatmentarm), y=titer))+
  geom_point()+
  # geom_smooth(aes(x=incidence_value, y=titer), method="lm")+
  geom_boxplot(aes(fill=factor(treatmentarm)), outliers = F)+
  stat_summary(geom="crossbar", colour="white", fun = "median", width=0.65, linewidth=0.19 )+
  scale_y_log10()+
  xlab(inci_type)+
  ylab("antibody concentration [MSD units]")+
  xlab("treatment arm")+
  facet_wrap(~antigen+timepoint)+
  ggpubr::stat_compare_means(comparisons = list(c("DP 1 year", "Placebo")))+
  theme_minimal()+
  scale_fill_manual(values = viridis::viridis(n=2))+
  theme(legend.position = "none")


antibodies_and_epi%>%
  filter(antigen%in%c("Varicella"))%>%
  ggplot(., aes(x=qPCRparsdens+0.001, y=titer, color=factor(treatmentarm)))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_y_log10()+
  scale_x_log10()+
  ylab("antibody concentration [MSD units]")+
  xlab("treatment arm")+
  facet_wrap(~timepoint)+
  ggpubr::stat_cor(method = "spearman")+
  theme_minimal()+
  scale_color_manual(values = viridis::viridis(n=3))+
  theme(legend.position = "none")
# Figure 5 viral infection is associated with vaccine titers ####

df_pred <- wide_long_msd%>%
  distinct(id, antigen, timepoint, titer) %>%
  filter(timepoint == "52 weeks", antigen %in% c(vaccines, vaccines_iga[c(4,5)])) %>%
  rename(
    vaccine_antigen = antigen,
    vaccine_titer = titer
  )

df_target <- wide_long_msd%>%
  filter(timepoint == "24 weeks" & fc_flavour=="log2fc_8_24" |
           timepoint == "52 weeks" & fc_flavour=="log2fc_24_52",
         antigen %in% viruses) %>%
  distinct(id, timepoint, antigen, conversion) %>%
  pivot_wider(values_from = conversion, names_from = timepoint)%>%
  mutate(sero_conversion=case_when(`24 weeks`=="nothing"&`52 weeks`=="nothing"~"nothing",
                                   `24 weeks`=="nothing"&`52 weeks`=="converts 24 to 52"~"converts 24 to 52",
                                   `24 weeks`=="converts 8 to 24"&`52 weeks`=="nothing"~"converts 8 to 24", .default = ))%>%
  distinct(id, antigen, sero_conversion)%>%
  rename(
    virus_antigen = antigen,
    virus_seroconversion = sero_conversion
  )%>%
  filter(!is.na(virus_seroconversion))



wide_long_msd_cross <- df_pred %>%
  inner_join(df_target, by = "id") %>%
  mutate(log2_titer = log2(vaccine_titer))

purf3 <- wide_long_msd_cross%>%
  filter(virus_antigen %notin%c("PIV.1", "PIV.2", "Flu.A.H3"),
         vaccine_antigen!="Pneumo.1.4.14",
         !(virus_antigen=="Mumps"&vaccine_antigen=="Diphtheria"))%>%
  filter()
  pivot_wider(names_from = virus_seroconversion, values_from = vaccine_titer, id_cols = c("id", "vaccine_antigen", "virus_antigen"))%>%
  group_by(vaccine_antigen, virus_antigen)%>%
  nest()%>%
  mutate(wilcox1p=map_dbl(data, ~wilcox.test(x = .$nothing, .$`converts 8 to 24`)$p.value))%>%
  mutate(wilcox2p=map_dbl(data, ~wilcox.test(x = .$nothing, .$`converts 24 to 52`)$p.value))%>%
  mutate(wilcox3p=map_dbl(data, ~wilcox.test(x = .$`converts 8 to 24`, .$`converts 24 to 52`)$p.value))%>%
  ungroup()%>%
  mutate(wilcox1padj=map_dbl(wilcox1p, ~p.adjust(., method="fdr")))%>%
  mutate(wilcox2padj=map_dbl(wilcox2p, ~p.adjust(., method="fdr")))%>%
  mutate(wilcox3padj=map_dbl(wilcox3p, ~p.adjust(., method="fdr")))


virus_conversion_vaccine_sigs <- purf3%>%
  filter(wilcox1padj < 0.05 | wilcox2padj < 0.05 | wilcox3padj < 0.05)



plot_list <- vector("list", length = nrow(virus_conversion_vaccine_sigs))

for(i in 1:nrow(virus_conversion_vaccine_sigs)){
  
  virus = virus_conversion_vaccine_sigs$virus_antigen[i]
  vaccine = virus_conversion_vaccine_sigs$vaccine_antigen[i]
  
  p = wide_long_msd_cross%>%
    mutate(virus_seroconversion=factor(virus_seroconversion, levels=c("converts 8 to 24",
                                                                      "converts 24 to 52",
                                                                      "nothing")))%>%
    filter(virus_antigen==virus, vaccine_antigen==vaccine)%>%
    ggplot(., aes(x=virus_seroconversion, y=vaccine_titer, fill=virus_seroconversion))+
    # geom_violin()+
    geom_boxplot(outliers = F, width = 0.3)+
    geom_point()+
    ggpubr::stat_compare_means(size=2, comparisons = list(
      c("converts 8 to 24", "converts 24 to 52"),
      c("nothing", "converts 24 to 52"),
      c("nothing", "converts 8 to 24")))+
    facet_wrap(~vaccine_antigen+virus_antigen, nrow=4, scales="free_y")+
    scale_y_log10(expand = expansion(mult = 0.1))+
    scale_x_discrete(labels = scales::label_wrap(width = 10))+
    scale_fill_manual(values = seroconversion_cols)+
    theme_minimal()+
    theme(axis.title.x = element_blank(),
          legend.position = "none",
          axis.title.y = element_blank())
  
  plot_list[[i]] = p
  
}

multiplot <- cowplot::plot_grid(plotlist = plot_list, ncol = 3)
ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/virus_conversion_vaccine_plot.png", multiplot, width=7, height=9, dpi=444, bg="white")



# S1 8 week titers correlate with 24 week titers but not 52 week titers ####

# correlate 8 weeks with 52 weeks
eight_52_corr_plot <- long_msd %>%
  filter(!is.na(timepoint), antigen %in% c(vaccines, vaccines_iga[c(4,5)]))%>%
  pivot_wider(names_from = "timepoint", values_from = "titer", id_cols = c("id", "antigen", "Ig_class"))%>%
  ggplot(., aes(x=`8 weeks`, y=`52 weeks`))+
  geom_point()+
  facet_wrap(~antigen, nrow=2)+
  scale_x_log10()+
  scale_y_log10()+
  xlab("antibody concentration at 8 weeks (MSD units)")+
  ylab("antibody concentration at 52 weeks (MSD units)")+
  ggpubr::stat_cor(method="spearman", color="red", size=3)+
  geom_smooth(method="lm")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/eight_52_corr_plot.png", eight_52_corr_plot, width=8, height=4, dpi=444, bg="white")

eight_24_corr_plot <- long_msd %>%
  filter(!is.na(timepoint), antigen %in% c(vaccines, vaccines_iga[c(4,5)]))%>%
  pivot_wider(names_from = "timepoint", values_from = "titer", id_cols = c("id", "antigen", "Ig_class"))%>%
  ggplot(., aes(x=`8 weeks`, y=`24 weeks`))+
  geom_point()+
  facet_wrap(~antigen, nrow=2)+
  scale_x_log10()+
  scale_y_log10()+
  xlab("antibody concentration at 8 weeks (MSD units)")+
  ylab("antibody concentration at 24 weeks (MSD units)")+
  ggpubr::stat_cor(method="spearman", color="red", size=3, vjust=0.5)+
  geom_smooth(method="lm")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/eight_24_corr_plot.png", eight_24_corr_plot, width=8, height=4, dpi=444, bg="white")



# S2 IgG~IgA correlation ####

# correlate 8 weeks with 52 weeks
iga_igg_corr_plot <- long_msd %>%
  mutate(base_antigen=ifelse(Ig_class=="IgA", substr(antigen, 1, nchar(antigen)-4), antigen))%>%
  filter(base_antigen %in% vaccines & base_antigen!="Pneumo.1.4.14")%>%
  pivot_wider(names_from = "Ig_class", values_from = "titer", id_cols = c("id", "base_antigen", "timepoint"))%>%
  ggplot(., aes(x=IgG, y=IgA, color=timepoint))+
  geom_point()+
  facet_wrap(~base_antigen, nrow=2)+
  scale_x_log10()+
  scale_y_log10()+
  xlab("IgG (MSD units)")+
  ylab("IgA (MSD units)")+
  ggpubr::stat_cor(method="spearman", size=3)+
  geom_smooth(method="lm")+
  scale_color_manual(values = time_cols)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/iga_igg_corr_plot.png", iga_igg_corr_plot, width=8, height=4, dpi=444, bg="white")

# S3 correlation between parasitemia and anitbodies


qpcr_vaccine_antibody_corr <- antibodies_and_epi %>%
  filter(antigen %in% c(vaccines, vaccines_iga[c(4,5)]), antigen!="Pneumo.1.4.14")%>%
  ggplot(., aes(x=qPCRparsdens+0.001 , y=titer, color=timepoint))+
  geom_point()+
  facet_wrap(~antigen, nrow=2, scale="free")+
  scale_x_log10()+
  scale_y_log10()+
  xlab("qPCR parasitemia")+
  ylab("antibody concentration (MSD units)")+
  ggpubr::stat_cor(method="spearman", size=2.5, r.accuracy = 0.01, p.accuracy = 0.01)+
  geom_smooth(method="lm")+
  scale_color_manual(values = time_cols)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/qpcr_vaccine_antibody_corr.png", qpcr_vaccine_antibody_corr, width=8, height=6, dpi=444, bg="white")



qpcr_virus_antibody_corr <- antibodies_and_epi %>%
  filter(antigen %in% c(viruses))%>%
  ggplot(., aes(x=qPCRparsdens+0.001 , y=titer, color=timepoint))+
  geom_point()+
  facet_wrap(~antigen, nrow=3, scales="free")+
  scale_x_log10()+
  scale_y_log10()+
  xlab("qPCR parasitemia")+
  ylab("antibody concentration (MSD units)")+
  ggpubr::stat_cor(method="spearman", size=2.5, r.accuracy = 0.01, p.accuracy = 0.01)+
  geom_smooth(method="lm")+
  scale_color_manual(values = time_cols)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/qpcr_virus_antibody_corr.png", qpcr_virus_antibody_corr, width=12, height=8, dpi=444, bg="white")




# S3 full viral seroconversion plot ####


full_virus_seroconversion_plot <- wide_long_msd%>%
  filter(antigen%in%viruses)%>%
  arrange(desc(conversion))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=titer))+
  geom_segment(data = filter(line_data, next_color=="nothing", antigen%in%viruses), aes(x = timepoint, y = titer, xend = next_x, yend = next_y, color = next_color)) +
  geom_segment(data = filter(line_data, next_color!="nothing", antigen%in%viruses), aes(x = timepoint, y = titer, xend = next_x, yend = next_y, color = next_color)) +
  geom_point(aes(color=conversion))+
  scale_y_continuous(trans='log10')+
  scale_color_manual(values = seroconversion_cols)+
  theme_minimal()+
  facet_wrap(~antigen, scales="free_y", ncol=3)+
  theme(axis.title = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(size=8, angle = 90, vjust=0.5)
        #axis.text.x = element_text(angle=30, vjust=1, hjust=1)
  )

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/full_virus_seroconversion_plot.png", full_virus_seroconversion_plot, width=6, height=8, dpi=444, bg="white")

#vaccine coverage plot ####

wide_long_msd %>%
  filter(antigen%in%c(vaccines, vaccines_iga[4:5]))%>%
  mutate(course_outcome=factor(course_outcome, levels=c("none", "partial", "no record", "complete")))%>%
  arrange(desc(course_outcome))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=titer))+
  # geom_segment(data = filter(line_data2, next_color=="nothing", antigen%in%c(vaccines, vaccines_iga[4:5])), aes(x = timepoint, y = titer, xend = next_x, yend = next_y, color = next_color)) +
  # geom_segment(data = filter(line_data2, next_color!="nothing", antigen%in%c(vaccines, vaccines_iga[4:5])), aes(x = timepoint, y = titer, xend = next_x, yend = next_y, color = next_color)) +
  geom_line(aes(color=course_outcome, group=id))+
  geom_point(aes(color=course_outcome))+
  scale_y_continuous(trans='log10')+
  scale_color_manual(values = c("red", "orange", "grey", "black"))+
  theme_minimal()+
  facet_wrap(~antigen, scales="free_y", nrow=6)+
  theme(axis.title = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        #axis.text.x = element_text(angle=30, vjust=1, hjust=1)
  )

# S4 treatmetn arm effects on viral seroconversion #####

seroconversion_table <- wide_long_msd%>%
  distinct(id.x, antigen, conversion, treatmentarm)%>%
  group_by(treatmentarm)%>%
  mutate(n_per_arm=n_distinct(id.x))%>%
  filter(conversion!="nothing", antigen %in% viruses)%>%
  group_by(antigen, conversion, treatmentarm, n_per_arm)%>%
  count()%>%
  arrange()%>%
  group_by(antigen, treatmentarm)%>%
  mutate(conversion_perc=n / n_per_arm)%>%
  arrange(desc(conversion_perc))


extra_rows <- seroconversion_table %>%
  group_by(antigen, treatmentarm)%>%
  summarise(conversion_perc=sum(conversion_perc))%>%
  mutate(
    conversion = "nothing",
    conversion_perc = 1 - conversion_perc,
    n = 89*conversion_perc,
    
  )%>%
  arrange(conversion_perc)

seroconversion_table2 <- seroconversion_table%>%
  bind_rows(., extra_rows)%>%
  arrange(antigen)

seroconvversion_barchart <- seroconversion_table2%>%
  ggplot(., aes(x=factor(antigen, levels=unique(extra_rows$antigen)), y=conversion_perc, fill=conversion))+
  geom_bar(stat="identity", position="stack")+
  scale_fill_manual(values = seroconversion_cols)+
  scale_y_continuous(labels = scales::label_percent())+
  ylab("Seroconversion")+
  facet_wrap(~treatmentarm)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=8, angle = 90, vjust=0.5),
        legend.position = "none")

# test for significanvr 
wide <- seroconversion_table2 %>%
  select(antigen, conversion, treatmentarm, n_per_arm, n) %>%
  pivot_wider(
    names_from = treatmentarm,
    values_from = c(n_per_arm, n)
  )

prop_test_fun <- function(x1, n1, x2, n2) {
  prop.test(
    x = c(x1, x2),
    n = c(n1, n2)
  )
}

safe_prop_test <- possibly(prop_test_fun, otherwise = NULL)


results <- wide %>%
  rowwise() %>%
  mutate(
    test = list(
      safe_prop_test(
        n_Placebo, n_per_arm_Placebo,
        `n_DP 1 year`, `n_per_arm_DP 1 year`
      )
    ),
    p_value = ifelse(is.null(test), NA, test$p.value),
    diff    = ifelse(is.null(test), NA, diff(test$estimate))
  ) %>%
  ungroup() %>%
  select(antigen, conversion, p_value, diff)


df <- read.csv("~/Downloads/df_coverage_with_recall.csv")

long_vaccine_coverage <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/long_vaccine_coverage_with_recall.csv")

vaccine_uptake_table <- long_vaccine_coverage%>%
  filter(id %in% msd_data$SubjectID)%>%
  group_by(vaccine_type, course_outcome)%>%
  count()%>%
  filter(!is.na(vaccine_type))%>%
  pivot_wider(names_from = course_outcome, values_from = n)%>%
  mutate(complete_perc=complete/n_distinct(msd_data$SubjectID, na.rm = T),
         none_perc=none/n_distinct(msd_data$SubjectID, na.rm = T),
         partial_perc=partial/n_distinct(msd_data$SubjectID, na.rm = T),
         some=sum(partial_perc, complete_perc, na.rm = T))%>%
  mutate(complete=paste(complete, "(", round(complete_perc*100, 1), "%)", sep = ""),
         partial=paste(partial, "(", round(partial_perc*100, 1), "%)", sep = ""),
         none=paste(none, "(", round(none_perc*100, 1), "%)", sep = ""))%>%
  tibble::column_to_rownames("vaccine_type")%>%
  select(complete, partial, none)

write.csv(vaccine_uptake_table, "~/postdoc/stanford/plasma_analytes/MICDROP/MSD/vaccine_uptake_table.csv", row.names = T)  

# S5 tetanus ####

tet52 <- antibodies_and_epi%>%
  filter(timepoint=="52 weeks", antigen=="Tetanus")%>%
  ggplot(., aes(x=factor(total_n_para_12), y=titer, fill=factor(total_n_para_12)))+
  geom_boxplot(outliers = F)+
  geom_point()+
  scale_y_log10()+
  scale_fill_manual(values = viridis::viridis(n=13))+
  theme_minimal()+
  ylab("concentration at 52 weeks [MSD units]")+
  xlab("number of parasitemic months in the first year of life")+
  theme(legend.position="none")

tet24 <- antibodies_and_epi%>%
  filter(timepoint=="24 weeks", antigen=="Tetanus")%>%
  ggplot(., aes(x=factor(total_n_para_6), y=titer, fill=factor(total_n_para_6)))+
  geom_boxplot(outliers = F)+
  geom_point()+
  scale_y_log10()+
  scale_fill_manual(values = viridis::viridis(n=13))+
  theme_minimal()+
  ylab("concentration at 24 weeks [MSD units]")+
  xlab("number of parasitemic months in the first 24 weeks of life")+
  theme(legend.position="none")

tet_combo <- tet24 / tet52
ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/S5_tet_combo.png", tet_combo, width=8, height=6, dpi=444, bg="white")


# SX paraseitemia~IgG sensitivity analysis ####


antibodies_and_epi %>%
  filter(antigen %in% c("PIV.1" ,    "PIV.2",      "Varicella", "Measles", "PIV.4", "Mumps",  "Rotavirus", "Rotavirus_IgA", "RV.C", "Diphtheria"),
         timepoint!="8 weeks")%>%
  filter(pardens!=0)%>%
  ggplot(., aes(x=qPCRparsdens+0.001 , y=titer+0.001, color=timepoint))+
  geom_smooth(method="lm", )+
  geom_point(alpha=0.68)+
  facet_wrap(~antigen, nrow=2)+
  scale_x_log10()+
  scale_y_log10()+
  xlab("qPCR parasitemia")+
  ylab("antibody concentration (MSD units)")+
  # ggpubr::stat_cor(aes(label = ..r.label..),method="spearman",
  #                  size=3.5, label.y=c(4.9, 5.25),label.x.npc=0.1, r.accuracy = 0.01, p.accuracy = 0.01)+
  ggpubr::stat_cor(method="spearman",
                   size=3.5, label.y=c(4.9, 5.25))+
  scale_color_manual(values = time_cols[2:3])+
  theme_minimal(base_size = 15)+
  theme(axis.text.x = element_text(angle=90, vjust=0.5),
        legend.position = "bottom")


