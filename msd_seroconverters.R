library(tidyr)
library(dplyr)
library(xlsx)
library(ggplot2)
library(purrr)
library(emmeans)

nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")
mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")
maternal_treatment_arms <- haven::read_dta("~/Library/CloudStorage/Box-Box/DP+SP study/Databases and preliminary findings/Final database used for analyses/DPSP treatment allocation_FINAL.dta")

vaccines = c("Diptheria",     "Measles" ,      "Mumps",         "Pertussis",     "Polio",
             "Rotavirus" ,    "Rubella",       "Tetanus", "Pneumo.1.4.14")


nulisa_data <- nulisa_data%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))

metadata_columns <- c("id", "anyDP", "treatmentarm",  "dob", "date", "ageinwks", "gender_categorical", "mstatus", "qPCRparsdens", "visittype", "fever", "febrile", "rogerson", "GAcomputed", "gi", "SGA", "qPCRdich", "mqPCRparsdens")

epi_data <- nulisa_data%>%
  distinct(sample, total_n_para_12, total_n_malaria_12, total_n_malaria_6, total_n_para_6,
           id, dob, date, ageinwks, gender_categorical, mstatus, qPCRparsdens, visittype, fever, febrile, rogerson, GAcomputed, gi, SGA, qPCRdich, mqPCRparsdens, anyHP)%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"),
         mom_rx=maternal_treatment_arms$treatmentarm[match(id-10000, maternal_treatment_arms$id)],
         mom_rx=case_match(mom_rx,
                           1~"SP",
                           2~"DP",
                           3~"DPSP"))

msd_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/batch_one.csv")

long_msd <- msd_data%>%
  mutate(sample=paste(SubjectID, "_", "tp", TimePt, sep=""))%>%
  mutate(id=SubjectID, timepoint=paste(TimePt, "weeks"))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))%>%
  select(-SubjectID, -TimePt)%>%
  pivot_longer(cols=-c(sample, id, timepoint), names_to = "antigen", values_to = "titer")

antibodies_and_epi <- epi_data%>%
  inner_join(., long_msd, by="sample")

antibodies_and_nulisa <-  nulisa_data%>%
  select(-id, -timepoint)%>%
  inner_join(., long_msd, by="sample")

# defining seroconversion ####
#  between 8 and 24: calculate average fold change from 8 to 24; assume 2 STDvs above mean = serconversion even if decrease; else if increase at all;
#  between 24 and 52: simple 10fold increase

wide_long_msd <-  long_msd %>%
  filter(!is.na(timepoint))%>%
  pivot_wider(names_from = "timepoint", values_from = "titer", id_cols = c("id", "antigen"))%>%
  mutate(log2fc_8_24=log2(`24 weeks`/`8 weeks`),
         log2fc_24_52=log2(`52 weeks`/`24 weeks`))%>%
  group_by(antigen)%>%
  mutate("z_fc_8_24"=scale(log2fc_8_24, center = T))%>%
  pivot_longer(cols = c("8 weeks", "24 weeks","52 weeks"), names_to = "timepoint", values_to = "titer")%>%
  pivot_longer(cols = c("log2fc_8_24", "log2fc_24_52"), names_to = "fc_flavour", values_to = "log2fc")%>%
  mutate(conversion=case_when(timepoint=="52 weeks" & fc_flavour=="log2fc_24_52" & log2fc > 2 ~ "converts 24 to 52",
                              timepoint=="24 weeks" & fc_flavour=="log2fc_8_24" & z_fc_8_24 > 2 | 
                              timepoint=="24 weeks" & fc_flavour=="log2fc_8_24" & log2fc > 1 ~ "converts 8 to 24",
                              .default="nothing"))%>%
  mutate(mom_rx=maternal_treatment_arms$treatmentarm[match(id-10000, maternal_treatment_arms$id)],
         mom_rx=case_match(mom_rx,
                           1~"SP",
                           2~"DP",
                           3~"DPSP"))
  
table(wide_long_msd$conversion)
# two_sigma_cutoff <- wide_msd%>%
#   group_by(antigen)%>%
#   filter(z_fc_8_24>2)

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


big_plot <- wide_long_msd%>%
  arrange(desc(conversion))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=titer))+
  geom_segment(data = filter(line_data, next_color=="nothing"), aes(x = timepoint, y = titer, xend = next_x, yend = next_y, color = next_color)) +
  geom_segment(data = filter(line_data, next_color!="nothing"), aes(x = timepoint, y = titer, xend = next_x, yend = next_y, color = next_color)) +
  geom_point(aes(color=conversion))+
  scale_y_continuous(trans='log10')+
  scale_color_manual(values = c("black", "red", "orange"))+
  theme_minimal()+
  facet_wrap(~antigen, scales="free")+
  theme(axis.title = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/big_seroconversion_plot.png", big_plot, width=15, height=9, dpi=444, bg="white")


converters_only824 <- wide_long_msd%>%
  group_by(id, antigen)%>%
  filter(any(conversion=="converts 8 to 24"))%>%
  arrange(desc(conversion))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=titer))+
  geom_line(aes(group=id), color="red")+
  scale_y_continuous(trans='log10')+
  theme_minimal()+
  facet_wrap(~antigen, scales="free")+
  theme(axis.title = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/big_seroconverters824_only_plot.png", converters_only824, width=15, height=9, dpi=444, bg="white")


converters_only2452 <- wide_long_msd%>%
  group_by(id, antigen)%>%
  filter(any(conversion=="converts 24 to 52"))%>%
  arrange(desc(conversion))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=titer))+
  geom_line(aes(group=id), color="orange")+
  scale_y_continuous(trans='log10')+
  theme_minimal()+
  facet_wrap(~antigen, scales="free")+
  theme(axis.title = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/big_seroconverters2452_only.png", converters_only2452, width=15, height=9, dpi=444, bg="white")


non_converters <- wide_long_msd%>%
  group_by(id, antigen)%>%
  filter(all(conversion=="nothing"))%>%
  arrange(desc(conversion))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=titer))+
  geom_line(aes(group=id), color="black")+
  scale_y_continuous(trans='log10')+
  theme_minimal()+
  facet_wrap(~antigen, scales="free")+
  theme(axis.title = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/big_serononconverters_only.png", non_converters, width=15, height=9, dpi=444, bg="white")


any_converters <- wide_long_msd%>%
  group_by(id, antigen)%>%
  filter(any(conversion!="nothing"))%>%
  arrange(desc(conversion))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=titer))+
  geom_line(aes(group=id, color=conversion))+
  scale_color_manual(values = c("black", "orange", "red"))+
  scale_y_continuous(trans='log10')+
  theme_minimal()+
  facet_wrap(~antigen, scales="free")+
  theme(axis.title = element_blank(),
        legend.position = "none")

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/big_any_serononconverters_only.png", any_converters, width=15, height=9, dpi=444, bg="white")




serconversion_df <- wide_long_msd%>%
  filter(if (n()>1) conversion!="nothing" else TRUE)%>%
  mutate(conversion=gsub("nothing", "no", conversion))%>%
  distinct(id, antigen, timepoint, conversion)%>%
  pivot_wider(names_from =antigen , values_from = conversion, values_fill = "no")


serconversion_df2 <- wide_long_msd%>%
  # filter(if (n()>1) conversion!="nothing" else TRUE)%>%
  mutate(conversion=gsub("nothing", "no", conversion))%>%
  distinct(id, antigen, conversion)%>%
  pivot_wider(names_from = antigen , values_from = conversion)

antigens <- colnames(serconversion_df)[3:ncol(serconversion_df)]
timepoints <- c("8 weeks", "24 weeks", "52 weeks")

nulisa_and_seroconversion <-  nulisa_data%>%
  left_join(., serconversion_df, by=c("id", "timepoint"))%>%
  pivot_longer(cols = all_of(antigens), names_to = "antigen", values_to = "conversion")%>%
  select(id, timepoint, targetName, conc, antigen, conversion)%>%
  filter(!is.na(conversion))%>%
  filter(timepoint %in% timepoints)


for(a in antigens){
  
  plt = nulisa_and_seroconversion%>%
    filter(antigen==a)%>%
    ggplot(., aes(x=timepoint, y=conc, fill=conversion))+
    geom_boxplot(outliers = F)+
    facet_wrap(~targetName, scales="free")+
    ggtitle(paste0(a, " at ", t)) +
    scale_fill_manual(values = c("orange", "red", "black"))+
    theme_minimal()+
    theme(legend.position = "none")
  
  ggsave(paste0("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/serconversion_nulisa_cor_",a, ".png"),plt, width = 40, height = 40, dpi = 250, bg="white", limitsize = F)
  
}

# defining seropositivity ####
# mixture models didn't converge, so we can define seropositivity as the lowest level of what we have called
# seroconverters; also, Nguyen-Tran et al. have defined thresholds in their paper

seropositivity_thresholds <- wide_long_msd%>%
  filter(conversion=="converts 8 to 24" & timepoint=="24 weeks" | conversion=="converts 24 to 52" & timepoint=="52 weeks")%>%
  filter(log2fc>=0.5)%>%
  group_by(antigen)%>%
  slice_min(n = 1, with_ties = F, order_by = titer)%>%
  mutate("method"="sero_conversion")

seropositivity_thresholds2 <- seropositivity_thresholds %>%
  mutate(titer=case_when(antigen=="EV.D68"~1219, 
            antigen=="EV.71"~5275, 
            antigen=="hMPV"~4096,
            antigen=="RSV"~2724,
            antigen=="Flu.Mix"~2070,
            antigen=="RV.C"~166,
            antigen=="PIV.Mix"~1257,
            antigen=="Flu.A.H1"~9490,
            antigen=="Flu.A.H3"~13498,
            antigen=="Flu.B"~8209,
            antigen=="PIV.1"~21190,
            antigen=="PIV.2"~7992,
            antigen=="PIV.3"~18371,
            antigen=="PIV.4"~26057, .default = NA))%>%
  mutate("method"="nguyen-tran_et_al")%>%
  filter(!is.na(titer))


seroprev_est_df <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/mixture_model_df.csv")

colnames(seroprev_est_df) <- c("antigen", "prevalence", "titer")
seroprev_est_df$method="mixture_model"
seroprev_est_df <- subset(seroprev_est_df, !is.na(titer))
seroprev_est_df$titer <- 10^seroprev_est_df$titer

wide_long_msd%>%
  arrange(desc(conversion))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=titer))+
  geom_hline(data=seropositivity_thresholds, aes(yintercept = titer), linetype="dashed", color="blue")+
  geom_hline(data=seropositivity_thresholds2, aes(yintercept = titer), linetype="dotted", color="blue")+
  geom_hline(data=seroprev_est_df, aes(yintercept = titer), linetype="solid", color="blue")+
  geom_segment(data = filter(line_data, next_color=="nothing"), aes(x = timepoint, y = titer, xend = next_x, yend = next_y, color = next_color)) +
  geom_segment(data = filter(line_data, next_color!="nothing"), aes(x = timepoint, y = titer, xend = next_x, yend = next_y, color = next_color)) +
  geom_point(aes(color=conversion))+
  scale_y_continuous(trans='log10')+
  scale_color_manual(values = c("black", "red", "orange"))+
  theme_minimal()+
  facet_wrap(~antigen, scales="free")+
  theme(axis.title = element_blank(),
        legend.title = element_blank())


wide_long_msd%>%
  arrange(desc(conversion))%>%
  ggplot(., aes(fill=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), x=titer))+
  geom_vline(data=seropositivity_thresholds, aes(xintercept = titer), linetype="dashed", color="red")+
  geom_vline(data=seropositivity_thresholds2, aes(xintercept = titer), linetype="dotted", color="red")+
  geom_vline(data=seroprev_est_df, aes(xintercept = titer), linetype="solid", color="red")+
  geom_histogram()+
  scale_fill_manual(values=colorspace::sequential_hcl(n=4, "RdPu")[1:3])+
  scale_x_continuous(trans="log10")+
  theme_minimal()+
  facet_wrap(~antigen, scales="free")+
  theme(axis.title = element_blank(),
        legend.title = element_blank())

final_cutoff_frame <- bind_rows(seropositivity_thresholds, seropositivity_thresholds2, seroprev_est_df)%>%
  select(antigen, titer, method)%>%
  group_by(antigen)%>%
  filter( min_rank(case_when(
    method == "mixture_model" ~ 1,
    method == "nguyen-tran_et_al" ~ 2,
    TRUE ~ 3
  )) == 1)

write.csv(final_cutoff_frame, "~/postdoc/stanford/plasma_analytes/MICDROP/MSD/final_cutoff_frame.csv", row.names = F)









wide_long_msd%>%
  arrange(desc(conversion))%>%
  ggplot(., aes(fill=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), x=titer))+
  geom_vline(data=seropositivity_thresholds, aes(xintercept = titer), linetype="dashed", color="red")+
  geom_vline(data=seropositivity_thresholds2, aes(xintercept = titer), linetype="dotted", color="red")+
  geom_vline(data=seroprev_est_df, aes(xintercept = titer), linetype="solid", color="red")+
  geom_histogram()+
  scale_fill_manual(values=colorspace::sequential_hcl(n=4, "RdPu")[1:3])+
  scale_x_continuous(trans="log10")+
  theme_minimal()+
  facet_wrap(~antigen, scales="free")+
  theme(axis.title = element_blank(),
        legend.title = element_blank())


# pairwise comparisons stats ####

purf <- nulisa_and_seroconversion%>%
  mutate(bino_convert=if_else(conversion=="no", as.numeric(0), as.numeric(1)))%>%
  group_by(targetName, antigen, timepoint)%>%
  nest()%>%
  mutate(model = purrr::map(data, ~glm(bino_convert~conc,  family = "binomial",  data=.)))%>%
  mutate(model_summary=purrr::map(model, ~summary(.)))%>%
  mutate(p=purrr::map(model_summary, ~coef(.)[8]))%>%
  group_by(timepoint, targetName)%>%
  mutate(padj=p.adjust(p, method="fdr"))

sig_purf_conversion_nulisa <- purf%>%
  filter(padj<0.15)

polio_analytes <- purf%>%
  filter(padj<0.15, antigen=="Polio")

seroconverter_hits <- nulisa_and_seroconversion%>%
  semi_join(., sig_purf_conversion_nulisa, by=c("targetName", "antigen"))%>%
  # filter(targetName %in% polio_analytes$targetName)%>%
  # filter(antigen=="Polio")%>%
  ggplot(., aes(x=timepoint, y=conc, fill=conversion))+
  geom_boxplot(outliers = F)+
  ggpubr::stat_compare_means(label="p.format", size=2, vjust=.3)+
  facet_wrap(~targetName+antigen, scales="free", nrow=2)+
  scale_fill_manual(values = c("orange", "red", "black"))+
  theme_minimal()+
  theme(legend.position = "bottom")
  
ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/seroconverter_nulisa_sigs.png", seroconverter_hits, width=8, height=4, bg="white", dpi=444)

# treatment effect

converter_table <- wide_long_msd%>%
  filter(if (n()>1) conversion!="nothing" else TRUE)%>%
  mutate(conversion=gsub("nothing", "no", conversion))%>%
  filter(antigen %in% vaccines)%>%
  distinct(id, antigen, conversion)%>%
  group_by(id)%>%
  summarise("8 to 24 weeks"=sum(conversion=="converts 8 to 24"),
            "24 to 52 weeks"=sum(conversion=="converts 24 to 52"))%>%
  pivot_longer(cols=c("8 to 24 weeks", "24 to 52 weeks"), names_to = "time_window", values_to = "seroconversions")%>%
  left_join(., mic_drop_key, by="id")%>%
  left_join(., df_clusters, by="id")%>%
  mutate(treatmentarm=case_match(treatmentarm,
                          1~"Placebo",
                          2~"DP 1 year",
                          3~"DP 2 years"))
  
converter_table%>%
  ggplot(., aes(x=time_window, y=seroconversions, fill=Cluster))+
  geom_violin(draw_quantiles = seq(0,1,0.25), color="white")+
  geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.height = 0.1))+
  theme_minimal()+
  ggpubr::stat_compare_means(label = "p.format")+
  ggtitle("seroconversion to any vaccine antigen(s)")+
  theme(legend.position = "right")
  


long_msd <- msd_data%>%
  mutate(sample=paste(SubjectID, "_", "tp", TimePt, sep=""))%>%
  mutate(id=SubjectID, timepoint=paste(TimePt, "weeks"))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))%>%
  select(-SubjectID, -TimePt)%>%
  pivot_longer(cols=-c(sample, id, timepoint), names_to = "antigen", values_to = "titer")%>%
  left_join(., df_clusters, by="id")

long_msd%>%
  filter(timepoint=="52 weeks")%>%
  # mutate(any_para = if_else(total_n_para_6 > 0, 1, 0),
  #        any_malaria = if_else(total_n_malaria_6 > 0, 1, 0))%>%
  ggplot(., aes(x=factor(timepoint), y=titer, fill=factor(Cluster)))+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  geom_boxplot()+
  scale_y_continuous(trans="log10")+
  ggpubr::stat_compare_means(vjust = 0.5)+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()


trajectory_antibodies <- long_msd%>%
  filter(timepoint=="52 weeks", antigen %in% c("Flu.A.H1", "Mumps", "PIV.4", "hMPV"))%>%
  # mutate(any_para = if_else(total_n_para_6 > 0, 1, 0),
  #        any_malaria = if_else(total_n_malaria_6 > 0, 1, 0))%>%
  ggplot(., aes(x=factor(timepoint), y=titer, fill=factor(Cluster)))+
  geom_violin(draw_quantiles = seq(0,1,0.25))+
  # geom_boxplot()+
  scale_y_continuous(trans="log10")+
  ggpubr::stat_compare_means(vjust = -0.2, label = "p.format")+
  facet_wrap(~antigen, scales="free", nrow=1)+
  theme_minimal()+
  theme(legend.title = element_blank(),
        axis.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/trajectory_antibodies.png",trajectory_antibodies, width=8, height=4 , bg="white")




serconversion_df3 <- wide_long_msd%>%
  filter(if (n()>1) conversion!="nothing" else TRUE)%>%
  left_join(., df_clusters, by="id")%>%
  mutate(conversion=gsub("nothing", "no", conversion))%>%
  distinct(id, antigen, timepoint, conversion, Cluster)%>%
  group_by(antigen, Cluster, conversion)%>%
  summarise("n"=n())
# pivot_wider(names_from = antigen , values_from = conversion, values_fill = "no")

ggplot(serconversion_df3, aes(x=factor(Cluster), y=n, fill=conversion))+
  geom_bar(stat="identity")+
  facet_wrap(~antigen)+
  scale_fill_manual(values=c("orange", "red", "black"))+
  theme_minimal()
  

# maternal rx ####
serconversion_df <- wide_long_msd%>%
  distinct(id, conversion, antigen, mom_rx)

serconversion_df%>%
  group_by(conversion, antigen, mom_rx)%>%
  summarise("n"=n(), "freq"=n/70)%>%
  ggplot(., aes(x=antigen, y=freq, fill=mom_rx))+
  geom_bar(stat="identity")+
  facet_wrap(~conversion)

#tbd
serocon_momrx_chisq <- serconversion_df %>%
  group_by(antigen) %>%
  nest() %>%
  mutate(
    test = map(
      data,
      ~ {
        tbl <- with(.x, table(conversion, mom_rx))
        chisq.test(tbl)
      }
    ),
    tidied = map(test, broom::tidy)
  ) %>%
  unnest(tidied) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "fdr"))


neuro_cog <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/NCT_infants_share032025.csv")
neuro_cog_edit <- neuro_cog%>%
  mutate(id=subjid, date=as.Date(vdate))%>%
  select(-gender, -GAcomputed, -SGA, -date, -dob)%>%
  pivot_longer(cols = ends_with("composite"), names_to = "composite_kind", values_to = "composite_score")

seroconv_plus_neuro <- left_join(serconversion_df, neuro_cog_edit, by="id")


seropos_neuro_cog_purrr <- seroconv_plus_neuro%>%
  filter(!is.na(composite_kind))%>%
  mutate(conversion=factor(conversion, levels=c("nothing", "converts 8 to 24", "converts 24 to 52")))%>%
  group_by(antigen, composite_kind)%>%
  nest()%>%
  mutate(model = purrr::map(data, ~glm(composite_score~conversion, data=.)))%>%
  mutate(model_summary=purrr::map(model, ~summary(.)))%>%
  mutate(p24=purrr::map(model_summary, ~coef(.)[11]))%>%
  mutate(p52=purrr::map(model_summary, ~coef(.)[12]))%>%
  ungroup(antigen)%>%
  mutate(padj24=p.adjust(p24, method="fdr"))%>%
  mutate(padj52=p.adjust(p52, method="fdr"))

seropos_neuro_cog_purrr%>%
  filter(padj24<0.05 | padj52<0.05)


# get all levels of conversion
df <- seroconv_plus_neuro%>%
  filter(!is.na(composite_kind))%>%
  filter(antigen=="PIV.Mix")%>%
  mutate(conversion=factor(conversion))
  
conv_levels <- levels(df$conversion)
# define comparisons vs "nothing"

comparisons <- lapply(setdiff(conv_levels, "nothing"), function(x) c("nothing", x))

seroconv_plus_neuro%>%
  filter(!is.na(composite_kind))%>%
  filter(antigen=="PIV.Mix")%>%
  ggplot(., aes(x=composite_kind, y=composite_score, fill=factor(conversion)))+
  geom_boxplot(outliers = F)+
  theme_minimal()+
  ggpubr::stat_compare_means()





