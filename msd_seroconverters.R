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
  mutate(conversion=case_when(timepoint=="52 weeks" & fc_flavour=="log2fc_24_52" & log2fc > 3.321928 ~ "converts 24 to 52",
                              timepoint=="24 weeks" & fc_flavour=="log2fc_8_24" & z_fc_8_24 > mean(z_fc_8_24)+2 | 
                              timepoint=="24 weeks" & fc_flavour=="log2fc_8_24" & log2fc > 1 ~ "converts 8 to 24",
                              .default="nothing"))


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


wide_long_msd%>%
  arrange(desc(conversion))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=titer))+
  geom_segment(data = filter(line_data, next_color=="nothing"), aes(x = timepoint, y = titer, xend = next_x, yend = next_y, color = next_color)) +
  geom_segment(data = filter(line_data, next_color!="nothing"), aes(x = timepoint, y = titer, xend = next_x, yend = next_y, color = next_color)) +
  geom_point(aes(color=conversion))+
  scale_y_continuous(trans='log10')+
  scale_color_manual(values = c("black", "red", "orange"))+
    theme_minimal()+
    facet_wrap(~antigen, scales="free")


serconversion_df <- wide_long_msd%>%
  filter(if (n()>1) conversion!="nothing" else TRUE)%>%
  mutate(conversion=gsub("nothing", "no", conversion))%>%
  distinct(id, antigen, timepoint, conversion)%>%
  pivot_wider(names_from =antigen , values_from = conversion, values_fill = "no")

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

purf%>%
  filter(padj<0.2)

nulisa_and_seroconversion%>%
  filter(targetName %in% "IL36B", timepoint=="52 weeks")%>%
  ggplot(., aes(x=timepoint, y=conc, fill=conversion))+
  geom_boxplot(outliers = F)+
  facet_wrap(~antigen, scales="free")+
  scale_fill_manual(values = c("orange", "red", "black"))+
  theme_minimal()
  

# treatment effect

converter_table <- wide_long_msd%>%
  filter(if (n()>1) conversion!="nothing" else TRUE)%>%
  mutate(conversion=gsub("nothing", "no", conversion))%>%
  filter(antigen%notin% vaccines)%>%
  distinct(id, antigen, conversion)%>%
  group_by(id)%>%
  summarise("8 to 24 weeks"=sum(conversion=="converts 8 to 24"),
            "24 to 52 weeks"=sum(conversion=="converts 24 to 52"))%>%
  pivot_longer(cols=c("8 to 24 weeks", "24 to 52 weeks"), names_to = "time_window", values_to = "seroconversions")%>%
  left_join(., mic_drop_key, by="id")%>%
  mutate(treatmentarm=case_match(treatmentarm,
                          1~"Placebo",
                          2~"DP 1 year",
                          3~"DP 2 years"))
  
converter_table%>%
  ggplot(., aes(x=time_window, y=seroconversions, fill=treatmentarm))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.height = 0.1))+
  theme_minimal()+
  theme(legend.position = "right")
  
  