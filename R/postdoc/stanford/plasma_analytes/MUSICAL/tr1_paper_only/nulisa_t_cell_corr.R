#preamble ####
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)


da_boxplot_theme <- theme(legend.position = "none",
                          axis.title = element_blank())

stim_palette <- c("darkred", "darkblue", "black")
names(stim_palette) <- c("iRBCs", "PMA", "unstim")

time_cols <- list("baseline"="#E4DEBD",
                  "day0" = "#C03F3E",
                  "day7" = "#D87E1F",
                  "day14" = "#E6B85F")

`%notin%` <- Negate(`%in%`)

# data generation ####
nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/tr1_paper/revised_baseline_clean_musical_combo_with_metadata.csv")
nulisa_data <- nulisa_data%>%
  mutate(infectiontype=substr(infectiontype, 1, 1))
  

slim_nulisa_data <- nulisa_data %>%
  mutate(id=as.character(id))%>%
  select(id, date, timepoint, timepoint_imm, infectiontype, targetName, concentration, log_qpcr, ageyrs, gender_categorical)%>%
  filter(timepoint!="bad_baseline")%>%
  mutate(timepoint_imm=if_else(timepoint=="baseline", -1, timepoint_imm))

# cell_count_data <- read.csv("~/postdoc/stanford/manuscripts/jason_tr1_2/t_cell_cytometry_count_data.csv")
slim_cell_count_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/df_jason_analysis.csv")

slim_cell_count_data <- slim_cell_count_data%>%
  pivot_longer(cols = ends_with("Frequency"), names_to = "gate", values_to = "freq")%>%
  mutate(timepoint_imm=timepoint)%>%
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

# write.csv(long_combo, "~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/cell_freqs_and_nulisa.csv", col.names = F)

# grand correlation 

grand_cor <- long_combo%>%
  filter(gate=="Tr1_Frequency", stim=="unstim", infectiontype%in% c("A", "S"), !is.na(targetName))%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(correlation=map(data, ~cor.test(.$concentration, .$freq, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  ungroup()%>%
  mutate(padj=p.adjust(p))


real_count_data <- read.csv("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/t_cell_cytometry_count_data.csv")
real_count_data <- real_count_data%>%
  mutate(infectiontype=toupper(substr(infectiontype, 1, 1)))
  
  
count_nulisa_combo <- slim_nulisa_data%>%
  mutate(id=as.integer(id))%>%
  full_join(., real_count_data, by = c("id", "timepoint", "infectiontype"))


grand_count_cor <- count_nulisa_combo%>%
  filter(gate=="Tr1_Frequency", stim=="unstim", infectiontype%in% c("A", "S"), !is.na(targetName))%>%
  distinct(id, timepoint, infectiontype, absolute_Tr1, targetName, concentration)%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(correlation=map(data, ~cor.test(.$concentration, .$absolute_Tr1, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  ungroup()%>%
  mutate(padj=p.adjust(p))

#noting
sig_count_cor <- grand_count_cor%>%
  filter(p<0.05)



# disjointed timepoints, absolute values ####
day0_concs <- long_combo %>%
  filter(infectiontype%in% c("S", "A"))%>%
  filter(timepoint=="day0")%>%
  distinct(targetName, concentration, id, timepoint, infectiontype)

day7_concs <- long_combo %>%
  filter(infectiontype%in% c("S", "A"))%>%
  filter(timepoint=="day7")%>%
  distinct(targetName, concentration, id, timepoint, infectiontype)

base_concs <- long_combo %>%
  filter(infectiontype%in% c("S", "A"))%>%
  filter(timepoint=="baseline")%>%
  distinct(targetName, concentration, id, timepoint, infectiontype)

base_freqs <- long_combo%>%
  filter(infectiontype%in% c("S", "A"))%>%
  filter(timepoint=="baseline")%>%
  filter(gate=="Tr1_Frequency", stim=="unstim")%>%
  distinct(gate, freq, stim, id, timepoint, infectiontype)

day0_freqs <- long_combo%>%
  filter(infectiontype%in% c("S", "A"))%>%
  filter(timepoint=="day0")%>%
  filter(gate=="Tr1_Frequency", stim=="unstim")%>%
  distinct(gate, freq, stim, id, timepoint, infectiontype)

day7_freqs <- long_combo%>%
  filter(infectiontype%in% c("S", "A"))%>%
  filter(infectiontype=="S", timepoint=="day7")%>%
  filter(gate=="Tr1_Frequency", stim=="unstim")%>%
  distinct(gate, freq, stim, id, timepoint, infectiontype)

day14_freqs <- long_combo%>%
  filter(infectiontype%in% c("S", "A"))%>%
  filter(timepoint=="day14")%>%
  filter(gate=="Tr1_Frequency", stim=="unstim")%>%
  distinct(gate, freq, stim, id, timepoint, infectiontype)

base_freq_base_conc <- inner_join(base_concs, base_freqs, by=c("id", "infectiontype"), relationship = "many-to-many")

base_freq_day0_conc <- inner_join(day0_concs, base_freqs, by=c("id", "infectiontype"), relationship = "many-to-many")

conc_freq_0 <- inner_join(day0_concs, day0_freqs, by=c("id", "infectiontype"), relationship = "many-to-many")

conc_freq_014 <- inner_join(day0_concs, day14_freqs, by=c("id", "infectiontype"), relationship = "many-to-many")


conc_freq_07 <- day0_concs%>%
  filter(infectiontype=="S")%>%
  inner_join(day7_freqs, by=c("id", "infectiontype"), relationship = "many-to-many")


conc_freq_77 <- day7_concs%>%
  filter(infectiontype=="S")%>%
  inner_join(day7_freqs, by=c("id", "infectiontype"), relationship = "many-to-many")

freq_freq_base_0 <- inner_join(base_freqs, day0_freqs, by=c("id", "gate", "stim", "infectiontype"), relationship = "many-to-many")


conc_freq_0_corr <- conc_freq_0%>%
  filter(!is.na(freq))%>%
  group_by(gate, stim, targetName, infectiontype)%>%
  nest()%>%
  mutate(nrows=map_dbl(data, ~nrow(.)))%>%
  filter(nrows>1)%>%
  mutate(correlation=map(data, ~cor.test(.$concentration, .$freq, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  ungroup()%>%
  group_by(infectiontype)%>%
  mutate(padj=p.adjust(p))


#CALCA & CXCL10 significnat!
conc_freq_07_corr <- conc_freq_07%>%
  filter(!is.na(freq), infectiontype=="S")%>%
  group_by(gate, stim, targetName)%>%
  nest()%>%
  mutate(nrows=map_dbl(data, ~nrow(.)))%>%
  filter(nrows>1)%>%
  mutate(correlation=map(data, ~cor.test(.$concentration, .$freq, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  ungroup()%>%
  mutate(padj=p.adjust(p))

conc_freq_77_corr <- conc_freq_77%>%
  filter(!is.na(freq), infectiontype=="S")%>%
  group_by(gate, stim, targetName)%>%
  nest()%>%
  mutate(nrows=map_dbl(data, ~nrow(.)))%>%
  filter(nrows>1)%>%
  mutate(correlation=map(data, ~cor.test(.$concentration, .$freq, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  ungroup()%>%
  mutate(padj=p.adjust(p))

#ANGPT maybe? correlates
base_freq_base_conc_corr <- base_freq_base_conc%>%
  filter(!is.na(freq))%>%
  group_by(gate, stim, targetName, infectiontype)%>%
  nest()%>%
  mutate(nrows=map_dbl(data, ~nrow(.)))%>%
  filter(nrows>1)%>%
  mutate(correlation=map(data, ~cor.test(.$concentration, .$freq, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  ungroup()%>%
  mutate(padj=p.adjust(p))

#nothing
base_freq_day0_conc_corr <- base_freq_day0_conc%>%
  filter(!is.na(freq))%>%
  group_by(gate, stim, targetName, infectiontype)%>%
  nest()%>%
  mutate(nrows=map_dbl(data, ~nrow(.)))%>%
  filter(nrows>1)%>%
  mutate(correlation=map(data, ~cor.test(.$concentration, .$freq, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  ungroup()%>%
  mutate(padj=p.adjust(p))

#nothing
conc_freq_014_corr <- conc_freq_014%>%
  filter(!is.na(freq))%>%
  group_by(gate, stim, targetName)%>%
  nest()%>%
  mutate(nrows=map_dbl(data, ~nrow(.)))%>%
  filter(nrows>1)%>%
  mutate(correlation=map(data, ~cor.test(.$concentration, .$freq, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  ungroup()%>%
  mutate(padj=p.adjust(p))
  

# tr1 freqs at baseline DO correlate with day0 freqs during symptomatic infection;
# more so during asymp, which makes sense, because they don't change much'
freq_freq_base_0_corr <- freq_freq_base_0%>%
  filter(!is.na(freq.x), !is.na(freq.y))%>%
  filter(gate %notin% c("Memory_CD4_T_Cell_Frequency",
                        "CD4_Lymphocyte_Frequency",
                        unique(grep("FOXP3", gate, value=TRUE))))%>%
  group_by(gate, stim, infectiontype)%>%
  nest()%>%
  mutate(nrows=map_dbl(data, ~nrow(.)))%>%
  filter(nrows>1)%>%
  mutate(correlation=map(data, ~cor.test(.$freq.x, .$freq.y, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  ungroup()%>%
  group_by()%>%
  mutate(padj=p.adjust(p))


# fold changes ####
# you need exactly the same number of observations between x and y for spearman;
# this is an issue with many missing samples in flow and some missing smaples in NULISA leading to infinites
# need to fix because otherwise p value calculation is garbage.

wide_concs <-  long_combo %>%
  select(targetName, concentration, id, infectiontype, timepoint)%>%
  filter(infectiontype!="nmf", timepoint %in% c("baseline", "day0", "day7", "day14"))%>%
  distinct(targetName, concentration, id, infectiontype, timepoint)%>%
  pivot_wider(names_from = timepoint, values_from = c(concentration), id_cols=c(targetName, id, infectiontype), names_prefix = "conc_")%>%
  group_by(infectiontype)%>%
  mutate(conc_base_d0_fc=conc_day0-conc_baseline,
         conc_base_d7_fc=conc_day7-conc_baseline,
         conc_base_d14_fc=conc_day14-conc_baseline)%>%
  distinct(targetName, id, infectiontype, conc_base_d0_fc, conc_base_d14_fc, conc_base_d7_fc)
  

wide_freqs <-  long_combo %>%
  # do some fancy subsetting so that phenotypic gates are only present from unstim samples and cytokine gates only from stim conditions
  # filter(!(stim=="PMA" & !grepl("_[^_]+_", gate)),
  #        !(stim=="iRBCs" & !grepl("_[^_]+_", gate)),
  #        #!(stim=="unstim" & grepl("_[^_]+_", gate)),
  #        #!(stim=="PMA" & gate=="CD4_T_Cell_Frequency")
  #        )%>%
  filter(gate%in%c("Tr1_Frequency", "Tr1_IL21_Frequency", "Tr1_IL10_Frequency", "Tr1_IFNg_Frequency", "IFNg_Frequency", "IL10_Frequency", "IL21_Frequency"))%>%
  select(gate, stim, freq, id, infectiontype, timepoint)%>%
  filter(infectiontype!="NM", timepoint %in% c("baseline", "day0", "day7", "day14"))%>%
  distinct(gate, stim, freq, id, infectiontype, timepoint)%>%
  pivot_wider(names_from = timepoint, values_from = c(freq), id_cols=c(gate, stim, id, infectiontype), names_prefix = "freq_")%>%
  group_by(infectiontype)%>%
  mutate(freq_baseline=freq_baseline+0.001)%>%
  mutate(freq_base_d0_fc=freq_day0/freq_baseline,
         freq_base_d7_fc=freq_day7/freq_baseline,
         freq_base_d14_fc=freq_day14/freq_baseline)%>%
  distinct(gate, stim, id, infectiontype, freq_base_d0_fc, freq_base_d7_fc, freq_base_d14_fc)


fc_combo_frame <- inner_join(wide_freqs, wide_concs, by=c("id", "infectiontype"))%>%
  filter(is.finite(freq_base_d0_fc),
         is.finite(freq_base_d14_fc),
         is.finite(conc_base_d0_fc),
         is.finite(conc_base_d14_fc))



fc_corr_purr_S <- fc_combo_frame %>%
  filter(targetName!="CTSS")%>%
  filter(infectiontype%in%c("S")
         # ,gate%notin%c("Memory_CD4_T_Cell_Frequency",
         #                                  "CD4_Lymphocyte_Frequency",
         #                                  "CCR4_P_CXCR6_N_IFNg_Frequency",
         #                                  "CCR4_P_CXCR6_N_IL21_Frequency",
         #                                  unique(grep("FOXP3", gate, value=TRUE)))
  )%>%
  # filter(!is.infinite(freq_base_d7_fc), !is.infinite(freq_base_d14_fc), !is.infinite(conc_base_d0_fc))%>%
  # filter(!is.na(freq_base_d7_fc), !is.na(freq_base_d14_fc), !is.na(conc_base_d0_fc))%>%
  group_by(infectiontype, stim, gate, targetName)%>%
  nest()%>%
  mutate(base77_correlation=map(data, ~cor.test(.$freq_base_d7_fc, .$conc_base_d7_fc, method = "spearman")))%>%
  mutate(base0_correlation=map(data, ~cor.test(.$freq_base_d0_fc, .$conc_base_d0_fc, method = "spearman")))%>%
  mutate(base7_correlation=map(data, ~cor.test(.$freq_base_d7_fc, .$conc_base_d0_fc, method = "spearman")))%>%
  mutate(base14_correlation=map(data, ~cor.test(.$freq_base_d14_fc, .$conc_base_d0_fc, method = "spearman")))%>%
  mutate(
    base77_p=map_dbl(base77_correlation, ~.$p.value),
    base77_rho=map_dbl(base77_correlation, ~.$estimate),
    
    base0_p=map_dbl(base0_correlation, ~.$p.value),
    base0_rho=map_dbl(base0_correlation, ~.$estimate),
    base7_p=map_dbl(base7_correlation, ~.$p.value),
    base7_rho=map_dbl(base7_correlation, ~.$estimate),
    base14_p=map_dbl(base14_correlation, ~.$p.value),
    base14_rho=map_dbl(base14_correlation, ~.$estimate))%>%
  group_by(infectiontype, stim, targetName)%>%
  mutate(base77_padj=p.adjust(base77_p),
         base14_padj=p.adjust(base14_p),
         base7_padj=p.adjust(base7_p),
         base0_padj=p.adjust(base0_p))%>%
  pivot_longer(cols=ends_with(c("padj", "rho")), names_to = c("Var", ".value"), names_sep  = "_")



  

fc_corr_purr_A <- fc_combo_frame %>%
  filter(infectiontype%in%c("A")) %>%
  group_by(infectiontype, stim, gate, targetName)%>%
  nest()%>%
  # mutate(base_d7_correlation=if_else(all(is.na(data$freq_day7)), map(data, ~cor.test(.$freq_base_d7_fc, .$conc_base_d0_fc, method = "spearman"))))%>%
  mutate(base0_correlation=map(data, ~cor.test(.$freq_base_d0_fc, .$conc_base_d0_fc, method = "spearman")))%>%
  mutate(base14freq_base0_conc_correlation=map(data, ~cor.test(.$freq_base_d14_fc, .$conc_base_d0_fc, method = "spearman")))%>%
  mutate(base14_correlation=map(data, ~cor.test(.$freq_base_d14_fc, .$conc_base_d14_fc, method = "spearman")))%>%
  mutate(
    base0_p=map_dbl(base0_correlation, ~.$p.value),
    base0_rho=map_dbl(base0_correlation, ~.$estimate),
    base14_p=map_dbl(base14_correlation, ~.$p.value),
    base14_rho=map_dbl(base14_correlation, ~.$estimate),
    base14freq_base0_conc_correlation_p=map_dbl(base14freq_base0_conc_correlation, ~.$p.value),
    base14freq_base0_conc_correlation_rho=map_dbl(base14freq_base0_conc_correlation, ~.$estimate))%>%
  group_by(infectiontype, stim, targetName)%>%
  mutate(base0_padj=p.adjust(base0_p),
         base14_padj=p.adjust(base14_p),
         base14freq_base0_conc_correlation_padj=p.adjust(base14freq_base0_conc_correlation_p))%>%
  pivot_longer(cols=ends_with(c("padj", "rho")), names_to = c("Var", ".value"), names_sep  = "_")

fc_corr_purr <- bind_rows(fc_corr_purr_A, fc_corr_purr_S)



sig_fc_corr_purr <- fc_corr_purr%>%
  select(gate, stim, infectiontype, targetName, Var, rho, padj)%>%
  # filter(!duplicated)
  filter(padj<0.1)

## day 0 fold change plots ####
### Tr1 frequency ####
#analytes where FC in concetration from base to day0 correlates with FC from base to day0 in Tr1 frequency
base_d0_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$stim=="unstim"&sig_fc_corr_purr$gate=="Tr1_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="S"]

base_d0_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("Tr1_Frequency"), stim=="unstim", targetName%in%base_d0_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d0_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in Tr1 frequency from baseline to day 0")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous(trans="log2")+
  facet_wrap(~targetName)+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d0_sigs_plot.png", base_d0_sigs_plot, height = 6, width=6, dpi=444, bg="white")


### il10 frequency ####
base_d0_il10_s_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IL10_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="PMA"]
base_d0_il10_s_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IL10_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="iRBCs"]

# sig_fc_corr_purr[sig_fc_corr_purr$stim=="PMA"&sig_fc_corr_purr$gate=="IL10_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="S",]

base_d0_il10_s_pma_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("IL10_Frequency"), stim=="PMA", targetName%in%base_d0_il10_s_pma_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d0_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman", size=3)+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IL10+ CD4 memory frequency from baseline to day 0")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous(trans="log2")+
  facet_wrap(~targetName, ncol=5)+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d0_il10_s_pma_sigs_plot.png", base_d0_il10_s_pma_sigs_plot, height = 8, width=8, dpi=444, bg="white")



base_d0_il10_s_irbc_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("IL10_Frequency"), stim=="iRBCs", targetName%in%base_d0_il10_s_irbc_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d0_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IL10+ CD4 memory frequency from baseline to day 0")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous(trans="log2")+
  facet_wrap(~targetName, nrow=2)+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d0_il10_s_irbc_sigs_plot.png", base_d0_il10_s_irbc_sigs_plot, height = 6, width=6, dpi=444, bg="white")


#noting
base_d0_il10_a_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IL10_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="PMA"]
base_d0_il10_a_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IL10_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="iRBCs"]


base_d0_il10_a_irbc_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("IL10_Frequency"), stim=="iRBCs", targetName%in%base_d0_il10_a_irbc_sigs, infectiontype=="A")%>%
  ggplot(., aes(x=freq_base_d0_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Asymptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IL10+ CD4 memory frequency from baseline to day 0")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous(trans="log2")+
  facet_wrap(~targetName, nrow=3, scales="free")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d0_il10_a_irbc_sigs_plot.png", base_d0_il10_a_irbc_sigs_plot, height = 6, width=6, dpi=444, bg="white")



### tr1_il10 ####

base_d0_tr1_il10_s_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IL10_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="PMA"]
base_d0_tr1_il10_s_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IL10_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="iRBCs"]

# sig_fc_corr_purr[sig_fc_corr_purr$stim=="PMA"&sig_fc_corr_purr$gate=="Tr1_IL10_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="S",]

base_d0_tr1_il10_s_pma_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("Tr1_IL10_Frequency"), stim=="PMA", targetName%in%base_d0_tr1_il10_s_pma_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d0_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IL10+ Tr1 frequency from baseline to day 0")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous(trans="log2")+
  facet_wrap(~targetName, nrow=3, scales = "free")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d0_tr1_il10_s_pma_sigs_plot.png", base_d0_tr1_il10_s_pma_sigs_plot, height = 8, width=16, dpi=444, bg="white")




base_d0_tr1_il10_a_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IL10_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="PMA"]
base_d0_tr1_il10_a_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IL10_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="iRBCs"]


base_d0_tr1_il10_a_irbc_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("Tr1_IL10_Frequency"), stim=="iRBCs", targetName%in%base_d0_tr1_il10_a_irbc_sigs, infectiontype=="A")%>%
  ggplot(., aes(x=freq_base_d0_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Asymptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IL10+ Tr1 frequency from baseline to day 0")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous(trans="log2")+
  facet_wrap(~targetName, nrow=1, scales="free")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d0_tr1_il10_a_irbc_sigs_plot.png", base_d0_tr1_il10_a_irbc_sigs_plot, height = 5, width=8, dpi=444, bg="white")

### ifng frequency ####

base_d0_ifng_s_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IFNg_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="PMA"]
base_d0_ifng_s_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IFNg_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="iRBCs"]

# sig_fc_corr_purr[sig_fc_corr_purr$stim=="PMA"&sig_fc_corr_purr$gate=="IFNg_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="S",]

base_d0_ifng_s_pma_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("IFNg_Frequency"), stim=="PMA", targetName%in%base_d0_ifng_s_pma_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d0_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IFNg+ CD4 memory frequency from baseline to day 0")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous(trans="log2")+
  facet_wrap(~targetName, nrow=2, scales = "free")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d0_ifng_s_pma_sigs_plot.png", base_d0_ifng_s_pma_sigs_plot, height = 6, width=9, dpi=444, bg="white")



base_d0_ifng_s_irbc_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("IFNg_Frequency"), stim=="iRBCs", targetName%in%base_d0_ifng_s_irbc_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d0_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IFNg+ CD4 memory frequency from baseline to day 0")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous(trans="log2")+
  facet_wrap(~targetName, nrow=1, scales="free")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d0_ifng_s_irbc_sigs_plot.png", base_d0_ifng_s_irbc_sigs_plot, height = 6, width=6, dpi=444, bg="white")


#noting
base_d0_ifng_a_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IFNg_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="PMA"]
base_d0_ifng_a_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IFNg_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="iRBCs"]


base_d0_ifng_a_irbc_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("IFNg_Frequency"), stim=="iRBCs", targetName%in%base_d0_ifng_a_irbc_sigs, infectiontype=="A")%>%
  ggplot(., aes(x=freq_base_d0_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Asymptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IFNg+ CD4 memory frequency from baseline to day 0")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous (trans="log2")+
  facet_wrap(~targetName, nrow=3, scales="free")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d0_ifng_a_irbc_sigs_plot.png", base_d0_ifng_a_irbc_sigs_plot, height = 8, width=12, dpi=444, bg="white")


### Tr1 IFNg frequency ####


base_d0_tr1_ifng_s_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IFNg_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="PMA"]
base_d0_tr1_ifng_s_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IFNg_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="iRBCs"]

# sig_fc_corr_purr[sig_fc_corr_purr$stim=="PMA"&sig_fc_corr_purr$gate=="Tr1_IFNg_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="S",]

base_d0_tr1_ifng_s_pma_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("Tr1_IFNg_Frequency"), stim=="PMA", targetName%in%base_d0_tr1_ifng_s_pma_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d0_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IFNg+ Tr1 frequency from baseline to day 0")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous (trans="log2")+
  facet_wrap(~targetName, nrow=5, scale="free")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d0_tr1_ifng_s_pma_sigs_plot.png", base_d0_tr1_ifng_s_pma_sigs_plot, height = 8, width=12, dpi=444, bg="white")



base_d0_tr1_ifng_s_irbc_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("Tr1_IFNg_Frequency"), stim=="iRBCs", targetName%in%base_d0_tr1_ifng_s_irbc_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d0_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IFNg+ Tr1 frequency from baseline to day 0")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous (trans="log2")+
  facet_wrap(~targetName, nrow=3, scales="free")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d0_tr1_ifng_s_irbc_sigs_plot.png", base_d0_tr1_ifng_s_irbc_sigs_plot, height = 6, width=6, dpi=444, bg="white")



base_d0_tr1_ifng_a_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IFNg_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="PMA"]
base_d0_tr1_ifng_a_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IFNg_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="iRBCs"]



base_d0_tr1_ifng_a_pma_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("Tr1_IFNg_Frequency"), stim=="PMA", targetName%in%base_d0_tr1_ifng_a_irbc_sigs, infectiontype=="A")%>%
  ggplot(., aes(x=freq_base_d0_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Asymptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IFNg+ Tr1 frequency from baseline to day 0")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous(trans="log2")+
  facet_wrap(~targetName, nrow=2, scales="free")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d0_tr1_ifng_a_pma_sigs_plot.png", base_d0_tr1_ifng_a_pma_sigs_plot, height = 6, width=6, dpi=444, bg="white")

base_d0_tr1_ifng_a_irbc_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("Tr1_IFNg_Frequency"), stim=="iRBCs", targetName%in%base_d0_tr1_ifng_a_irbc_sigs, infectiontype=="A")%>%
  ggplot(., aes(x=freq_base_d0_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Asymptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IFNg+ Tr1 frequency from baseline to day 0")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous(trans="log2")+
  facet_wrap(~targetName, nrow=2, scales="free")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d0_tr1_ifng_a_irbc_sigs_plot.png", base_d0_tr1_ifng_a_irbc_sigs_plot, height = 6, width=6, dpi=444, bg="white")



base_d0_tr1_ifng_s_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IFNg_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="PMA"]
base_d0_tr1_ifng_s_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IFNg_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="iRBCs"]

# sig_fc_corr_purr[sig_fc_corr_purr$stim=="PMA"&sig_fc_corr_purr$gate=="Tr1_IFNg_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="S",]

base_d0_tr1_ifng_s_pma_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("Tr1_IFNg_Frequency"), stim=="PMA", targetName%in%base_d0_tr1_ifng_s_pma_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d0_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IFNg+ Tr1 frequency from baseline to day 0")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous(trans="log2")+
  facet_wrap(~targetName, nrow=5, scales = "free")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d0_tr1_ifng_s_pma_sigs_plot.png", base_d0_tr1_ifng_s_pma_sigs_plot, height = 16, width=8, dpi=444, bg="white")



base_d0_tr1_ifng_s_irbc_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("Tr1_IFNg_Frequency"), stim=="iRBCs", targetName%in%base_d0_tr1_ifng_s_irbc_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d0_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IFNg+ Tr1 frequency from baseline to day 0")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous(trans="log2")+
  facet_wrap(~targetName, nrow=3)+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d0_tr1_ifng_s_irbc_sigs_plot.png", base_d0_tr1_ifng_s_irbc_sigs_plot, height = 6, width=6, dpi=444, bg="white")



base_d0_tr1_ifng_a_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IFNg_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="PMA"]
base_d0_tr1_ifng_a_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IFNg_Frequency"&sig_fc_corr_purr$Var=="base0"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="iRBCs"]


base_d0_tr1_ifng_a_irbc_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("Tr1_IFNg_Frequency"), stim=="iRBCs", targetName%in%base_d0_tr1_ifng_a_irbc_sigs, infectiontype=="A")%>%
  ggplot(., aes(x=freq_base_d0_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Asymptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IFNg+ Tr1 frequency from baseline to day 0")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous(trans="log2")+
  facet_wrap(~targetName, nrow=2, scales="free")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d0_tr1_ifng_a_irbc_sigs_plot.png", base_d0_tr1_ifng_a_irbc_sigs_plot, height = 6, width=6, dpi=444, bg="white")



## day 14 fold change plots ####

### Tr1 frequency ####
#analytes where FC frome base to day0 correlates with FC from base to day0 in Tr1 frequency
base_d14_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$stim=="unstim"&sig_fc_corr_purr$gate=="Tr1_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="S"]

base_d14_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("Tr1_Frequency"), stim=="unstim", targetName%in%base_d14_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d14_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in Tr1 frequency from baseline to day 14")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  facet_wrap(~targetName)+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d14_sigs_plot.png", base_d14_sigs_plot, height = 6, width=6, dpi=444, bg="white")


### il10 frequency ####
base_d14_il10_s_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IL10_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="PMA"]
base_d14_il10_s_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IL10_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="iRBCs"]

# sig_fc_corr_purr[sig_fc_corr_purr$stim=="PMA"&sig_fc_corr_purr$gate=="IL10_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="S",]

base_d14_il10_s_pma_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("IL10_Frequency"), stim=="iRBCs", targetName%in%base_d14_il10_s_irbc_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d14_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IL10+ CD4 memory frequency from baseline to day 14")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values = stim_palette)+
  scale_x_continuous(trans="log2")+
  facet_wrap(~targetName, scales = "free")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d14_il10_s_irbc_sigs_plot.png", base_d14_il10_s_pma_sigs_plot, height = 16, width=8, dpi=444, bg="white")



base_d14_il10_a_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IL10_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="PMA"]
#nothing
base_d14_il10_a_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IL10_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="iRBCs"]

# both nothing
# base_d14_il10_a_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IL10_Frequency"&sig_fc_corr_purr$Var=="base14freq"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="PMA"]
# base_d14_il10_a_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IL10_Frequency"&sig_fc_corr_purr$Var=="base14freq"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="iRBCs"]


base_d14_il10_a_irbc_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("IL10_Frequency"), stim=="PMA", targetName%in%base_d14_il10_a_pma_sigs, infectiontype=="A")%>%
  ggplot(., aes(x=freq_base_d14_fc, y=conc_base_d14_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Asymptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IL10+ CD4 memory frequency from baseline to day 14")+
  ylab("fold change in plasma protein from baseline to day 14")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous(trans="log2")+
  facet_wrap(~targetName, nrow=4, scales="free")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d14_il10_a_pma_sigs_plot.png", base_d14_il10_a_irbc_sigs_plot, height = 8, width=12, dpi=444, bg="white")



### tr1_il10 ####

#nothing
base_d14_tr1_il10_s_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IL10_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="PMA"]
base_d14_tr1_il10_s_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IL10_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="iRBCs"]

# sig_fc_corr_purr[sig_fc_corr_purr$stim=="PMA"&sig_fc_corr_purr$gate=="Tr1_IL10_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="S",]

base_d14_tr1_il10_s_irbc_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("Tr1_IL10_Frequency"), stim=="iRBCs", targetName%in%base_d14_tr1_il10_s_irbc_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d14_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IL10+ Tr1 frequency from baseline to day 14")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  facet_wrap(~targetName, scales="free")+
  scale_x_continuous(trans="log2")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d14_tr1_il10_s_irbc_sigs_plot.png", base_d14_tr1_il10_s_irbc_sigs_plot, height = 16, width=8, dpi=444, bg="white")



base_d14_tr1_il10_a_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IL10_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="PMA"]
#nothing
base_d14_tr1_il10_a_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IL10_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="iRBCs"]


base_d14_tr1_il10_a_pma_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("Tr1_IL10_Frequency"), stim=="PMA", targetName%in%base_d14_tr1_il10_a_pma_sigs[1:18], infectiontype=="A")%>%
  ggplot(., aes(x=freq_base_d14_fc, y=conc_base_d14_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  # ggtitle("Asymptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IL10+ Tr1 frequency from baseline to day 14")+
  ylab("fold change in plasma protein from baseline to day 14")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous(trans="log2")+
  facet_wrap(~targetName, ncol=9, scales="free")+
  theme_minimal(13)+
  theme(legend.title = element_blank())

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d14_tr1_il10_a_pma_sigs_plot.png", base_d14_tr1_il10_a_pma_sigs_plot, height = 6, width=18, dpi=444, bg="white")

### ifng frequency ####

base_d14_ifng_s_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IFNg_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="PMA"]
base_d14_ifng_s_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IFNg_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="iRBCs"]

# sig_fc_corr_purr[sig_fc_corr_purr$stim=="PMA"&sig_fc_corr_purr$gate=="IFNg_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="S",]


base_d14_ifng_s_irbc_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("IFNg_Frequency"), stim=="iRBCs", targetName%in%base_d14_ifng_s_irbc_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d14_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IFNg+ CD4 memory frequency from baseline to day 14")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  facet_wrap(~targetName, ncol=3, scales="free")+
  scale_x_continuous(trans="log2")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d14_ifng_s_irbc_sigs_plot.png", base_d14_ifng_s_irbc_sigs_plot, height = 6, width=8, dpi=444, bg="white")



base_d14_ifng_a_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IFNg_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="PMA"]
#nothing
base_d14_ifng_a_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IFNg_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="iRBCs"]


base_d14_ifng_a_pma_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("IFNg_Frequency"), stim=="PMA", targetName%in%base_d14_ifng_a_pma_sigs, infectiontype=="A")%>%
  ggplot(., aes(x=freq_base_d14_fc, y=conc_base_d14_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Asymptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IFNg+ CD4 memory frequency from baseline to day 14")+
  ylab("fold change in plasma protein from baseline to day day 14")+
  scale_color_manual(values =stim_palette)+
  facet_wrap(~targetName, ncol=3, scales="free")+
  scale_x_continuous(trans="log2")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d14_ifng_a_pma_sigs_plot.png", base_d14_ifng_a_pma_sigs_plot, height = 6, width=8, dpi=444, bg="white")

### Tr1 IFNg frequency ####

#nothing both
base_d14_tr1_ifng_s_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IFNg_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="PMA"]
base_d14_tr1_ifng_s_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IFNg_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="iRBCs"]


base_d14_tr1_ifng_a_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IFNg_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="PMA"]
#nothing
base_d14_tr1_ifng_a_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IFNg_Frequency"&sig_fc_corr_purr$Var=="base14"&sig_fc_corr_purr$infectiontype=="A"&sig_fc_corr_purr$stim=="iRBCs"]


base_d14_tr1_ifng_a_pma_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("Tr1_IFNg_Frequency"), stim=="PMA", targetName%in%base_d14_tr1_ifng_a_pma_sigs, infectiontype=="A")%>%
  ggplot(., aes(x=freq_base_d14_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Asymptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IFNg+ Tr1 frequency from baseline to day 14")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous(trans="log2", limits = c(0.5, 2))+
  facet_wrap(~targetName, scales="free")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d14_tr1_ifng_a_pma_sigs_plot.png", base_d14_tr1_ifng_a_pma_sigs_plot, height = 8, width=10, dpi=444, bg="white")




# day 7 fold changes ####
### Tr1 frequency ####
#analytes where FC frome base to day0 correlates with FC from base to day0 in Tr1 frequency
base_d7_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$stim=="unstim"&sig_fc_corr_purr$gate=="Tr1_Frequency"&sig_fc_corr_purr$Var=="base7"&sig_fc_corr_purr$infectiontype=="S"]

base_d7_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("Tr1_Frequency"), stim=="unstim", targetName%in%base_d7_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d7_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in Tr1 frequency from baseline to day 7")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  facet_wrap(~targetName)+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d7_sigs_plot.png", base_d7_sigs_plot, height = 6, width=6, dpi=444, bg="white")


### il10 frequency ####
base_d7_il10_s_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IL10_Frequency"&sig_fc_corr_purr$Var=="base7"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="PMA"]
base_d7_il10_s_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IL10_Frequency"&sig_fc_corr_purr$Var=="base7"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="iRBCs"]

# sig_fc_corr_purr[sig_fc_corr_purr$stim=="PMA"&sig_fc_corr_purr$gate=="IL10_Frequency"&sig_fc_corr_purr$Var=="base7"&sig_fc_corr_purr$infectiontype=="S",]

base_d7_il10_s_pma_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("IL10_Frequency"), stim=="PMA", targetName%in%base_d7_il10_s_pma_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d7_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IL10+ CD4 memory frequency from baseline to day 7")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  facet_wrap(~targetName, nrow=2, scales = "free")+
  scale_x_continuous(trans="log2")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d7_il10_s_pma_sigs_plot.png", base_d7_il10_s_pma_sigs_plot, height = 8, width=8, dpi=444, bg="white")



base_d7_il10_s_irbc_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("IL10_Frequency"), stim=="iRBCs", targetName%in%base_d7_il10_s_irbc_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d7_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IL10+ CD4 memory frequency from baseline to day 7")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  facet_wrap(~targetName, nrow=3)+
  scale_x_continuous(trans="log2")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d7_il10_s_irbc_sigs_plot.png", base_d7_il10_s_irbc_sigs_plot, height = 8, width=6, dpi=444, bg="white")



### tr1_il10 ####

#nothing
base_d7_tr1_il10_s_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IL10_Frequency"&sig_fc_corr_purr$Var=="base7"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="PMA"]
base_d7_tr1_il10_s_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IL10_Frequency"&sig_fc_corr_purr$Var=="base7"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="iRBCs"]

### ifng frequency ####

base_d7_ifng_s_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IFNg_Frequency"&sig_fc_corr_purr$Var=="base7"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="PMA"]
base_d7_ifng_s_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="IFNg_Frequency"&sig_fc_corr_purr$Var=="base7"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="iRBCs"]

# sig_fc_corr_purr[sig_fc_corr_purr$stim=="PMA"&sig_fc_corr_purr$gate=="IFNg_Frequency"&sig_fc_corr_purr$Var=="base7"&sig_fc_corr_purr$infectiontype=="S",]


base_d7_ifng_s_pma_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("IFNg_Frequency"), stim=="PMA", targetName%in%base_d7_ifng_s_pma_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d7_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IFNg+ CD4 memory frequency from baseline to day 7")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  facet_wrap(~targetName, nrow=3)+
  scale_x_continuous(trans="log2")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d7_ifng_s_pma_sigs_plot.png", base_d7_ifng_s_pma_sigs_plot, height = 8, width=6, dpi=444, bg="white")


base_d7_ifng_s_irbc_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("IFNg_Frequency"), stim=="iRBCs", targetName%in%base_d7_ifng_s_irbc_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d7_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IFNg+ CD4 memory frequency from baseline to day 7")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  facet_wrap(~targetName, nrow=3)+
  scale_x_continuous(trans="log2")+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d7_ifng_s_irbc_sigs_plot.png", base_d7_ifng_s_irbc_sigs_plot, height = 8, width=6, dpi=444, bg="white")



### Tr1 IFNg frequency ####
base_d7_tr1_ifng_s_pma_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IFNg_Frequency"&sig_fc_corr_purr$Var=="base7"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="PMA"]
base_d7_tr1_ifng_s_irbc_sigs <- sig_fc_corr_purr$targetName[sig_fc_corr_purr$gate=="Tr1_IFNg_Frequency"&sig_fc_corr_purr$Var=="base7"&sig_fc_corr_purr$infectiontype=="S"&sig_fc_corr_purr$stim=="iRBCs"]


base_d7_tr1_ifng_s_pma_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("Tr1_IFNg_Frequency"), stim=="PMA", targetName%in%base_d7_tr1_ifng_s_pma_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d7_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IFNg+ Tr1 frequency from baseline to day 7")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  facet_wrap(~targetName, nrow=1)+
  scale_x_continuous()+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d7_tr1_ifng_s_pma_sigs_plot.png", base_d7_tr1_ifng_s_pma_sigs_plot, height = 16, width=8, dpi=444, bg="white")


base_d7_tr1_ifng_s_irbc_sigs_plot <- fc_combo_frame %>%
  filter(gate %in% c("Tr1_IFNg_Frequency"), stim=="iRBCs", targetName%in%base_d7_tr1_ifng_s_irbc_sigs, infectiontype=="S")%>%
  ggplot(., aes(x=freq_base_d7_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("Symptomatic")+
  geom_smooth(method="lm")+
  xlab("fold change in IFNg+ Tr1 frequency from baseline to day 7")+
  ylab("fold change in plasma protein from baseline to day 0")+
  scale_color_manual(values =stim_palette)+
  scale_x_continuous(trans="log2")+
  facet_wrap(~targetName, nrow=1)+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_d7_tr1_ifng_s_irbc_sigs_plot.png", base_d7_tr1_ifng_s_irbc_sigs_plot, height = 16, width=8, dpi=444, bg="white")



# emmeans ####
library(emmeans)

combo_as_results <- nulisa_data %>%
  filter(infectiontype%in%c("A", "S"), timepoint %notin% c("day28"), timepoint!="bad_baseline")%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  group_by(targetName)%>%
  nest() %>%
  #mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id_cat), data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id_cat), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  # mutate(emm=map(model, ~emmeans(., specs = ~ timepoint:infectiontype)))%>%
  mutate(emm=map(model, ~emmeans(., specs = pairwise ~ timepoint | infectiontype)))%>%
  mutate(emm2=map(model, ~emmeans(., specs = pairwise ~ infectiontype | timepoint)))%>%
  
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise")))%>%
  
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  
  mutate("baseline A - baseline S"=map_dbl(emm_contrast_summary2, ~.$p.value[1])) %>%
  mutate("baseline A - day0 A"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("baseline A - day14 A"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  mutate("baseline S - day0 S"=map_dbl(emm_contrast_summary, ~.$p.value[7])) %>%
  mutate("baseline S - day7 S"=map_dbl(emm_contrast_summary, ~.$p.value[8])) %>%
  mutate("baseline S - day14 S"=map_dbl(emm_contrast_summary, ~.$p.value[9]))%>%
  mutate("coef baseline A - baseline S"=map_dbl(emm_contrast_summary2, ~.$estimate[1])) %>%
  mutate("coef baseline A - day0 A"=map_dbl(emm_contrast_summary, ~.$estimate[1])) %>%
  mutate("coef baseline A - day14 A"=map_dbl(emm_contrast_summary, ~.$estimate[3])) %>%
  mutate("coef baseline S - day0 S"=map_dbl(emm_contrast_summary, ~.$estimate[7])) %>%
  mutate("coef baseline S - day7 S"=map_dbl(emm_contrast_summary, ~.$estimate[8])) %>%
  mutate("coef baseline S - day14 S"=map_dbl(emm_contrast_summary, ~.$estimate[9]))%>%
  pivot_longer(cols=starts_with("baseline"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))%>%
  ungroup()%>%
  pivot_longer(cols=starts_with("coef"), names_to = "coef_name", values_to = "coef")%>%
  rowwise()%>%
  filter(grepl(contrast, coef_name))

fc_data <- nulisa_data %>%
  pivot_wider(names_from = timepoint, values_from = c(concentration), id_cols=c(targetName, id, infectiontype), names_prefix = "conc_")%>%
  group_by(infectiontype)%>%
  mutate(conc_base_d0_fc=conc_day0-conc_baseline,
         conc_base_d14_fc=conc_day14-conc_baseline)%>%
  distinct(targetName, id, infectiontype, conc_base_d0_fc, conc_base_d14_fc)%>%
  group_by(targetName, infectiontype)%>%
  summarise("mean_conc_base_d0_fc"=mean(conc_base_d0_fc, na.rm = T), "mean_conc_base_d14_fc"=mean(conc_base_d14_fc, na.rm = T))

analytes_of_interest <- c("IL10", "LAG3", "GZMA", "IL27", "IL6", "TNF", "IFNG", "CTLA4", "LILRB2", "CRP", "CCL24", "IL7R", "CX3CL1", "KDR", "TNFSF11")

(base_day0_s_volcano <- combo_as_results %>%
    filter(contrast == "baseline S - day0 S")%>%
    left_join(., fc_data, by="targetName")%>%
    distinct(targetName, coef, padj, infectiontype, mean_conc_base_d0_fc)%>%
    filter(infectiontype=="S", is.finite(mean_conc_base_d0_fc))%>%
    mutate("label2" = if_else(targetName %in% analytes_of_interest, targetName, NA))%>%
    ggplot(., aes(x=mean_conc_base_d0_fc, y=-log10(padj+10^-14), alpha=padj<0.05, color=coef>0))+
    geom_point()+
    # geom_label(aes(label=label2), nudge_x = 0.1, )+
    ggrepel::geom_text_repel(aes(label=label2), nudge_x = 0)+
    ggtitle("n = 145")+
    geom_hline(yintercept = -log10(0.05+10^-14), linetype="dashed")+
    scale_alpha_manual(values=c(0.5, 1))+
    scale_color_manual(values=c("darkred", "darkblue"))+
    scale_x_continuous(limits=c(-5, 5))+
    xlab("log2 fold change")+
    ylab("-log10 q value")+
    theme_minimal()+
    theme(legend.position="none"))

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/revised_baselines/base_day0_s_volcano.png", base_day0_s_volcano, width = 5, height = 5, dpi=444, bg="white")
# ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/revised_baselines/base_day0_s_volcano_no_label.png", base_day0_s_volcano, width = 5, height = 5, dpi=444, bg="white")


# figures for jason's paper ####
tr1_genes <- c("GZMA", "IL10", "CTLA4", "LAG3", "IL7R", "CXCL6")

jason_box <- clean_data%>%
    filter(targetName %in% tr1_genes)%>%
    filter(timepoint%notin%c("day28", "bad_baseline"))%>%
    filter(infectiontype %in% c("A", "S"))%>%
    ggplot(., aes(x=infectiontype, y=concentration, fill = timepoint))+
    geom_boxplot(outliers = F)+
    theme_minimal()+
    facet_wrap(~targetName, scales = "free")+
    scale_fill_manual(values=time_cols)+
    da_boxplot_theme+
  theme(legend.position = "right",
        legend.title = element_blank())

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/musical_tr1_genes_boxplot.png", jason_box, width =8, height=6, dpi=444, bg="white")


tr1_freq_plasma_protein_cor <- long_combo %>%
  filter(timepoint=="baseline", infectiontype %in% c("A", "S"), stim!="unstim")%>%
  filter(gate %in% c("Tr1_IL10_Frequency", "IL10_Frequency"),targetName %in% c("IL10", "LAG3", "GZMA"))%>%
  ggplot(., aes(x=freq, y=concentration, color=stim))+
  geom_point()+
  facet_wrap(~targetName+gate+stim, scales="free", ncol=4)+
  geom_smooth(method="lm")+
  ylab("plasma concentration at baseline")+
  xlab("% Cytokine+ Cells")+
  scale_y_continuous(limits = c(NA, 18))+
  ggpubr::stat_cor(method="spearman", label.y = 17.2)+
  scale_color_manual(values=stim_palette)+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/tr1_freq_plasma_protein_cor.png", tr1_freq_plasma_protein_cor, width = 9, height=9, dpi=444, bg="white")




tr1_freq_base_plasma_protein_cor <- long_combo %>%
  filter(timepoint=="baseline", infectiontype %in% c("A", "S"), stim=="unstim")%>%
  filter(gate %in% c("Tr1_Frequency"),targetName %in% c("IL10", "LAG3", "GZMA"))%>%
  ggplot(., aes(x=freq, y=concentration, color=stim))+
  geom_point()+
  facet_wrap(~gate+targetName+stim, scales="free")+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman", label.y = 17.5)+
  ylab("IL10 plasma concentration at baseline")+
  xlab("Frequency of Tr1 cells")+
  scale_y_continuous(limits = c(NA, 18))+
  scale_color_manual(values=stim_palette)+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/tr1_freq_base_plasma_protein_cor.png", tr1_freq_base_plasma_protein_cor, width = 8, height=3, dpi=444, bg="white")

base_freq_day0_conc <- inner_join(day0_concs, base_freqs, by=c("id", "infectiontype"), relationship = "many-to-many")

 base_freq_day0_conc %>%
  filter(infectiontype%in%c("A", "S"))%>%
  filter(gate %in% c("Tr1_IL10_Frequency", "IL10_Frequency"),targetName %in% c("IL10", "LAG3", "GZMA"), stim!="unstim")%>%
  ggplot(., aes(x=freq, y=concentration, color=stim))+
  geom_point()+
  facet_wrap(~targetName+gate+stim, scales="free", ncol=4)+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman", label.y = 17.5)+
  ggtitle("All Infections")+
  ylab("plasma concentration at day0")+
  xlab("Frequency of Tr1 cells at baseline")+
  scale_y_continuous(limits = c(NA, 18))+
  scale_color_manual(values=stim_palette)+
  theme_minimal()

 ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/all_cell_freq_base_day0_conc_tr1_genes.png", width = 9, height=9, dpi=444, bg="white")
 
 base_freq_day0_conc %>%
 filter(infectiontype %in% c("S"), stim=="unstim")%>%
   filter(gate %in% c("Tr1_Frequency"),targetName %in% c("IL10", "LAG3", "GZMA"))%>%
   ggplot(., aes(x=freq, y=concentration, color=stim))+
   geom_point()+
   facet_wrap(~gate+targetName+stim, scales="free")+
   geom_smooth(method="lm")+
   ggpubr::stat_cor(method="spearman", label.y = 17.5)+
   ylab("IL10 plasma concentration at baseline")+
   xlab("Frequency of Tr1 cells")+
   scale_y_continuous(limits = c(NA, 18))+
   scale_color_manual(values=stim_palette)+
   theme_minimal()
 
 ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/all_tr1_freq_base_day0_conc.png", width = 9, height=4.5, dpi=444, bg="white")

 
 individuals_with_cleaner_baseline <- nulisa_data%>%
   filter(timepoint=="bad_baseline")%>%
   distinct(id)
   
 
 nulisa_data%>%
   filter(id %in% individuals_with_cleaner_baseline$id)%>%
   distinct(id, timepoint, log_qpcr, parasitedensity, infectiontype, new_qpcr)%>%
   mutate(timepoint2=if_else(timepoint=="bad_baseline", "cleaner baseline", "most recent baseline"))%>%
   filter(infectiontype%in%c("S", "A"), timepoint%in%c("bad_baseline", "baseline"))%>%
   ggplot(., aes(x=timepoint2, y=new_qpcr+0.001, fill=timepoint2))+
   geom_label(aes(label=id))+
   scale_y_log10()+
   ylab("qPCR parasites / μL")+
   xlab("")+
   geom_point(position = position_jitter(width = 0.2))+
   facet_wrap(~infectiontype, scales="free")+
   ggpubr::stat_cor(method="spearman", label.y = 17.5)+
   theme_minimal()+
   theme(legend.position = "none")

 
 