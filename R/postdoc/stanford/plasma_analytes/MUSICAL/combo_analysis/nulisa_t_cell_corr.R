#preamble ####
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)

stim_palette <- c("darkred", "darkblue", "black")
names(stim_palette) <- c("iRBCs", "PMA", "unstim")

`%notin%` <- Negate(`%in%`)

# data generation ####
nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")
slim_nulisa_data <- nulisa_data %>%
  mutate(id=as.character(id))%>%
  select(id, timepoint, timepoint_imm, infectiontype, targetName, concentration)

# cell_count_data <- read.csv("~/postdoc/stanford/manuscripts/jason_tr1_2/t_cell_cytometry_count_data.csv")
cell_count_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/df_jason_analysis.csv")

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
  select(id, timepoint, timepoint_imm, infectiontype, gate, stim, freq)


long_combo <- slim_nulisa_data%>%
  mutate(id=as.integer(id))%>%
  full_join(., slim_cell_count_data, by = c("id", "timepoint_imm", "infectiontype"))%>%
  mutate(timepoint=timepoint.x)%>%
  select(-timepoint.x, -timepoint.y)
#   
# long_combo%>%
#   filter(infectiontype=="S",
#          timepoint_imm==7,
#          stim=="unstim",
#          gate=="Tr1_Frequency",
#          targetName=="AGER")


# disjointed timepoints, absolute values ####
day0_concs <- long_combo %>%
  filter(infectiontype%in% c("S", "A"))%>%
  filter(timepoint=="day0")%>%
  distinct(targetName, concentration, id, timepoint, infectiontype)
  
base_concs <- long_combo %>%
  filter(infectiontype%in% c("S", "A"))%>%
  filter(timepoint=="baseline")%>%
  distinct(targetName, concentration, id, timepoint, infectiontype)

base_freqs <- long_combo%>%
  filter(infectiontype%in% c("S", "A"))%>%
  filter(timepoint=="baseline")%>%
  distinct(gate, freq, stim, id, timepoint, infectiontype)

day0_freqs <- long_combo%>%
  filter(infectiontype%in% c("S", "A"))%>%
  filter(timepoint=="day0")%>%
  distinct(gate, freq, stim, id, timepoint, infectiontype)

day7_freqs <- long_combo%>%
  filter(infectiontype%in% c("S", "A"))%>%
  filter(infectiontype=="S", timepoint=="day7")%>%
  distinct(gate, freq, stim, id, timepoint, infectiontype)

day14_freqs <- long_combo%>%
  filter(infectiontype%in% c("S", "A"))%>%
  filter(timepoint=="day14")%>%
  distinct(gate, freq, stim, id, timepoint, infectiontype)


# 
# base_counts <- long_combo%>%
#   filter(infectiontype%in% c("S", "A"))%>%
#   filter(timepoint=="baseline")%>%
#   distinct(gate, count, stim, id, timepoint, infectiontype)
# 
# day0_counts <- long_combo%>%
#   filter(infectiontype%in% c("S", "A"))%>%
#   filter(timepoint=="day0")%>%
#   distinct(gate, count, stim, id, timepoint, infectiontype)
# 
# day7_counts <- long_combo%>%
#   filter(infectiontype%in% c("S", "A"))%>%
#   filter(infectiontype=="S", timepoint=="day7")%>%
#   distinct(gate, count, stim, id, timepoint, infectiontype)
# 
# day14_counts <- long_combo%>%
#   filter(infectiontype%in% c("S", "A"))%>%
#   filter(timepoint=="day14")%>%
#   distinct(gate, count, stim, id, timepoint, infectiontype)
# 


base_freq_base_conc <- inner_join(base_concs, base_freqs, by=c("id", "infectiontype"), relationship = "many-to-many")

base_freq_day0_conc <- inner_join(day0_concs, base_freqs, by=c("id", "infectiontype"), relationship = "many-to-many")

conc_freq_0 <- inner_join(day0_concs, day0_freqs, by=c("id", "infectiontype"), relationship = "many-to-many")

conc_freq_014 <- inner_join(day0_concs, day14_freqs, by=c("id", "infectiontype"), relationship = "many-to-many")



# base_count_base_conc <- inner_join(base_concs, base_counts, by=c("id", "infectiontype"), relationship = "many-to-many")
# 
# base_count_day0_conc <- inner_join(day0_concs, base_counts, by=c("id", "infectiontype"), relationship = "many-to-many")
# 
# conc_count_0 <- inner_join(day0_concs, day0_counts, by=c("id", "infectiontype"), relationship = "many-to-many")
# 
# conc_count_014 <- inner_join(day0_concs, day14_counts, by=c("id", "infectiontype"), relationship = "many-to-many")





conc_freq_07 <- day0_concs%>%
  filter(infectiontype=="S")%>%
  inner_join(day7_freqs, by=c("id", "infectiontype"), relationship = "many-to-many")
# 
# conc_count_07 <- day0_concs%>%
#   filter(infectiontype=="S")%>%
#   inner_join(day7_counts, by=c("id", "infectiontype"), relationship = "many-to-many")

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


conc_freq_07_corr <- conc_freq_07%>%
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


# 
# 
# conc_count_0_corr <- conc_count_0%>%
#   filter(!is.na(count))%>%
#   group_by(gate, stim, targetName, infectiontype)%>%
#   nest()%>% 
#   mutate(nrows=map_dbl(data, ~nrow(.)))%>%
#   filter(nrows>1)%>%
#   mutate(correlation=map(data, ~cor.test(.$concentration, .$count, method = "spearman")))%>%
#   mutate(p=map_dbl(correlation, ~.$p.value),
#          rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$count, method="spearman")))%>%
#   ungroup()%>%
#   group_by(infectiontype)%>%
#   mutate(padj=p.adjust(p))
# 
# 
# 
# conc_count_07_corr <- conc_count_07%>%
#   filter(!is.na(count))%>%
#   group_by(gate, stim, targetName)%>%
#   nest()%>%
#   mutate(nrows=map_dbl(data, ~nrow(.)))%>%
#   filter(nrows>1)%>%
#   mutate(correlation=map(data, ~cor.test(.$concentration, .$count, method = "spearman")))%>%
#   mutate(p=map_dbl(correlation, ~.$p.value),
#          rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$count, method="spearman")))%>%
#   ungroup()%>%
#   mutate(padj=p.adjust(p))
# 



remove_from_conc_freq_014 <- conc_freq_014 %>%
  filter(gate == "CD4_Lymphocyte_Frequency", stim == "PMA", infectiontype == "A")

conc_freq_014_corr <- conc_freq_014%>%
  filter(gate %notin% c("Memory_CD4_T_Cell_Frequency",
                      "CD4_Lymphocyte_Frequency",
                      unique(grep("FOXP3", gate, value=TRUE))))%>%
  filter(!is.na(gate))%>%
  group_by(gate, stim, targetName, infectiontype)%>%
  nest()%>%
  mutate(nrows=map_dbl(data, ~nrow(.)))%>%
  filter(nrows>1)%>%
  mutate(correlation=map(data, ~cor.test(.$concentration, .$freq, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  ungroup()%>%
  group_by()%>%
  mutate(padj=p.adjust(p))


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

base_freq_day0_conc_cor <- base_freq_day0_conc%>%
  filter(!is.na(concentration), !is.na(freq))%>%
  filter(gate %notin% c("Memory_CD4_T_Cell_Frequency",
                        "CD4_Lymphocyte_Frequency",
                        unique(grep("FOXP3", gate, value=TRUE))))%>%
  group_by(gate, stim, infectiontype)%>%
  nest()%>%
  mutate(nrows=map_dbl(data, ~nrow(.)))%>%
  filter(nrows>1)%>%
  mutate(correlation=map(data, ~cor.test(.$concentration, .$freq, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  ungroup()%>%
  group_by()%>%
  mutate(padj=p.adjust(p))

base_freq_base_conc_cor <- base_freq_base_conc%>%
  filter(!is.na(concentration), !is.na(freq))%>%
  filter(gate %notin% c("Memory_CD4_T_Cell_Frequency",
                        "CD4_Lymphocyte_Frequency",
                        unique(grep("FOXP3", gate, value=TRUE))))%>%
  group_by(gate, stim, targetName)%>%
  nest()%>%
  mutate(nrows=map_dbl(data, ~nrow(.)))%>%
  filter(nrows>1)%>%
  mutate(correlation=map(data, ~cor.test(.$concentration, .$freq, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%
  ungroup()%>%
  mutate(padj=p.adjust(p))


fdr_cutoff=0.05


sig_conc_freq_0_corr <- conc_freq_0_corr%>%
  filter(padj<fdr_cutoff)

sig_conc_freq_07_corr <- conc_freq_07_corr%>%
  filter(padj<fdr_cutoff)

sig_conc_freq_014_corr <- conc_freq_014_corr%>%
  filter(padj<fdr_cutoff)

sig_freq_freq_base_0 <- freq_freq_base_0_corr%>%
  filter(padj<fdr_cutoff)

sig_base_freq_day0_conc_cor <- base_freq_day0_conc_cor%>%
  filter(padj<fdr_cutoff)


sig_conc_count_0_corr <-  conc_count_0_corr%>%
  filter(padj<fdr_cutoff)

sig_conc_count_07_corr <- conc_count_07_corr%>%
  filter(padj<fdr_cutoff)


sig_base_freq_base_conc_cor <- base_freq_base_conc_cor%>%
  filter(padj<fdr_cutoff)
# 
# for_jason <- base_freq_base_conc_cor%>%
#   select(-data, -correlation)
# 
# write.csv(for_jason, "~/Downloads/base_freq_base_conc_cor.csv", row.names = F)


sig_conc_freq_07_corr_plot <- conc_freq_07_corr %>%
  # filter(targetName %in% c("GZMA", "LAG3", "IL10", "CTLA4", "LILRB2", "IL6"), gate %in% c("Tr1_Frequency", "IL10_Frequency"))%>%
  filter(targetName %in% "GZMA", gate %in% c("IL10_Frequency"))%>%
  filter(stim %in% c("unstim", "iRBCs"))%>%
  unnest(data)%>%
  mutate(gate=gsub("_", " ", gate, fixed=TRUE))%>%
  ggplot(., aes(x=freq, y=concentration, color=stim))+
  geom_point()+
  facet_wrap(~gate+stim+targetName, ncol=4, scales="free")+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  ggtitle("GMZMA concentration at day 0 correlates with\n IL10 secretion cells at day7")+
  geom_smooth(method="lm")+
  xlab("Concentration")+
  ylab("Percentage")+
  scale_color_manual(values=stim_palette)+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/GZMA_IL10_corr.png", sig_conc_freq_07_corr_plot, width=5.3333, height=5.333, bg="white")



cxcl11_angpt2_treg_corr_plot <- conc_freq_014_corr %>%
  filter(targetName %in% c("CXCL11", "ANGPT2"), gate %in% c("Tregs_CD25_Frequency"), infectiontype=="S")%>%
  unnest(data)%>%
  filter(stim =="unstim")%>%
  ggplot(., aes(x=freq, y=concentration))+
  geom_point()+
  facet_wrap(~gate+targetName, scales="free")+
  ggpubr::stat_cor(method="spearman")+
  ggtitle("plasma concentration of CXCL11 and ANGPT2 at day0 correlates with\nCD25 expression on Tregs at day14")+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("darkred", "darkblue", "white"))+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/cxcl11_angpt2_treg_corr.png", cxcl11_angpt2_treg_corr_plot, width=5.3333, height=5.333, bg="white")



TNFSF10_tr1_il21_corr_plot <- conc_freq_014_corr %>%
  filter(targetName %in% c("TNFSF10"), gate %in% c("Tr1_IL21_Frequency"), infectiontype=="A")%>%
  unnest(data)%>%
  filter(stim =="PMA")%>%
  ggplot(., aes(x=freq, y=concentration))+
  geom_point()+
  facet_wrap(~gate+targetName, scales="free")+
  ggpubr::stat_cor(method="spearman")+
  ggtitle("plasma concentration of TNSF10 day0 correlates with\nIL12 expression by Tr1 at day14")+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("darkred", "darkblue", "white"))+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/cxcl11_angpt2_treg_corr.png", cxcl11_angpt2_treg_corr_plot, width=5.3333, height=5.333, bg="white")


# freq_freq_base_0 %>%
#   filter(gate %in% sig_freq_freq_base_0$gate)%>%
#   distinct(freq.x, freq.y, stim, gate, id, timepoint.x, infectiontype)%>%
#   ggplot(., aes(x=freq.x, y=freq.y, color=stim))+
#   geom_point()+
#   geom_smooth()+
#   facet_wrap(~gate+infectiontype)+
#   ggpubr::stat_cor(method="spearman")+
#   geom_smooth(method="lm")+
#   scale_color_manual(values=c("darkred", "darkblue", "white"))+
#   theme_minimal()
  

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
  mutate(conc_base_d0_fc=conc_day0/conc_baseline,
         conc_base_d14_fc=conc_day14/conc_baseline)%>%
  distinct(targetName, id, infectiontype, conc_base_d0_fc, conc_base_d14_fc)
  

wide_freqs <-  long_combo %>%
  # do some fancy subsetting so that phenotypic gates are only present from unstim samples and cytokine gates only from stim conditions
  # filter(!(stim=="PMA" & !grepl("_[^_]+_", gate)),
  #        !(stim=="iRBCs" & !grepl("_[^_]+_", gate)),
  #        #!(stim=="unstim" & grepl("_[^_]+_", gate)),
  #        #!(stim=="PMA" & gate=="CD4_T_Cell_Frequency")
  #        )%>%
  select(gate, stim, freq, id, infectiontype, timepoint)%>%
  filter(infectiontype!="NM", timepoint %in% c("baseline", "day0", "day7", "day14"))%>%
  distinct(gate, stim, freq, id, infectiontype, timepoint)%>%
  pivot_wider(names_from = timepoint, values_from = c(freq), id_cols=c(gate, stim, id, infectiontype), names_prefix = "freq_")%>%
  group_by(infectiontype)%>%
  mutate(freq_base_d0_fc=freq_day0/freq_baseline,
         freq_base_d7_fc=freq_day7/freq_baseline,
         freq_base_d14_fc=freq_day14/freq_baseline)%>%
  distinct(gate, stim, id, infectiontype, freq_base_d0_fc, freq_base_d7_fc, freq_base_d14_fc)


# wide_counts <-  long_combo %>%
#   select(gate, stim, count, id, infectiontype, timepoint)%>%
#   filter(infectiontype!="NM", timepoint %in% c("baseline", "day0", "day7", "day14"))%>%
#   distinct(gate, stim, count, id, infectiontype, timepoint)%>%
#   pivot_wider(names_from = timepoint, values_from = c(count), id_cols=c(gate, stim, id, infectiontype), names_prefix = "count_")%>%
#   group_by(infectiontype)%>%
#   mutate(count_base_d0_fc=count_day0/count_baseline,
#          count_base_d7_fc=count_day7/count_baseline,
#          count_base_d14_fc=count_day14/count_baseline)%>%
#   distinct(gate, stim, id, infectiontype, count_base_d0_fc, count_base_d7_fc, count_base_d14_fc)
# 


fc_combo_frame <- inner_join(wide_freqs, wide_concs, by=c("id", "infectiontype"))%>%
  filter(is.finite(freq_base_d0_fc),
         is.finite(freq_base_d14_fc),
         is.finite(conc_base_d0_fc),
         is.finite(conc_base_d14_fc))


fc_corr_purr_S <- fc_combo_frame %>%
  filter(targetName!="CTSS")%>%
  filter(infectiontype%in%c("S"), gate%notin%c("Memory_CD4_T_Cell_Frequency",
                                          "CD4_Lymphocyte_Frequency",
                                          "CCR4_P_CXCR6_N_IFNg_Frequency",
                                          "CCR4_P_CXCR6_N_IL21_Frequency",
                                          unique(grep("FOXP3", gate, value=TRUE))))%>%
  # filter(!is.infinite(freq_base_d7_fc), !is.infinite(freq_base_d14_fc), !is.infinite(conc_base_d0_fc))%>%
  # filter(!is.na(freq_base_d7_fc), !is.na(freq_base_d14_fc), !is.na(conc_base_d0_fc))%>%
  group_by(infectiontype, stim, gate, targetName)%>%
  nest()%>%
  # mutate(base_d7_correlation=if_else(all(is.na(data$freq_day7)), map(data, ~cor.test(.$freq_base_d7_fc, .$conc_base_d0_fc, method = "spearman"))))%>%
  mutate(base_d0_correlation=map(data, ~cor.test(.$freq_base_d0_fc, .$conc_base_d0_fc, method = "spearman")))%>%
  mutate(base_d7_correlation=map(data, ~cor.test(.$freq_base_d7_fc, .$conc_base_d0_fc, method = "spearman")))%>%
  mutate(base_d14_correlation=map(data, ~cor.test(.$freq_base_d14_fc, .$conc_base_d0_fc, method = "spearman")))%>%
  mutate(
         base_d0_p=map_dbl(base_d0_correlation, ~.$p.value),
         base_d0_rho=map_dbl(base_d0_correlation, ~.$estimate),
         base_d7_p=map_dbl(base_d7_correlation, ~.$p.value),
         base_d7_rho=map_dbl(base_d7_correlation, ~.$estimate),
         base_d14_p=map_dbl(base_d14_correlation, ~.$p.value),
         base_d14_rho=map_dbl(base_d14_correlation, ~.$estimate))%>%
  group_by(infectiontype, stim, targetName)%>%
  mutate(base_d14_padj=p.adjust(base_d14_p),
         base_d7_padj=p.adjust(base_d7_p),
         base_d0_padj=p.adjust(base_d0_p))
  

fc_corr_purr_A <- fc_combo_frame %>%
  filter(infectiontype%in%c("A"), gate %notin% c("Memory_CD4_T_Cell_Frequency",
                                               "CD4_Lymphocyte_Frequency",
                                               "CCR4_P_CXCR6_N_IFNg_Frequency",
                                               unique(grep("FOXP3", gate, value=TRUE))))%>%
  # filter(!is.infinite(freq_base_d7_fc), !is.infinite(freq_base_d14_fc), !is.infinite(conc_base_d0_fc))%>%
  # filter(!is.na(freq_base_d7_fc), !is.na(freq_base_d14_fc), !is.na(conc_base_d0_fc))%>%
  group_by(infectiontype, stim, gate, targetName)%>%
  nest()%>%
  # mutate(base_d7_correlation=if_else(all(is.na(data$freq_day7)), map(data, ~cor.test(.$freq_base_d7_fc, .$conc_base_d0_fc, method = "spearman"))))%>%
  mutate(base_d0_correlation=map(data, ~cor.test(.$freq_base_d0_fc, .$conc_base_d0_fc, method = "spearman")))%>%
  mutate(base_d14_correlation=map(data, ~cor.test(.$freq_base_d14_fc, .$conc_base_d0_fc, method = "spearman")))%>%
  mutate(base_d14_correlation2=map(data, ~cor.test(.$freq_base_d14_fc, .$conc_base_d14_fc, method = "spearman")))%>%
  mutate(
    base_d0_p=map_dbl(base_d0_correlation, ~.$p.value),
    base_d0_rho=map_dbl(base_d0_correlation, ~.$estimate),
    base_d14_p=map_dbl(base_d14_correlation, ~.$p.value),
    base_d14_rho=map_dbl(base_d14_correlation, ~.$estimate),
    base_d14_p2=map_dbl(base_d14_correlation2, ~.$p.value),
    base_d14_rho2=map_dbl(base_d14_correlation2, ~.$estimate))%>%
  group_by(infectiontype, stim, targetName)%>%
  mutate(base_d14_padj=p.adjust(base_d14_p),
         base_d0_padj=p.adjust(base_d0_p),
         base_d14_padj2=p.adjust(base_d14_p2))# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%

fc_corr_purr <- bind_rows(fc_corr_purr_A, fc_corr_purr_S)



sig_fc_corr_purr <- fc_corr_purr%>%
  filter(base_d14_padj2<0.05 | base_d14_padj < 0.05 | base_d7_padj<0.05 | base_d0_padj < 0.05)
  


tr1_gzma <- fc_combo_frame %>%
  filter(gate%in%c("Tr1_IL10_Frequency"), stim!="unstim", targetName%in%c("GZMA", "LAG3", "IL10", "IFNG", "TNF", "IL6"))%>%
  ggplot(., aes(x=freq_base_d7_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("fold change in Tr1 Frequency & PD1 expression vs. fold change in Tr1 proteins")+
  geom_smooth(method="lm")+
  # scale_color_manual(values =stim_palette)+
  facet_wrap(gate~targetName)+
  theme_minimal()

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/tr1_gzma_fc_corr_day0_day14.png", tr1_gzma, height = 6, width=6, dpi=444, bg="white")

tr1_sigs <- sig_fc_corr_purr %>%
  filter(grepl("Tr1", gate), base_d7_padj<0.05)%>%
  select(gate, stim, targetName, base_d7_padj, base_d7_rho)


tr1_il10_freq_cor_plot <- fc_combo_frame %>%
  filter(gate %in% tr1_sigs$gate[tr1_sigs$gate=="Tr1_IL10_Frequency"],
         targetName %in% tr1_sigs$targetName[tr1_sigs$gate=="Tr1_IL10_Frequency"],
         infectiontype=="S",
         stim=="iRBCs")%>%
  ggplot(., aes(x=freq_base_d7_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman", label.y = 0.7)+
  xlab("fold change in Tr1 proteins baseline vs day 0")+
  ylab("fold change in Tr1 IL10 expression baseline vs day 7")+
  geom_smooth(method="lm")+
  ylim(c(0.65, NA))+
  scale_color_manual(values =stim_palette)+
  facet_wrap(~targetName, scales="free", nrow=2)+
  theme_minimal()+
  theme(legend.position = "bottom")

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/tr1_il10_freq_cor_plot.png", tr1_il10_freq_cor_plot, width=8, height=5.33, bg="white", dpi=444)

tr1_ifng_freq_cor_plot <- fc_combo_frame %>%
  filter(gate == "Tr1_IFNg_Frequency",
         targetName %in% tr1_sigs$targetName[tr1_sigs$gate=="Tr1_IFNg_Frequency" &tr1_sigs$stim=="PMA"],
         infectiontype=="S",
         stim=="PMA")%>%
  ggplot(., aes(x=freq_base_d7_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman", label.y = 0.7)+
  xlab("fold change in Tr1 proteins baseline vs day 0")+
  ylab(expression(paste("fold change in Tr1 IFN", gamma, sep="")))+
  geom_smooth(method="lm")+
  ylim(c(0.65, NA))+
  scale_color_manual(values =stim_palette)+
  facet_wrap(~targetName, scales="free")+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/tr1_ifng_freq_cor_plot.png", tr1_ifng_freq_cor_plot, width=5.33, height=3, bg="white", dpi=444)


tr1_freq_cor_plot <- fc_combo_frame %>%
  filter(gate %in% tr1_sigs$gate[tr1_sigs$gate=="Tr1_Frequency"],
         targetName %in% tr1_sigs$targetName[tr1_sigs$gate=="Tr1_Frequency"],
         infectiontype=="S",
         stim=="unstim")%>%
  ggplot(., aes(x=freq_base_d7_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman", label.y = 0.7)+
  xlab("fold change in Tr1 proteins baseline vs day 0")+
  ylab("fold change in Tr1 % baseline vs day 7")+
  geom_smooth(method="lm")+
  ylim(c(0.65, NA))+
  scale_color_manual(values =stim_palette)+
  facet_wrap(~targetName, scales="free")+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/tr1_freq_cor_plot.png", tr1_freq_cor_plot, width=5.33, height=3, bg="white", dpi=444)




sig_fc_corr_purr7 <- fc_corr_purr%>%
  filter(base_d7_padj<0.05)
  

list_of_plots <- list(matrix(nrow = nrow(sig_fc_corr_purr7)))

for(i in 1:nrow(sig_fc_corr_purr7)){
  
  plot_data <- fc_combo_frame %>%
    filter(gate==sig_fc_corr_purr7$gate[i],
           stim==sig_fc_corr_purr7$stim[i],
           infectiontype==sig_fc_corr_purr7$infectiontype[i],
           targetName==sig_fc_corr_purr7$targetName[i])
  
  plot <- ggplot(plot_data, aes(x=freq_base_d7_fc, y=conc_base_d0_fc, color=stim))+
    geom_point()+
    geom_smooth(method="lm")+
    ggpubr::stat_cor(method = "spearman", size=2)+
    scale_color_manual(values = stim_palette)+
    xlab(paste(sig_fc_corr_purr7$gate[i]))+
    ylab(paste(sig_fc_corr_purr7$targetName[i]))+
    theme_minimal()+
    theme(legend.position = "none")
  
  list_of_plots[[i]] <- plot
  
}

big_plot <- cowplot::plot_grid(plotlist = list_of_plots, nrow =  round(nrow(sig_fc_corr_purr7)/10))

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/sig_fc_corr_day0_day7.png", big_plot, height = 32, width=32, dpi=444, bg="white")






sig_fc_corr_purr14 <- fc_corr_purr%>%
  filter(base_d14_padj<0.05)


list_of_plots <- list(matrix(nrow = nrow(sig_fc_corr_purr14)))

for(i in 1:nrow(sig_fc_corr_purr14)){
  
  plot_data <- fc_combo_frame %>%
    filter(gate==sig_fc_corr_purr14$gate[i],
           stim==sig_fc_corr_purr14$stim[i],
           infectiontype==sig_fc_corr_purr14$infectiontype[i],
           targetName==sig_fc_corr_purr14$targetName[i])
  
  plot <- ggplot(plot_data, aes(x=freq_base_d14_fc, y=conc_base_d0_fc))+
    geom_point(aes(color=stim))+
    geom_smooth(method="lm")+
    ggpubr::stat_cor(method = "spearman", size=2)+
    scale_color_manual(values = stim_palette)+
    xlab(paste(sig_fc_corr_purr14$gate[i]))+
    ylab(paste(sig_fc_corr_purr14$targetName[i]))+
    theme_minimal()
  
  list_of_plots[[i]] <- plot
  
}

big_plot <- cowplot::plot_grid(plotlist = list_of_plots, nrow =  round(nrow(sig_fc_corr_purr14)/8))

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/sig_fc_corr_day0_day14.png", big_plot, height = 32, width=32, dpi=444, bg="white")




sig_fc_corr_purr0 <- fc_corr_purr%>%
  filter(base_d0_padj<0.05)


list_of_plots <- list(matrix(nrow = nrow(sig_fc_corr_purr0)))

for(i in 1:nrow(sig_fc_corr_purr0)){
  
  plot_data <- fc_combo_frame %>%
    filter(gate==sig_fc_corr_purr0$gate[i],
           stim==sig_fc_corr_purr0$stim[i],
           infectiontype==sig_fc_corr_purr0$infectiontype[i],
           targetName==sig_fc_corr_purr0$targetName[i])
  
  plot <- ggplot(plot_data, aes(x=freq_base_d0_fc, y=conc_base_d0_fc))+
    geom_point(aes(color=stim))+
    geom_smooth(method="lm")+
    ggpubr::stat_cor(method = "spearman", size=2)+
    scale_color_manual(values = stim_palette)+
    xlab(paste(sig_fc_corr_purr0$gate[i]))+
    ylab(paste(sig_fc_corr_purr0$targetName[i]))+
    theme_minimal()
  
  list_of_plots[[i]] <- plot
  
}

big_plot <- cowplot::plot_grid(plotlist = list_of_plots, nrow =  round(nrow(sig_fc_corr_purr0)/6))

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/sig_fc_corr_day0_day0.png", big_plot, height = 32, width=32, dpi=444, bg="white")










# sandbox####



fc_combo_frame %>%
  filter(gate%in%c("IFNg_Frequency", "IL10_Frequency"),
         stim=="iRBCs",
         targetName %in% c("TNFRSF17", "FLT1", "LILRB2"),
         infectiontype=="A")%>%
  ggplot(., aes(x=freq_base_d0_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  ggtitle("fold change in Tr1 Frequency & PD1 expression vs. fold change in Tr1 proteins")+
  geom_smooth(method="lm")+
  # scale_color_manual(values =stim_palette)+
  facet_wrap(gate~targetName, scales = "free")+
  theme_minimal()


cytometry_plot <- long_combo%>%
  filter(gate %in% c("IFNg_Frequency", "IL10_Frequency", "Tr1_IL10_Frequency", "Tr1_IFNg_Frequency"),
         infectiontype %in% c("A", "S"),
         #stim=="PMA",
         timepoint %in% c("baseline", "day0", "day7", "day14"))%>%
  distinct(timepoint, freq, gate, stim, infectiontype)%>%
  ggplot(., aes(x=timepoint, y=freq, fill=infectiontype))+
  geom_boxplot(outliers = F)+
  facet_wrap(~gate+stim, scales="free")+
  theme_minimal()

ggsave("~/Downloads/cytometry_plot.png", width=12, height=12, bg="white")
combo_data %>%
  filter(gate=="CD4_T_Cell_Frequency", infectiontype!="nmf", targetName=="CXCL10")%>%
  group_by(infectiontype, timepoint.x)%>%
  filter(!duplicated(id), !is.na(timepoint.x), timepoint.x%in% c("baseline", "day0", "day14"))%>%
  mutate(timepoint.x = factor(timepoint.x, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(., aes(x=qpcr, y=unstim))+
  geom_point(aes(color=timepoint.x))+
  scale_x_log10()+
  # geom_line(aes(group=id))+
  # geom_violin(aes(fill=infectiontype), draw_quantiles = quantile(seq(0,1,0.25)))+
  ggpubr::stat_cor(method = "spearman")+
  # ggtitle(paste(gate, stim, infectiontype, targetName))+
  geom_smooth(method="lm")+
  # scale_color_manual(values =stim_palette)+
  facet_wrap(~infectiontype)+
  theme_minimal()


# combo_data%>%
#   filter(targetName %in% c("IL10", "CXCL10", "TNF"))%>%
#   group_by(id, timepoint)%>%
#   ggplot(aes(x=timepoint.x, y=concentration, fill=timepoint.x))+
#   geom_point()+#
#   geom_line(aes(group=id))+
#   # geom_boxplot(outliers = FALSE)+
#   # geom_violin(draw_quantiles = seq(0,1,0.25))+
#   # ggtitle("regulated during asymptomatic parasitemia")+
#   facet_wrap(~targetName+infectiontype, scales = "free")+
#   scale_fill_manual(values=viridis::magma(4))+
#   theme_minimal()+
#   da_boxplot_theme


cd4_t_plot <- combo_data%>%
  filter(gate %in% c("CD4_T_Cell_Frequency"))%>%
  filter(infectiontype!="nmf", !is.na(timepoint.x))%>%
  ggplot(aes(x=factor(timepoint.x, levels=c("baseline", "day0", "day7", "day14")), y=unstim, fill=timepoint.x))+
  geom_line(aes(group=id))+
  geom_boxplot(outliers = FALSE)+
  geom_point(color="grey")+#
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  # ggtitle("regulated during asymptomatic parasitemia")+
  facet_wrap(gate~infectiontype, scales = "free")+
  scale_color_manual(values=viridis::magma(3))+
  scale_fill_manual(values=viridis::magma(4))+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/cd4_freqs.png", cd4_t_plot, width=4, height=4, bg="white", dpi=444)


# old luis data stuff ####


pilot_data <- read.csv("~/Library/CloudStorage/Box-Box/Jagannathan_Lab_Folder/PROJECTS/Tr1/MUSICAL_Tr1_Flow/Luis_Analysis/cell_freqs_pilot.csv")

list.files("~/Library/CloudStorage/Box-Box/Jagannathan_Lab_Folder/PROJECTS/Tr1/MUSICAL_Tr1_Flow/Luis_Analysis/CSV_Files/")

metadata_files <- list.files("~/Library/CloudStorage/Box-Box/Jagannathan_Lab_Folder/PROJECTS/Tr1/MUSICAL_Tr1_Flow/Luis_Analysis/CSV_Files/", pattern = "Metadata_MUS[0-9].csv", full.names = T)
data_files <-     list.files("~/Library/CloudStorage/Box-Box/Jagannathan_Lab_Folder/PROJECTS/Tr1/MUSICAL_Tr1_Flow/Luis_Analysis/CSV_Files/", pattern = "^FCdata_MUS[0-9].csv", full.names = T)

metadata_list <- list()
data_list <- list()

for(files in 1:length(metadata_files)){
  
  tmp <- read.csv(metadata_files[files])
  tmp$MusicalID <- paste("MUS", files, sep="")
  metadata_list[[files]] <- tmp[!is.na(tmp$plate),]
}



for(files in 1:length(data_files)){
  
  tmp <- read.csv(data_files[files])
  tmp$study <- "big batch"
  tmp$MusicalID <- paste("MUS", files, sep="")
  #ignore letter before 0, find 0, ignore number after 0
  tmp$well <- stringr::str_remove(tmp$well, "(?<=[A-Z])0(?=[0-9])")
  tmp <- tmp[!is.na(tmp$plate),]
  tmp <- tmp[tmp$X%notin%c("SD", "Mean"),]
  tmp <- tmp[,colnames(tmp)!="X.1"]
  
  if(any(grepl("Drop", colnames(tmp))))
    tmp <- tmp %>%
    filter(Drop == "")%>%
    select(-Drop)
  
  data_list[[files]] <- tmp
}

big_metadata <- do.call(rbind, metadata_list)
big_data <- do.call(rbind, data_list)

big_metadata$well_id <- paste(big_metadata$MusicalID, big_metadata$plate, big_metadata$well)
big_data$well_id <- paste(big_data$MusicalID, big_data$plate, big_data$well)

big_data <- big_data %>%
  mutate(plate=as.numeric(plate))

luis_data <- big_data%>%
  left_join(., big_metadata, by=c("MusicalID", "plate", "well"))%>%
  filter(!is.na(stim), timepoint!="na")%>%
  pivot_longer(cols=ends_with("Frequency"), names_to = "gate", values_to = "Freq")%>%
  pivot_wider(names_from = stim, values_from = Freq, id_cols = c(id, timepoint, inf_type, gate))%>%
  mutate(sample_id2=paste(id, inf_type, timepoint, sep="_"))



# create variable that can be matched to flow cytyometry metadata
slim_nulisa_data <- nulisa_data %>%
  mutate(sample_id2 = paste(substr(id, nchar(id)-2, nchar(id)),
                            case_when(class=="A"~"asymp",
                                      class=="S"~"symp"),
                            case_when(timepoint=="baseline"~"-1",
                                      timepoint=="day0"~"0",
                                      timepoint=="day7"~"7",
                                      timepoint=="day14"~"14",
                            ), sep="_"))%>%
  select(targetName, concentration, qpcr, temperature, sample_id2, timepoint)


combo_data <- slim_nulisa_data %>%
  right_join(., luis_data, by = "sample_id2", relationship = "many-to-many")%>%
  filter(!is.na(timepoint.x))



long_combo <-  combo_data %>%
  pivot_longer(cols=c(iRBC, unstim, PMA), names_to = "stim", values_to = "freq")%>%
  filter(!(stim=="PMA" & !grepl("_[^_]+_", gate)),
         !(stim=="iRBC" & !grepl("_[^_]+_", gate)),
         !(stim=="unstim" & grepl("_[^_]+_", gate)),
         !(stim=="PMA" & gate=="CD4_T_Cell_Frequency"))


fdr_cutoff <- 0.1

unstim_corrs <- combo_data %>%
  group_by(targetName, gate, timepoint.x, inf_type)%>%
  nest()%>%
  mutate(correlation=map(data, ~cor.test(.$concentration, .$unstim, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  ungroup()%>%
  mutate(padj=p.adjust(p))

sig_unstim <- unstim_corrs %>%
  filter(padj<fdr_cutoff)



PMA_corrs <- combo_data %>%
  group_by(targetName, gate, timepoint.x, inf_type)%>%
  nest()%>%
  mutate(correlation=map(data, ~cor.test(.$concentration, .$PMA, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  ungroup()%>%
  mutate(padj=p.adjust(p))


sig_PMA <- PMA_corrs %>%
  filter(padj<fdr_cutoff)



iRBC_corrs <- combo_data %>%
  group_by(targetName, gate, timepoint.x, inf_type)%>%
  nest()%>%
  mutate(correlation=map(data, ~cor.test(.$concentration, .$iRBC, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  ungroup()%>%
  mutate(padj=p.adjust(p))

sig_iRBC<- iRBC_corrs %>%
  filter(padj<fdr_cutoff)



combo_data %>%
  filter(grepl("^Tr1_", gate), targetName%in%c("LILRB2"))%>%
  pivot_longer(cols=c(iRBC, unstim, PMA), names_to = "stim", values_to = "freq")%>%
  ggplot(., aes(x=freq, y=concentration))+
  geom_point(aes(color=inf_type))+
  facet_wrap(~gate+stim+targetName, scales = "free")+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("orange", "darkblue"))+
  theme_minimal()

combo_data %>%
  filter(grepl("^Treg", gate), targetName=="IL10")%>%
  ggplot(., aes(x=PMA, y=concentration))+
  geom_point(aes(color=inf_type))+
  facet_wrap(~gate)+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("orange", "darkblue"))+
  theme_minimal()




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
    scale_fill_manual(values=viridis::magma(n=5))+
    da_boxplot_theme+
  theme(legend.position = "right",
        legend.title = element_blank())

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/musical_tr1_genes_boxplot.png", jason_box, width =8, height=6, dpi=444, bg="white")

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

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/tr1_freq_plasma_protein_cor.png", tr1_freq_plasma_protein_cor, width = 9, height=9, dpi=444, bg="white")




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

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/tr1_freq_base_plasma_protein_cor.png", tr1_freq_base_plasma_protein_cor, width = 8, height=3, dpi=444, bg="white")

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

 ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/all_cell_freq_base_day0_conc_tr1_genes.png", width = 9, height=9, dpi=444, bg="white")
 
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
 
 ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/all_tr1_freq_base_day0_conc.png", width = 9, height=4.5, dpi=444, bg="white")
 
 