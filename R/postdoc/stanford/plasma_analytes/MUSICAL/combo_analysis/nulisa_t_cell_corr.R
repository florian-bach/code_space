#preamble ####

library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)

stim_palette <- c("darkred", "darkblue", "black")
names(stim_palette) <- c("iRBC", "PMA", "unstim")

`%notin%` <- Negate(`%in%`)

# data generation ####
nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")

## read in cytometry data ####
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
  pivot_wider(names_from = stim, values_from = Freq, id_cols = c(cohortid, timepoint, inf_type, gate))%>%
  mutate(sample_id2=paste(cohortid, inf_type, timepoint, sep="_"))
  


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

## putting it all together ####
combo_data <- slim_nulisa_data %>%
  right_join(., luis_data, by = "sample_id2", relationship = "many-to-many")%>%
  filter(!is.na(timepoint.x))

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

combo_data %>%
  filter(gate=="Treg2_Frequency", targetName=="CD276")%>%
  ggplot(., aes(x=iRBC, y=concentration, color=inf_type))+
  geom_point()+
  facet_wrap(~gate)+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("orange", "darkblue"))+
  theme_minimal()
  

# disjointed timepoints, absolute values ####
long_combo <-  combo_data %>%
  pivot_longer(cols=c(iRBC, unstim, PMA), names_to = "stim", values_to = "freq")%>%
  filter(!(stim=="PMA" & !grepl("_[^_]+_", gate)),
         !(stim=="iRBC" & !grepl("_[^_]+_", gate)),
         !(stim=="unstim" & grepl("_[^_]+_", gate)),
         !(stim=="PMA" & gate=="CD4_T_Cell_Frequency"))

day0_concs <- long_combo %>%
  filter(inf_type=="symp", timepoint.x=="baseline")%>%
  distinct(targetName, concentration, cohortid, timepoint.x, inf_type)
  
day7_freqs <- long_combo%>%
  filter(inf_type=="symp", timepoint.x=="day7")%>%
  distinct(gate, freq, stim, cohortid, timepoint.x, inf_type)

day14_freqs <- long_combo%>%
  filter(inf_type=="symp", timepoint.x=="day14")%>%
  distinct(gate, freq, stim, cohortid, timepoint.x, inf_type)


conc_freq_07 <- inner_join(day0_concs, day7_freqs, by=c("cohortid", "inf_type"), relationship = "many-to-many")
conc_freq_014 <- inner_join(day0_concs, day14_freqs, by=c("cohortid", "inf_type"), relationship = "many-to-many")

conc_freq_07_corr <- conc_freq_07%>%
  group_by(gate, stim, targetName)%>%
  nest()%>%
  mutate(correlation=map(data, ~cor.test(.$concentration, .$freq, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  ungroup()%>%
  mutate(padj=p.adjust(p))

conc_freq_014_corr <- conc_freq_014%>%
  group_by(gate, stim, targetName)%>%
  nest()%>%
  mutate(correlation=map(data, ~cor.test(.$concentration, .$freq, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  ungroup()%>%
  mutate(padj=p.adjust(p))

sig_conc_freq_07_corr <- conc_freq_07_corr%>%
  filter(padj<fdr_cutoff)

sig_conc_freq_014_corr <- conc_freq_014_corr%>%
  filter(padj<fdr_cutoff)

long_combo %>%
  filter(inf_type=="symp", targetName %in% c("CRP", "IL10"), grepl("^Treg2", gate))%>%
  ggplot(., aes(x=freq, y=concentration))+
  geom_point()+
  facet_wrap(~gate+stim)+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("darkred", "darkblue", "white"))+
  theme_minimal()
  

pilot_data <- read.csv("~/Library/CloudStorage/Box-Box/Jagannathan_Lab_Folder/PROJECTS/Tr1/MUSICAL_Tr1_Flow/Luis_Analysis/cell_freqs_pilot.csv")

# fold changes ####

wide_concs <-  combo_data %>%
  pivot_longer(cols=c(iRBC, unstim, PMA), names_to = "stim", values_to = "freq")%>%
  filter(!(stim=="PMA" & !grepl("_[^_]+_", gate)),
         !(stim=="iRBC" & !grepl("_[^_]+_", gate)),
         !(stim=="unstim" & grepl("_[^_]+_", gate)),
         !(stim=="PMA" & gate=="CD4_T_Cell_Frequency"))%>%
  select(targetName, concentration, cohortid, inf_type, timepoint.x, sample_id2)%>%
  filter(inf_type!="nmf", timepoint.x %in% c("baseline", "day0", "day7", "day14"))%>%
  distinct(targetName, concentration, cohortid, inf_type, timepoint.x)%>%
  # distinct(c(targetName, sample_id2))%>%
  # pivot_longer(cols=c(iRBC, unstim, PMA), names_to = "stim", values_to = "freq")%>%
  pivot_wider(names_from = timepoint.x, values_from = c(concentration), id_cols=c(targetName, cohortid, inf_type), names_prefix = "conc_")%>%
  group_by(inf_type)%>%
  mutate(conc_base_d0_fc=conc_day0/conc_baseline)


wide_freqs <-  combo_data %>%
  pivot_longer(cols=c(iRBC, unstim, PMA), names_to = "stim", values_to = "freq")%>%
  filter(!(stim=="PMA" & !grepl("_[^_]+_", gate)),
         !(stim=="iRBC" & !grepl("_[^_]+_", gate)),
         !(stim=="unstim" & grepl("_[^_]+_", gate)),
         !(stim=="PMA" & gate=="CD4_T_Cell_Frequency"))%>%
  select(gate, stim, freq, cohortid, inf_type, timepoint.x, sample_id2)%>%
  filter(inf_type!="nmf", timepoint.x %in% c("baseline", "day0", "day7", "day14"))%>%
  distinct(gate, stim, freq, cohortid, inf_type, timepoint.x)%>%
  # distinct(c(targetName, sample_id2))%>%
  # pivot_longer(cols=c(iRBC, unstim, PMA), names_to = "stim", values_to = "freq")%>%
  pivot_wider(names_from = timepoint.x, values_from = c(freq), id_cols=c(gate, stim, cohortid, inf_type), names_prefix = "freq_")%>%
  group_by(inf_type)%>%
  mutate(freq_base_d0_fc=freq_day0/freq_baseline,
         freq_base_d7_fc=freq_day7/freq_baseline,
         freq_base_d14_fc=freq_day14/freq_baseline)

fc_combo_frame <- left_join(wide_freqs, wide_concs, by=c("cohortid", "inf_type"))

fc_corr_purr <- fc_combo_frame %>%
  filter(inf_type=="symp")%>%
  group_by(inf_type, stim, gate, targetName)%>%
  nest()%>%
  # mutate(base_d7_correlation=if_else(all(is.na(data$freq_day7)), map(data, ~cor.test(.$freq_base_d7_fc, .$conc_base_d0_fc, method = "spearman"))))%>%
  mutate(base_d7_correlation=map(data, ~cor.test(.$freq_base_d7_fc, .$conc_base_d0_fc, method = "spearman")))%>%
  mutate(base_d14_correlation=map(data, ~cor.test(.$freq_base_d14_fc, .$conc_base_d0_fc, method = "spearman")))%>%
  mutate(
         base_d7_p=map_dbl(base_d7_correlation, ~.$p.value),
         base_d7_rho=map_dbl(base_d7_correlation, ~.$estimate),
         base_d14_p=map_dbl(base_d14_correlation, ~.$p.value),
         base_d14_rho=map_dbl(base_d14_correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  group_by(inf_type, stim, targetName)%>%
  mutate(base_d14_padj=p.adjust(base_d14_p),
         base_d7_padj=p.adjust(base_d7_p)
         )
  
sig_fc_corr_purr <- fc_corr_purr%>%
  filter(base_d14_padj<0.05 | base_d7_padj<0.05)


tr1_gzma <- fc_combo_frame %>%
  filter(gate=="Tr1_Frequency", stim!="PMA", targetName%in%c("GZMA"))%>%
  ggplot(., aes(x=freq_base_d7_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor()+
  ggtitle("Tr1 Frequency change vs. GZMA change")+
  geom_smooth(method="lm")+
  scale_color_manual(values =stim_palette)+
  facet_wrap(~targetName)+
  theme_minimal()

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/tr1_gzma_fc_corr_day0_day14.png", tr1_gzma, height = 6, width=6, dpi=444, bg="white")

ctfh_cxcl10 <- fc_combo_frame %>%
  filter(gate=="cTfh_Frequency", stim!="PMA", targetName%in%c("CXCL10"))%>%
  ggplot(., aes(x=freq_base_d14_fc, y=conc_base_d0_fc, color=stim))+
  geom_point()+
  ggpubr::stat_cor()+
  ggtitle("cTfh change vs. CXCL10 change")+
  geom_smooth(method="lm")+
  scale_color_manual(values =stim_palette)+
  facet_wrap(~targetName)+
  theme_minimal()

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/ctfh_cxcl10_fc_corr_day0_day14.png", ctfh_cxcl10, height = 6, width=6, dpi=444, bg="white")



list_of_plots <- list(matrix(nrow = nrow(sig_fc_corr_purr)))


for(i in 1:nrow(sig_fc_corr_purr)){
  
  plot_data <- fc_combo_frame %>%
    filter(gate==sig_fc_corr_purr$gate[i],
           stim==sig_fc_corr_purr$stim[i],
           inf_type==sig_fc_corr_purr$inf_type[i],
           targetName==sig_fc_corr_purr$targetName[i])
  
  plot <- ggplot(plot_data, aes(x=freq_base_d14_fc, y=conc_base_d0_fc))+
    geom_point(aes(color=stim))+
    geom_smooth(method="lm")+
    ggpubr::stat_cor(method = "spearman", na.rm = TRUE, size=2)+
    scale_color_manual(values = stim_palette)+
    xlab(paste(sig_fc_corr_purr$gate[i]))+
    ylab(paste(sig_fc_corr_purr$targetName[i]))+
    theme_minimal()
  
  list_of_plots[[i]] <- plot
  
}

big_plot <- cowplot::plot_grid(plotlist = list_of_plots, nrow =  round(nrow(sig_fc_corr_purr)/5))

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/sig_fc_corr_day0_day14.png", big_plot, height = 32, width=16, dpi=444, bg="white")




fc_corr_purr7 <- fc_combo_frame%>%
  filter(inf_type=="symp")%>%
  group_by(stim, gate, targetName)%>%
  nest()%>%
  # mutate(base_d7_correlation=if_else(all(is.na(freq_day7)), map(data, ~cor.test(.$freq_base_d7_fc, .$conc_base_d0_fc, method = "spearman")))%>%
  mutate(base_d7_correlation=map(data, ~cor.test(.$freq_base_d7_fc, .$conc_base_d0_fc, method = "spearman", )))%>%
  mutate(
    # base_d7_p=map_dbl(base_d7_correlation, ~.$p.value),
    # base_d7_rho=map_dbl(base_d7_correlation, ~.$estimate),
    base_d7_p=map_dbl(base_d7_correlation, ~.$p.value),
    base_d7_rho=map_dbl(base_d7_correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  group_by(inf_type, stim, targetName)%>%
  mutate(padj=p.adjust(base_d7_p))

sig_fc_corr_purr7 <- fc_corr_purr7%>%
  filter(padj<0.05)
  

list_of_plots <- list(matrix(nrow = nrow(sig_fc_corr_purr)))

for(i in 1:nrow(sig_fc_corr_purr7)){
  
  plot_data <- fc_combo_frame %>%
    filter(gate==sig_fc_corr_purr7$gate[i],
           stim==sig_fc_corr_purr7$stim[i],
           inf_type==sig_fc_corr_purr7$inf_type[i],
           targetName==sig_fc_corr_purr7$targetName[i])
  
  plot <- ggplot(plot_data, aes(x=freq_base_d7_fc, y=conc_base_d0_fc))+
    geom_point(aes(color=stim))+
    geom_smooth(method="lm")+
    ggpubr::stat_cor(method = "spearman", size=2)+
    scale_color_manual(values = stim_palette)+
    xlab(paste(sig_fc_corr_purr$gate[i]))+
    ylab(paste(sig_fc_corr_purr$targetName[i]))+
    theme_minimal()
  
  list_of_plots[[i]] <- plot
  
}

big_plot <- cowplot::plot_grid(plotlist = list_of_plots, nrow =  round(nrow(sig_fc_corr_purr7)/5))

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/sig_fc_corr_day0_day7.png", big_plot, height = 32, width=16, dpi=444, bg="white")






fc_corr_purr0 <- fc_combo_frame%>%
  # filter(inf_type=="symp")%>%
  group_by(inf_type, stim, gate, targetName)%>%
  nest()%>%
  # mutate(base_d7_correlation=if_else(all(is.na(freq_day7)), map(data, ~cor.test(.$freq_base_d7_fc, .$conc_base_d0_fc, method = "spearman")))%>%
  mutate(base_d0_correlation=map(data, ~cor.test(.$freq_base_d0_fc, .$conc_base_d0_fc, method = "spearman", )))%>%
  mutate(
    # base_d0_p=map_dbl(base_d0_correlation, ~.$p.value),
    # base_d0_rho=map_dbl(base_d0_correlation, ~.$estimate),
    base_d0_p=map_dbl(base_d0_correlation, ~.$p.value),
    base_d0_rho=map_dbl(base_d0_correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  group_by(inf_type, stim, targetName)%>%
  mutate(padj=p.adjust(base_d0_p))

sig_fc_corr_purr0 <- fc_corr_purr0%>%
  filter(padj<0.05)


list_of_plots <- list(matrix(nrow = nrow(sig_fc_corr_purr)))

for(i in 1:nrow(sig_fc_corr_purr7)){
  
  plot_data <- fc_combo_frame %>%
    filter(gate==sig_fc_corr_purr7$gate[i],
           stim==sig_fc_corr_purr7$stim[i],
           inf_type==sig_fc_corr_purr7$inf_type[i],
           targetName==sig_fc_corr_purr7$targetName[i])
  
  plot <- ggplot(plot_data, aes(x=freq_base_d7_fc, y=conc_base_d0_fc))+
    geom_point(aes(color=stim))+
    geom_smooth(method="lm")+
    ggpubr::stat_cor(method = "spearman", size=2)+
    scale_color_manual(values = stim_palette)+
    xlab(paste(sig_fc_corr_purr$gate[i]))+
    ylab(paste(sig_fc_corr_purr$targetName[i]))+
    theme_minimal()
  
  list_of_plots[[i]] <- plot
  
}

big_plot <- cowplot::plot_grid(plotlist = list_of_plots, nrow =  round(nrow(sig_fc_corr_purr7)/5))

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/sig_fc_corr_day0_day7.png", big_plot, height = 32, width=16, dpi=444, bg="white")





# sandbox

combo_data %>%
  filter(gate=="CD4_T_Cell_Frequency", inf_type!="nmf", targetName=="CXCL10")%>%
  group_by(inf_type, timepoint.x)%>%
  filter(!duplicated(cohortid), !is.na(timepoint.x), timepoint.x%in% c("baseline", "day0", "day14"))%>%
  mutate(timepoint.x = factor(timepoint.x, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(., aes(x=qpcr, y=unstim))+
  geom_point(aes(color=timepoint.x))+
  scale_x_log10()+
  # geom_line(aes(group=cohortid))+
  # geom_violin(aes(fill=inf_type), draw_quantiles = quantile(seq(0,1,0.25)))+
  ggpubr::stat_cor()+
  # ggtitle(paste(gate, stim, inf_type, targetName))+
  geom_smooth(method="lm")+
  # scale_color_manual(values =stim_palette)+
  facet_wrap(~inf_type)+
  theme_minimal()


# combo_data%>%
#   filter(targetName %in% c("IL10", "CXCL10", "TNF"))%>%
#   group_by(cohortid, timepoint)%>%
#   ggplot(aes(x=timepoint.x, y=concentration, fill=timepoint.x))+
#   geom_point()+#
#   geom_line(aes(group=cohortid))+
#   # geom_boxplot(outliers = FALSE)+
#   # geom_violin(draw_quantiles = seq(0,1,0.25))+
#   # ggtitle("regulated during asymptomatic parasitemia")+
#   facet_wrap(~targetName+inf_type, scales = "free")+
#   scale_fill_manual(values=viridis::magma(4))+
#   theme_minimal()+
#   da_boxplot_theme


cd4_t_plot <- combo_data%>%
  filter(gate %in% c("CD4_T_Cell_Frequency"))%>%
  filter(inf_type!="nmf", !is.na(timepoint.x))%>%
  ggplot(aes(x=factor(timepoint.x, levels=c("baseline", "day0", "day7", "day14")), y=unstim, fill=timepoint.x))+
  geom_line(aes(group=cohortid))+
  geom_boxplot(outliers = FALSE)+
  geom_point(color="grey")+#
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  # ggtitle("regulated during asymptomatic parasitemia")+
  facet_wrap(gate~inf_type, scales = "free")+
  scale_color_manual(values=viridis::magma(3))+
  scale_fill_manual(values=viridis::magma(4))+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/cd4_freqs.png", cd4_t_plot, width=4, height=4, bg="white", dpi=444)
