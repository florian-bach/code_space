library(emmeans)
library(tidyr)
library(dplyr)
library(ggplot2)

# read in data

nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/tr1_paper/revised_baseline_clean_musical_combo_with_metadata.csv")

nulisa_data <- nulisa_data%>%
  mutate(timepoint=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))


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
  group_by(targetName, infectiontype)%>%
  mutate(conc_base_d0_fc=conc_day0-conc_baseline,
         conc_base_d14_fc=conc_day14-conc_baseline)%>%
  distinct(targetName, id, infectiontype, conc_base_d0_fc, conc_base_d14_fc)%>%
  group_by(targetName, infectiontype)%>%
  summarise("mean_conc_base_d0_fc"=mean(conc_base_d0_fc, na.rm = T), "mean_conc_base_d14_fc"=mean(conc_base_d14_fc, na.rm = T))

combo_as_results_with_fc <- combo_as_results%>%
  left_join(., fc_data, by="targetName")%>%
  distinct(targetName, contrast,padj, infectiontype, mean_conc_base_d0_fc, mean_conc_base_d14_fc )
  

# Panel C ####
analytes_of_interest <- c("IL10", "LAG3", "GZMA", "IL27", "IL6", "TNF", "IFNG", "CTLA4", "LILRB2", "CRP", "CCL24", "IL7R", "CX3CL1", "KDR", "TNFSF11")

(base_day0_s_volcano <- combo_as_results_with_fc %>%
    filter(contrast == "baseline S - day0 S")%>%
    distinct(targetName, padj, infectiontype, mean_conc_base_d0_fc)%>%
    filter(infectiontype=="S", is.finite(mean_conc_base_d0_fc))%>%
    mutate("label2" = if_else(targetName %in% analytes_of_interest, targetName, NA))%>%
    ggplot(., aes(x=mean_conc_base_d0_fc, y=-log10(padj+10^-14), alpha=padj<0.05&abs(mean_conc_base_d0_fc)>0.5, color=mean_conc_base_d0_fc<0))+
    geom_point()+
    ggrepel::geom_text_repel(aes(label=label2, alpha=padj<0.05&abs(mean_conc_base_d0_fc)>0.5),
                            size = 5,min.segment.length = 0,
                            position=ggpp::position_nudge_center(center_x = 0, x = 2, y=0.00000000001, center_y = 5))+
    ggtitle("n = 112")+
    geom_hline(yintercept = -log10(0.05+10^-14), linetype="dashed")+
    geom_vline(xintercept = -0.5, linetype="dashed")+
    geom_vline(xintercept = 0.5, linetype="dashed")+
    scale_alpha_manual(values=c(0.5, 1))+
    scale_color_manual(values=c("darkred", "darkblue"))+
    scale_x_continuous(limits=c(-5, 5))+
    xlab("log2 fold change")+
    ylab("-log10 q value")+
    theme_minimal()+
    theme(legend.position="none"))

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_day0_s_volcano.png", base_day0_s_volcano, width = 5, height = 5, dpi=444, bg="white")


# Panel D ####

tr1_proteins_symp_plot <- nulisa_data %>%
  filter(targetName %in% c("IL10", "LAG3", "GZMA", "IFNG"),
         infectiontype=="S",
         timepoint %in% names(time_cols))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
  geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetName, scales = "free")+
  scale_fill_manual(values=time_cols)+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/tr1_proteins_symp_plot.png", tr1_proteins_symp_plot, width = 5.33333, height = 5.3333, dpi=444, bg="white")


# Panel E ####
cell_count_data <- read.csv("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/t_cell_cytometry_count_data.csv")

tr1_freq_plot <- cell_count_data %>%
  filter(infectiontype %in% c("S", "A"), timepoint!="bad_baseline", stim=="unstim")%>%
  distinct(id, infectiontype, timepoint, absolute_Tr1)%>%
  ggplot(., aes(x=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")), y=absolute_Tr1, fill=timepoint))+
  geom_boxplot(outliers = F)+
  # scale_y_log10()+
  ylab("Tr1 cells / μL")+
  facet_wrap(~infectiontype, scales = "free_x")+
  scale_fill_manual(values=time_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        legend.position = "none")

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/tr1_freq_plot.png", tr1_freq_plot, width =5.3333333, height=3, dpi=444, bg="white")



tr1_freq_plot <- cell_count_data %>%
  mutate(day14_para=if_else(timepoint=="day14" & parasitedensity.x > 10 &infectiontype=="A", "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id, infectiontype)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), ">10 parasites / μL on day 14", "no parasites day14"))%>%
  filter(infectiontype %in% c("S", "A"), timepoint!="bad_baseline", stim=="unstim", !is.na(class2))%>%
  distinct(id, infectiontype, timepoint, absolute_Tr1, class2)%>%
  ggplot(., aes(x=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")), y=absolute_Tr1, fill=timepoint, color=class2))+
  geom_boxplot(outliers = F)+
  # scale_y_log10()+
  ylab("Tr1 cells / μL")+
  facet_wrap(~infectiontype, scales = "free_x")+
  scale_fill_manual(values=time_cols)+
  scale_color_manual(values=c("darkred", "black"))+
  guides(fill=guide_none(), color=guide_legend(title = NULL))+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/tr1_freq_plot_ctrl.png", tr1_freq_plot, width =5.3333333, height=3, dpi=444, bg="white")

# Panel F ####

long_combo <- read.csv("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/cell_freqs_and_nulisa.csv")


grand_cor <- long_combo%>%
  filter(gate=="Tr1_Frequency", stim=="unstim", infectiontype%in% c("A", "S"), !is.na(targetName))%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(correlation=map(data, ~cor.test(.$concentration, .$freq, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  ungroup()%>%
  mutate(padj=p.adjust(p))

sig_grand_cor <- grand_cor%>%
  filter(padj<0.05)


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

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/big_correlation_tr1_proteins.png", target_cor_plot, width=5.33333, height=5.33333, bg="white")


base_target_cor_plot <- long_combo%>%
  filter(gate=="Tr1_Frequency", stim=="unstim", infectiontype%in%c("A", "S"), targetName %in%  c("IL10", "LAG3", "GZMA", "IFNG"), timepoint=="baseline")%>%
  distinct(id, timepoint, infectiontype, targetName, gate, freq, stim, concentration)%>%
  ggplot(., aes(x=freq, y=concentration))+
  geom_point(aes(fill=timepoint, shape=infectiontype), stroke=0.1)+
  ggpubr::stat_cor(method="spearman")+
  geom_smooth(aes(color=infectiontype), method="lm", color="black")+
  xlab("Tr1% of non-naive CD4")+
  scale_fill_manual(values=time_cols)+
  scale_shape_manual(values=c(21, 25))+
  facet_wrap(~targetName, scales = "free")+
  guides(fill=guide_none())+
  theme_minimal()

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_target_cor_plots.png", base_target_cor_plot, width=5.33333, height=5.33333, bg="white")

# Panel G (part of E?)####
# Panel H ####
tr1_proteins_asymp_plot <- nulisa_data %>%
  filter(targetName %in% c("IL10", "LAG3", "GZMA", "IFNG"),
         infectiontype=="A",
         timepoint %in% names(time_cols))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
  geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetName, scales = "free")+
  scale_fill_manual(values=time_cols)+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/tr1_proteins_asymp_plot.png", tr1_proteins_asymp_plot, width = 5.33333, height = 5.3333, dpi=444, bg="white")


# Panel I ####
day_a_sigs <- combo_as_results_with_fc %>%
  filter(contrast == "baseline A - day0 A", infectiontype=="A", padj<0.05, abs(mean_conc_base_d0_fc)>0.5)

(base_day0_a_volcano <- combo_as_results_with_fc %>%
   filter(contrast == "baseline A - day0 A")%>%
   distinct(targetName, padj, infectiontype, mean_conc_base_d0_fc)%>%
   filter(infectiontype=="A", is.finite(mean_conc_base_d0_fc))%>%
   mutate("label2" = if_else(targetName %in% c("CLEC4A", "CSF3", analytes_of_interest), targetName, NA))%>%
   ggplot(., aes(x=mean_conc_base_d0_fc, y=-log10(padj+10^-14), alpha=padj<0.05&abs(mean_conc_base_d0_fc)>0.5, color=mean_conc_base_d0_fc<0))+
   geom_point()+
   ggrepel::geom_text_repel(aes(label=label2, alpha=padj<0.05&abs(mean_conc_base_d0_fc)>0.5),
                            size = 5,min.segment.length = 0,
                            position=ggpp::position_nudge_center(center_x = 0, x = 2, y=0.00000000001, center_y = 5))+
   ggtitle("n = 9")+
   geom_hline(yintercept = -log10(0.05+10^-14), linetype="dashed")+
   geom_vline(xintercept = -0.5, linetype="dashed")+
   geom_vline(xintercept = 0.5, linetype="dashed")+
   scale_alpha_manual(values=c(0.5, 1))+
   scale_color_manual(values=c("darkred", "darkblue"))+
   scale_x_continuous(limits=c(-5, 5))+
   xlab("log2 fold change")+
   ylab("-log10 q value")+
   theme_minimal()+
   theme(legend.position="none"))

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/base_day0_a_volcano.png", base_day0_a_volcano, width = 5, height = 5, dpi=444, bg="white")

# panel J ####
ctrl_tr1_plot <- nulisa_data %>%
  filter(infectiontype %in% c("A"), timepoint!="day28", timepoint!="bad_baseline")%>%
  mutate(day14_para=if_else(timepoint=="day14" & qpcr > 10, "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), ">10 parasites / μL on day 14", "no parasites day14"))%>%
  # filter(targetName %in% sig_ctrl$targetName)%>%
  filter(targetName %in% c("IL10", "LAG3", "GZMA", "IFNG"))%>%#[seq((i-1)*16+1, i*16)])%>%
  filter(!is.na(class2))%>%
  # mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=timepoint, y=concentration, group=interaction(class2, timepoint)))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(aes(fill=timepoint, color=class2), outliers = FALSE)+
  facet_wrap(~targetName, scales = "free_y", nrow=2)+
  scale_fill_manual(values=time_cols)+
  scale_color_manual(values=c("#8B0000", "black"))+
  theme_minimal()+
  guides(fill=guide_none())+
  theme(legend.title = element_blank(),
         legend.position = "bottom",
         axis.title = element_blank())

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/ctrl_tr1_plot.png", ctrl_tr1_plot, width = 5.33333, height = 5.3333, dpi=444, bg="white")

