library(tidyr)
library(dplyr)
library(ggplot2)

da_boxplot_theme <- theme(legend.position = "none",
                          axis.title = element_blank())

time_cols <- list("baseline"="#E4DEBD",
                  "day0" = "#C03F3E",
                  "day7" = "#D87E1F",
                  "day14" = "#E6B85F")

combo_as_results <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/differential_abundance/combo_as_purff.csv", check.names = FALSE)

clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")

fc_data <- clean_data %>%
  pivot_wider(names_from = timepoint, values_from = c(concentration), id_cols=c(targetName, id, infectiontype), names_prefix = "conc_")%>%
  group_by(infectiontype)%>%
  mutate(conc_base_d0_fc=conc_day0/conc_baseline,
         conc_base_d14_fc=conc_day14/conc_baseline)%>%
  distinct(targetName, id, infectiontype, conc_base_d0_fc, conc_base_d14_fc)%>%
  group_by(targetName, infectiontype)%>%
  summarise("mean_conc_base_d0_fc"=mean(conc_base_d0_fc, na.rm = T), "mean_conc_base_d14_fc"=mean(conc_base_d14_fc, na.rm = T))

combo_as_purff%>%
  group_by(contrast)%>%
  summarise("n"=sum(padj<0.05))

# symptomatic day 0 ####
increase_at_day0_S <- combo_as_results %>%
  filter(padj<0.05, contrast=="baseline S - day0 S", coef<0)%>%
  arrange(padj)

decrease_at_day0_S <- combo_as_results %>%
  filter(padj<0.05, contrast=="baseline S - day0 S", coef>0)%>%
  arrange(padj)

top9_s_up <- clean_data %>%
  filter(infectiontype %in% c("S"), timepoint!="day28", timepoint!="bad_baseline")%>%
  filter(targetName %in% increase_at_day0_S$targetName[1:9])%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
  geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~factor(targetName, levels=increase_at_day0_S$targetName), scales = "free")+
  scale_fill_manual(values=time_cols)+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/top9_s_up.png", top9_s_up, width = 8, height = 8, dpi=444, bg="white")

top9_s_down <- clean_data %>%
  filter(infectiontype %in% c("S"), timepoint!="day28", timepoint!="bad_baseline")%>%
  filter(targetName %in% decrease_at_day0_S$targetName[1:9])%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
  geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~factor(targetName, levels=decrease_at_day0_S$targetName), scales = "free")+
  scale_fill_manual(values=time_cols)+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/top9_s_down.png", top9_s_down, width = 8, height = 8, dpi=444, bg="white")

analytes_of_interest <- c("IL10", "IL1RN", "LAG3", "GZMA", "IL27", "IL6", "TNF", "IFNG", "CTLA4", "LILRB2", "CRP")
(base_day0_s_volcano <- combo_as_results %>%
  filter(contrast == "baseline S - day0 S")%>%
  left_join(., fc_data, by="targetName")%>%
  distinct(targetName, coef, padj, infectiontype, mean_conc_base_d0_fc)%>%
  filter(infectiontype=="S", is.finite(mean_conc_base_d0_fc))%>%
  mutate("label2" = if_else(targetName %in% analytes_of_interest, targetName, NA))%>%
  ggplot(., aes(x=log2(mean_conc_base_d0_fc), y=-log10(padj+10^-14), alpha=padj<0.05, color=coef>0))+
  geom_point()+
  # geom_label(aes(label=label2), nudge_x = 0.1, )+
  ggrepel::geom_text_repel(aes(label=label2), nudge_x = 0.12)+
  ggtitle("n = 153")+
  geom_hline(yintercept = -log10(0.05+10^-14), linetype="dashed")+
  scale_alpha_manual(values=c(0.5, 1))+
  scale_color_manual(values=c("darkred", "darkblue"))+
  scale_x_continuous(limits=c(-0.75, 0.75))+
  xlab("log2 fold change")+
  ylab("-log10 P value")+
  theme_minimal()+
  theme(legend.position="none"))

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/base_day0_s_volcano.png", base_day0_s_volcano, width = 5, height = 5, dpi=444, bg="white")
ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/base_day0_s_volcano.png", base_day0_s_volcano, width = 5, height = 5, dpi=444, bg="white")


# asyptomatic day 0 ####

(base_day0_a_volcano <- combo_as_results %>%
  filter(contrast == "baseline A - day0 A")%>%
  left_join(., fc_data, by="targetName")%>%
  distinct(targetName, coef, padj, infectiontype, mean_conc_base_d0_fc)%>%
  filter(infectiontype=="A", is.finite(mean_conc_base_d0_fc))%>%
  mutate("label2" = if_else(targetName %in% analytes_of_interest, targetName, NA))%>%
  ggplot(., aes(x=log2(mean_conc_base_d0_fc), y=-log10(padj+10^-14), alpha=padj<0.05, color=coef>0))+
  geom_point()+
  # geom_label(aes(label=label2), nudge_x = 0.1, )+
  ggrepel::geom_text_repel(aes(label=label2), nudge_x = 0.12)+
  ggtitle("n = 36")+
  geom_hline(yintercept = -log10(0.05+10^-14), linetype="dashed")+
  scale_alpha_manual(values=c(0.5, 1))+
  scale_color_manual(values=c("darkred", "darkblue"))+
  scale_x_continuous(limits=c(-0.75, 0.75))+
  xlab("log2 fold change")+
  ylab("-log10 P value")+
  theme_minimal()+
  theme(legend.position="none"))

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/base_day0_a_volcano.png", base_day0_a_volcano, width = 5, height = 5, dpi=444, bg="white")
ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/base_day0_a_volcano.png", base_day0_a_volcano, width = 5, height = 5, dpi=444, bg="white")


increase_at_day0_a <- combo_as_results %>%
  filter(padj<0.05, contrast=="baseline A - day0 A", coef<0)%>%
  arrange(padj)

decrease_at_day0_a <- combo_as_results %>%
  filter(padj<0.05, contrast=="baseline A - day0 A", coef>0)%>%
  arrange(padj)

top9_a_up <- clean_data %>%
  filter(infectiontype %in% c("A"), timepoint!="day28", timepoint!="bad_baseline")%>%
  filter(targetName %in% increase_at_day0_a$targetName[1:9])%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
  geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~factor(targetName, levels=increase_at_day0_a$targetName), scales = "free")+
  scale_fill_manual(values=time_cols)+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/top9_a_up.png", top9_a_up, width = 5.33, height = 5.33, dpi=444, bg="white")

top9_a_down <- clean_data %>%
  filter(infectiontype %in% c("A"), timepoint!="day28", timepoint!="bad_baseline")%>%
  filter(targetName %in% decrease_at_day0_a$targetName[1:9])%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
  geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~factor(targetName, levels=decrease_at_day0_a$targetName), scales = "free")+
  scale_fill_manual(values=time_cols)+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/top9_a_down.png", top9_a_down, width = 8, height = 8, dpi=444, bg="white")


# asyptomatic day 14 ####

base_day14_a_volcano <- combo_as_results %>%
  filter(contrast == "baseline A - day14 A")%>%
  left_join(., fc_data, by="targetName")%>%
  distinct(targetName, coef, padj, infectiontype, mean_conc_base_d0_fc)%>%
  filter(infectiontype=="A", is.finite(mean_conc_base_d0_fc))%>%
  mutate("label2" = if_else(targetName %in% analytes_of_interest, targetName, NA))%>%
  ggplot(., aes(x=log2(mean_conc_base_d0_fc), y=-log10(padj+10^-14), alpha=padj<0.05, color=coef>0))+
  geom_point()+
  # geom_label(aes(label=label2), nudge_x = 0.1, )+
  ggrepel::geom_text_repel(aes(label=label2), nudge_x = 0.12)+
  ggtitle("n = 35")+
  geom_hline(yintercept = -log10(0.05+10^-14), linetype="dashed")+
  scale_alpha_manual(values=c(0.5, 1))+
  scale_color_manual(values=c("darkred", "darkblue"))+
  scale_x_continuous(limits=c(-0.75, 0.75))+
  xlab("log2 fold change")+
  ylab("-log10 P value")+
  theme_minimal()+
  theme(legend.position="none")

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/base_day14_a_volcano.png", base_day14_a_volcano, width = 5, height = 5, dpi=444, bg="white")
ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/base_day14_a_volcano.png", base_day14_a_volcano, width = 5, height = 5, dpi=444, bg="white")


ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/base_day14_a_volcano.png", base_day14_a_volcano, width = 5, height = 5, dpi=444, bg="white")


increase_at_day14_a <- combo_as_results %>%
  filter(padj<0.05, contrast=="baseline A - day14 A", coef<0)%>%
  arrange(padj)

decrease_at_day14_a <- combo_as_results %>%
  filter(padj<0.05, contrast=="baseline A - day14 A", coef>0)%>%
  arrange(padj)

top9_a_up <- clean_data %>%
  filter(infectiontype %in% c("A"), timepoint!="day28", timepoint!="bad_baseline")%>%
  filter(targetName %in% increase_at_day14_a$targetName[1:9])%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
  geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~factor(targetName, levels=increase_at_day14_a$targetName), scales = "free")+
  scale_fill_manual(values=time_cols)+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/top9_day14_a_up.png", top9_a_up, width = 8, height = 8, dpi=444, bg="white")

top9_a_down <- clean_data %>%
  filter(infectiontype %in% c("A"), timepoint!="day28", timepoint!="bad_baseline")%>%
  filter(targetName %in% decrease_at_day14_a$targetName[1:9])%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
  geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~factor(targetName, levels=decrease_at_day14_a$targetName), scales = "free")+
  scale_fill_manual(values=time_cols)+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/top9_day14_a_down.png", top9_a_down, width = 8, height = 8, dpi=444, bg="white")


symp_target <- combo_as_results %>%
  filter(padj<0.05, contrast=="baseline S - day0 S")

asymp_target2 <- combo_as_results %>%
  filter(padj<0.05, contrast=="baseline A - day14 A")

asymp_target <- combo_as_results %>%
  filter(padj<0.05, contrast=="baseline A - day0 A")

a_only <- asymp_target2$targetName[asymp_target2$targetName %notin% symp_target$targetName]


a_only_plot <- clean_data %>%
  filter(infectiontype %in% c("A"), timepoint!="day28", timepoint!="bad_baseline")%>%
  filter(targetName %in% a_only)%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
  geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetName, scales = "free")+
  scale_fill_manual(values=time_cols)+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/a_only_plot.png", a_only_plot, width = 8, height = 5.3333, dpi=444, bg="white")


tr1_proteins_symp_plot <- clean_data %>%
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

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/tr1_proteins_symp_plot.png", tr1_proteins_symp_plot, width = 5.33333, height = 5.3333, dpi=444, bg="white")
ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/tr1_proteins_symp_plot.png", tr1_proteins_symp_plot, width =5.3333333, height=5.333333, dpi=444, bg="white")

tr1_proteins_asymp_plot <- clean_data %>%
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

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/tr1_proteins_asymp_plot.png", tr1_proteins_asymp_plot, width =5.3333333, height=5.333333, dpi=444, bg="white")




(ctrl_tr1_plot <- clean_data %>%
  filter(infectiontype %in% c("A"), timepoint!="day28", timepoint!="bad_baseline")%>%
  mutate(day14_para=if_else(timepoint=="day14" & qpcr > 1, "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), "TBS positive\n at day 14", "TBS negative\n at day 14"))%>%
  # filter(targetName %in% sig_ctrl$targetName)%>%
  filter(targetName %in% c("IL10", "LAG3", "GZMA", "IFNG"))%>%#[seq((i-1)*16+1, i*16)])%>%
  filter(!is.na(class2))%>%
  # mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=timepoint, y=concentration, group=interaction(class2, timepoint)))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(aes(fill=timepoint, color=class2), outliers = FALSE)+
  facet_wrap(~targetName, scales = "free_y", nrow=2)+
  scale_fill_manual(values=time_cols)+
  scale_color_manual(values=c("black", "#8B0000"))+
  theme_minimal()+
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_blank()))


ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/tr1_proteins_parasite_ctrl_plot.png", ctrl_tr1_plot, height=5.33, width=5.33, dpi=444, bg="white")


brain_and_bone_stuff <- clean_data %>%
  filter(infectiontype %in% c("A", "S"), timepoint%notin%c("day28","bad_baseline", "day7"))%>%
  mutate(day14_para=if_else(timepoint=="day14" & qpcr > 1, "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), "non_controller", "controller"))%>%
  # filter(targetName %in% sig_ctrl$targetName)%>%
  filter(targetName %in% c("VSNL1", "GFAP", "GDF2", "SPP1"))%>%#[seq((i-1)*16+1, i*16)])%>%
  filter(class2 %in% c("non_controller", "controller"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=factor(timepoint)))+
  geom_line(aes(group=id), alpha=0.2)+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetName, scales = "free", nrow=2)+
  scale_fill_manual(values=time_cols)+
  theme_minimal()+
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.title = element_blank())
da_boxplot_theme

ggsave("~/Downloads/brain_and_bone_stuff.png", width=5.3333, height=5.3333, dpi=444, bg="white")



controller_stuff <- clean_data %>%
  filter(infectiontype %in% c("A", "S"), timepoint%notin%c("day28","bad_baseline", "day7"))%>%
  mutate(day14_para=if_else(timepoint=="day14" & qpcr > 1, "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), "non_controller", "controller"))%>%
  # filter(targetName %in% sig_ctrl$targetName)%>%
  filter(targetName %in% c("IL10", "CRP", "IL15", "CCL22"))%>%#[seq((i-1)*16+1, i*16)])%>%
  filter(class2 %in% c("non_controller", "controller"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=factor(class2)))+
  geom_line(aes(group=id), alpha=0.2)+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetName, scales = "free", nrow=2)+
  # scale_fill_manual(values=time_cols)+
  theme_minimal()+
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.title = element_blank())
da_boxplot_theme

ggsave("~/Downloads/brain_and_bone_stuff.png", width=5.3333, height=5.3333, dpi=444, bg="white")




infection_overview <- clean_data%>%
  distinct(id, timepoint, infectiontype, parasitedensity, new_qpcr)%>%
  filter(timepoint %in% c("baseline", "day0", "day14"),
         infectiontype %in% c("A", "S"))%>%
  mutate(day14_para=if_else(timepoint=="day14" & parasitedensity > 10 &infectiontype=="A", "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), "non_controller", "controller"))



infection_overview%>%
  ggplot(., aes(x=timepoint, y=new_qpcr+0.01, fill=day14_para))+
  geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  theme_minimal()+
  facet_wrap(~infectiontype)




para_plot <- infection_overview%>%
  ggplot(., aes(x=timepoint, y=parasitedensity+0.01, fill=interaction(day14_para, infectiontype)))+
  geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = F)+
  # annotation_logticks(sides="l")+
  facet_wrap(~infectiontype, scales = 'fixed')+
  scale_fill_manual(values=c("lightgrey", "darkgrey", "darkred"))+
  # theme_minimal()+
  ylab("parasites /  μL")+
  scale_y_log10(guide = "axis_logticks")+
  theme(legend.position = "none",
        panel.grid = element_line(color="#f0f0f0"),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank()
  )

ggsave("~/Downloads/para_plot.png", para_plot, width=4.3333, height=4.3333, dpi=444, bg="white")
