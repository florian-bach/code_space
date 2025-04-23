library(tidyr)
library(dplyr)
library(ggplot2)

# read in data

nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/tr1_paper/revised_baseline_clean_musical_combo_with_metadata.csv")

nulisa_data <- nulisa_data%>%
  mutate(timepoint=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))

time_cols <- list("baseline"="#E4DEBD",
                  "day0" = "#C03F3E",
                  "day7" = "#D87E1F",
                  "day14" = "#E6B85F")
# 
infectiontype_cols <- list("A"="#36454F",
                           "S" = "darkred")
long_combo <- read.csv("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/cell_freqs_and_nulisa.csv")
long_combo <-long_combo%>%
  mutate(timepoint=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))




individuals_with_parasites_at_day14 <- long_combo%>%
  distinct(id, infectiontype, timepoint, log_qpcr)%>%
  filter(infectiontype %in% c("A"), timepoint!="day28", timepoint!="bad_baseline")%>%
  mutate(day14_para=if_else(timepoint=="day14" & log_qpcr > 0, "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), "non_controller", "controller"))%>%
  filter(class2=="non_controller")%>%
  distinct(id, infectiontype, class2)

parasitemia_plot <- nulisa_data%>%
  filter(infectiontype%in%c("A", "S"), timepoint!="bad_baseline")%>%
  distinct(id, timepoint, log_qpcr, parasitedensity, infectiontype, new_qpcr)%>%
  mutate(day14_para=if_else(timepoint=="day14" & parasitedensity > 10 &infectiontype=="A", "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), "non_controller", "controller"))%>%
  ggplot(., aes(x=timepoint, y=parasitedensity+0.01, fill=interaction(day14_para, infectiontype)))+
  geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = F)+
  facet_wrap(~infectiontype, scales = 'fixed')+
  scale_fill_manual(values=c("lightgrey", "darkgrey", "darkred"))+
  ylab("parasites /  μL")+
  scale_y_log10(guide = "axis_logticks")+
  theme(legend.position = "none",
        panel.grid = element_line(color="#f0f0f0"),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank()
  )

ggsave("~/postdoc/stanford/abstracts/keystone25/poster_figures/parasitemia_plot.png", parasitemia_plot, height=3, width=6, bg="white", dpi=444)






nulisa_data <- nulisa_data%>%
  mutate(timepoint=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  mutate(infectiontype=substr(infectiontype, 1, 1))

combo_as_results_with_fc <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/tr1_paper/combo_as_results_with_fc.csv")

time_cols <- list("baseline"="#E4DEBD",
                  "day0" = "#C03F3E",
                  "day7" = "#D87E1F",
                  "day14" = "#E6B85F")


base_d0_s_genes <- combo_as_results_with_fc%>%
  filter(padj<0.05, contrast=="baseline S - day0 S", abs(mean_conc_base_d0_fc)>0.5)

base_d7_s_genes <- combo_as_results_with_fc%>%
  filter(padj<0.05, contrast=="baseline S - day7 S", abs(mean_conc_base_d7_fc)>0.5)

base_d0_a_genes <- combo_as_results_with_fc%>%
  filter(padj<0.05, contrast=="baseline A - day0 A", abs(mean_conc_base_d0_fc)>0.5)

base_d14_a_genes <- combo_as_results_with_fc%>%
  filter(padj<0.05, contrast=="baseline A - day14 A", abs(mean_conc_base_d14_fc)>0.5)

base_d7_s_genes$targetName[base_d7_s_genes$targetName%notin% base_d0_s_genes$targetName]
# "CXCL11"  "IL12p70" "IL15"    "SDC1"    "TNFSF9" 
base_d0_a_genes$targetName[base_d0_a_genes$targetName%notin% base_d0_s_genes$targetName]
# nothing
base_d0_a_genes$targetName[base_d0_a_genes$targetName%notin% base_d14_a_genes$targetName]
#"KITLG" "TREM2"
base_d14_a_genes$targetName[base_d14_a_genes$targetName%notin% base_d0_a_genes$targetName]
#"GZMA"   "IL15"   "IL2"    "IL7R"   "LILRB2"
base_d14_a_genes$targetName[base_d14_a_genes$targetName%notin% base_d0_s_genes$targetName]

library(ggvenn)

s_deg_data <- list("base_d0_s"=base_d0_s_genes$targetName,
                   #"base_d0_a"=base_d0_a_genes$targetName,
                   "base_d7_s"=base_d7_s_genes$targetName
                   #"base_d14_a"=base_d14_a_genes$targetName
)
# make data a list (of assyas) of lists (of individuals) 
s_venn_diagram <- ggvenn(
  s_deg_data, 
  set_name_color = NA,
  show_stats = "c",
  text_color = "black",
  text_size = 6,
  fill_color = c("darkred", "orange"),
  stroke_size = 0.5, set_name_size = 4
)

ggsave("~/postdoc/stanford/abstracts/keystone25/poster_figures/day0_s_venn.png", s_venn_diagram, width=5, height=5, bg='white')


a_deg_data <- list(#"base_d0_s"=base_d0_s_genes$targetName,
  "base_d0_a"=base_d0_a_genes$targetName,
  #"base_d7_s"=base_d7_s_genes$targetName
  "base_d14_a"=base_d14_a_genes$targetName
)
# make data a list (of assyas) of lists (of individuals)
a_venn_diagram <- ggvenn(
  a_deg_data, set_name_color = NA,
  show_stats = "c",
  text_color = "black",
  text_size = 6,
  fill_color = unname(time_cols[c(1,4)]),
  stroke_size = 0.5, set_name_size = 4
)
ggsave("~/postdoc/stanford/abstracts/keystone25/poster_figures/day14_a_venn.png", a_venn_diagram, width=5, height=5, bg='white')



unique_day14_plot <- nulisa_data %>%
  mutate(day14_para=if_else(timepoint=="day14" & new_parasitedensity > 10 & infectiontype=="A", "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id, infectiontype)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), ">10 parasites / μL on day 14", "no parasites day14"))%>%
  filter(targetName %in% base_d14_a_genes$targetName[base_d14_a_genes$targetName%notin% base_d0_a_genes$targetName],
         infectiontype%in%c("A"),
         timepoint %in% names(time_cols),
         !is.na(class2))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, color=class2))+
  geom_boxplot(aes(fill=timepoint), outliers = FALSE)+
  # ggpubr::stat_compare_means()+
  facet_wrap(~targetName+infectiontype, scales = "free", nrow=1)+
  scale_fill_manual(values=time_cols)+
  scale_color_manual(values=c("#8B0000", "black"))+
  # scale_fill_manual(values=c("grey", "darkred"))+
  theme_minimal()+
  guides(fill=guide_none())+
  theme(legend.position = "bottom",
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/abstracts/keystone25/poster_figures/unique_day14_plot.png", unique_day14_plot, width=12, height=4, bg='white')


unique_day7_plot <- nulisa_data %>%
  mutate(day14_para=if_else(timepoint=="day14" & new_parasitedensity > 1 & infectiontype=="A", "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id, infectiontype)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), ">1 parasites / μL on day 14", "no parasites day14"))%>%
  # filter(targetName %in% base_d7_s_genes$targetName,
  filter(targetName %in% base_d7_s_genes$targetName[base_d7_s_genes$targetName%notin% base_d0_s_genes$targetName],
         infectiontype%in%c("S"),
         timepoint %in% names(time_cols))%>%
  ggplot(aes(x=factor(timepoint), y=concentration))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(aes(fill=timepoint), outliers = FALSE)+
  # ggpubr::stat_compare_means()+
  facet_wrap(~targetName+infectiontype, scales = "free", nrow=1)+
  scale_fill_manual(values=time_cols)+
  # scale_fill_manual(values=c("grey", "darkred"))+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        legend.position = "none")

ggsave("~/postdoc/stanford/abstracts/keystone25/poster_figures/unique_day7_plot.png", unique_day7_plot, width=12, height=4, bg='white')



nulisa_data %>%
  # mutate(day14_para=if_else(timepoint=="day14" & new_parasitedensity > 1 & infectiontype=="A", "parasitemic_day14", "no_parasites_day14"))%>%
  # group_by(id, infectiontype)%>%
  # mutate(class2= if_else(any(day14_para=="parasitemic_day14"), ">1 parasites / μL on day 14", "no parasites day14"))%>%
  filter(targetName %in% c("CXCL9", "CXCL11"),
         infectiontype%in%c("S", "A"),
         timepoint %in% names(time_cols))%>%
  pivot_wider(names_from = targetName, values_from = concentration, id_cols = c(id, timepoint, infectiontype))%>%
  ggplot(., aes(x=CXCL9, y=CXCL11))+
  geom_point(aes(color=timepoint), outliers = FALSE)+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method = "spearman")+
  facet_wrap(~infectiontype, scales = "fixed")+
  scale_color_manual(values=time_cols)+
  # scale_fill_manual(values=c("grey", "darkred"))+
  theme_minimal()

# plots about parasite controllers ####

ctrl_fc_data <- nulisa_data %>%
  filter(infectiontype%in%c("A"))%>%
  mutate(day14_para=if_else(timepoint=="day14" & new_parasitedensity > 10 & infectiontype=="A", "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id, infectiontype)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), ">10 parasites / μL on day 14", "no parasites day14"))%>%
  filter(class2==">10 parasites / μL on day 14")%>%
  pivot_wider(names_from = timepoint, values_from = c(concentration), id_cols=c(targetName, id, infectiontype, class2), names_prefix = "conc_")%>%
  group_by(targetName)%>%
  mutate(conc_base_d0_fc=conc_day0-conc_baseline,
         conc_base_d14_fc=conc_day14-conc_baseline)%>%
  distinct(targetName, id, infectiontype, conc_base_d0_fc, conc_base_d14_fc)%>%
  group_by(targetName, infectiontype)%>%
  summarise("mean_conc_base_d0_fc"=mean(conc_base_d0_fc, na.rm = T),
            "mean_conc_base_d14_fc"=mean(conc_base_d14_fc, na.rm = T))


ctrl_a_results_with_fc <- day14_a_purff%>%
  left_join(., ctrl_fc_data, by="targetName")%>%
  distinct(targetName, contrast, padj, infectiontype, mean_conc_base_d0_fc, mean_conc_base_d14_fc)

(day14_ctrl_volcano <- ctrl_a_results_with_fc %>%
    filter(contrast == "baseline controller - day14 controller p")%>%
    distinct(targetName, padj, infectiontype, mean_conc_base_d14_fc)%>%
    filter(infectiontype=="A", is.finite(mean_conc_base_d14_fc))%>%
    mutate("label2" = if_else(targetName %in% c(sig_non_ctrl_base_14$targetName), targetName, NA))%>%
    ggplot(., aes(x=mean_conc_base_d14_fc, y=-log10(padj+10^-14), alpha=padj<0.05&abs(mean_conc_base_d14_fc)>0.5, color=mean_conc_base_d14_fc<0))+
    geom_point()+
    ggrepel::geom_text_repel(aes(label=label2, alpha=padj<0.05&abs(mean_conc_base_d14_fc)>0.5),
                             size = 5,min.segment.length = 0, ,
                             position=ggpp::position_nudge_center(center_x = 0, x = 1,
                                                                  center_y = 2.5, y=0))+
    ggtitle("n = 0")+
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

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/day14_ctrl_volcano.png", day14_ctrl_volcano, width = 5, height = 5, dpi=444, bg="white")



(day0_ctrl_volcano <- no_ctrl_a_results_with_fc %>%
    filter(contrast == "baseline controller - day0 controller p")%>%
    distinct(targetName, padj, infectiontype, mean_conc_base_d0_fc)%>%
    filter(infectiontype=="A", is.finite(mean_conc_base_d0_fc))%>%
    mutate("label2" = if_else(targetName %in% c(sig_ctrl_base_0$targetName), targetName, NA))%>%
    ggplot(., aes(x=mean_conc_base_d0_fc, y=-log10(padj+10^-14), alpha=padj<0.05&abs(mean_conc_base_d0_fc)>0.5, color=mean_conc_base_d0_fc<0))+
    geom_point()+
    ggrepel::geom_text_repel(aes(label=label2, alpha=padj<0.05&abs(mean_conc_base_d0_fc)>0.5),
                             size = 5,min.segment.length = 0, ,
                             position=ggpp::position_nudge_center(center_x = 0, x = 1,
                                                                  center_y = 2.5, y=0))+
    ggtitle("n = 20")+
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

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/day0_ctrl_volcano.png", day0_ctrl_volcano, width = 5, height = 5, dpi=444, bg="white")


# plots about parasite non controllers ####


no_ctrl_fc_data <- nulisa_data %>%
  filter(infectiontype%in%c("A"))%>%
  mutate(day14_para=if_else(timepoint=="day14" & new_parasitedensity > 10 & infectiontype=="A", "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id, infectiontype)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), ">10 parasites / μL on day 14", "no parasites day14"))%>%
  filter(class2=="no parasites day14")%>%
  pivot_wider(names_from = timepoint, values_from = c(concentration), id_cols=c(targetName, id, infectiontype, class2), names_prefix = "conc_")%>%
  group_by(targetName)%>%
  mutate(conc_base_d0_fc=conc_day0-conc_baseline,
         conc_base_d14_fc=conc_day14-conc_baseline)%>%
  distinct(targetName, id, infectiontype, conc_base_d0_fc, conc_base_d14_fc)%>%
  group_by(targetName, infectiontype)%>%
  summarise("mean_conc_base_d0_fc"=mean(conc_base_d0_fc, na.rm = T),
            "mean_conc_base_d14_fc"=mean(conc_base_d14_fc, na.rm = T))


no_ctrl_a_results_with_fc <- day14_a_purff%>%
  left_join(., no_ctrl_fc_data, by="targetName")%>%
  distinct(targetName, contrast, padj, infectiontype, mean_conc_base_d0_fc, mean_conc_base_d14_fc)


(day14_no_ctrl_volcano <- no_ctrl_a_results_with_fc %>%
    filter(contrast == "baseline non_controller - day14 non_controller p")%>%
    distinct(targetName, padj, infectiontype, mean_conc_base_d14_fc)%>%
    filter(infectiontype=="A", is.finite(mean_conc_base_d14_fc))%>%
    mutate("label2" = if_else(targetName %in% c(sig_non_ctrl_base_14$targetName), targetName, NA))%>%
    ggplot(., aes(x=mean_conc_base_d14_fc, y=-log10(padj+10^-14), alpha=padj<0.05&abs(mean_conc_base_d14_fc)>0.5, color=mean_conc_base_d14_fc<0))+
    geom_point()+
    ggrepel::geom_text_repel(aes(label=label2, alpha=padj<0.05&abs(mean_conc_base_d14_fc)>0.5),
                             size = 5,min.segment.length = 0, ,
                             position=ggpp::position_nudge_center(center_x = 0, x = 1,
                                                                  center_y = 2.5, y=0))+
    ggtitle("n = 20")+
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

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/day14_ctrl_volcano.png", day14_ctrl_volcano, width = 5, height = 5, dpi=444, bg="white")



(day0_no_ctrl_volcano <- no_ctrl_a_results_with_fc %>%
    filter(contrast == "baseline non_controller - day0 non_controller p")%>%
    distinct(targetName, padj, infectiontype, mean_conc_base_d0_fc)%>%
    filter(infectiontype=="A", is.finite(mean_conc_base_d0_fc))%>%
    mutate("label2" = if_else(targetName %in% c(sig_non_ctrl_base_0$targetName), targetName, NA))%>%
    ggplot(., aes(x=mean_conc_base_d0_fc, y=-log10(padj+10^-14), alpha=padj<0.05&abs(mean_conc_base_d0_fc)>0.5, color=mean_conc_base_d0_fc<0))+
    geom_point()+
    ggrepel::geom_text_repel(aes(label=label2, alpha=padj<0.05&abs(mean_conc_base_d0_fc)>0.5),
                             size = 5,min.segment.length = 0, ,
                             position=ggpp::position_nudge_center(center_x = 0, x = 1,
                                                                  center_y = 2.5, y=0))+
    ggtitle("n = 20")+
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

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/day0_no_ctrl_volcano.png", day0_no_ctrl_volcano, width = 5, height = 5, dpi=444, bg="white")





# plots about parasite controllers vs non_ontrollers ####

ctrl_vs_nonctrl_base14_fc_data <- nulisa_data %>%
  filter(infectiontype%in%c("A"), timepoint=="day14")%>%
  mutate(day14_para=if_else(timepoint=="day14" & new_parasitedensity > 10 & infectiontype=="A", "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id, infectiontype)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), ">10 parasites / μL on day 14", "no parasites day14"))%>%
  pivot_wider(names_from = class2, values_from = c(concentration), id_cols=c(targetName, id), names_prefix = "conc_")%>%
  group_by(targetName)%>%
  mutate(ctrl_vs_no_fc=mean(`conc_>10 parasites / μL on day 14`, na.rm=T)-mean(`conc_no parasites day14`, na.rm=T))%>%
  distinct(targetName, ctrl_vs_no_fc)


ctrl_vs_nonctrl_base14_results_with_fc <- day14_a_purff%>%
  left_join(., ctrl_vs_nonctrl_base14_fc_data, by="targetName")%>%
  distinct(targetName, contrast, padj, ctrl_vs_no_fc)


(day14_ctrl_no_ctrl_volcano <- ctrl_vs_nonctrl_base14_results_with_fc %>%
    filter(contrast == "controller - non_controller day14 p")%>%
    distinct(targetName, padj, ctrl_vs_no_fc)%>%
    filter(is.finite(ctrl_vs_no_fc))%>%
    mutate("label2" = if_else(targetName %in% c(sig_day14_control$targetName), targetName, NA))%>%
    ggplot(., aes(x=ctrl_vs_no_fc, y=-log10(padj+10^-14), alpha=padj<0.05&abs(ctrl_vs_no_fc)>0.5, color=ctrl_vs_no_fc<0))+
    geom_point()+
    ggrepel::geom_text_repel(aes(label=label2, alpha=padj<0.05&abs(ctrl_vs_no_fc)>0.5),
                             size = 5,min.segment.length = 0, ,
                             position=ggpp::position_nudge_center(center_x = 0, x = 1,
                                                                  center_y = 2.5, y=0))+
    ggtitle("n = 31")+
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

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/day14_ctrl_volcano.png", day14_ctrl_no_ctrl_volcano, width = 5, height = 5, dpi=444, bg="white")

# other plots ####

unique_day14_ctrl_plot <- nulisa_data %>%
  mutate(day14_para=if_else(timepoint=="day14" & new_parasitedensity > 10 & infectiontype=="A", "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id, infectiontype)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), ">10 parasites / μL on day 14", "no parasites day14"))%>%
  filter(targetName %in% c("LAG3", "GZMA", "IL10", "TNF", "IL27"),
         infectiontype%in%c("A"),
         timepoint %in% names(time_cols),
         !is.na(class2))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, color=class2))+
  geom_boxplot(aes(fill=timepoint), outliers = FALSE)+
  # ggpubr::stat_compare_means()+
  facet_wrap(~targetName, scales = "free", nrow=1)+
  scale_fill_manual(values=time_cols)+
  scale_color_manual(values=c("#8B0000", "black"))+
  # scale_fill_manual(values=c("grey", "darkred"))+
  theme_minimal()+
  guides(fill=guide_none())+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/abstracts/keystone25/poster_figures/unique_day14_ctrl_plot.png", unique_day14_ctrl_plot, width=12, height=4, bg='white', dpi=444)


# freq and count plots ####

library(patchwork)
cell_count_data <- read.csv("~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/t_cell_cytometry_count_data.csv")

s_tr1_count_plot <- cell_count_data %>%
  filter(infectiontype %in% c("S"), timepoint!="bad_baseline", stim=="unstim")%>%
  distinct(id, infectiontype, timepoint, absolute_Tr1)%>%
  ggplot(., aes(x=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")), y=absolute_Tr1, fill=timepoint))+
  geom_boxplot(outliers = F)+
  # scale_y_log10()+
  ylab("Tr1 cells / μL")+
  scale_fill_manual(values=time_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")


s_tr1_freq_plot <- cell_count_data %>%
  filter(infectiontype %in% c("S"), timepoint!="bad_baseline", gate=="Tr1_Frequency", stim=="unstim")%>%
  distinct(id, infectiontype, timepoint, freq)%>%
  ggplot(., aes(x=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")), y=freq/100, fill=timepoint))+
  geom_boxplot(outliers = F)+
  scale_y_continuous(labels = scales::label_percent())+
  ylab("Tr1 cells of non-naive CD4")+
  scale_fill_manual(values=time_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")



s_tr1_freq_and_count <- s_tr1_freq_plot + s_tr1_count_plot

ggsave("~/postdoc/stanford/abstracts/keystone25/poster_figures/s_tr1_freq_and_count.png", s_tr1_freq_and_count, width = 6, height = 4, dpi=444, bg="white")



a_tr1_count_plot <- cell_count_data %>%
  filter(infectiontype %in% c("A"), timepoint!="bad_baseline", stim=="unstim")%>%
  mutate(class2=case_when(infectiontype=="A"&id%in%individuals_with_parasites_at_day14$id~"parasites at day 14",
                          infectiontype=="A"&id%notin%individuals_with_parasites_at_day14$id~"no parasites at day 14", .default = NULL))%>%
  distinct(id, infectiontype, timepoint, absolute_Tr1, class2)%>%
  ggplot(., aes(x=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")), y=absolute_Tr1, fill=timepoint, color=class2))+
  geom_boxplot(outliers = F)+
  # scale_y_log10()+
  ylab("Tr1 cells / μL")+
  scale_fill_manual(values=time_cols)+
  scale_color_manual(values=c("black", "darkred"))+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")


a_tr1_freq_plot <- cell_count_data %>%
  filter(infectiontype %in% c("A"), timepoint!="bad_baseline", gate=="Tr1_Frequency", stim=="unstim")%>%
  mutate(class2=case_when(infectiontype=="A"&id%in%individuals_with_parasites_at_day14$id~"parasites at day 14",
                          infectiontype=="A"&id%notin%individuals_with_parasites_at_day14$id~"no parasites at day 14", .default = NULL))%>%
  distinct(id, infectiontype, timepoint, freq, class2)%>%
  ggplot(., aes(x=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")), y=freq/100, fill=timepoint, color=class2))+
  geom_boxplot(outliers = F)+
  scale_y_continuous(labels = scales::label_percent())+
  ylab("Tr1 cells of non-naive CD4")+
  scale_fill_manual(values=time_cols)+
  scale_color_manual(values=c("black", "darkred"))+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")



a_tr1_freq_and_count <- a_tr1_freq_plot + a_tr1_count_plot

ggsave("~/postdoc/stanford/abstracts/keystone25/poster_figures/a_tr1_freq_and_count.png", a_tr1_freq_and_count, width = 6, height = 4, dpi=444, bg="white")
