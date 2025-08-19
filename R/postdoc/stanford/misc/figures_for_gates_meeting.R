library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(emmeans)

`%notin%` <- Negate(`%in%`)
treatment_palette <- c("Placebo"="darkred", "DP 1 year"="#00555A", "DP 2 years"="orange")
# var gene stuff ####
long_luminex <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/long_luminex.csv")%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks", "104 weeks")))

luminex_fc <- long_luminex %>%
  group_by(antigen, timepoint, treatmentarm)%>%
  summarise("mean_conc"=mean(MFI, na.rm = T))%>%
  pivot_wider(names_from = treatmentarm, values_from = mean_conc)%>%
  mutate(fc=log2(Placebo/`DP 1 year`))

var_gene_treatment_purf <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/var_gene_treatment_regression.csv")

luminex_purf_fc <- var_gene_treatment_purf%>%
  mutate("timepoint"=contrast)%>%
  left_join(., luminex_fc, by=c("antigen", "timepoint"))%>%
  mutate(contrast=factor(contrast, levels=c("8 weeks", "24 weeks", "52 weeks", "104 weeks")))


sig_treatment_varying_vars <- luminex_purf_fc %>% 
  filter(padj<0.1, abs(fc)>=0.5)%>%
  arrange(padj)


# treatment_varying_vars <- c("schizont", "CIDRa1.7 2083-1", "CIDRg3 IT4var08", "CIDRa2.2 IT4var24", "CIDRa1.4 HB3var03", 
# "IT4var02", "CIDRg3", "CIDRa2.7 D3 IT4var61", "CIDRa1.8a 1965-3 BT", "CIDRa1.6b IT4var18")  

## var gene volcano ####


(var_treatment_volcano <- luminex_purf_fc %>%
   filter(contrast %in% c("8 weeks", "24 weeks", "52 weeks", "104 weeks"))%>%
   arrange(padj)%>%
   mutate("label2" = if_else(antigen%in%c("schizont", "CIDRa1.7 2083-1", "CIDRg3 IT4var08") & timepoint %in% c("52 weeks"), antigen,
                             if_else(antigen=="CIDRa6"&timepoint=="24 weeks", "CIDRa6", NA)))%>%
   # mutate("label2" = if_else(.$targetName, .$contrast)== cbind(treatment_nulisa_sigs$targetName, treatment_nulisa_sigs$contrast)), targetName, NA))%>%
   ggplot(., aes(x=-fc, y=padj, alpha=padj<0.1&abs(fc)>=0.5, color=fc<0))+
   geom_hline(yintercept = 0.1, linetype="dashed", alpha = 0.5)+
   geom_vline(xintercept = -0.5, linetype="dashed", alpha = 0.5)+
   geom_vline(xintercept = 0.5, linetype="dashed", alpha = 0.5)+
   geom_point()+
   ggrepel::geom_text_repel(aes(label=label2, alpha=padj<0.1&abs(fc)>=0.5),
                            size = 4, force = 20,
                            position=ggpp::position_nudge_center(center_x = 0, x = 0.5,
                                                                 center_y = 0.01, y=0.5)
   )+
   scale_alpha_manual(values=c(0.5, 1))+
   scale_color_manual(values=rev(c("#ff8c00", "darkblue")))+
   scale_x_continuous(limits=c(-2.5, 2.5))+
   scale_y_continuous(trans=c("log10", "reverse"))+
   # ggtitle("differentially abundant plasma proteins in placebo vs. DP")+
   facet_wrap(~contrast, nrow=1)+
   xlab("log2 fold change")+
   ylab("padj")+
   theme_minimal()+
   theme(legend.position="none"))


ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/sig_var_treatment_volcano_plot.png", var_treatment_volcano, height=4, width=8, dpi=400, bg="white")


## var gene boxplots ####

treatment_vars52 <- long_luminex %>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  filter(antigen %in% sig_treatment_varying_vars$antigen[1:3])%>% 
  filter(timepoint %in% c("52 weeks"))%>%
  ggplot(., aes(x=factor(timepoint), y=MFI, fill=treatmentarm))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  facet_wrap(~factor(antigen, levels=unique(sig_treatment_varying_vars$antigen)), scales="free", labeller = label_wrap_gen(width = 10), nrow=1)+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  # viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  scale_fill_manual(values=treatment_palette)+
  ggpubr::stat_compare_means(size=2, vjust = 0.2)+
  theme_minimal()+
  ylab("MFI at 52 weeks")+
  theme(legend.position="bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text.x = element_text(size=7),
        axis.title.x=element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/sig_var_gene_treatment_vars52.png", treatment_vars52, height=4, width=8, dpi=400, bg="white")  



treatment_vars104 <- long_luminex %>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  filter(antigen %in%c("schizont", "CIDRa1.7 2083-1", "CIDRg3 IT4var08"))%>% 
  filter(timepoint %in% c("104 weeks"))%>%
  ggplot(., aes(x=factor(timepoint), y=MFI, fill=treatmentarm))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  facet_wrap(~factor(antigen, levels=c("schizont", "CIDRa1.7 2083-1", "CIDRg3 IT4var08")), scales="free", labeller = label_wrap_gen(width = 10), nrow=1)+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  # ggpubr::stat_compare_means(size=2, vjust = 0.2)+
  # viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  scale_fill_manual(values=treatment_palette)+
  theme_minimal()+
  ylab("MFI at 104 weeks")+
  theme(legend.position="bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text.x = element_text(size=7),
        axis.title.x=element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/sig_var_gene_treatment_vars104.png", treatment_vars104, height=4, width=8, dpi=400, bg="white")  

inf_sigs_52 <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/inf_sigs_52.csv")

inf_varying_vars <- inf_sigs_52%>%
  slice_min(n_para_padj, n=4, with_ties = F)%>%
  filter(antigen !="CIDRa1.4 HB3var03")
  

malaria_vars_placebo <- long_luminex %>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  # filter(antigen %in% treatment_varying_vars[c(1:3, 5)])%>% 
  filter(antigen %in% inf_varying_vars$antigen)%>% 
  filter(timepoint %in% c("52 weeks"))%>%
  mutate(total_n_malaria_12=ifelse(total_n_malaria_12>=2, "2+", total_n_malaria_12))%>%
  ggplot(., aes(x=factor(total_n_malaria_12), y=MFI, fill=treatmentarm))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  # geom_smooth(aes(x=total_n_para_12, y=MFI), method="lm")+
  facet_wrap(~factor(antigen, levels=inf_varying_vars$antigen), scales="free", labeller = label_wrap_gen(width = 10), nrow=1)+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  # viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  scale_fill_manual(values=treatment_palette)+
  xlab("number of malaria episodes in the first year of life")+
  ylab("MFI at 52 weeks")+
  theme_minimal()+
  theme(legend.position="none",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text.x = element_text(size=7))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/sig_malaria_vars.png", malaria_vars_placebo, height=4, width=8, dpi=400, bg="white")  



malaria_vars_placebo2 <- long_luminex %>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  # filter(antigen %in% treatment_varying_vars[c(1:3, 5)])%>% 
  filter(antigen %in% inf_varying_vars$antigen[c(1:8)])%>% 
  filter(timepoint %in% c("104 weeks"))%>%
  mutate(total_n_malaria_12_24=ifelse(total_n_malaria_12_24>=2, "2+", total_n_malaria_12_24))%>%
  ggplot(., aes(x=factor(total_n_malaria_12_24), y=MFI, fill=treatmentarm))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  # geom_smooth(aes(x=total_n_para_12, y=MFI), method="lm")+
  facet_wrap(~factor(antigen, levels=inf_varying_vars$antigen), scales="free", labeller = label_wrap_gen(width = 10), nrow=2)+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  # viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  scale_fill_manual(values=treatment_palette)+
  xlab("number of malaria episodes in the second year of life")+
  ylab("MFI at 104 weeks")+
  theme_minimal()+
  theme(legend.position="none",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text.x = element_text(size=7))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/sig_malaria_vars2.png", malaria_vars_placebo2, height=4, width=8, dpi=400, bg="white")  




parasitemia_vars_placebo <- long_luminex %>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  mutate(total_n_para_12=ifelse(total_n_para_12<5, total_n_para_12, "5+"))%>%
  filter(antigen %in% inf_varying_vars$antigen)%>% 
  filter(timepoint %in% c("52 weeks"))%>%
  ggplot(., aes(x=factor(total_n_para_12), y=MFI, fill=treatmentarm))+
  stat_boxplot(geom = 'errorbar', color="black", width = 0,
               position = position_dodge(0.75), linewidth = 0.5) +
  geom_boxplot(outliers = F, coef=0, color="white", lwd=0.1, size = 4)+
  scale_y_log10()+
  # geom_smooth(aes(x=total_n_para_12, y=MFI), method="lm")+
  facet_wrap(~factor(antigen, levels=inf_varying_vars$antigen), scales="free", labeller = label_wrap_gen(width = 10), nrow=1)+
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  # viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  scale_fill_manual(values=treatment_palette)+
  xlab("number of parasitemic episodes in the first year of life")+
  ylab("MFI at 52 weeks")+
  theme_minimal()+
  theme(legend.position="none",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text.x = element_text(size=7))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/parasitemia_vars.png", parasitemia_vars_placebo, height=4, width=8, dpi=400, bg="white")  



n_para_var_select104 <- long_luminex %>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  mutate(total_n_para_12_24=ifelse(total_n_para_12_24<4, total_n_para_12_24, "4+"))%>%
  filter(antigen %in% c("schizont", "IT4var02", "CIDRa6"))%>% 
  filter(timepoint %in% c("104 weeks"))%>%
  ggplot(., aes(x=factor(total_n_para_12_24), y=MFI, fill=treatmentarm))+
  stat_boxplot(geom = 'errorbar', color="black", width = 0,
               position = position_dodge(0.75), linewidth = 0.5) +
  geom_boxplot(outliers = F, coef=0, color="white", lwd=0.1, size = 4)+
  scale_y_log10()+
  # geom_smooth(aes(x=total_n_para_12, y=MFI), method="lm")+
  facet_wrap(~factor(antigen,c("CIDRa6", "schizont", "IT4var02")), scales="free", labeller = label_wrap_gen(width = 10), nrow=1)+
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  # viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  scale_fill_manual(values=treatment_palette)+
  xlab("number of parasitemic episodes in the second year of life")+
  ylab("MFI at 104 weeks")+
  theme_minimal()+
  theme(legend.position="none",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text.x = element_text(size=7))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/n_para_var_select104.png", n_para_var_select104, height=4, width=8, dpi=400, bg="white")  



time_var_select52_104 <- long_luminex %>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  mutate(total_n_para_12_24=ifelse(total_n_para_12_24<4, total_n_para_12_24, "4+"))%>%
  filter(antigen %in% c("schizont", "IT4var02", "CIDRa6"))%>% 
  filter(timepoint %in% c("52 weeks", "104 weeks"))%>%
  ggplot(., aes(x=factor(timepoint), y=MFI, fill=treatmentarm))+
  stat_boxplot(geom = 'errorbar', color="black", width = 0,
               position = position_dodge(0.75), linewidth = 0.5) +
  geom_boxplot(outliers = F, coef=0, color="white", lwd=0.1, size = 4)+
  scale_y_log10()+
  # geom_smooth(aes(x=total_n_para_12, y=MFI), method="lm")+
  facet_wrap(~factor(antigen,c("CIDRa6", "schizont", "IT4var02")), scales="free", labeller = label_wrap_gen(width = 10), nrow=1)+
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  # viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  scale_fill_manual(values=treatment_palette)+
  xlab("number of parasitemic episodes in the second year of life")+
  ylab("MFI at 104 weeks")+
  theme_minimal()+
  theme(legend.position="none",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text.x = element_text(size=7))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/time_var_select52_104.png", time_var_select52_104, height=4, width=8, dpi=400, bg="white")  

# 
# malaria_vars_dp <- long_luminex %>%
#   mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
#   filter(antigen %in% treatment_varying_vars[c(1:3, 5)])%>% 
#   filter(timepoint %in% c("52 weeks"))%>%
#   ggplot(., aes(x=factor(total_n_malaria_12), y=MFI, fill=treatmentarm))+
#   geom_boxplot(outliers = F)+
#   scale_y_log10()+
#   # geom_smooth(aes(x=total_n_para_12, y=MFI), method="lm")+
#   facet_wrap(~factor(antigen, levels=treatment_varying_vars), scales="free", labeller = label_wrap_gen(width = 10), nrow=1)+
#   stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
#                geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
#   # viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
#   scale_fill_manual(values=treatment_palette)+
#   xlab("number of malaria episodes")+
#   theme_minimal()+
#   theme(legend.position="none",
#         legend.direction = "horizontal",
#         legend.title = element_blank(),
#         axis.text.x = element_text(size=7))
# 
# ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/malaria_vars_dp.png", malaria_vars_dp, height=2, width=8, dpi=400, bg="white")  
# 
# 
# 
# parasitemia_vars_dp <- long_luminex %>%
#   mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
#   # filter(antigen %in% treatment_varying_vars[c(1:3, 5)])%>% 
#   filter(timepoint %in% c("104 weeks"), treatmentarm=="DP 1 year")%>%
#   ggplot(., aes(x=factor(total_n_para_12_24), y=MFI, fill=treatmentarm))+
#   geom_boxplot(outliers = F)+
#   scale_y_log10()+
#   # geom_smooth(aes(x=total_n_para_12, y=MFI), method="lm")+
#   facet_wrap(~factor(antigen), scales="free", labeller = label_wrap_gen(width = 10))+
#   stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
#                geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
#   # viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
#   scale_fill_manual(values=treatment_palette)+
#   xlab("number of parasitemic episodes")+
#   theme_minimal()+
#   theme(legend.position="none",
#         legend.direction = "horizontal",
#         legend.title = element_blank(),
#         axis.text.x = element_text(size=7))
# 
# ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/parasitemia_vars_dp.png", parasitemia_vars_dp, height=2, width=8, dpi=400, bg="white")  
# 


# nulisa stuff ####
clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks", "68 weeks")))%>%
  filter(targetName %notin% c("CTSS", "LTA|LTB", "IFNA2"))

mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")

clean_data <- clean_data%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))

## treatment effect ####
treatment_nulisa_purf <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/nulisa_treatment_regression.csv")

treatment_nulisa_fc_data <- clean_data %>%
  filter(treatmentarm != "DP 2 years", timepoint !=c("68 weeks"))%>%
  group_by(targetName, timepoint, treatmentarm)%>%
  summarise("mean_conc"=mean(conc, na.rm = T))%>%
  pivot_wider(names_from = treatmentarm, values_from = mean_conc)%>%
  mutate(fc=Placebo-`DP 1 year`)


treatment_nulisa_purf_fc <- treatment_nulisa_purf%>%
  mutate("timepoint"=contrast)%>%
  left_join(., treatment_nulisa_fc_data, by=c("targetName", "timepoint"))%>%
  mutate(contrast=factor(contrast, levels=c("8 weeks", "24 weeks", "52 weeks")))

treatment_nulisa_sigs <- treatment_nulisa_purf_fc %>%
  filter(padj<0.1 & abs(fc)>=0.5 )

treatment_nulisa_kinda_sigs <- treatment_nulisa_purf_fc%>%
  filter(padj<0.05 )



## volcano ####
(treatment_volcano <- treatment_nulisa_purf_fc %>%
    mutate("label2" = if_else(targetName %in% c(treatment_nulisa_sigs$targetName) & timepoint != "8 weeks", targetName, NA))%>%
    # mutate("label2" = if_else(.$targetName, .$contrast)== cbind(treatment_nulisa_sigs$targetName, treatment_nulisa_sigs$contrast)), targetName, NA))%>%
    ggplot(., aes(x=-fc, y=padj, alpha=padj<0.1&abs(fc)>=0.5, color=fc<0))+
    geom_hline(yintercept = 0.1, linetype="dashed", alpha = 0.5)+
    geom_vline(xintercept = -0.5, linetype="dashed", alpha = 0.5)+
    geom_vline(xintercept = 0.5, linetype="dashed", alpha = 0.5)+
   
    geom_point()+
    ggrepel::geom_text_repel(aes(label=label2, alpha=padj<0.1&abs(fc)>=0.5),
                             size = 4,min.segment.length = 0,
                             position=ggpp::position_nudge_center(center_x = 0, x = 0.8,
                                                                  center_y = 0.1, y=2.5)
    )+
    scale_alpha_manual(values=c(0.5, 1))+
    scale_color_manual(values=rev(c("#FF8C00", "darkblue")))+
    scale_x_continuous(limits=c(-2, 2))+
    scale_y_continuous(trans=c("log10", "reverse"))+
    # ggtitle("differentially abundant plasma proteins in placebo vs. DP")+
    facet_wrap(~contrast)+
    xlab("log2 fold change")+
    ylab("padj")+
    theme_minimal()+
    theme(legend.position="none"))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/sig_nulisa_treatment_volcano_plot.png", treatment_volcano, height=4, width=10, dpi=400, bg="white")


treatment_nulisa_purf <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/nulisa_treatment_regression.csv")

treatment_nulisa_sigs <- treatment_nulisa_purf%>%
  filter(padj<0.1)%>%
  filter(contrast=="52 weeks")

treatment_plot <- clean_data%>%
  filter(treatmentarm!="DP 2 years", timepoint=="52 weeks", mstatus==0)%>%
  filter(targetName %in% treatment_nulisa_sigs$targetName)%>%
  ggplot(aes(x=factor(timepoint), y=conc, fill=treatmentarm))+
  # geom_violin(draw_quantiles = seq(0,1,0.25), color="white")+
  # ggpubr::stat_compare_means(size=2, )+
  geom_boxplot(outliers=F)+
  # # stat_summary(aes(group=c(treatmentarm)), geom="cro", fun.y=median, colour="white")+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  facet_wrap(~targetName, scales="free", nrow=2)+
  scale_fill_manual(values=treatment_palette)+
  theme_minimal()+
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")


ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/sig_nulisa_treatment_plot.png", treatment_plot, height=4, width=8, dpi=400, bg="white")



## infection history ####
n_infection_purf <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/n_infection_purf.csv")

sig_mala_12 <- n_infection_purf%>%
  filter(n_malaria_padj<0.1)

sig_para_12 <- n_infection_purf%>%
  filter(n_para_padj<0.1)

n_para_nulisa <- clean_data%>%
  filter(treatmentarm!="DP 2 years", timepoint=="52 weeks", mstatus==0)%>%
  filter(targetName %in% sig_para_12$targetName, !(total_n_para_12==3&treatmentarm=="DP 1 year"))%>%
  mutate(total_n_para_12=ifelse(total_n_para_12<7, total_n_para_12, "7+"))%>%
  ggplot(aes(x=factor(total_n_para_12), y=conc, fill=factor(treatmentarm)))+
  # geom_violin(draw_quantiles = seq(0,1,0.25), color="white")+
  # ggpubr::stat_compare_means(size=2, )+
  geom_boxplot(outliers=F)+
  # # stat_summary(aes(group=c(treatmentarm)), geom="cro", fun.y=median, colour="white")+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  facet_wrap(~targetName, scales="free", nrow=2)+
  scale_fill_manual(values=treatment_palette)+
  theme_minimal()+
  xlab("number of parasitemic episodes in first year of life")+
  theme(legend.title = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.direction = "horizontal")

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/n_para_nulisa_plot.png", n_para_nulisa, height=4, width=8, dpi=400, bg="white")


n_para_nulisa_select <- clean_data%>%
  filter(treatmentarm!="DP 2 years", timepoint=="52 weeks", mstatus==0)%>%
  filter(targetName %in% c("IL10", "TLR3", "IL15"), !(total_n_para_12==3&treatmentarm=="DP 1 year"))%>%
  mutate(total_n_para_12=ifelse(total_n_para_12<7, total_n_para_12, "7+"))%>%
  ggplot(aes(x=factor(total_n_para_12), y=conc, fill=factor(treatmentarm)))+
  # geom_violin(draw_quantiles = seq(0,1,0.25), color="white")+
  # ggpubr::stat_compare_means(size=2, )+
  geom_boxplot(outliers=F)+
  # # stat_summary(aes(group=c(treatmentarm)), geom="cro", fun.y=median, colour="white")+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  facet_wrap(~targetName, scales="free", nrow=1)+
  scale_fill_manual(values=treatment_palette)+
  theme_minimal()+
  xlab("number of parasitemic episodes in first year of life")+
  theme(legend.title = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.direction = "horizontal")

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/n_para_nulisa_plot_select2.png", n_para_nulisa_select, height=4, width=8, dpi=400, bg="white")



n_malaria_nulisa <- clean_data%>%
  filter(treatmentarm!="DP 2 years", timepoint=="52 weeks", mstatus==0)%>%
  filter(targetName %in% sig_mala_12$targetName)%>%
  mutate(total_n_malaria_12=ifelse(total_n_malaria_12>=2, "2+", total_n_malaria_12))%>%
  ggplot(aes(x=factor(total_n_malaria_12), y=conc, fill=factor(treatmentarm)))+
  # geom_violin(draw_quantiles = seq(0,1,0.25), color="white")+
  # ggpubr::stat_compare_means(size=2, )+
  geom_boxplot(outliers=F)+
  # # stat_summary(aes(group=c(treatmentarm)), geom="cro", fun.y=median, colour="white")+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  facet_wrap(~targetName, scales="free", nrow=2)+
  scale_fill_manual(values=treatment_palette)+
  theme_minimal()+
  xlab("number of malaria episodes in first year of life")+
  theme(legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=10),
        legend.position = "none",
        legend.direction = "horizontal")

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/n_malaria_nulisa_plot.png", n_malaria_nulisa, height=3, width=3, dpi=400, bg="white")


# msd stuff ####
mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")
maternal_treatment_arms <- haven::read_dta("~/Library/CloudStorage/Box-Box/DP+SP study/Databases and preliminary findings/Final database used for analyses/DPSP treatment allocation_FINAL.dta")

nulisa_data <- clean_data%>%
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


treatment_msd_purf <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/msd_antibodies_treatment_regression.csv")

treatment_msd_sigs <- treatment_msd_purf%>%
  filter(padj<0.15)

vaccines = c("Diptheria",     "Measles" ,      "Mumps",         "Pertussis",     "Polio",
             "Rotavirus" ,    "Rubella",       "Tetanus", "Pneumo.1.4.14")

(msd_treatment_plot <- antibodies_and_epi%>%
    filter(treatmentarm!="DP 2 years")%>%
    # filter(antigen %in% c("Diptheria", "Tetanus", "Pneumo.1.4.14", "Rotavirus", "Measles", "PIV.1"))%>%
    filter(antigen %in% treatment_msd_sigs$antigen)%>%
    filter(timepoint=="52 weeks")%>%
    ggplot(aes(x=factor(timepoint), y=titer, fill=treatmentarm))+
    geom_violin(draw_quantiles = seq(0,1,0.25), color="white")+
    # geom_boxplot(outliers = F)+
    scale_y_log10()+
    # stat_summary(aes(group=c(treatmentarm)), geom="cro", fun.y=median, colour="white")+
    # stat_summary(aes(group=treatmentarm), fun.y = median, fun.ymin = median, fun.ymax = median,
    #              geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
    facet_wrap(~antigen, scales="free", nrow=1)+
    scale_fill_manual(values=treatment_palette)+
    theme_minimal()+
    theme(legend.title = element_blank(),
          axis.title = element_blank()
          #axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
          ))+
  scale_fill_manual(values=c("darkred", "#00555A"))
  
ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/msd_treatment_plot.png", msd_treatment_plot, height=3, width=9, dpi=400, bg="white")
  
## infection history ####
any_malaria_purf <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/msd_antibodies_any_malaria_regression.csv")

any_para_purf <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/msd_antibodies_any_para_regression.csv")

any_malaria_purf_sigs <- any_malaria_purf %>%
  filter(padj<0.1)

any_para_kinda_sigs <- any_para_purf %>%
  filter(padj<0.1)



any_para_06_msd <- antibodies_and_epi%>%
  filter(treatmentarm!="DP 2 years")%>%
  filter(antigen %in% any_para_kinda_sigs$antigen)%>%
  # filter(timepoint=="52 weeks")%>%
  mutate(any_para = if_else(total_n_para_6 > 0, "some", "none"),
         any_malaria = if_else(total_n_malaria_6 > 0,"some", "none"))%>%
  ggplot(aes(x=factor(timepoint), y=titer, fill=factor(any_para)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ggpubr::stat_compare_means(size=3, label = "p.format")+
  # stat_summary(aes(group=c(treatmentarm)), geom="cro", fun.y=median, colour="white")+
  stat_summary(aes(group=any_para), fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  facet_wrap(~antigen, scales="free", nrow=1)+
  scale_fill_manual(values=c("grey", "red"))+
  theme_minimal()+
  xlab("any parasitemia in the first six months of life")+
  theme(legend.title = element_blank(),
        axis.title.y = element_blank()
        #axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
  )



ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/any_para_06_msd.png", any_para_06_msd, height=3, width=9, dpi=400, bg="white")



any_mala_12_msd <- antibodies_and_epi%>%
  filter(treatmentarm!="DP 2 years")%>%
  filter(antigen %in% any_malaria_purf_sigs$antigen)%>%
  # filter(timepoint=="52 weeks")%>%
  mutate(any_para = if_else(total_n_para_6 > 0, "some", "none"),
         any_malaria = if_else(total_n_malaria_6 > 0,"some", "none"))%>%
  ggplot(aes(x=factor(timepoint), y=titer, fill=factor(any_malaria)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ggpubr::stat_compare_means(size=3, label = "p.format")+
  # stat_summary(aes(group=c(treatmentarm)), geom="cro", fun.y=median, colour="white")+
  # stat_summary(aes(group=any_para), fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  facet_wrap(~antigen, scales="free", nrow=1)+
  scale_fill_manual(values=c("grey", "red"))+
  theme_minimal()+
  xlab("any malaria in the first six months of life")+
  theme(legend.title = element_blank(),
        axis.title.y = element_blank()
        #axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
  )



ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/any_para_06_msd.png", any_para_06_msd, height=3, width=9, dpi=400, bg="white")

# blood counts ####

mic_drop <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/MICDroP Data/MICDROP all visit database through June 30th 2025.dta")

blood_counts <- mic_drop %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7,
         "flow_age_in_integer_years"=(floor(flo_age_in_wks/52))+1)%>%
  mutate("Timepoint_in_weeks"=if_else(
    flo_age_in_wks %in% sample_ages, flo_age_in_wks, ifelse(
      flo_age_in_wks %in% sample_ages_minus, flo_age_in_wks+1, if_else(
        flo_age_in_wks %in% sample_ages_plus, flo_age_in_wks-1, 999)))
  )%>%
  mutate(disease=case_when(mstatus==2 ~ "complicated",
                           mstatus==1 ~ "uncomplicated",
                           mstatus==0 ~ "no malaria"),
         ever_comp=if_else(id %in% kids_with_comp$id, "complicated", "never complicated"))%>%
  mutate(hb_less_than_eight =ifelse(hb<8, TRUE, FALSE))%>%
  pivot_longer(cols=c(plt, wbc, hb, eosino, mono, lymph, neutro), names_to = "cell_type", values_to = "cell_freq")%>%
  filter(!is.na(cell_freq))%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(id, mic_drop_key$id)])%>%
  mutate(treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))


blood_counts_treatment <- blood_counts%>%
  filter(cell_type%in%c("mono", "plt", "hb"), Timepoint_in_weeks<53, mstatus==0)%>%
  filter(cell_type!="mono"|cell_freq<2500)%>%
  filter(cell_type!="plt"|cell_freq<10^6)%>%
  
  ggplot(., aes(x=factor(Timepoint_in_weeks), y=cell_freq,  fill = treatmentarm))+
  geom_boxplot(outliers = F)+
  stat_summary(aes(group=treatmentarm), fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  facet_wrap(~cell_type, scales="free")+
  # scale_y_log10()+
  ggpubr::stat_compare_means(label="p.signif")+
  xlab("age in weeks")+
  ylab("blood count")+
  # scale_fill_manual(values=colorspace::sequential_hcl(n = 9, palette="LaJolla"))+
  scale_fill_manual(values=treatment_palette)+
  theme_minimal()


ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/blood_counts_treatment.png", blood_counts_treatment, height=3, width=9, dpi=400, bg="white")

# extra stuff ####
# asymp specific plot ####
mstatus_purf <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/mstatus_purf.csv")


sig_symp <- mstatus_purf%>%
  filter(padj < 0.05, grepl("^symp", contrast))

sig_asymp <- mstatus_purf%>%
  filter(padj < 0.05, grepl("^asymp vs uninfected", contrast))

top_overall_genes <- mstatus_purf%>%
  ungroup()%>%
  slice_min(padj, n=20)%>%
  distinct(targetName)

top_symp_genes2 <- mstatus_purf%>%
  ungroup()%>%
  filter(grepl("^symp", contrast))%>%
  slice_min(padj, n=30)%>%
  distinct(targetName)


top_asymp_genes <- mstatus_purf%>%
  ungroup()%>%
  filter(grepl("^asymp vs uninfected", contrast))%>%
  slice_min(padj, n=20)%>%
  distinct(targetName)

top_asymp_symp_genes <- mstatus_purf%>%
  ungroup()%>%
  filter(grepl("^asymp vs symp", contrast))%>%
  filter(padj>0.5)%>%
  filter(targetName %in% top_symp_genes2$targetName)



# c("WNT16", "CSF3", "TNRSF21")
asymp_specific <- unique(sig_asymp$targetName)[unique(sig_asymp$targetName) %notin% unique(sig_malaria$targetName)]

symp_specific <- unique(sig_symp$targetName)[unique(sig_symp$targetName) %notin% unique(sig_asymp$targetName)]


clean_data %>%
  mutate(mstatusf=factor(case_when(mstatus==0&qPCRparsdens>0.5~"asymptomatic",
                                   mstatus==0&qPCRparsdens==0 ~"uninfected",
                                   mstatus==1~"symptomatic", .default = "unsure"
  )))%>%
  filter(mstatusf !="unsure")%>%
  filter(targetName %in% top_overall_genes$targetName,
         timepoint!="68 weeks")%>%
  ggplot(aes(x=factor(timepoint_num), y=conc, fill=factor(mstatusf)))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~factor(targetName, levels = top_overall_genes$targetName), scales = "free", nrow=2)+
  scale_fill_manual(values=c("orange", "darkred", "darkgrey"))+
  xlab("age in weeks")+
  theme_minimal()+
  theme(legend.position = "bottom",
        axis.title.y = element_blank())


clean_data %>%
  mutate(mstatusf=factor(case_when(mstatus==0&qPCRparsdens>0.5~"asymptomatic",
                                   mstatus==0&qPCRparsdens==0 ~"uninfected",
                                   mstatus==1~"symptomatic", .default = "unsure"
  )))%>%
  filter(mstatusf !="unsure")%>%
  filter(targetName %in% top_overall_genes$targetName,
         timepoint!="68 weeks")%>%
  ggplot(aes(x=factor(n), y=conc, fill=factor(mstatusf)))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~factor(targetName, levels = top_overall_genes$targetName), scales = "free", nrow=2)+
  scale_fill_manual(values=c("orange", "darkred", "darkgrey"))+
  xlab("age in weeks")+
  theme_minimal()+
  theme(legend.position = "bottom",
        axis.title.y = element_blank())


clean_data %>%
  mutate(mstatusf=factor(case_when(mstatus==0&qPCRparsdens>0.5~"asymptomatic",
                                   mstatus==0&qPCRparsdens==0 ~"uninfected",
                                   mstatus==1~"symptomatic", .default = "unsure"
  )))%>%
  filter(mstatusf !="unsure")%>%
  filter(targetName %in% top_overall_genes$targetName,
         timepoint!="68 weeks")%>%
  ggplot(aes(x=factor(timepoint_num), y=conc, fill=factor(mstatusf)))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~factor(targetName, levels = top_overall_genes$targetName), scales = "free", nrow=2)+
  scale_fill_manual(values=c("orange", "darkred", "darkgrey"))+
  xlab("age in weeks")+
  theme_minimal()+
  theme(legend.position = "bottom",
        axis.title.y = element_blank())


clean_data %>%
  mutate(mstatusf=factor(case_when(mstatus==0&qPCRparsdens>0.5~"asymptomatic",
                                   mstatus==0&qPCRparsdens==0 ~"uninfected",
                                   mstatus==1~"symptomatic", .default = "unsure"
  )))%>%
  filter(mstatusf !="unsure")%>%
  filter(targetName %in% top_asymp_genes$targetName,
         timepoint!="68 weeks")%>%
  ggplot(aes(x=factor(timepoint_num), y=conc, fill=factor(mstatusf)))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~factor(targetName, levels = top_asymp_genes$targetName), scales = "free", nrow=3)+
  scale_fill_manual(values=c("orange", "darkred", "darkgrey"))+
  xlab("age in weeks")+
  theme_minimal()+
  theme(legend.position = "bottom",
        axis.title.y = element_blank())

# n_infection vs age plots


mic_drop <-  haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/Jun25/MICDSpecimenBoxJun25_withclinical.dta")

mic_drop_data <- mic_drop %>%
  filter(mstatus != 0, !is.na(mstatus), incidentmalaria==1 | is.na(incidentmalaria) & mstatus==2 | incidentmalaria==0 & mstatus==2)%>%
  group_by(id) %>%
  add_count(name="total_n_infection") %>%
  mutate(n_infection = seq(1, max(total_n_infection)),
         mstatus = case_match(mstatus,
                              0~"no malaria",
                              1~"uncomplicated",
                              2~"complicated",
                              3~"quinine for AL failure",
                              4~"Q/AS failure"))%>%
  mutate(date=as.character(date))%>%
  select(id, date, n_infection)

infections <- clean_data%>%
  filter(mstatus==1)%>%
  distinct(id, date, ageinwks)


n_infections <- left_join(infections, mic_drop_data, by=c("id", "date"))


analytes <- unique(clean_data$targetName)

clean_data %>%
  left_join(., n_infections, by=c("id", "date"))%>%
  mutate(n_infection=case_when(n_infection>=4~"4+",
                               is.na(n_infection)~"no malaria",
                               .default = as.character(n_infection)))%>%
  filter(targetName %in% top_overall_genes$targetName)%>%
  mutate(targetNamef=factor(targetName, levels= top_overall_genes$targetName))%>%
  ggplot(aes(x=factor(n_infection, levels=c("no malaria", "1", "2", "3", "4+")), y=conc, fill=factor(n_infection)))+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetNamef, scales = "free", nrow=2)+
  scale_fill_manual(values=viridis::rocket(n = 5))+
  xlab("Order of infection")+
  ylab("Concentration in Blood")+
  theme_minimal(base_size = 16)+
  theme(legend.position = "none")


# mstatust age / n_infection dive deeper  / create narrative####

## growth factors ####
growth_factors_mstatus <- clean_data %>%
  # filter(timepoint_num=="52")%>%
  left_join(., n_infections, by=c("id", "date"))%>%
  mutate(n_infection=case_when(n_infection>=4~"4+",
                               is.na(n_infection)~"healthy",
                               .default = as.character(n_infection)))%>%
  filter(targetName %in% c("KDR", "ANGPT1", "EPO", "ANGPT2", "PGF"))%>%
  # mutate(targetNamef=factor(targetName, levels= top_overall_genes$targetName))%>%
  ggplot(aes(x=factor(n_infection, levels=c("healthy", "1", "2", "3", "4+")), y=conc, fill=factor(n_infection)))+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetName, scales = "free", nrow=1)+
  scale_fill_manual(values=viridis::rocket(n = 5))+
  xlab("order of malaria episode")+
  ylab("concentration in blood")+
  theme_minimal(base_size = 16)+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/growth_factors_mstatus.png", growth_factors_mstatus, height=4, width=10, dpi=400, bg="white")  

growth_factors_age <- clean_data %>%
  mutate(mstatusf=factor(case_when(mstatus==0&qPCRparsdens>0.5~"asymptomatic",
                                   mstatus==0&qPCRparsdens==0 ~"uninfected",
                                   mstatus==1~"symptomatic", .default = "unsure"
  )))%>%
  filter(mstatusf !="unsure")%>%
  filter(targetName %in% c("KDR", "ANGPT1", "EPO", "ANGPT2", "PGF"),
         timepoint!="68 weeks")%>%
  ggplot(aes(x=factor(timepoint_num), y=conc, fill=factor(mstatusf)))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetName, scales = "free", nrow=1)+
  scale_fill_manual(values=c("orange", "darkred", "darkgrey"))+
  xlab("age in weeks")+
  theme_minimal()+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/figures_for_gates_meeting/growth_factors_age.png", growth_factors_age, height=4, width=10, dpi=400, bg="white")  

## immunoregulatory proteins
(regulatory_mstatus <- clean_data %>%
    # filter(timepoint_num=="52")%>%
    left_join(., n_infections, by=c("id", "date"))%>%
    mutate(n_infection=case_when(n_infection>=4~"4+",
                                 is.na(n_infection)~"healthy",
                                 .default = as.character(n_infection)))%>%
    filter(targetName %in% c("CTLA4", "CD274", "IL1RN", "FASLG"))%>%
    # mutate(targetNamef=factor(targetName, levels= top_overall_genes$targetName))%>%
    ggplot(aes(x=factor(n_infection, levels=c("healthy", "1", "2", "3", "4+")), y=conc, fill=factor(n_infection)))+
    geom_boxplot(outliers = FALSE)+
    facet_wrap(~targetName, scales = "free", nrow=1)+
    scale_fill_manual(values=viridis::rocket(n = 5))+
    xlab("order of malaria episode")+
    ylab("concentration in blood")+
    theme_minimal(base_size = 16)+
    theme(legend.position = "none",
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=8)))


# up or across
(up_activatory_mstatus <- clean_data %>%
    # filter(timepoint_num=="52")%>%
    left_join(., n_infections, by=c("id", "date"))%>%
    mutate(n_infection=case_when(n_infection>=4~"4+",
                                 is.na(n_infection)~"healthy",
                                 .default = as.character(n_infection)))%>%
    filter(targetName %in% c("TNF", "CCL3", "CCL4"))%>%
    # mutate(targetNamef=factor(targetName, levels= top_overall_genes$targetName))%>%
    ggplot(aes(x=factor(n_infection, levels=c("healthy", "1", "2", "3", "4+")), y=conc, fill=factor(n_infection)))+
    geom_boxplot(outliers = FALSE)+
    facet_wrap(~targetName, scales = "free", nrow=1)+
    scale_fill_manual(values=viridis::rocket(n = 5))+
    xlab("order of malaria episode")+
    ylab("concentration in blood")+
    theme_minimal(base_size = 16)+
    theme(legend.position = "none",
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=8)))

# down
(down_activatory_mstatus <- clean_data %>%
    # filter(timepoint_num=="52")%>%
    left_join(., n_infections, by=c("id", "date"))%>%
    mutate(n_infection=case_when(n_infection>=4~"4+",
                                 is.na(n_infection)~"healthy",
                                 .default = as.character(n_infection)))%>%
    filter(targetName %in% c("CRP", "CX3CL1", "IFNG", "IL4"))%>%
    # mutate(targetNamef=factor(targetName, levels= top_overall_genes$targetName))%>%
    ggplot(aes(x=factor(n_infection, levels=c("healthy", "1", "2", "3", "4+")), y=conc, fill=factor(n_infection)))+
    geom_boxplot(outliers = FALSE)+
    facet_wrap(~targetName, scales = "free", nrow=1)+
    scale_fill_manual(values=viridis::rocket(n = 5))+
    xlab("order of malaria episode")+
    ylab("concentration in blood")+
    theme_minimal(base_size = 16)+
    theme(legend.position = "none",
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=8)))
# clinical immunity plot? ####


all_para_symp_summary <- all_para%>%
  filter(treatmentarm!= "DP 2 years", study=="micdrop")%>%
  group_by(bino_symp, n_any_para, treatmentarm)%>%
  summarise("n_cases"=n())%>%
  pivot_wider(names_from = bino_symp, values_from = n_cases, names_prefix = "symp_")%>%
  mutate(total_infections=symp_1+symp_0,
         risk=symp_1/total_infections)%>%
  arrange(n_any_para)

all_mala_symp_summary <- all_para%>%
  filter(treatmentarm!= "DP 2 years", study=="micdrop")%>%
  group_by(bino_symp, n_any_malaria, treatmentarm)%>%
  summarise("n_cases"=n())%>%
  pivot_wider(names_from = bino_symp, values_from = n_cases, names_prefix = "symp_")%>%
  mutate(total_infections=symp_1+symp_0,
         risk=symp_1/total_infections)%>%
  arrange(n_any_malaria)


(all_para_symp_summary_plot <- all_para_symp_summary%>%
    filter(treatmentarm=="No DP")%>%
    ggplot(., aes(x=n_any_para, y=risk))+
    theme_minimal()+
    geom_smooth(method = "lm", color="black")+
    geom_point(color="darkred")+
    geom_text(aes(y=0.8, label= paste0("frac(",symp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
    scale_x_continuous(breaks = 1:50, limits=c(1,18))+
    scale_y_continuous(labels = scales::label_percent())+
    ggtitle("")+
    xlab("")+
    ylab("P(symptoms | parsitemia)"))

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/all_para_symp_summary_plot.png", all_para_symp_summary_plot, width=5, height=5, bg="white")


(all_mala_symp_summary <- ggplot(all_mala_symp_summary, aes(x=n_any_malaria, y=risk))+
    theme_minimal()+
    geom_smooth(method = "lm", color="black")+
    geom_point(color="darkred")+
    # geom_ribbon(data=para_summary_model_prediction$se, aes(x=n_any_para, ymin = exp(lci), ymax = exp(uci)),
    #             alpha = 0.2, inherit.aes = FALSE)+
    # geom_function(fun = para_summary_model_prediction$fun, colour="black")+
    geom_text(aes(y=0.8, label= paste0("frac(",symp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
    scale_x_continuous(breaks = 1:50, limits=c(1,12))+
    scale_y_continuous(limits=c(0.3, 0.85), labels = scales::label_percent())+
    ggtitle("")+
    facet_wrap(~treatmentarm)+
    xlab("Months with any parasitemia")+
    ylab("Probability of Symptoms"))

