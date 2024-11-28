library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)



`%notin%` <- Negate(`%in%`)

da_boxplot_theme <- theme(legend.position = "none",
                          axis.title = element_blank())

# read clean data ####
clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")

clean_data <- clean_data %>%
  mutate(timepoint = factor(timepoint, levels=c("bad_baseline", "baseline", "day0", "day7", "day14", "day28")))

# linear regression ####
# the plan is to have two models for two nuances: one with timepoint, and the other with parasitemia;
# the two are heavily correlated because timepoints were chosen to be smear negative and day0 by definition have to be parasitemic.
# parasitemia has explanatory power that can sometimes rival / outrank the impact of timepoint, so a model with both may have small 
# contributions from either, meaning neither may be significant when indeed there is patent change. so both can't be in the same model.
# a model without timepoint resolves the question to what extent parasitemia has an impact on inflammation and how the magnitude of the effect
# varies between asymptomatic and symptomatic infection

## regressing on timepoint ####
fdr_cutoff = 0.05

base_zero_contrast <- t(matrix(c(0, 1, 0, 0, 0, 0, 0, 0)))
base_14_contrast <- t(matrix(c(0, 0, 1, 0, 0, 0, 0, 0)))
base_zero_infectiontype_contrast <- t(matrix(c(rep(0, 6), 1, 0)))
base_infectiontype_contrast <- t(matrix(c(0, 0, 0, 1, rep(0, 4))))
base_14_infectiontype_contrast <- t(matrix(c(rep(0, 6), 0, 1)))

combo_as_purff <- clean_data %>%
  filter(infectiontype%in%c("A", "S"), timepoint %notin% c("day7", "day28"), timepoint!="bad_baseline")%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  group_by(targetName)%>%
  nest() %>%
  #mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_zero=map(model, ~multcomp::glht(., base_zero_contrast)),
         base_zero_p=map_dbl(base_zero, ~summary(.)$test$pvalues)) %>%
  mutate(base_14=map(model, ~multcomp::glht(., base_14_contrast)),
         base_14_p=map_dbl(base_14, ~summary(.)$test$pvalues)) %>%
  mutate(base_14_infectiontype=map(model, ~multcomp::glht(., base_14_infectiontype_contrast)),
         base_14_infectiontype_p=map_dbl(base_14_infectiontype, ~summary(.)$test$pvalues)) %>%
  mutate(base_zero_infectiontype=map(model, ~multcomp::glht(., base_zero_infectiontype_contrast)),
         base_zero_infectiontype_p=map_dbl(base_zero_infectiontype, ~summary(.)$test$pvalues)) %>%
  mutate(base_infectiontype=map(model, ~multcomp::glht(., base_infectiontype_contrast)),
         base_infectiontype_p=map_dbl(base_infectiontype, ~summary(.)$test$pvalues))%>%
  ungroup()%>%
  mutate(base_zero_padj=p.adjust(base_zero_p, method="BH"),
         base_14_padj=p.adjust(base_14_p, method="BH"),
         base_14_infectiontype_padj=p.adjust(base_14_infectiontype_p, method="BH"),
         base_infectiontype_padj=p.adjust(base_infectiontype_p, method="BH"),
         base_zero_infectiontype_padj=p.adjust(base_zero_infectiontype_p, method="BH")
  )


results_table <- combo_as_purff %>%
  dplyr::select(targetName,
                base_zero_padj,
                base_14_padj,
                base_14_infectiontype_padj,
                base_infectiontype_padj,
                base_zero_infectiontype_padj)%>%
  ungroup()

# table(results_table$base_zero_padj<fdr_cutoff)# 63 
# table(results_table$base_14_padj<fdr_cutoff)# 0 this is different when you model A on their own
# table(results_table$base_infectiontype_padj<fdr_cutoff)# 15 
# table(results_table$base_zero_infectiontype_padj<fdr_cutoff)# 89 --> ~110 vary by timepoint 
# table(results_table$base_para_padj<fdr_cutoff)# 94 --> ~110 vary by timepoint 

#104
sig_base_zero_infectiontype <- results_table%>%
  filter(base_zero_infectiontype_padj<fdr_cutoff)%>%
  arrange(base_zero_infectiontype_padj)

write.csv(sig_base_zero_infectiontype, "~/postdoc/stanford/plasma_analytes/MUSICAL/combo/differential_abundance/sig_base_zero_infectiontype.csv", row.names = F)

#39
sig_base_zero <- results_table%>%
  filter(base_zero_padj<fdr_cutoff)%>%
  arrange(base_zero_padj)

#42
sig_base_14 <- results_table%>%
  filter(base_14_padj<fdr_cutoff)%>%
  arrange(base_14_padj)


for(i in 1:ceiling(nrow(sig_base_zero)/16)){
  
  plt <- clean_data %>%
    filter(infectiontype %in% c("A"), timepoint!="day28", timepoint!="bad_baseline")%>%
    filter(targetName %in% sig_base_zero$targetName[seq((i-1)*16+1, i*16)])%>%
    mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
    ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
    # geom_point()+#
    # geom_line(aes(group=id))+
    geom_boxplot(outliers = FALSE)+
    # geom_violin(draw_quantiles = seq(0,1,0.25))+
    # ggtitle("regulated during asymptomatic parasitemia")+
    facet_wrap(~targetName, scales = "free")+
    scale_fill_manual(values=viridis::magma(5))+
    theme_minimal()+
    da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/combo_as_base_zero", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}



for(i in 1:ceiling(nrow(sig_base_14)/16)){
  
  plt <- clean_data %>%
    filter(infectiontype %in% c("A", "S"), timepoint=="baseline", timepoint!="bad_baseline")%>%
    filter(targetName %in% sig_base_zero_14$targetName[seq((i-1)*16+1, i*16)])%>%
    mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
    ggplot(aes(x=factor(timepoint), y=concentration, fill=infectiontype))+
    # geom_point()+#
    # geom_line(aes(group=id))+
    geom_boxplot(outliers = FALSE)+
    # geom_violin(draw_quantiles = seq(0,1,0.25))+
    # ggtitle("regulated during asymptomatic parasitemia")+
    facet_wrap(~targetName, scales = "free")+
    scale_fill_manual(values=viridis::magma(5))+
    theme_minimal()+
    da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/combo_as_base_14", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}


for(i in 1:ceiling(nrow(sig_base_zero_infectiontype)/16)){
  
  plt <- clean_data %>%
    filter(infectiontype %in% c("S"), timepoint!="day28", timepoint!="bad_baseline")%>%
    filter(targetName %in% sig_base_zero_infectiontype$targetName[seq((i-1)*16+1, i*16)])%>%
    mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
    ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
    # geom_point()+#
    # geom_line(aes(group=id))+
    geom_boxplot(outliers = FALSE)+
    # geom_violin(draw_quantiles = seq(0,1,0.25))+
    # ggtitle("regulated during asymptomatic parasitemia")+
    facet_wrap(~targetName, scales = "free")+
    scale_fill_manual(values=viridis::magma(5))+
    theme_minimal()+
    da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/combo_as_base_zero_infectiontype", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}


## regressing on parasitemia, not timepoint ####
qpcr_infectiontype_contrast <- t(matrix(c(0, 0, 0, 0, 0, 1)))
qpcr_contrast <- t(matrix(c(0, 1, 0, 0, 0, 0)))
infectiontype_contrast <- t(matrix(c(0, 0, 1, 0, 0, 0)))


combo_as_purff_para <- clean_data %>%
  filter(infectiontype%in%c("A", "S"), timepoint %notin% c("day7", "day28"), timepoint!="bad_baseline")%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")),
         log_qpcr = log10(qpcr+0.1))%>%
  group_by(targetName)%>%
  nest() %>%
  # mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~log10(parasitedensity+0.1)*infectiontype+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(qpcr_infectiontype=map(model, ~multcomp::glht(., qpcr_infectiontype_contrast)),
         qpcr_infectiontype_p=map_dbl(qpcr_infectiontype, ~summary(.)$test$pvalues)) %>%
  mutate(qpcr=map(model, ~multcomp::glht(., qpcr_contrast)),
         qpcr_p=map_dbl(qpcr, ~summary(.)$test$pvalues)) %>%
  mutate(infectiontype=map(model, ~multcomp::glht(., infectiontype_contrast)),
         infectiontype_p=map_dbl(infectiontype, ~summary(.)$test$pvalues)) %>%
  mutate(qpcr_infectiontype_coef=map_dbl(summary, ~coef(.)[6]),
         qpcr_coef=map_dbl(summary, ~coef(.)[2]),
         infectiontype_coef=map_dbl(summary, ~coef(.)[3]))%>%
  ungroup()%>%
  # group_by(timepoint)%>%
  mutate(qpcr_infectiontype_padj=p.adjust(qpcr_infectiontype_p, method="BH"),
         qpcr_padj=p.adjust(qpcr_p, method="BH"),
         infectiontype_padj=p.adjust(infectiontype_p, method="BH")
  )


para_results_table <- combo_as_purff_para %>%
  dplyr::select(targetName,
                qpcr_infectiontype_coef,
                qpcr_infectiontype_padj,
                qpcr_coef,
                qpcr_padj, 
                infectiontype_coef,
                infectiontype_padj,
  )%>%
  ungroup()

# table(para_results_table$qpcr_infectiontype_padj<fdr_cutoff)# 63 
# table(para_results_table$qpcr_padj<fdr_cutoff)# 0 this is different when you model A on their own
# table(para_results_table$base_infectiontype_padj<fdr_cutoff)# 15 
# table(para_results_table$qpcr_infectiontype_infectiontype_padj<fdr_cutoff)# 89 --> ~110 vary by timepoint 
# table(para_results_table$base_para_padj<fdr_cutoff)# 94 --> ~110 vary by timepoint 

#62 with qPCR; 79 with BS
sig_qpcr_infectiontype <- para_results_table%>%
  filter(qpcr_infectiontype_padj<fdr_cutoff)%>%
  arrange(qpcr_infectiontype_padj)

#57 with qPCR; 63 with BS
sig_qpcr <- para_results_table%>%
  filter(qpcr_padj<fdr_cutoff)%>%
  arrange(qpcr_padj)
#18 with qPCR

sig_infectiontype <- para_results_table%>%
  filter(infectiontype_padj<fdr_cutoff)%>%
  arrange(infectiontype_padj)

sig_qpcr$targetName[sig_qpcr$targetName %notin% sig_infectiontype$targetName]

### visualise correlation between parasitemia and inflammation####

extreme_coefs <- sig_qpcr_infectiontype %>%
  slice_max(n=10, order_by = abs(.$qpcr_infectiontype_coef))

clean_data %>%
  filter(infectiontype%in%c("A", "S"), timepoint %notin% c("day7", "day28"))%>%
  filter(targetName %in% extreme_coefs$targetName[1:10])%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=parasitedensity+0.1, y=concentration, color=infectiontype))+
  geom_point()+#
  geom_smooth(method="lm")+
  ggpubr::stat_cor()+
  scale_x_log10()+
  facet_wrap(~targetName, scales = "free")+
  scale_color_manual(values=viridis::magma(3))+
  theme_minimal()



clean_data %>%
  filter(infectiontype%in%c("A", "S"), timepoint %notin% c("day7", "day28"), timepoint!="bad_baseline")%>%
  filter(targetName %in% sig_qpcr_infectiontype$targetName)%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=qpcr+0.1, y=concentration, color=infectiontype))+
  geom_smooth(method = "lm")+
  scale_x_log10()+
  facet_wrap(~targetName, scales = "fixed")+
  scale_color_manual(values=viridis::magma(3))+
  theme_minimal()




clean_data %>%
  filter(infectiontype %in% c("A", "S"), timepoint=="baseline")%>%
  filter(targetName %in% sig_infectiontype$targetName)%>%
  filter(log_qpcr<1)%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=infectiontype, y=concentration, fill=infectiontype))+
  # geom_point()+#
  # geom_line(aes(group=id))+
  geom_boxplot(outliers = FALSE)+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  # ggtitle("regulated during asymptomatic parasitemia")+
  facet_wrap(~targetName, scales = "free")+
  scale_fill_manual(values=viridis::magma(5))+
  theme_minimal()+
  da_boxplot_theme





for(i in 1:ceiling(nrow(sig_qpcr_infectiontype)/16)){
  
  plotable_proteins <- factor(sig_qpcr_infectiontype$targetName[seq((i-1)*16+1, i*16)], levels = sig_qpcr_infectiontype$targetName)
  
  plt <- clean_data %>%
    filter(infectiontype %in% c("A"), timepoint!="day28", timepoint!="bad_baseline")%>%
    filter(targetName %in% plotable_proteins)%>%
    mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
    ggplot(aes(x=infectiontype, y=concentration, fill=timepoint))+
    # geom_point()+#
    # geom_line(aes(group=id))+
    geom_boxplot(outliers = FALSE)+
    # geom_violin(draw_quantiles = seq(0,1,0.25))+
    # ggtitle("regulated during asymptomatic parasitemia")+
    facet_wrap(~factor(targetName, levels = plotable_proteins), scales = "free")+
    scale_fill_manual(values=viridis::magma(5))+
    theme_minimal()+
    da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/combo_as_qpcr_infectiontype", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}



for(i in 1:ceiling(nrow(sig_qpcr)/16)){
  
  plotable_proteins <- factor(sig_qpcr$targetName[seq((i-1)*16+1, i*16)], levels = sig_qpcr$targetName)
  
  plt <- clean_data %>%
    filter(infectiontype %in% c("S"), timepoint!="day28", timepoint!="bad_baseline")%>%
    filter(targetName %in% plotable_proteins)%>%
    mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
    ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
    # geom_point()+#
    # geom_line(aes(group=id))+
    geom_boxplot(outliers = FALSE)+
    # geom_violin(draw_quantiles = seq(0,1,0.25))+
    # ggtitle("regulated during asymptomatic parasitemia")+
    facet_wrap(~factor(targetName, levels=plotable_proteins), scales = "free")+
    scale_fill_manual(values=viridis::magma(5))+
    theme_minimal()+
    da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/combo_as_qpcr", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}



for(i in 1:ceiling(nrow(sig_infectiontype)/16)){
  
  plotable_proteins <- factor(sig_infectiontype$targetName[seq((i-1)*16+1, i*16)])
  
  plt <- clean_data %>%
    filter(infectiontype %in% c("A", "S"), timepoint!="day28", timepoint=="baseline")%>%
    filter(targetName %in% plotable_proteins)%>%
    mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
    ggplot(aes(x=infectiontype, y=concentration, fill=infectiontype))+
    # geom_point()+#
    # geom_line(aes(group=id))+
    geom_boxplot(outliers = FALSE)+
    # geom_violin(draw_quantiles = seq(0,1,0.25))+
    # ggtitle("regulated during asymptomatic parasitemia")+
    facet_wrap(~factor(targetName, levels=plotable_proteins), scales = "free")+
    scale_fill_manual(values=viridis::magma(3))+
    theme_minimal()+
    da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/combo_as_infectiontype", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}






## baselines only ####
base_base_contrast <- t(matrix(c(0,1)))


base_only_purff <- clean_data %>%
  filter(infectiontype %in% c("S", "A"), timepoint=="baseline")%>%
  group_by(targetName)%>%
  mutate(bino_infectiontype=if_else(infectiontype=="S", 1, 0))%>%
  nest() %>%
  # mutate(model=map(data, ~lm(concentration~infectiontype+id, data=.))) %>%
  # mutate(model=map(data, ~lme4::lmer(concentration~infectiontype+(1|id), data=.))) %>%
  # this one "works"
  mutate(model=map(data, ~glm(bino_infectiontype~concentration+id, data=., family="binomial"))) %>%
  # mutate(model=map(data, ~lme4::glmer(bino_infectiontype~concentration+(1|id), data=., family="binomial")))%>%
  # mutate(model=map(data, ~lme4::lmer(concentration~infectiontype+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  # mutate(base_base=map(model, ~multcomp::glht(., base_base_contrast)),
  #        base_base_p=map_dbl(summary, ~coef(.)[8])) %>%
  # this one "works" too
  mutate(base_base_p=map_dbl(summary, ~coef(.)[140])) %>%
  # mutate(base_14=map(model, ~multcomp::glht(., base_14_contrast)),
  #        base_14_p=map_dbl(base_14, ~summary(.)$test$pvalues)) %>%
  # mutate(zero_14=map(model, ~multcomp::glht(., zero_14_contrast)),
  #        zero_14_p=map_dbl(zero_14, ~summary(.)$test$pvalues))%>%
  # mutate(zero_28=map(model, ~multcomp::glht(., zero_28_contrast)),
  #        zero_28_p=map_dbl(zero_28, ~summary(.)$test$pvalues))%>%
  ungroup()%>%
  mutate(base_base_padj=p.adjust(base_base_p, method="BH"),
         # base_14_padj=p.adjust(base_14_p, method="BH"),
         # zero_14_padj=p.adjust(zero_14_p, method="BH"),
         # zero_28_padj=p.adjust(zero_28_p, method="BH")
  )

base_only_results_table <- base_only_purff %>%
  dplyr::select(targetName,
                base_base_p,
                base_base_padj)%>%
  ungroup()

#8 with glm
base_only_sig <- base_only_results_table %>%
  filter(base_base_padj<fdr_cutoff)%>%
  select(targetName)


base_only_plot <- clean_data %>%
  filter(timepoint=="baseline", infectiontype %in% c("S", "A"))%>%
  # filter(targetName %in% c("CXCL8", "MMP8", "MMP9", "MP9", "OSM", "PTX3", "TNFSF11"))%>%
  filter(targetName %in% base_only_sig$targetName)%>%
  ggplot(aes(x=infectiontype, y=concentration, fill=targetName))+
  # geom_violin(draw_quantiles = c(seq(0,1,0.25)), colour = "darkgrey")+
  geom_boxplot(outliers = FALSE)+
  ggtitle("baselines before parasitemia")+
  facet_wrap(~targetName, scales = "free", ncol=6)+
  scale_fill_manual(values=c(viridis::magma(18)))+
  theme_minimal()+
  theme(legend.position = "none")+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/base_only_bino_glm.png", base_only_plot, width=8, height=6, bg="white", dpi=444)


neutrophil_base_only_plot <- clean_data %>%
  filter(timepoint=="baseline", infectiontype %in% c("S", "A"))%>%
  filter(targetName %in% c("CXCL8", "MMP8", "MMP9", "MP9", "OSM", "PTX3", "TNFSF11"))%>%
  # filter(targetName %in% base_only_sig$targetName)%>%
  ggplot(aes(x=infectiontype, y=concentration, fill=infectiontype))+
  # geom_violin(draw_quantiles = c(seq(0,1,0.25)), colour = "darkgrey")+
  geom_boxplot(outliers = FALSE)+
  ggtitle("baselines before parasitemia")+
  facet_wrap(~targetName, scales = "free", ncol=6)+
  scale_fill_manual(values=c(viridis::magma(5)))+
  theme_minimal()+
  theme(legend.position = "none")+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/neutrophil_base_only.png", neutrophil_base_only_plot, width=8, height=6, bg="white", dpi=444)



## nmf ####
base_zero_contrast <- t(matrix(c(0,1,0,0,0,0)))
base_7_contrast <- t(matrix(c(0,0,1,0,0,0)))
base_14_contrast <- t(matrix(c(0,0,0,1,0,0)))
zero_14_contrast <- t(matrix(c(0,-1,0,1,0,0)))


nmf_only_purff <- clean_data %>%
  filter(infectiontype=="NM", timepoint!="day28")%>%
  group_by(targetName)%>%
  nest() %>%
  mutate(model=map(data, ~lme4::lmer(concentration~timepoint+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_zero=map(model, ~multcomp::glht(., base_zero_contrast)),
         base_zero_p=map_dbl(base_zero, ~summary(.)$test$pvalues)) %>%
  mutate(base_14=map(model, ~multcomp::glht(., base_14_contrast)),
         base_14_p=map_dbl(base_14, ~summary(.)$test$pvalues)) %>%
  mutate(base_7=map(model, ~multcomp::glht(., base_7_contrast)),
         base_7_p=map_dbl(base_7, ~summary(.)$test$pvalues))%>%
  mutate(zero_14=map(model, ~multcomp::glht(., zero_14_contrast)),
         zero_14_p=map_dbl(zero_14, ~summary(.)$test$pvalues))%>%
  ungroup()%>%
  mutate(base_zero_padj=p.adjust(base_zero_p, method="BH"),
         base_7_padj=p.adjust(base_7_p, method="BH"),
         base_14_padj=p.adjust(base_14_p, method="BH"),
         zero_14_padj=p.adjust(zero_14_p, method="BH")
  )


nmf_results_table <- nmf_only_purff %>%
  dplyr::select(targetName,
                base_zero_padj,
                base_7_padj,
                base_14_padj,
                zero_14_padj
                )%>%
  ungroup()

nmf_sig_base_zero <- nmf_results_table %>%
  filter(base_zero_padj<fdr_cutoff)%>%
  arrange(desc(base_zero_padj))%>%
  select(targetName)



nmf_plot <- clean_data %>%
  filter(targetName %in% c("CRP", "IL6", "IL1RL1", "CCL2", "TNF", "IFNG"), infectiontype %in% c("NM"), timepoint!="day28")%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=interaction(infectiontype, timepoint)))+
  geom_boxplot()+
  geom_point(position = position_dodge(width=0.75))+# geom_boxplot(outliers = FALSE)+
  ggtitle("non-malarial fever")+
  facet_wrap(~targetName, scales = "free")+
  scale_fill_manual(values=c(viridis::magma(4)))+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/nmf_plot.png", nmf_plot, width=4, height=4, bg="white", dpi=444)




# post hoc clean up ####
## regrssing on timepoint ####

posthoc_clean_data <- clean_data %>%
  filter(parasitedensity<500 & infectiontype=="S" & timepoint=="day0" | parasitedensity>500 & timepoint=="baseline")%>%
  anti_join(clean_data, ., by="sample_id")
  
fdr_cutoff = 0.05

base_zero_contrast <- t(matrix(c(0, 1, 0, 0, 0, 0, 0, 0)))
base_14_contrast <- t(matrix(c(0, 0, 1, 0, 0, 0, 0, 0)))
base_zero_infectiontype_contrast <- t(matrix(c(rep(0, 6), 1, 0)))
base_infectiontype_contrast <- t(matrix(c(0, 0, 0, 1, rep(0, 4))))
base_14_infectiontype_contrast <- t(matrix(c(rep(0, 6), 0, 1)))

posthoc_combo_as_purff <- posthoc_clean_data %>%
  filter(infectiontype%in%c("A", "S"), timepoint %notin% c("day7", "day28"), timepoint!="bad_baseline")%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  group_by(targetName)%>%
  nest() %>%
  #mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_zero=map(model, ~multcomp::glht(., base_zero_contrast)),
         base_zero_p=map_dbl(base_zero, ~summary(.)$test$pvalues)) %>%
  mutate(base_14=map(model, ~multcomp::glht(., base_14_contrast)),
         base_14_p=map_dbl(base_14, ~summary(.)$test$pvalues)) %>%
  mutate(base_14_infectiontype=map(model, ~multcomp::glht(., base_14_infectiontype_contrast)),
         base_14_infectiontype_p=map_dbl(base_14_infectiontype, ~summary(.)$test$pvalues)) %>%
  mutate(base_zero_infectiontype=map(model, ~multcomp::glht(., base_zero_infectiontype_contrast)),
         base_zero_infectiontype_p=map_dbl(base_zero_infectiontype, ~summary(.)$test$pvalues)) %>%
  mutate(base_infectiontype=map(model, ~multcomp::glht(., base_infectiontype_contrast)),
         base_infectiontype_p=map_dbl(base_infectiontype, ~summary(.)$test$pvalues))%>%
  ungroup()%>%
  mutate(base_zero_padj=p.adjust(base_zero_p, method="BH"),
         base_14_padj=p.adjust(base_14_p, method="BH"),
         base_14_infectiontype_padj=p.adjust(base_14_infectiontype_p, method="BH"),
         base_infectiontype_padj=p.adjust(base_infectiontype_p, method="BH"),
         base_zero_infectiontype_padj=p.adjust(base_zero_infectiontype_p, method="BH")
  )


posthoc_results_table <- posthoc_combo_as_purff %>%
  dplyr::select(targetName,
                base_zero_padj,
                base_14_padj,
                base_14_infectiontype_padj,
                base_infectiontype_padj,
                base_zero_infectiontype_padj)%>%
  ungroup()

# table(results_table$base_zero_padj<fdr_cutoff)# 63 
# table(results_table$base_14_padj<fdr_cutoff)# 0 this is different when you model A on their own
# table(results_table$base_infectiontype_padj<fdr_cutoff)# 15 
# table(results_table$base_zero_infectiontype_padj<fdr_cutoff)# 89 --> ~110 vary by timepoint 
# table(results_table$base_para_padj<fdr_cutoff)# 94 --> ~110 vary by timepoint 

#104
sig_base_zero_infectiontype <- results_table%>%
  filter(base_zero_infectiontype_padj<fdr_cutoff)%>%
  arrange(base_zero_infectiontype_padj)

#39
sig_base_zero <- results_table%>%
  filter(base_zero_padj<fdr_cutoff)%>%
  arrange(base_zero_padj)

#42
sig_base_14 <- results_table%>%
  filter(base_14_padj<fdr_cutoff)%>%
  arrange(base_14_padj)



## regressing on parasitemia ####

ble(para_results_table$base_para_padj<fdr_cutoff)# 94 --> ~110 vary by timepoint 

#62 with qPCR; 79 with BS
sig_qpcr_infectiontype <- para_results_table%>%
  filter(qpcr_infectiontype_padj<fdr_cutoff)%>%
  arrange(qpcr_infectiontype_padj)

#57 with qPCR; 63 with BS
sig_qpcr <- para_results_table%>%
  filter(qpcr_padj<fdr_cutoff)%>%
  arrange(qpcr_padj)
#18 with qPCR

sig_infectiontype <- para_results_table%>%
  filter(infectiontype_padj<fdr_cutoff)%>%
  arrange(infectiontype_padj)

qpcr_infectiontype_contrast <- t(matrix(c(0, 0, 0, 0, 0, 1)))
qpcr_contrast <- t(matrix(c(0, 1, 0, 0, 0, 0)))
infectiontype_contrast <- t(matrix(c(0, 0, 1, 0, 0, 0)))


posthoc_clean_data_as_purff_para <- posthoc_clean_data %>%
  filter(infectiontype%in%c("A", "S"), timepoint %notin% c("day7", "day28"), timepoint!="bad_baseline")%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")),
         log_qpcr = log10(qpcr+0.1))%>%
  group_by(targetName)%>%
  nest() %>%
  # mutate(model=map(data, ~lme4::lmer(concentration~timepoint*infectiontype+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~log10(parasitedensity+0.1)*infectiontype+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(qpcr_infectiontype=map(model, ~multcomp::glht(., qpcr_infectiontype_contrast)),
         qpcr_infectiontype_p=map_dbl(qpcr_infectiontype, ~summary(.)$test$pvalues)) %>%
  mutate(qpcr=map(model, ~multcomp::glht(., qpcr_contrast)),
         qpcr_p=map_dbl(qpcr, ~summary(.)$test$pvalues)) %>%
  mutate(infectiontype=map(model, ~multcomp::glht(., infectiontype_contrast)),
         infectiontype_p=map_dbl(infectiontype, ~summary(.)$test$pvalues)) %>%
  mutate(qpcr_infectiontype_coef=map_dbl(summary, ~coef(.)[6]),
         qpcr_coef=map_dbl(summary, ~coef(.)[2]),
         infectiontype_coef=map_dbl(summary, ~coef(.)[3]))%>%
  ungroup()%>%
  # group_by(timepoint)%>%
  mutate(qpcr_infectiontype_padj=p.adjust(qpcr_infectiontype_p, method="BH"),
         qpcr_padj=p.adjust(qpcr_p, method="BH"),
         infectiontype_padj=p.adjust(infectiontype_p, method="BH")
  )


posthoc_clean_data_para <- combo_as_purff_para %>%
  dplyr::select(targetName,
                qpcr_infectiontype_coef,
                qpcr_infectiontype_padj,
                qpcr_coef,
                qpcr_padj, 
                infectiontype_coef,
                infectiontype_padj,
  )%>%
  ungroup()

# table(para_results_table$qpcr_infectiontype_padj<fdr_cutoff)# 63 
# table(para_results_table$qpcr_padj<fdr_cutoff)# 0 this is different when you model A on their own
# table(para_results_table$base_infectiontype_padj<fdr_cutoff)# 15 
# table(para_results_table$qpcr_infectiontype_infectiontype_padj<fdr_cutoff)# 89 --> ~110 vary by timepoint 
# table(para_results_table$base_para_padj<fdr_cutoff)# 94 --> ~110 vary by timepoint 

#62 with qPCR; 79 with BS
posthoc_sig_qpcr_infectiontype <- para_results_table%>%
  filter(qpcr_infectiontype_padj<fdr_cutoff)%>%
  arrange(qpcr_infectiontype_padj)

#57 with qPCR; 63 with BS
posthoc_sig_qpcr <- para_results_table%>%
  filter(qpcr_padj<fdr_cutoff)%>%
  arrange(qpcr_padj)

#18 with qPCR
posthoc_sig_infectiontype <- para_results_table%>%
  filter(infectiontype_padj<fdr_cutoff)%>%
  arrange(infectiontype_padj)
