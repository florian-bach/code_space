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
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))

# linear regression ####

## combined modelling of both A & S infections ####

fdr_cutoff = 0.05

base_zero_contrast <- t(matrix(c(0, 1, 0, 0, 0, 0, 0, 0, 0)))
base_zero_class_contrast <- t(matrix(c(rep(0, 7), 1, 0)))
base_class_contrast <- t(matrix(c(0, 0, 0, 1, rep(0, 5))))
base_14_contrast <- t(matrix(c(rep(0, 7), 0, 1)))
base_para_contrast <- t(matrix(c(rep(0, 4), 1, rep(0, 4))))

#short contrasts
# base_zero_contrast <- t(matrix(c(0, 1, rep(0, 6))))
# base_zero_class_contrast <- t(matrix(c(rep(0, 7), 1)))
# base_class_contrast <- t(matrix(c(0, 0, 0, 1, rep(0, 4))))
# base_14_contrast <- t(matrix(c(rep(0, 6), 0, 1)))
# base_para_contrast <- t(matrix(c(rep(0, 4), 1, rep(0, 3))))

combo_as_purff <- clean_data %>%
  filter(class%in%c("A", "S"), timepoint %notin% c("day7", "day28"), timepoint!="bad_baseline")%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")),
         log_qpcr = log10(qpcr+0.1))%>%
  group_by(targetName)%>%
  nest() %>%
  # mutate(model=map(data, ~lme4::lmer(concentration~timepoint*class+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~log_qpcr*class+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_zero=map(model, ~multcomp::glht(., base_zero_contrast)),
         base_zero_p=map_dbl(base_zero, ~summary(.)$test$pvalues)) %>%
  mutate(base_14=map(model, ~multcomp::glht(., base_14_contrast)),
         base_14_p=map_dbl(base_14, ~summary(.)$test$pvalues)) %>%
  mutate(base_zero_class=map(model, ~multcomp::glht(., base_zero_class_contrast)),
         base_zero_class_p=map_dbl(base_zero_class, ~summary(.)$test$pvalues)) %>%
  mutate(base_class=map(model, ~multcomp::glht(., base_class_contrast)),
         base_class_p=map_dbl(base_class, ~summary(.)$test$pvalues))%>%
  mutate(base_para=map(model, ~multcomp::glht(., base_para_contrast)),
         base_para_p=map_dbl(base_para, ~summary(.)$test$pvalues))%>%
  ungroup()%>%
  mutate(base_zero_padj=p.adjust(base_zero_p, method="BH"),
         base_14_padj=p.adjust(base_14_p, method="BH"),
         base_class_padj=p.adjust(base_class_p, method="BH"),
         base_para_padj=p.adjust(base_para_p, method="BH"),
         base_zero_class_padj=p.adjust(base_zero_class_p, method="BH")
  )


results_table <- combo_as_purff %>%
  dplyr::select(targetName,
                base_zero_padj,
                base_14_padj,
                base_class_padj,
                base_para_padj,
                base_zero_class_padj)%>%
  ungroup()

# table(results_table$base_zero_padj<fdr_cutoff)# 63 
# table(results_table$base_14_padj<fdr_cutoff)# 0 this is different when you model A on their own
# table(results_table$base_class_padj<fdr_cutoff)# 15 
# table(results_table$base_zero_class_padj<fdr_cutoff)# 89 --> ~110 vary by timepoint 
# table(results_table$base_para_padj<fdr_cutoff)# 94 --> ~110 vary by timepoint 

sig_base_zero_class <- results_table%>%
  filter(base_zero_class_padj<fdr_cutoff)%>%
  arrange(base_zero_class_padj)

sig_base_zero <- results_table%>%
  filter(base_zero_padj<fdr_cutoff)%>%
  arrange(base_zero_padj)

sig_base_class <- results_table%>%
  filter(base_class_padj<fdr_cutoff)%>%
  arrange(base_class_padj)

# model with age and sex
# table(results_table$base_zero_padj<0.05)# 89 with cont. pcr
# table(results_table$base_class_padj<0.05)# 0 with cont. pcr
# table(results_table$base_zero_para_padj<0.05)# 0 with cont. pcr
# 
# table(results_table$base_zero_padj<0.05)# 97 with day0 pcr
# table(results_table$base_class_padj<0.05)# 15 with day0 pcr
# table(results_table$base_14_padj<0.05)# 0 with cont. pcr
# 
# table(results_table$base_zero_padj<0.05)# 95 without pcr
# table(results_table$base_class_padj<0.05)# 0 without pcr
# table(results_table$base_14_padj<0.05)# 0 without pcr
# 
# # model without age and sex
# table(results_table$base_zero_padj<0.05)# 95 without pcr
# table(results_table$base_class_padj<0.05)# 0 without pcr
# table(results_table$base_14_padj<0.05)# 0 without pcr
# 
# table(results_table$base_zero_padj<0.05)# 98 with day0 pcr
# table(results_table$base_class_padj<0.05)# 15 with day0 pcr
# table(results_table$base_14_padj<0.05)# 0 with day0 pcr

for(i in 1:ceiling(nrow(sig_base_zero)/16)){
  
  plt <- clean_data %>%
    filter(class %in% c("A"), timepoint!="day28", timepoint!="bad_baseline")%>%
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



for(i in 1:ceiling(nrow(sig_base_zero_class)/16)){
  
  plt <- clean_data %>%
    filter(class %in% c("S"), timepoint!="day28", timepoint!="bad_baseline")%>%
    filter(targetName %in% sig_base_zero_class$targetName[seq((i-1)*16+1, i*16)])%>%
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
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/combo_as_base_zero_class", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}



  for(i in 1:ceiling(nrow(sig_base_class)/16)){
  
  plt <- clean_data %>%
    filter(class %in% c("A", "S"), timepoint=="baseline", timepoint!="bad_baseline")%>%
    filter(targetName %in% sig_base_zero_class$targetName[seq((i-1)*16+1, i*16)])%>%
    mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
    ggplot(aes(x=factor(timepoint), y=concentration, fill=class))+
    # geom_point()+#
    # geom_line(aes(group=id))+
    geom_boxplot(outliers = FALSE)+
    # geom_violin(draw_quantiles = seq(0,1,0.25))+
    # ggtitle("regulated during asymptomatic parasitemia")+
    facet_wrap(~targetName, scales = "free")+
    scale_fill_manual(values=viridis::magma(5))+
    theme_minimal()+
    da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/combo_as_base_class", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}

clean_data %>%
  filter(class%in%c("A", "S"), timepoint %notin% c("day7", "day28"), timepoint!="bad_baseline")%>%
  filter(targetName %in% results_table$targetName[results_table$base_class_padj<0.1])%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  filter(timepoint=="baseline")%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=class))+
  # geom_point()+#
  # geom_line(aes(group=id))+
  geom_boxplot(outliers = FALSE)+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  ggtitle("differentially abundant at baseline")+
  facet_wrap(~targetName, scales = "free")+
  scale_fill_manual(values=viridis::magma(4))+
  theme_minimal()+
  da_boxplot_theme+
  theme(legend.position = "right")




clean_data %>%
  filter(class%in%c("A", "S"), timepoint %notin% c("day7", "day28"), timepoint!="bad_baseline")%>%
  filter(targetName %in% c("CXCL8", "MMP8", "MMP9", "MPO", "PTX3", "OSM", "TNFSF11"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  filter(timepoint=="baseline")%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=class))+
  # geom_point()+#
  # geom_line(aes(group=id))+
  geom_boxplot(outliers = FALSE)+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  ggtitle("differentially abundant at baseline in pilot")+
  facet_wrap(~targetName, scales = "free")+
  scale_fill_manual(values=viridis::magma(4))+
  theme_minimal()+
  da_boxplot_theme+
  theme(legend.position = "right")


asymp_pcr_plot <- clean_data %>%
  filter(class%in%c("A"), timepoint %notin% c("day28"), timepoint!="bad_baseline")%>%
  distinct(id, timepoint, qpcr)%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=qpcr+0.1, color=id))+
  geom_point()+
  geom_line(aes(group=id))+
  ggtitle("asymptomatic parasitemia")+
  scale_color_manual(values=viridis::magma(60))+
  scale_y_log10(limits=c(0.1, 10^6))+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/asymp_pcr_plot.png", asymp_pcr_plot, width=5, height=5, dpi=444, bg="white")


symp_pcr_plot <- clean_data %>%
  filter(class%in%c("S"), timepoint %notin% c("day28"), timepoint!="bad_baseline")%>%
  distinct(id, timepoint, qpcr)%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=qpcr+0.1, color=id))+
  geom_point()+
  geom_line(aes(group=id))+
  ggtitle("symptomatic parasitemia")+
  scale_color_manual(values=viridis::magma(60))+
  scale_y_log10(limits=c(0.1, 10^6))+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/symp_pcr_plot.png", symp_pcr_plot, width=5, height=5, dpi=444, bg="white")

## symptomatic only ####
fdr_cutoff = 0.05

base_zero_contrast <- t(matrix(c(0,1,0, rep(0, 2), 0)))
base_7_contrast <- t(matrix(c(0,0,1, rep(0, 3))))
base_14_contrast <- t(matrix(c(0,0,0, rep(0, 3))))
zero_14_contrast <- t(matrix(c(0,-1,1, rep(0, 3))))
# base_zero_contrast <- t(matrix(c(0,0,0,0,rep(0, 2), 1)))

symp_only_purff <- clean_data %>%
  filter(class=="S", timepoint!="day28", timepoint!="bad_baseline")%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  mutate(log_qpcr=log10(qpcr+0.1))%>%
  group_by(targetName)%>%
  nest() %>%
  # mutate(model=map(data, ~lm(concentration~timepoint+ageyrs+id, data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~timepoint+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_zero=map(model, ~multcomp::glht(., base_zero_contrast)),
         base_zero_p=map_dbl(base_zero, ~summary(.)$test$pvalues)) %>%
  mutate(base_14=map(model, ~multcomp::glht(., base_14_contrast)),
         base_14_p=map_dbl(base_14, ~summary(.)$test$pvalues)) %>%
  mutate(base_7=map(model, ~multcomp::glht(., base_7_contrast)),
         base_7_p=map_dbl(base_7, ~summary(.)$test$pvalues))%>%
  ungroup()%>%
  mutate(base_zero_padj=p.adjust(base_zero_p, method="BH"),
         base_7_padj=p.adjust(base_7_p, method="BH"),
         base_14_padj=p.adjust(base_14_p, method="BH")
  )


symp_results_table <- symp_only_purff %>%
  dplyr::select(targetName,
                base_zero_padj,
                base_7_padj,
                base_14_padj)%>%
  ungroup()

#79 with lm(); 159 with lme4;
symp_sig_base_zero <- symp_results_table %>%
  filter(base_zero_padj<fdr_cutoff)%>%
  arrange(base_zero_padj)%>%
  select(targetName)

#0 with lm(); 36 with lme4
symp_sig_base_7 <- symp_results_table %>%
  filter(base_7_padj<fdr_cutoff)%>%
  arrange(base_7_padj)%>%
  select(targetName)

#2 with lme4
symp_sig_base_14 <- symp_results_table %>%
  filter(base_14_padj<fdr_cutoff)%>%
  select(targetName)



clean_data %>%
  filter(class=="S", timepoint%notin%c("day28"), timepoint!="bad_baseline")%>%
  filter(targetName %in% c("IL21", "IFNG"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
  # geom_line(aes(group=id))+
  # geom_point()+#
  geom_boxplot(outliers = FALSE)+
  # geom_violin(draw_quantiles = c(0,0.25,0.5,0.75), color="white")+
  ggtitle("symptomatic malaria")+
  facet_wrap(~targetName, scales = "free")+
  scale_fill_manual(values=viridis::magma(5))+
  theme_minimal()


for(i in 1:ceiling(nrow(symp_sig_base_zero)/16)){
  
plt <- clean_data %>%
  filter(class=="S", timepoint!="day28", timepoint!="bad_baseline")%>%
  filter(targetName %in% symp_sig_base_zero$targetName[seq((i-1)*16+1, i*16)])%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
  # geom_point()+#
  # geom_line(aes(group=id))+
  geom_boxplot(outliers = FALSE)+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  ggtitle("regulated during malaria")+
  facet_wrap(~targetName, scales = "free")+
  scale_fill_manual(values=viridis::magma(5))+
  theme_minimal()+
  da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/symp_base_zero", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
  }


for(i in 1:ceiling(nrow(symp_sig_base_7)/16)){
  
  plt <- clean_data %>%
    filter(class=="S", timepoint!="day28", timepoint!="bad_baseline")%>%
    filter(targetName %in% symp_sig_base_7$targetName[seq((i-1)*16+1, i*16)])%>%
    mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
    ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
    # geom_point()+#
    # geom_line(aes(group=id))+
    geom_boxplot(outliers = FALSE)+
    # geom_violin(draw_quantiles = seq(0,1,0.25))+
    ggtitle("regulated on day 7 post malaria")+
    facet_wrap(~targetName, scales = "free")+
    scale_fill_manual(values=viridis::magma(5))+
    theme_minimal()+
    da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/symp_base_7", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}

analytes_da_7_not_0 <- symp_sig_base_7$targetName[(symp_sig_base_7$targetName %notin% symp_sig_base_zero$targetName)]

analytes_da_7_not_0_plot <- clean_data %>%
  filter(class=="S", timepoint!="day28", timepoint!="bad_baseline")%>%
  filter(targetName %in% analytes_da_7_not_0)%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
  # geom_point()+#
  # geom_line(aes(group=id))+
  geom_boxplot(outliers = FALSE)+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  ggtitle("regulated on day 7 post malaria")+
  facet_wrap(~targetName, scales = "free", ncol=4)+
  scale_fill_manual(values=viridis::magma(5))+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/symp_base_7_not_0.png", analytes_da_7_not_0_plot, height=4, width=8, dpi=444, bg="white")


## asymptomatic only ####
# base_zero_contrast <- t(matrix(c(0,1,0,0, rep(0,44))))
# base_14_contrast <- t(matrix(c(0,0,1,0, rep(0,44))))
# # base_28_contrast <- t(matrix(c(0,0,0)))
# zero_14_contrast <- t(matrix(c(0,-1,1,0, rep(0,44))))

base_zero_contrast <- t(matrix(c(0,1,0, rep(1,3))))
base_14_contrast <- t(matrix(c(0,0,1, rep(1,3))))
# base_28_contrast <- t(matrix(c(0,0,0)))
zero_14_contrast <- t(matrix(c(0,-1,1, rep(0,3))))

asymp_only_purff <- clean_data %>%
  filter(class=="A", timepoint!="day28", timepoint!="bad_baseline")%>%
  mutate("log_qpcr"=log10(qpcr+0.1))%>%
  group_by(targetName)%>%
  nest() %>%
  mutate(model=map(data, ~lme4::lmer(concentration~timepoint+ageyrs+gender_categorical+log_qpcr+(1|id), data=.))) %>%
  # mutate(model=map(data, ~lm(concentration~timepoint+ageyrs+id, data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_zero=map(model, ~multcomp::glht(., base_zero_contrast)),
         base_zero_p=map_dbl(base_zero, ~summary(.)$test$pvalues)) %>%
  mutate(base_14=map(model, ~multcomp::glht(., base_14_contrast)),
         base_14_p=map_dbl(base_14, ~summary(.)$test$pvalues)) %>%
  mutate(zero_14=map(model, ~multcomp::glht(., zero_14_contrast)),
         zero_14_p=map_dbl(zero_14, ~summary(.)$test$pvalues))%>%
  # mutate(base_28=map(model, ~multcomp::glht(., base_28_contrast)),
  #        base_28_p=map_dbl(base_28, ~summary(.)$test$pvalues))%>%
  ungroup()%>%
  mutate(base_zero_padj=p.adjust(base_zero_p, method="BH"),
         base_14_padj=p.adjust(base_14_p, method="BH"),
         zero_14_padj=p.adjust(zero_14_p, method="BH"),
         # base_28_padj=p.adjust(base_28_p, method="BH")
  )


asymp_results_table <- asymp_only_purff %>%
  dplyr::select(targetName,
                base_zero_padj,
                base_14_padj,
                zero_14_padj
  )%>%
  ungroup()

#0 28
asymp_sig_base_zero <- asymp_results_table %>%
  filter(base_zero_padj<fdr_cutoff)%>%
  arrange(base_zero_padj)%>%
  select(targetName)

#26
asymp_sig_base_14 <- asymp_results_table %>%
  filter(base_14_padj<fdr_cutoff)%>%
  arrange(base_14_padj)%>%
  select(targetName)

#0 
asymp_sig_zero_14 <- asymp_results_table %>%
  filter(zero_14_padj<fdr_cutoff)%>%
  arrange(zero_14_padj)%>%
  select(targetName)


for(i in 1:ceiling(nrow(asymp_sig_base_zero)/16)){
  
  plt <- clean_data %>%
    filter(class=="A", timepoint!="day28", timepoint!="bad_baseline")%>%
    filter(targetName %in% asymp_sig_base_zero$targetName[seq((i-1)*16+1, i*16)])%>%
    mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
    ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
    # geom_point()+#
    # geom_line(aes(group=id))+
    geom_boxplot(outliers = FALSE)+
    # geom_violin(draw_quantiles = seq(0,1,0.25))+
    ggtitle("regulated during asymptomatic parasitemia")+
    facet_wrap(~targetName, scales = "free")+
    scale_fill_manual(values=viridis::magma(5))+
    theme_minimal()+
    da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/asymp_base_zero", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}


for(i in 1:ceiling(nrow(asymp_sig_base_14)/16)){
  
  plt <- clean_data %>%
    filter(class=="A", timepoint!="day28", timepoint!="bad_baseline")%>%
    filter(targetName %in% asymp_sig_base_14$targetName[seq((i-1)*16+1, i*16)])%>%
    mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
    ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
    # geom_point()+#
    # geom_line(aes(group=id))+
    geom_boxplot(outliers = FALSE)+
    # geom_violin(draw_quantiles = seq(0,1,0.25))+
    ggtitle("regulated on day 14 of asymptomatic parasitemia")+
    facet_wrap(~targetName, scales = "free")+
    scale_fill_manual(values=viridis::magma(5))+
    theme_minimal()+
    da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/asymp_base_14", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}

# 
base_zero_asymp_plot <- clean_data %>%
  filter(class=="A", timepoint!="day28", timepoint!="bad_baseline")%>%
  filter(targetName %in% c("IL10", "LAG3", "GZMA", "PDCD1", "IL2RA", "CTLA4"))%>%
  ggplot(., aes(x=timepoint, y=concentration, fill=timepoint))+
  geom_boxplot()+
  ggtitle("regulated during asymptomatic parasitemia")+
  facet_wrap(~targetName, scales = "free", ncol = 3)+
  scale_fill_manual(values=c(viridis::magma(4)))+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/select_asymp_base_zero.png", base_zero_asymp_plot, width=8, height=6, bg="white", dpi=444)

# base_14_asymp_plot <- clean_data %>%
#   filter(class=="A", timepoint!="day28", timepoint!="bad_baseline")%>%
#   filter(targetName %in% asymp_sig_base_14$targetName)%>%
#   ggplot(., aes(x=timepoint, y=concentration, fill=timepoint))+
#   geom_boxplot()+
#   ggtitle("upregulated at day 14 of asymptomatic parasitemia")+
#   facet_wrap(~targetName, scales = "free", ncol = 4)+
#   scale_fill_manual(values=c(viridis::magma(4)))+
#   theme_minimal()+
#   da_boxplot_theme
# 
# ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/asymp_base_14.png", base_14_asymp_plot, width=8, height=4, bg="white", dpi=444)

## baselines only ####
base_base_contrast <- t(matrix(c(0,1)))


base_only_purff <- clean_data %>%
  filter(class %in% c("S", "A"), timepoint=="baseline")%>%
  group_by(targetName)%>%
  mutate(bino_class=if_else(class=="S", 1, 0))%>%
  nest() %>%
  # mutate(model=map(data, ~lm(concentration~class+id, data=.))) %>%
  # mutate(model=map(data, ~lme4::lmer(concentration~class+(1|id), data=.))) %>%
  # this one "works"
  mutate(model=map(data, ~glm(bino_class~concentration+id, data=., family="binomial"))) %>%
  # mutate(model=map(data, ~lme4::lmer(concentration~class+(1|id), data=.)))%>%
  # mutate(model=map(data, ~lme4::lmer(concentration~class+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  # mutate(base_base=map(model, ~multcomp::glht(., base_base_contrast)),
  #        base_base_p=map_dbl(base_base, ~summary(.)$test$pvalues)) %>%
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

#18 with glm
base_only_sig <- base_only_results_table %>%
  filter(base_base_padj<fdr_cutoff)%>%
  select(targetName)


base_only_plot <- clean_data %>%
  filter(timepoint=="baseline", class %in% c("S", "A"))%>%
  # filter(targetName %in% c("CXCL8", "MMP8", "MMP9", "MP9", "OSM", "PTX3", "TNFSF11"))%>%
  filter(targetName %in% base_only_sig$targetName)%>%
  ggplot(aes(x=class, y=concentration, fill=targetName))+
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
  filter(timepoint=="baseline", class %in% c("S", "A"))%>%
  filter(targetName %in% c("CXCL8", "MMP8", "MMP9", "MP9", "OSM", "PTX3", "TNFSF11"))%>%
  # filter(targetName %in% base_only_sig$targetName)%>%
  ggplot(aes(x=class, y=concentration, fill=class))+
  # geom_violin(draw_quantiles = c(seq(0,1,0.25)), colour = "darkgrey")+
  geom_boxplot(outliers = FALSE)+
  ggtitle("baselines before parasitemia")+
  facet_wrap(~targetName, scales = "free", ncol=6)+
  scale_fill_manual(values=c(viridis::magma(5)))+
  theme_minimal()+
  theme(legend.position = "none")+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/neutrophil_base_only.png", neutrophil_base_only_plot, width=8, height=6, bg="white", dpi=444)


clean_data %>%
  filter(timepoint=="baseline", class %in% c("S", "A"))%>%
  # filter(targetName %in% c("CXCL8", "MMP8", "MMP9", "MP9", "OSM", "PTX3", "TNFSF11"))%>%
  filter(targetName %in% base_only_sig$targetName)%>%
  ggplot(aes(x=class, y=concentration, fill=gender_categorical))+
  geom_boxplot()+
  # geom_violin(draw_quantiles = c(seq(0,1,0.25)), colour = "darkgrey")+
  # geom_point()+
  # geom_line(aes(group=id))+
  ggtitle("baselines before parasitemia")+
  facet_wrap(~targetName, scales = "free")+
  # scale_fill_manual(values=c(viridis::magma(10)))+
  theme_minimal()+
  theme()


## nmf ####
base_zero_contrast <- t(matrix(c(0,1,0,0)))
base_7_contrast <- t(matrix(c(0,0,1,0)))
base_14_contrast <- t(matrix(c(0,0,0,1)))
zero_14_contrast <- t(matrix(c(0,-1,0,1)))



nmf_only_purff <- clean_data %>%
  filter(class=="NM", timepoint!="day28")%>%
  group_by(targetName)%>%
  nest() %>%
  mutate(model=map(data, ~lme4::lmer(concentration~timepoint+(1|id), data=.))) %>%
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
                zero_14_padj)%>%
  ungroup()

nmf_sig_base_zero <- nmf_results_table %>%
  filter(base_zero_padj<fdr_cutoff)%>%
  arrange(desc(base_zero_padj))%>%
  select(targetName)

nmf_sig_base_7 <- nmf_results_table %>%
  filter(base_7_padj<fdr_cutoff)%>%
  select(targetName)

nmf_sig_base_14 <- nmf_results_table %>%
  filter(base_14_padj<fdr_cutoff)%>%
  select(targetName)

nmf_sig_zero_14 <- nmf_results_table %>%
  filter(zero_14_padj<fdr_cutoff)%>%
  select(targetName)




nmf_plot <- clean_data %>%
  filter(targetName %in% c("CX3CL1", "CCL2", "CXCL10", "IL10"), class %in% c("NM"), timepoint!="day28")%>%
  ggplot(aes(x=timepoint, y=concentration, fill=interaction(class, timepoint)))+
  geom_boxplot()+
  geom_point(position = position_dodge(width=0.75))+# geom_boxplot(outliers = FALSE)+
  ggtitle("non-malarial fever")+
  facet_wrap(~targetName, scales = "free")+
  scale_fill_manual(values=c(viridis::magma(4)))+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/nmf_plot.png", nmf_plot, width=4, height=4, bg="white", dpi=444)



# predicting parasitemia ####

base_zero_contrast <- t(matrix(c(0,1,0,0)))
base_7_contrast <- t(matrix(c(0,0,1,0)))
base_14_contrast <- t(matrix(c(0,0,0,1)))
zero_14_contrast <- t(matrix(c(0,-1,0,1)))



para_only_purff <- clean_data %>%
  filter(class%in%c("S"), timepoint=="baseline", timepoint!="bad_baseline", !is.na(qpcr_cat))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  group_by(targetName)%>%
  nest() %>%
  mutate(model=map(data, ~lm(day0_qpcr~concentration, data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(para_p=map_dbl(summary, ~coef(.)[8])) %>%
  ungroup()%>%
  mutate(para_padj=p.adjust(para_p, method="BH")
  )


para_results_table <- para_only_purff %>%
  dplyr::select(targetName,
                para_padj)%>%
  ungroup()

sig_para <- para_results_table %>%
  # filter(para_padj<fdr_cutoff)%>%
  arrange(para_padj)%>%
  select(targetName)

# model selection ####
  ## lm() covariate inclusion####
model_selection <- clean_data %>%
  filter(timepoint!="day28", timepoint!="bad_baseline")%>%
  mutate(log_qpcr=log10(qpcr+0.1))%>%
  group_by(targetName)%>%
  nest() %>%
  mutate(simple_model=map(data,  ~lm(concentration~timepoint+class+id, data=.)), simple_AIC=map_dbl(simple_model, ~AIC(.)))%>%
  mutate(qpcr_model=map(data,    ~lm(concentration~timepoint+class+id+log_qpcr, data=.)), qpcr_AIC=map_dbl(qpcr_model, ~AIC(.)))%>%
  mutate(age_model=map(data,     ~lm(concentration~timepoint+class+id+ageyrs, data=.)), age_AIC=map_dbl(age_model, ~AIC(.)))%>%
  mutate(sex_model=map(data,     ~lm(concentration~timepoint+class+id+gender_categorical, data=.)), sex_AIC=map_dbl(sex_model, ~AIC(.)))%>%
  mutate(age_sex_model=map(data, ~lm(concentration~timepoint+class+id+gender_categorical+ageyrs, data=.)), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))%>%
  mutate(simple_model_summary=map(simple_model, ~summary(.))) %>%
  mutate(sex_model_summary=map(sex_model, ~summary(.))) %>%
  rowwise()%>%
  mutate("lowest"=min(simple_AIC, age_AIC, sex_AIC, age_sex_AIC, qpcr_AIC))%>%
  mutate("best_model"=list("simple_AIC", "age_AIC", "sex_AIC", "age_sex_AIC", "qpcr_AIC")[which(c(simple_AIC, age_AIC, sex_AIC, age_sex_AIC, qpcr_AIC)==lowest)])%>%
  mutate("better_than_base_model"=(lowest+4)<simple_AIC)%>%
  muate("aic_diff"=)



table(mixed_model_selection$best_model, mixed_model_selection$better_than_base_model==TRUE) #101
age

mixed_model_selection$targetName[mixed_model_selection$best_model=="qpcr_AIC" & mixed_model_selection$better_than_base_model==TRUE]

# --> inlcude age as covariate

## lm() model structure####

model_selection <- clean_data %>%
  filter(timepoint!="day28", timepoint!="bad_baseline")%>%
  group_by(targetName)%>%
  nest() %>%
  mutate(additive_model=map(data,  ~lm(concentration~timepoint+class+ageyrs+id, data=.)), additive_model_AIC=map_dbl(additive_model, ~AIC(.)))%>%
  mutate(age_interact=map(data,    ~lm(concentration~timepoint+class*ageyrs+id, data=.)), age_interact_AIC=map_dbl(age_interact, ~AIC(.)))%>%
  mutate(class_interact=map(data,     ~lm(concentration~timepoint*class+ageyrs+id, data=.)), class_interact_AIC=map_dbl(class_interact, ~AIC(.)))%>%
  # mutate(age_sex_model=map(data, ~lm(concentration~timepoint+class+ageyrs+gender_categorical+id, data=.)), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))%>%
  rowwise()%>%
  mutate(lowest=min(additive_model_AIC, age_interact_AIC, class_interact_AIC))%>%
  ungroup()%>%
  mutate(across(c(additive_model_AIC, age_interact_AIC, class_interact_AIC), ~ .x - lowest))

table(model_selection$additive_model_AIC==0)#129
table(model_selection$age_interact_AIC==0 & model_selection$additive_model_AIC>4)#4
table(model_selection$class_interact_AIC==0 & model_selection$additive_model_AIC>4)#66
table(model_selection$age_sex_AIC==0 & model_selection$age_AIC>4)#47

# --> simple additive model is best


## lmer () ####
mixed_model_selection <- clean_data %>%
  filter(timepoint!="day28", timepoint!="bad_baseline", class%in% c("S", "A"))%>%
  group_by(targetName)%>%
  mutate(log_qpcr=log10(qpcr+0.1))%>%
  nest() %>%
  mutate(simple_model=map(data,  ~lme4::lmer(concentration~timepoint+class+(1|id), data=.)), simple_AIC=map_dbl(simple_model, ~AIC(.)))%>%
  mutate(qpcr_model=map(data,     ~lme4::lmer(concentration~timepoint+class+log_qpcr+(1|id), data=.)), qpcr_AIC=map_dbl(qpcr_model, ~AIC(.)))%>%
  mutate(age_model=map(data,     ~lme4::lmer(concentration~timepoint+class+ageyrs+(1|id), data=.)), age_AIC=map_dbl(age_model, ~AIC(.)))%>%
  mutate(sex_model=map(data,     ~lme4::lmer(concentration~timepoint+class+gender_categorical+(1|id), data=.)), sex_AIC=map_dbl(sex_model, ~AIC(.)))%>%
  mutate(age_sex_model=map(data, ~lme4::lmer(concentration~timepoint+class+ageyrs+gender_categorical+(1|id), data=.)), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))%>%
  rowwise()%>%
  mutate("lowest"=min(simple_AIC, age_AIC, sex_AIC, age_sex_AIC, qpcr_AIC))%>%
  mutate("best_model"=c("simple_AIC", "age_AIC", "sex_AIC", "age_sex_AIC", "qpcr_AIC")[which(c(simple_AIC, age_AIC, sex_AIC, age_sex_AIC, qpcr_AIC)==lowest)])%>%
  mutate("better_than_base_model"=(lowest+4)<simple_AIC)
  # ungroup()%>%
  # mutate(across(c(simple_AIC, age_AIC, sex_AIC, age_sex_AIC, qpcr_AIC), ~ .x - lowest))
  
table(mixed_model_selection$best_model, mixed_model_selection$better_than_base_model==TRUE) #101
#               FALSE TRUE
# age_AIC        15   19
# age_sex_AIC     0    1
# qpcr_AIC       22   79
# sex_AIC        14    3
# simple_AIC     97    0

for(i in 1:nrow(mixed_model_selection)){
  
  models <- list(
    "simple model" = mixed_model_selection$simple_model[[i]],
    "qpcr model" = mixed_model_selection$qpcr_model[[i]],
    "age model" = mixed_model_selection$age_model[[i]],
    "sex model" = mixed_model_selection$sex_model[[i]],
    "age sex model" = mixed_model_selection$age_sex_model[[i]]
  )
  
  plt <- ggcoef_compare(models, type = "faceted")
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/model_plots/mixed_model_selection", mixed_model_selection$targetName[[i]], ".png", sep=""), plt, height=6, width=18, dpi=444, bg="white")
  
}

# --> use simple mixed model
# > table(mixed_model_selection$sex_AIC==0 & mixed_model_selection$simple_AIC>4)
# FALSE  TRUE 
# 237    13 
# [1] "ANGPT1"  "CD40LG"  "CRP"     "CTF1"    "GZMB"    "IL12p70" "IL17F"   "IL1RN"   "IL27"    "IRAK4"   "MIF"     "PTX3"    "TNFSF14"

# > table(mixed_model_selection$age_AIC==0 & mixed_model_selection$simple_AIC>4)
# FALSE  TRUE
# 213    37 
# [1] "AGRP"          "BST2"          "CCL17"         "CCL28"         "CD70"          "CD80"          "CSF2"          "CTLA4"        
# [9] "CXCL2"         "EGF"           "FGF2"          "GDF15"         "GFAP"          "IFNA1; IFNA13" "IFNG"          "IKBKG"        
# [17] "IL17RB"        "IL22"          "IL2RA"         "IL33"          "IL7"           "KLRK1"         "LGALS9"        "MUC16"        
# [25] "NCR1"          "NGF"           "S100A9"        "SCG2"          "SELE"          "SELP"          "TEK"           "TGFB1"        
# [33] "THBS2"         "TIMP1"         "TNFSF15"       "VEGFD"         "VSNL1" 


# sandbox ####

# ditching timepoint####
qpcr_class_contrast <- t(matrix(c(0, 0, 0, 0, 0, 1)))
qpcr_contrast <- t(matrix(c(0, 1, 0, 0, 0, 0)))
class_contrast <- t(matrix(c(0, 0, 1, 0, 0, 0)))


combo_as_purff <- clean_data %>%
  filter(class%in%c("A", "S"), timepoint %notin% c("day7", "day28"), timepoint!="bad_baseline")%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")),
         log_qpcr = log10(qpcr+0.1))%>%
  group_by(targetName)%>%
  nest() %>%
  # mutate(model=map(data, ~lme4::lmer(concentration~timepoint*class+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~log_qpcr*class+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(qpcr_class=map(model, ~multcomp::glht(., qpcr_class_contrast)),
         qpcr_class_p=map_dbl(qpcr_class, ~summary(.)$test$pvalues)) %>%
  mutate(qpcr=map(model, ~multcomp::glht(., qpcr_contrast)),
         qpcr_p=map_dbl(qpcr, ~summary(.)$test$pvalues)) %>%
  mutate(class=map(model, ~multcomp::glht(., class_contrast)),
         class_p=map_dbl(class, ~summary(.)$test$pvalues)) %>%
  ungroup()%>%
  mutate(qpcr_class_padj=p.adjust(qpcr_class_p, method="BH"),
         qpcr_padj=p.adjust(qpcr_p, method="BH"),
         class_padj=p.adjust(class_p, method="BH")
  )


results_table <- combo_as_purff %>%
  dplyr::select(targetName,
                qpcr_class_padj,
                qpcr_padj, 
                class_padj)%>%
  ungroup()

# table(results_table$qpcr_class_padj<fdr_cutoff)# 63 
# table(results_table$qpcr_padj<fdr_cutoff)# 0 this is different when you model A on their own
# table(results_table$base_class_padj<fdr_cutoff)# 15 
# table(results_table$qpcr_class_class_padj<fdr_cutoff)# 89 --> ~110 vary by timepoint 
# table(results_table$base_para_padj<fdr_cutoff)# 94 --> ~110 vary by timepoint 

sig_qpcr_class <- results_table%>%
  filter(qpcr_class_padj<fdr_cutoff)%>%
  arrange(qpcr_class_padj)

sig_qpcr <- results_table%>%
  filter(qpcr_padj<fdr_cutoff)%>%
  arrange(qpcr_padj)

sig_class <- results_table%>%
  filter(class_padj<fdr_cutoff)%>%
  arrange(class_padj)


clean_data %>%
  filter(class%in%c("A", "S"), timepoint %notin% c("day7", "day28"), timepoint!="bad_baseline")%>%
  filter(targetName %in% sig_class$targetName)%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=class, y=concentration, fill=class))+
  # geom_point()+#
  # geom_line(aes(group=id))+
  geom_boxplot(outliers = FALSE)+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  # ggtitle("regulated during asymptomatic parasitemia")+
  facet_wrap(~targetName, scales = "free")+
  scale_fill_manual(values=viridis::magma(3))+
  theme_minimal()



clean_data %>%
  filter(class%in%c("A", "S"), timepoint %notin% c("day7", "day28"), timepoint!="bad_baseline")%>%
  filter(targetName %in% sig_qpcr_class$targetName)%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=qpcr+0.1, y=concentration, color=class))+
  # geom_point()+#
  geom_smooth(method = "lm")+
  # geom_line(aes(group=id))+
  # geom_boxplot(outliers = FALSE)+
  scale_x_log10()+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  # ggtitle("regulated during asymptomatic parasitemia")+
  facet_wrap(~targetName, scales = "free")+
  scale_color_manual(values=viridis::magma(3))+
  theme_minimal()




clean_data %>%
  filter(class%in%c("A", "S"), timepoint %notin% c("day7", "day28"), timepoint!="bad_baseline")%>%
  filter(targetName %in% c("IL4", "IL5", "IL9", "IL13", "IL3RA", "CCL17", "CCL22"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=qpcr+0.1, y=concentration, color=class))+
  # geom_point()+#
  geom_smooth(method = "lm")+
  # geom_line(aes(group=id))+
  # geom_boxplot(outliers = FALSE)+
  scale_x_log10()+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  # ggtitle("regulated during asymptomatic parasitemia")+
  facet_wrap(~targetName, scales = "fixed")+
  scale_color_manual(values=viridis::magma(3))+
  theme_minimal()

# misc ####



target_vars <- clean_data %>%
  filter(targetName%notin%c("CTSS"))%>%
  group_by(targetName)%>%
  summarise("mean"=mean(concentration), "var"=var(concentration), "ratio"=var/mean)
 
most_variable <- slice_max(target_vars, ratio, prop = 0.5)
 
clean_data%>%
  filter(targetName %in% most_variable_features$targetName,
         class %notin% c("PS", "A2", "NM"))%>%
  ggplot(aes(x=timepoint, y=concentration, fill=class))+
  geom_violin()+
  facet_wrap(~targetName)+
  theme_minimal()
  



model_selection <- clean_data %>%
  filter(timepoint!="day28", timepoint!="bad_baseline")%>%
  group_by(targetName)%>%
  nest() %>%
  mutate(simple_model=map(data,  ~lm(concentration~timepoint*class+id, data=.)), simple_AIC=map_dbl(simple_model, ~AIC(.)))%>%
  mutate(age_model=map(data,     ~lm(concentration~timepoint*class+ageyrs+id, data=.)), age_AIC=map_dbl(age_model, ~AIC(.)))%>%
  mutate(sex_model=map(data,     ~lm(concentration~timepoint*class+gender_categorical+id, data=.)), sex_AIC=map_dbl(sex_model, ~AIC(.)))%>%
  mutate(age_sex_model=map(data, ~lm(concentration~timepoint*class+ageyrs+gender_categorical+id, data=.)), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))%>%
  rowwise()%>%
  mutate(lowest=min(simple_AIC, age_AIC, sex_AIC, age_sex_AIC))%>%
  ungroup()%>%
  mutate(simplest_is_best=lowest==simple_AIC,
         age_is_best=lowest==age_AIC,
         sex_is_best=lowest==sex_AIC,
         age_sex_is_best=lowest==age_sex_AIC)

model_selection$targetName[model_selection$simplest_is_best]
model_selection$targetName[model_selection$age_is_best]
model_selection$targetName[model_selection$sex_is_best]
model_selection$targetName[model_selection$age_sex_is_best]






#most variable features

topVarProts <- head(order(rowVars(cpm(d_norm)), decreasing = TRUE))



# clean_data <- big_df %>%
#   mutate(id=substr(sample_id, 2,4),
#          class=substr(sample_id, 6,6),
#          barcode= substr(sample_id, nchar(sample_id)-5, nchar(sample_id)))%>%
#   rowwise()%>%
#   mutate("sample_type"=substr(sample_id,
#                               stringi::stri_locate_all(pattern = "_", sample_id, fixed=TRUE)[[1]][1]+1,
#                               stringi::stri_locate_all(pattern = "_", sample_id, fixed=TRUE)[[1]][2]-1))%>%
#   ungroup()%>%
#   mutate("timepoint"=case_when(sample_type=="A.1"~"baseline",
#                                sample_type=="A.2"~"drop",
#                                sample_type=="A0"~"day0",
#                                sample_type=="A14"~"day14",
#                                sample_type=="A2.1"~"drop",
#                                sample_type=="A20"~"day0",
#                                sample_type=="A214"~"day14",
#                                sample_type=="A28"~"day28",
#                                sample_type=="NM.1"~"baseline",
#                                sample_type=="NM0"~"day0",
#                                sample_type=="NM14"~"day14",
#                                sample_type=="NM7"~"day7",
#                                sample_type=="PS0"~"PS0",
#                                sample_type=="S.1"~"baseline",
#                                sample_type=="S.2"~"drop",
#                                sample_type=="S0"~"day0",
#                                sample_type=="S14"~"day14",
#                                sample_type=="S28"~"day28",
#                                sample_type=="S7"~"day7")
#   )%>%
#   mutate(timepoint=factor(timepoint, levels=c("baseline", "day0", "day7", "day14", "day28")))%>%
#   group_by(targetName) %>%
#   mutate(z_conc=scale(concentration, center = TRUE, scale = TRUE))%>%
#   ungroup()%>%
#   group_by(id)%>%
#   mutate("mean_z_conc"=median(z_conc))%>%
#   ungroup()%>%
#   filter(mean_z_conc > -0.35)
# parasitemia inclussion?


base_zero_contrast <- t(matrix(c(0,1,0,0)))
base_14_contrast <- t(matrix(c(0,0,1,0)))
# base_28_contrast <- t(matrix(c(0,0,0)))
zero_14_contrast <- t(matrix(c(0,-1,1,0)))

base_zero_para_contrast <- t(matrix(c(0,1,0,1)))
base_14_para_contrast <- t(matrix(c(0,0,1,1)))

asymp_only_purff <- clean_data %>%
  filter(class=="A", timepoint!="day28", timepoint!="bad_baseline")%>%
  group_by(targetName)%>%
  mutate(log_pcr=log10(qpcr+0.1),
         log_pcr_dich=if_else(log_pcr>=2.69, "over 500k/ml", "under 500k/ml"))%>%
  nest() %>%
  # mutate(model=map(data, ~lm(concentration~timepoint+ageyrs+id, data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~timepoint+log_pcr+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_zero=map(model, ~multcomp::glht(., base_zero_contrast)),
         base_zero_p=map_dbl(base_zero, ~summary(.)$test$pvalues)) %>%
  mutate(base_zero_para=map(model, ~multcomp::glht(., base_zero_para_contrast)),
         base_zero_para_p=map_dbl(base_zero_para, ~summary(.)$test$pvalues),
         base_zero_para_coef=map_dbl(summary, ~coef(.)[5]),
         base_zero_para_coef_p=map_dbl(summary, ~coef(.)[5]))%>%
  mutate(base_14_para=map(model, ~multcomp::glht(., base_14_para_contrast)),
         base_14_para_p=map_dbl(base_14_para, ~summary(.)$test$pvalues),
         base_14_para_coef=map_dbl(summary, ~coef(.)[5]),
         base_14_para_coef_p=map_dbl(summary, ~coef(.)[5]))%>%
  ungroup()%>%
  mutate(base_zero_padj=p.adjust(base_zero_p, method="BH"),
         base_14_para_padj=p.adjust(base_14_para_p, method="BH"),
         base_zero_para_padj=p.adjust(base_zero_para_p, method="BH")
  )



asymp_results_table <- asymp_only_purff %>%
  dplyr::select(targetName,
                base_zero_padj,
                base_zero_para_padj,
                base_zero_para_coef,
                base_14_para_padj)%>%
  ungroup()

#79 with lm(); 159 with lme4;
asymp_sig_base_zero <- asymp_results_table %>%
  filter(base_zero_padj<fdr_cutoff)


#0 with lm(); 36 with lme4
asymp_sig_base_0_para <- asymp_only_purff %>%
  filter(base_zero_para_padj<fdr_cutoff)%>%
  arrange(desc(base_zero_para_coef))

asymp_sig_base_14_para <- asymp_only_purff %>%
  filter(base_14_para_padj<fdr_cutoff)%>%
  arrange(desc(base_14_para_padj))


#2 with lme4
symp_sig_base_14 <- symp_results_table %>%
  filter(base_14_padj<fdr_cutoff)%>%
  select(targetName)




clean_data %>%
  filter(class=="A", timepoint!="day28", timepoint!="bad_baseline")%>%
  filter(targetName %in% head(asymp_sig_base_14_para$targetName, n=29))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  mutate(log_pcr=log10(qpcr+0.1),
         log_pcr_dich=if_else(log_pcr>=2.69, "over 500k/ml", "under 500k/ml"))%>%
  ggplot(aes(x=log_pcr, y=concentration))+
  geom_point(aes(color=timepoint))+#
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman", color="red")+
  # geom_line(aes(group=id, color=qpcr_cat))+
  # geom_boxplot(aes(fill=log_pcr_dich), outliers = FALSE)+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  # ggtitle("regulated during asymptomatic parasitemia")+
  facet_wrap(~targetName, scales = "free")+
  scale_color_manual(values=viridis::magma(6))+
  theme_minimal()




# examining effect of covariates on timepoint ####
## model selection ####

mixed_model_selection <- clean_data %>%
  filter(timepoint!="day28", timepoint!="bad_baseline", class%in% c("S"))%>%
  group_by(targetName)%>%
  mutate(log_qpcr=log10(qpcr+0.1))%>%
  nest() %>%
  mutate(simple_model=map(data,  ~lme4::lmer(concentration~timepoint+(1|id), data=.)), simple_AIC=map_dbl(simple_model, ~AIC(.)))%>%
  # mutate(qpcr_model=map(data,     ~lme4::lmer(concentration~timepoint+class+log_qpcr+(1|id), data=.)), qpcr_AIC=map_dbl(qpcr_model, ~AIC(.)))%>%
  mutate(age_model=map(data,     ~lme4::lmer(concentration~timepoint+ageyrs+(1|id), data=.)), age_AIC=map_dbl(age_model, ~AIC(.)))%>%
  mutate(sex_model=map(data,     ~lme4::lmer(concentration~timepoint+gender_categorical+(1|id), data=.)), sex_AIC=map_dbl(sex_model, ~AIC(.)))%>%
  mutate(age_sex_model=map(data, ~lme4::lmer(concentration~timepoint+ageyrs+gender_categorical+(1|id), data=.)), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))%>%
  mutate(summary=map(simple_model, ~summary(.))) %>%
  mutate(simple_day0_coef=map(simple_model,  ~coef(summary(.))[2]), simple_AIC=map_dbl(simple_model, ~AIC(.)))%>%
  # mutate(qpcr_day0_coef=map(data,     ~lme4::lmer(concentration~timepoint+class+log_qpcr+(1|id), data=.)), qpcr_AIC=map_dbl(qpcr_model, ~AIC(.)))%>%
  mutate(age_day0_coef=map(age_model,     ~coef(summary(.))[2]), age_AIC=map_dbl(age_model, ~AIC(.)))%>%
  mutate(sex_day0_coef=map(sex_model,     ~coef(summary(.))[2]), sex_AIC=map_dbl(sex_model, ~AIC(.)))%>%
  mutate(age_sex_day0_coef=map(age_sex_model, ~coef(summary(.))[2]), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))%>%
  rowwise()%>%
  mutate("lowest"=min(simple_AIC, age_AIC, sex_AIC, age_sex_AIC))%>%
  mutate("best_model"=c("simple_AIC", "age_AIC", "sex_AIC", "age_sex_AIC")[which(c(simple_AIC, age_AIC, sex_AIC, age_sex_AIC)==lowest)])%>%
  mutate("better_than_base_model"=(lowest+4)<simple_AIC)



slim_mixed_model_selection <- mixed_model_selection %>%
  select(ends_with("_AIC"), ends_with("coef"))%>%
  rowwise()%>%
  mutate("lowest_coef"=min(simple_day0_coef, age_day0_coef, sex_day0_coef, age_sex_day0_coef))%>%
  mutate("highest_coef"=max(simple_day0_coef, age_day0_coef, sex_day0_coef, age_sex_day0_coef))%>%
  mutate("range"=lowest_coef/highest_coef)%>%
  filter()

# for S0, there are only ~10 that vary by 20% or more and the coefs are < | 0.1 |
# for S7 there are only ~10 that vary by 20% or more and the coefs are < | 0.1 |
# for A0 there are ~30 that vary by 20% or more, ~5 that vary by more than 2 fold  
# for A14 there are ~60 that vary by 20% or more, ~10 2 fold, couple of sign changes

# ungroup()%>%
# mutate(across(c(simple_AIC, age_AIC, sex_AIC, age_sex_AIC, qpcr_AIC), ~ .x - lowest))
# 
# table(mixed_model_selection$best_model, mixed_model_selection$better_than_base_model==TRUE) #101



mixed_model_selection <- clean_data %>%
  filter(timepoint!="day28", timepoint!="bad_baseline", class%in% c("S"))%>%
  group_by(targetName)%>%
  mutate(log_qpcr=log10(qpcr+0.1))%>%
  nest() %>%
  mutate(simple_model=map(data,  ~lme4::lmer(concentration~timepoint+(1|id), data=.)), simple_AIC=map_dbl(simple_model, ~AIC(.)))%>%
  mutate(qpcr_model=map(data,     ~lme4::lmer(concentration~timepoint+day0_qpcr+(1|id), data=.)), qpcr_AIC=map_dbl(qpcr_model, ~AIC(.)))%>%
  mutate(qpcr_i_model=map(data,     ~lme4::lmer(concentration~timepoint*day0_qpcr+(1|id), data=.)), qpcr_AIC=map_dbl(qpcr_model, ~AIC(.)))%>%
  # mutate(age_model=map(data,     ~lme4::lmer(concentration~timepoint+ageyrs+(1|id), data=.)), age_AIC=map_dbl(age_model, ~AIC(.)))%>%
  # mutate(sex_model=map(data,     ~lme4::lmer(concentration~timepoint+gender_categorical+(1|id), data=.)), sex_AIC=map_dbl(sex_model, ~AIC(.)))%>%
  # mutate(age_sex_model=map(data, ~lme4::lmer(concentration~timepoint+ageyrs+gender_categorical+(1|id), data=.)), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))%>%
  mutate(summary=map(qpcr_model, ~summary(.))) %>%
  mutate(simple_day0_coef=map(simple_model,  ~coef(summary(.))[2]), simple_AIC=map_dbl(simple_model, ~AIC(.)))%>%
  mutate(qpcr_day0_coef=map(qpcr_model,     ~coef(summary(.))[2]), qpcr_AIC=map_dbl(qpcr_model, ~AIC(.)))%>%
  mutate(qpcr_day0_i_coef=map(qpcr_i_model,     ~coef(summary(.))[2]), qpcr_i_AIC=map_dbl(qpcr_i_model, ~AIC(.)))%>%
  rowwise()%>%
  mutate("lowest"=min(simple_AIC, qpcr_AIC, qpcr_i_AIC))%>%
  mutate("best_model"=c("simple_AIC", "qpcr_AIC", "qpcr_i_AIC")[which(c(simple_AIC, qpcr_AIC, qpcr_i_AIC)==lowest)])%>%
  mutate("better_than_base_model"=(lowest+4)<simple_AIC)



slim_mixed_model_selection <- mixed_model_selection %>%
  select(ends_with("_AIC"), ends_with("coef"))%>%
  rowwise()%>%
  mutate("lowest_coef"=min(simple_day0_coef, qpcr_day0_coef, qpcr_day0_i_coef))%>%
  mutate("highest_coef"=max(simple_day0_coef, qpcr_day0_coef, qpcr_day0_i_coef))%>%
  mutate("log2FC"=log2(abs(lowest_coef/highest_coef)))%>%
  mutate("lowest"=min(simple_AIC, qpcr_AIC, qpcr_i_AIC))%>%
  mutate("best_model"=list("simple_AIC", "qpcr_AIC", "qpcr_i_AIC")[which(c(simple_AIC, qpcr_AIC, qpcr_i_AIC)==lowest)])%>%
  mutate("better_than_base_model"=(lowest+4)<simple_AIC)


table(abs(slim_mixed_model_selection$log2FC)>log2(1.2)) 
table(mixed_model_selection$best_model, mixed_model_selection$better_than_base_model==TRUE) #101
# FALSE TRUE
# qpcr_AIC      40  160
# qpcr_i_AIC     4   20
# simple_AIC    26    0

# mixed_model_selection$targetName[mixed_model_selection$best_model=="qpcr_i_AIC"]
# [1] "CCL2"      "CCL20"     "CCL3"      "CCL4"      "CD80"      "CSF1R"     "CXCL11"    "CXCL16"   
# [9] "FLT1"      "ICAM1"     "IL10"      "IL1B"      "IL1RN"     "IL34"      "IL6"       "PDGFB"    
# [17] "PTX3"      "TNF"       "TNFRSF11A" "TNFRSF1A"  "TNFRSF1B"  "TREM1"     "TSLP"      "VSTM1" 

#substantial change in day0 coef: log qpcr  #240
#substantial change in day0 coef: day0 qpcr #220 

