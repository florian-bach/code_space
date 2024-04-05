modelable_antigens <- c("TT", "SBP1", "Rh5", "SEA", "AMA1", "Hyp2", "HSP40", "GST", "GEXP", "CSP")



(birth_prot <- combo_data %>%
  filter(antigen %in% c("SEA", "HSP40"), is.finite(conc), timepoint==1)%>%
  ggplot(., aes(x=factor(inf_0_12), y=conc))+
  geom_hline(data=filter(cutoff_df, antigen %in% c("SEA", "HSP40")), aes(yintercept = max_below_standard), linetype="dashed")+
  geom_point(alpha=0.2, shape=21, position = position_dodge(width=0.75))+
  geom_boxplot(aes(fill=factor(inf_0_12)), outlier.shape = NA)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow = 1)+
  ylab("Concentration")+
  scale_y_log10(labels=scales::label_log(), breaks=10^seq(-5, 2), limits=10^c(-5, 2.5))+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1),
    axis.title.x = element_blank(),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols))
  
 combo_data %>%
  filter(antigen %in% antigens_exposure_6)%>%
  ggplot(., aes(x=factor(timepointf), y=conc))+
  geom_hline(data=filter(cutoff_df, antigen %in% antigens_exposure_6), aes(yintercept = log10(max_below_standard)), linetype="dashed")+
  geom_point(alpha=0.2, shape=21, position = position_dodge(width=0.75))+
  geom_boxplot(aes(fill=factor(antigen)), outlier.shape = NA)+
  facet_grid(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  ggtitle("Parasitaemia in The First Six Months of Life Leads to The Production of\nAntimalarial Antibodies in Infants")+
  xlab("Number of Parasitaemic Episodes Months 0-6")+
  ylab("Concentration")+
  scale_y_log10(labels=scales::label_log(), breaks=10^seq(-5, 2), limits=10^c(-5, 2.5))+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    axis.text.x = element_text(),
    strip.text = element_text())+
  scale_fill_manual(values=pc1_cols)


ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/yes_inf_six.png", six_inf, height=3, width=8, bg="white", dpi=444)


twelve_inf <- combo_data %>%
  filter(timepointf == "12 Months", antigen %in% placental_antigens)%>%
  ggplot(., aes(x=factor(inf_6_12), y=conc))+
  geom_hline(data=filter(cutoff_df, antigen %in% placental_antigens), aes(yintercept = log10(max_below_standard)), linetype="dashed")+
  geom_point(alpha=0.2, shape=21, position = position_dodge(width=0.75))+
  geom_boxplot(aes(fill=factor(inf_6_12)), outlier.shape = NA)+
  facet_grid(timepointf~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  ggtitle("Parasitaemia in The First Six Months of Life Leads to The Production of\nAntimalarial Antibodies in Infants")+
  xlab("Number of Parasitaemic Episodes Months 6 to 12")+
  ylab("Concentration")+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    axis.text.x = element_text(),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)


ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/yes_inf_twelve.png", twelve_inf, height=3, width=8, bg="white", dpi=444)







a <- combo_data %>%
  filter(timepointf == "6 Months", antigen %in% fave_antigens)%>%
  ggplot(., aes(x=factor(inf_0_6), y=conc))+
  geom_hline(data=filter(cutoff_df, antigen %in% fave_antigens), aes(yintercept = log10(max_below_standard)), linetype="dashed")+
  geom_point(alpha=0.2, shape=21, position = position_dodge(width=0.75))+
  geom_boxplot(aes(fill=factor(inf_0_6)), outlier.shape = NA)+
  facet_grid(timepointf~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  # ggtitle("Parasitaemia in The First Six Months of Life Leads to The Production of\nAntimalarial Antibodies in Infants")+
  xlab("Number of Parasitaemic Episodes Months 0-6")+
  ylab("Concentration")+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    axis.text.x = element_text(),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)


ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/yes_inf_six.png", six_inf, height=3, width=8, bg="white", dpi=444)


b <- combo_data %>%
  filter(timepointf == "12 Months", antigen %in% fave_antigens)%>%
  ggplot(., aes(x=factor(inf_6_12), y=conc))+
  geom_hline(data=filter(cutoff_df, antigen %in% fave_antigens), aes(yintercept = log10(max_below_standard)), linetype="dashed")+
  geom_point(alpha=0.2, shape=21, position = position_dodge(width=0.75))+
  geom_boxplot(aes(fill=factor(inf_6_12)), outlier.shape = NA)+
  facet_grid(timepointf~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  # ggtitle("Parasitaemia in The First Six Months of Life Leads to The Production of\nAntimalarial Antibodies in Infants")+
  xlab("Number of Parasitaemic Episodes Months 6 to 12")+
  ylab("Concentration")+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    axis.text.x = element_text(),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)


c <- combo_data %>%
  filter(timepointf == "12 Months", inf_6_12==0, antigen %in% fave_antigens)%>%
  mutate(any_inf_0_6=case_when(any_inf_0_6=="0"~"none",
                               any_inf_0_6=="1"~"one or more"))%>%
  ggplot(., aes(x=factor(any_inf_0_6), y=conc))+
  geom_hline(data=filter(cutoff_df, antigen %in% fave_antigens), aes(yintercept = log10(max_below_standard)), linetype="dashed")+
  geom_point(alpha=0.2, shape=21, position = position_dodge(width=0.75))+
  geom_boxplot(aes(fill=factor(any_inf_0_6, levels = c("none", "one or more"))), outlier.shape = NA)+
  facet_grid(timepointf~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  # ggtitle("Parasitaemia in The First Six Months of a Defect in\nAntimalarial Antibodies in Infants at 12 Months of Age")+
  xlab("Number of Parasitaemic Episodes Months 0 to 6")+
  ylab("Concentration")+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_text(angle = 90, hjust=1),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)

fig1= a / b / c

fig1=fig1+plot_annotation(
  title = 'Infants\' Antibody Responses to Infection',
  subtitle = 'The number of infections in the previous six months influences antibody secretion, but infections in the first six\nmonths of life lead to decreased antiobdy responses at 12 months (though mostly below the limit of detection)'
  # caption = 'Disclaimer: None of these plots are insightful'
)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/fig1_patch.png", fig1, height=8, width=12, bg="white", dpi=444)












# incidence tfh / abc stuff ####


tfh_12_18 <- long_tfh_clinab %>%
  filter(!is.na(symp_12_18), cell_pop %in% c("Th_memory","Th_naive"))%>%
  filter(timepointf=="12 Months")%>%
  ggplot(., aes(x=factor(symp_12_18), y=cell_freq/100))+
  geom_point(alpha=0.2, shape=21)+
  geom_boxplot(aes(fill=factor(symp_12_18)), outlier.shape = NA)+
  facet_wrap(~cell_popf, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow=1)+
  ggtitle("lower T cell memory differentiation is\nassociated with increased future malaria")+
  ylab("% of CD4 T cells")+
  xlab("Number of symptomatic episodes in months 12-18")+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_text(angle = 90, hjust=1),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/t_cell_incidence_12_18.png", tfh_12_18, height=3, width=4, bg="white", dpi=444)





abc_12_18 <- long_abc_clinab %>%
  filter(!is.na(symp_12_18), cell_pop %in% c("naive_b", "memory_b", "immature_b_perc"))%>%
  filter(timepointf=="12 Months")%>%
  ggplot(., aes(x=factor(symp_12_18), y=cell_freq/100))+
  geom_point(alpha=0.2, shape=21)+
  geom_boxplot(aes(fill=factor(symp_12_18)), outlier.shape = NA)+
  facet_wrap(~factor(cell_popf, levels=c("naive B Cells", "memory B Cells", "CD10+ B Cells")), labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow=1)+
  ggtitle("higher memory differentiation in B cells is associated\nwith increased future malaria")+
  ylab("% of parent population")+
  xlab("Number of symptomatic episodes in months 12-18")+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_text(angle = 90, hjust=1),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/b_cell_incidence_12_18.png", abc_12_18, height=3, width=6, bg="white", dpi=444)





list_of_abc_plots <- list()

for(i in 1:nrow(abc_incidence_sigs)){
  
  plot_data <- abc_clin %>%
    # filter(!is.na(antigen), timepoint==3) %>%
    filter(cell_pop==paste(abc_incidence_sigs[i,1]) & incidence_type==paste(abc_incidence_sigs[i,2]))%>%
    filter(!duplicated(id))
  
  plot <- ggplot(plot_data, aes(x=factor(incidence_value), y=cell_freq, fill=factor(incidence_value)))+
    geom_boxplot()+
    geom_point(shape=21)+
    # ggtitle(paste(unique(sigs[i,2])))+
    scale_fill_manual(values=incidence_cols)+
    # geom_point(fill=pc1_cols[i], alpha=1, shape=21)+
    # geom_smooth(method="lm")+
    # scale_y_log10()+
    # ggpubr::stat_cor(method = "spearman", label.y = -1.5, na.rm = TRUE, size=2)+
    ylab(abc_incidence_sigs[i,1])+
    xlab(abc_incidence_sigs[i,2])+
    theme_minimal()+
    theme(legend.position = "none")
  list_of_abc_plots[[i]] <- plot
}

sig_abc_plot <- cowplot::plot_grid(plotlist = list_of_abc_plots, nrow = 3)

