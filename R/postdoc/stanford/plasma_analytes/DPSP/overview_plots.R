library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(emmeans)


dpsp_nulisa <- read.csv("~/postdoc/stanford/plasma_analytes/DPSP/dpsp_nulisa_data.csv")

fc_data <- dpsp_nulisa %>%
  filter(sample!="1094_pregnant_1")%>%
  pivot_wider(names_from = timepoint2, values_from = c(conc), id_cols=c(targetName, id), names_prefix = "conc_")%>%
  group_by(targetName)%>%
  mutate(preg_fc=conc_delivery-conc_pregnant)%>%
  distinct(targetName, id, preg_fc)%>%
  group_by(targetName)%>%
  summarise("mean_preg_fc"=mean(preg_fc, na.rm = T))

neuro_cog <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/NCT_infants_share032025.csv")

maternal_neuro_cog_edit <- neuro_cog%>%
  mutate(id=subjid-10000)

nulisa_plus_neuro <- left_join(maternal_data, maternal_neuro_cog_edit, by="id")



dpsp_results <- dpsp_nulisa %>%
  filter(timepoint2!="post")%>%
  mutate(id_cat=as.character(id))%>%
  group_by(targetName)%>%
  nest() %>%
  mutate(model=map(data, ~lme4::lmer(conc~timepoint2+(1|id_cat), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(emm=map(model, ~emmeans(., specs = pairwise ~ timepoint2)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate("delivery - pregnant p"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("coef delivery - pregnant"=map_dbl(emm_contrast_summary, ~.$estimate[1])) %>%
  ungroup()%>%
  mutate("delivery - pregnant padj" = p.adjust(`delivery - pregnant p`, method="fdr"))

combo_dpsp_results <- left_join(dpsp_results, fc_data, by="targetName")

sig_preg <- combo_dpsp_results %>%
  filter(`delivery - pregnant padj`<0.05 & abs(mean_preg_fc)>0.5)

top_16 <- sig_preg%>%
  slice_min(order_by = `delivery - pregnant padj`, n = 16)



dpsp_nulisa%>%
  filter(timepoint2!="post")%>%
  mutate(id_cat=as.character(id))%>%
  filter(targetName%in%top_16$targetName)%>%
  ggplot(., aes(x=factor(timepoint2), y=conc, fill=timepoint2))+
  geom_line(aes(group=id_cat), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetName, scales = "free")+
  theme_minimal()


combo_dpsp_results %>%
  mutate("label2" = if_else(targetName %in% top_16$targetName, targetName, NA))%>%
  ggplot(., aes(x=mean_preg_fc, y=-log10(`delivery - pregnant padj`), alpha=`delivery - pregnant padj`<0.05&abs(mean_preg_fc)>0.5, color=mean_preg_fc<0))+
  geom_point()+
  ggrepel::geom_text_repel(aes(label=label2, alpha=`delivery - pregnant padj`<0.05&abs(mean_preg_fc)>0.5),
                           size = 5,min.segment.length = 0,
                           position=ggpp::position_nudge_center(center_x = 0, x = 2, y=0.00000000001, center_y = 5))+
  ggtitle("n = 92, red means higher at delivery")+
  geom_hline(yintercept = -log10(0.05+10^-14), linetype="dashed")+
  geom_vline(xintercept = -0.5, linetype="dashed")+
  geom_vline(xintercept = 0.5, linetype="dashed")+
  scale_alpha_manual(values=c(0.5, 1))+
  scale_color_manual(values=c("darkred", "darkblue"))+
  scale_x_continuous(limits=c(-5, 5))+
  xlab("log2 fold change")+
  ylab("-log10 q value")+
  theme_minimal()+
  theme(legend.position="none")

id_columns <- c("id", "timepoint2", "plate", "sample")

wide_df2 <- dpsp_nulisa %>%
  filter(sample!="1094_pregnant_1")%>%
  filter(targetName %notin% c("IFNA2", "CTSS"))%>%
  # filter(sample_id %notin% c("X384_S.1_D1V93U", "X744_A.1_D1KT5Z", "X323_A14_D1DZ2D", "X496_NM7_D1CAYS", "X316_S.1_D12FNR", "X667_NM7_D1GBYA", "X 176 S_t14 D1FXRJ", "X164_NM0_D1WA6Q"))%>%
  pivot_wider(names_from = targetName, values_from = conc, id_cols = all_of(id_columns))


rownames(wide_df2) <- wide_df2$sample

# each row = measurement; each column = feature
big_pca <-  prcomp(wide_df2[,(length(id_columns)+1):ncol(wide_df2)], center = T)
pca_plot_data <- as.data.frame(cbind(wide_df2, big_pca$x))

loadings_df <- data.frame(big_pca$rotation)
loadings_df$targetName <- rownames(loadings_df)
loadings_df$targetName <- factor(loadings_df$targetName, levels = loadings_df$targetName[order(loadings_df$PC1)])

pc1_plot <- ggplot(loadings_df, aes(x=factor(targetName, levels = targetName[order(loadings_df$PC1)]), y=PC1, fill=targetName))+
  geom_bar(stat = "identity")+
  # scale_fill_manual(values = pc2_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")

plate_pca12 <- pca_plot_data %>%
  # mutate(qc_failed=as.factor(qc_failed))%>%
  ggplot(., aes(x=PC1, y=PC2, color=plate))+
  geom_point()+
  stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  # viridis::scale_color_viridis(discrete = T)+
  # viridis::scale_color_viridis(discrete = T)+
  theme(axis.text = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/DPSP/big_experiment/figures/plate_pca12.png", plate_pca12, height=4, width = 6, bg="white")



# neurocog linear regression



dpsp_neuro_results <- nulisa_plus_neuro %>%
  filter(timepoint2!="post")%>%
  mutate(id_cat=as.character(id))%>%
  pivot_longer(cols = ends_with("composite"), names_to = "composite_kind", values_to = "composite_score")%>%
  group_by(targetName, timepoint2, composite_kind)%>%
  nest() %>%
  mutate(model=map(data, ~lm(conc~composite_score, data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(p=map_dbl(summary, ~coef(.)[8])) %>%
  group_by(timepoint2)%>%
  mutate("padj" = p.adjust(p, method="fdr"))


dpsp_neuro_results <- nulisa_plus_neuro %>%
  filter(timepoint2!="post")%>%
  mutate(id_cat=as.character(id))%>%
  pivot_longer(cols = ends_with("composite"), names_to = "composite_kind", values_to = "composite_score")%>%
  group_by(targetName, timepoint2, composite_kind)%>%
  nest() %>%
  mutate(model=map(data, ~lm(conc~composite_score, data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(p=map_dbl(summary, ~coef(.)[8])) %>%
  group_by(timepoint2)%>%
  mutate("padj" = p.adjust(p, method="fdr"))



#did a bunch of simple logistic regressions
#no relationship with evermalaria
#no relationship with everfever
#no relationship with preterm
# p=0.1 for motor skills and anyhp, both timepoints
#no relationship with LBWdich
# p=0.16 pregnant, lang composite
#no relationship with gender

#significant interaction between language skills and maternal education

#simple lms
# no relationship with GA computed
bino_lms <- nulisa_plus_neuro %>%
  filter(timepoint2!="post")%>%
  mutate(id_cat=as.character(id))%>%
  mutate("multi"=if_else(gravidcat==3, 1, 0))%>%
  pivot_longer(cols = ends_with("composite"), names_to = "composite_kind", values_to = "composite_score")%>%
  pivot_longer(cols=c(evermalaria, everfever, preterm, anyHP, LBWdich, Olevel, SGA, multi), names_to = "covariate", values_to = "covarite_value")%>%
  distinct(id, composite_kind, composite_score, gender, covariate, covarite_value)%>%
  group_by(covariate, composite_kind)%>%
  nest() %>%
  mutate(model=map(data, ~glm(covarite_value~composite_score+gender, data=., family = "binomial"))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(score_p=map_dbl(summary, ~coef(.)[11]),
         gender_p=map_dbl(summary, ~coef(.)[12]))


#maternal age
linear_lms <- nulisa_plus_neuro %>%
  filter(timepoint2!="post")%>%
  mutate(id_cat=as.character(id))%>%
  pivot_longer(cols = ends_with("composite"), names_to = "composite_kind", values_to = "composite_score")%>%
  pivot_longer(cols=c(GAcomputed, totalfever, BMIenrol, AGE), names_to = "covariate", values_to = "covariate_value")%>%
  distinct(id, composite_kind, composite_score, gender, covariate, covariate_value)%>%
  group_by(covariate, composite_kind)%>%
  nest() %>%
  mutate(model=map(data, ~lm(covariate_value~composite_score+gender, data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(score_p=map_dbl(summary, ~coef(.)[11]),
         gender_p=map_dbl(summary, ~coef(.)[12]))


#nothing
nb_lms <- nulisa_plus_neuro %>%
  filter(timepoint2!="post")%>%
  mutate(id_cat=as.character(id))%>%
  pivot_longer(cols = ends_with("composite"), names_to = "composite_kind", values_to = "composite_score")%>%
  pivot_longer(cols=c(gravidcat, totalvisits, totalmalaria), names_to = "covariate", values_to = "covarite_value")%>%
  distinct(id, composite_kind, composite_score, gender, covariate, covarite_value)%>%
  group_by(covariate, composite_kind)%>%
  nest() %>%
  mutate(model=map(data, ~MASS::glm.nb(covarite_value~composite_score+gender, data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(score_p=map_dbl(summary, ~coef(.)[11]),
         gender_p=map_dbl(summary, ~coef(.)[12]))



nulisa_plus_neuro %>%
  filter(timepoint2!="post")%>%
  mutate(id_cat=as.character(id))%>%
  pivot_longer(cols = ends_with("composite"), names_to = "composite_kind", values_to = "composite_score")%>%
  distinct(id, composite_kind, composite_score, gender, AGE)%>%
  ggplot(., aes(x=AGE, y=composite_score, colour = factor(gender)))+
  geom_point()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  facet_wrap(~composite_kind)



nulisa_plus_neuro %>%
  filter(timepoint2!="post")%>%
  mutate(id_cat=as.character(id))%>%
  pivot_longer(cols = ends_with("composite"), names_to = "composite_kind", values_to = "composite_score")%>%
  distinct(id, composite_kind, composite_score, gender, Olevel)%>%
  ggplot(., aes(x=factor(Olevel), y=composite_score, color=factor(gender)))+
  geom_boxplot()+
  ggpubr::stat_compare_means()+
  facet_wrap(~composite_kind)



# maternal analytes and neurocog ####
linear_lms <- nulisa_plus_neuro %>%
  filter(timepoint2!="post")%>%
  mutate(id_cat=as.character(id))%>%
  pivot_longer(cols = ends_with("composite"), names_to = "composite_kind", values_to = "composite_score")%>%
  group_by(timepoint2, targetName, composite_kind)%>%
  nest() %>%
  mutate(model=map(data, ~lm(conc~composite_score+factor(Olevel), data=.))) %>%
  mutate(summary=map(model, ~summary(.)))%>%
  mutate(p=map_dbl(summary, ~coef(.)[8]))%>%
  group_by(composite_kind,timepoint2)%>%
  mutate(padj=p.adjust(p, method="fdr"))





