library(tidyr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(patchwork)
library(emmeans)

# CSP ELISA figures ####
norm_df_with_meta <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/csp_elisa_norm_df_with_meta.csv")
random_200 <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/Samples selected for Florian Oct 30 2024.dta")


(norm_csp_time_plot <- norm_df_with_meta%>%
  filter(old_id %notin% c("donor3", "pool"), timepoint=="52 weeks")%>%
  ggplot(., aes(x=treatmentarm, y=norm_OD, fill=treatmentarm))+
  geom_violin()+
  geom_boxplot(color="white", width=0.2, outliers = F, position = position_dodge(width=0.75))+
  geom_point(shape=21, size=1.5, color="white", position=position_jitter(width=0.02), stroke = 0.5)+
  ggpubr::stat_compare_means(
    aes(group = treatmentarm),
    label = "p.format",        # shows formatted p-values (e.g., 0.001, 0.05)
    method = "wilcox.test",    # or "t.test", depending on your data
    vjust=0.51, label.x = 1.35
  )+
  geom_hline(yintercept = 0.31)+
  # geom_violin(draw_quantiles = seq(0,1,0.25), aes(fill=treatmentarm))+
  ylab("normalized intensity (AU)")+
  theme_minimal(base_size = 15)+
  scale_fill_manual(values=treatment_palette)+
  theme(legend.position = "none",
        legend.title = element_blank(),
        #axis.text.y = element_blank(),
        axis.text = element_text(size=14),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank()))

ggsave("~/postdoc/stanford/manuscripts/MICDROP_clinial/norm_csp_time_plot.png", norm_csp_time_plot, width=5, height = 4, dpi=444, bg="white")
ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/ELISA/normalized_csp_52weeks.png", norm_csp_time_plot, width=5, height = 4, dpi=444, bg="white")
ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/ELISA/normalized_csp_52weeks.pdf", norm_csp_time_plot, width=5, height = 4, dpi=444, bg="white")


para12_plot <- norm_df_with_meta%>%
  filter(old_id %notin% c("donor3", "pool"))%>%
  filter(!is.na(total_n_para_12), timepoint=="52 weeks")%>%
  ggplot(., aes(x=factor(total_n_para_12), y=norm_OD, fill=factor(total_n_para_12)))+
  geom_boxplot(, outliers = F)+
  geom_point(shape=21, stroke=0.5, position=position_jitterdodge(jitter.width =  0.1, dodge.width = 0.75))+
  theme_minimal(base_size = 15)+
  ggpubr::stat_cor(method="spearman", mapping = aes(x=total_n_para_12, y=norm_OD, group=treatmentarm), inherit.aes = F)+
  facet_wrap(~treatmentarm)+
  xlab("\nparasitemic months in the first year of life")+
  ylab("normalized intensity (AU) at 52 weeks\n")+
  scale_fill_manual(values=viridis::magma(n=13))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(size=11),
        axis.title = element_text(size=14.5))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/figures/norm_para12_plot.png", para12_plot, width=6, height = 4, dpi=444, bg="white")
ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/ELISA/normalized_csp_parasite_incidence.png", para12_plot, width=5, height = 4, dpi=444, bg="white")
ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/ELISA/normalized_csp_parasite_incidence.pdf", para12_plot, width=5, height = 4, dpi=444, bg="white")



(parasitemia_plot <- norm_df_with_meta%>%
  filter(old_id %notin% c("donor3", "pool"), timepoint=="52 weeks", treatmentarm=="Placebo")%>%
  ggplot(., aes(x=10^log_qpcr, y=norm_OD, color=treatmentarm))+
  geom_smooth(method="lm")+
  geom_point(position = position_jitter(width=0.1))+
  # geom_hline(yintercept =0.31)+
  # geom_violin(draw_quantiles = seq(0,1,0.25), aes(fill=treatmentarm))+
  xlab("parasites / μL")+
  ylab("normalized intensity (AU)")+
  ggpubr::stat_cor(method="spearman", label.y = c(1.5, 1.5))+
  facet_wrap(~treatmentarm)+
  theme_minimal(base_size = 15)+
  scale_color_manual(values=treatment_palette)+
  scale_x_log10(labels=scales::label_log(), breaks = c(10^-3, 10^-1, 10^1, 10^3, 10^5))+
  scale_y_continuous(limits = c(-0.2, NA))+
  guides(color = guide_legend(override.aes = list(label = "")))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        #axis.text.y = element_blank(),
        axis.text = element_text(size=14),
        axis.title.y = element_text(size=16),
        strip.text = element_text(size=16)))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/figures/norm_parasitemia_plot.png", parasitemia_plot, width=6, height = 4, dpi=444, bg="white")
ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/ELISA/normalized_csp_parasitemia_correlation.png", parasitemia_plot, width=5, height = 4, dpi=444, bg="white")
ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/ELISA/normalized_csp_parasitemia_correlation.pdf", parasitemia_plot, width=5, height = 4, dpi=444, bg="white")


# PfEMP1 heatmap ####
long_luminex <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/long_luminex.csv")

# 1. create matrix for clustering
df <- long_luminex %>%
  filter(
    antigen != "tetanus",
    antigen != "var2csa",
    timepoint %in% c("8 weeks", "24 weeks", "52 weeks")
  ) %>%
  mutate(
    log_mfi = log10(MFI + 1),
    timepoint = factor(timepoint, levels = c("8 weeks", "24 weeks", "52 weeks"))
  )


mat <- df %>%
  arrange(treatmentarm)%>%
  filter(timepoint!= "8 weeks")%>%
  unite(sample_id, id, timepoint, treatmentarm, remove = FALSE) %>%
  mutate(log2_mfi=log2(MFI+1))%>%
  select(sample_id, antigen, log2_mfi) %>%
  pivot_wider(names_from = sample_id, values_from = log2_mfi) %>%
  tibble::column_to_rownames("antigen") %>%
  as.matrix()

meta <- df %>%
  filter(timepoint!= "8 weeks")%>%
  arrange(treatmentarm)%>%
  distinct(sample = paste(id, timepoint, treatmentarm, sep="_"),
           id, timepoint, treatmentarm)

library(grid)

ha <- HeatmapAnnotation(
  treatmentarm = meta$treatmentarm,
  col = list(
    treatmentarm = c(
      "Placebo" = "darkred",
      "DP 1 year" = "darkblue"
    )
  ),
  show_annotation_name = c(treatmentarm = FALSE),
  show_legend = c(treatmentarm = FALSE)
  )

mat_z <- t(scale(t(mat)))

q <- quantile(mat_z, c(0.02, 0.5, 0.98), na.rm = TRUE)


var_gene_heatmap <- Heatmap(
  mat_z,
  name = "log2 MFI\nz score",
  top_annotation = ha,
  col = circlize::colorRamp2(q, c("#2C3E8F", "white", "#E66101")),
  # breaks=c(min(mat, na.rm = TRUE), 9.7, max(mat, na.rm = TRUE)),
  # colors=c("#2C3E8F", "#F7F7F7", "#E66101")),
# ),
#   col = circlize::colorRamp2(
#     c(min(mat, na.rm = TRUE), 2.89, max(mat, na.rm = TRUE)),
#     c("#4C72B0", "#FAFAFA", "#DD8452")
#   ),
  column_split = meta$timepoint,
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  cluster_column_slices = TRUE,   # independent clustering
  cluster_rows = TRUE,
  row_names_side = "left",
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  rect_gp = gpar(col = "white", lwd = 0.3)
)


arm <- meta$treatmentarm

add_arm_labels <- function(y = 0.62) {
  n <- length(meta$treatmentarm)
  
  r <- rle(as.character(meta$treatmentarm))
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  mids <- (starts + ends) / 2
  
  for (i in seq_along(r$values)) {
    grid.text(
      label = r$values[i],
      x = unit(mids[i] / n, "npc"),
      y = unit(y, "npc"),
      gp = gpar(col = "white", fontsize = 10, fontface = "bold")
    )
  }
}

pdf("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/var_genes/var_gene_heatmap.pdf", width = 10, height=6)
draw(var_gene_heatmap)
decorate_annotation("treatmentarm", slice = 1, {
  add_arm_labels(y = 0.55)
})
# add to slice 2 (the other column slice)
decorate_annotation("treatmentarm", slice = 2, {
  add_arm_labels(y = 0.55)
})
dev.off()

png("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/var_genes/var_gene_heatmap.png", width = 10, height=8, units = "in", res = 444)
draw(var_gene_heatmap)
dev.off()



# 
# (p1 <- long_luminex %>%
#     filter(antigen!="tetanus", antigen!="var2csa", timepoint %in% c("8 weeks", "24 weeks", "52 weeks")) %>%
#     mutate(log_mfi = log10(MFI + 1), timepoint=factor(timepoint, levels = c("8 weeks", "24 weeks", "52 weeks"))) %>%
#     mutate(sample_id = paste(id, timepoint, sep="_"))%>%
#     ggplot(aes(x = factor(sample_id, levels=id_order), y = antigen, fill = log_mfi)) +
#     geom_tile() +
#     # scale_fill_gradientn(colors=viridis::mako(n=5))+
#     scale_fill_gradient2(
#       low = "#233875",
#       mid = "white",
#       high = "#FF5A00",
#       midpoint = 2.89
#     ) +
#     facet_wrap(~timepoint)+
#     theme_minimal()+
#     theme(axis.text.x = element_blank())
# )
# 

# full length vs. NANP ####

full_length <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/csp_elisa_norm_df_with_meta.csv")
nanp <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/nanp_repeat_elisa_alex.csv")

nanp_merge <- nanp%>%
  mutate(id=as.numeric(ID),
         CSP_NANP=Average)%>%
  filter(Dup_remove=="no")

full_length_merge <- full_length%>%
  filter(timepoint=="52 weeks")%>%
  mutate(CSP_full=norm_OD)

combo <- full_length_merge%>%
  inner_join(., nanp_merge, by="id")


combo %>%
  ggplot(., aes(x=CSP_NANP, y=CSP_full))+
  geom_point(aes( color=treatmentarm))+
  theme_minimal()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method = "spearman", label.y = 2)+
  scale_x_continuous(trans="log2", expand=expansion(mult=0.15))+
  scale_y_continuous(trans="log2")+
  theme(axis.text = element_blank())+
  scale_color_manual(values=rev(c("red", "blue")))

full_length_plot <- full_length%>%
  filter(old_id %notin% c("donor3", "pool"), timepoint=="52 weeks", id %in% random_200$id)%>%
  ggplot(., aes(x=treatmentarm, y=norm_OD, fill=treatmentarm))+
  geom_violin()+
  # geom_boxplot(color="white", width=0.2, outliers = F, position = position_dodge(width=0.75))+
  geom_point(shape=21, size=1.5, color="white", position=position_jitter(width=0.04), stroke = 0.5)+
  ggpubr::stat_compare_means(
    aes(group = treatmentarm),
    label = "p.format",        # shows formatted p-values (e.g., 0.001, 0.05)
    method = "wilcox.test",    # or "t.test", depending on your data
    vjust=0.51, label.x = 1.35
  )+
  ggtitle("full length CSP")+
  ylab("normalized intensity (AU)")+
  # scale_y_continuous(trans="log2")+
  theme_minimal(base_size = 13)+
  scale_fill_manual(values=rev(c("darkred", "darkblue")))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        # axis.text.y = element_blank(),
        axis.text = element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.title.x = element_blank())

nanplot <- nanp_merge%>%
  left_join(full_length_merge, by="id")%>%
  filter(!is.na(treatmentarm), id %in% random_200$id)%>%
  ggplot(., aes(x=treatmentarm, y=Average, fill=treatmentarm))+
  geom_violin()+
  # geom_boxplot(color="white", width=0.2, outliers = F, position = position_dodge(width=0.75))+
  geom_point(shape=21, size=1.5, color="white", position=position_jitter(width=0.04), stroke = 0.5)+
  ggpubr::stat_compare_means(
    aes(group = treatmentarm),
    label = "p.format",        # shows formatted p-values (e.g., 0.001, 0.05)
    method = "wilcox.test",    # or "t.test", depending on your data
    vjust=0.51, label.x = 1.35
  )+
  # scale_y_continuous(trans="log2")+
  ggtitle("NANP repeats")+
  ylab("normalized intensity (AU)")+
  theme_minimal(base_size = 13)+
  scale_fill_manual(values=rev(c("darkred", "darkblue")))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        # axis.text.y = element_blank(),
        axis.text = element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.title.x = element_blank())

combo_plot <- full_length_plot + nanplot + plot_annotation(tag_levels = 'A')

ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/Pf_antibodies/CSP/nanp_and_csp_by_treatmentarm.png", combo_plot, width=6, height = 3, dpi=444, bg="white")
ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/Pf_antibodies/CSP/nanp_and_csp_by_treatmentarm.pdf", combo_plot, width=6, height = 3, dpi=444, bg="white")

### log2 version ####
log_full_length_plot <- full_length%>%
  filter(old_id %notin% c("donor3", "pool"), timepoint=="52 weeks")%>%
  ggplot(., aes(x=treatmentarm, y=norm_OD, fill=treatmentarm))+
  geom_violin()+
  geom_boxplot(color="white", width=0.2, outliers = F, position = position_dodge(width=0.75))+
  geom_point(shape=21, size=1.5, color="white", position=position_jitter(width=0.02), stroke = 0.5)+
  ggpubr::stat_compare_means(
    aes(group = treatmentarm),
    label = "p.format",        # shows formatted p-values (e.g., 0.001, 0.05)
    method = "wilcox.test",    # or "t.test", depending on your data
    vjust=0.51, label.x = 1.35, label.y = 40
  )+
  ggtitle("full length CSP")+
  ylab("log2 normalized intensity (AU)")+
  scale_y_continuous(trans="log2")+
  scale_fill_manual(values=rev(c("darkred", "darkblue")))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.title.x = element_blank())

log_nanplot <- nanp_merge%>%
  left_join(full_length_merge, by="id")%>%
  filter(!is.na(treatmentarm))%>%
  ggplot(., aes(x=treatmentarm, y=Average, fill=treatmentarm))+
  geom_violin()+
  geom_boxplot(color="white", width=0.2, outliers = F, position = position_dodge(width=0.75))+
  geom_point(shape=21, size=1.5, color="white", position=position_jitter(width=0.02), stroke = 0.5)+
  ggpubr::stat_compare_means(
    aes(group = treatmentarm),
    label = "p.format",        # shows formatted p-values (e.g., 0.001, 0.05)
    method = "wilcox.test",    # or "t.test", depending on your data
    vjust=0.51, label.x = 1.35
  )+
  scale_y_continuous(trans="log2")+
  ggtitle("NANP repeats")+
  ylab("log2 normalized intensity (AU)")+
  scale_fill_manual(values=rev(c("darkred", "darkblue")))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.title.x = element_blank())

combo_plot2 <- log_full_length_plot + log_nanplot + plot_annotation(tag_levels = 'A')

(combo_plot3 <-  full_length_plot + log_nanplot + plot_annotation(tag_levels = 'A'))
ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/ELISA/nanp_and_csp.png", combo_plot3, width=6, height = 3, dpi=444, bg="white")
ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/ELISA/nanp_and_csp_.pdf", combo_plot3, width=6, height = 3, dpi=444, bg="white")

# PfEMP1 cox regression ####

seroprev_est_df <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/mixture_model_df.csv")
colnames(seroprev_est_df) <- c("antigen", "estimated_seroprevalence", "cutoff")

short_luminex <- long_luminex%>%
  group_by(antigen)%>%
  mutate(seroprev_cutoff=if_else(log_mfi>seroprev_est_df$cutoff[seroprev_est_df$antigen==antigen], 1, 0))%>%
  filter(timepoint=="52 weeks")

# NULISA volcano ####

clean_data <- clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks", "68 weeks")))%>%
  filter(targetName %notin% c("CTSS", "LTA|LTB", "IFNA2"))

random_200 <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/Samples selected for Florian Oct 30 2024.dta")

mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")

clean_data <- clean_data%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))

treatment_purf <- clean_data%>%
  filter(id %in% random_200$id)%>%
  filter(mstatus==0, treatmentarm!="DP 2 years")%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(conc~timepoint*treatmentarm+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ treatmentarm | timepoint)))%>%
  mutate(emm2=map(time_model, ~emmeans(., specs = pairwise ~ timepoint | treatmentarm)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  mutate("8 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("24 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("52 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  pivot_longer(cols=ends_with("weeks"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

sigs <- treatment_purf%>%
  filter(padj<0.05)%>%
  filter(contrast=="52 weeks")

treatment_purf%>%
  select(targetName, contrast, p, padj)%>%
  write.csv(., "~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/nulisa_treatment_regression.csv")



treatment_nulisa_fc_data <- clean_data %>%
  filter(id %in% random_200$id)%>%
  filter(treatmentarm != "DP 2 years", timepoint !=c("68 weeks"))%>%
  group_by(targetName, timepoint, treatmentarm)%>%
  summarise("mean_conc"=mean(conc, na.rm = T))%>%
  pivot_wider(names_from = treatmentarm, values_from = mean_conc)%>%
  mutate(fc=Placebo-`DP 1 year`)


treatment_nulisa_purf_fc <- treatment_purf%>%
  mutate("timepoint"=contrast)%>%
  left_join(., treatment_nulisa_fc_data, by=c("targetName", "timepoint"))%>%
  mutate(contrast=factor(contrast, levels=c("8 weeks", "24 weeks", "52 weeks")))

treatment_nulisa_sigs <- treatment_nulisa_purf_fc %>%
  filter(padj<0.05)




## volcano ####
(treatment_volcano <- treatment_nulisa_purf_fc %>%
   mutate("label2" = if_else(targetName %in% c(treatment_nulisa_sigs$targetName) & timepoint != "8 weeks", targetName, NA))%>%
   # mutate("label2" = if_else(.$targetName, .$contrast)== cbind(treatment_nulisa_sigs$targetName, treatment_nulisa_sigs$contrast)), targetName, NA))%>%
   ggplot(., aes(x=-fc, y=padj, alpha=padj<0.05&abs(fc)>=0.5, color=fc<0))+
   geom_hline(yintercept = 0.05, linetype="dashed", alpha = 0.5)+
   geom_vline(xintercept = -0.5, linetype="dashed", alpha = 0.5)+
   geom_vline(xintercept = 0.5, linetype="dashed", alpha = 0.5)+
   
   geom_point()+
   ggrepel::geom_text_repel(aes(label=label2, alpha=padj<0.05&abs(fc)>=0.5),force = 3,
                            size = 4, min.segment.length = 0,
                            position=ggpp::position_nudge_center(center_x = 0, x = 0.8,
                                                                 center_y = 0.05, y=2.5)
   )+
   scale_alpha_manual(values=c(0.5, 1))+
   scale_color_manual(values=rev(c("#DC2626", "#2563EB")))+
   scale_x_continuous(limits=c(-2, 2))+
   scale_y_continuous(trans=c("log10", "reverse"))+
   labs(
     title = "**Differentially Abundant Plasma Proteins**  
    <span style='font-size:11pt'>
    <span style='color:#FF8C00;'>up in DP group </span> vs.
    <span style='color:darkblue;'>up in placebo group </span>",
   ) +
   # ggtitle("differentially abundant plasma proteins in placebo vs. DP")+
   facet_wrap(~contrast)+
   xlab("log2 fold change")+
   ylab("padj")+
   theme_minimal()+
   theme(legend.position="none",
         plot.title=ggtext::element_markdown()))

ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/NULISA/sig_nulisa_treatment_volcano_plot.png", treatment_volcano, height=4, width=10, dpi=400, bg="white")
ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/NULISA/sig_nulisa_treatment_volcano_plot.pdf", treatment_volcano, height=4, width=10, dpi=400, bg="white")


