library(dplyr)
library(tidyr)
library(sva)

#--------------------------------------------------
# 0) optional but strongly recommended: log transform
#--------------------------------------------------
clean_data <- read.csv("~/Library/CloudStorage/Box-Box/Florian Bach's Externally Shareable Files/ELIA/raw_prosynk_micdrop_combined.csv")
clean_data2 <- clean_data %>%
  filter(!is.na(conc))   # small offset avoids log(0)

#--------------------------------------------------
# 1) PLATE NORMALIZATION (manufacturer option 2)
#    remove within-study plate drift
#--------------------------------------------------
plate_norm <- clean_data2 %>%
  # global per-protein median across ALL plates
  group_by(targetName) %>%
  mutate(global_median = median(conc, na.rm = TRUE)) %>%
  ungroup() %>%
  # plate-specific median
  group_by(targetName, plate) %>%
  mutate(plate_median = median(conc, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(conc_norm = conc - plate_median + global_median) %>%
  select(-plate_median, -global_median)

#--------------------------------------------------
# 2) CREATE SAMPLE METADATA TABLE
#--------------------------------------------------
meta <- plate_norm %>%
  mutate(sample_id=paste(id, timepoint, study, sep="_"))%>%
  filter(sample_id!="micdrop_11660_24 weeks_MICDROP")%>%
  select(id, study, treatmentarm, timepoint, sample_id) %>%
  distinct()

#--------------------------------------------------
# 3) LONG → WIDE (samples x proteins)
#--------------------------------------------------



wide <- plate_norm %>%
  mutate(sample_id=paste(id, timepoint, study, sep="_"))%>%
  filter(sample_id!="micdrop_11660_24 weeks_MICDROP")%>%
  select(sample_id, study, treatmentarm, timepoint, targetName, conc_norm) %>%
  pivot_wider(names_from = targetName, values_from = conc_norm)

#--------------------------------------------------
# 4) BUILD EXPRESSION MATRIX (proteins x samples)
#--------------------------------------------------
mat <- wide %>%
  tibble::column_to_rownames("sample_id") %>%
  select(-study, -treatmentarm, -timepoint) %>%  # drop metadata columns
  as.matrix() %>%                                 # explicitly cast to matrix
  t()

#--------------------------------------------------
# 5) FILTER LOW-VARIANCE PROTEINS (prevents ComBat NA/Inf)
#--------------------------------------------------
protein_sd <- apply(mat, 1, sd, na.rm = TRUE)
mat <- mat[protein_sd > 0.05, ]

#--------------------------------------------------
# 6) DESIGN MATRIX (protect biology!)
#--------------------------------------------------
meta <- meta %>% filter(sample_id %in% colnames(mat))

bio_model <- model.matrix(
  ~ treatmentarm + timepoint,
  data = meta
)

#--------------------------------------------------
# 7) COMBAT — STUDY / KIT LOT HARMONIZATION
#--------------------------------------------------
combat_mat <- ComBat(
  dat   = mat,
  batch = meta$study,
  mod   = bio_model,
  par.prior = TRUE
)

# combat_mat is now your harmonized expression matrix

long_combat_mat <- combat_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("targetName") %>%
  pivot_longer(
    cols = -targetName,
    names_to = "sample_id",
    values_to = "conc_combat"
  )%>%
  left_join(., meta, by="sample_id")
  
  
  



batchy_analytes <- clean_data%>%
  group_by(targetName, study)%>%
  filter(!is.na(conc))%>%
  summarise("mean_conc"=mean(conc))%>%
  pivot_wider(names_from = study, values_from = mean_conc)%>%
  mutate(diff=MICDROP-PROSYNK)%>%
  arrange(desc(abs(diff)))%>%
  pull(targetName)




exhibition_analytes <- c("IL10", "TLR3", "LAG3", "IFNW1")
exhibition_analytes1 <-  batchy_analytes[1:4]
exhibition_analytes2 <-  batchy_analytes[5:8]
exhibition_analytes3 <-  batchy_analytes[9:12]
exhibition_analytes4 <- batchy_analytes[13:16]
exhibition_analytes5 <- batchy_analytes[17:20]
exhibition_analytes6 <- batchy_analytes[21:24]
exhibition_analytes7 <- batchy_analytes[235:238]
exhibition_analytes8 <- batchy_analytes[239:242]
exhibition_analytes9 <- batchy_analytes[243:246]
exhibition_analytes10 <- c("GDF15")

list_of_exhibitions <- list(exhibition_analytes,
                            exhibition_analytes1,
                            exhibition_analytes2,
                            exhibition_analytes3,
                            exhibition_analytes4,
                            exhibition_analytes5,
                            exhibition_analytes6,
                            exhibition_analytes7,
                            exhibition_analytes8,
                            exhibition_analytes9,
                            exhibition_analytes10)

study_and_arm_cols <- c("darkred", "#233875", "darkmagenta", "orchid", "brown1", "blueviolet")
names(study_and_arm_cols) <- c("MICDROP.Placebo",   "MICDROP.DP 1 year", "PROSYNK.Labinic",   "PROSYNK.Probiotic", "PROSYNK.Placebo",   "PROSYNK.Lab4b")


corrected_long <- read.csv("~/Library/CloudStorage/Box-Box/Florian Bach's Externally Shareable Files/ELIA/limma_corrected_long_format.csv")



for(i in 1:length(list_of_exhibitions)){
pre_batch_correction <- clean_data %>%
  filter(timepoint %in% c("8 weeks", "24 weeks", "52 weeks"),
         #treatmentarm == "Placebo"
  )%>%
  filter(targetName %in% list_of_exhibitions[[i]])%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=conc, fill=interaction(study,treatmentarm)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName, scales="free_y",nrow=1)+
  scale_fill_manual(values=study_and_arm_cols)+
  theme_minimal()+
  ylab("NPQ")+
  theme(legend.title = element_blank(),
        axis.title.x = element_blank())

post_z_score_correction <- clean_data %>%
  filter(timepoint %in% c("8 weeks", "24 weeks", "52 weeks"),
         #treatmentarm == "Placebo"
  )%>%
  filter(targetName %in% list_of_exhibitions[[i]])%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=normalized_conc, fill=interaction(study,treatmentarm)))+
  geom_boxplot(outliers = F)+
  ylab("z score")+
  facet_wrap(~targetName, scales="free_y",nrow=1)+
  scale_fill_manual(values=study_and_arm_cols)+
  theme_minimal()+
  theme(legend.title = element_blank(),
        axis.title.x = element_blank())

post_combat_correction <- long_combat_mat %>%
  filter(timepoint %in% c("8 weeks", "24 weeks", "52 weeks"),
         #treatmentarm == "Placebo"
  )%>%
  filter(targetName %in% list_of_exhibitions[[i]])%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=conc_combat, fill=interaction(study,treatmentarm)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName, scales="free_y",nrow=1)+
  scale_fill_manual(values=study_and_arm_cols)+
  theme_minimal()+
  ylab("combat")+
  theme(legend.title = element_blank(),
        axis.title.x = element_blank())

post_limma_correction <- corrected_long %>%
  filter(timepoint %in% c("8 weeks", "24 weeks", "52 weeks"),
         #treatmentarm == "Placebo"
  )%>%
  filter(targetName %in% list_of_exhibitions[[i]])%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=corrected_npq, fill=interaction(study,treatmentarm)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName, scales="free_y",nrow=1)+
  scale_fill_manual(values=study_and_arm_cols)+
  theme_minimal()+
  ylab("limma")+
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank())

combo_plot <- pre_batch_correction /
  post_z_score_correction /
  post_combat_correction  +patchwork::plot_layout(axes = "collect")

ggsave(paste("~/postdoc/stanford/ELIA/nulisa_combo_figures/batch_integration_comparison_plot",i, ".png", sep=""), combo_plot, width=12, height=12, dpi=444)

} 





#normal
clean_data
#combat
long_combat_mat
#limma
corrected_long


big_merge <- long_combat_mat%>%
  select(id, timepoint, study, targetName, conc_combat)%>%
  right_join(., clean_data)
  
bigger_merge <- corrected_long%>%
  mutate(conc_limma=corrected_npq)%>%
  select(id, timepoint, study, targetName, conc_limma)%>%
  right_join(., big_merge)

write.csv(bigger_merge, "~/Library/CloudStorage/Box-Box/Florian Bach's Externally Shareable Files/ELIA/corrected_prosynk_micdrop_combined.csv", row.names = F)
