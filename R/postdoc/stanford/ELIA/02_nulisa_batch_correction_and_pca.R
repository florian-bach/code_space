library(tidyverse)
library(ggplot2)
library(limma)

clean_data <- read.csv("~/Library/CloudStorage/Box-Box/Florian Bach's Externally Shareable Files/ELIA/corrected_prosynk_micdrop_combined.csv")

# some problematic samples:
outlier_samples <- scores_df%>%
  filter(PC1>=25)%>%
  select(id, study, treatmentarm, timepoint)%>%
  pull(id)
# id   study treatmentarm timepoint
# 1  micdrop_10951 MICDROP    DP 1 year  52 weeks
# 2  micdrop_11030 MICDROP      Placebo  52 weeks
# 3  micdrop_11110 MICDROP      Placebo  52 weeks
# 4  micdrop_11639 MICDROP    DP 1 year  24 weeks
# 5      PSK-10222 PROSYNK        Lab4b   8 weeks
# 6      PSK-10236 PROSYNK    Probiotic   8 weeks
# 7      PSK-10269 PROSYNK      Placebo   8 weeks
# 8      PSK-10280 PROSYNK      Labinic   8 weeks
# 9      PSK-10300 PROSYNK        Lab4b  24 weeks
# 10     PSK-10407 PROSYNK        Lab4b  24 weeks
# 11     PSK-10492 PROSYNK      Labinic  24 weeks
# PCA ####
#--- prepare data ---
pca_input <- clean_data %>%
  # filter(!id %in% c("micdrop_11660", "micdrop_10951", "PSK-10222", "PSK-10269", "PSK-10280", "PSK-10407"))%>%
  filter(!id %in% c("micdrop_11660", outlier_samples))%>%
  select(study, id, targetName, treatmentarm, timepoint, conc, gender_categorical) %>%
  drop_na()%>%
  pivot_wider(names_from = targetName, values_from = conc)

# needs features as rows, samples as columns
mat <- pca_input %>%
  select(where(is.numeric)) %>%
  t()

# protect biological variables so they aren't removed along with batch
bio_model <- model.matrix(~ timepoint + treatmentarm, data = pca_input)

mat_corrected <- removeBatchEffect(
  mat,
  batch = pca_input$study,
  design = bio_model
)

# check how many features have NA/Inf after correction
sum(is.infinite(mat_corrected))
sum(is.na(mat_corrected))


# are the NAs concentrated in specific features or spread across samples?
na_per_feature <- rowSums(is.na(mat_corrected))
na_per_sample <- colSums(is.na(mat_corrected))

summary(na_per_feature[na_per_feature > 0])
summary(na_per_sample[na_per_sample > 0])

bad_features <- names(na_per_feature[na_per_feature > 4])
bad_features  # check which proteins you're dropping

mat_clean <- mat_corrected[!rownames(mat_corrected) %in% bad_features, ]

sum(is.na(mat_clean))  # should be 0 or close to it

# impute remaining NAs with feature mean (sufficient for PCA)
mat_clean <- apply(mat_clean, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}) %>% t()

sum(is.na(mat_clean))  # should be 0 now

pca_corrected <- prcomp(t(mat_clean), scale. = TRUE)


# --- scores dataframe ---
scores_df <- as.data.frame(pca_corrected$x) %>%
  bind_cols(pca_input %>% select(id, study, treatmentarm, timepoint, gender_categorical))
# --- PCA scatter plot ---

ggplot(scores_df, aes(x = PC1, y = PC2, colour = study, shape = treatmentarm)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained["PC1"]),
    y = sprintf("PC2 (%.1f%%)", var_explained["PC2"]),
    title = "PCA scores"
  ) +
  geom_density2d()+
  facet_wrap(~factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))+
  theme_minimal()

ggplot(scores_df, aes(x = PC1, y = PC2, colour = gender_categorical, shape = treatmentarm)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained["PC1"]),
    y = sprintf("PC2 (%.1f%%)", var_explained["PC2"]),
    title = "PCA scores"
  ) +
  geom_density2d()+
  facet_wrap(~factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))+
  theme_minimal()

# --- loadings: top contributors to PC1 and PC2 ---
loadings_df <- as.data.frame(pca_corrected$rotation) %>%
  tibble::rownames_to_column("feature") %>%
  select(feature, PC1, PC2, PC3, PC4) %>%
  pivot_longer(-feature, names_to = "PC", values_to = "loading")

# bar plot of top N loadings per PC
top_n <- 20

library(dplyr)
library(ggplot2)

plot_df <- loadings_df %>%
  group_by(PC) %>%
  slice_max(abs(loading), n = top_n) %>%
  ungroup() %>%
  # create per-facet ordering variable
  arrange(PC, desc(abs(loading))) %>%
  group_by(PC) %>%
  mutate(feature_order = factor(
    paste(PC, feature, sep = "___"),
    levels = paste(PC, feature, sep = "___")
  )) %>%
  ungroup()

ggplot(plot_df, aes(x = feature_order, y = loading, fill = loading > 0)) +
  geom_col(show.legend = FALSE) +
  scale_x_discrete(labels = function(x) sub(".*___", "", x)) +  # remove PC prefix
  scale_fill_manual(values = c("tomato", "steelblue")) +
  facet_wrap(~ PC, scales = "free_y") +
  coord_flip() +
  labs(x = NULL, y = "Loading", title = paste("Top", top_n, "loadings per PC")) +
  theme_minimal()


# back_to_long_format linear regression ####

corrected_long <- t(mat_clean) %>%
  as.data.frame() %>%
  bind_cols(pca_input %>% select(id, study, treatmentarm, timepoint, gender_categorical)) %>%
  pivot_longer(
    cols = -c(id, study, treatmentarm, timepoint, gender_categorical),
    names_to = "targetName",
    values_to = "corrected_npq"
  )

write.csv(corrected_long, "~/Library/CloudStorage/Box-Box/Florian Bach's Externally Shareable Files/ELIA/limma_corrected_long_format.csv", row.names = F)

sutdy_and_arm_cols <- c("darkred", "#233875", "darkmagenta", "orchid", "brown1", "blueviolet")
names(sutdy_and_arm_cols) <- c("MICDROP.Placebo",   "MICDROP.DP 1 year", "PROSYNK.Labinic",   "PROSYNK.Probiotic", "PROSYNK.Placebo",   "PROSYNK.Lab4b")

exhibition_analytes <- c("IL10", "TLR3", "LAG3", "IFNW1")

pre_batch_correction <- clean_data %>%
  filter(timepoint %in% c("8 weeks", "24 weeks", "52 weeks"),
         #treatmentarm == "Placebo"
  )%>%
  filter(targetName %in% exhibition_analytes)%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=conc, fill=interaction(study,treatmentarm)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName, scales="free_y",nrow=1)+
  scale_fill_manual(values=sutdy_and_arm_cols)+
  theme_minimal()+
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank())

post_batch_correction <- corrected_long %>%
  filter(timepoint %in% c("8 weeks", "24 weeks", "52 weeks"),
         #treatmentarm == "Placebo"
  )%>%
  filter(targetName %in% exhibition_analytes)%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=corrected_npq, fill=interaction(study,treatmentarm)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName, scales="free_y",nrow=1)+
  scale_fill_manual(values=sutdy_and_arm_cols)+
  theme_minimal()+
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank())

pre_batch_correction /
  post_batch_correction
