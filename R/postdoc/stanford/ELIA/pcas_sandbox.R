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


# check how many features have NA/Inf after correction
sum(is.infinite(mat))
sum(is.na(mat))


# are the NAs concentrated in specific features or spread across samples?
na_per_feature <- rowSums(is.na(mat))
na_per_sample <- colSums(is.na(mat))

summary(na_per_feature[na_per_feature > 0])
summary(na_per_sample[na_per_sample > 0])

bad_features <- names(na_per_feature[na_per_feature > 4])
bad_features  # check which proteins you're dropping

mat_clean <- mat[!rownames(mat) %in% bad_features, ]

sum(is.na(mat_clean))  # should be 0 or close to it

# impute remaining NAs with feature mean (sufficient for PCA)
mat_clean <- apply(mat_clean, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}) %>% t()

sum(is.na(mat_clean))  # should be 0 now

pca_corrected <- prcomp(t(mat_clean), scale. = TRUE)

# Calculate variance explained
var_explained <- summary(pca_corrected)$importance["Proportion of Variance", ] * 100
# Or equivalently:
var_explained <- (pca_corrected$sdev^2 / sum(pca_corrected$sdev^2)) * 100
names(var_explained) <- colnames(pca_corrected$x)  # ensures PC1, PC2, etc. as names

# --- scores dataframe ---
scores_df <- as.data.frame(pca_corrected$x) %>%
  bind_cols(pca_input %>% select(id, study, treatmentarm, timepoint, gender_categorical))
# --- PCA scatter plot ---

raw_pca <- ggplot(scores_df, aes(x = PC1, y = PC2, colour = study, shape = treatmentarm)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained["PC1"]),
    y = sprintf("PC2 (%.1f%%)", var_explained["PC2"]),
    title = "raw PCA scores"
  ) +
  geom_density2d()+
  facet_wrap(~factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))+
  theme_minimal(base_size = 12)

ggsave("~/postdoc/stanford/ELIA/nulisa_combo_figures/raw_pca.png", width=8, height = 5, dpi=444)




combat_pca_input <- clean_data %>%
  # filter(!id %in% c("micdrop_11660", "micdrop_10951", "PSK-10222", "PSK-10269", "PSK-10280", "PSK-10407"))%>%
  filter(!id %in% c("micdrop_11660", outlier_samples))%>%
  select(study, id, targetName, treatmentarm, timepoint, conc_combat, gender_categorical) %>%
  drop_na()%>%
  pivot_wider(names_from = targetName, values_from = conc_combat)

# needs features as rows, samples as columns
combat_mat <- combat_pca_input %>%
  select(where(is.numeric)) %>%
  t()


# check how many features have NA/Inf after correction
sum(is.infinite(combat_mat))
sum(is.na(combat_mat))


# are the NAs concentrated in specific features or spread across samples?
na_per_feature <- rowSums(is.na(combat_mat))
na_per_sample <- colSums(is.na(combat_mat))

summary(na_per_feature[na_per_feature > 0])
summary(na_per_sample[na_per_sample > 0])

bad_features <- names(na_per_feature[na_per_feature > 4])
bad_features  # check which proteins you're dropping

combat_mat_clean <- combat_mat[!rownames(combat_mat) %in% bad_features, ]

sum(is.na(combat_mat_clean))  # should be 0 or close to it

# impute remaining NAs with feature mean (sufficient for PCA)
combat_mat_clean <- apply(combat_mat_clean, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}) %>% t()

sum(is.na(combat_mat_clean))  # should be 0 now

pca_corrected <- prcomp(t(combat_mat_clean), scale. = TRUE)

# Calculate variance explained
var_explained <- summary(pca_corrected)$importance["Proportion of Variance", ] * 100
# Or equivalently:
var_explained <- (pca_corrected$sdev^2 / sum(pca_corrected$sdev^2)) * 100
names(var_explained) <- colnames(pca_corrected$x)  # ensures PC1, PC2, etc. as names

# --- scores dataframe ---
scores_df <- as.data.frame(pca_corrected$x) %>%
  bind_cols(pca_input %>% select(id, study, treatmentarm, timepoint, gender_categorical))
# --- PCA scatter plot ---

combat_pca <- ggplot(scores_df, aes(x = PC1, y = PC2, colour = study, shape = treatmentarm)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained["PC1"]),
    y = sprintf("PC2 (%.1f%%)", var_explained["PC2"]),
    title = "combat PCA scores"
  ) +
  geom_density2d()+
  facet_wrap(~factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))+
  theme_minimal(base_size = 12)

ggsave("~/postdoc/stanford/ELIA/nulisa_combo_figures/combat_pca.png", combat_pca, width=8, height = 5, dpi=444)



treat_combat_pca <- ggplot(scores_df, aes(x = PC1, y = PC2, colour = treatmentarm, shape = treatmentarm)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained["PC1"]),
    y = sprintf("PC2 (%.1f%%)", var_explained["PC2"]),
    title = "combat PCA scores"
  ) +
  geom_density2d()+
  facet_wrap(~factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))+
  theme_minimal(base_size = 12)

ggsave("~/postdoc/stanford/ELIA/nulisa_combo_figures/treat_combat_pca.png", treat_combat_pca, width=8, height = 5, dpi=444)
