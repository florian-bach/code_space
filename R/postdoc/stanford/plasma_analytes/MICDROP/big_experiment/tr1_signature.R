library(dplyr)
library(tidyr)
library(glmnet)
library(ggplot2)

`%notin%`<- Negate(`%in%`)
# WGCNA try ####

rnaseq_networks <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/wgcna_RNAseq_MEs.csv")

nulisa_networks <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/wgcna_nulisa_MEs.csv")

t_cell_freqs <- read.csv("~/Downloads/df_jason_analysis (1).csv")

rnaseq_networks <- rnaseq_networks%>%
  mutate(timepoint_imm=case_when(timepoint_category=="Baseline" ~ -1,
                                 timepoint_category=="Day 0" ~ 0,
                                 timepoint_category=="Day 14" ~ 14,
                                 timepoint_category=="Day 28" ~ 28,
                                 timepoint_category=="Day 7" ~ 7,
                                 timepoint_category=="Older baseline" ~ -2,
                                 timepoint_category=="Day 84"~ -1),
         inf_type=case_when(infection=="A2"~"A",
                            infection=="asymptomatic"~"A",
                            infection=="symptomatic"~"S", .default = infection))

rnaseq_networks$sample_name <- paste(rnaseq_networks$cohortid,
                                     rnaseq_networks$inf_type,
                                     rnaseq_networks$timepoint_imm, sep="_")

t_cell_freqs$sample_name <- paste(t_cell_freqs$cohortid, substr(t_cell_freqs$infectiontype, 1, 1), t_cell_freqs$timepoint, sep="_")


nulisa_networks$crappy_names <- substr(nulisa_networks$sample_id, 5, 9)
nulisa_networks <-nulisa_networks%>%
  mutate(crappy_names2=case_when(nulisa_networks$crappy_names=="bad_b" ~ -2,
                                 nulisa_networks$crappy_names=="basel" ~ -1,
                                 nulisa_networks$crappy_names=="day0 " ~ 0,
                                 nulisa_networks$crappy_names=="day14" ~ 14,
                                 nulisa_networks$crappy_names=="day28"  ~ 28))

nulisa_networks$sample_name <- paste(substr(nulisa_networks$sample_id, 1, 3),
                                     substr(nulisa_networks$sample_id, nchar(nulisa_networks$sample_id), nchar(nulisa_networks$sample_id)),
                                     nulisa_networks$crappy_names2,
                                     sep = "_")


nulisa_networks_for_merge <- nulisa_networks%>%
  select(sample_name, MEthistle3, MEtomato3, MEturquoise3, MEwhite)
colnames(nulisa_networks_for_merge) <- c("sample_name", "protein_MEthistle3", "protein_MEtomato3", "protein_MEturquoise3", 'protein_MEwhite')


rna_networks_for_merge <- rnaseq_networks%>%
  select(sample_name, starts_with("ME", ignore.case = F))
colnames(rna_networks_for_merge)[2:15] <- paste("rna", colnames(rna_networks_for_merge)[2:15], sep="_")

combo_df <- t_cell_freqs%>%
  left_join(., rna_networks_for_merge, by="sample_name")%>%
  left_join(., nulisa_networks_for_merge, by="sample_name")

tr1_combo_df <- combo_df%>%
  pivot_longer(cols = ends_with("Frequency"), names_to = "gate", values_to = "freq")%>%
  pivot_longer(cols = contains("_ME", ignore.case = F), names_to = "network", values_to = "network_value")%>%
  filter(gate=="Tr1_Frequency", stim=="unstim")

tr1_combo_corr <- tr1_combo_df%>%
  group_by(network)%>%
  nest()%>%
  mutate(spearman_cor=map(data, ~cor.test(.$network_value, .$freq, method="spearman")),
         spearman_rho=map_dbl(spearman_cor, ~.$estimate),
         spearman_p=map_dbl(spearman_cor, ~.$p.value))

tr1_combo_df%>%
  # filter(network %in% c(rna_MEsalmon, protein_MEtomato3))%>%
  ggplot(., aes(x=network_value, y=freq))+
  geom_smooth(method="lm")+
  geom_point()+
  ggpubr::stat_cor(method = "spearman", size=4, color="red")+
  facet_wrap(~network)+
  theme_minimal()


# table(t_cell_freqs$sample_name %in% nulisa_networks$sample_name)
# table(t_cell_freqs$sample_name %in% rnaseq_networks$sample_name)
# 
# table(nulisa_networks$sample_name %in% t_cell_freqs$sample_name)
# table(nulisa_networks$sample_name %in% rnaseq_networks$sample_name)
# 
# table(rnaseq_networks$sample_name %in% nulisa_networks$sample_name)
# table(rnaseq_networks$sample_name %in% t_cell_freqs$sample_name)

# lasso ####
library(glmnet)


# Install if needed
# install.packages("glmnet")

nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/tr1_paper/revised_baseline_clean_musical_combo_with_metadata.csv")
nulisa_data <- nulisa_data%>%
  mutate(infectiontype=substr(infectiontype, 1, 1))

wide_nulisa <- nulisa_data%>%
  filter(targetName %notin% c("CTSS", "LTA|LTB", "IFNA2"))%>%
  pivot_wider(names_from = targetName, names_prefix = "p_", values_from = concentration, id_cols = sample_id)

slim_nulisa_data <- nulisa_data %>%
  mutate(id=as.character(id))%>%
  select(id, date, timepoint, timepoint_imm, infectiontype, targetName, concentration, log_qpcr, ageyrs, gender_categorical)%>%
  filter(timepoint!="bad_baseline")%>%
  mutate(timepoint_imm=if_else(timepoint=="baseline", -1, timepoint_imm))

slim_cell_count_data <- read.csv("~/Downloads/df_jason_analysis (1).csv")

slim_cell_count_data <- slim_cell_count_data%>%
  pivot_longer(cols = ends_with("Frequency"), names_to = "gate", values_to = "freq")%>%
  mutate(id=cohortid, timepoint_imm=timepoint)%>%
  mutate(infectiontype=substr(infectiontype, 1, 1))%>%
  mutate("timepoint"=case_when(
    timepoint_imm==-1~"baseline",
    timepoint_imm==0~"day0",
    timepoint_imm==7~"day7",
    timepoint_imm==14~"day14",
    timepoint_imm==28~"day28"))%>%
  select(id, date, timepoint, timepoint_imm, infectiontype, gate, stim, freq)


long_combo <- slim_nulisa_data%>%
  mutate(id=as.integer(id))%>%
  inner_join(., slim_cell_count_data, by = c("id", "infectiontype", "timepoint_imm"))%>%
  mutate(timepoint=timepoint.x)%>%
  select(-timepoint.x, -timepoint.y)%>%
  filter(!is.na(freq))

# -----------------------------
# Prepare outcome vector 
# -----------------------------

y_df <- long_combo %>%
  filter(gate=="Tr1_Frequency", stim=="iRBCs", !is.na(timepoint)) %>%
  mutate(sample_name = paste(id, timepoint, substr(infectiontype,1,1), sep=" ")) %>%
  filter(sample_name %in% wide_nulisa$sample_id) %>%
  distinct(sample_name, id, freq) %>%
  mutate(z_freq=scale(freq, center = TRUE, scale = TRUE))%>%
  arrange(sample_name)

y <- y_df$freq
sample_ids <- y_df$sample_name
subject_ids <- y_df$id

# -----------------------------
# Predictor matrix
# -----------------------------

X <- wide_nulisa %>%
  filter(sample_id %in% sample_ids) %>%
  arrange(sample_id) %>%
  tibble::column_to_rownames("sample_id") %>%
  select(starts_with("p_")) %>%
  as.matrix()

# -----------------------------
# Scale predictors explicitly
# -----------------------------

X_scaled <- scale(X)

train_means <- attr(X_scaled,"scaled:center")
train_sds   <- attr(X_scaled,"scaled:scale")

# -----------------------------
# Create LOIO folds
# -----------------------------

# each individual becomes its own fold
foldid <- as.numeric(factor(subject_ids))

# -----------------------------
# Fit elastic net with grouped CV
# -----------------------------

set.seed(123)

cvfit <- cv.glmnet(
  x = X_scaled,
  y = y,
  alpha = 0.5,           # elastic net
  family = "gaussian",
  foldid = foldid,       # <- key change
  keep = TRUE
)

plot(cvfit)

# -----------------------------
# Extract model
# -----------------------------
lambda_use <- cvfit$lambda.min
lambda_1se <- cvfit$lambda.1se

coef_lasso <- coef(cvfit, s = lambda_use)

selected <- rownames(coef_lasso)[coef_lasso[,1] != 0]
selected <- selected[selected != "(Intercept)"]

cat("Selected proteins:\n")
print(selected)

 # -----------------------------
# Training predictions
# -----------------------------

signature_score <- predict(cvfit, newx=X_scaled, s=lambda_use)

cor(signature_score, y, method="spearman")

# -----------------------------
# LOIO cross-validated predictions
# -----------------------------

# Replace the lambda_index extraction with this:
lambda_index <- which(cvfit$lambda == lambda_use)
cv_pred      <- cvfit$fit.preval[, lambda_index]

cor(cv_pred, y, method="spearman")
cor(cv_pred, y, method="pearson")


plot(cv_pred, y,
     xlab="Predicted frequency (CV)",
     ylab="Observed frequency")

abline(lm(y ~ cv_pred), col="red", lwd=2)

pred_train <- as.numeric(predict(cvfit, newx = X_scaled, s="lambda.min"))

df_plot_data <- data.frame(sample_id = rownames(X_scaled), pred = pred_train) |>
  dplyr::inner_join(
    data.frame(sample_id = sample_ids, obs = as.numeric(y)),
    by="sample_id"
  )

(df_plot=ggplot(df_plot_data, aes(pred, obs)) +
    geom_point() +
    geom_smooth(method="lm") +
    labs(x="Predicted frequency (CV)",
         y="Observed frequency")+
    ggpubr::stat_cor(method="spearman"))



# apply to MICDROP data ####
clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv") %>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks","24 weeks","52 weeks","68 weeks"))) %>%
  filter(targetName %notin% c("CTSS", "LTA|LTB", "IFNA2")) %>%
  mutate(sample_id = sample)

wide_new <- clean_data %>%
  pivot_wider(names_from = targetName, names_prefix = "p_", values_from = conc, id_cols = sample_id)
# ------------------------------
# 2. Extract predictor matrix
# ------------------------------
train_features <- colnames(X)

X_new <- wide_new %>%
  select(sample_id, all_of(train_features)) %>%
  arrange(sample_id)

sample_ids <- X_new$sample_id

X_new <- X_new %>%
  select(-sample_id) %>%
  as.matrix()




# ------------------------------
# 3. Apply training scaling
# ------------------------------

X_new_scaled <- sweep(X_new, 2, train_means[train_features], "-")
X_new_scaled <- sweep(X_new_scaled, 2, train_sds[train_features], "/")
# ------------------------------
# 4. Predict cell frequencies
# ------------------------------

predicted_freq <- predict(cvfit, newx = X_new_scaled, s = cvfit$lambda.min)

# ------------------------------
# 5. Store predictions
# ------------------------------

prediction_df <- data.frame(
  sample_id = sample_ids,
  predicted_tr1_freq = predicted_freq[,1]
)

clean_data_with_tr1_signature <- clean_data%>%
  left_join(., prediction_df, by="sample_id")

# write.csv(clean_data_with_tr1_signature, "~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_tr1_freq.csv")
write.csv(prediction_df, "~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/tr1_freq.csv")

#for pras
clean_data_with_tr1_signature%>%
  distinct(id, timepoint_num, predicted_tr1_freq)%>%
  filter(!is.na(timepoint_num))%>%
  write.csv(., "~/Downloads/tr1_freq.csv")



(p <- clean_data_with_tr1_signature%>%
    filter(mstatus==0, treatmentarm!="DP 2 years")%>%
    distinct(id, timepoint, treatmentarm, predicted_tr1_freq)%>%
    ggplot(., aes(x=timepoint, y=predicted_tr1_freq, fill=treatmentarm))+
    geom_violin()+
    geom_boxplot(color="white",position = position_dodge(width=0.9), width=0.3)+
    ggpubr::stat_compare_means(method = "wilcox.test")+
    scale_fill_manual(values = c('darkred', "darkblue"))+
    theme_minimal()
)


tr1_symp_prob_purf <- clean_data_with_tr1_signature%>%
  filter(mstatus==0, timepoint=="52 weeks")%>%
  mutate(symp_prob_12_24=total_n_malaria_12_24/total_n_para_12_24)%>%
  filter(!is.na(symp_prob_12_24))%>%
  distinct(id, symp_prob_12_24, predicted_tr1_freq, gender_categorical, treatmentarm, total_n_para_12_24)%>%
  nest()%>%
  # mutate(n_malaria_model=map(data, ~glm(symp_prob_12_24~conc+gender_categorical+log_qpcr+total_n_para_12, data=., family = "binomial" , weights = total_n_para_12_24))) %>%
  mutate(n_malaria_model=map(data, ~glm(symp_prob_12_24~predicted_tr1_freq*treatmentarm, data=., family = "binomial"))) %>%
  mutate(tr1_slopes = map(n_malaria_model, ~emmeans::emtrends(.x, ~ treatmentarm, var = "predicted_tr1_freq")))%>%
  mutate(tr1_treatment_slopes=map(n_malaria_model, ~emmeans::emtrends(.x, pairwise ~ treatmentarm, var = "predicted_tr1_freq")))



(symp_prob_12_24_tr1_plot <- clean_data_with_tr1_signature %>%
    filter(
           timepoint=="52 weeks",
           treatmentarm!="DP 2 years",
           !is.na(total_n_para_12_24)
           #log_qpcr<0.001
    )%>%
    mutate(symp_prob_12_24=total_n_malaria_12_24/total_n_para_12_24)%>%
    distinct(id, symp_prob_12_24, predicted_tr1_freq, gender_categorical, treatmentarm, total_n_para_12_24)%>%
    # mutate(symp_prob_12_24f=case_when(symp_prob_12_24<0.25~"< 25%",
    #                                   symp_prob_12_24>=0.25&symp_prob_12_24<0.5~"25-50%",
    #                                   symp_prob_12_24>=0.5&symp_prob_12_24<0.75~"50-75%",
    #                                   symp_prob_12_24>=0.75~"> 75%",
    #                                   total_n_para_12_24==0~"no parasitemia"))%>%
    # mutate(symp_prob_12_24f=factor(symp_prob_12_24f, levels=c("< 25%", "25-50%", "50-75%", "> 75%")))%>%
    mutate(symp_prob_12_24f=case_when(symp_prob_12_24<0.33~"< 33%",
                                      symp_prob_12_24>=0.33&symp_prob_12_24<0.66~"33-66%",
                                      symp_prob_12_24>=0.66~"> 66%",
                                      total_n_para_12_24==0~"no parasitemia"))%>%
    mutate(symp_prob_12_24f=factor(symp_prob_12_24f, levels=c("< 33%", "33-66%", "> 66%")))%>%
    
    # mutate(symp_prob_12_24f=case_when(symp_prob_12_24==0~"2% ",
    #                                   symp_prob_12_24>=0.25&symp_prob_12_24<0.5~"25%-50%",
    #                                   symp_prob_12_24>=0.5&symp_prob_12_24<0.75~"50%-75%",
    #                                   symp_prob_12_24>=0.75~"75% or more",
    #                                   total_n_para_12_24==0~"no parasitemia"))%>%
    filter(symp_prob_12_24f!="no parasitemia")%>%
    ggplot(aes(y=predicted_tr1_freq, x="", fill=symp_prob_12_24f))+
    # geom_violin()+
    geom_boxplot(outliers = F)+
    facet_wrap(~treatmentarm, scales = "fixed", ncol=2)+
    # geom_smooth(method="gam")+
    # ggpubr::stat_compare_means(size=2, hide.ns = T, ref.group = "33% or less")+
    viridis::scale_fill_viridis(discrete = T)+
    ylab("concentration at 52 weeks")+
    xlab("probability of symptoms, given parasitemia in months 12-24")+
    theme_minimal()+
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)
    ))


clean_data_with_tr1_signature %>%
  filter(
    timepoint=="52 weeks",
    treatmentarm!="DP 2 years",
    !is.na(total_n_para_12_24)
    #log_qpcr<0.001
  )%>%
  mutate(symp_prob_12_24=total_n_malaria_12_24/total_n_para_12_24)%>%
  distinct(id, symp_prob_12_24, predicted_tr1_freq, gender_categorical, treatmentarm, total_n_para_12_24)%>%
  ggplot(., aes(y=predicted_tr1_freq, x=symp_prob_12_24))+
  # geom_violin()+
  geom_point(shape=21)+
  geom_smooth(method="gam")+
  # ggpubr::stat_compare_means(size=2, hide.ns = T, ref.group = "33% or less")+
  viridis::scale_fill_viridis(discrete = T)+
  ylab("concentration at 52 weeks")+
  xlab("probability of symptoms, given parasitemia in months 12-24")+
  theme_minimal()+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)
  )
