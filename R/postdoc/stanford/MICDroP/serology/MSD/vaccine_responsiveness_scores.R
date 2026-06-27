## ============================================================
## Vaccine responsiveness summary statistics
## 1. Z-score sum within antigen class (IgG / IgA separately)
## 2. PCA-derived composite score
## ============================================================

library(tidyverse)
library(lme4)
library(emmeans)
library(ggplot2)


mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")
maternal_treatment_arms <- haven::read_dta("~/Library/CloudStorage/Box-Box/DP+SP study/Databases and preliminary findings/Final database used for analyses/DPSP treatment allocation_FINAL.dta")
epi_data <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/micdrop_metadata_for_immunology.csv")
nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")
msd_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/msd_results_v2.csv")
vaccine_coverage <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/long_vaccine_coverage_with_recall.csv")
IU_cutoffs <- readxl::read_excel("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/iu_conversion.xlsx")

diagnoses_of_kids_in_msd <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/diagnoses_of_kids_in_msd.csv")


vaccines = c("Diphtheria",     "Measles" ,         "Pertussis",     "Polio",
             "Rotavirus" ,    "Rubella",       "Tetanus", "Pneumo.1.4.14")
vaccines_iga <- c("Diphtheria_IgA", "Measles_IgA",  "Pertussis_IgA", "Polio_IgA",  "Rotavirus_IgA", "Rubella_IgA",  "Tetanus_IgA",  "Varicella_IgA")
viruses <- c("Mumps",    
             "Varicella",
             "CoV.2",    
             "EV.71",    
             "EV.D68",   
             "Flu.Mix",  
             'hMPV',     
             "PIV.Mix",  
             "RSV",      
             "RV.C",     
             "Flu.A.H1", 
             "Flu.A.H3", 
             "Flu.B",    
             "PIV.1",    
             "PIV.2",    
             'PIV.3',    
             "PIV.4")
cmv_data <-  read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/CMV_ELISA_Uganda.csv", skip = 1)



long_msd <- msd_data%>%
  mutate(sample=paste(SubjectID, "_", "tp", TimePt, sep=""))%>%
  mutate(id=SubjectID, timepoint=paste(TimePt, "weeks"))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))%>%
  select(-SubjectID, -TimePt)%>%
  pivot_longer(cols=-c(sample, id, timepoint), names_to = "antigen", values_to = "titer")%>%
  mutate(Ig_class=ifelse(grepl("*IgA$", antigen), "IgA", "IgG"))%>%
  filter(!is.na(titer))%>%
  mutate(mom_rx=maternal_treatment_arms$treatmentarm[match(id-10000, maternal_treatment_arms$id)],
         CMV_at_52_weeks=cmv_data$Diagnosis[match(id, cmv_data$id)],
         antigen=ifelse(antigen=="Diptheria", "Diphtheria", antigen),
         vaccine_type=case_match(antigen,
                                 "Pneumo.1.4.14"~"PCV10",
                                 "Diphtheria"~"DTP_Hib",
                                 "Diphtheria_IgA"~"DTP_Hib",
                                 "Pertussis"~"DTP_Hib",
                                 "Pertussis_IgA"~"DTP_Hib",
                                 "Tetanus"~"DTP_Hib",
                                 "Measles"~"Measles Rubella",
                                 "Rubella"~"Measles Rubella",
                                 "Polio"~"Sabin Polio",
                                 "Polio_IgA"~"Oral Polio",
                                 "Rotavirus"~"Rotarix",
                                 "Rotavirus_IgA"~"Rotarix", .default = "no vaccine"
         ))%>%
  mutate(mom_rx=case_match(mom_rx, 1~"SP", 2~"DP", 3~"DPSP"))%>%
  left_join(., vaccine_coverage, by=c("id", "vaccine_type"))


# additional data from the moment of sampling
additional_clinical_data <- nulisa_data%>%
  distinct(id, sample, date, ageinwks, gender_categorical, mstatus, qPCRparsdens, pardens)

antibodies_and_epi <- epi_data%>%
  inner_join(., long_msd, by="id")%>%
  left_join(additional_clinical_data, by=c("sample", "id"))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))


final_cutoff_frame <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/final_cutoff_frame.csv")
colnames(final_cutoff_frame)[2]="cut_off"

wide_long_msd <-  long_msd %>%
  filter(!is.na(timepoint))%>%
  pivot_wider(names_from = "timepoint", values_from = "titer", id_cols = c("id", "antigen", "Ig_class"))%>%
  mutate(log2fc_8_24=log2(`24 weeks`/`8 weeks`),
         log2fc_24_52=log2(`52 weeks`/`24 weeks`))%>%
  group_by(antigen)%>%
  mutate("z_fc_8_24"=scale(log2fc_8_24, center = T))%>%
  pivot_longer(cols = c("8 weeks", "24 weeks","52 weeks"), names_to = "timepoint", values_to = "titer")%>%
  pivot_longer(cols = c("log2fc_8_24", "log2fc_24_52"), names_to = "fc_flavour", values_to = "log2fc")%>%
  mutate(conversion=case_when(timepoint=="52 weeks" & fc_flavour=="log2fc_24_52" & log2fc > 2 ~ "converts 24 to 52",
                              timepoint=="24 weeks" & fc_flavour=="log2fc_8_24" & z_fc_8_24 > 2 | 
                                timepoint=="24 weeks" & fc_flavour=="log2fc_8_24" & log2fc > 1 ~ "converts 8 to 24",
                              .default="nothing"))%>%
  left_join(., final_cutoff_frame, by="antigen")%>%
  mutate(sero_positive=if_else(titer>=cut_off, 1, 0))%>%
  mutate(conversion=factor(conversion, levels=c("nothing", "converts 8 to 24", "converts 24 to 52")))

df_pred <- antibodies_and_epi%>%
  distinct(id, antigen, timepoint, titer, qPCRparsdens, course_outcome) %>%
  filter(timepoint %in% c("24 weeks", "52 weeks"),
         antigen %in% c(vaccines, vaccines_iga)) %>%
  rename(
    vaccine_antigen = antigen,
    vaccine_titer = titer
  )

df_target <- wide_long_msd%>%
  filter(timepoint == "24 weeks" & fc_flavour=="log2fc_8_24" |
          timepoint == "52 weeks" & fc_flavour=="log2fc_24_52",
          antigen %in% viruses) %>%
  distinct(id, timepoint, antigen, conversion) %>%
  pivot_wider(values_from = conversion, names_from = timepoint)%>%
  mutate(sero_conversion=case_when(`24 weeks`=="nothing"&`52 weeks`=="nothing"~"nothing",
                                   `24 weeks`=="nothing"&`52 weeks`=="converts 24 to 52"~"converts 24 to 52",
                                   `24 weeks`=="converts 8 to 24"&`52 weeks`=="nothing"~"converts 8 to 24", .default = ))%>%
  distinct(id, antigen, sero_conversion)%>%
  rename(
    virus_antigen = antigen,
    virus_seroconversion = sero_conversion
  )%>%
  filter(!is.na(virus_seroconversion))



wide_long_msd_cross <- df_pred %>%
  inner_join(df_target, by = "id") %>%
  mutate(log2_titer = log2(vaccine_titer))%>%
  mutate(id_cat=factor(as.character(id)))

personal_seroconversions <- wide_long_msd_cross %>%
  distinct(id, virus_antigen, virus_seroconversion) %>%
  group_by(id)%>%
  summarise(n_seroconversions = sum(virus_seroconversion != "nothing"),
            early_seroconversions = sum(virus_seroconversion == "converts 8 to 24"),
            late_seroconversions = sum(virus_seroconversion == "converts 24 to 52"))%>%
  arrange(desc(n_seroconversions))

wide_long_msd_cross <- wide_long_msd_cross%>%
  left_join(., personal_seroconversions, by="id")
## ---- 0. Define antigen classes ------------------------------
## Adjust these vectors to match your exact antigen names

antigens_igg <- c("Diphtheria", "Measles", "Pertussis",
                  "Pneumo.1.4.14", "Polio", "Rotavirus", "Rubella")

antigens_iga <- c("Pertussis_IgA", "Polio_IgA", "Rotavirus_IgA")

## ---- 1. Z-score sum -----------------------------------------
## Normalise each antigen to z-score within timepoint and cohort,
## then sum within IgG and IgA classes per individual

zscore_summary <- wide_long_msd_cross %>%
  distinct(id, timepoint, vaccine_antigen, vaccine_titer) %>%
  group_by(timepoint, vaccine_antigen) %>%
  mutate(z = scale(log2(vaccine_titer), center = TRUE, scale = TRUE)[,1]) %>%
  ungroup() %>%
  mutate(antigen_class = case_when(
    vaccine_antigen %in% antigens_igg ~ "IgG",
    vaccine_antigen %in% antigens_iga ~ "IgA",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(antigen_class)) %>%
  group_by(id, timepoint, antigen_class) %>%
  summarise(
    z_sum   = sum(z, na.rm = TRUE),
    n_antigens = sum(!is.na(z)),
    .groups = "drop"
  ) %>%
  ## Only keep individuals with complete data across antigens
  filter(n_antigens == case_when(
    antigen_class == "IgG" ~ length(antigens_igg),
    antigen_class == "IgA" ~ length(antigens_iga)
  ))

## Audit
zscore_summary %>%
  count(timepoint, antigen_class) %>%
  print()

## ---- 2. PCA composite score ---------------------------------
## Run PCA on log2 titer matrix (individuals x antigens)
## Extract PC1 as composite responsiveness axis
## Run separately per timepoint

all_antigens <- c(antigens_igg, antigens_iga)

## Build wide matrix: one row per individual per timepoint
wide_vaccine <- wide_long_msd_cross %>%
  distinct(id, timepoint, vaccine_antigen, log2_titer) %>%
  filter(vaccine_antigen %in% all_antigens) %>%
  pivot_wider(
    id_cols     = c(id, timepoint),
    names_from  = vaccine_antigen,
    values_from = log2_titer
  ) %>%
  drop_na()   # PCA requires complete cases

## Run PCA per timepoint and extract PC1
pca_scores <- wide_vaccine %>%
  group_by(timepoint) %>%
  group_modify(function(dat, grp) {

    mat <- dat %>%
      select(all_of(all_antigens)) %>%
      as.matrix()

    pca <- prcomp(mat, center = TRUE, scale. = TRUE)

    ## Variance explained
    var_exp <- summary(pca)$importance[2, 1:3]
    message(sprintf(
      "Timepoint %s — PC1: %.1f%%  PC2: %.1f%%  PC3: %.1f%%",
      grp$timepoint, var_exp[1]*100, var_exp[2]*100, var_exp[3]*100
    ))

    ## PC1 loadings — check sign convention:
    ## positive loading = higher titer -> higher PC1 (expected)
    ## if majority of loadings are negative, flip PC1
    loadings <- pca$rotation[, 1]
    if (mean(loadings) < 0) {
      message("Flipping PC1 sign so higher score = better response")
      pca$x[, 1] <- -pca$x[, 1]
    }

    data.frame(
      id      = dat$id,
      pc1     = pca$x[, 1],
      pc2     = pca$x[, 2],
      pct_var = var_exp[1] * 100
    )
  }) %>%
  ungroup()

## Loadings for interpretation
pca_loadings <- wide_vaccine %>%
  group_by(timepoint) %>%
  group_modify(function(dat, grp) {
    mat <- dat %>% select(all_of(all_antigens)) %>% as.matrix()
    pca <- prcomp(mat, center = TRUE, scale. = TRUE)
    loadings <- pca$rotation[, 1]
    if (mean(loadings) < 0) loadings <- -loadings
    data.frame(antigen = names(loadings), loading = loadings)
  }) %>%
  ungroup()

## ---- 3. Combine into one summary dataframe ------------------

vaccine_responsiveness <- zscore_summary %>%
  pivot_wider(
    id_cols     = c(id, timepoint),
    names_from  = antigen_class,
    values_from = z_sum,
    names_prefix = "zscore_"
  ) %>%
  full_join(
    pca_scores %>% select(id, timepoint, pc1, pct_var),
    by = c("id", "timepoint")
  )

## ---- 4. Join with seroconversion burden ---------------------
## Assumes n_seroconversions column exists in wide_long_msd_cross
## from option 1 analysis; adjust if named differently

burden <- wide_long_msd_cross %>%
  distinct(id, timepoint, n_seroconversions)

vaccine_responsiveness <- vaccine_responsiveness %>%
  left_join(burden, by = c("id", "timepoint"))

## ---- 5. Model: seroconversion burden -> vaccine responsiveness
## One model per summary score per timepoint

outcomes <- c("zscore_IgG", "zscore_IgA", "pc1")

burden_models <- expand_grid(
  outcome  = outcomes,
  tp       = unique(vaccine_responsiveness$timepoint)
) %>%
  mutate(
    data = map2(outcome, tp, ~ {
      vaccine_responsiveness %>%
        filter(timepoint == .y) %>%
        select(id, y = all_of(.x), n_seroconversions) %>%
        drop_na()
    }),
    n      = map_int(data, nrow),
    model  = map2(data, outcome, ~ tryCatch(
      lm(y ~ n_seroconversions, data = .x),
      error = function(e) NULL
    )),
    tidy   = map(model, ~ if (inherits(.x, "lm"))
      broom::tidy(.x, conf.int = TRUE) else NULL)
  ) %>%
  filter(!map_lgl(tidy, is.null)) %>%
  select(outcome, tp, n, tidy) %>%
  unnest(tidy) %>%
  filter(term == "n_seroconversions") %>%
  group_by(outcome) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

print(burden_models %>%
  select(outcome, tp, n, estimate, std.error, p.value, padj))

## ---- 6. Visualise -------------------------------------------

## 6a. PCA loadings heatmap — what drives PC1?
ggplot(pca_loadings, aes(x = timepoint, y = antigen, fill = loading)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(loading, 2)), size = 3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, name = "PC1 loading") +
  labs(title = "PC1 loadings by timepoint",
       subtitle = "Positive = higher titer contributes to higher responsiveness score",
       x = "Timepoint", y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

## 6b. Seroconversion burden vs PC1 score per timepoint
vaccine_responsiveness %>%
  drop_na(pc1, n_seroconversions) %>%
  ggplot(aes(x = n_seroconversions, y = pc1)) +
  geom_point(alpha = 0.6, size = 1.5, position = position_jitter(width=0.2)) +
  geom_smooth(method = "lm", colour = "firebrick", se = TRUE) +
  ggpubr::stat_cor(method = "spearman", label.x.npc = 0.05) +
  facet_wrap(~ timepoint) +
  labs(x = "Number of viral seroconversions",
       y = "PC1 vaccine responsiveness score",
       title = "Viral exposure burden vs. vaccine responsiveness") +
  theme_classic(base_size = 11)


## 6c. Same for z-score sum
vaccine_responsiveness %>%
  drop_na(zscore_IgG, n_seroconversions) %>%
  ggplot(aes(x = n_seroconversions, y = zscore_IgG)) +
  geom_point(alpha = 0.6, size = 1.5, position = position_jitter(width=0.2)) +
  geom_smooth(method = "lm", colour = "steelblue", se = TRUE) +
  ggpubr::stat_cor(method = "spearman", label.x.npc = 0.05) +
  facet_wrap(~ timepoint) +
  labs(x = "Number of viral seroconversions",
       y = "IgG z-score sum",
       title = "Viral exposure burden vs. IgG vaccine responsiveness") +
  theme_classic(base_size = 11)

vaccine_responsiveness %>%
  filter(timepoint == "52 weeks") %>%
  drop_na(zscore_IgG, n_seroconversions) %>%
  ggplot(aes(x = n_seroconversions, y = zscore_IgG)) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  geom_smooth(method = "lm", colour = "firebrick") +
  geom_smooth(method = "loess", colour = "steelblue", linetype = "dashed") +
  theme_classic()

## 6d. Biaxial PC1 vs PC2 coloured by seroconversion burden
pca_scores %>%
  left_join(burden, by = c("id", "timepoint")) %>%
  drop_na(n_seroconversions) %>%
  ggplot(aes(x = pc1, y = pc2, colour = n_seroconversions)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_viridis_c(name = "N seroconversions") +
  facet_wrap(~ timepoint) +
  labs(
    x        = "PC1 (vaccine responsiveness)",
    y        = "PC2",
    title    = "PCA of vaccine antibody titers",
    subtitle = sprintf("PC1: %.1f%% variance explained", 
                       mean(pca_scores$pct_var))
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")
