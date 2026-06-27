library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(emmeans)

`%notin%` <- Negate(`%in%`)

clean_data <- read.csv("~/Library/CloudStorage/Box-Box/Florian Bach's Externally Shareable Files/ELIA/corrected_prosynk_micdrop_combined.csv")


# MICDROP ####
raw_micdrop_treatment_purf <-  clean_data %>%
  filter(study == "MICDROP", !is.na(conc_limma)) %>%
  group_by(targetName) %>%
  mutate(timepoint = factor(timepoint, levels = c("8 weeks", "24 weeks", "52 weeks"))) %>%
  nest() %>%
  mutate(time_model = map(data, ~lme4::lmer(conc ~ timepoint * treatmentarm + gender_categorical + (1|id), data = .))) %>%
  mutate(emm = map(time_model, ~emmeans(., specs = pairwise ~ treatmentarm | timepoint))) %>%
  mutate(emm_contrast_summary = map(emm, ~contrast(., "pairwise", adjust = "none") %>% summary())) %>%
  # p-values
  mutate("8 weeks_p"  = map_dbl(emm_contrast_summary, ~.$p.value[1]),
         "24 weeks_p" = map_dbl(emm_contrast_summary, ~.$p.value[2]),
         "52 weeks_p" = map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  # estimates
  mutate("8 weeks_est"  = map_dbl(emm_contrast_summary, ~.$estimate[1]),
         "24 weeks_est" = map_dbl(emm_contrast_summary, ~.$estimate[2]),
         "52 weeks_est" = map_dbl(emm_contrast_summary, ~.$estimate[3])) %>%
  pivot_longer(
    cols = ends_with(c("_p", "_est")),
    names_to = c("contrast", ".value"),
    names_pattern = "(.+)_(p|est)"
  ) %>%
  group_by(contrast) %>%
  mutate(padj = p.adjust(p, method = "fdr")) %>%
  ungroup()

raw_sigs <- raw_micdrop_treatment_purf%>%
  filter(padj<0.05, contrast=="52 weeks")
  

raw_micdrop_treatment_purf %>%
  mutate(
    neg_log10_p = -log10(p),
    sig = case_when(
      padj < 0.05 & est > 0 ~ "Up",
      padj < 0.05 & est < 0 ~ "Down",
      TRUE                   ~ "NS"
    ),
    label = ifelse(padj < 0.05, targetName, NA)
  ) %>%
  ggplot(aes(x = est, y = neg_log10_p, colour = sig, label = label)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_text_repel(size = 2.5, max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(
    values = c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey70"),
    name = NULL
  ) +
  facet_wrap(~ contrast) +
  labs(
    x = "Estimated difference (treatment vs control)",
    y = expression(-log[10](p)),
    title = "MICDROP — treatment arm contrasts"
  ) +
  theme_minimal(base_size = 11)








limma_micdrop_treatment_purf <- clean_data %>%
  filter(study == "MICDROP", !is.na(conc_limma)) %>%
  group_by(targetName) %>%
  mutate(timepoint = factor(timepoint, levels = c("8 weeks", "24 weeks", "52 weeks"))) %>%
  nest() %>%
  mutate(time_model = map(data, ~lme4::lmer(conc_limma ~ timepoint * treatmentarm + gender_categorical + (1|id), data = .))) %>%
  mutate(emm = map(time_model, ~emmeans(., specs = pairwise ~ treatmentarm | timepoint))) %>%
  mutate(emm_contrast_summary = map(emm, ~contrast(., "pairwise", adjust = "none") %>% summary())) %>%
  # p-values
  mutate("8 weeks_p"  = map_dbl(emm_contrast_summary, ~.$p.value[1]),
         "24 weeks_p" = map_dbl(emm_contrast_summary, ~.$p.value[2]),
         "52 weeks_p" = map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  # estimates
  mutate("8 weeks_est"  = map_dbl(emm_contrast_summary, ~.$estimate[1]),
         "24 weeks_est" = map_dbl(emm_contrast_summary, ~.$estimate[2]),
         "52 weeks_est" = map_dbl(emm_contrast_summary, ~.$estimate[3])) %>%
  pivot_longer(
    cols = ends_with(c("_p", "_est")),
    names_to = c("contrast", ".value"),
    names_pattern = "(.+)_(p|est)"
  ) %>%
  group_by(contrast) %>%
  mutate(padj = p.adjust(p, method = "fdr")) %>%
  ungroup()


limma_sigs <- limma_micdrop_treatment_purf%>%
  filter(padj<0.05, contrast=="52 weeks")





limma_micdrop_treatment_purf %>%
  mutate(
    neg_log10_p = -log10(p),
    sig = case_when(
      padj < 0.05 & est > 0 ~ "Up",
      padj < 0.05 & est < 0 ~ "Down",
      TRUE                   ~ "NS"
    ),
    label = ifelse(padj < 0.05, targetName, NA)
  ) %>%
  ggplot(aes(x = est, y = neg_log10_p, colour = sig, label = label)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_text_repel(size = 2.5, max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(
    values = c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey70"),
    name = NULL
  ) +
  facet_wrap(~ contrast) +
  labs(
    x = "Estimated difference (treatment vs control)",
    y = expression(-log[10](p)),
    title = "MICDROP — treatment arm contrasts"
  ) +
  theme_minimal(base_size = 11)




combat_micdrop_treatment_purf <- clean_data%>%
  filter(study=="MICDROP", !is.na(conc_combat))%>%
  group_by(targetName)%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(conc_combat~timepoint*treatmentarm+gender_categorical+(1|id), data=.))) %>%
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

combat_sigs <- combat_micdrop_treatment_purf%>%
  filter(padj<0.05, contrast=="52 weeks")


# PROSYNK individual treatments ####
raw_PROSYNK_treatment <- clean_data %>%
  filter(study == "PROSYNK") %>%
  mutate(timepoint = factor(timepoint, levels = c("8 weeks", "24 weeks", "52 weeks"))) %>%
  group_by(targetName) %>%
  nest() %>%
  mutate(
    model = map(data, ~ purrr::quietly(lme4::lmer)(
      conc ~ timepoint * treatmentarm + gender_categorical + (1 | id),
      data = .x
    )$result),
    
    # treatmentarm comparisons within each timepoint
    arm_contrasts = map(model, ~ emmeans(
      .x, ~ treatmentarm | timepoint
    ) %>%
      contrast("trt.vs.ctrl", ref = "Placebo", adjust = "none") %>%
      as.data.frame()),
    
    # timepoint comparisons within each arm
    time_contrasts = map(model, ~ emmeans(
      .x, ~ timepoint | treatmentarm
    ) %>%
      contrast("pairwise", adjust = "none") %>%
      as.data.frame())
  ) %>%
  select(targetName, arm_contrasts, time_contrasts) %>%
  pivot_longer(
    cols = c(arm_contrasts, time_contrasts),
    names_to = "contrast_type",
    values_to = "df"
  ) %>%
  group_by(contrast_type) %>%
  mutate(df = map(df, ~ mutate(.x, p_adj = p.adjust(p.value, method = "BH")))) %>%
  ungroup()

# access results — contrast labels are automatic from emmeans
raw_PROSYNK_treatment %>%
  filter(contrast_type == "arm_contrasts") %>%
  unnest(df) %>%
  filter(p_adj < 0.1) %>%
  arrange(p_adj)



volcano_df <- raw_PROSYNK_treatment %>%
  unnest(df) %>%
  mutate(
    neg_log10_p = -log10(p.value),
    sig = case_when(
      p_adj < 0.05 & estimate > 0  ~ "Up",
      p_adj < 0.05 & estimate < 0  ~ "Down",
      TRUE                          ~ "NS"
    ),
    label = ifelse(p_adj < 0.05, targetName, NA)
  )

# --- arm contrasts volcano (one panel per timepoint) ---
volcano_df %>%
  filter(contrast_type == "arm_contrasts") %>%
  ggplot(aes(x = estimate, y = neg_log10_p, colour = sig, label = label)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "solid", colour = "grey80") +
  geom_text_repel(
    size = 2.5, max.overlaps = 20,
    box.padding = 0.3, show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey70"),
    name = NULL
  ) +
  facet_wrap(~ paste(contrast, timepoint, sep = "\n"), scales = "free_y") +
  labs(
    x = "Estimated difference (vs Placebo)",
    y = expression(-log[10](p)),
    title = "PROSYNK — treatment arm contrasts"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(size = 8))








limma_PROSYNK_treatment <- clean_data %>%
  filter(study == "PROSYNK", !is.na(conc_limma)) %>%
  mutate(timepoint = factor(timepoint, levels = c("8 weeks", "24 weeks", "52 weeks"))) %>%
  group_by(targetName) %>%
  nest() %>%
  mutate(
    model = map(data, ~ purrr::quietly(lme4::lmer)(
      conc_limma ~ timepoint * treatmentarm + gender_categorical + (1 | id),
      data = .x
    )$result),
    
    # treatmentarm comparisons within each timepoint
    arm_contrasts = map(model, ~ emmeans(
      .x, ~ treatmentarm | timepoint
    ) %>%
      contrast("trt.vs.ctrl", ref = "Placebo", adjust = "none") %>%
      as.data.frame()),
    
    # timepoint comparisons within each arm
    time_contrasts = map(model, ~ emmeans(
      .x, ~ timepoint | treatmentarm
    ) %>%
      contrast("pairwise", adjust = "none") %>%
      as.data.frame())
  ) %>%
  select(targetName, arm_contrasts, time_contrasts) %>%
  pivot_longer(
    cols = c(arm_contrasts, time_contrasts),
    names_to = "contrast_type",
    values_to = "df"
  ) %>%
  group_by(contrast_type) %>%
  mutate(df = map(df, ~ mutate(.x, p_adj = p.adjust(p.value, method = "BH")))) %>%
  ungroup()

# access results — contrast labels are automatic from emmeans
limma_PROSYNK_treatment %>%
  filter(contrast_type == "arm_contrasts") %>%
  unnest(df) %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj)%>%
  View()


combat_PROSYNK_treatment <- clean_data %>%
  filter(study == "PROSYNK", !is.na(conc_combat)) %>%
  mutate(timepoint = factor(timepoint, levels = c("8 weeks", "24 weeks", "52 weeks"))) %>%
  group_by(targetName) %>%
  nest() %>%
  mutate(
    model = map(data, ~ purrr::quietly(lme4::lmer)(
      conc_combat ~ timepoint * treatmentarm + gender_categorical + (1 | id),
      data = .x
    )$result),
    
    # treatmentarm comparisons within each timepoint
    arm_contrasts = map(model, ~ emmeans(
      .x, ~ treatmentarm | timepoint
    ) %>%
      contrast("trt.vs.ctrl", ref = "Placebo", adjust = "none") %>%
      as.data.frame()),
    
    # timepoint comparisons within each arm
    time_contrasts = map(model, ~ emmeans(
      .x, ~ timepoint | treatmentarm
    ) %>%
      contrast("pairwise", adjust = "none") %>%
      as.data.frame())
  ) %>%
  select(targetName, arm_contrasts, time_contrasts) %>%
  pivot_longer(
    cols = c(arm_contrasts, time_contrasts),
    names_to = "contrast_type",
    values_to = "df"
  ) %>%
  group_by(contrast_type) %>%
  mutate(df = map(df, ~ mutate(.x, p_adj = p.adjust(p.value, method = "BH")))) %>%
  ungroup()

# access results — contrast labels are automatic from emmeans
combat_PROSYNK_treatment %>%
  filter(contrast_type == "arm_contrasts") %>%
  unnest(df) %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj)





volcano_df <- combat_PROSYNK_treatment %>%
  unnest(df) %>%
  mutate(
    neg_log10_p = -log10(p.value),
    sig = case_when(
      p_adj < 0.05 & estimate > 0  ~ "Up",
      p_adj < 0.05 & estimate < 0  ~ "Down",
      TRUE                          ~ "NS"
    ),
    label = ifelse(p_adj < 0.05, targetName, NA)
  )

# --- arm contrasts volcano (one panel per timepoint) ---
volcano_df %>%
  filter(contrast_type == "arm_contrasts") %>%
  ggplot(aes(x = estimate, y = neg_log10_p, colour = sig, label = label)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "solid", colour = "grey80") +
  geom_text_repel(
    size = 2.5, max.overlaps = 20,
    box.padding = 0.3, show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey70"),
    name = NULL
  ) +
  facet_wrap(~ paste(contrast, timepoint, sep = "\n"), scales = "free_y") +
  labs(
    x = "Estimated difference (vs Placebo)",
    y = expression(-log[10](p)),
    title = "PROSYNK — treatment arm contrasts (combat)"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(size = 8))









# PROSYNK all treatment groups collapsed ####
raw_PROSYNK_treatment_purf <- clean_data%>%
  filter(study=="PROSYNK")%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))%>%
  mutate(any_probiotic=ifelse(treatmentarm=="Placebo"&study=="PROSYNK", "no probiotic", "probiotic"))%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(conc~timepoint*any_probiotic+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ any_probiotic | timepoint)))%>%
  mutate(emm2=map(time_model, ~emmeans(., specs = pairwise ~ timepoint | any_probiotic)))%>%
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

raw_sigs <- raw_PROSYNK_treatment_purf%>%
  filter(padj<0.05)




limma_PROSYNK_treatment_purf <- clean_data%>%
  filter(study=="PROSYNK", !is.na(conc_limma))%>%
  group_by(targetName)%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))%>%
  mutate(any_probiotic=ifelse(treatmentarm=="Placebo"&study=="PROSYNK", "no probiotic", "probiotic"))%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(conc_limma~timepoint*any_probiotic+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ any_probiotic | timepoint)))%>%
  mutate(emm2=map(time_model, ~emmeans(., specs = pairwise ~ timepoint | any_probiotic)))%>%
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

limma_sigs <- limma_PROSYNK_treatment_purf%>%
  filter(padj<0.05)




combat_PROSYNK_treatment_purf <- clean_data%>%
  filter(study=="PROSYNK", !is.na(conc_combat))%>%
  group_by(targetName)%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))%>%
  mutate(any_probiotic=ifelse(treatmentarm=="Placebo"&study=="PROSYNK", "no probiotic", "probiotic"))%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(conc_combat~timepoint*any_probiotic+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ any_probiotic | timepoint)))%>%
  mutate(emm2=map(time_model, ~emmeans(., specs = pairwise ~ timepoint | any_probiotic)))%>%
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

combat_sigs <- combat_PROSYNK_treatment_purf%>%
  filter(padj<0.05)
