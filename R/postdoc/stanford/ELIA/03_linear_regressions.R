# modelz ####
## change through time ####
library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)
library(purrr)
library(tidyr)

corrected_long <- read.csv("~/Library/CloudStorage/Box-Box/Florian Bach's Externally Shareable Files/ELIA/limma_corrected_long_format.csv")

time_results <- corrected_long %>%
  group_by(targetName) %>%
  nest() %>%
  mutate(
    model = map(data, ~ purrr::quietly(lme4::lmer)(
      corrected_npq ~ timepoint * study * treatmentarm + gender_categorical + (1 | id),
      data = .x
    )$result),
    
    # MICDROP: DP 1 year vs Placebo at each timepoint
    micdrop_arm = map(model, ~ purrr::quietly(function(m) {
      emmeans(m, ~ treatmentarm | timepoint,
              at = list(study = "MICDROP",
                        treatmentarm = c("DP 1 year", "Placebo"))) %>%
        contrast(method = "pairwise", adjust = "none") %>%
        as.data.frame()
    })(.x)$result),
    
    # PROSYNK: each active arm vs Placebo at each timepoint
    prosynk_arm = map(model, ~ purrr::quietly(function(m) {
      emm <- emmeans(m, ~ treatmentarm | timepoint,
                     at = list(study = "PROSYNK",
                               treatmentarm = c("Lab4b", "Labinic", "Probiotic", "Placebo")))
      contrast(emm, method = "trt.vs.ctrl", ref = "Placebo", adjust = "none") %>%
        as.data.frame()
    })(.x)$result),
    
    # Within-arm timepoint comparisons (adjacent only): both studies
    time_contrasts = map(model, ~ purrr::quietly(function(m) {
      emm <- emmeans(m, ~ timepoint | study * treatmentarm)
      contrast(emm, method = list(
        "8w vs 24w" = c(-1, 1, 0),
        "24w vs 52w" = c(0, -1, 1)
      ), adjust = "none") %>%
        as.data.frame()
    })(.x)$result),
    
    # Cross-study placebo: adjacent slopes only
    placebo_slopes = map(model, ~ purrr::quietly(function(m) {
      emm <- emmeans(m, ~ study * timepoint,
                     at = list(treatmentarm = "Placebo"))
      contrast(emm, interaction = "pairwise", adjust = "none") %>%
        as.data.frame() %>%
        filter(timepoint_pairwise %in% c("8 weeks - 24 weeks", "24 weeks - 52 weeks"))
    })(.x)$result)
    
  ) %>%
  select(targetName, micdrop_arm, prosynk_arm, time_contrasts, placebo_slopes) %>%
  pivot_longer(
    cols = c(micdrop_arm, prosynk_arm, time_contrasts, placebo_slopes),
    names_to = "contrast_type",
    values_to = "df"
  ) %>%
  group_by(contrast_type) %>%
  mutate(df = map(df, ~ mutate(.x, p_adj = p.adjust(p.value, method = "BH")))) %>%
  ungroup()


study_differences <- time_results%>%
  filter(contrast_type == "placebo_slopes") %>%
  unnest(df)%>%
  filter(p_adj<0.05)%>%
  arrange(p_adj)

arm_differences <- time_results %>%
  filter(contrast_type %in% c("micdrop_arm", "prosynk_arm")) %>%
  unnest(df)%>%
  filter(p_adj<0.05)%>%
  arrange(p_adj)



top9_placebo_study_differences <- study_differences%>%
  slice_min(n = 9, with_ties = F, order_by = p_adj)%>%
  pull(targetName)



# placebo / no treatment only ####


placebo_time_results <- clean_data %>%
  filter(timepoint %in% c("24 weeks", "52 weeks"),
         treatmentarm == "Placebo") %>%
  group_by(targetName) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lme4::lmer(normalized_conc ~ timepoint * study + gender_categorical + (1 | id), data = .x)),
    emm = map(model, ~ emmeans(
      .x,
      ~ timepoint * study
    )),
    time_contrasts = map(emm, ~ contrast(
      .x,
      method = "pairwise",
      by = "study",
      adjust = "none"
    ) %>% as.data.frame()),
    study_contrasts = map(model, ~ emmeans(
      .x,
      ~ study | timepoint
    ) %>%
      contrast(method = "pairwise") %>%
      as.data.frame())
  ) %>%
  select(targetName, time_contrasts, study_contrasts) %>%
  pivot_longer(
    cols = c(time_contrasts, study_contrasts),
    names_to = "contrast_type",
    values_to = "df"
  ) %>%
  unnest(df) %>%
  group_by(contrast_type) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

placebo_study_differences <- placebo_time_results%>%
  filter(p_adj<0.05, contrast_type=="study_contrasts")%>%
  arrange(p_adj)

top9_placebo_study_differences <- placebo_study_differences%>%
  slice_min(n = 9, with_ties = F, order_by = p_adj)%>%
  pull(targetName)


clean_data %>%
  filter(timepoint %in% c("8 weeks", "24 weeks", "52 weeks"),
         treatmentarm == "Placebo")%>%
  filter(targetName %in% top9_placebo_study_differences)%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=normalized_conc, fill=study))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName)+
  theme_minimal()

