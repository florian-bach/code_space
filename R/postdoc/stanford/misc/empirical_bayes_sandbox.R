library(lme4)
library(emmeans)
library(tidyverse)

#fold change cutoff test function
# moderate_treat <- function(df, lfc = 0.5){
# 
#   # Empirical Bayes shrinkage of variances
#   squeeze <- squeezeVar(df$SE^2, df$df)
#   se_post <- sqrt(squeeze$var.post)
#   df_post <- squeeze$df.prior + df$df
#   
#   # TREAT statistic (tests |effect| > lfc)
#   t_treat <- (abs(df$estimate) - lfc) / se_post
#   
#   # one-sided p-value because we already use abs()
#   p_treat <- 2 * pt(-abs(t_treat), df = df_post)
#   
#   df %>%
#     mutate(
#       t_treat = t_treat,
#       p_treat = p_treat
#     )
# }

moderate_treat <- function(tbl, lfc = 0.25) {
  tbl %>%
    mutate(
      t_treat = (abs(estimate) - lfc) / se_post,
      ## Two-sided: abs() folds both tails; multiply by 2
      p_treat = 2 * pt(-abs(t_treat), df = df_post),
      ## Clip to [0,1] — numerically t_treat can be very negative
      ## giving pt() > 0.5, doubling past 1
      p_treat = pmin(p_treat, 1)
    )
}

clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")

clean_data <- clean_data %>%
  mutate(timepoint = factor(timepoint, levels=c("bad_baseline", "baseline", "day0", "day7", "day14", "day28")))

# build model 
fits <- clean_data %>%
  filter(infectiontype %in% c("A","S"),
         !timepoint %in% c("day28","bad_baseline")) %>%
  mutate(timepoint = factor(timepoint,
                            levels=c("baseline","day0","day7","day14"))) %>%
  group_by(targetName) %>%
  nest() %>%
  mutate(model = map(data, ~ lmer(
    concentration ~ timepoint*infectiontype +
      ageyrs + gender_categorical + (1|id_cat),
    data = .
  )))

custom_contrasts <- list(
  "baseline A - baseline S" = c( 1,0,0,0, -1,0,0,0),
  "baseline A - day0 A"     = c( 1,-1,0,0,  0,0,0,0),
  "baseline A - day14 A"    = c( 1,0,0,-1,  0,0,0,0),
  "baseline S - day0 S"     = c( 0,0,0,0,  1,-1,0,0),
  "baseline S - day7 S"     = c( 0,0,0,0,  1,0,-1,0),
  "baseline S - day14 S"    = c( 0,0,0,0,  1,0,0,-1)
)

#extract emmeans constasts, estimates, se df
fits <- fits %>%
  mutate(
    emm = map(model, ~ emmeans(.x, ~ timepoint * infectiontype)),
    cont = map(emm, ~ contrast(.x, custom_contrasts) %>% summary(infer = TRUE))
  )

# unnnest into tidy table 
contrast_tbl <- fits %>%
  select(targetName, cont) %>%
  unnest(cont) %>%
  select(targetName, contrast, estimate, SE, df)

eb_results <- contrast_tbl %>%
  group_by(contrast) %>%
  group_modify(~ moderate_treat(.x,)) %>%
  ungroup()

final_results <- eb_results %>%
  group_by(contrast) %>%
  mutate(padj = p.adjust(p_treat, method="fdr")) %>%
  ungroup()

final_results%>%
  filter(padj<0.05)%>%
  group_by(contrast)%>%
  count()


limma_sigs <- final_results%>%
  filter(padj<0.05)%>%
  filter(contrast=="baseline A - day0 A")%>%
  arrange(padj)%>%
  pull(targetName)
  

emmeans_sigs <- combo_as_purff%>%
  filter(padj<0.05)%>%
  filter(contrast=="baseline A - day0 A")%>%
  arrange(padj)%>%
  pull(targetName)



# this tests for significance WITHOUT an effect size a cutoff
moderate_contrast <- function(df) {
  squeeze <- squeezeVar(df$SE^2, df$df)
  t_mod <- df$estimate / sqrt(squeeze$var.post)
  df_mod <- squeeze$df.prior + df$df
  p_mod <- 2 * pt(-abs(t_mod), df = df_mod)
  
  df %>%
    mutate(
      t_moderated = t_mod,
      df_moderated = df_mod,
      p_moderated = p_mod
    )
}
