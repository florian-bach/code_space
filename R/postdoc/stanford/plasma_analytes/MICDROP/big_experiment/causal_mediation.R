# mediation try ####
library(mediation)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks", "68 weeks")))%>%
  mutate(wealthcat=factor(wealthcat),
         educ=factor(educ))%>%
  filter(!targetName %in% c("CTSS", "LTA|LTB", "IFNA2"), 
         timepoint!="68 weeks", treatmentarm=="DP 1 year")

run_med_long <- function(df) {
  
  # mediator model (same for every analyte, but that's fine)
  med.fit <- lm(log_qpcr ~ total_n_para_12, data=df)
  
  # outcome model
  out.fit <- lm(conc~total_n_para_12+log_qpcr+gender_categorical+wealthcat+educ, data=df)
  
  med.out <- mediate(
    med.fit,
    out.fit,
    treat = "total_n_para_12",
    mediator = "log_qpcr",
    sims = 1000,
    boot = FALSE
  )
  
  tibble(
    ACME  = med.out$d.avg,
    ADE   = med.out$z.avg,
    TE    = med.out$tau.coef,
    PM    = med.out$n.avg,
    p_ACME = med.out$d.avg.p,
    p_ADE  = med.out$z.avg.p,
    p_TE   = med.out$tau.p
  )
}

system.time(
  results <- clean_data %>%
    filter(timepoint=="52 weeks", mstatus==0, !is.na(wealthcat))%>%
    group_by(targetName) %>%
    nest() %>%
    mutate(med = map(data, run_med_long)) %>%
    unnest(med) %>%
    ungroup()
)


results <- results %>%
  mutate(
    q_ACME = p.adjust(p_ACME, method = "fdr"),
    q_ADE  = p.adjust(p_ADE,  method = "fdr"),
    q_TE   = p.adjust(p_TE,   method = "fdr"),
    prop_mediated = ACME / TE
  )

results <- results %>%
  mutate(
    mediation_class = case_when(
      q_ACME < 0.05 & q_ADE >= 0.05 ~ "Mediated",
      q_ADE  < 0.05 & q_ACME >= 0.05 ~ "Direct",
      q_ACME < 0.05 & q_ADE  < 0.05 ~ "Both",
      q_TE   < 0.05 ~ "Total only",
      TRUE ~ "None"
    )
  )
# write.csv(results[,-2], "~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/causal_mediation_results_para_prev_12_with_confounders.csv", row.names = F)
results%>%
  group_by(mediation_class)%>%
  count()

## with dp
# prev:
#   Both     Direct   Mediated       None Total only 
# 32         31         34        148          2
# n_para:
#   Both     Direct   Mediated       None Total only 
# 23         27         44        151          2 
# 

# without dp
# prev:
#   Both     Direct   Mediated       None Total only 
# 8         17         64        154          4
# n_para:
#   Both     Direct   Mediated       None Total only 
# 5         9         68        162          3 
# n_malaria:
#   Both     Direct   Mediated       None Total only 
# 0         0         40        207          0 

# Category	Criteria	Interpretation
# Mediated	q_ACME < 0.05	Driven by pathogen load
# Direct effect	q_ADE < 0.05	Evidence for immune imprinting
# Total effect only	q_TE < 0.05 but q_ACME ≥ 0.05	Overall association but unclear pathway
# No association	none significant	No evidence

# ACME > 0: infections increase immune markers via pathogen load
# ADE < 0 or small: infections may reduce immune markers when controlling for current load


# Step 1
# Children with more infections end up with higher pathogen load at 12 months.
# 
# Step 2
# Current pathogen load strongly drives immune activation.
# 
# Step 3
# Once current pathogen load is accounted for, infection history itself has little independent effect.

# 5. The proportion mediated (> 1 !!)
# Prop mediated = 1.21 (121%)
# 
# This is the most confusing line.
# 
# How can mediation be >100%?
#   
#   This happens when:
#   
#   ACME and ADE have opposite signs.
# 
# Check the math:
#   
#   ACME  = +2.10
# ADE   = −0.36
# Total =  1.74
# 
# So pathogen load explains more than the total effect, and the direct pathway slightly cancels it.
# 
# This is called inconsistent mediation or suppression.
# 
# Interpretation:
#   
#   Infection history increases immune activation only because it increases pathogen load.
# Any direct effect of past infections may even slightly reduce immune activation.
# 
# This is actually biologically fascinating.
direct_effect_sigs <- results%>%
  filter(mediation_class=="Direct")%>%
  pull(targetName)

clean_data %>%
  filter(timepoint=="52 weeks", mstatus==0, !is.na(wealthcat))%>%
  filter(targetName %in% direct_effect_sigs)%>%
  ggplot(., aes(x=para_prev_12, y=conc))+
    geom_smooth(method="lm")+
    geom_point()+
  ggpubr::stat_cor(method="spearman", color="red")+
    facet_wrap(~targetName, nrow=3, scales="free")+
    theme_minimal()


clean_data %>%
  filter(timepoint=="52 weeks", mstatus==0, !is.na(wealthcat))%>%
  filter(targetName %in% direct_effect_sigs)%>%
  ggplot(., aes(x=total_n_para_12, y=conc))+
  geom_smooth(method="lm")+
  geom_point()+
  ggpubr::stat_cor(method="spearman", color="red")+
  facet_wrap(~targetName, nrow=3, scales="free")+
  theme_minimal()

results%>%
  group_by(mediation_class)%>%
  count()
