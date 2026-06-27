library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(emmeans)

`%notin%` <- Negate(`%in%`)

micdrop_nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks", "68 weeks")))%>%
  filter(targetName %notin% c("CTSS", "LTA|LTB", "IFNA2"))

slim_micdrop_data <- micdrop_nulisa_data%>%
  filter(mstatus==0, treatmentarm!="DP 2 years")%>%
  mutate(id=paste("micdrop_", id, sep=''))%>%
  select(id, timepoint, anyDP, treatmentarm, gender_categorical, targetName, conc)%>%
  mutate(study="MICDROP")
  
prosynk_nulisa_data <- readRDS("~/Library/CloudStorage/Box-Box/Florian Bach's Externally Shareable Files/ELIA/PROSYNK_data_for_Stanford_ELIA/Data/rds/NULISAseq_Data.rds")
prosynk_enrolment_data <- readRDS("~/Library/CloudStorage/Box-Box/Florian Bach's Externally Shareable Files/ELIA/PROSYNK_data_for_Stanford_ELIA/Data/rds/Enrolment_Data.rds")

slim_prosynk_data <- prosynk_nulisa_data%>%
  left_join(., prosynk_enrolment_data, by="IDNumber")%>%
  mutate(id=IDNumber, timepoint=case_match(VisitID,
                                               1~"8 weeks",
                                               2~"12 weeks",
                                               3~"24 weeks",
                                               4~"52 weeks"
                                              ),
         treatmentarm=case_match(Arm,
                                 1~"Labinic",
                                 2~"Lab4b",
                                 3~"Probiotic",
                                 4~"Placebo"),
         gender_categorical=case_match(enr_Sex, 1~"male", 2~"female", 3~"undetermined"),
         conc=targetName_Value,
         anyDP=case_match(Arm,
                          1~"yes",
                          2~"yes",
                          3~"yes",
                          4~"no"))%>%
  select(id, timepoint, treatmentarm, gender_categorical, targetName, conc)%>%
  mutate(study="PROSYNK")


clean_data <- bind_rows(slim_micdrop_data, slim_prosynk_data)%>%
  filter(targetName %in% slim_micdrop_data$targetName & targetName %in% slim_prosynk_data$targetName)

# modelz ####
## change through time ####
library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)
library(purrr)
library(tidyr)

time_results <- clean_data %>%
  filter(timepoint %in% c("24 weeks", "52 weeks")) %>%
  group_by(targetName) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lme4::lmer(conc ~ timepoint * study * treatmentarm + gender_categorical + (1 | id), data = .x)),
    emm = map(model, ~ emmeans(
      .x,
      ~ timepoint * study * treatmentarm
    )),
    time_contrasts = map(emm, ~ contrast(
      .x,
      method = "pairwise",
      by = c("study", "treatmentarm"),
      adjust = "none"
    ) %>% as.data.frame()),
    treat_contrasts = map(emm, ~ contrast(
      .x,
      method = "pairwise",
      by = c("timepoint", "study"),
      adjust = "none"
    ) %>% as.data.frame()),
    study_contrasts = map(model, ~ emmeans(
      .x,
      ~ study | treatmentarm * timepoint
    ) %>%
      contrast(method = "pairwise") %>%
      as.data.frame())
  ) %>%
  select(targetName, time_contrasts, treat_contrasts, study_contrasts) %>%
  pivot_longer(
    cols = c(time_contrasts, treat_contrasts, study_contrasts),
    names_to = "contrast_type",
    values_to = "df"
  ) %>%
  unnest(df) %>%
  group_by(contrast_type) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

study_differences <- time_results%>%
  filter(p_adj<0.05, contrast_type=="study_contrasts")




top_study_differences <- study_differences%>%
  slice_min(n = 100, order_by = p_adj)%>%
  distinct(targetName)

clean_data%>%
  filter(timepoint %in%c("8 weeks", "24 weeks", "52 weeks"))%>%
  filter(targetName %in% top_study_differences$targetName[1:9])%>%
  mutate(treatmentarm=factor(treatmentarm, levels=c("DP 1 year", "Placebo", "Labinic", "Lab4b", "Probiotic")))%>%
  ggplot(., aes(x=factor(timepoint, levels = c("8 weeks", "24 weeks", "52 weeks")),
                y=conc,
                fill=study))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName)+
  scale_fill_viridis_d()+
  theme_minimal()





