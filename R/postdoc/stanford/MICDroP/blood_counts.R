library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

comp_pal <- c("no malaria"="lightgrey",
              "uncomplicated"="black",
              "complicated"="orange",
              "quinine for AL failure"="purple",
              "Q/AS failure"="purple")


mic_drop <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2024_07/MICDROP expanded database through July 31st 2024.dta")


sample_ages <- c(8, 24, 52, 68, 84, 104, 120)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1

sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))

kids_with_comp <- mic_drop %>%
  filter(mstatus==2)%>%
  group_by(id)%>%
  filter(!duplicated(date))


blood_counts <- mic_drop %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7)%>%
  mutate("Timepoint_in_weeks"=if_else(
    flo_age_in_wks %in% sample_ages, flo_age_in_wks, ifelse(
      flo_age_in_wks %in% sample_ages_minus, flo_age_in_wks+1, if_else(
        flo_age_in_wks %in% sample_ages_plus, flo_age_in_wks-1, 999)))
    )%>%
  mutate(disease=case_when(mstatus==2 ~ "complicated",
                           mstatus==1 ~ "uncomplicated",
                           mstatus==0 ~ "no malaria"),
         ever_comp=if_else(id %in% kids_with_comp$id, "complicated", "never complicated"))%>%
  pivot_longer(cols=c(plt, wbc, hb, eosino, mono, lymph, neutro), names_to = "cell_type", values_to = "cell_freq")%>%
  filter(!is.na(cell_freq), Timepoint_in_weeks %in% sample_ages)
  

 
blood_counts%>%
  filter(cell_type!="hb")%>%
  arrange(ever_comp)%>%
  group_by(id) %>%
  add_count(name="total_n_infection") %>%
  arrange(AGE) %>%
  mutate(n_infection = seq(1, max(total_n_infection)))%>%
  ggplot(., aes(x=factor(Timepoint_in_weeks), y=cell_freq))+
    geom_violin(aes(fill=factor(Timepoint_in_weeks)), draw_quantiles = c(0.25, 0.5, 0.75))+
    # geom_point(aes(alpha=factor(ever_comp), colour = disease))+
    # geom_boxplot(aes(fill=factor(Timepoint_in_weeks)), outliers = FALSE)+
    # geom_line(aes(group=id))+
    facet_wrap(~cell_type, scales="free")+
    scale_y_log10()+
    scale_fill_manual(values=colorspace::sequential_hcl(n = 7, palette="LaJolla"))+
    scale_alpha_manual(values = c(1,0))+
    scale_color_manual(values=comp_pal)+
    theme_minimal()



blood_counts%>%
  filter(cell_type!="hb")%>%
  arrange(disease)%>%
  ggplot(., aes(x=qPCRparsdens+0.01, y=cell_freq))+
  geom_point(aes(color=factor(disease)))+
  # geom_line(aes(group=id))+
  facet_wrap(~cell_type, scales="free")+
  scale_y_log10()+
  scale_x_log10()+
  scale_fill_manual(values=colorspace::sequential_hcl(n = 7, palette="LaJolla"))+
  scale_color_manual(values=comp_pal)+
  theme_minimal()


blood_counts %>%
  filter(cell_type=="hb")%>%
  arrange(disease)%>%
  ggplot(., aes(x=qPCRparsdens+0.01, y=cell_freq))+
  geom_point(aes(color=factor(disease)))+
  geom_boxplot(aes(fill=factor(Timepoint_in_weeks)), outliers = FALSE)+
  # geom_line(aes(group=id))+
  # facet_wrap(~cell_type, scales="free")+
  scale_x_log10()+
  scale_fill_manual(values=colorspace::sequential_hcl(n = 7, palette="LaJolla"))+
  scale_color_manual(values=comp_pal)+
  theme_minimal()

blood_count_qpcr_purrr <- blood_counts %>%
  mutate(log_pd=log10(qPCRparsdens+0.01), 
         log_cell_freq=log10(cell_freq+0.01),
         subject_id=factor(id))%>%
  filter(cell_type!="hb", !is.na(log_cell_freq), !is.na(log_pd))%>%
  group_by(cell_type)%>%
  nest()%>%
  mutate(model= map(data, ~lme4::lmer(log_pd~log_cell_freq+(1|subject_id), data=.)))%>%
  mutate(summary=map(model, ~summary(.))) 


blood_count_age_purrr <- blood_counts %>%
  mutate(log_pd=log10(qPCRparsdens+0.01), 
         log_cell_freq=log10(cell_freq+0.01),
         subject_id=factor(id))%>%
  filter(cell_type!="hb", !is.na(log_cell_freq), !is.na(log_pd))%>%
  group_by(cell_type)%>%
  nest()%>%
  mutate(model= map(data, ~lme4::lmer(log_cell_freq~Timepoint_in_weeks+(1|subject_id), data=.)))%>%
  mutate(summary=map(model, ~summary(.))) 
             
