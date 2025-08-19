library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

comp_pal <- c("no malaria"="lightgrey",
              "uncomplicated"="black",
              "complicated"="orange",
              "quinine for AL failure"="purple",
              "Q/AS failure"="purple")


mic_drop <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/MICDroP Data/MICDROP all visit database through June 30th 2025.dta")
mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")


sample_ages <- c(8, 24, 52, 68, 84, 104, 120, 136, 156)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1

sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))

kids_with_comp <- mic_drop %>%
  filter(mstatus==2)%>%
  group_by(id)%>%
  filter(!duplicated(date))


blood_counts <- mic_drop %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7,
         "flow_age_in_integer_years"=(floor(flo_age_in_wks/52))+1)%>%
  mutate("Timepoint_in_weeks"=if_else(
    flo_age_in_wks %in% sample_ages, flo_age_in_wks, ifelse(
      flo_age_in_wks %in% sample_ages_minus, flo_age_in_wks+1, if_else(
        flo_age_in_wks %in% sample_ages_plus, flo_age_in_wks-1, 999)))
    )%>%
  mutate(disease=case_when(mstatus==2 ~ "complicated",
                           mstatus==1 ~ "uncomplicated",
                           mstatus==0 ~ "no malaria"),
         ever_comp=if_else(id %in% kids_with_comp$id, "complicated", "never complicated"))%>%
  mutate(hb_less_than_eight =ifelse(hb<8, TRUE, FALSE))%>%
  pivot_longer(cols=c(plt, wbc, hb, eosino, mono, lymph, neutro), names_to = "cell_type", values_to = "cell_freq")%>%
  filter(!is.na(cell_freq))%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(id, mic_drop_key$id)])%>%
  mutate(treatmentarm=case_match(treatmentarm,
                                 1~"No DP",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))
  

  

 
blood_counts%>%
  filter(cell_type%in%c("mono", "plt", "hb"), Timepoint_in_weeks<53, mstatus==0)%>%
  mutate(anyDP=if_else(treatmentarm=="No DP", "no DP", "DP"))%>%
  ggplot(., aes(x=factor(Timepoint_in_weeks), y=cell_freq,  fill = anyDP))+
    geom_boxplot(outliers = F)+
    # geom_point(aes(alpha=factor(ever_comp), colour = disease))+
    # geom_boxplot(aes(fill=factor(Timepoint_in_weeks)), outliers = FALSE)+
    # geom_line(aes(group=id))+
    facet_wrap(~cell_type, scales="free")+
    scale_y_log10()+
  ggpubr::stat_compare_means(size=3, label = "p.format")+
    # scale_fill_manual(values=colorspace::sequential_hcl(n = 9, palette="LaJolla"))+
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
  ggplot(., aes(x=factor(mstatus), y=cell_freq))+
  # geom_point(aes(color=factor(disease)))+
  geom_boxplot(aes(fill=factor(mstatus)), outliers = FALSE)+
  # geom_line(aes(group=id))+
  # facet_wrap(~cell_type, scales="free")+
  # scale_x_log10()+
  scale_fill_manual(values=colorspace::sequential_hcl(n = 7, palette="LaJolla"))+
  scale_color_manual(values=comp_pal)+
  theme_minimal()



blood_counts %>%
  filter(cell_type=="hb")%>%
  filter(cell_freq<100)%>%
  arrange(mstatus)%>%
  ggplot(., aes(x=flo_age_in_wks, y=cell_freq))+
  geom_point(aes(color=factor(disease)))+
  # geom_boxplot(aes(fill=factor(mstatus)), outliers = FALSE)+
  # geom_line(aes(group=id))+
  facet_wrap(~treatmentarm, scales="fixed")+
  # scale_x_log10()+
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
             
