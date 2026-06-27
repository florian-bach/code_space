library(tidyr)
library(purrr)
library(dplyr)
library(ggplot2)

`%notin%`=Negate(`%in%`)

clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks", "68 weeks")))%>%
  filter(targetName %notin% c("CTSS", "LTA|LTB", "IFNA2"))


# negative binomial means coefficient is on natural log scale.
#so coef of -0.3 -> each 1-unit increase in conc is associated with ~26.6% decrease in the expected number of infections
# higher conc at month 12 is associated with fewer cumulative infections over the first year.



# symp prob ####


## symp prob 12_24 ####
symp_prob_purf_12_24 <- clean_data%>%
  filter(mstatus==0, timepoint=="52 weeks")%>%
  mutate(symp_prob_12_24=total_n_malaria_12_24/total_n_para_12_24)%>%
  filter(!is.na(symp_prob_12_24), symp_prob_12_24<=1)%>%
  group_by(targetName)%>%
  nest()%>%
  # mutate(n_malaria_model=map(data, ~glm(symp_prob_12_24~conc+gender_categorical+log_qpcr+total_n_para_12, data=., family = "binomial" , weights = total_n_para_12_24))) %>%
  mutate(n_malaria_model=map(data, ~glm(symp_prob_12_24~conc+gender_categorical+treatmentarm, data=., family = "binomial" , weights = total_n_para_12_24))) %>%
  # mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(total_n_malaria_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  # mutate(n_para_model=map(data, ~glmmTMB::glmmTMB(total_n_para_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[14]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"))

symp_mala_12_24 <- symp_prob_purf_12_24%>%
  filter(n_malaria_padj<0.1)


(symp_prob_12_24_plot <- clean_data %>%
  filter(targetName %in% symp_mala_12_24$targetName,
         timepoint=="52 weeks",
         treatmentarm!="DP 2 years",
         !is.na(total_n_para_12_24),
         #log_qpcr<0.001
         )%>%
  mutate(pardens_dich = if_else(log_qpcr>0.001, "some", "none"))%>%
  mutate(symp_prob_12_24=total_n_malaria_12_24/total_n_para_12_24)%>%
  mutate(symp_prob_12_24f=case_when(symp_prob_12_24<0.25~"< 25%",
                                    symp_prob_12_24>=0.25&symp_prob_12_24<0.5~"25-50%",
                                    symp_prob_12_24>=0.5&symp_prob_12_24<0.75~"50-75%",
                                    symp_prob_12_24>=0.75~"> 75%",
                                    total_n_para_12_24==0~"no parasitemia"))%>%
  mutate(symp_prob_12_24f=factor(symp_prob_12_24f, levels=c("< 25%", "25-50%", "50-75%", "> 75%")))%>%
  # mutate(symp_prob_12_24f=case_when(symp_prob_12_24<0.33~"< 33%",
  #                                   symp_prob_12_24>=0.33&symp_prob_12_24<0.66~"33-66%",
  #                                   symp_prob_12_24>=0.66~"> 66%",
  #                                   total_n_para_12_24==0~"no parasitemia"))%>%
  # mutate(symp_prob_12_24f=factor(symp_prob_12_24f, levels=c("< 33%", "33-66%", "> 66%")))%>%
  
  # mutate(symp_prob_12_24f=case_when(symp_prob_12_24==0~"2% ",
  #                                   symp_prob_12_24>=0.25&symp_prob_12_24<0.5~"25%-50%",
  #                                   symp_prob_12_24>=0.5&symp_prob_12_24<0.75~"50%-75%",
  #                                   symp_prob_12_24>=0.75~"75% or more",
  #                                   total_n_para_12_24==0~"no parasitemia"))%>%
  filter(symp_prob_12_24f!="no parasitemia")%>%
  ggplot(aes(y=conc, x="", fill=symp_prob_12_24f))+
  # geom_violin()+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName*treatmentarm, scales = "free", ncol=2)+
  geom_smooth(method="gam")+
  # ggpubr::stat_compare_means(size=2, hide.ns = T, ref.group = "33% or less")+
  viridis::scale_fill_viridis(discrete = T)+
  ylab("concentration at 52 weeks")+
  xlab("probability of symptoms, given parasitemia in months 12-24")+
  theme_minimal()+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)
        ))

ggsave("~/Downloads/symp_prob_12_24_plot.png", symp_prob_12_24_plot, width=6, height = 6, bg="white", dpi=444)


n_malaria_12_24_plot <- clean_data %>%
  filter(targetName %in% symp_mala_12_24$targetName,
         timepoint=="52 weeks",
         treatmentarm!="DP 2 years",
         !is.na(total_n_malaria_12_24),
         #log_qpcr<0.001
  )%>%
  mutate(total_n_malaria_24f=case_when(total_n_malaria_12_24==0~"no malaria",
                                       total_n_malaria_12_24>=1&total_n_malaria_12_24<=2~"1-2",
                                       total_n_malaria_12_24>=3&total_n_malaria_12_24<=4~"3-4",
                                       total_n_malaria_12_24>=5~">5",
                                       ))%>%
  mutate(total_n_malaria_24f=factor(total_n_malaria_24f, levels = c("no malaria","1-2", "3-4", ">5")))%>%
  ggplot(aes(y=conc, x=factor(total_n_malaria_24f), fill=factor(total_n_malaria_24f)))+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  geom_boxplot(outliers = F)+
  facet_wrap(~treatmentarm+targetName, scales = "fixed", ncol=3)+
  # ggpubr::stat_compare_means(size=2, hide.ns = T, ref.group = "33% or less")+
  viridis::scale_fill_viridis(discrete = T)+
  ylab("concentration at 52 weeks")+
  xlab("number of malaria episodes months 12-24")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)
  )

ggsave("~/Downloads/symp_prob_12_24_plot2.png", n_malaria_12_24_plot, width=6, height = 6, bg="white", dpi=444)




clean_data %>%
  filter(targetName %in% symp_mala_12_24$targetName,
         timepoint=="52 weeks",
         treatmentarm!="DP 2 years",
         !is.na(total_n_para_12_24),
         log_qpcr<0.001)%>%
  mutate(pardens_dich = if_else(log_qpcr>0.001, "some", "none"))%>%
  mutate(symp_prob_12_24=total_n_malaria_12_24/total_n_para_12_24)%>%
  mutate(symp_prob_12_24f=case_when(symp_prob_12_24<0.25~"25% or less",
                                    symp_prob_12_24>=0.25&symp_prob_12_24<0.5~"25%-50%",
                                    symp_prob_12_24>=0.5&symp_prob_12_24<0.75~"50%-75%",
                                    symp_prob_12_24>=0.75~"75% or more",
                                    total_n_para_12_24==0~"symp_prob_12_24f"))%>%
  filter(symp_prob_12_24f!="symp_prob_12_24f")%>%
  ggplot(aes(y=conc, x=symp_prob_12_24))+
  geom_smooth(method="lm")+
  geom_point(shape=21, aes(fill=factor(total_n_malaria_12_24)))+
  facet_wrap(~targetName+treatmentarm, scales = "free")+
  ggpubr::stat_cor(method = "spearman")+
  viridis::scale_fill_viridis(discrete = T)+
  ylab("conc")+
  xlab("probability of symptoms, given infection")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_blank())





# sandbox LOOK! ####

future_purf_52 <- clean_data%>%
  filter(mstatus==0, timepoint=="52 weeks")%>%
  group_by(targetName)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~MASS::glm.nb(total_n_malaria_12_24~conc+gender_categorical+log_qpcr+total_n_para_12, data=.))) %>%
  mutate(n_para_model=map(data, ~MASS::glm.nb(total_n_para_12_24~conc+gender_categorical+log_qpcr+total_n_para_12, data=.))) %>%
  # mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(total_n_malaria_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  # mutate(n_para_model=map(data, ~glmmTMB::glmmTMB(total_n_para_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  mutate(n_para_model_summary=map(n_para_model, ~summary(.)))%>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[17]))%>%
  mutate(n_para_p=map_dbl(n_para_model_summary, ~coef(.)[17]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"),
         n_para_padj=p.adjust(n_para_p, method="BH"))

# negative binomial means coefficient is on natural log scale.
#so coef of -0.3 -> each 1-unit increase in conc is associated with ~26.6% decrease in the expected number of infections
# higher conc at month 12 is associated with fewer cumulative infections over the first year.


sig_mala_12 <- future_purf_52%>%
  filter(n_malaria_padj<0.1)

sig_para_12 <- future_purf_52%>%
  filter(n_para_padj<0.1)




clean_data %>%
  filter(targetName %in% unique(sig_mala_12$targetName),
         timepoint=="52 weeks",
         !is.na(total_n_para_12_18),
         log_qpcr<0.001,
         treatmentarm!="DP 2 years"
  )%>%
  mutate(pardens_dich = if_else(log_qpcr>0.001, "some", "none"))%>%
  mutate(targetNamef=factor(targetName, levels=sig_mala_12$targetName))%>%
  ggplot(., aes(x=factor(total_n_malaria_12_18), y=conc, fill=factor(total_n_malaria_12_18)))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetName+treatmentarm, scales = "free")+
  viridis::scale_fill_viridis(discrete = T)+
  xlab("number of malaria episodes in months 12-18")+
  ylab("concentration at 12 months")+
  theme_minimal()+
  theme(legend.position = "none")




clean_data %>%
  filter(#targetName %in% unique(symp_mala_12_18$targetName),
         targetName %in% c("FURIN",   "IL4R" ,   "IL6ST" ,  "SELE"),
         timepoint=="52 weeks",
         !is.na(total_n_para_12_18),
         log_qpcr<0.001,
         treatmentarm!="DP 2 years"
  )%>%
  mutate(pardens_dich = if_else(log_qpcr>0.001, "some", "none"))%>%
  mutate(targetNamef=factor(targetName, levels=sig_mala_12$targetName))%>%
  mutate(symp_prob_12_18=total_n_malaria_12_18/total_n_para_12_18)%>%
  mutate(symp_prob_12_18f=case_when(symp_prob_12_18<0.25~"25% or less",
                                    symp_prob_12_18>=0.25&symp_prob_12_18<0.5~"25%-50%",
                                    symp_prob_12_18>=0.5&symp_prob_12_18<0.75~"50%-75%",
                                    symp_prob_12_18>=0.75~"75% or more",
                                    total_n_para_12_18==0~"never parasitaemic"))%>%
  ggplot(., aes(x=symp_prob_12_18f, y=conc, fill=factor(symp_prob_12_18f)))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~treatmentarm+targetName, scales = "free")+
  viridis::scale_fill_viridis(discrete = T)+
  xlab("number of malaria episodes idvided by prevalent infectionsin months 12-18")+
  ylab("concentration at 12 months")+
  theme_minimal()+
  theme(legend.position = "none")





clean_data %>%
  filter(targetName %in% unique(sig_para_12$targetName),
         timepoint=="52 weeks",
         !is.na(total_n_para_12_18),
         log_qpcr<0.001
  )%>%
  mutate(pardens_dich = if_else(log_qpcr>0.001, "some", "none"))%>%
  mutate(targetNamef=factor(targetName, levels=sig_mala_12$targetName))%>%
  ggplot(aes(x=factor(total_n_para_12_18), y=conc, fill=factor(total_n_para_12_18)))+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_violin(outliers = FALSE)+
  facet_wrap(~targetNamef, scales = "free")+
  viridis::scale_fill_viridis(discrete = T)+
  xlab("number of parasitemic episodes in months 12-18")+
  ylab("concentration at 12 months")+
  facet_wrap(~targetName)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_blank())

# the relationship between tr1 freq and parasitemia at first malaria ####

# read tr1 predictions ####

tr1_freqs <- read.csv("~/Downloads/tr1_freq.csv")%>%
  filter(timepoint_num==52)

# make file containing info on the first malaria infection after 1st birthday
raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/MICDroP Data/MICDROP all visit database through February 28th 2026.dta")

first_malaria_after_one <- raw_data %>%
  mutate("flo_age_in_wks"=compage,
         "flo_age_in_months"=as.numeric(date-dob)%/%30.5)%>%
  filter(mstatus%in%c(1,2) & flo_age_in_wks>52 & flo_age_in_wks < 104)%>%
  group_by(id)%>%
  slice_min(flo_age_in_wks)%>%
  select(id, flo_age_in_wks, qPCRparsdens, pardens, temp)%>%
  inner_join(., tr1_freqs, by="id")

first_malaria_after_one%>%
  ggplot(., aes(x=predicted_tr1_freq, y=pardens))+
  geom_point()+
  scale_y_log10()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor()+
  theme_minimal()

first_malaria_after_one%>%
  ggplot(., aes(x=predicted_tr1_freq, y=temp))+
  geom_point()+
  # scale_y_log10()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor()+
  theme_minimal()

temp_model <- lm(temp~log10(pardens)+predicted_tr1_freq, data=first_malaria_after_one)

