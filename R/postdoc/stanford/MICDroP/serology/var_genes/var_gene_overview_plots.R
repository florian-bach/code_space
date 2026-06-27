# preamble ####
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(emmeans)

treatment_palette <-  c("darkred", "darkblue")
names(treatment_palette) <-c("Placebo", "DP 1 year")

var_luminex <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/slim_luminex_results.csv", check.names = F)
colnames(var_luminex)[12:13] <- c("CIDRg3",	"CIDRa6")
colnames(var_luminex)[31:32] <- c("schizont",	"tetanus")

plate_map <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/final_plate_map.csv")
mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")
mic_drop_hbs <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/MICDROP SickleTr final.dta")
maternal_treatment_arms <- haven::read_dta("~/Library/CloudStorage/Box-Box/DP+SP study/Databases and preliminary findings/Final database used for analyses/DPSP treatment allocation_FINAL.dta")

raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/Oct25/MICDSpecimenBoxOct25_withclinical.dta")

plate_map <- plate_map%>%
  mutate(well=paste(plate.row, plate.column, sep=""),
         location=paste(dane.plate, well, sep=""))

slim_plate_map <- plate_map%>%
  select(location, id, date, ageinwks)%>%
  mutate("date"=as.Date(lubridate::parse_date_time(.$date, orders="%m/%d/%y")))

all_samples <- distinct(plate_map, id, date, ageinwks)%>%
  mutate("date"=as.Date(lubridate::parse_date_time(.$date, orders="%m/%d/%y")))

# add epi data ####
## add visit-level metadata ####
metadata_columns <- c("id", "date", "mstatus", "pardens", "qPCRparsdens", "fever", "febrile", "gender", "mom_rx")

dobs <- raw_data%>%
  distinct(id, dob)%>%
  filter(!is.na(dob))

visit_metadata <- raw_data%>%
  mutate(date=as.Date(lubridate::parse_date_time(date, "%y/%m/%d")))%>%
  select(all_of(metadata_columns))%>%
  right_join(., all_samples, by=c("id", "date"))%>%
  mutate(dob2=dobs$dob[match(id, dobs$id)])%>%
  mutate(real_dob=as.Date(lubridate::parse_date_time(dob2, "%y/%m/%d")))%>%
  # mutate(ageinwks=round((date-real_dob)/7))%>%
  mutate(timepoint_num=as.numeric(ageinwks))%>%
  mutate(timepoint_num=case_when(timepoint_num==9~8,
                                 timepoint_num==53~52,
                                 timepoint_num==105~104,
                                 .default=timepoint_num))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.001))%>%
  mutate(gender_categorical=if_else(gender==1, "male", "female"))


## static meta data (infection, birth outcomes) ####
static_metadata <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/micdrop_metadata_for_immunology.csv")

long_luminex <- var_luminex%>%
  mutate(location=paste(Plate, Well, sep=""))%>%
  pivot_longer(cols=c(colnames(var_luminex)[3:32]), names_to = "antigen", values_to = "MFI")%>%
  left_join(., slim_plate_map, by="location")%>%
  left_join(., visit_metadata, by=c("id", "date"))%>%
  mutate(date=as.Date(lubridate::parse_date_time(date, "%y/%m/%d")))%>%
  mutate(log_mfi=log10(MFI+1),
         timepointw=paste(timepoint_num, "weeks"),
         timepoint=factor(timepointw, levels=c("8 weeks", "24 weeks", "52 weeks", "104 weeks")))%>%
  filter(!is.na(log_mfi), is.finite(log_mfi))%>%
  left_join(., static_metadata, by=c("id"))%>%
  mutate(antigen=gsub("_", " ", antigen),
         id_cat=as.character(id))%>%
  # there's a couple of mispipetted samples that have NAs for timepoint
  filter(!is.na(timepoint_num))%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))

## add binder information to same table ####
var_domains <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/var_domain_dictionary.csv")
colnames(var_domains)[1] <- "antigen"
slim_var_domains <- var_domains[,1:3]

long_luminex <- long_luminex%>%
  left_join(., slim_var_domains, by="antigen")

write.csv(long_luminex, "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/long_luminex.csv", row.names = F)

# overview plot ####
large_plot <- long_luminex %>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=timepoint, y=MFI, fill=timepoint))+
  geom_boxplot(outliers = F, color="grey")+
  scale_y_log10()+
  facet_wrap(~antigen, scales="free", labeller = label_wrap_gen(width = 10))+
  viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  theme_minimal()+
  theme(legend.position="bottom",
        legend.direction = "horizontal",
        axis.text.x = element_blank(),
        axis.title.x=element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/large_overview_plot.png", large_plot, height=9, width=15, dpi=400, bg="white")  

large_plot2 <- long_luminex %>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=timepoint, y=MFI, fill=treatmentarm))+
  geom_boxplot(outliers = F, color="grey")+
  scale_y_log10()+
  facet_wrap(~antigen, scales="free", labeller = label_wrap_gen(width = 10))+
  scale_fill_manual(values=c("darkred", "#00555A"))+
  theme_minimal()+
  theme(legend.position="bottom",
        legend.direction = "horizontal",
        axis.text.x = element_blank(),
        axis.title.x=element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/large_overview_plot2.png", large_plot2, height=9, width=15, dpi=400, bg="white")  



 # linear regression ####
## kids treatment ####
treatment_purf <- long_luminex%>%
  filter(MFI>100)%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(log_mfi~timepoint*treatmentarm+log_qpcr+(1|id_cat), data=.))) %>%
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
  mutate("104 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[4])) %>%
  mutate("DP boost 52 - 104 weeks"=map_dbl(emm_contrast_summary2, ~.$p.value[6])) %>%
  mutate("Placebo boost 52 - 104 weeks"=map_dbl(emm_contrast_summary2, ~.$p.value[12])) %>%
  pivot_longer(cols=ends_with("weeks"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="BH"))


#nothing
sigs_8 <- treatment_purf%>%
  filter(padj<0.1, contrast=="8 weeks")%>%
  select(antigen, contrast, p, padj)
#nothing
sigs_24 <- treatment_purf%>%
  filter(padj<0.1, contrast=="24 weeks")%>%
  select(antigen, contrast, p, padj)

# 10 without pcr; 15 with pcr
sigs_52 <- treatment_purf%>%
  filter(padj<0.01, contrast=="52 weeks")%>%
  select(antigen, contrast, p, padj)%>%
  arrange(padj)

# nothing
sigs_104 <- treatment_purf%>%
  filter(padj<0.1, contrast=="104 weeks")%>%
  select(antigen, contrast, p, padj)

#6 without pcr, 5 with
sigs_52_104 <- treatment_purf%>%
  filter(padj<0.1, contrast%in%c("Placebo boost 52 - 104 weeks", "DP boost 52 - 104 weeks"))%>%
  select(antigen, contrast, p, padj)





treatment_purf%>%
  select(antigen, contrast, p, padj)%>%
  write.csv(., "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/var_gene_treatment_regression.csv")


sig_52_plot1 <- long_luminex %>%
  # mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  filter(antigen %in% sigs_52$antigen[1:8])%>%
  mutate(antigen=factor(antigen, levels=sigs_52$antigen[1:8]))%>%
  ggplot(., aes(x=timepoint, y=MFI, fill=treatmentarm))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  facet_wrap(~antigen, scales="free", labeller = label_wrap_gen(width = 10), nrow=2)+
  scale_fill_manual(values=treatment_palette)+
  theme_minimal()+
  theme(legend.position="bottom",
        legend.direction = "horizontal",
        axis.text.x = element_blank(),
        axis.title.x=element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/sig_52_overview_plot1.png", sig_52_plot1, height=4, width=8, dpi=400, bg="white")  


sig_52_plot2 <- long_luminex %>%
  # mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  filter(antigen %in% sigs_52$antigen[9:16])%>%
  mutate(antigen=factor(antigen, levels=sigs_52$antigen[9:16]))%>%
  ggplot(., aes(x=timepoint, y=MFI, fill=treatmentarm))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  facet_wrap(~antigen, scales="free", labeller = label_wrap_gen(width = 10), nrow=2)+
  scale_fill_manual(values=c("darkred", "#00555A"))+
  theme_minimal()+
  theme(legend.position="bottom",
        legend.direction = "horizontal",
        axis.text.x = element_blank(),
        axis.title.x=element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/sig_52_overview_plot2.png", sig_52_plot2, height=4, width=8, dpi=400, bg="white")  

sig_104_plot <- long_luminex %>%
  # mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  filter(antigen %in% sigs_104$antigen)%>%
  ggplot(., aes(x=timepoint, y=MFI, fill=treatmentarm))+
  geom_boxplot(outliers = F)+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  scale_y_log10()+
  facet_wrap(~antigen, scales="free", labeller = label_wrap_gen(width = 10))+
  scale_fill_manual(values=c("darkred", "#00555A"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/sig_104_overview_plot.png", sig_104_plot, height=3, width=6, dpi=400, bg="white")  



sig_52_104_plot <- long_luminex %>%
  filter(MFI>100)%>%
  filter(antigen %in% sigs_52_104$antigen)%>%
  ggplot(., aes(x=timepoint, y=MFI, fill=treatmentarm))+
  geom_boxplot(outliers = F)+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  scale_y_log10()+
  facet_grid(treatmentarm~antigen, scales="free", labeller = label_wrap_gen(width = 10))+
  scale_fill_manual(values=treatment_palette)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        axis.title.x = element_blank(),
        legend.position = "none"
          )

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/sig_52_104_overview_plot.png", sig_52_104_plot, height=5, width=8, dpi=400, bg="white")  


## mom treatment ####
mom_rx_purf <- long_luminex%>%
  group_by(antigen)%>%
  # filter(MFI>100)%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(log_mfi~timepoint*mom_rx+(1|id_cat), data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ mom_rx | timepoint)))%>%
  mutate(emm2=map(time_model, ~emmeans(., specs = pairwise ~ timepoint | mom_rx)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  mutate("8 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("24 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[5])) %>%
  mutate("52 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[8])) %>%
  mutate("104 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[11])) %>%
  pivot_longer(cols=ends_with("weeks"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

# DP vs DPSP: var2csa @ 8 weeks
# DP vs SP: var2csa @ 8 weeks
# DPSP vs SP: no diff
sigs <- mom_rx_purf%>%
  filter(padj<0.1)%>%
  select(antigen, contrast, p, padj)


var2csa_time <- long_luminex %>%
  filter(antigen=="var2csa")%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=timepoint, y=MFI, fill=mom_rx))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  scale_fill_manual(values=c("darkred", "#00555A", "#FEC111"))+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/var2csa_time.png", var2csa_time, height=3, width=6, dpi=400, bg="white")  



var2csa_gravid <- long_luminex %>%
  mutate(gravidcat=case_when(gravidcat==1~"Primigravid",
                             gravidcat==2~"Secundigravid",
                             gravidcat==3~"Multigravid"),
         gravidcat=factor(gravidcat, levels=c("Primigravid", "Secundigravid", "Multigravid")))%>%
  filter(antigen=="var2csa", timepoint=="8 weeks")%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=factor(gravidcat), y=MFI, fill=factor(mom_rx)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  scale_fill_manual(values=c("darkred", "#00555A", "#FEC111"))+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/var2csa_gravid.png", var2csa_gravid, height=3, width=6, dpi=400, bg="white")  

## mom gravidity ####

mom_rx_gravid_purf <- long_luminex%>%
  mutate(gravidcat=case_when(gravidcat==1~"Primigravid",
                             gravidcat==2~"Secundigravid",
                             gravidcat==3~"Multigravid"),
         gravidcat=factor(gravidcat, levels=c("Primigravid", "Secundigravid", "Multigravid")))%>%
  filter(timepoint=="8 weeks")%>%
  group_by(antigen)%>%
  # filter(MFI>100)%>%
  nest()%>%
  mutate(time_model=map(data, ~lm(log_mfi~mom_rx*gravidcat, data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ mom_rx | gravidcat)))%>%
  mutate(emm2=map(time_model, ~emmeans(., specs = pairwise ~ gravidcat | mom_rx)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  #DP vs SP comparison only
  mutate("primigravid"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("secundigravid"=map_dbl(emm_contrast_summary, ~.$p.value[5])) %>%
  mutate("multiigravid"=map_dbl(emm_contrast_summary, ~.$p.value[8])) %>%
  pivot_longer(cols=ends_with("gravid"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

gravid_rx_sigs <- mom_rx_gravid_purf%>%
  filter(padj<0.1)%>%
  select(antigen, contrast, padj)


sigs_gravid <- long_luminex %>%
  mutate(gravidcat=case_when(gravidcat==1~"Primigravid",
                             gravidcat==2~"Secundigravid",
                             gravidcat==3~"Multigravid"),
         gravidcat=factor(gravidcat, levels=c("Primigravid", "Secundigravid", "Multigravid")))%>%
  filter(antigen%in%c(gravid_rx_sigs$antigen), timepoint=="8 weeks")%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=factor(gravidcat), y=MFI, fill=factor(mom_rx)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  facet_wrap(~antigen)+
  scale_fill_manual(values=c("darkred", "#00555A", "#FEC111"))+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/sigs_rx_gravid.png", sigs_gravid, height=4, width=6, dpi=400, bg="white")  

## placental malaria####

any_hp_8_weeks <- long_luminex %>%
  mutate(gravidcat=case_when(gravidcat==1~"Primigravid",
                             gravidcat==2~"Secundigravid",
                             gravidcat==3~"Multigravid"),
         gravidcat=factor(gravidcat, levels=c("Primigravid", "Secundigravid", "Multigravid")))%>%
  filter(timepoint=="8 weeks", !is.na(anyHP))%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=factor(anyHP), y=MFI, fill=factor(anyHP)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ggpubr::stat_compare_means(size=2.2)+
  facet_wrap(~antigen)+
  scale_fill_manual(values=c("darkred", "#00555A", "#FEC111"))+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/any_hp_8_weeks.png", any_hp_8_weeks, height=9, width=14, dpi=400, bg="white")  

####
## infection incidence ####
### 52 weeks ####
inf_purf_52 <- long_luminex%>%
  filter(timepoint=="52 weeks", MFI>100)%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~MASS::glm.nb(total_n_malaria_12~log_mfi, data=.))) %>%
  mutate(n_para_model=map(data, ~MASS::glm.nb(total_n_para_12~log_mfi, data=.))) %>%
  # mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(total_n_malaria_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  # mutate(n_para_model=map(data, ~glmmTMB::glmmTMB(total_n_para_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  mutate(n_para_model_summary=map(n_para_model, ~summary(.)))%>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[8]))%>%
  mutate(n_para_p=map_dbl(n_para_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"),
         n_para_padj=p.adjust(n_para_p, method="BH"))

inf_sigs_52 <- inf_purf_52%>%
  filter(n_para_padj<0.05 |n_malaria_padj <0.05 )%>%
  select(antigen,n_para_padj, n_malaria_padj )%>%
  arrange(n_para_padj)

inf_purf_52%>%
  select(antigen,n_para_padj, n_malaria_padj )%>%
  write.csv(., "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/inf_purf_52.csv")


n_para_sigs_52_1 <- long_luminex %>%
  filter(antigen %in% inf_sigs_52$antigen[1:9])%>%
  mutate(antigen=factor(antigen, levels=inf_sigs_52$antigen[1:9]))%>%
  filter(timepoint =="52 weeks", treatmentarm=="Placebo")%>%
  mutate(total_n_para_12f=if_else(total_n_para_12>6, "7+", as.character(total_n_para_12)))%>%
  ggplot(., aes(x=factor(total_n_para_12f), y=MFI, fill=factor(total_n_para_12f)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ylab("MFI at 52 weeks")+
  facet_wrap(~antigen, scales="free", nrow=3)+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/n_para_sigs_52_1.png", n_para_sigs_52_1, height=8, width=8, dpi=400, bg="white")  

n_para_sigs_52_2 <- long_luminex %>%
  filter(antigen %in% inf_sigs_52$antigen[10:18])%>%
  mutate(antigen=factor(antigen, levels=inf_sigs_52$antigen[10:18]))%>%
  filter(timepoint =="52 weeks", treatmentarm=="Placebo")%>%
  mutate(total_n_para_12f=if_else(total_n_para_12>6, "7+", as.character(total_n_para_12)))%>%
  ggplot(., aes(x=factor(total_n_para_12f), y=MFI, fill=factor(total_n_para_12f)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ylab("MFI at 52 weeks")+
  facet_wrap(~antigen, scales="free", nrow=3)+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_blank())
ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/n_para_sigs_52_2.png", n_para_sigs_52_2, height=8, width=8, dpi=400, bg="white")  

n_malaria_sigs_52 <- long_luminex %>%
  filter(antigen %in% inf_sigs_52$antigen[inf_sigs_52$n_malaria_padj<0.05])%>%
  mutate(antigen=factor(antigen, levels=inf_sigs_52$antigen[inf_sigs_52$n_malaria_padj<0.05]))%>%
  filter(timepoint =="52 weeks", treatmentarm=="Placebo")%>%
  ggplot(., aes(x=factor(total_n_malaria_12), y=MFI, fill=factor(total_n_malaria_12)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ylab("MFI at 52 weeks")+
  facet_wrap(~antigen, scales="free", nrow=2)+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/n_malaria_sigs_52.png", n_malaria_sigs_52, height=4, width=8, dpi=400, bg="white")  

### 104 weeks ####

inf_purf_104 <- long_luminex%>%
  filter(timepoint=="104 weeks", MFI>100)%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~MASS::glm.nb(total_n_malaria_12_24~log_mfi*treatmentarm, data=.))) %>%
  mutate(n_para_model=map(data, ~MASS::glm.nb(total_n_para_12_24~log_mfi*treatmentarm, data=.))) %>%
  mutate(emm=map(n_malaria_model, ~emmeans(., specs = pairwise ~ treatmentarm)))%>%
  mutate(emm2=map(n_para_model, ~emmeans(., specs = pairwise ~ treatmentarm)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  mutate("8 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("24 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("52 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  mutate("104 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[4])) %>%
  
  # mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(total_n_malaria_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  # mutate(n_para_model=map(data, ~glmmTMB::glmmTMB(total_n_para_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  mutate(n_para_model_summary=map(n_para_model, ~summary(.)))%>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[11]))%>%
  mutate(n_para_p=map_dbl(n_para_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"),
         n_para_padj=p.adjust(n_para_p, method="BH"))

inf_sigs_104 <- inf_purf_104%>%
  filter(n_malaria_padj<0.05|
         n_para_padj<0.05)%>%
  select(antigen,n_para_padj, n_malaria_padj )%>%
  arrange(n_para_padj)


n_malaria_sigs_104 <- long_luminex %>%
  filter(antigen %in% inf_sigs_104$antigen[1:9])%>%
  mutate(antigen=factor(antigen, levels=inf_sigs_104$antigen[1:9]))%>%
  filter(timepoint =="104 weeks")%>%
  mutate(total_n_para_24f=if_else(total_n_para_24>6, "7+", as.character(total_n_para_24)))%>%
  ggplot(., aes(x=factor(total_n_para_24f), y=MFI, fill=factor(total_n_para_24f)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ylab("MFI at 104 weeks")+
  facet_wrap(~antigen, scales="free", nrow=3)+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/n_malaria_sigs_104.png", n_malaria_sigs_104, height=8, width=8, dpi=400, bg="white")  



n_malaria_sigs_104_dp <- long_luminex %>%
  filter(antigen %in% inf_sigs_52$antigen[1:9])%>%
  mutate(antigen=factor(antigen, levels=inf_sigs_52$antigen[1:9]))%>%
  filter(timepoint =="104 weeks", treatmentarm=="DP 1 year")%>%
  mutate(total_n_para_24f=if_else(total_n_para_24>6, "7+", as.character(total_n_para_24)))%>%
  ggplot(., aes(x=total_n_para_24f, y=MFI, fill=total_n_para_24f))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ylab("MFI at 104 weeks (DP)")+
  facet_wrap(~antigen, scales="free", nrow=3)+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/n_malaria_sigs_104_dp.png", n_malaria_sigs_104_dp, height=8, width=8, dpi=400, bg="white")  


### 52 - 104 weeks ####
inf_purf_52_104 <- long_luminex%>%
  filter(timepoint=="52 weeks", MFI>100)%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~lm(log_mfi~total_n_malaria_12_24+log_qpcr, data=.))) %>%
  mutate(n_para_model=map(data, ~lm(log_mfi~total_n_para_12_24+log_qpcr, data=.))) %>%
  # mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(total_n_malaria_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  # mutate(n_para_model=map(data, ~glmmTMB::glmmTMB(total_n_para_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  mutate(n_para_model_summary=map(n_para_model, ~summary(.)))%>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[11]))%>%
  mutate(n_para_p=map_dbl(n_para_model_summary, ~coef(.)[11]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"),
         n_para_padj=p.adjust(n_para_p, method="BH"))

inf_sigs_52_104 <- inf_purf_52_104%>%
  filter(n_malaria_p<0.05,
         n_para_p<0.05)


n_malaria_sigs_52_104 <- long_luminex %>%
  filter(!is.na(total_n_malaria_12_24))%>%
  filter(antigen %in% c(inf_sigs_52_104$antigen[1:3], "schizont"))%>%
  # mutate(antigen=factor(antigen, levels=inf_sigs_52_104$antigen[1:9]))%>%
  filter(timepoint %in% c("52 weeks"))%>%
  # mutate(total_n_malaria_12_24f=if_else(total_n_malaria_12_24>3, "4+", as.character(total_n_malaria_12_24)))%>%
  ggplot(., aes(x=factor(total_n_malaria_12_24), y=MFI, fill=factor(total_n_malaria_12_24)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ylab("MFI at 52 weeks")+
  xlab("number malaria episodes in second year of life")+
  facet_wrap(~antigen+treatmentarm, scales="fixed")+
  scale_fill_manual(values=viridis::magma(n=7))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/n_malaria_sigs_52_104.png", n_malaria_sigs_52_104, height=8, width=8, dpi=400, bg="white")  


n_para_non_sigs_52_104 <- long_luminex %>%
  filter(!is.na(total_n_malaria_12_24))%>%
  filter(antigen %in% c("CIDRa6",
                        "CIDRa1.7 2083-1",
                        "CIDRg3 IT4var08"
  ))%>%
  # mutate(antigen=factor(antigen, levels=inf_sigs_52_104$antigen[1:9]))%>%
  filter(timepoint %in% c("52 weeks", "104 weeks"))%>%
  # mutate(total_n_malaria_12_24f=if_else(total_n_malaria_12_24>3, "4+", as.character(total_n_malaria_12_24)))%>%
  ggplot(., aes(x=factor(total_n_malaria_12_24), y=MFI, fill=factor(total_n_malaria_12_24)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ylab("MFI at 52 weeks")+
  xlab("number malaria episodes in second year of life")+
  facet_wrap(~antigen+treatmentarm, scales="fixed")+
  scale_fill_manual(values=viridis::magma(n=7))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/n_para_non_sigs_52_104.png", n_para_non_sigs_52_104, height=8, width=12, dpi=400, bg="white")  






n_malaria_sigs_52_104_DP <- long_luminex %>%
  filter(antigen %in% inf_sigs_52$antigen[1:9])%>%
  mutate(antigen=factor(antigen, levels=inf_sigs_52$antigen[1:9]))%>%
  filter(timepoint =="104 weeks", treatmentarm=="DP 1 year")%>%
  mutate(total_n_para_12_24f=if_else(total_n_para_12_24>3, "4+", as.character(total_n_para_12_24)))%>%
  ggplot(., aes(x=factor(total_n_para_12_24f), y=MFI, fill=factor(total_n_para_12_24f)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  facet_wrap(~antigen, scales="free", nrow=3)+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/n_malaria_sigs_52_104_DP.png", n_malaria_sigs_52_104_DP, height=8, width=8, dpi=400, bg="white")  

### 8 weeks ####
inf_purf_8 <- long_luminex%>%
  filter(timepoint=="8 weeks", MFI>100)%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~glm(any_malar_6~log_mfi+gender_categorical+log_qpcr, family = "binomial", data=.))) %>%
  mutate(n_para_model=map(data, ~glm(any_para_6~log_mfi+gender_categorical+log_qpcr,family = "binomial",  data=.))) %>%
  # mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(total_n_malaria_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  # mutate(n_para_model=map(data, ~glmmTMB::glmmTMB(total_n_para_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  mutate(n_para_model_summary=map(n_para_model, ~summary(.)))%>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[14]))%>%
  mutate(n_para_p=map_dbl(n_para_model_summary, ~coef(.)[14]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"),
         n_para_padj=p.adjust(n_para_p, method="BH"))


inf_sigs_8 <- inf_purf_8%>%
  filter(n_para_padj<0.1 |n_malaria_padj <0.1 )%>%
  select(antigen,n_para_padj, n_malaria_padj )%>%
  arrange(n_para_padj)

inf_purf_8%>%
  select(antigen,n_para_padj, n_malaria_padj )%>%
  write.csv(., "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/inf_purf_8.csv")


n_malaria_sigs_8a <- long_luminex %>%
  filter(antigen %in% c("CIDRa5 IT4var14"))%>%
  filter(timepoint =="8 weeks", !is.na(total_n_para_6))%>%
  ggplot(., aes(x=factor(total_n_para_6), y=MFI, fill=factor(total_n_para_6)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ylab("MFI at 8 weeks")+
  xlab("number of parasitemic months in the first six months of life")+
  facet_wrap(~antigen, scales="free", nrow=3)+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/var_para8_plot.png", n_malaria_sigs_8a, height=5, width=8, dpi=400, bg="white")  

n_malaria_sigs_8b <- long_luminex %>%
  filter(antigen %in% c("CIDRa6",
                        "CIDRa1.7 2083-1",
                        "CIDRg3 IT4var08",
                        "schizont"))%>%
  filter(timepoint =="8 weeks", !is.na(total_n_malaria_6))%>%
  ggplot(., aes(x=factor(total_n_malaria_6), y=MFI, fill=factor(total_n_malaria_6)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ylab("MFI at 8 weeks")+
  xlab("number of malaria episodesin the first six months of life")+
  facet_wrap(~antigen, scales="free", nrow=3)+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(legend.position = "none")
ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/var_malaria8_plot.png", n_malaria_sigs_8b, height=5, width=8, dpi=400, bg="white")  

#### parasitemia ####

var_parasitemia_plot <- long_luminex %>%
  filter(antigen %in% c("CIDRa6",
                        "CIDRa1.7 2083-1",
                        "CIDRg3 IT4var08",
                        "schizont"))%>%
  ggplot(., aes(x=log_qpcr, y=log_mfi, color=antigen))+
  geom_smooth(method="lm")+
  geom_point()+
  ylab("log10 MFI")+
  xlab("log10 qPCR")+
  ggpubr::stat_cor(method = "spearman", label.y = c(0.8, 0.6, 0.4, 0.2))+
  facet_wrap(~timepoint, scales="fixed", nrow=1)+
  scale_color_manual(values=viridis::magma(n=7))+
  theme_minimal()+
  guides(color = guide_legend(override.aes = list(label = "")))+
  theme(legend.position = "bottom",
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key = element_rect(fill = "white", colour = NA),        axis.text = element_text(size=14),
        strip.text = element_text(size=16),
        axis.title = element_text(size=16))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/var_parasitemia_plot.png", var_parasitemia_plot, height=5, width=8, dpi=400, bg="white")  

# probability of symptoms, given infection ####

symp_prob_purf_52 <- long_luminex%>%
  filter(timepoint=="52 weeks")%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~glm(symp_prop_12_24~log_mfi+gender_categorical+treatmentarm, data=., family = "binomial" , weights = total_n_para_12_24))) %>%
  # mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(total_n_malaria_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  # mutate(n_para_model=map(data, ~glmmTMB::glmmTMB(total_n_para_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[14]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"))

symp_prob_sigs_52 <- symp_prob_purf_52%>%
  filter(n_malaria_padj < 0.1 )%>%
  select(antigen,n_malaria_p, n_malaria_padj )%>%
  arrange(n_malaria_padj)


long_luminex %>%
  mutate(symp_prob_12_24f=case_when(symp_prop_12_24<0.25~"< 25%",
                                    symp_prop_12_24>=0.25&symp_prop_12_24<0.5~"25-50%",
                                    symp_prop_12_24>=0.5&symp_prop_12_24<0.75~"50-75%",
                                    symp_prop_12_24>=0.75~"> 75%",
                                    total_n_para_12_24==0~"no parasitemia"))%>%
  mutate(symp_prob_12_24f=factor(symp_prob_12_24f, levels=c("< 25%", "25-50%", "50-75%", "> 75%")))%>%
  filter(antigen %in% symp_prob_sigs_52$antigen, !is.na(symp_prob_12_24f))%>%
  filter(timepoint =="52 weeks")%>%
  ggplot(., aes(x=symp_prob_12_24f, y=MFI, fill=treatmentarm))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  facet_wrap(~treatmentarm+antigen, nrow=2)+
  ylab("MFI at 52 weeks")+
  scale_fill_manual(values=treatment_palette)+
  theme_minimal()+
  theme(legend.position = "bottom",
        axis.title.x = element_blank())

symp_prob_purf_104 <- long_luminex%>%
  filter(timepoint=="104 weeks", symp_prop_24_36<=1)%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~glm(symp_prop_24_36~log_mfi+gender_categorical+treatmentarm, data=., family = "binomial" , weights = total_n_para_12_24))) %>%
  # mutate(n_malaria_model=map(data, ~glmmTMB::glmmTMB(total_n_malaria_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  # mutate(n_para_model=map(data, ~glmmTMB::glmmTMB(total_n_para_12~conc+log_qpcr, data=., family=nbinom2))) %>%
  mutate(n_malaria_model_summary=map(n_malaria_model, ~summary(.))) %>%
  #11 when additional covariate is included
  mutate(n_malaria_p=map_dbl(n_malaria_model_summary, ~coef(.)[14]))%>%
  ungroup()%>%
  mutate(n_malaria_padj=p.adjust(n_malaria_p, method="BH"))

symp_prob_sigs_104 <- symp_prob_purf_104%>%
  filter(n_malaria_padj < 0.03 )%>%
  select(antigen,n_malaria_p, n_malaria_padj )%>%
  arrange(n_malaria_padj)



long_luminex %>%
  mutate(symp_prob_24_36f=case_when(symp_prop_24_36<0.25~"< 25%",
                                    symp_prop_24_36>=0.25&symp_prop_24_36<0.5~"25-50%",
                                    symp_prop_24_36>=0.5&symp_prop_24_36<0.75~"50-75%",
                                    symp_prop_24_36>=0.75~"> 75%",
                                    total_n_para_24_36==0~"no parasitemia"))%>%
  mutate(symp_prob_24_36f=factor(symp_prob_24_36f, levels=c("< 25%", "25-50%", "50-75%", "> 75%")))%>%
  filter(antigen %in% symp_prob_sigs_104$antigen, !is.na(symp_prob_24_36f))%>%
  filter(timepoint =="104 weeks")%>%
  ggplot(., aes(x=symp_prob_24_36f, y=MFI, fill=treatmentarm))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  facet_wrap(~treatmentarm+antigen, nrow=2)+
  ylab("MFI at 104 weeks")+
  scale_fill_manual(values=treatment_palette)+
  theme_minimal()+
  theme(legend.position = "bottom",
        axis.title.x = element_blank())


# pca plots ####

id_columns <- c("id", "Plate", "location", "treatmentarm", "timepoint", "gender")

wide_df2 <- long_luminex %>%
  select(-treatmentarm)%>%
  pivot_wider(names_from = antigen, values_from = log_mfi)%>%
  na.omit()
rownames(wide_df2) <- paste(wide_df2$location)

# data_for_louise <- long_luminex %>%
#   pivot_wider(names_from = antigen, values_from = MFI, id_cols = c("Plate", "Well", "id", "timepoint"))
# 
# write.csv(data_for_louise, "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/compiled_data_with_ids.csv", row.names = F)

# each row = measurement; each column = feature
big_pca <-  prcomp(wide_df2[,(length(id_columns)+1):ncol(wide_df2)], center = T, scale. = T)
pca_plot_data <- as.data.frame(cbind(wide_df2, big_pca$x))

loadings_df <- data.frame(big_pca$rotation)
loadings_df$targetName <- rownames(loadings_df)
loadings_df$targetName <- factor(loadings_df$targetName, levels = loadings_df$targetName[order(loadings_df$PC1)])


PC1_plot <- ggplot(loadings_df, aes(x=factor(targetName, levels = targetName[order(loadings_df$PC1)]), y=-PC1, fill=targetName))+
  geom_bar(stat = "identity")+
  scale_fill_viridis_d(name = "mako")+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")



plate_pca_plot <- pca_plot_data %>%
  ggplot(., aes(x=PC1, y=PC2, color=factor(Plate)))+
  geom_point()+
  stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_viridis_d()+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/plate_pca_plot.png", plate_pca_plot, width=5, height=4, bg="white", dpi=444)



dp_pca_plot <- pca_plot_data %>%
  ggplot(., aes(x=PC1, y=PC2, color=factor(treatmentarm)))+
  geom_point()+
  stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=c("darkred", "#00555A"))+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/dp_pca_plot.png", dp_pca_plot, width=5, height=4, bg="white", dpi=444)



dp_pca_plot2 <- pca_plot_data %>%
  ggplot(., aes(x=PC1, y=PC2, color=factor(treatmentarm)))+
  geom_point()+
  stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  facet_wrap(~timepoint)+
  theme_minimal()+
  scale_color_manual(values=c("darkred", "#00555A"))+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/dp_pca_plot2.png", dp_pca_plot2, width=8, height=8, bg="white", dpi=444)



time_pca_plot <- pca_plot_data %>%
  ggplot(., aes(x=PC1, y=PC2, color=factor(timepoint)))+
  geom_point()+
  stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  viridis::scale_color_viridis(option = "turbo", discrete = T, direction = -1)+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/time_pca_plot.png", time_pca_plot, width=5, height=4, bg="white", dpi=444)



gender_pca_plot <- pca_plot_data %>%
  ggplot(., aes(x=PC1, y=PC2, color=factor(gender)))+
  geom_point()+
  stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  facet_wrap(~timepoint)+
  viridis::scale_color_viridis(option = "turbo", discrete = T, direction = -1)+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/time_pca_plot.png", time_pca_plot, width=5, height=4, bg="white", dpi=444)


# individual plots ####

for(i in(unique(long_luminex$antigen))){

file_name <-gsub(" ", "_", i)

plt <- long_luminex %>%
  filter(antigen==i)%>%
  filter(MFI>100)%>%
  ggplot(., aes(x=timepoint, y=MFI, color=timepoint))+
  geom_line(aes(group=id_cat), alpha=0.15)+
  geom_violin()+
  geom_point()+
  scale_y_log10()+
  ggtitle(paste(i))+
  facet_wrap(~treatmentarm, scales="fixed", labeller = label_wrap_gen(width = 10))+
  viridis::scale_color_viridis(option = "turbo", discrete = T, direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        axis.text.x = element_blank(),
        axis.title.x=element_blank())
  
ggsave(paste("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/individual_plots/indie_", file_name, ".png",sep=""), plt, width=6, height=4, bg="white", dpi=444)

  }


wide_df3 <- long_luminex %>%
  select(-treatmentarm, -MFI)%>%
  pivot_wider(names_from = antigen, values_from = log_mfi)

write.csv(wide_df3, "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/wide_micdrop_luminex_df.csv", row.names = F)



long_luminex%>%
  filter(#antigen %in% c("schizont"),
    timepoint=="52 weeks", 
    treatmentarm!="DP 2 years")%>%
  ggplot(aes(x=log_qpcr, y=log_mfi, color=antigen))+
  ggpubr::stat_cor(method="spearman", na.rm = T)+
  # geom_line(aes(group=id), alpha=0.2)+
  geom_smooth(method="lm")+
  geom_point(alpha=0.5)+
  facet_wrap(~targetName, scales = "free")+
  viridis::scale_fill_viridis(discrete = T)+
  xlab("log10 parasites per microlitre")+
  facet_wrap(~antigen, scales = "free")+
  scale_color_viridis_d(option = "D")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_blank())
)

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/pardens_nulisa_cor_plot.png", pardens_nulisa_cor_plot, width=8, height=4, dpi=444, bg="white")



# heatmap for paper ####

(p1 <- long_luminex %>%
  filter(antigen!="tetanus", antigen!="var2csa", timepoint %in% c("8 weeks", "24 weeks", "52 weeks")) %>%
  mutate(log_mfi = log10(MFI + 1)) %>%
  # group_by(antigen, timepoint, treatmentarm) %>%
  # summarise(mean_log_mfi = mean(log_mfi, na.rm = TRUE)) %>%
  ggplot(aes(x = factor(id), y = antigen, fill = log_mfi)) +
  geom_tile() +
  # scale_fill_gradientn(colors=viridis::mako(n=5))+
  scale_fill_gradient2(
    low = "#233875",
    mid = "white",
    high = "#FF5A00",
    midpoint = 2.89
  ) +
  facet_wrap(~treatmentarm+timepoint, scales="free_x") +
  theme_minimal()+
   theme(axis.text.x = element_blank())
)

