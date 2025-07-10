# preamble ####
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(emmeans)

var_luminex <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/slim_luminex_results.csv", check.names = F)
plate_map <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/final_plate_map.csv")
mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")

raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/May25/MICDSpecimenBoxMay25_withclinical.dta")


plate_map <- plate_map%>%
  mutate(well=paste(plate.row, plate.column, sep=""),
         location=paste(dane.plate, well, sep=""))
         
slim_plate_map <- plate_map%>%
  select(location, id, date, ageinwks)%>%
  mutate("date"=as.Date(lubridate::parse_date_time(.$date, orders="%m/%d/%y")))


colnames(var_luminex)[12:13] <- c("CIDRg3",	"CIDRa6")
colnames(var_luminex)[31:32] <- c("schizont",	"tetanus")
maternal_treatment_arms <- haven::read_dta("~/Library/CloudStorage/Box-Box/DP+SP study/Databases and preliminary findings/Final database used for analyses/DPSP treatment allocation_FINAL.dta")

maternal_gravid_cat <-  haven::read_dta("~/Library/CloudStorage/Box-Box/DP+SP study/Databases and preliminary findings/Final database used for analyses/DPSP individual level analysis database_FINAL.dta")%>%
  distinct(id, gravidcat)

maternal_malaria_cat <-  haven::read_dta("~/Library/CloudStorage/Box-Box/DP+SP study/Databases and preliminary findings/Final database used for analyses/DPSP individual level analysis database_FINAL.dta")%>%
  distinct(id, gravidcat)

raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/May25/MICDSpecimenBoxMay25_withclinical.dta")

# write.csv(all_samples, "~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/all_samples.csv", row.names = F)

# add epi data ####
metadata_columns <- c("id", "dob", "date", "ageinwks", "gender", "mstatus", "qPCRparsdens", "visittype", "fever", "febrile", "rogerson", "anyHP", "GAcomputed", "gi", "SGA", "qPCRdich", "mqPCRparsdens")

all_samples <- distinct(plate_map, id, date)%>%
  mutate("date"=as.Date(lubridate::parse_date_time(.$date, orders="%m/%d/%y")))

#merge based on id and date
metadata <- raw_data%>%
  select(all_of(metadata_columns))%>%
  right_join(., all_samples, by=c("id", "date"))%>%
  mutate(timepoint_num=as.numeric(ageinwks), id=as.numeric(id))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.01))%>%
  mutate(mom_rx=maternal_treatment_arms$treatmentarm[match(id-10000, maternal_treatment_arms$id)],
         mom_rx=case_match(mom_rx,
                           1~"SP",
                           2~"DP",
                           3~"DPSP"),
         gravidcat=maternal_gravid_cat$gravidcat[match(id-10000, maternal_gravid_cat$id)],)

## infection data ####
infs_and_meta <- raw_data%>%
  # filter(ageinwks<54)%>%
  group_by(id) %>%
  mutate("total_n_para_12"=sum(pardens>1&ageinwks<53, na.rm=TRUE),
         "total_n_malaria_12"=sum(mstatus!=0&ageinwks<53, na.rm=TRUE),
         "total_n_para_24"=sum(pardens>1&ageinwks<105, na.rm=TRUE),
         "total_n_malaria_24"=sum(mstatus!=0&ageinwks<105, na.rm=TRUE),
         "total_n_para_12_24"=sum(pardens>1&ageinwks<105&ageinwks>52, na.rm=TRUE),
         "total_n_malaria_12_24"=sum(mstatus!=0&ageinwks<105&ageinwks>52, na.rm=TRUE),
         "total_n_para_6"=sum(qPCRparsdens>1&ageinwks<27, na.rm=TRUE),
         "total_n_malaria_6"=sum(mstatus!=0&ageinwks<27, na.rm=TRUE),
         "any_malar_6"=if_else(total_n_malaria_6==0, FALSE, TRUE),
         "any_malar_12"=if_else(total_n_malaria_12==0, FALSE, TRUE))%>%
  select(id, date, total_n_para_6, total_n_malaria_6, total_n_para_12, total_n_malaria_12, total_n_para_24, total_n_malaria_24, total_n_para_12_24, total_n_malaria_12_24, -gender)%>%
  right_join(., metadata, by=c("id", "date"))%>%
  group_by(id)%>%
  mutate(total_n_para=max(total_n_para_12, na.rm = T),
         total_n_malaria=max(total_n_malaria_12, na.rm = T))%>%
  mutate(timepoint_num=case_when(timepoint_num==9~8,
                                 timepoint_num==25~24,
                                 timepoint_num==53~52,
                                 .default=timepoint_num))

long_luminex <- var_luminex%>%
  mutate(location=paste(Plate, Well, sep=""))%>%
  pivot_longer(cols=c(colnames(var_luminex)[3:32]), names_to = "antigen", values_to = "MFI")%>%
  left_join(., slim_plate_map, by="location")%>%
  mutate(log_mfi=log10(MFI+1))%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(id, mic_drop_key$id)])%>%
  mutate(timepoint=case_when(ageinwks<10~"8 weeks",
                             ageinwks==24~"24 weeks",
                             ageinwks%in%c(52, 53)~"52 weeks",
                             ageinwks%in%c(104, 105)~"104 weeks"),
         timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks", "104 weeks")),
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))%>%
  filter(!is.na(treatmentarm), !is.na(log_mfi), is.finite(log_mfi))%>%
  left_join(., infs_and_meta, by=c("id", "date"))%>%
  mutate(antigen=gsub("_", " ", antigen),
         id_cat=as.character(id))


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
  mutate(time_model=map(data, ~lme4::lmer(log_mfi~timepoint*treatmentarm+(1|id_cat), data=.))) %>%
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
  filter(padj<0.05, contrast=="8 weeks")%>%
  select(antigen, contrast, p, padj)
#nothing
sigs_24 <- treatment_purf%>%
  filter(padj<0.05, contrast=="24 weeks")%>%
  select(antigen, contrast, p, padj)

sigs_52 <- treatment_purf%>%
  filter(padj<0.05, contrast=="52 weeks")%>%
  select(antigen, contrast, p, padj)%>%
  arrange(padj)

sigs_104 <- treatment_purf%>%
  filter(padj<0.05, contrast=="104 weeks")%>%
  select(antigen, contrast, p, padj)

sigs_52_104 <- treatment_purf%>%
  filter(padj<0.05, contrast%in%c("Placebo boost 52 - 104 weeks", "DP boost 52 - 104 weeks"))%>%
  select(antigen, contrast, p, padj)



sig_52_plot1 <- long_luminex %>%
  # mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  filter(antigen %in% sigs_52$antigen[1:8])%>%
  mutate(antigen=factor(antigen, levels=sigs_52$antigen[1:8]))%>%
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
  scale_fill_manual(values=c("darkred", "#00555A"))+
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



n_para_sigs_52_1 <- long_luminex %>%
  filter(antigen %in% inf_sigs_52$antigen[1:9])%>%
  mutate(antigen=factor(antigen, levels=inf_sigs_52$antigen[1:9]))%>%
  filter(timepoint =="52 weeks", treatmentarm=="Placebo")%>%
  mutate(total_n_para_12f=if_else(total_n_para_12>6, "7+", as.character(total_n_para_12)))%>%
  ggplot(., aes(x=factor(total_n_para_12f), y=MFI, fill=factor(total_n_para_12f)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
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
  facet_wrap(~antigen, scales="free", nrow=3)+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_blank())
ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/n_para_sigs_52_2.png", n_para_sigs_52_2, height=8, width=8, dpi=400, bg="white")  

n_para_sigs_52_3 <- long_luminex %>%
  filter(antigen %in% inf_sigs_52$antigen[19:20])%>%
  mutate(antigen=factor(antigen, levels=inf_sigs_52$antigen[19:20]))%>%
  filter(timepoint =="52 weeks", treatmentarm=="Placebo")%>%
  mutate(total_n_para_12f=if_else(total_n_para_12>6, "7+", as.character(total_n_para_12)))%>%
  ggplot(., aes(x=factor(total_n_para_12f), y=MFI, fill=factor(total_n_para_12f)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  facet_wrap(~antigen, scales="free", nrow=3)+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/n_para_sigs_52_3.png", n_para_sigs_52_3, height=8, width=8, dpi=400, bg="white")  

n_malaria_sigs_52 <- long_luminex %>%
  filter(antigen %in% inf_sigs_52$antigen[inf_sigs_52$n_malaria_padj<0.05])%>%
  mutate(antigen=factor(antigen, levels=inf_sigs_52$antigen[inf_sigs_52$n_malaria_padj<0.05]))%>%
  filter(timepoint =="52 weeks", treatmentarm=="Placebo")%>%
  ggplot(., aes(x=factor(total_n_malaria_12), y=MFI, fill=factor(total_n_malaria_12)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  facet_wrap(~antigen, scales="free", nrow=3)+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/n_malaria_sigs_52.png", n_malaria_sigs_52, height=3, width=3, dpi=400, bg="white")  

### 104 weeks ####

inf_purf_104 <- long_luminex%>%
  filter(timepoint=="104 weeks", MFI>100)%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~MASS::glm.nb(total_n_malaria_24~log_mfi, data=.))) %>%
  mutate(n_para_model=map(data, ~MASS::glm.nb(total_n_para_24~log_mfi, data=.))) %>%
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

inf_sigs_104 <- inf_purf_104%>%
  filter(n_malaria_p<0.05,
         n_para_p<0.05)


n_malaria_sigs_104 <- long_luminex %>%
  filter(antigen %in% inf_sigs_52$antigen[1:9])%>%
  mutate(antigen=factor(antigen, levels=inf_sigs_52$antigen[1:9]))%>%
  filter(timepoint =="104 weeks")%>%
  mutate(total_n_para_24f=if_else(total_n_para_24>6, "7+", as.character(total_n_para_24)))%>%
  ggplot(., aes(x=factor(total_n_para_24f), y=MFI, fill=factor(total_n_para_24f)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
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
  facet_wrap(~antigen, scales="free", nrow=3)+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/n_malaria_sigs_104_dp.png", n_malaria_sigs_104_dp, height=8, width=8, dpi=400, bg="white")  


### 52 - 104 weeks ####

inf_purf_52_104 <- long_luminex%>%
  filter(timepoint=="104 weeks", MFI>100)%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(n_malaria_model=map(data, ~MASS::glm.nb(total_n_malaria_12_24~log_mfi, data=.))) %>%
  mutate(n_para_model=map(data, ~MASS::glm.nb(total_n_para_12_24~log_mfi, data=.))) %>%
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

inf_sigs_52_104 <- inf_purf_52_104%>%
  filter(n_malaria_p<0.05,
         n_para_p<0.05)


n_malaria_sigs_52_104 <- long_luminex %>%
  filter(antigen %in% inf_sigs_52$antigen[1:9])%>%
  mutate(antigen=factor(antigen, levels=inf_sigs_52$antigen[1:9]))%>%
  filter(timepoint =="104 weeks")%>%
  mutate(total_n_para_12_24f=if_else(total_n_para_12_24>3, "4+", as.character(total_n_para_12_24)))%>%
  ggplot(., aes(x=factor(total_n_para_12_24f), y=MFI, fill=factor(total_n_para_12_24f)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  facet_wrap(~antigen, scales="free", nrow=3)+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/n_malaria_sigs_52_104.png", n_malaria_sigs_52_104, height=8, width=8, dpi=400, bg="white")  



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

# pca plots ####

id_columns <- c("id", "Plate", "location", "treatmentarm", "timepoint", "gender")

wide_df2 <- long_luminex %>%
  pivot_wider(names_from = antigen, values_from = log_mfi, id_cols = all_of(id_columns))%>%
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
