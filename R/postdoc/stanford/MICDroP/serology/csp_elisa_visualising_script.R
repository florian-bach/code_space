# read in the data ####
library(tidyr)
library(dplyr)
library(ggplot2)
`%notin%` <- Negate(`%in%`)
treatment_palette <- c("Placebo"="darkred", "DP 1 year"="darkblue", "DP 2 years"="orange")

mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")
mic_drop_hbs <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/MICDROP SickleTr final.dta")

elisa_data1 <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/plate_1_and_2_data.csv")
elisa_data2 <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/plate_3_and_4_data.csv")
elisa_data3 <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/plate_5_and_6_data.csv")
elisa_data4 <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/plate_7_and_8_data.csv")%>%
  select(sample_id, plate, plate.row, plate.column, OD450.540)

standards1 <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/plate_1_and_2_standards.csv")
standards2 <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/plate_3_and_4_standards.csv")
standards3 <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/plate_5678_standards.csv")
two_year_samples <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/104_week_timepoints_to_pick.csv")

cleaned_standards <- bind_rows(standards1, standards2, standards3)%>%
  group_by(Plate, Dilution)%>%
  summarise("mean_OD"=mean(Read1, Read2))

master_micdrop_metadata <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/micdrop_metadata_for_immunology.csv")
master_micdrop_metadata <- master_micdrop_metadata%>%
  mutate(id=as.character(id))

elisa_data <- elisa_data1%>%
  bind_rows(elisa_data2, elisa_data3, elisa_data4)

cleaned_elisa_data <- elisa_data%>%
  mutate(id=ifelse(plate<7, substr(sample_id, 1,5), sample_id))%>%
  mutate(date=as.Date(ifelse(plate<7,
                     as.Date(lubridate::parse_date_time(substr(sample_id, 7,15), "%m/%d/%y")),
                     as.Date(two_year_samples$date[match(id, two_year_samples$id)]))))%>%
  group_by(id, plate, date)%>%
  summarise(mean_OD=mean(`OD450.540`))%>%
  left_join(., master_micdrop_metadata, by="id")%>%
  mutate(timepoint=case_when(plate<7~"52 weeks",
                             plate>=7~"104 weeks"))%>%
  ungroup()%>%
  group_by(plate)%>%
  mutate(sample_number=seq(1, n()))

# normalize ####
list_of_data <- split(cleaned_elisa_data, cleaned_elisa_data$plate)
list_of_standards <- split(cleaned_standards, cleaned_standards$Plate)
list_of_std_conc <- c(1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 0)


analyze_elisa_base <- function(std_conc, std_signal, sample_signal, offset = 1e-6) {
  # Add offset to avoid log(0) or div-by-zero
  std_conc_adj <- std_conc + offset
  
  # Define the 4PL model
  model_formula <- std_signal ~ (A - D) / (1 + (std_conc_adj / C)^B) + D
  
  # Reasonable starting values
  start_vals <- list(
    A = max(std_signal),
    D = min(std_signal),
    C = median(std_conc_adj),
    B = 1
  )
  
  # Try to fit the model
  fit <- try(
    nls(model_formula, start = start_vals,
        control = nls.control(maxiter = 200, warnOnly = TRUE)),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    warning("⚠️ 4PL model failed to converge. Using linear fit instead.")
    lin_fit <- lm(std_signal ~ std_conc_adj)
    pred <- (sample_signal - coef(lin_fit)[1]) / coef(lin_fit)[2]
    return(pred)
  }
  
  coefs <- coef(fit)
  A <- coefs["A"]; D <- coefs["D"]; B <- coefs["B"]; C <- coefs["C"]
  
  # Inverse 4PL: absorbance -> concentration
  predict_conc <- function(y) {
    (((A - D) / (y - D) - 1)^(1 / B)) * C - offset
  }
  
  pred <- sapply(sample_signal, predict_conc)
  return(pred)
}

list_of_normalised_results <- lapply(seq(1:8), function(x)analyze_elisa_base(list_of_std_conc[1:6],
                                                                             list_of_standards[[x]]$mean_OD[1:6],
                                                                             list_of_data[[x]]$mean_OD,
                                                                             offset = 1e-6))

list_of_normalised_dfs <- lapply(seq(1:8), function(x){
  data.frame("norm_OD"=list_of_normalised_results[[x]],
             "plate"=x,
             "sample_number"=seq(1, length(list_of_normalised_results[[x]])))
}
)

norm_df <- do.call(rbind, list_of_normalised_dfs)


elisa_metadata <- cleaned_elisa_data%>%
  select(-mean_OD)

# add static metadata ####
norm_df_with_meta <- norm_df%>%
  left_join(., elisa_metadata, by=c("plate", "sample_number"))%>%
  mutate(seropositive=ifelse(norm_OD>0.3154542, 1, 0))


# add visit-level metadata
raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/Oct25/MICDSpecimenBoxOct25_withclinical.dta")
metadata_columns <- c("id", "date", "mstatus", "ageinwks", "pardens", "qPCRparsdens", "fever", "febrile", "gender")

dobs <- raw_data%>%
  distinct(id, dob)%>%
  filter(!is.na(dob))

visit_metadata <- raw_data%>%
  select(all_of(metadata_columns))%>%
  # right_join(., all_samples, by=c("id", "date"))%>%
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

norm_df_with_meta <- norm_df_with_meta%>%
  mutate(old_id=id)%>%
  mutate(id=as.numeric(id))%>%
  left_join(., visit_metadata, by=c("id", "date"))
  
write.csv(norm_df_with_meta, "~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/csp_elisa_norm_df_with_meta.csv", row.names = F)
# figures ####
## overview figures visualise batch data ####
raw_csp_time_plot <- cleaned_elisa_data%>%
  filter(!is.na(treatmentarm))%>%
  ggplot(., aes(x=factor(timepoint, levels = c("52 weeks", "104 weeks")), y=mean_OD, fill=treatmentarm))+
  geom_boxplot(, outliers = F)+
  geom_point(shape=21, position=position_jitterdodge(jitter.width =  0.1, dodge.width = 0.75))+
  ggpubr::stat_compare_means(
    aes(group = treatmentarm),
    label = "p.format",        # shows formatted p-values (e.g., 0.001, 0.05)
    method = "wilcox.test",    # or "t.test", depending on your data
    comparisons = NULL,        # auto by group
    ref.group = NULL, vjust=-0.31
  )+
  # geom_hline(yintercept =0.31)+
  # geom_violin(draw_quantiles = seq(0,1,0.25), aes(fill=treatmentarm))+
  ylab("raw OD (AU)")+
  theme_minimal()+
  scale_fill_manual(values=treatment_palette)+
  theme(legend.position = "right",
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(size=14),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/figures/raw_csp_time_plot.png", raw_csp_time_plot, width=5, height = 4, dpi=444, bg="white")


(norm_csp_time_plot <- norm_df_with_meta%>%
  filter(old_id %notin% c("donor3", "pool"), timepoint=="52 weeks")%>%
  ggplot(., aes(x=treatmentarm, y=norm_OD, fill=treatmentarm))+
  geom_violin()+
  geom_boxplot(color="white", width=0.2, outliers = F, position = position_dodge(width=0.75))+
  geom_point(shape=21, size=1.5, color="white", position=position_jitter(width=0.02), stroke = 0.5)+
  ggpubr::stat_compare_means(
    aes(group = treatmentarm),
    label = "p.format",        # shows formatted p-values (e.g., 0.001, 0.05)
    method = "wilcox.test",    # or "t.test", depending on your data
    vjust=0.51, label.x = 1.35
  )+
  # geom_hline(yintercept =0.31)+
  # geom_violin(draw_quantiles = seq(0,1,0.25), aes(fill=treatmentarm))+
  ylab("normalized intensity (AU)")+
  theme_minimal(base_size = 15)+
  scale_fill_manual(values=treatment_palette)+
  theme(legend.position = "none",
        legend.title = element_blank(),
        #axis.text.y = element_blank(),
        axis.text = element_text(size=14),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank()))

ggsave("~/postdoc/stanford/manuscripts/MICDROP_clinial/norm_csp_time_plot.png", norm_csp_time_plot, width=5, height = 4, dpi=444, bg="white")
ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/ELISA/normalized_csp_52weeks.png", norm_csp_time_plot, width=5, height = 4, dpi=444, bg="white")
ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/ELISA/normalized_csp_52weeks.pdf", norm_csp_time_plot, width=5, height = 4, dpi=444, bg="white")

## infections first year of life ####
para12_plot <- norm_df_with_meta%>%
  filter(old_id %notin% c("donor3", "pool"))%>%
  filter(!is.na(total_n_para_12), timepoint=="52 weeks")%>%
  ggplot(., aes(x=factor(total_n_para_12), y=norm_OD, fill=factor(total_n_para_12)))+
  geom_boxplot(, outliers = F)+
  geom_point(shape=21, stroke=0.5, position=position_jitterdodge(jitter.width =  0.1, dodge.width = 0.75))+
  theme_minimal(base_size = 15)+
  ggpubr::stat_cor(method="spearman", mapping = aes(x=total_n_para_12, y=norm_OD, group=treatmentarm), inherit.aes = F)+
  facet_wrap(~treatmentarm)+
  xlab("\nparasitemic months in the first year of life")+
  ylab("normalized intensity (AU) at 52 weeks\n")+
  scale_fill_manual(values=viridis::magma(n=13))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(size=11),
        axis.title = element_text(size=14.5))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/figures/norm_para12_plot.png", para12_plot, width=6, height = 4, dpi=444, bg="white")
ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/ELISA/normalized_csp_parasite_incidence.png", para12_plot, width=5, height = 4, dpi=444, bg="white")
ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/ELISA/normalized_csp_parasite_incidence.pdf", para12_plot, width=5, height = 4, dpi=444, bg="white")



mala12_plot <- norm_df_with_meta%>%
  filter(old_id %notin% c("donor3", "pool"))%>%
  filter(!is.na(total_n_malaria_12), timepoint=="52 weeks")%>%
  ggplot(., aes(x=factor(total_n_malaria_12), y=norm_OD, fill=factor(total_n_malaria_12)))+
  geom_boxplot(, outliers = F)+
  geom_point(shape=21, position=position_jitterdodge(jitter.width =  0.1, dodge.width = 0.75))+
  theme_minimal()+
  ggpubr::stat_cor(method="spearman", mapping = aes(x=total_n_malaria_12, y=norm_OD, group=treatmentarm), inherit.aes = F)+
  facet_wrap(~treatmentarm, scales="free_x")+
  xlab("\nmalaria episodes in the first year of life")+
  ylab("normalized intensity (AU) at 52 weeks\n")+
  scale_fill_manual(values=viridis::magma(n=8)[-c(1)])+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(size=11),
        axis.title = element_text(size=14.5))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/figures/norm_mala12_plot.png", mala12_plot, width=6, height = 4, dpi=444, bg="white")


## infections in the second year of life ####
(norm52_para12_24plot <- norm_df_with_meta%>%
  filter(id %notin% c("donor3", "pool"))%>%
  filter(!is.na(total_n_para_12_24), timepoint=="52 weeks")%>%
  ggplot(.)+
  # geom_smooth(aes(x=total_n_para_12_24, y=norm_OD), inherit.aes = F)+
  geom_boxplot(aes(x=factor(total_n_para_12_24), y=norm_OD, fill=factor(total_n_para_12_24)), inherit.aes = F, outliers = F)+
  geom_point(aes(x=factor(total_n_para_12_24), y=norm_OD, fill=factor(total_n_para_12_24)), shape=21, position=position_jitterdodge(jitter.width =  0.1, dodge.width = 0.75))+
  ggpubr::stat_cor(method="spearman", mapping = aes(x=total_n_para_12_24, y=norm_OD), inherit.aes = F)+
  geom_smooth(aes(x=total_n_para_12_24, y=norm_OD), inherit.aes = F, method = "lm", color="black")+
  geom_boxplot(aes(x=factor(total_n_para_12_24), y=norm_OD, fill=factor(total_n_para_12_24)), inherit.aes = F, outliers = F)+
  # facet_wrap(~treatmentarm, scales="free_x")+
  xlab("\nmonths with any parasitemia in year 1-2")+
  ylab("normalized intensity (AU) at 1 year\n")+
  scale_fill_manual(values=viridis::magma(n=14))+
  theme_minimal(base_size = 20)+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text.y = element_blank()))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/figures/norm52_para12_24plot.png", norm52_para12_24plot, width=8, height = 6, dpi=444, bg="white")

para12_24plot <- norm_df_with_meta%>%
  filter(id %notin% c("donor3", "pool"))%>%
  filter(!is.na(total_n_para_12_24), timepoint=="104 weeks")%>%
  ggplot(., aes(x=factor(total_n_para_12_24), y=norm_OD, fill=factor(total_n_para_12_24)))+
  geom_boxplot(, outliers = F)+
  geom_point(shape=21, position=position_jitterdodge(jitter.width =  0.1, dodge.width = 0.75))+
  theme_minimal()+
  ggpubr::stat_cor(method="spearman", mapping = aes(x=total_n_para_12_24, y=norm_OD, group=treatmentarm), inherit.aes = F)+
  facet_wrap(~treatmentarm, scales="free_x")+
  xlab("\nmonths with any parasitemia in year 1-2")+
  ylab("normalized intensity (AU) at 2 years\n")+
  scale_fill_manual(values=viridis::magma(n=13))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(size=11),
        axis.title = element_text(size=14.5))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/figures/norm_para12_24plot.png", para12_24plot, width=6, height = 4, dpi=444, bg="white")



mala12_24plot <- norm_df_with_meta%>%
  filter(old_id %notin% c("donor3", "pool"))%>%
  filter(!is.na(total_n_malaria_12_24), timepoint=="104 weeks")%>%
  ggplot(., aes(x=factor(total_n_malaria_12_24), y=norm_OD, fill=factor(total_n_malaria_12_24)))+
  geom_boxplot(, outliers = F)+
  geom_point(shape=21, position=position_jitterdodge(jitter.width =  0.1, dodge.width = 0.75))+
  theme_minimal()+
  ggpubr::stat_cor(method="spearman", mapping = aes(x=total_n_malaria_12_24, y=norm_OD, group=treatmentarm), inherit.aes = F)+
  facet_wrap(~treatmentarm, scales="free_x")+
  xlab("\nmalaria episodes in the second year of life")+
  ylab("normalized intensity (AU) at 104 weeks\n")+
  scale_fill_manual(values=viridis::magma(n=9))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(size=11),
        axis.title = element_text(size=14.5))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/figures/norm_mala12_24plot.png", mala12_24plot, width=6, height = 4, dpi=444, bg="white")


## parasitemia ####

(parasitemia_plot <- norm_df_with_meta%>%
  filter(old_id %notin% c("donor3", "pool"), timepoint=="52 weeks", treatmentarm=="Placebo")%>%
  ggplot(., aes(x=10^log_qpcr, y=norm_OD, color=treatmentarm))+
  geom_smooth(method="lm")+
  geom_point(position = position_jitter(width=0.1))+
  # geom_hline(yintercept =0.31)+
  # geom_violin(draw_quantiles = seq(0,1,0.25), aes(fill=treatmentarm))+
  xlab("parasites / μL")+
  ylab("normalized intensity (AU)")+
  ggpubr::stat_cor(method="spearman", label.y = c(1.5, 1.5))+
  facet_wrap(~treatmentarm)+
  theme_minimal(base_size = 15)+
  scale_color_manual(values=treatment_palette)+
  scale_x_log10(labels=scales::label_log(), breaks = c(10^-3, 10^-1, 10^1, 10^3, 10^5))+
  scale_y_continuous(limits = c(-0.2, NA))+
  guides(color = guide_legend(override.aes = list(label = "")))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        #axis.text.y = element_blank(),
        axis.text = element_text(size=14),
        axis.title.y = element_text(size=16),
        strip.text = element_text(size=16)))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/figures/norm_parasitemia_plot.png", parasitemia_plot, width=6, height = 4, dpi=444, bg="white")
ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/ELISA/normalized_csp_parasitemia_correlation.png", parasitemia_plot, width=5, height = 4, dpi=444, bg="white")
ggsave("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/ELISA/normalized_csp_parasitemia_correlation.pdf", parasitemia_plot, width=5, height = 4, dpi=444, bg="white")


para12_12_24_plot <- master_micdrop_metadata%>%
  distinct(id, treatmentarm, total_n_para_12, total_n_para_12_24)%>%
  filter(!is.na(total_n_para_12_24))%>%
  ggplot(., aes(x=total_n_para_12, y=total_n_para_12_24))+
  geom_point(position = position_jitter(height=0.3, width=0.1))+
  geom_smooth(method = "lm")+
  scale_y_continuous(breaks = seq(1,12))+
  scale_x_continuous(breaks = seq(1,12))+
  xlab("number of parasitemic months in year 1")+
  ylab('number of parasitemic months in year 2')+
  ggpubr::stat_cor()+
  # facet_wrap(~treatmentarm)+
  theme_minimal(base_size = 20)

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/CSP_elisas/figures/para12_12_24_plot.png", para12_12_24_plot, width=12, height = 6, dpi=444, bg="white")
