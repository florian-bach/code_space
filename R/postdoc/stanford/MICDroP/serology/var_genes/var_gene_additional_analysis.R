# complicated malaria ####
kids_with_comp <- raw_data%>%
  filter(mstatus==2)%>%
  select(id, date, dob, rogerson)%>%
  mutate(flo_age_in_wks=as.numeric(floor((date-dob)/7)))

kids_with_comp_before_130 <- kids_with_comp%>%
  filter(flo_age_in_wks>104 & flo_age_in_wks<130)

kids_with_comp_anytime_before_3 <- kids_with_comp%>%
  filter(flo_age_in_wks & flo_age_in_wks<156)

comp_0_36_per_antigen <- long_luminex %>%
  mutate(comp130 = if_else(id %in% c(kids_with_comp_anytime_before_3$id), "has comp", "has not"))%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=factor(timepoint_num), y=MFI))+
  # geom_point()+
  # geom_line(aes(group=id))+
  ggtitle("complicated malaria before 36 months months (n=10)")+
  geom_boxplot(aes(fill=factor(comp130)), outliers = F)+
  scale_y_log10()+
  facet_wrap(~antigen, scales="free", labeller = label_wrap_gen(width = 10))+
  #viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  scale_fill_manual(values=c("darkred", "darkgrey"))+
  ggpubr::stat_compare_means(size=5, vjust=1,aes(group=comp130), hide.ns = T, label="p.signif")+
  xlab("age in weeks")+
  theme_minimal()+
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.title.x=element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/comp_0_36_per_antigen.png", comp_0_36_per_antigen, width=14, height=10, dpi=444, bg="white")


comp_0_36_per_receptor <- long_luminex %>%
  mutate(comp130 = if_else(id %in% c(kids_with_comp_anytime_before_3$id), "has comp", "has not"))%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=factor(timepoint_num), y=MFI))+
  # geom_point()+
  # geom_line(aes(group=id))+
  ggtitle("complicated malaria before 36 months (n=10)")+
  geom_boxplot(aes(fill=factor(comp130)), outliers = F)+
  scale_y_log10()+
  facet_wrap(~receptor, scales="free", labeller = label_wrap_gen(width = 10))+
  #viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  scale_fill_manual(values=c("darkred", "darkgrey"))+
  ggpubr::stat_compare_means(size=5, vjust=1,aes(group=comp130), hide.ns = T, label="p.signif")+
  theme_minimal()+
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.title.x=element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/comp_0_36_per_receptor.png", comp_0_36_per_receptor, width=8, height=6, dpi=444, bg="white")

# clinical malaria future ####
## year 2 ####
year2_n_malaria_antigen <- long_luminex %>%
  filter(timepoint_num==52, !is.na(total_n_malaria_12_24))%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=total_n_malaria_12_24, y=MFI))+
  geom_point(position = position_jitter(width=0.1))+
  geom_smooth(method="lm")+
  scale_y_log10()+
  facet_wrap(~antigen, scales="free", labeller = label_wrap_gen(width = 10))+
  viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  ggpubr::stat_cor(color="red")+
  ylab("MFI at 52 weeks")+
  xlab("number of malaria episodes in months 12-24")+
  theme_minimal()+
  theme(legend.position="bottom",
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/year2_n_malaria_antigen.png", year2_n_malaria_antigen, width=14, height=10, dpi=444, bg="white")

year2_symp_prob_antigen <- long_luminex %>%
  filter(timepoint_num==52, !is.na(total_n_malaria_12_24), !is.na(total_n_para_12_24))%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=total_n_malaria_12_24/total_n_para_12_24, y=MFI))+
  geom_point(position = position_jitter(width=0.1))+
  geom_smooth(method="lm")+
  ylab("MFI at 52 weeks")+
  scale_y_log10()+
  facet_wrap(~antigen, scales="free", labeller = label_wrap_gen(width = 10))+
  viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  ggpubr::stat_cor(color="red")+
  xlab("number of malaria episodes\ndivided by parasitemic months 12-24")+
  theme_minimal()+
  theme(legend.position="bottom",
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/year2_symp_prob_antigen.png", year2_symp_prob_antigen, width=14, height=10, dpi=444, bg="white")





year2_n_malaria_receptor <- long_luminex %>%
  filter(timepoint_num==52, !is.na(total_n_malaria_12_24))%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=total_n_malaria_12_24, y=MFI))+
  geom_point(position = position_jitter(width=0.2))+
  geom_smooth(method="lm")+
  scale_y_log10()+
  facet_wrap(~receptor, scales="free", labeller = label_wrap_gen(width = 10))+
  viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  ggpubr::stat_cor(color="red")+
  ylab("MFI at 52 weeks")+
  xlab("number of malaria episodes in months 12-24")+
  theme_minimal()+
  theme(legend.position="bottom",
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/year2_n_malaria_receptor.png", year2_n_malaria_receptor, width=8, height=6, dpi=444, bg="white")

year2_symp_prob_receptor <- long_luminex %>%
  filter(timepoint_num==52, !is.na(total_n_malaria_12_24), !is.na(total_n_para_12_24))%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=total_n_malaria_12_24/total_n_para_12_24, y=MFI))+
  geom_point(position = position_jitter(width=0.1))+
  geom_smooth(method="lm")+
  ylab("MFI at 52 weeks")+
  scale_y_log10()+
  facet_wrap(~receptor, scales="free", labeller = label_wrap_gen(width = 10))+
  viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  ggpubr::stat_cor(color="red")+
  xlab("number of malaria episodes\ndivided by parasitemic months 12-24")+
  theme_minimal()+
  theme(legend.position="bottom",
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/year2_symp_prob_receptor.png", year2_symp_prob_receptor, width=8, height=6, dpi=444, bg="white")



## year 3 ####
year3_n_malaria_antigen <- long_luminex %>%
  filter(timepoint_num==104, !is.na(total_n_malaria_24_36))%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=total_n_malaria_12_24, y=MFI))+
  geom_point(position = position_jitter(width=0.1))+
  geom_smooth(method="lm")+
  facet_wrap(~antigen, scales="free", labeller = label_wrap_gen(width = 10))+
  viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = 1)+
  # scale_y_log10(limits = c(NA, 70000))+
  ggpubr::stat_cor(color="red")+
  ylab("MFI at 104 weeks")+
  xlab("number of malaria episodes in months 24-36")+
  theme_minimal()+
  theme(legend.position="bottom",
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/year3_n_malaria_antigen.png", year3_n_malaria_antigen, width=14, height=10, dpi=444, bg="white")

year3_symp_prob_antigen <- long_luminex %>%
  filter(timepoint_num==104, !is.na(total_n_malaria_24_36), !is.na(total_n_para_24_36))%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI), frac_24_36=total_n_malaria_24_36/total_n_para_24_36)%>%
  ggplot(., aes(x=frac_24_36, y=MFI))+
  geom_point(position = position_jitter(width=0.1))+
  geom_smooth(method="lm")+
    ylab("MFI at 104 weeks")+
# scale_y_log10()+
  ggpubr::stat_cor(method="spearman", color="red")+
  facet_wrap(~antigen, scales="free", labeller = label_wrap_gen(width = 10))+
  viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  xlab("number of malaria episodes\ndivided by parasitemic months 24-36")+
  ylab("MFI at 104 weeks")+
  theme_minimal()+
  theme(legend.position="bottom",
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/year3_symp_prob_antigen.png", year3_symp_prob_antigen, width=14, height=10, dpi=444, bg="white")





year3_n_malaria_receptor <- long_luminex %>%
  filter(timepoint_num==104, !is.na(total_n_malaria_24_36))%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=total_n_malaria_24_36, y=MFI))+
  geom_point(position = position_jitter(width=0.2))+
  geom_smooth(method="lm")+
  scale_y_log10()+
  facet_wrap(~receptor, scales="free", labeller = label_wrap_gen(width = 10))+
  viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  ggpubr::stat_cor(color="red")+
  xlab("number of malaria episodes in months 24-36")+
  ylab("MFI at 104 weeks")+
  theme_minimal()+
  theme(legend.position="bottom",
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/year3_n_malaria_receptor.png", year3_n_malaria_receptor, width=8, height=6, dpi=444, bg="white")

year3_symp_prob_receptor <- long_luminex %>%
  filter(timepoint_num==104, !is.na(total_n_malaria_24_36), !is.na(total_n_para_24_36))%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=total_n_malaria_24_36/total_n_para_24_36, y=MFI))+
  geom_point(position = position_jitter(width=0.1))+
  geom_smooth(method="lm")+
  scale_y_log10()+
  facet_wrap(~receptor, scales="free", labeller = label_wrap_gen(width = 10))+
  viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  ggpubr::stat_cor(color="red", vjust=11)+
  ylab("MFI at 104 weeks")+
  xlab("number of malaria episodes\ndivided by parasitemic months 24-36")+
  theme_minimal()+
  theme(legend.position="bottom",
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/year3_symp_prob_receptor.png", year3_symp_prob_receptor, width=8, height=6, dpi=444, bg="white")







# binder timecourse ####

binder_timecourse <- long_luminex%>%
  filter(antigen!="tetanus")%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=factor(timepoint), y=MFI))+
  geom_point(aes(color=antigen))+
  geom_line(aes(group=id, color=antigen))+
  geom_boxplot(outliers = F)+
  # geom_boxplot(outliers = F, color="grey")+
  scale_y_log10()+
  facet_wrap(~receptor, scales="free", labeller = label_wrap_gen(width = 10), nrow=2)+
  viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        #legend.direction = "horizontal",
        # axis.text.x = element_blank(),
        axis.title.x=element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/binder_timecourse.png", binder_timecourse, width=10, height=8, dpi=444, bg="white")

binder_by_timepoint <- long_luminex%>%
  filter(antigen!="tetanus")%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=factor(receptor), y=MFI))+
  geom_point(aes(color=antigen))+
  geom_boxplot(aes(fill=receptor), outliers = F)+
  # geom_boxplot(outliers = F, color="grey")+
  scale_y_log10()+
  facet_wrap(~timepoint, scales="free", labeller = label_wrap_gen(width = 10), nrow=2)+
  viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        #legend.direction = "horizontal",
        # axis.text.x = element_blank(),
        axis.title.x=element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/binder_by_timepoint.png", binder_by_timepoint, width=10, height=8, dpi=444, bg="white")


# var2csa ####
wide_luminex <- long_luminex%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  pivot_wider(names_from = antigen, values_from = MFI, id_cols = c(id, timepoint))

## per timepoint correlation ####
var2csa_corr <- wide_luminex%>%
  pivot_longer(cols=colnames(wide_luminex)[c(3:10, 12:ncol(wide_luminex))], names_to = "antigen", values_to = "MFI")%>%
  filter(antigen!="tetanus")%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=var2csa, y=MFI))+
  geom_point(aes(color=antigen))+
  geom_smooth(method="lm")+
  ggpubr::stat_cor()+
  xlab("var2csa MFI")+
  ylab("other MFI")+
  scale_y_log10()+
  scale_x_log10()+
  facet_wrap(~timepoint, scales="free", labeller = label_wrap_gen(width = 10), nrow=2)+
  viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  theme_minimal()+
  theme(legend.position="none")
ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/var2csa_corr.png", var2csa_corr, width=8, height=8, dpi=444, bg="white")


## 8 week var2csa ####

var2csa_luminex <- wide_luminex%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  filter(timepoint=="8 weeks")%>%
  select(id, var2csa)

var2csa_corr8 <- wide_luminex%>%
  select(-var2csa)%>%
  pivot_longer(cols=colnames(wide_luminex)[c(3:10, 12:ncol(wide_luminex))], names_to = "antigen", values_to = "MFI")%>%
  filter(antigen!="tetanus")%>%
  left_join(., var2csa_luminex, by="id")%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  ggplot(., aes(x=var2csa, y=MFI))+
  geom_point(aes(color=antigen))+
  geom_smooth(method="lm")+
  ggpubr::stat_cor()+
  xlab("var2csa MFI @ 8 weeks")+
  ylab("other MFI")+
  scale_y_log10()+
  scale_x_log10()+
  facet_wrap(~timepoint, scales="free", labeller = label_wrap_gen(width = 10), nrow=2)+
  viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  theme_minimal()+
  theme(legend.position="none")

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/var2csa_corr8.png", var2csa_corr8, width=8, height=8, dpi=444, bg="white")


## per receptor ####
var2csa_corr_receptor <- wide_luminex%>%
  pivot_longer(cols=colnames(wide_luminex)[c(3:10, 12:ncol(wide_luminex))], names_to = "antigen", values_to = "MFI")%>%
  filter(antigen!="tetanus")%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  left_join(., slim_var_domains, by="antigen")%>%
  ggplot(., aes(x=var2csa, y=MFI))+
  geom_point(aes(color=antigen))+
  geom_smooth(method="lm")+
  ggpubr::stat_cor()+
  xlab("var2csa MFI")+
  ylab("other MFI")+
  scale_y_log10()+
  scale_x_log10()+
  facet_wrap(~timepoint+receptor, scales="free", labeller = label_wrap_gen(width = 10), nrow=4)+
  viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  theme_minimal()+
  theme(legend.position="none")

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/var2csa_corr_receptor.png", var2csa_corr_receptor, width=8, height=8, dpi=444, bg="white")

var2csa_corr8_receptor <- wide_luminex%>%
  select(-var2csa)%>%
  pivot_longer(cols=colnames(wide_luminex)[c(3:10, 12:ncol(wide_luminex))], names_to = "antigen", values_to = "MFI")%>%
  filter(antigen!="tetanus")%>%
  left_join(., var2csa_luminex, by="id")%>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  left_join(., slim_var_domains, by="antigen")%>%
  ggplot(., aes(x=var2csa, y=MFI))+
  geom_point(aes(color=antigen))+
  geom_smooth(method="lm")+
  ggpubr::stat_cor()+
  xlab("var2csa MFI @ 8 weeks")+
  ylab("other MFI")+
  scale_y_log10()+
  scale_x_log10()+
  facet_wrap(~timepoint+receptor, scales="free", labeller = label_wrap_gen(width = 10), nrow=4)+
  viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  theme_minimal()+
  theme(legend.position="none")

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/var2csa_corr8_receptor.png", var2csa_corr8_receptor, width=8, height=8, dpi=444, bg="white")

# parasitemia correlation####
parasite_receptor_cor52 <- long_luminex %>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  filter(timepoint=="52 weeks")%>%
  ggplot(., aes(x=log_qpcr, y=MFI))+
  geom_point()+
  scale_y_log10()+
  facet_wrap(~receptor+timepoint, scales="free", labeller = label_wrap_gen(width = 10))+
  #viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  scale_fill_manual(values=c("darkred", "darkgrey"))+
  ggpubr::stat_cor(color="red")+
  ggtitle("")+
  xlab("age in weeks")+
  theme_minimal()+
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.title.x=element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/parasite_receptor_cor52.png", parasite_receptor_cor52, width=8, height=8, dpi=444, bg="white")


parasite_antigen_cor <- long_luminex %>%
  mutate(MFI=ifelse(MFI<1, 1, MFI))%>%
  filter(timepoint!="8 weeks")%>%
  ggplot(., aes(x=qPCRparsdens+0.001, y=MFI))+
  geom_point()+
  geom_smooth(method="lm")+
  scale_y_log10()+
  scale_x_log10()+
  facet_wrap(~antigen, scales="free", labeller = label_wrap_gen(width = 10))+
  #viridis::scale_fill_viridis(option = "turbo", discrete = T, direction = -1)+
  scale_fill_manual(values=c("darkred", "darkgrey"))+
  ggpubr::stat_cor(color="red", vjust=5.7)+
  xlab("qPCR parasitemia")+
  theme_minimal()+
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/parasite_antigen_cor.png", parasite_antigen_cor, width=12, height=10, dpi=444, bg="white")

# seropositivity ####

library(mclust)
library(data.table)

seroprev_est_df <- data.table()

for(v in unique(long_luminex$antigen)){
  fmm_model=NULL
  log_serotype_vector <- long_luminex$log_mfi[long_luminex$antigen==v&!is.na(long_luminex$log_mfi)&long_luminex$timepoint_num==52]
  print(v)
  ### fit initial mixutre model with 2 components
  k <- 2 #number of components 
  fmm_model <- mclust::Mclust(log_serotype_vector, G = k, modelNames = "V") 
  
  if(is.null(fmm_model)) print(paste(v, "failed")) 
  if(is.null(fmm_model)) next
  # Extract parameters
  means <- fmm_model$parameters$mean
  sds <- sqrt(fmm_model$parameters$variance$sigmasq)
  props <- fmm_model$parameters$pro
  
  # Create density curves for each component
  x_vals <- seq(min(log_serotype_vector), max(log_serotype_vector), length.out = 1000)
  
  dens_df <- data.frame(
    x = rep(x_vals, 2),
    density = c(
      dnorm(x_vals, mean = means[1], sd = sds[1]) * props[1],
      dnorm(x_vals, mean = means[2], sd = sds[2]) * props[2]
    ),
    component = factor(rep(1:2, each = length(x_vals)))
  )
  
  seropos_threshold <- dens_df%>%
    pivot_wider(names_from = component, names_prefix = "component_", values_from = density)%>%
    filter(component_1/component_2>0.99 &component_1/component_2<1.01)
  
  # Plot histogram + density curves
  gg_density <- ggplot() +
    geom_histogram(aes(x = log_serotype_vector, y = after_stat(density)), 
                   bins = 50, fill = "gray80", color = "white") +
    geom_line(data = dens_df, aes(x = x, y = density, color = component), size = 1.2) +
    geom_vline(xintercept = seropos_threshold$x)+
    scale_x_continuous(limits = c(0, NA))+
    labs(title = "FMM: Mixture of 2 Components",
         subtitle=paste0("Serotype: ", v),
         x = "log(concentration)", y = "Density") +
    scale_color_manual(values = c("blue", "red")) +
    theme_minimal()+
    theme(legend.position = "none")
  
  ggsave(paste0("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/fmm_2components_",v, ".png", sep=""),gg_density, width = 8, height = 4, dpi = 444, bg="white")
  
  # Step 1: Out of the 2 components, we identify which distribution has a higher mean, and assume that is the seropositive component
  component_means <- fmm_model$parameters$mean
  seropositive_component <- which.max(component_means)
  
  # Step 4: Get the posterior probabilities for each sample
  posteriors <- fmm_model$z  #For each sample, this is the probability of the component being in the first component or second component
  seropositive_probs <- posteriors[, seropositive_component] #this pulls out just the probability of each sample being seropositive
  
  # Step 5: Estimate seroprevalence
  estimated_seroprevalence <- mean(seropositive_probs) #the mean of all the probabilities will be equal to the overall population seroprevalence
  
  # Step 6: Report as percentage
  to_add <- cbind(v, estimated_seroprevalence, seropos_threshold$x[1])
  
  seroprev_est_df <- rbind(seroprev_est_df, to_add)
  
}

write.csv(seroprev_est_df, "~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/mixture_model_df.csv", row.names = F)

View(seroprev_est_df)



long_luminex %>%
  filter(timepoint =="8 weeks", !is.na(total_n_para_6), !is.na(total_n_malaria_6))%>%
  mutate(symp_porb_6=total_n_malaria_6/total_n_para_6)%>%
  filter(is.finite(symp_porb_6))%>%
  ggplot(., aes(x=symp_porb_6, y=MFI, ))+
  geom_point()+
  scale_y_log10()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(color="red", vjust=8)+
  ylab("MFI at 8 weeks")+
  xlab("number of malaria episodes\n divided by parasitemic months in the first six months of life")+
  facet_wrap(~receptor, scales="free", nrow=3)+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(legend.position = "none")
