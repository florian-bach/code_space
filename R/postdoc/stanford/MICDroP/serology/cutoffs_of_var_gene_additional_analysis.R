censored_long_luminex <- long_luminex%>%
  filter(MFI>500)


# complicated malaria ####
kids_with_comp <- raw_data%>%
  filter(mstatus==2)%>%
  select(id, date, dob, rogerson)%>%
  mutate(flo_age_in_wks=as.numeric(floor((date-dob)/7)))

kids_with_comp_before_130 <- kids_with_comp%>%
  filter(flo_age_in_wks>104 & flo_age_in_wks<130)

kids_with_comp_anytime_before_3 <- kids_with_comp%>%
  filter(flo_age_in_wks & flo_age_in_wks<156)

comp_0_36_per_antigen <- censored_long_luminex %>%
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

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/comp_0_36_per_antigen_500_plus.png", comp_0_36_per_antigen, width=14, height=10, dpi=444, bg="white")


comp_0_36_per_receptor <- censored_long_luminex %>%
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

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/comp_0_36_per_receptor_500_plus.png", comp_0_36_per_receptor, width=8, height=6, dpi=444, bg="white")

# clinical malaria future ####
## year 2 ####
year2_n_malaria_antigen <- censored_long_luminex %>%
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

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/year2_n_malaria_antigen_500_plus.png", year2_n_malaria_antigen, width=14, height=10, dpi=444, bg="white")

year2_symp_prob_antigen <- censored_long_luminex %>%
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

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/year2_symp_prob_antigen_500_plus.png", year2_symp_prob_antigen, width=14, height=10, dpi=444, bg="white")





year2_n_malaria_receptor <- censored_long_luminex %>%
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

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/year2_n_malaria_receptor_500_plus.png", year2_n_malaria_receptor, width=8, height=6, dpi=444, bg="white")

year2_symp_prob_receptor <- censored_long_luminex %>%
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

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/year2_symp_prob_receptor_500_plus.png", year2_symp_prob_receptor, width=8, height=6, dpi=444, bg="white")



## year 3 ####
year3_n_malaria_antigen <- censored_long_luminex %>%
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

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/year3_n_malaria_antigen_500_plus.png", year3_n_malaria_antigen, width=14, height=10, dpi=444, bg="white")

year3_symp_prob_antigen <- censored_long_luminex %>%
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

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/year3_symp_prob_antigen_500_plus.png", year3_symp_prob_antigen, width=14, height=10, dpi=444, bg="white")





year3_n_malaria_receptor <- censored_long_luminex %>%
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

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/year3_n_malaria_receptor_500_plus.png", year3_n_malaria_receptor, width=8, height=6, dpi=444, bg="white")

year3_symp_prob_receptor <- censored_long_luminex %>%
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

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/year3_symp_prob_receptor_500_plus.png", year3_symp_prob_receptor, width=8, height=6, dpi=444, bg="white")







# binder timecourse ####

binder_timecourse <- censored_long_luminex%>%
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

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/binder_timecourse_500_plus.png", binder_timecourse, width=10, height=8, dpi=444, bg="white")

binder_by_timepoint <- censored_long_luminex%>%
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

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/binder_by_timepoint_500_plus.png", binder_by_timepoint, width=10, height=8, dpi=444, bg="white")


# var2csa ####
wide_luminex <- censored_long_luminex%>%
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
ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/var2csa_corr_500_plus.png", var2csa_corr, width=8, height=8, dpi=444, bg="white")


## 8 week var2csa ####

var2csa_luminex <- wide_luminex%>%
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

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/var2csa_corr8_500_plus.png", var2csa_corr8, width=8, height=8, dpi=444, bg="white")


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

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/var2csa_corr_receptor_500_plus.png", var2csa_corr_receptor, width=8, height=8, dpi=444, bg="white")

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

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/var2csa_corr8_receptor_500_plus.png", var2csa_corr8_receptor, width=8, height=8, dpi=444, bg="white")

# parasitemia correlation####
parasite_receptor_cor52 <- censored_long_luminex %>%
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

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/parasite_receptor_cor52_500_plus.png", parasite_receptor_cor52, width=8, height=8, dpi=444, bg="white")


parasite_antigen_cor <- censored_long_luminex %>%
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

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/figures/parasite_antigen_cor_500_plus.png", parasite_antigen_cor, width=12, height=10, dpi=444, bg="white")


censored_long_luminex %>%
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
