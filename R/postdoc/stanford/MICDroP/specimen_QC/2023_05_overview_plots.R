library(dplyr)
library(tidyr)
library(ggplot2)


mic_drop <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2023_05/MICDROP all visit database through May 31st 2023.dta")


smear_positive <- mic_drop %>%
  filter(BSdich==1)%>%
  mutate(id=factor(id))%>%
  group_by(id) %>%
  add_count(name="total_positive_smears") %>%
  arrange(compage) %>%
  mutate(n_smear = seq(1, max(total_positive_smears)))


n_infection_cols <- c("white", colorspace::sequential_hcl(n=max(smear_positive$n_smear), palette = "Lajolla"))


# overview plots ####

micro_dens <- ggplot(smear_positive, aes(x=factor(n_smear), y=pardens, fill=factor(n_smear)))+
  geom_point(alpha=0.3)+
  geom_line(alpha=0.3, aes(group=id))+
  geom_boxplot(aes(group=factor(n_smear)))+
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", width = 0.8, color="white")+
  ylab("Microscopic Parasite Density")+
  xlab("Order of Positive Smear")+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_manual(values=n_infection_cols[-1])


pcr_dens <- ggplot(smear_positive, aes(x=factor(n_smear), y=qPCRparsdens, fill=factor(n_smear)))+
  geom_point(alpha=0.3)+
  geom_line(alpha=0.3, aes(group=id))+
  geom_boxplot(aes(group=factor(n_smear)))+
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", width = 0.8, color="white")+
  ylab("qPCR Parasite Density")+
  xlab("Order of Positive Smear")+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_manual(values=n_infection_cols[-1])

dens_n_smear_combo <- cowplot::plot_grid(micro_dens, pcr_dens,  nrow = 1)






micro_age <- ggplot(smear_positive, aes(x=AGE, y=pardens, fill=factor(AGE)))+
  # geom_point(alpha=0.2)+
  # geom_line(alpha=0.2, aes(group=id))+
  geom_boxplot(aes(group=factor(AGE)))+
  # stat_summary(fun.y = median, fatten=0, fun.min = median, fun.max = median,
  #              geom = "crossbar", width = 0.8, color="white")+
  ylab("Microscopic Parasite Density")+
  xlab("Age (Weeks)")+
  geom_smooth()+
  scale_y_log10(limits=c(0.01, 10^6))+
  theme_minimal()+
  theme(legend.position = "none")+
  viridis::scale_fill_viridis(discrete = TRUE)


pcr_age <- ggplot(smear_positive, aes(x=AGE, y=qPCRparsdens+0.001, fill=factor(AGE)))+
  geom_point(alpha=0.3)+
  geom_line(alpha=0.3, aes(group=id))+
  geom_boxplot(aes(group=factor(AGE)))+
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", width = 0.8, color="white")+
  ylab("qPCR Parasite Density")+
  xlab("Age (Weeks)")+
  scale_y_log10(limits=c(0.01, 10^6))+
  theme_minimal()+
  theme(legend.position = "none")+
  viridis::scale_fill_viridis(discrete = TRUE)

dens_age_combo <- cowplot::plot_grid(micro_age, pcr_age,  nrow = 1)


