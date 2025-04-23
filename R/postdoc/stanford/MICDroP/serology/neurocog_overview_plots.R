# relationships between maternal and infant characteristics
library(patchwork)
library(tidyr)
library(dplyr)
library(ggplot2)

nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")

kid_epi_data <- nulisa_data%>%
  distinct(id, timepoint_num, gender_categorical, log_qpcr, total_n_para_6, total_n_malaria_6, total_n_para_12, total_n_malaria_12, total_n_para_24, total_n_malaria_24)%>%
  mutate(age=timepoint_num,
         child_id=id)


neuro_cog <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/NCT_infants_share032025.csv")
maternal_neuro_cog_edit <- neuro_cog%>%
  mutate(child_id=subjid,
         mother_id=subjid-10000)

combo_data <- full_join(kid_epi_data, maternal_neuro_cog_edit, by="child_id")%>%
  pivot_longer(cols = ends_with("composite"), names_to = "composite_kind", values_to = "composite_score")%>%
  mutate(any_malaria_12=if_else(total_n_malaria_12>0, 1, 0),
         any_para_12=if_else(total_n_para_12>0, 1, 0),
         composite_score2=composite_score/5)%>%
  group_by(id)%>%
  filter(!duplicated(composite_kind))
  


#maternal characteristics ####

mat_age_plot1<- combo_data%>%
  ggplot(., aes(x=AGE, y=composite_score))+
  geom_point(alpha=0.3)+
  geom_smooth(method="lm", linewidth=0.31)+
  ggpubr::stat_cor(method="spearman", size=2)+
  # ggpubr::stat_compare_means(label.y = 26, size=2)+
  facet_wrap(~composite_kind)+
  scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  ylab("composite score")+
  xlab("maternal age")+
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))

mat_age_plot2<- combo_data%>%
  ggplot(., aes(x=AGE, y=composite_score, color = gender_categorical))+
  geom_point(alpha=0.3)+
  geom_smooth(method="lm", linewidth=0.31)+
  ggpubr::stat_cor(method="spearman", size=2)+
  # ggpubr::stat_compare_means(label.y = 26, size=2)+
  facet_wrap(~composite_kind)+
  scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  ylab("composite score")+
  xlab("maternal age")+
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))
mat_age_plot <- mat_age_plot1 + mat_age_plot2

ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/mat_age_plot.png", mat_age_plot, width = 8, height=4, bg="white", dpi=444)


mat_age_grav_plot <- combo_data%>%
  distinct(id, gravidcat, AGE)%>%
  ggplot(., aes(x=AGE, color=factor(gravidcat)))+
  geom_density()+
  # ggpubr::stat_cor(method="spearman", size=2)+
  # ggpubr::stat_compare_means(label.y = 130, size=2, label = "p.signif", comparisons = list( c("0", "1")))+
  # scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  scale_color_manual(values=c("#F89F5B", "#DA70D6",  "#9D0759"))+
  ggtitle("gravidity")+
  xlab("maternal age")+
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))


mat_bmi_plot1<- combo_data%>%
  ggplot(., aes(x=BMIenrol, y=composite_score))+
  geom_point(alpha=0.3)+
  geom_smooth(method="lm", linewidth=0.31)+
  ggpubr::stat_cor(method="spearman", size=2)+
  # ggpubr::stat_compare_means(label.y = 26, size=2)+
  facet_wrap(~composite_kind)+
  scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  ylab("composite score")+
  xlab("maternal BMI at enrolment")+
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))

mat_bmi_plot2<- combo_data%>%
  ggplot(., aes(x=BMIenrol, y=composite_score, color=gender_categorical))+
  geom_point(alpha=0.3)+
  geom_smooth(method="lm", linewidth=0.31)+
  ggpubr::stat_cor(method="spearman", size=2)+
  # ggpubr::stat_compare_means(label.y = 26, size=2)+
  facet_wrap(~composite_kind)+
  scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  ylab("composite score")+
  xlab("maternal BMI at enrolment")+
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))

mat_bmi_plot <- mat_bmi_plot1 +mat_bmi_plot2
ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/mat_bmi_plot.png", mat_bmi_plot, width = 8, height=4, bg="white", dpi=444)



mat_mal_plot<- combo_data%>%
  ggplot(., aes(x=factor(evermalaria), y=composite_score, fill = factor(evermalaria)))+
  geom_violin(draw_quantiles=seq(0,1,0.25), color="white")+
  # ggpubr::stat_cor(method="spearman", size=2)+
  ggpubr::stat_compare_means(label.y = 130, size=2, label = "p.signif", comparisons = list( c("0", "1")))+
  facet_wrap(~composite_kind)+
  scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  ylab("composite score")+
  xlab("maternal ever malaria")+
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/mat_mal_plot.png", mat_mal_plot, width = 4, height=4, bg="white", dpi=444)


mat_hp_plot<- combo_data%>%
  filter(!is.na(anyHP))%>%
  ggplot(., aes(x=factor(anyHP), y=composite_score, fill = factor(anyHP)))+
  geom_violin(draw_quantiles=seq(0,1,0.25), color="white")+
  # ggpubr::stat_cor(method="spearman", size=2)+
  ggpubr::stat_compare_means(label.y = 130, size=2, label = "p.signif", comparisons = list( c("0", "1")))+
  facet_wrap(~composite_kind)+
  scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  scale_fill_manual(values=c("grey", "darkred"))+
  ylab("composite score")+
  xlab("any hp")+
  theme(legend.position = "none",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/mat_hp_plot.png", mat_hp_plot, width = 4, height=4, bg="white", dpi=444)



gravidity_plot<- combo_data%>%
  filter(!is.na(gravidcat))%>%
  ggplot(., aes(x=factor(gravidcat), y=composite_score, fill = factor(gravidcat)))+
  geom_violin(draw_quantiles=seq(0,1,0.25), color="white")+
  # ggpubr::stat_cor(method="spearman", size=2)+
  ggpubr::stat_compare_means(label.y = 130, size=2, label = "p.signif",
                             comparisons = list( c("1", "2"), c("2", "3"), c("1", "3")))+
  facet_wrap(~composite_kind)+
  scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  ylab("composite score")+
  xlab("gravidity")+
  scale_fill_manual(values=c("#F89F5B", "#DA70D6",  "#9D0759"))+
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/mat_gravidity_plot.png", gravidity_plot, width = 5, height=4, bg="white", dpi=444)

n_fever_plot <- combo_data%>%
  ggplot(., aes(x=factor(totalfever), y=composite_score, fill = factor(totalfever)))+
  geom_boxplot()+
  facet_wrap(~composite_kind, scales="free_x")+
  scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  ylab("composite score")+
  xlab("fevers")+
  viridis::scale_fill_viridis(discrete = T)+
  theme(legend.position = "none",
        legend.title = element_blank())+
  guides(fill=guide_legend(nrow = 1))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/mat_n_fever_plot.png", n_fever_plot, width = 5, height=4, bg="white", dpi=444)


mat_fever_plot<- combo_data%>%
  ggplot(., aes(x=factor(everfever), y=composite_score, fill = factor(everfever)))+
  geom_violin(draw_quantiles=seq(0,1,0.25), color="white")+
  # ggpubr::stat_cor(method="spearman", size=2)+
  ggpubr::stat_compare_means(label.y = 130, size=2, label = "p.signif", comparisons = list( c("0", "1")))+
  facet_wrap(~composite_kind)+
  scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  ylab("composite score")+
  xlab("maternal ever fever")+
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/mat_fever_plot.png", mat_fever_plot, width = 4, height=4, bg="white", dpi=444)

mat_Olevel_plot<- combo_data%>%
  ggplot(., aes(x=factor(Olevel), y=composite_score, fill = factor(Olevel)))+
  geom_violin(draw_quantiles=seq(0,1,0.25), color="white")+
  # ggpubr::stat_cor(method="spearman", size=2)+
  ggpubr::stat_compare_means(label.y = 130, size=2, label = "p.signif", comparisons = list( c("0", "1")))+
  facet_wrap(~composite_kind)+
  scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  ylab("composite score")+
  xlab("maternal O level education")+
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/mat_Olevel_plot.png", mat_Olevel_plot, width = 4, height=4, bg="white", dpi=444)



# child characteristics ####

sex_plot<- combo_data%>%
  ggplot(., aes(x=gender_categorical, y=composite_score, fill = gender_categorical))+
  geom_violin(draw_quantiles=seq(0,1,0.25), color="white")+
  ggpubr::stat_compare_means(label.y = 130, size=2)+
  facet_wrap(~composite_kind)+
  theme_minimal()+
  ylab("composite score")+
  xlab("")+
  scale_y_continuous(breaks = seq(0,150,5))+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/kid_sex_plot.png", sex_plot, width = 4, height=4, bg="white", dpi=444)


gestage_plot<- combo_data%>%
  filter(!is.na(gravidcat))%>%
  ggplot(., aes(x=GAcomputed, y=composite_score, color=gender_categorical))+
  geom_point(alpha=0.3)+
  geom_smooth(method="lm", linewidth=0.31)+
  ggpubr::stat_cor(method="spearman", size=2)+
 facet_wrap(~composite_kind)+
  scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  ylab("composite score")+
  xlab("gestational age")+
  scale_fill_manual(values=c("#F89F5B", "#9D0759", "#DA70D6"))+
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/kid_gestage_plot_plot.png", gestage_plot, width = 5, height=4, bg="white", dpi=444)



sga_plot<- combo_data%>%
  filter(!is.na(anyHP))%>%
  ggplot(., aes(x=factor(SGA), y=composite_score, fill = factor(SGA)))+
  geom_violin(draw_quantiles=seq(0,1,0.25), color="white")+
  # ggpubr::stat_cor(method="spearman", size=2)+
  ggpubr::stat_compare_means(label.y = 130, size=2, label = "p.signif", comparisons = list( c("0", "1")))+
  facet_wrap(~composite_kind)+
  scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  scale_fill_manual(values=c("grey", "darkred"))+
  ylab("composite score")+
  xlab("small for gestational age")+
  theme(legend.position = "none",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/kid_sga_plot.png", sga_plot, width = 4, height=4, bg="white", dpi=444)


preterm_plot<- combo_data%>%
  filter(!is.na(anyHP))%>%
  ggplot(., aes(x=factor(preterm), y=composite_score, fill = factor(preterm)))+
  geom_violin(draw_quantiles=seq(0,1,0.25), color="white")+
  # ggpubr::stat_cor(method="spearman", size=2)+
  ggpubr::stat_compare_means(label.y = 130, size=2, label = "p.signif", comparisons = list( c("0", "1")))+
  facet_wrap(~composite_kind)+
  scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  scale_fill_manual(values=c("grey", "darkred"))+
  ylab("composite score")+
  xlab("small for gestational age")+
  theme(legend.position = "none",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/kid_preterm_plot.png", preterm_plot, width = 4, height=4, bg="white", dpi=444)



kid_mal_plot<- combo_data%>%
  filter(!is.na(anyHP))%>%
  ggplot(., aes(x=factor(any_malaria_12), y=composite_score, fill = factor(any_malaria_12)))+
  geom_violin(draw_quantiles=seq(0,1,0.25), color="white")+
  # ggpubr::stat_cor(method="spearman", size=2)+
  ggpubr::stat_compare_means(label.y = 130, size=2, label = "p.signif", comparisons = list( c("0", "1")))+
  facet_wrap(~composite_kind)+
  scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  scale_fill_manual(values=c("grey", "darkred"))+
  ylab("composite score")+
  xlab("any child malaria")+
  theme(legend.position = "none",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/kid_mal_plot.png", kid_mal_plot, width = 4, height=4, bg="white", dpi=444)



kid_mal_plot2 <- combo_data%>%
  distinct(id, total_n_para_24, composite_score, composite_kind)%>%
  ggplot(., aes(x=total_n_para_24, y=composite_score))+
  geom_point()+
  geom_smooth(method="lm")+
  # geom_violin(draw_quantiles=seq(0,1,0.25), color="white")+
  # ggpubr::stat_cor(method="spearman", size=2)+
  # ggpubr::stat_compare_means(label.y = 130, size=2, label = "p.signif", comparisons = list( c("0", "2")))+
  facet_wrap(~composite_kind)+
  scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  # scale_fill_manual(values=c("grey", "darkred"))+
  ylab("composite score")+
  xlab("any child malaria")+
  theme(legend.position = "none",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/kid_mal_plot.png", kid_mal_plot, width = 4, height=4, bg="white", dpi=444)


kid_para_plot<- combo_data%>%
  filter(!is.na(anyHP))%>%
  ggplot(., aes(x=factor(any_para_12), y=composite_score, fill = factor(any_para_12)))+
  geom_violin(draw_quantiles=seq(0,1,0.25), color="white")+
  # ggpubr::stat_cor(method="spearman", size=2)+
  ggpubr::stat_compare_means(label.y = 130, size=2, label = "p.signif", comparisons = list( c("0", "1")))+
  facet_wrap(~composite_kind)+
  scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  scale_fill_manual(values=c("grey", "darkred"))+
  ylab("composite score")+
  xlab("any child parasitemia")+
  theme(legend.position = "none",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/kid_para_plot.png", kid_para_plot, width = 4, height=4, bg="white", dpi=444)


# child vs maternal characteristics ####



Olevel_kid_para<- combo_data%>%
  distinct(id, total_n_para_12, Olevel)%>%
  ggplot(., aes(x=total_n_para_12, color=factor(Olevel)))+
  geom_density()+
  # ggpubr::stat_cor(method="spearman", size=2)+
  # ggpubr::stat_compare_means(label.y = 130, size=2, label = "p.signif", comparisons = list( c("0", "1")))+
  # scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  ggtitle("O level education")+
  xlab("number of parasitemic months in first year of life")+
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/Olevel_kid_para_plot.png", Olevel_kid_para, width = 4, height=4, bg="white", dpi=444)


Olevel_kid_malaria<- combo_data%>%
  distinct(id, total_n_malaria_12, Olevel)%>%
  ggplot(., aes(x=total_n_malaria_12, color=factor(Olevel)))+
  geom_density()+
  # ggpubr::stat_cor(method="spearman", size=2)+
  # ggpubr::stat_compare_means(label.y = 130, size=2, label = "p.signif", comparisons = list( c("0", "1")))+
  # scale_y_continuous(breaks = seq(0,150,5))+
  theme_minimal()+
  ggtitle("O level education")+
  xlab("number of malaria episodes in first year of life")+
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/Olevel_kid_malaria.png", Olevel_kid_malaria, width = 4, height=4, bg="white", dpi=444)


check <- combo_data %>%
  distinct(id, total_n_para_12, Olevel, AGE, gender_categorical)

#nothing
olevel_model <- glm(Olevel~total_n_para_12+gender_categorical, data=check, family = "binomial")
#nothing
olevel_model2 <- MASS::glm.nb(total_n_para_12~Olevel+gender_categorical, data=check)
#nothing
mat_age_model <- MASS::glm.nb(total_n_para_12~AGE+gender_categorical, data=check)


#
mat_age_kid_malaria_plot<- combo_data%>%
  distinct(id, any_malaria_12, AGE)%>%
  ggplot(., aes(x=AGE, color=factor(any_malaria_12)))+
  geom_density()+
  theme_minimal()+
  ggtitle("malaria in the first year of life")+
  xlab("maternal age")+
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/mat_age_kid_malaria_plot.png", mat_age_kid_malaria_plot, width = 4, height=4, bg="white", dpi=444)

mat_grav_kid_malaria_plot <- combo_data%>%
  distinct(id, gravidcat, any_malaria_12, AGE)%>%
  ggplot(., aes(x=factor(any_malaria_12), fill=factor(gravidcat)))+
  geom_bar(position = "fill")+
  theme_minimal()+
  scale_fill_manual(values=c("#F89F5B", "#DA70D6",  "#9D0759"))+
  ggtitle("gravidity")+
  scale_y_continuous(labels = scales::label_percent())+
  xlab("any malaria in first year of life")+
  theme(legend.position = "bottom",axis.title.y = element_blank())+
  guides(color=guide_legend(override.aes = list(alpha=0)),
         fill=guide_legend(title="gravidity"))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/figures/mat_grav_kid_malaria_plot.png", mat_grav_kid_malaria_plot, width = 4, height=4, bg="white", dpi=444)

