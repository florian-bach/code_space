
wide_flagless <- flagless_df %>%
  # filter(antigen %in% c("CSP GENOVA", "GEXP", "PfSEA", "Tet Tox ", "Rh5"))%>%
  dplyr::select(-flag)%>%
  filter(!is.na(conc))%>%
  pivot_wider(names_from = antigen, values_from = conc)

# no_na_flagless <- na.omit(wide_flagless)

trans_wide_flagless <- cbind(wide_flagless[,1:11], apply(log10(wide_flagless[,12:33]), 2, function(x) scale(x, center = TRUE, scale = TRUE)))

long_trans_flagless <- trans_wide_flagless %>%
  pivot_longer(colnames(trans_wide_flagless)[12:33], names_to = "antigen", values_to = "zscore")

visits_with_positive_lamp <- bc1_all%>%
  mutate("age_at_lamp"=date-dob)%>%
  mutate("age_cat_at_lamp"=case_when(age_at_lamp<180 ~ "2",
                                     age_at_lamp>180 & age_at_lamp<365 ~"3",
                                     age_at_lamp > 365 ~ "4"
  ))%>%
  # filter(age_cat_at_lamp %in% c("2", "3"))%>%
  filter(lamp_result==1)%>%
  mutate("sample_id"=paste(id, age_cat_at_lamp, sep="_"))%>%
  group_by(sample_id)%>%
  summarise("n_positives"=n())







malaria_age_12_24 <- bc1_all%>%
  mutate("age_at_visit"=date-dob)%>%
  filter(MSPinfectiondich==1, age_at_visit<730 & age_at_visit>365)%>%
  group_by(id)%>%
  summarise("n_positives"=n())

lamp_age_12_24 <- bc1_all%>%
  mutate("age_at_visit"=date-dob)%>%
  filter(lamp_result==1, age_at_visit<730 & age_at_visit>365)%>%
  group_by(id)%>%
  summarise("n_positives"=n())

flagless_df$n_malaria_12_24 <- malaria_age_12_24$n_positives[match(flagless_df$id, malaria_age_12_24$id)]
flagless_df$n_lamp_12_24 <- lamp_age_12_24$n_positives[match(flagless_df$id, lamp_age_12_24$id)]

flagless_df$n_malaria_12_24 <- ifelse(is.na(flagless_df$n_malaria_12_24), 0, flagless_df$n_malaria_12_24)
flagless_df$n_lamp_12_24 <- ifelse(is.na(flagless_df$n_lamp_12_24), 0, flagless_df$n_lamp_12_24)




malaria_age_12_15 <- bc1_all%>%
  mutate("age_at_visit"=date-dob)%>%
  filter(MSPinfectiondich==1, age_at_visit<455 & age_at_visit>365)%>%
  group_by(id)%>%
  summarise("n_positives"=n())

lamp_age_12_15 <- bc1_all%>%
  mutate("age_at_visit"=date-dob)%>%
  filter(lamp_result==1, age_at_visit<455 & age_at_visit>365)%>%
  group_by(id)%>%
  summarise("n_positives"=n())


flagless_df$n_malaria_12_15 <- malaria_age_12_15$n_positives[match(flagless_df$id, malaria_age_12_15$id)]
flagless_df$n_lamp_12_15 <- lamp_age_12_15$n_positives[match(flagless_df$id, lamp_age_12_15$id)]

flagless_df$n_malaria_12_15 <- ifelse(is.na(flagless_df$n_malaria_12_15), 0, flagless_df$n_malaria_12_15)
flagless_df$n_lamp_12_15 <- ifelse(is.na(flagless_df$n_lamp_12_15), 0, flagless_df$n_lamp_12_15)


malaria_age_12_13 <- bc1_all%>%
  mutate("age_at_visit"=date-dob)%>%
  filter(MSPinfectiondich==1, age_at_visit<395 & age_at_visit>365)%>%
  group_by(id)%>%
  summarise("n_positives"=n())

lamp_age_12_13 <- bc1_all%>%
  mutate("age_at_visit"=date-dob)%>%
  filter(lamp_result==1, age_at_visit<395 & age_at_visit>365)%>%
  group_by(id)%>%
  summarise("n_positives"=n())

flagless_df$n_malaria_12_13 <- malaria_age_12_13$n_positives[match(flagless_df$id, malaria_age_12_13$id)]
flagless_df$n_lamp_12_13 <- lamp_age_12_13$n_positives[match(flagless_df$id, lamp_age_12_13$id, )]

flagless_df$n_malaria_12_13 <- ifelse(is.na(flagless_df$n_malaria_12_13), 0, flagless_df$n_malaria_12_13)
flagless_df$n_lamp_12_13 <- ifelse(is.na(flagless_df$n_lamp_12_13), 0, flagless_df$n_lamp_12_13)







episodes_12_24_plot <- flagless_df %>%
  filter(antigen %in% modelable_antigens)%>%
  arrange(n_malaria_12_24)%>%
  ggplot(aes(x=n_malaria_12_24, y=conc, color=n_malaria_12_24))+
  geom_point()+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)))+
  geom_smooth(formula = "y~x", method="lm")+
  scale_x_continuous(breaks=seq(0,10))+
  xlab("number of malaria episodes in months 12-24")+
  ylab("log concentration")+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/malaria_episodes_vs_antibodies_12_24.png", episodes_12_24_plot, width=8, height=6, bg="white")




lamp_12_24_plot <- flagless_df %>%
  filter(antigen %in% modelable_antigens)%>%
  arrange(n_lamp_12_24)%>%
  ggplot(aes(x=n_lamp_12_24, y=conc, color=n_lamp_12_24))+
  geom_point()+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)))+
  geom_smooth(formula = "y~x", method="lm")+
  scale_x_continuous(breaks=seq(0,10))+
  xlab("number of positive lamp in months 12-24")+
  ylab("log concentration")+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/positive_lamp_vs_antibodies_12_24.png", lamp_12_24_plot, width=8, height=6, bg="white")






episodes_12_15_plot <- flagless_df %>%
  filter(antigen %in% modelable_antigens)%>%
  arrange(n_malaria_12_15)%>%
  ggplot(aes(x=n_malaria_12_15, y=conc, color=n_malaria_12_15))+
  geom_point()+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)))+
  geom_smooth(formula = "y~x", method="lm")+
  scale_x_continuous(breaks=seq(1,10))+
  xlab("number of malaria episodes in months 12-15")+
  ylab("log concentration")+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/malaria_episodes_vs_antibodies_12_15.png", episodes_12_15_plot, width=8, height=6, bg="white")




lamp_12_15_plot <- flagless_df %>%
  filter(antigen %in% modelable_antigens)%>%
  arrange(n_lamp_12_15)%>%
  ggplot(aes(x=n_lamp_12_15, y=conc, color=n_lamp_12_15))+
  geom_point()+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)))+
  geom_smooth(formula = "y~x", method="lm")+
  scale_x_continuous(breaks=seq(1,10))+
  xlab("number of positive lamp in months 12-15")+
  ylab("log concentration")+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/positive_lamp_vs_antibodies_12_15.png", lamp_12_15_plot, width=8, height=6, bg="white")







episodes_12_13_plot <- flagless_df %>%
  filter(antigen %in% modelable_antigens)%>%
  arrange(n_malaria_12_13)%>%
  ggplot(aes(x=n_malaria_12_13, y=conc, color=n_malaria_12_13))+
  geom_point()+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)))+
  geom_smooth(formula = "y~x", method="lm")+
  scale_x_continuous(breaks=seq(1,10))+
  xlab("number of malaria episodes in months 12-13")+
  ylab("log concentration")+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/malaria_episodes_vs_antibodies_12_13.png", episodes_12_13_plot, width=8, height=6, bg="white")




lamp_12_13_plot <- flagless_df %>%
  filter(antigen %in% modelable_antigens)%>%
  arrange(n_lamp_12_13)%>%
  ggplot(aes(x=n_lamp_12_13, y=conc, color=n_lamp_12_13))+
  geom_point()+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)))+
  geom_smooth(formula = "y~x", method="lm")+
  scale_x_continuous(breaks=seq(1,10))+
  xlab("number of positive lamp in months 12-13")+
  ylab("log concentration")+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/positive_lamp_vs_antibodies_12_13.png", lamp_12_13_plot, width=8, height=6, bg="white")







post12_purf <- flagless_df %>%
  filter(antigen %in% modelable_antigens)%>%
  group_by(antigen) %>%
  nest() %>%
  mutate(dich_model=map(data, ~lm(log(conc)~n_malaria_12_13, data=.))) %>%
  mutate(dich_summary=map(dich_model, ~summary(.)))%>%
  mutate(lamp_model=map(data, ~lm(log(conc)~n_lamp_12_13, data=.))) %>%
  mutate(lamp_summary=map(lamp_model, ~summary(.)))%>%
  mutate(lamp_p=map_dbl(lamp_model, ~summary(.)$coefficients[2,"Pr(>|t|)"]))%>%
  mutate(dich_p=map_dbl(dich_model, ~summary(.)$coefficients[2,"Pr(>|t|)"]))

# no significance whatsoever when predicting the number of lamp OR malaria episodes across
# months 12-24 of life by modelable antibodies, with or without including all the zeroes

# no significance in months 12-15
# Hyp2 "significant" in first month (positive coefficient), but only two observations




malaria_age_0_6 <- bc1_all%>%
  mutate("age_at_visit"=date-dob)%>%
  filter(MSPinfectiondich==1, age_at_visit<181 & age_at_visit>0)%>%
  group_by(id)%>%
  summarise("n_positives"=n())

lamp_age_0_6 <- bc1_all%>%
  mutate("age_at_visit"=date-dob)%>%
  filter(lamp_result==1, age_at_visit<181 & age_at_visit>0)%>%
  group_by(id)%>%
  summarise("n_positives"=n())


flagless_df$n_malaria_0_6 <- malaria_age_0_6$n_positives[match(flagless_df$id, malaria_age_0_6$id)]
flagless_df$n_lamp_0_6 <- lamp_age_0_6$n_positives[match(flagless_df$id, lamp_age_0_6$id)]

flagless_df$n_malaria_0_6 <- ifelse(is.na(flagless_df$n_malaria_0_6), 0, flagless_df$n_malaria_0_6)
flagless_df$n_lamp_0_6 <- ifelse(is.na(flagless_df$n_lamp_0_6), 0, flagless_df$n_lamp_0_6)



t2_flagless_df <- filter(flagless_df, timepoint==2)
t3_flagless_df <- filter(flagless_df, timepoint==3)


t23_flagless_df <- filter(flagless_df, timepoint %in% c(2,3))

t23_flagless_df$conc_2 <- t2_flagless_df$conc[match(t23_flagless_df$id, t2_flagless_df$id)]
t23_flagless_df$conc_3 <- t3_flagless_df$conc[match(t23_flagless_df$id, t3_flagless_df$id)]

t23_flagless_df$conc_2 <- ifelse(is.na(t23_flagless_df$conc_2), 0, t23_flagless_df$conc_2)
t23_flagless_df$conc_3 <- ifelse(is.na(t23_flagless_df$conc_3), 0, t23_flagless_df$conc_3)




croos_plot <- t23_flagless_df%>%
  filter(antigen %in% modelable_antigens)%>%
  arrange(n_malaria_0_6)%>%
  ggplot(., aes(x=conc_2, y=conc_3, color=factor(n_malaria_0_6)))+
  geom_point()+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)))+
  theme_minimal()+
  scale_y_log10()+
  scale_x_log10()+
  xlab("[Ab] at 6 months")+
  ylab("[Ab] at 12 months")+
  scale_color_manual(values=c("grey", "#FFA500", "#FC6A03", "darkred"))+
  guides(color=guide_legend(title="# of smear+ in months 0-6", title.position = "right", override.aes = list(size=3)),
         # shape=guide_legend(title="# symptomatic malaria"),
         size=guide_none(),
         alpha=guide_none())+
  theme(legend.key.height = unit(2.5, "cm"),
        legend.title = element_text(size = 12, angle = 90),
        legend.title.align = 0.5,
        legend.direction = "vertical",
        panel.grid.major.x = element_blank()
  )


ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/612_dich_plot.png", croos_plot, width=8, height=6, bg="white")





croos2_plot <- t23_flagless_df%>%
  filter(antigen %in% modelable_antigens)%>%
  arrange(n_lamp_0_6)%>%
  ggplot(., aes(x=conc_2, y=conc_3, color=factor(n_lamp_0_6)))+
  geom_point()+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)))+
  theme_minimal()+
  scale_y_log10()+
  scale_x_log10()+
  xlab("[Ab] at 6 months")+
  ylab("[Ab] at 12 months")+
  scale_color_manual(values=c("grey", "#FFA500", "#FC6A03", "darkred"))+
  guides(color=guide_legend(title="# of LAMP+ in months 0-6", title.position = "right", override.aes = list(size=3)),
         # shape=guide_legend(title="# symptomatic malaria"),
         size=guide_none(),
         alpha=guide_none())+
  theme(legend.key.height = unit(2.5, "cm"),
        legend.title = element_text(size = 12, angle = 90),
        legend.title.align = 0.5,
        legend.direction = "vertical",
        panel.grid.major.x = element_blank()
  )


ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/612_lamp_plot.png", croos2_plot, width=8, height=6, bg="white")
