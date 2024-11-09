histopath_lm <- combo_data %>%
  filter(!is.na(matmal), !is.na(conc), timepoint==1)%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(model=map(data, ~lm(log10(conc)~matmal, data = .)))%>%
  mutate("no_vs_non_placental" = map_dbl(model, ~summary(.)$coefficients[11]))%>%
  mutate("no_vs_placental" = map_dbl(model, ~summary(.)$coefficients[12]))%>%
  ungroup()%>%
  mutate("no_vs_non_placental_adj" = p.adjust(no_vs_non_placental, method="BH"))%>%
  mutate("no_vs_placental_adj" = p.adjust(no_vs_placental, method="BH"))

 
histopath_glm <- combo_data %>%
  filter(!is.na(anyHPfinal), anyHPfinal!=9, !is.na(conc), timepoint==1)%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(model=map(data, ~glm(anyHPfinal~log10(conc), data = ., family="binomial")))%>%
  mutate("any_vs_placental" = map_dbl(model, ~summary(.)$coefficients[8]))%>%
  ungroup()%>%
  mutate("any_vs_placental_adj" = p.adjust(any_vs_placental, method="fdr"))%>%
  filter(any_vs_placental_adj<0.1)





histopath_wilcox <- combo_data %>%
  filter(!is.na(matmal), !is.na(conc), timepoint==1)%>%
  group_by(antigen)%>%
  pivot_wider(names_from = anyHPfinalx, values_from = conc)%>%
  nest()%>%
  mutate(wilcox=map(data, ~wilcox.test(.$`No Pathology`, .$`Placental Malaria`)))%>%
  mutate(wilcox2=map(data, ~wilcox.test(.$`No Pathology`, .$`Placental Malaria`)))%>%
  mutate(raw_p=map_dbl(wilcox, ~.$p.value))%>%
  ungroup()%>%
  mutate(padj=p.adjust(raw_p, method="BH"))


# child malaria incidence ~ maternal malaria outcome ####
inf_012_matmal_frac <- combo_data %>%
  filter(!duplicated(id), !is.na(matmal))%>%
  group_by(inf_0_12, matmal)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = matmal, values_from = n)%>%
  mutate(sum=sum(`no malaria`, `non-placental malaria`, `placental malaria`, na.rm = TRUE))%>%
  pivot_longer(c(`no malaria`, `non-placental malaria`, `placental malaria`), names_to = "matmal", values_to = "n")%>%
  mutate(n=if_else(is.na(n), 0, n), "frac"=n/sum)%>%
  ggplot(., aes(x=inf_0_12, y=frac, fill=matmal))+
  geom_bar(stat="identity")+
  geom_text(aes(x=inf_0_12, label=paste0("n=", sum, sep=""), y=0.75))+
  xlab("Number of parasitaemic episodes in the first year of life")+
  ylab("Fraction of maternal malaria status")+
  scale_fill_manual(values=disease_pal)+
  theme_minimal()+
  theme(legend.title = element_blank())

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/inf_012_matmal_frac.png", inf_012_matmal_frac, height=3, width=6, bg="white", dpi=444)



any_inf_012_matmal_frac <- combo_data %>%
  filter(!duplicated(id), !is.na(matmal))%>%
  mutate("any_inf_0_12"=if_else(inf_0_12==0, "none", "some"))%>%
  group_by(any_inf_0_12, matmal)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = matmal, values_from = n)%>%
  mutate(sum=sum(`no malaria`, `non-placental malaria`, `placental malaria`, na.rm = TRUE))%>%
  pivot_longer(c(`no malaria`, `non-placental malaria`, `placental malaria`), names_to = "matmal", values_to = "n")%>%
  mutate(n=if_else(is.na(n), 0, n), "frac"=n/sum)%>%
  ggplot(., aes(x=any_inf_0_12, y=frac, fill=matmal))+
  geom_bar(stat="identity")+
  geom_text(aes(label=paste("n = ", sum, sep="")), y=0.75)+
  xlab("Number of parasitaemic episodes in the first year of life")+
  ylab("Fraction of maternal malaria status")+
  scale_fill_manual(values=rev(time_palette))+
  theme_minimal()+
  theme(legend.title = element_blank())

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/any_inf_012_matmal_frac.png", any_inf_012_matmal_frac, height=3, width=6, bg="white", dpi=444)



inf_06_matmal_frac <- combo_data %>%
  filter(!duplicated(id), !is.na(matmal))%>%
  group_by(inf_0_6, matmal)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = matmal, values_from = n)%>%
  mutate(sum=sum(`no malaria`, `non-placental malaria`, `placental malaria`, na.rm = TRUE))%>%
  pivot_longer(c(`no malaria`, `non-placental malaria`, `placental malaria`), names_to = "matmal", values_to = "n")%>%
  mutate(n=if_else(is.na(n), 0, n), "frac"=n/sum)%>%
  ggplot(., aes(x=inf_0_6, y=frac, fill=matmal))+
  geom_bar(stat="identity")+
  xlab("Number of parasitaemic episodes in the first year of life")+
  ylab("Fraction of maternal malaria status")+
  scale_fill_manual(values=rev(time_palette))+
  theme_minimal()


inf_012_matmal_hist <- combo_data %>%
  filter(!duplicated(id), !is.na(matmal))%>%
  group_by(inf_0_12, matmal)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = matmal, values_from = n)%>%
  mutate(sum=sum(`no malaria`, `non-placental malaria`, `placental malaria`, na.rm = TRUE))%>%
  pivot_longer(c(`no malaria`, `non-placental malaria`, `placental malaria`), names_to = "matmal", values_to = "n")%>%
  mutate(n=if_else(is.na(n), 0, n), "frac"=n/sum)%>%
  ggplot(., aes(x=inf_0_12, y=n, fill=matmal))+
  geom_bar(stat="identity", position = position_dodge())+
  xlab("Number of parasitaemic episodes in the first year of life")+
  ylab("Fraction of maternal malaria status")+
  scale_fill_manual(values=rev(time_palette))+
  theme_minimal()



matmal_df <- combo_data %>%
  filter(!duplicated(id), !is.na(matmal))%>%
  mutate(any_inf_0_12=if_else(inf_0_12==0, 0, 1))

# significant
matmal_glm <- glm(data=matmal_df, inf_0_12~matmal, family=poisson)
#significant
matmal_glm2 <- glm(data=matmal_df, inf_0_12~matmal+id, family=poisson)
# p=0.115
matmal_glm3 <- glm(data=matmal_df, any_inf_0_12~matmal, family=binomial)
matmal_glm4 <- glm(data=matmal_df, inf_0_12~matmal)

#significant
matmal_glmm1 <- MASS::glmmPQL(data=matmal_df, inf_0_12~matmal , random = ~1 | id,family=poisson)
#not significant
matmal_glmm2 <- lme4::glmer(inf_0_12~matmal+(1|id), data=matmal_df, family=poisson)


# any maternal malaria
anymalariapreg_inf012 <- combo_data %>%
  filter(!duplicated(id), !is.na(matmal))%>%
  group_by(inf_0_12, anymalariapreg)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = anymalariapreg, values_from = n)%>%
  mutate(sum=sum(`0`, `1`, na.rm = TRUE))%>%
  pivot_longer(c(`0`, `1`,), names_to = "matmal", values_to = "n")%>%
  mutate(n=if_else(is.na(n), 0, n), "frac"=n/sum)%>%
  ggplot(., aes(x=inf_0_12, y=frac, fill=matmal))+
  geom_bar(stat="identity")+
  geom_text(aes(x=inf_0_12, label=sum), y=0.75)+
  xlab("Number of parasitaemic episodes in the first year of life")+
  ylab("Fraction of any maternal malaria")+
  scale_fill_manual(values=rev(time_palette))+
  theme_minimal()




# maternal treatment ~ maternal malaria ####
momrx_matmal_frac <- combo_data %>%
  filter(!duplicated(id), !is.na(matmal), !is.na(MomFinalRx))%>%
  group_by(MomFinalRx, matmal)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = matmal, values_from = n)%>%
  mutate(sum=sum(`no malaria`, `non-placental malaria`, `placental malaria`, na.rm = TRUE))%>%
  pivot_longer(c(`no malaria`, `non-placental malaria`, `placental malaria`), names_to = "matmal", values_to = "n")%>%
  mutate(n=if_else(is.na(n), 0, n), "frac"=n/sum)%>%
  ggplot(., aes(x=MomFinalRx, y=frac, fill=matmal))+
  geom_bar(stat="identity")+
  geom_text(aes(label=paste("n = ", sum, sep="")), y=0.75)+
  ggtitle("DP lower")
  xlab("Number of parasitaemic episodes in the first year of life")+
  ylab("Maternal Treatment Arm")+
  scale_fill_manual(values=rev(time_palette))+
  theme_minimal()

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/momrx_matmal_frac.png", momrx_matmal_frac, height=3, width=6, bg="white", dpi=444)





momrx_matmal_hist <- combo_data %>%
  filter(!duplicated(id), !is.na(matmal), !is.na(MomFinalRx))%>%
  group_by(MomFinalRx, matmal)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = matmal, values_from = n)%>%
  mutate(sum=sum(`no malaria`, `non-placental malaria`, `placental malaria`, na.rm = TRUE))%>%
  pivot_longer(c(`no malaria`, `non-placental malaria`, `placental malaria`), names_to = "matmal", values_to = "n")%>%
  mutate(n=if_else(is.na(n), 0, n), "frac"=n/sum)%>%
  ggplot(., aes(x=MomFinalRx, y=n, fill=matmal))+
  geom_bar(stat="identity", position = position_dodge())+
  xlab("Number of parasitaemic episodes in the first year of life")+
  ylab("Maternal Treatment Arm")+
  scale_fill_manual(values=rev(time_palette))+
  theme_minimal()


# child malaria incidence ~ maternal treatment
inf_012_matmal_frac <- combo_data %>%
  filter(!duplicated(id), !is.na(matmal))%>%
  group_by(inf_0_12, MomFinalRx)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = MomFinalRx, values_from = n)%>%
  mutate(sum=sum(c(`3 Dose SP`, `3 Dose DP`, `Monthly DP`), na.rm = TRUE))%>%
  pivot_longer(c(`3 Dose SP`, `3 Dose DP`, `Monthly DP`), names_to = "MomFinalRx", values_to = "n")%>%
  mutate(n=if_else(is.na(n), 0, n), "frac"=n/sum)%>%
  ggplot(., aes(x=inf_0_12, y=frac, fill=MomFinalRx))+
  geom_bar(stat="identity")+
  xlab("Number of parasitaemic episodes in the first year of life")+
  ylab("Fraction of maternal treatment")+
  scale_fill_manual(values=rev(time_palette))+
  theme_minimal()
