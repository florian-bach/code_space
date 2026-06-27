ab_incidence_purf <- ab_clin %>%
  group_by(antigen, incidence_type, timepoint)%>%
  filter(!duplicated(id), incidence_type!="symp_6_12")%>%
  nest() %>%
  # negative binomial model
  # mutate(nb_cell_incidence_model=map(data, ~lme4::glmer.nb(incidence_value ~ conc + (1|id), data=.)))%>%
  # mutate(nb_cell_incidence_model_summary=map(nb_cell_incidence_model, ~summary(.)))%>%
  # mutate(nb_cell_incidence_model_summary_p=map_dbl(nb_cell_incidence_model_summary, ~coef(.)[11]))%>%
  # 
  # MASS poisson model
  # mutate(cell_incidence_model=map(data, ~MASS::glmmPQL(data=., incidence_value ~ log10(conc), random = ~1 | id, family=poisson)))%>%
  # mutate(cell_incidence_model=map(data, ~MASS::glm.nb(data=., incidence_value ~ log10(conc))))%>%
  mutate(cell_incidence_model=map(data, ~glm(data=., incidence_value ~ log10(conc), family="poisson")))%>%
  # mutate(cell_incidence_model=map(data, ~lm(data=., log10(conc)~incidence_value)))%>%
  # lme4 poisson model
  # mutate(cell_incidence_model=map(data, ~lme4::glmer(data=., incidence_value ~ log10(conc) + (1 | id), family="poisson")))%>%
  
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  #MASS poisson p
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  #lme4 poisson p
  # mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  group_by(timepoint)%>%
  # mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p, method="BH"))
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p, method="BH"))



ab_poisson_incidence_sigs <- ab_incidence_purf %>%
  dplyr::select(antigen, incidence_type, timepoint, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.05)



ab_incidence_purf2 <- ab_clin %>%
  group_by(antigen, incidence_type, timepoint)%>%
  filter(!duplicated(id), incidence_type!="symp_6_12")%>%
  nest() %>%
  # negative binomial model
  # mutate(nb_cell_incidence_model=map(data, ~lme4::glmer.nb(incidence_value ~ conc + (1|id), data=.)))%>%
  # mutate(nb_cell_incidence_model_summary=map(nb_cell_incidence_model, ~summary(.)))%>%
  # mutate(nb_cell_incidence_model_summary_p=map_dbl(nb_cell_incidence_model_summary, ~coef(.)[11]))%>%
  # 
  # MASS poisson model
  # mutate(cell_incidence_model=map(data, ~MASS::glmmPQL(data=., incidence_value ~ log10(conc), random = ~1 | id, family=poisson)))%>%
  # mutate(cell_incidence_model=map(data, ~MASS::glm.nb(data=., incidence_value ~ log10(conc))))%>%
  mutate(cell_incidence_model=map(data, ~glm(data=., log10(conc)~incidence_value, family="gaussian")))%>%
  # mutate(cell_incidence_model=map(data, ~lm(data=., log10(conc)~incidence_value)))%>%
  # lme4 poisson model
  # mutate(cell_incidence_model=map(data, ~lme4::glmer(data=., incidence_value ~ log10(conc) + (1 | id), family="poisson")))%>%
  
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  #MASS poisson p
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  #lme4 poisson p
  # mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  group_by(timepoint)%>%
  # mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p, method="BH"))
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p, method="BH"))



ab_poisson_incidence_sigs2 <- ab_incidence_purf2 %>%
  dplyr::select(antigen, incidence_type, timepoint, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.1)













df_pred <- wide_long_msd%>%
  distinct(id, antigen, timepoint, titer) %>%
  filter(timepoint == "52 weeks", antigen %in% c(vaccines, vaccines_iga[c(4,5)])) %>%
  rename(
    vaccine_antigen = antigen,
    vaccine_titer = titer
  )

df_target <- wide_long_msd%>%
  filter(timepoint == "24 weeks" & fc_flavour=="log2fc_8_24" |
           timepoint == "52 weeks" & fc_flavour=="log2fc_24_52",
         antigen %in% viruses) %>%
  distinct(id, timepoint, antigen, conversion) %>%
  pivot_wider(values_from = conversion, names_from = timepoint)%>%
  mutate(sero_conversion=case_when(`24 weeks`=="nothing"&`52 weeks`=="nothing"~"nothing",
                                   `24 weeks`=="nothing"&`52 weeks`=="converts 24 to 52"~"converts 24 to 52",
                                   `24 weeks`=="converts 8 to 24"&`52 weeks`=="nothing"~"converts 8 to 24", .default = ))%>%
  distinct(id, antigen, sero_conversion)%>%
  rename(
    virus_antigen = antigen,
    virus_seroconversion = sero_conversion
  )%>%
  filter(!is.na(virus_seroconversion))



wide_long_msd_cross <- df_pred %>%
  inner_join(df_target, by = "id") %>%
  mutate(log2_titer = log2(vaccine_titer))

purf3 <- wide_long_msd_cross%>%
  filter(virus_antigen %notin%c("PIV.1", "PIV.2", "Flu.A.H3"),
         vaccine_antigen!="Pneumo.1.4.14",
         !(virus_antigen=="Mumps"&vaccine_antigen=="Diphtheria"))%>%
  pivot_wider(names_from = virus_seroconversion, values_from = vaccine_titer,
              id_cols = c("id", "vaccine_antigen", "virus_antigen")) %>%
  group_by(vaccine_antigen, virus_antigen) %>%
  filter(sum(!is.na(`converts 8 to 24`)) >= 12,   # set your threshold here
         sum(!is.na(`converts 24 to 52`)) >= 12) %>% 
  nest()%>%
  mutate(wilcox1p=map_dbl(data, ~wilcox.test(x = .$nothing, .$`converts 8 to 24`)$p.value))%>%
  mutate(wilcox2p=map_dbl(data, ~wilcox.test(x = .$nothing, .$`converts 24 to 52`)$p.value))%>%
  mutate(wilcox3p=map_dbl(data, ~wilcox.test(x = .$`converts 8 to 24`, .$`converts 24 to 52`)$p.value))%>%
  ungroup()%>%
  mutate(wilcox1padj = p.adjust(wilcox1p, method="fdr"))%>%
  mutate(wilcox2padj = p.adjust(wilcox2p, method="fdr"))%>%
  mutate(wilcox3padj = p.adjust(wilcox3p, method="fdr"))


virus_conversion_vaccine_sigs <- purf3%>%
  filter(wilcox1padj < 0.1 | wilcox2padj < 0.1 | wilcox3padj < 0.1)%>%
  select(vaccine_antigen, virus_antigen, wilcox1padj, wilcox2padj, wilcox3padj)


# 52 weeks
# virus_conversion_vaccine_sigs
# # A tibble: 4 × 5
# vaccine_antigen virus_antigen wilcox1padj wilcox2padj wilcox3padj
# <chr>           <chr>               <dbl>       <dbl>       <dbl>
#   1 Polio           RV.C              0.00876      0.407       0.158 
# 2 Rotavirus       RV.C              0.00513      0.201       0.0284
# 3 Tetanus         RV.C              0.0906       0.439       0.153 
# 4 Rotavirus_IgA   RV.C              0.593        0.0974      0.0549

library(pwr)

# Wilcoxon power ≈ t-test power / 0.955 (ARE correction)
# Assume a large effect (Cohen's d = 0.8) as the most optimistic scenario
pwr.t2n.test(n1 = 14, n2 = 75,   # 8 converters vs ~81 non-converters
             d = 0.8, sig.level = 0.05)

# Try varying n1 to find where power drops below 0.8