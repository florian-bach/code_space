library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

`%notin%` <- Negate(`%in%`)

clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")

clean_data%>%
  distinct(id, timepoint, total_n_malaria_12)%>%
  group_by(total_n_malaria_12)%>%
  summarise(n())

neuro_cog <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/NCT_infants_share032025.csv")
neuro_cog_edit <- neuro_cog%>%
  mutate(id=subjid, date=as.Date(vdate))%>%
  select(-gender, -GAcomputed, -SGA, -date, -dob)

nulisa_plus_neuro <- left_join(clean_data, neuro_cog_edit, by="id")


big_palette_hex <- c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                     "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                     "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                     "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                     "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                     "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                     "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                     "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                     "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                     "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                     "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                     "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                     "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")

# big_palette <- c("240,163,255","0,117,220","153,63,0","76,0,92","25,25,25","0,92,49","43,206,72","255,204,153","128,128,128","148,255,181","143,124,0","157,204,0","194,0,136","0,51,128","255,164,5","255,168,187","66,102,0","255,0,16","94,241,242","0,153,143","224,255,102","116,10,255","153,0,0","255,255,128","255,255,0","255,80,5")
# 
# big_palette_hex <- sapply(strsplit(big_palette, ","), function(x)
#   rgb(x[1], x[2], x[3], maxColorValue=255))


id_columns <- c("sample", "qc_failed","plate", "id", "gender", "ageinwks", "timepoint", "log_qpcr",
"total_n_malaria_12", "mstatus", "febrile", "fever", "rogerson", "SGA", "cog_scale", "rlan_scale",
"elan_scaled", "fmot_scaled", "gmot_scaled", "lang_composite","mot_composite", 
"cog_composite", "weeks_age", "mother_id", "dob", "AGE", "gravidcat", "LBWdich", 
"preterm", "anyHP", "totalvisits", "totalfever","everfever", "totalmalaria","evermalaria", "BMIenrol","Olevel")

wide_df2 <- nulisa_plus_neuro %>%
  filter(targetName %notin% c("IFNA2", "CTSS", "LTA|LTB"))%>%
  # filter(sample_id %notin% c("X384_S.1_D1V93U", "X744_A.1_D1KT5Z", "X323_A14_D1DZ2D", "X496_NM7_D1CAYS", "X316_S.1_D12FNR", "X667_NM7_D1GBYA", "X 176 S_t14 D1FXRJ", "X164_NM0_D1WA6Q"))%>%
  pivot_wider(names_from = targetName, values_from = conc, id_cols = all_of(id_columns))


rownames(wide_df2) <- wide_df2$sample

# each row = measurement; each column = feature
big_pca <-  prcomp(wide_df2[,(length(id_columns)+1):ncol(wide_df2)], center = T)
pca_plot_data <- as.data.frame(cbind(wide_df2, big_pca$x))

# [1] "sample"             "qc_failed"          "plate"              "id"                 "gender"             "ageinwks"           "timepoint"         
# [8] "log_qpcr"           "total_n_malaria_12" "mstatus"            "febrile"            "fever"


set.seed(123)

large_pca12 <- pca_plot_data %>%
  mutate(pcr_cat=round(log_qpcr))%>%
  # mutate(qc_failed=as.factor(qc_failed))%>%
  pivot_longer(cols=c(id_columns[-c(1, 4, 7, 8)], pcr_cat), names_to = "meta_column", values_to = "meta_value", values_transform = as.character)%>%
  arrange(meta_value)%>%
  ggplot(., aes(x=PC1, y=PC2, color=factor(meta_value)))+
  geom_point()+
  stat_ellipse()+
  facet_wrap(~meta_column)+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=(sample(big_palette_hex, size = 48)))+
  # viridis::scale_color_viridis(discrete = T)+
  theme(legend.position="none",
        axis.text = element_blank())



# linear regression ####
neuro_cog_purff52 <- nulisa_plus_neuro %>%
  filter(timepoint=="52 weeks")%>%
  filter(targetName %notin% c("IFNA2", "CTSS", "LTA|LTB"))%>%
  mutate(gender_categorical=if_else(gender==1, "male", "female"))%>%
  group_by(targetName)%>%
  nest() %>%
  mutate(cog_model=map(data, ~lm(conc~cog_composite+gender_categorical, data=., na))) %>%
  mutate(mot_model=map(data, ~lm(conc~mot_composite+gender_categorical, data=.))) %>%
  mutate(lang_model=map(data, ~lm(conc~lang_composite+gender_categorical, data=.))) %>%
  mutate(cog_model_summary=map(cog_model, ~summary(.))) %>%
  mutate(mot_model_summary=map(mot_model, ~summary(.))) %>%
  mutate(lang_model_summary=map(lang_model, ~summary(.))) %>%
  mutate(cog_model_summary_p=map(cog_model_summary, ~coef(.)[11])) %>%
  mutate(mot_model_summary_p=map(mot_model_summary, ~coef(.)[11])) %>%
  mutate(lang_model_summary_p=map(lang_model_summary, ~coef(.)[11])) %>%
  pivot_longer(cols=ends_with("_p"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))


neuro_cog_plot <- nulisa_plus_neuro%>%
  filter(timepoint=="52 weeks")%>%
  filter(targetName %notin% c("IFNA2", "CTSS", "LTA|LTB"))%>%
  mutate(gender_categorical=if_else(gender==1, "male", "female"))%>%
  # select(id, targetName, conc, timepoint,  GAcomputed, SGA, preterm, lang_composite, mot_composite, cog_composite, totalfever, totalmalaria, rogerson)%>%
  pivot_longer(cols = ends_with("composite"), names_to = "composite_kind", values_to = "composite_score")%>%
  ggplot(., aes(x=factor(gender_categorical), y=conc, fill=factor(gender_categorical)))+
  geom_violin(draw_quantiles = seq(0, 1, by=0.25))+
  facet_wrap(~composite_kind+timepoint)


nulisa_plus_neuro%>%
  filter(targetName %in% c("IL3RA", "TNFRSF14", "IFNB1", "SDC1"))%>%
  mutate(gender_categorical=if_else(gender==1, "male", "female"))%>%
  mutate(timepoint = factor(timepoint, levels = c("8 weeks", "24 weeks", "52 weeks")))%>%
  ggplot(., aes(x=timepoint, y=conc, fill=factor(gender_categorical)))+
  geom_violin(draw_quantiles = seq(0, 1, by=0.25))+
  facet_wrap(~targetName)

nulisa_plus_neuro%>%
  filter(targetName %in% c("ANGPT1", "BDNF", "PDGFA", "TGFB1",  "PDGFB"))%>%
  mutate(gender_categorical=if_else(gender==1, "male", "female"))%>%
  mutate(timepoint = factor(timepoint, levels = c("8 weeks", "24 weeks", "52 weeks")))%>%
  pivot_longer(cols = ends_with("composite"), names_to = "composite_kind", values_to = "composite_score")%>%
  ggplot(., aes(x=factor(composite_score), y=conc))+
  geom_violin(draw_quantiles = seq(0,1,0.25))+
  # ggpubr::stat_cor(method="spearman")+
  facet_wrap(~targetName+composite_kind, scales="free")+
  theme_minimal()



# maternal nulisa 
maternal_data <- read.csv("~/postdoc/stanford/plasma_analytes/DPSP/dpsp_nulisa_data.csv")

neuro_cog <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/NCT_infants_share032025.csv")
maternal_neuro_cog_edit <- neuro_cog%>%
  mutate(id=subjid-10000)

nulisa_plus_neuro <- left_join(maternal_data, maternal_neuro_cog_edit, by="id")



# other modelling approach ####

infant_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")

neuro_cog <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/NCT_infants_share032025.csv")

maternal_neuro_cog_edit <- neuro_cog%>%
  mutate(id=subjid)

nulisa_plus_neuro <- left_join(infant_data, maternal_neuro_cog_edit, by="id")
## complex modelling ####
linear_lms <- nulisa_plus_neuro %>%
  filter(targetName !="LTA|LTB", mstatus!=1, timepoint %notin% c("68 weeks"), sample%notin%c("11660_tp24_Repeat"))%>%
  mutate(id_cat=paste("id_", id, sep=""))%>%
  pivot_longer(cols = ends_with("composite"), names_to = "composite_kind", values_to = "composite_score")%>%
  # mutate(comp_score_7 = composite_score/7)%>%
  # pivot_wider(names_from = timepoint, values_from = conc, id_cols = c(id_cat, targetName, composite_score, composite_kind))%>%
  group_by(targetName, composite_kind)%>%
  nest() %>%
  mutate(model=map(data, ~lme4::lmer(conc~composite_score*timepoint+(1|id_cat), data=.))) %>%
  mutate(summary=map(model, ~summary(.)))%>%
  #The output will show the effect (slope) of composite_score on conc at each of the levels of timepoint.
  mutate(emm=map(model, ~emmeans(., ~ composite_score * timepoint)))%>%
  #The output will show the average values of conc at each timepoint for a given level of composite_score
  # mutate(emm2=map(model, ~emmeans(., ~ timepoint | composite_score)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "revpairwise")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate("52 weeks - 24 weeks p"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("slope difference 52 weeks - 24 weeks"=map_dbl(emm_contrast_summary, ~.$estimate[1])) %>%
  mutate("8 weeks - 24 weeks p"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("slope difference 8 weeks - 24 weeks"=map_dbl(emm_contrast_summary, ~.$estimate[2])) %>%
  mutate("8 weeks - 52 weeks p"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  mutate("slope difference 8 weeks - 52 weeks"=map_dbl(emm_contrast_summary, ~.$estimate[3])) %>%
  pivot_longer(cols=ends_with("p"), names_to = "contrast", values_to = "raw_p")%>%
  group_by(contrast, composite_kind)%>%
  mutate("padj" = p.adjust(raw_p, method="fdr"))

linear_lms%>%
  select(targetName, starts_with("slope"), composite_kind, raw_p, padj)%>%
  arrange(padj)
# p <- nulisa_plus_neuro %>%
#   filter(targetName %in%c("PDGFB", "AGER", "ANGPT2", "IL15", "IL10", "TLR3"), mstatus!=1, timepoint %notin% c("68 weeks"), sample%notin%c("11660_tp24_Repeat"))%>%
#   mutate(id_cat=paste("id_", id, sep=""))%>%
#   pivot_longer(cols = ends_with("composite"), names_to = "composite_kind", values_to = "composite_score")%>%
#   # pivot_wider(names_from = timepoint, values_from = conc, id_cols = c(id_cat, targetName, composite_score, composite_kind))%>%
#   mutate(comp_score_7 = composite_score/7)%>%
#   ggplot(., aes(x=comp_score_7, y=conc))+
#   geom_point()+
#   geom_smooth()+
#   geom_rug(alpha=0.2)+
#   facet_wrap(~targetName)+
#   theme_minimal()
# 
# ggExtra::ggMarginal(p)  

nulisa_plus_neuro %>%
  filter(targetName %in%c("FGF21"), mstatus!=1, timepoint %notin% c("68 weeks"), sample%notin%c("11660_tp24_Repeat"))%>%
  mutate(id_cat=paste("id_", id, sep=""))%>%
  pivot_longer(cols = ends_with("composite"), names_to = "composite_kind", values_to = "composite_score")%>%
  # pivot_wider(names_from = timepoint, values_from = conc, id_cols = c(id_cat, targetName, composite_score, composite_kind))%>%
  mutate(comp_score_7 = composite_score/7)%>%
  ggplot(., aes(x=comp_score_7, y=conc))+
  geom_point()+
  geom_smooth()+
  geom_rug(alpha=0.2)+
  facet_wrap(~targetName+timepoint+composite_kind)+
  theme_minimal()

## simpler modelling ####

simple_linear_lms <- nulisa_plus_neuro %>%
  filter(targetName !="LTA|LTB", mstatus!=1, timepoint %in% c("52 weeks"))%>%
  pivot_longer(cols = ends_with("composite"), names_to = "composite_kind", values_to = "composite_score")%>%
  mutate(Olevel_categorical = factor(Olevel))%>%
  # pivot_wider(names_from = timepoint, values_from = conc, id_cols = c(id_cat, targetName, composite_score, composite_kind))%>%
  group_by(targetName, composite_kind)%>%
  nest() %>%
  mutate(model=map(data, ~lm(conc~composite_score+gender_categorical+Olevel_categorical+AGE, data=.))) %>%
  mutate(summary=map(model, ~summary(.)))%>%
  mutate(raw_p=map_dbl(summary, ~coef(.)[17]))%>%
  group_by(composite_kind)%>%
  mutate("padj" = p.adjust(raw_p, method="fdr"))

simple_linear_lms%>%
  select(targetName, composite_kind, raw_p, padj)%>%
  arrange(padj)


