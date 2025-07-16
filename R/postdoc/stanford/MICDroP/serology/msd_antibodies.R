library(tidyr)
library(dplyr)
library(xlsx)
library(ggplot2)
library(purrr)
library(data.table)
library(mclust)
library(emmeans)

nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")
mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")
maternal_treatment_arms <- haven::read_dta("~/Library/CloudStorage/Box-Box/DP+SP study/Databases and preliminary findings/Final database used for analyses/DPSP treatment allocation_FINAL.dta")

vaccines = c("Diptheria",     "Measles" ,      "Mumps",         "Pertussis",     "Polio",
             "Rotavirus" ,    "Rubella",       "Tetanus", "Pneumo.1.4.14")


nulisa_data <- nulisa_data%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))

metadata_columns <- c("id", "anyDP", "treatmentarm",  "dob", "date", "ageinwks", "gender_categorical", "mstatus", "qPCRparsdens", "visittype", "fever", "febrile", "rogerson", "GAcomputed", "gi", "SGA", "qPCRdich", "mqPCRparsdens")

epi_data <- nulisa_data%>%
  distinct(sample, total_n_para_12, total_n_malaria_12, total_n_malaria_6, total_n_para_6,
           id, dob, date, ageinwks, gender_categorical, mstatus, qPCRparsdens, visittype, fever, febrile, rogerson, GAcomputed, gi, SGA, qPCRdich, mqPCRparsdens, anyHP)%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"),
         mom_rx=maternal_treatment_arms$treatmentarm[match(id-10000, maternal_treatment_arms$id)],
         mom_rx=case_match(mom_rx,
                           1~"SP",
                           2~"DP",
                           3~"DPSP"))

msd_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/batch_one.csv")

long_msd <- msd_data%>%
  mutate(sample=paste(SubjectID, "_", "tp", TimePt, sep=""))%>%
  mutate(id=SubjectID, timepoint=paste(TimePt, "weeks"))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))%>%
  select(-SubjectID, -TimePt)%>%
  pivot_longer(cols=-c(sample, id, timepoint), names_to = "antigen", values_to = "titer")

antibodies_and_epi <- epi_data%>%
  inner_join(., long_msd, by="sample")

antibodies_and_nulisa <-  nulisa_data%>%
  select(-id, -timepoint)%>%
  inner_join(., long_msd, by="sample")


antibodies_and_epi%>%
  filter(timepoint=="8 weeks", !is.na(anyHP))%>%
  mutate(any_para = if_else(total_n_para_6 > 0, 1, 0),
         any_malaria = if_else(total_n_malaria_6 > 0, 1, 0))%>%
  ggplot(., aes(x=factor(anyHP), y=titer))+
  geom_boxplot(outliers = F)+
  ggpubr::stat_compare_means(vjust = 0.5)+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()


antibodies_and_epi%>%
  filter(timepoint=="24 weeks")%>%
  mutate(any_para = if_else(total_n_para_6 > 0, 1, 0),
         any_malaria = if_else(total_n_malaria_6 > 0, 1, 0))%>%
  filter(antigen %in% c("RSV", "Pneumo.1.4.14", "PIV.1", "PIV.2"))%>%
  ggplot(., aes(x=factor(any_para), y=titer))+
  geom_boxplot()+
  scale_y_log10()+
  ggpubr::stat_compare_means()+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()


w52_plopt <- antibodies_and_epi%>%
  filter(timepoint=="52 weeks")%>%
  mutate(any_para = if_else(total_n_para_12 > 0, 1, 0),
         any_malaria = if_else(total_n_malaria_12 > 0, 1, 0))%>%
  filter(antigen %in% c("RSV", "Pneumo.1.4.14", "PIV.1", "PIV.2"))%>%
  ggplot(., aes(x=factor(any_para), y=titer, fill=factor(any_para)))+
  geom_boxplot()+
  scale_y_log10()+
  xlab("any parasitemia")+
  ylab("")+
  scale_fill_manual(values=c("grey", "darkred"))+
  ggpubr::stat_compare_means()+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()+
  theme(legend.position = "none")


ggsave("~/Downloads/52weekplot.png", w52_plopt, height=6, width = 6, dpi=444, bg="white")


big_purrf <- antibodies_and_nulisa%>%
  group_by(antigen, targetName, timepoint)%>%
  nest()%>%
  mutate(correlation=map(data, ~cor.test(.$conc, .$titer, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%
  group_by(antigen, timepoint)%>%
  mutate(padj=p.adjust(p, method="fdr"))

top_100 <- big_purrf%>%
  filter(padj<0.1)

big_purrf2 <- antibodies_and_nulisa%>%
  mutate(log_titer=log10(titer))%>%
  group_by(antigen, targetName, timepoint)%>%
  nest()%>%
  mutate(model=map(data, ~lm(conc~log_titer+gender_categorical, data=.)))%>%
  mutate(summary=map(model, ~summary(.)))%>%
  mutate(p=map_dbl(summary, ~coef(.)[11]),
         coef=map_dbl(summary, ~coef(.)[2]))%>%
  group_by(antigen, timepoint)%>%
  mutate(padj=p.adjust(p, method="fdr"))

topn <- big_purrf2%>%
  filter(padj<0.1)

# most common hits; probably not significant
# 1 CXCL12    
# 2 IFNL1     
# 3 VEGFA 

sigs_only <- big_purrf%>%
  filter(padj<0.05)
# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  # mutate(simple_lm=map(data, ~glm(conc~titer, data=.)))%>%
  # mutate(simple_lm_summary=map(simple_lm, ~summary(.))) %>%
  # mutate(simple_lm_summary_p=map_dbl(simple_lm_summary, ~coef(.)[11])) %>%
  # mutate(para_bino_model_summary_p=map_dbl(para_bino_model_summary, ~coef(.)[11])) %>%
  # pivot_longer(cols=ends_with("_p"), names_to = "contrast", values_to = "p")%>%
  # group_by(contrast)%>%
  # mutate(padj = p.adjust(p, method="fdr"))


big_plot <- antibodies_and_nulisa%>%
  mutate(any_para = if_else(total_n_para_12 > 0, 1, 0),
         any_malaria = if_else(total_n_malaria_12 > 0, 1, 0))%>%
  # filter(antigen %in% c("RSV", "Pneumo.1.4.14", "PIV.1", "PIV.2"))%>%
  distinct(sample, timepoint, titer, any_para, antigen)%>%
  ggplot(., aes(x=timepoint, y=titer, fill=factor(any_para)))+
  geom_boxplot(outliers = F)+
  scale_y_log10()+
  ggpubr::stat_compare_means()+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/big_any_para.png", big_plot, width=20, height = 10, dpi=444, bg="white")



# fold change analysis ####

antibody_fc <- antibodies_and_nulisa %>%
  distinct(id, timepoint, antigen, titer)%>%
  pivot_wider(names_from = timepoint, values_from = titer, id_cols = c("id", "antigen"))%>%
  mutate("8_to_24_fc"=`24 weeks`/`8 weeks`,
         "8_to_52_fc"=`52 weeks`/`8 weeks`,
         "24_to_52_fc"=`52 weeks`/`24 weeks`)%>%
  select(id, antigen, "8_to_24_fc", "8_to_52_fc", "24_to_52_fc")


nulisa_with_antibody_fc <- nulisa_data%>%
  inner_join(., antibody_fc, by="id")

slim_nulisa_with_antibody_fc <- nulisa_with_antibody_fc%>%
  select(id, timepoint, targetName, conc, antigen, ends_with("_fc"))

slim_nulisa_with_antibody_fc_purrf <- slim_nulisa_with_antibody_fc %>%
  group_by(antigen, targetName, timepoint)%>%
  nest()%>%
  mutate(correlation_8_24=map(data, ~cor.test(.$conc, .$`8_to_24_fc`, method = "spearman")))%>%
  mutate(correlation_8_52=map(data, ~cor.test(.$conc, .$`8_to_52_fc`, method = "spearman")))%>%
  mutate(correlation_24_52=map(data, ~cor.test(.$conc, .$`24_to_52_fc`, method = "spearman")))%>%
  mutate(p_8_24=map_dbl(correlation_8_24, ~.$p.value),
         rho_8_24=map_dbl(correlation_8_24, ~.$estimate))%>%
  mutate(p_8_52=map_dbl(correlation_8_52, ~.$p.value),
         rho_8_52=map_dbl(correlation_8_52, ~.$estimate))%>%
  mutate(p_24_52=map_dbl(correlation_24_52, ~.$p.value),
         rho_24_52=map_dbl(correlation_24_52, ~.$estimate))%>%
  pivot_longer(cols=starts_with("p_"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast, timepoint)%>%
  mutate(padj = p.adjust(p, method="fdr"))

kinda_sigs <- filter(slim_nulisa_with_antibody_fc_purrf, padj<0.1)

slim_nulisa_with_antibody_fc%>%
  filter(targetName %in% kinda_sigs$targetName,
         antigen %in% "Flu.A.H3",
         timepoint == "8 weeks")%>%
  distinct(id, targetName, antigen, conc, `8_to_24_fc`)%>%
  ggplot(., aes(x=log2(`8_to_24_fc`), y=conc))+
  geom_point()+
  scale_x_continuous()+
  geom_smooth(method="lm")+
  xlab("fold change in titer from week 8 to 24")+
  ggpubr::stat_cor(method="spearman")+
  theme_minimal()+
  facet_wrap(~targetName+antigen
  )

slim_nulisa_with_antibody_fc%>%
  filter(targetName %in% c("BST2", "CCL7", "CTLA4", "GFAP", "HAVCR1", "MMP12", "TNFSF8"),
         antigen %in% c("Flu.A.H1", "Flu.A.H3"),
         timepoint == "8 weeks")%>%
  distinct(id, targetName, antigen, conc, `8_to_24_fc`)%>%
  ggplot(., aes(x=log2(`8_to_24_fc`), y=conc))+
  geom_point()+
  scale_x_continuous()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  theme_minimal()+
  facet_wrap(~targetName+antigen
  )



slim_nulisa_with_antibody_fc%>%
  filter(targetName %in% c("CXCL12", "IL27", "IL2RB", "TSLP", "LTA"),
         antigen %in% c("Measles"),
         timepoint == "24 weeks")%>%
  distinct(id, targetName, antigen, conc, `24_to_52_fc`)%>%
  ggplot(., aes(x=log2(`24_to_52_fc`), y=conc))+
  geom_point()+
  scale_x_continuous()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  theme_minimal()+
  facet_wrap(~targetName+antigen)+
  xlab("log 2 fold change in antibody from week 24 to 52")+
  ylab("plasma protein at 24 weeks")

long_msd%>%
  filter(antigen %in% c("Flu.A.H1", "Flu.A.H3"))%>%
  ggplot(., aes(x=factor(timepoint), y=titer, fill=factor(timepoint)))+
  geom_point()+
  geom_line(aes(group=id))+
  geom_violin()+
  scale_y_continuous(trans='log10')+
  theme_minimal()+
  facet_wrap(~antigen, scales="free")



slim_nulisa_with_antibody_fc_purrf%>%
  filter(contrast != "p_8_24", p<0.01)%>%
  group_by(antigen, contrast)%>%
  summarise(n=n())%>%
  arrange(desc(n))

slim_nulisa_with_antibody_fc_purrf%>%
  filter(contrast != "p_8_24", p<0.01)%>%
  group_by(targetName, contrast)%>%
  summarise(n=n())%>%
  arrange(desc(n))
 

# linear regression ####


treatment_purf <- antibodies_and_epi%>%
  filter(mstatus==0, treatmentarm!="DP 2 years")%>%
  mutate(id=factor(id.x))%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(titer~timepoint*treatmentarm+(1|id), data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ treatmentarm | timepoint)))%>%
  mutate(emm2=map(time_model, ~emmeans(., specs = pairwise ~ timepoint | treatmentarm)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  mutate("8 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("24 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("52 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  pivot_longer(cols=ends_with("weeks"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))


sigs <- treatment_purf%>%
  filter(padj<0.15)

# (time_treatment_plot <- antibodies_and_epi%>%
#   filter(treatmentarm!="DP 2 years", timepoint=="52 weeks")%>%
#   filter(antigen %in% sigs$antigen)%>%
#   ggplot(aes(x=factor(timepoint), y=titer, color=treatmentarm))+
#   geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1))+
#   scale_y_log10()+
#   # stat_summary(aes(group=c(treatmentarm)), geom="cro", fun.y=median, colour="white")+
#   stat_summary(aes(group=treatmentarm), fun.y = median, fun.ymin = median, fun.ymax = median,
#                geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="black")+
#   ggpubr::stat_compare_means()+
#   facet_wrap(~antigen, scales="free")+
#   scale_fill_manual(values=c("darkred", "#00555A"))+
#   theme_minimal()+
#   theme(legend.title = element_blank(),
#         axis.title = element_blank()))

(time_treatment_plot <- antibodies_and_epi%>%
    filter(treatmentarm!="DP 2 years")%>%
    filter(antigen %in% sigs$antigen)%>%
    ggplot(aes(x=factor(timepoint), y=titer, fill=treatmentarm))+
    geom_boxplot(outliers = F)+
    scale_y_log10()+
    # stat_summary(aes(group=c(treatmentarm)), geom="cro", fun.y=median, colour="white")+
    stat_summary(aes(group=treatmentarm), fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
    facet_wrap(~antigen, scales="free", nrow=1)+
    scale_fill_manual(values=c("darkred", "#00555A"))+
    theme_minimal()+
    theme(legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/titer_time_treatment_plot.png", width=9, height=3, dpi=444, bg="white")

# sandbox ####
# 
# targets <- unique(big_purrf$targetName)
# 
# results <- data.frame(matrix(ncol = 3))
# colnames(results) <- c("Var1","Freq", "test_iter")
# 
# filler <-  data.frame(matrix(ncol = 3))
# colnames(filler) <- c("Var1","Freq", "test_iter")
# 
# filler_item=NULL
# 
# for(i in 1:1000){
#   
#   try <- sample(targets, size = 1002, replace = T)
#   item <- data.frame(table(table(try)), "test_iter"=i)
#   
#   filler_item <- data.frame(Var1=factor((nrow(item)+1):20), Freq=0, test_iter=i)
#   filled_item <- rbind(item, filler_item)
#   
#   results <- rbind(results, filled_item)
#   
# }
# 
# results <- results[2:nrow(results),]
# results$Var1 <- as.numeric(results$Var1)
# 
# results%>%
#   group_by(Var1)%>%
#   summarise(mean(Freq))
# 
# sigs_only%>%
#   group_by(targetName)%>%
#   add_count()%>%
#   filter(n>=12)%>%
#   distinct(targetName)


#slides for George PCA ####
id_columns <- c("sample", "timepoint", "gender_categorical", "ageinwks", "log_qpcr","mstatus", "febrile","fever",
                "total_n_malaria_12", "total_n_para_12",  "total_n_malaria_6", "total_n_para_6",
                "treatmentarm", "anyDP", "hbs", 
                "rogerson", "SGA", "anyHP")

wide_df2 <- long_msd %>%
  mutate(log_titer=log10(titer))%>%
  distinct(sample, antigen, log_titer)%>%
  # filter(sample_id %notin% c("X384_S.1_D1V93U", "X744_A.1_D1KT5Z", "X323_A14_D1DZ2D", "X496_NM7_D1CAYS", "X316_S.1_D12FNR", "X667_NM7_D1GBYA", "X 176 S_t14 D1FXRJ", "X164_NM0_D1WA6Q"))%>%
  pivot_wider(names_from = antigen, values_from = log_titer)%>%
  filter(sample!="NA_tpNA")

rownames(wide_df2) <- wide_df2$sample
id_data <- antibodies_and_nulisa%>%
  select(all_of(id_columns))%>%
  distinct()

wider_df2 <- left_join(wide_df2, id_data, by="sample")
# each row = measurement; each column = feature
big_pca <-  prcomp(wider_df2[,2:26], center = T)
pca_plot_data <- as.data.frame(cbind(wider_df2, big_pca$x))


pca_plot_data %>%
  ggplot(., aes(x=PC1, y=PC2, color=factor(hbs)))+
  geom_point()+
  stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  # ggrepel::geom_label_repel(aes(label=sample))+
  theme_minimal()+
  # geom_vline(xintercept = -2e+05)+
  # scale_color_viridis_d()+
  scale_fill_manual(values=c("darkred", "#00555A"))+
  facet_wrap(~timepoint)+
  theme(axis.text = element_blank(),
        legend.title=element_blank())

## PC loadings df ####
loadings_df <- data.frame(big_pca$rotation)
loadings_df$antibody <- rownames(loadings_df)
loadings_df$antibody <- factor(loadings_df$antibody, levels = loadings_df$antibody[order(loadings_df$PC1)])

pc1_cols <- colorspace::sequential_hcl(nrow(loadings_df), palette = "Purple Yellow")
names(pc1_cols) <- loadings_df$antibody[order(loadings_df$PC1)]

PC1_plot <- ggplot(loadings_df, aes(x=factor(antibody, levels = antibody[order(loadings_df$PC1)]), y=PC1, fill=antibody))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = pc1_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")

pc2_cols <- colorspace::sequential_hcl(nrow(loadings_df), palette = "Purple Yellow")
names(pc2_cols) <- loadings_df$antibody[order(loadings_df$PC2)]

PC2_plot <- ggplot(loadings_df, aes(x=factor(antibody, levels = antibody[order(loadings_df$PC2)]), y=PC2, fill=antibody))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = pc2_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")



# RSV, hMPV drive variance in both PC1 and PC 2
## outlier PCA ####
outliers <- filter(pca_plot_data, PC2< 0)


antibodies_and_epi%>%
  filter(antigen %in% c("Measles", "Rubella"))%>%
  mutate(outlier=ifelse(sample %in% outliers$sample & timepoint=="52 weeks", "yes", "no"))%>%
  arrange(outlier)%>%
  ggplot(., aes(x=factor(timepoint), y=titer))+
  geom_violin()+
  geom_point(aes(color=outlier), position = position_jitter(width=0.15))+
  scale_y_log10()+
  # ggrepel::geom_label_repel(aes(label=sample_label))+
  # ggpubr::stat_compare_means()+
  scale_color_manual(values=c("black", "red"))+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()+
  theme(axis.title=element_blank())


antibodies_and_epi%>%
  filter(antigen %in% vaccines)%>%
  ggplot(., aes(x=factor(timepoint), y=titer, fill=treatmentarm))+
  geom_violin()+
  scale_y_log10()+
  # ggrepel::geom_label_repel(aes(label=sample_label))+
  # ggpubr::stat_compare_means()+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()+
  theme(axis.title=element_blank())

# treatmentarm analysis ####
## child treatment arm ####

treatment_purf <- antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(id_cat=factor(id.x))%>%
  # filter(mstatus==0, treatmentarm!="DP 2 years")%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(log_titer~timepoint*treatmentarm+gender_categorical+(1|id_cat), data=.))) %>%
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
  pivot_longer(cols=ends_with("weeks"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

kinda_sigs <- treatment_purf %>%
  filter(padj<0.2)


antibodies_and_epi%>%
  filter(!antigen%in%vaccines)%>%
  ggplot(., aes(x=factor(timepoint), y=titer, fill=treatmentarm))+
  geom_point(position=position_dodge(width=0.75))+
  # geom_violin(draw_quantiles = c(seq(0,1,by=0.25)), color="white", )+
  geom_boxplot(outliers=F)+
  scale_y_log10()+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  
  # ggrepel::geom_label_repel(aes(label=sample_label))+
  ggpubr::stat_compare_means()+
  scale_fill_manual(values=c("darkred", "#00555A"))+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()+
  theme(axis.title=element_blank(),
        legend.title = element_blank(),
  )

## maternal treatment arm ####
mom_treatment_purf <- antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(id_cat=factor(id.x))%>%
  # filter(mstatus==0, treatmentarm!="DP 2 years")%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(log_titer~timepoint*mom_rx+(1|id_cat), data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ mom_rx | timepoint)))%>%
  mutate(emm2=map(time_model, ~emmeans(., specs = pairwise ~ timepoint | mom_rx)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  mutate("8 weeks DP - DPSP"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("8 weeks DP - SP"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("8 weeks DPSP - SP"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  mutate("24 weeks DP - DPSP"=map_dbl(emm_contrast_summary, ~.$p.value[4])) %>%
  mutate("24 weeks DP - SP"=map_dbl(emm_contrast_summary, ~.$p.value[5])) %>%
  mutate("24 weeks DPSP - SP"=map_dbl(emm_contrast_summary, ~.$p.value[6])) %>%
  mutate("52 weeks DP - DPSP"=map_dbl(emm_contrast_summary, ~.$p.value[7])) %>%
  mutate("52 weeks DP - SP"=map_dbl(emm_contrast_summary, ~.$p.value[8])) %>%
  mutate("52 weeks DPSP - SP"=map_dbl(emm_contrast_summary, ~.$p.value[9])) %>%
  pivot_longer(cols=ends_with("SP"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

kinda_sigs <- mom_treatment_purf %>%
  filter(padj<0.1)


antibodies_and_epi%>%
  # filter(antigen%in%vaccines)%>%
  ggplot(., aes(x=factor(timepoint), y=titer, fill=mom_rx))+
  # geom_point(position=position_dodge(width=0.75))+
  # geom_violin(draw_quantiles = c(seq(0,1,by=0.25)), color="white", )+
  geom_boxplot(outliers=F)+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  scale_y_log10()+
  facet_wrap(~antigen)+
  viridis::scale_fill_viridis(option = "cividis", discrete = T, direction = 1)+
  theme_minimal()
  
momrx_8week_plot <- antibodies_and_epi%>%
  filter(timepoint=="8 weeks")%>%
  ggplot(., aes(x=factor(timepoint), y=titer, fill=mom_rx))+
  # geom_point(position=position_dodge(width=0.75))+
  # geom_violin(draw_quantiles = c(seq(0,1,by=0.25)), color="white", )+
  geom_boxplot(outliers=F)+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  scale_y_log10()+
  facet_wrap(~antigen)+
  viridis::scale_fill_viridis(option = "cividis", discrete = T, direction = 1)+
  ggtitle("8 weeks")+
  theme_minimal()+
  theme(axis.title = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/momrx_8week_plot.png", momrx_8week_plot, width=8, height=8, dpi=444, bg="white")

## anyhp ####
anyhp_purf <- antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(id_cat=factor(id.x))%>%
  # filter(mstatus==0, treatmentarm!="DP 2 years")%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(log_titer~timepoint*anyHP+(1|id_cat), data=.))) %>%
  mutate(summary=map(time_model, ~summary(.))) %>%
  mutate(emm=map(time_model, ~emmeans(., specs = pairwise ~ anyHP | timepoint)))%>%
  mutate(emm2=map(time_model, ~emmeans(., specs = pairwise ~ timepoint | anyHP)))%>%
  mutate(emm_contrast=map(emm, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast2=map(emm2, ~contrast(., "pairwise", adjust="none")))%>%
  mutate(emm_contrast_summary=map(emm_contrast, ~summary(.)))%>%
  mutate(emm_contrast_summary2=map(emm_contrast2, ~summary(.)))%>%
  mutate("8 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[1])) %>%
  mutate("24 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[2])) %>%
  mutate("52 weeks"=map_dbl(emm_contrast_summary, ~.$p.value[3])) %>%
  pivot_longer(cols=ends_with("weeks"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

kinda_sigs <- anyhp_purf %>%
  filter(padj<0.1)

anyhp_8week_plot <- antibodies_and_epi%>%
  filter(timepoint=="8 weeks", !is.na(anyHP))%>%
  mutate(any_para = if_else(total_n_para_6 > 0, 1, 0),
         any_malaria = if_else(total_n_malaria_6 > 0, 1, 0))%>%
  ggplot(., aes(x=factor(anyHP), y=titer, fill=factor(anyHP)))+
  geom_violin(draw_quantiles = seq(0,1,0.25))+
  scale_y_log10()+
  facet_wrap(~antigen)+
  viridis::scale_fill_viridis(option = "cividis", discrete = T, direction = 1)+
  ggtitle("8 weeks")+
  theme_minimal()+
  theme(axis.title = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/anyhp_8week_plot.png", anyhp_8week_plot, width=8, height=8, dpi=444, bg="white")

## child treatment arm ####
treatment_purf <- antibodies_and_epi%>%
  mutate(log_titer=log10(titer))%>%
  mutate(id_cat=factor(id.x))%>%
  # filter(mstatus==0, treatmentarm!="DP 2 years")%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(time_model=map(data, ~lme4::lmer(log_titer~timepoint*treatmentarm+gender_categorical+(1|id_cat), data=.))) %>%
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
  pivot_longer(cols=ends_with("weeks"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))

kinda_sigs <- treatment_purf %>%
  filter(padj<0.2)


antibodies_and_epi%>%
  filter(!antigen%in%vaccines)%>%
  ggplot(., aes(x=factor(timepoint), y=titer, fill=treatmentarm))+
  geom_point(position=position_dodge(width=0.75))+
  # geom_violin(draw_quantiles = c(seq(0,1,by=0.25)), color="white", )+
  geom_boxplot(outliers=F)+
  scale_y_log10()+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", position = position_dodge(width = 0.75), width = 0.65, fatten=0.25, color="white")+
  
  # ggrepel::geom_label_repel(aes(label=sample_label))+
  # ggpubr::stat_compare_means()+
  scale_fill_manual(values=c("darkred", "#00555A"))+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()+
  theme(axis.title=element_blank(),
        legend.title = element_blank(),
  )

# high responder vs low responder


scaled_msd <- antibodies_and_epi %>%
  filter(sample!="NA_tpNA", antigen %in% vaccines)%>%
  group_by(antigen, timepoint)%>%
  mutate(scaled_titer=scales::rescale(titer, to = c(0, 1)))%>%
  mutate(quintile=case_when(scaled_titer<0.2~"lowest_quintile",
                            scaled_titer>0.2&scaled_titer<0.4~"second_quintile",
                            scaled_titer>0.4&scaled_titer<0.6~"third_quintile",
                            scaled_titer>0.6&scaled_titer<0.8~"fourth_quintile",
                            scaled_titer>0.8~"highest_quintile"))%>%
  mutate(quintile=factor(quintile, levels=c("lowest_quintile", "second_quintile", "third_quintile", "fourth_quintile", "highest_quintile")))

scaled_msd%>%
  ggplot(., aes(x=quintile, y=""))+
  geom_bar(stat="identity", aes(fill=treatmentarm))+
  facet_wrap(~antigen+timepoint)



# anitbody nulisa stuff

antibodies_and_nulisa%>%
  filter(grepl("^IFN*", targetName),
         antigen=="RV.C",
         timepoint=="52 weeks")%>%
  ggplot(., aes(x=titer, y=conc))+
  geom_point()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  facet_wrap(~targetName, scales="free")+
  theme_minimal()

antibodies_and_nulisa%>%
  filter(targetName %in% c("CCL26", "CCL28"),
         antigen%in%c("Flu.A.H3", "EV.71"), ,
         timepoint=="24 weeks")%>%
  ggplot(., aes(x=titer, y=conc))+
  geom_point()+
  scale_x_log10()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  facet_wrap(~targetName+antigen, scales="free")+
  theme_minimal()



antibodies_and_nulisa%>%
  filter(targetName %in% c("CSF3R", "CD40LG"),
         antigen%in%c("Flu.A.H3"), ,
         timepoint=="52 weeks")%>%
  ggplot(., aes(x=titer, y=conc))+
  geom_point()+
  scale_x_log10()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  facet_wrap(~targetName, scales="free")+
  theme_minimal()




# attempt at mixture models ####

# credit Kenneth & Alyssa for code below this line
timepoints = c("8 weeks", "24 weeks", "52 weeks")
  
seroprev_est_df <- data.table()

for(v in vaccines){
  for(t in timepoints){
  fmm_model=NULL
  concentration_values <- long_msd$titer[long_msd$antigen==v&long_msd$timepoint==t&!is.na(long_msd$titer)]
  log_serotype_vector <- log10(concentration_values)
  
  ### fit initial mixutre model with 2 components
  k <- 2 #number of components 
  fmm_model <- Mclust(log_serotype_vector, G = k, modelNames = "V") 
  
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
  
  # Plot histogram + density curves
  gg_density <- ggplot() +
    geom_histogram(aes(x = log_serotype_vector, y = after_stat(density)), 
                   bins = 50, fill = "gray80", color = "white") +
    geom_line(data = dens_df, aes(x = x, y = density, color = component), size = 1.2) +
    scale_x_continuous(limits = c(0, NA))+
    labs(title = "FMM: Mixture of 2 Components",
         subtitle=paste0("Serotype: ", v, " at ", t),
         x = "log(concentration)", y = "Density") +
    scale_color_manual(values = c("blue", "red")) +
    theme_minimal()+
    theme(legend.position = "none")
  
  ggsave(paste0("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/fmm_2components_",v, "_", t, ".png", sep=""),gg_density, width = 8, height = 4, dpi = 444, bg="white")
  
  # Step 1: Out of the 2 components, we identify which distribution has a higher mean, and assume that is the seropositive component
  component_means <- fmm_model$parameters$mean
  seropositive_component <- which.max(component_means)
  
  # Step 4: Get the posterior probabilities for each sample
  posteriors <- fmm_model$z  #For each sample, this is the probability of the component being in the first component or second component
  seropositive_probs <- posteriors[, seropositive_component] #this pulls out just the probability of each sample being seropositive
  
  # Step 5: Estimate seroprevalence
  estimated_seroprevalence <- mean(seropositive_probs) #the mean of all the probabilities will be equal to the overall population seroprevalence
  
  # Step 6: Report as percentage
  to_add <- cbind(v, estimated_seroprevalence, t)
  
  seroprev_est_df <- rbind(seroprev_est_df, to_add)
  
  }}

View(seroprev_est_df)
