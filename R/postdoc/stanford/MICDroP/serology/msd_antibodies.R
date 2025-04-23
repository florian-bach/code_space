library(tidyr)
library(dplyr)
library(xlsx)
library(ggplot2)
library(purrr)

nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")

metadata_columns <- c("id", "dob", "date", "ageinwks", "gender_categorical", "mstatus", "qPCRparsdens", "visittype", "fever", "febrile", "rogerson", "GAcomputed", "gi", "SGA", "qPCRdich", "mqPCRparsdens")

epi_data <- nulisa_data%>%
  distinct(sample, total_n_para_12, total_n_malaria_12, total_n_malaria_6, total_n_para_6,
           id, dob, date, ageinwks, gender_categorical, mstatus, qPCRparsdens, visittype, fever, febrile, rogerson, GAcomputed, gi, SGA, qPCRdich, mqPCRparsdens, anyHP)

msd_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/batch_one.csv")

long_msd <- msd_data%>%
  mutate(sample=paste(SubjectID, "_", "tp", TimePt, sep=""))%>%
  mutate(id=SubjectID, timepoint=paste(TimePt, "weeks"))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))%>%
  select(-SubjectID, -TimePt)%>%
  pivot_longer(cols=-c(sample, id, timepoint), names_to = "antigen", values_to = "titer")

antibodies_and_epi <- epi_data%>%
  inner_join(., long_msd, by="sample")


antibodies_and_epi%>%
  filter(timepoint=="24 weeks", !is.na(anyHP))%>%
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
  ggpubr::stat_compare_means()+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()

antibodies_and_epi%>%
  filter(timepoint=="52 weeks")%>%
  mutate(any_para = if_else(total_n_para_12 > 0, 1, 0),
         any_malaria = if_else(total_n_malaria_12 > 0, 1, 0))%>%
  filter(antigen %in% c("RSV", "Pneumo.1.4.14", "PIV.1", "PIV.2"))%>%
  ggplot(., aes(x=factor(any_para), y=titer))+
  geom_boxplot()+
  scale_y_log10()+
  ggpubr::stat_compare_means()+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()

antibodies_and_nulisa <-  nulisa_data%>%
  select(-id, -timepoint)%>%
  inner_join(., long_msd, by="sample")

big_purrf <- antibodies_and_nulisa%>%
  group_by(antigen, targetName)%>%
  nest()%>%
  mutate(correlation=map(data, ~cor.test(.$conc, .$titer, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%
  ungroup()%>%
  mutate(padj=p.adjust(p, method="fdr"))

top_100 <- big_purrf%>%
  arrange(p)%>%
  slice_head(n = 100)

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
