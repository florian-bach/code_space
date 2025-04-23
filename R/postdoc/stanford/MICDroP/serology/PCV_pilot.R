library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)



flexift_dir <- "/Users/fbach/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/PCV experiments/standardization//"

# pcv_data <- read.csv("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/PCV experiments/"
alyssa_file_path <- '~/'

library(data.table)

plate_list <- c(1:8)
serotype_list <- c('1', '14', '18C', '19F', '23F', '4', '5', '6B', '7F', '9V')

df_all_plate_data <- data.table()
for(p in plate_list){
  message('working on plate ', p)
  for(s in serotype_list){
    
    df_temp <- fread(paste0(flexift_dir,'Plate',p,'_serotype',s,'_flexfit_conc.csv'))
    
    df_temp$plate <- p
    df_temp$serotype <- s
    df_all_plate_data <- rbind(df_all_plate_data, df_temp)
  }
}

df_all_plate_data$id_clean <- sub(" .*", "", df_all_plate_data$id)
df_all_plate_data$age <- as.numeric(sub("^[^ ]+ wk", "", df_all_plate_data$id))

standard_removed <- df_all_plate_data%>%
  filter(!grepl("^STD", id))


nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")

epi_data <- nulisa_data%>%
  distinct(id, gender_categorical, total_n_para_12, total_n_malaria_12, total_n_malaria_6, total_n_para_6)%>%
  mutate(id_clean=as.character(id))

pcv_and_nulisa <- left_join(standard_removed, epi_data, by=c("id_clean"), relationship = "many-to-one")%>%
  mutate(any_malaria=ifelse(total_n_malaria_12==0, 0, 1))%>%
  mutate(any_para=ifelse(total_n_para_12==0, 0, 1))
  



#nothing
pcv_purrf <- pcv_and_nulisa %>%
  mutate(any_malaria=ifelse(total_n_malaria_12==0, 0, 1))%>%
  mutate(any_para=ifelse(total_n_para_12==0, 0, 1))%>%
  filter(!is.na(age))%>%
  group_by(age, serotype)%>%
  nest() %>%
  mutate(malaria_bino_model=map(data, ~glm(any_malaria~concentration+gender_categorical, data=., family = "binomial"))) %>%
  mutate(para_bino_model=map(data, ~glm(any_para~concentration+gender_categorical, data=., family = "binomial"))) %>%
  mutate(malaria_bino_model_summary=map(malaria_bino_model, ~summary(.))) %>%
  mutate(para_bino_model_summary=map(para_bino_model, ~summary(.))) %>%
  mutate(malaria_bino_model_summary_p=map_dbl(malaria_bino_model_summary, ~coef(.)[11])) %>%
  mutate(para_bino_model_summary_p=map_dbl(para_bino_model_summary, ~coef(.)[11])) %>%
  pivot_longer(cols=ends_with("_p"), names_to = "contrast", values_to = "p")%>%
  group_by(contrast)%>%
  mutate(padj = p.adjust(p, method="fdr"))


#nothing for 6 or 12; some spurious relationships because the few individuals with >1 malaria episode score high
nb_lms <- pcv_and_nulisa %>%
  filter(!is.na(age))%>%
  group_by(age, serotype)%>%
  nest() %>%
  mutate(malaria_model=map(data, ~MASS::glm.nb(total_n_malaria_12~concentration+gender_categorical, data=.))) %>%
  mutate(para_model=map(data, ~MASS::glm.nb(total_n_para_12~concentration+gender_categorical, data=.))) %>%
  mutate(malaria_summary=map(malaria_model, ~summary(.))) %>%
  mutate(para_summary=map(para_model, ~summary(.))) %>%
  mutate(malaria_p=map_dbl(malaria_summary, ~coef(.)[11]),
         para_p=map_dbl(para_summary, ~coef(.)[11]))%>%
  pivot_longer(cols=ends_with("_p"), names_to = "contrast", values_to = "p")%>%
  group_by(age)%>%
  mutate(padj=p.adjust(p, method="fdr"))



pcv_and_nulisa%>%
  filter(age==52, serotype=="6B")%>%
  ggplot(., aes(x=factor(total_n_malaria_12), y=concentration))+
  geom_boxplot()+
  geom_point()+
  scale_y_continuous(trans="log2")+
  facet_wrap(~serotype, scales = "free")+
  theme_minimal()

pcv_and_nulisa%>%
  filter(age==24)%>%
  ggplot(., aes(x=factor(total_n_para_12), y=concentration))+
  geom_boxplot()+
  scale_y_continuous(trans="log2")+
  facet_wrap(~serotype, scales = "free")+
  theme_minimal()

any_malar_plot <- pcv_and_nulisa%>%
  filter(age==52)%>%
  ggplot(., aes(x=factor(any_malaria), y=concentration+0.0001))+
  geom_boxplot(aes(fill=serotype))+
  scale_y_log10()+
  facet_wrap(~serotype, scales = "free", nrow=2)+
  ylab("Concentration (AU)")+
  xlab("any malaria")+
  ggpubr::stat_compare_means(vjust=0.5, size=3)+
  viridis::scale_fill_viridis(discrete = T)+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/Downloads/any_malar_plot.png", any_malar_plot, width=8, height=4, dpi=444, bg="white")

any_para_plot <- pcv_and_nulisa%>%
  filter(age==52)%>%
  ggplot(., aes(x=factor(any_para), y=concentration+0.0001))+
  geom_boxplot(outliers=F, aes(fill=serotype))+
  ggpubr::stat_compare_means(vjust=0.5, size=3)+
  scale_y_log10()+
  facet_wrap(~serotype, scales = "free", nrow=2)+
  ylab("Concentration (AU)")+
  xlab("any parasitemia")+
  viridis::scale_fill_viridis(discrete = T)+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/Downloads/any_para_plot.png", any_para_plot, width=8, height=4, dpi=444, bg="white")

