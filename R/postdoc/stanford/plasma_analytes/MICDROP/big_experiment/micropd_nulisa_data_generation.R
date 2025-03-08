library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

`%notin%` <- Negate(`%in%`)


# read in pilots NULISA data ####

pilot_nulisa <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/extra_24/combined_pilot_data.csv")

# get data in same shape as big experiemnt

pilot_nulisa_edit <- pilot_nulisa%>%
  mutate("sample"=gsub("w", "", sample_id))%>%
  mutate("sample"=gsub("_", "_tp", sample))%>%
  mutate("file_name"="doesnt_matter")%>%
  mutate(timepoint2=stringr::str_split(sample, pattern = "_", simplify = T)[,2])%>%
  mutate(plate=study)%>%
  mutate(conc=concentration)%>%
  mutate(study="MICDROP")%>%
  select(targetName, sample, conc, file_name, id, timepoint2, study, plate)
  
  

# read in big NULISA data ####
excel_files <- list.files(path="~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/Data Files/", pattern = "*Report.xlsx", full.names = T)

big_data <- data.frame()
for(file in excel_files){
  
  sheet <- readxl::read_excel(file, sheet=1)
  
  long_sheet <- sheet %>%
    pivot_longer(cols=colnames(sheet)[-1], names_to = "sample", values_to = "conc")%>%
    mutate("file_name"=file)%>%
    mutate(id=stringr::str_split(sample, pattern = "_", simplify = T)[,1],
           timepoint2=stringr::str_split(sample, pattern = "_", simplify = T)[,2])%>%
    mutate(id=as.numeric(id))%>%
    filter(!is.na(id))%>%
    mutate(study=if_else(timepoint2 %in% c("pregnant", "post", "delivery"), "DPSP", "MICDROP"))%>%
    mutate(plate=unique(substr(file_name, 126, 132)))
  
  big_data <- rbind(big_data, long_sheet)
  
}

very_big_data <- rbind(big_data, pilot_nulisa_edit)

micdrop_nulisa <- very_big_data%>%
  filter(study=="MICDROP")%>%
  filter(sample %notin% c("10720_tp24_1", "10640_tp24_1"))%>%
  mutate(timepoint=case_when(timepoint2=="tp7"~"8 weeks",
                              timepoint2=="tp8"~"8 weeks",
                              timepoint2=="tp23"~"24 weeks",
                              timepoint2=="tp24"~"24 weeks",
                              timepoint2=="tp51"~"52 weeks",
                              timepoint2=="tp52"~"52 weeks",
                              timepoint2=="tp68"~"68 weeks"))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks", "68 weeks")))%>%
  mutate(timepoint_num=as.numeric(substr(timepoint, 1, 2)))


dpsp_nulisa <- very_big_data%>%
  filter(study=="DPSP")

write.csv(dpsp_nulisa, "~/postdoc/stanford/plasma_analytes/DPSP/dpsp_nulisa_data.csv", row.names = F)

# read in metadata####
## clinical data####
raw_data <- haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/Nov24/MICDSpecimenBoxNov24_withclinical.dta")

master_plate_map <- readxl::read_excel("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/master_plate_map.xlsx", col_types = c("numeric", "date", "numeric", "text", "guess", "guess", "guess"))
  
big_experiment_samples <- master_plate_map%>%
  select(id, date)%>%
  filter(!is.na(id))%>%
  mutate(date=as.Date(date))

pilot_samples <- raw_data%>%
  filter(id %in% c(unique(pilot_nulisa_edit$id)))%>%
  filter(Plasma!="")%>%
  select(id, date, ageinwks, Plasma)%>%
  mutate(sample=paste(id, "_tp", ageinwks, sep=""))%>%
  filter(sample %in% unique(pilot_nulisa_edit$sample))%>%
  add_row(id=10766, date=as.Date("2022-11-04"), ageinwks=25, Plasma="10766-AB", sample="10766_tp24")

all_samples <- bind_rows(pilot_samples, big_experiment_samples)%>%
  select(id, date)

metadata_columns <- c("id", "dob", "date", "ageinwks", "gender", "mstatus", "qPCRparsdens", "visittype", "fever", "febrile")

#merge based on id and date
metadata <- raw_data%>%
  select(all_of(metadata_columns))%>%
  right_join(., all_samples, by=c("id", "date"))%>%
  mutate(timepoint_num=as.numeric(ageinwks), id=as.numeric(id))%>%
  mutate(log_qpcr=log10(qPCRparsdens+0.01))
    
## infection data ####
infs_and_meta <- raw_data%>%
# filter(ageinwks<54)%>%
  group_by(id) %>%
  mutate("total_n_para_12"=sum(qPCRparsdens>1&ageinwks<53, na.rm=TRUE),
         "total_n_malaria_12"=sum(mstatus!=0&ageinwks<53, na.rm=TRUE),
         "total_n_para_6"=sum(qPCRparsdens>1&ageinwks<27, na.rm=TRUE),
         "total_n_malaria_6"=sum(mstatus!=0&ageinwks<27, na.rm=TRUE),
         "any_malar_6"=if_else(total_n_malaria_6==0, FALSE, TRUE),
         "any_malar_12"=if_else(total_n_malaria_12==0, FALSE, TRUE))%>%
  select(id, date, total_n_para_6, total_n_malaria_6, total_n_para_12, total_n_malaria_12, -gender)%>%
  right_join(., metadata, by=c("id", "date"))%>%
  group_by(id)%>%
  mutate(total_n_para=max(total_n_para_12, na.rm = T),
         total_n_malaria=max(total_n_malaria_12, na.rm = T))%>%
  mutate(timepoint_num=case_when(timepoint_num==9~8,
                                 timepoint_num==25~24,
                                 timepoint_num==53~52,
                                 .default=timepoint_num))


## apply QC ####
qc_results <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/NULISA_QC_summary.csv")

## put it together####
#duplciation happens here?
micdrop_nulisa_with_meta <- micdrop_nulisa%>%
  inner_join(., infs_and_meta, by=c("id", "timepoint_num"), multiple = "first")%>%
  mutate(qc_failed=if_else(sample %in% qc_results$sample_id, TRUE, FALSE))

#13% dropout; 642 - 86 samples; mostly plate 8 (37)

write.csv(micdrop_nulisa_with_meta, "~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv", row.names = F)


# bigass figure ####

micdrop_nulisa_with_meta%>%
  filter(mstatus==0)%>%
  filter(targetName%in%c("IL10", "IFNG", "TNF"))%>%
  ggplot(., aes(x=timepoint, y=conc, fill=factor(total_n_para)))+
  geom_boxplot()+
  facet_wrap(~targetName, nrow = 16, scales="free")+
  # geom_smooth(method="lm")+
  viridis::scale_fill_viridis(discrete = T)+
  theme_minimal()

micdrop_nulisa_with_meta%>%
  filter(qPCRparsdens==0)%>%
  filter(targetName%in%c("IL10", "IFNG", "TNF"))%>%
  ggplot(., aes(x=timepoint, y=conc, fill=factor(total_n_para)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName, nrow = 16, scales="free")+
  # geom_smooth(method="lm")+
  viridis::scale_fill_viridis(discrete = T)+
  theme_minimal()

micdrop_nulisa_with_meta%>%
  # filter(mstatus==0)%>%
  filter(targetName%in%c("IL10", "IFNG", "TNF"))%>%
  ggplot(., aes(x=qPCRparsdens+0.01, y=conc, color=total_n_para))+
  geom_point(aes(shape=factor(mstatus)))+
  facet_wrap(~targetName, nrow = 16, scales="free")+
  # geom_smooth(method="lm")+
  scale_x_log10()+
  viridis::scale_color_viridis(discrete = F)+
  theme_minimal()

big_ass_figure2 <- micdrop_nulisa_with_meta%>%
  filter(qPCRparsdens==0, timepoint=="52 weeks")%>%
  # filter(targetName%in%c("IL10", "IFNG", "TNF", "GZMA", "LAG3", "LILRB2"))%>%
  ggplot(., aes(x=total_n_para, y=conc))+
  facet_wrap(~targetName+timepoint, scales="free")+
  geom_point()+
  ggpubr::stat_cor(method="spearman", size=2, color="red")+
  geom_smooth(method="lm")+
  viridis::scale_fill_viridis(discrete = T)+
  theme_minimal()
# 
ggsave("~/Downloads/big_ass_figure2.png", big_ass_figure2, width=24, height=24, dpi=444, bg="white")

big_ass_figure3 <- micdrop_nulisa_with_meta%>%
  filter(qPCRparsdens==0, timepoint=="52 weeks")%>%
  # filter(targetName%in%c("IL10", "IFNG", "TNF", "GZMA", "LAG3", "LILRB2"))%>%
  ggplot(., aes(x=total_n_malaria, y=conc))+
  facet_wrap(~targetName+timepoint, scales="free")+
  geom_point()+
  ggpubr::stat_cor(method="spearman", size=2, color="red")+
  geom_smooth(method="lm")+
  viridis::scale_fill_viridis(discrete = T)+
  theme_minimal()
# 
ggsave("~/Downloads/big_ass_figure3.png", big_ass_figure3, width=24, height=24, dpi=444, bg="white")

micdrop_nulisa_with_meta%>%
  # filter(qPCRparsdens==0, timepoint=="52 weeks")%>%
  filter(mstatus==0, timepoint=="52 weeks")%>%
  filter(targetName%in%c("CXCL1"))%>%
  ggplot(., aes(x=total_n_malar_12, y=conc))+
  facet_wrap(~targetName+timepoint, scales="free")+
  geom_point()+
  ggpubr::stat_cor(method="spearman", size=2, color="red")+
  geom_smooth(method="lm")+
  viridis::scale_fill_viridis(discrete = T)+
  theme_minimal()



clean_data %>%
  filter(targetName %in% c("IL10", "LAG3", "IL2", "GZMA"),
         timepoint!="68 weeks")%>%
  ggplot(aes(x=factor(timepoint_num), y=conc, fill=timepoint))+
  geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetName, scales = "free")+
  scale_fill_manual(values=viridis::magma(3))+
  xlab("age in weeks")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_blank())

# 
# micdrop_nulisa%>%
#   select(-targetName, -conc)%>%
#   distinct() -> micdrop_nulisa_metadata
# 
# write.csv(micdrop_nulisa_metadata, "~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/NULISA/micdrop_nulisa_metadata.csv")
# 
