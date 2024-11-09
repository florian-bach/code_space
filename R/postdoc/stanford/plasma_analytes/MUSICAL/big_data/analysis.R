library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)


`%notin%` <- Negate(`%in%`)


# data generation ####
  ## reading in pilot data ####


musical_metadata <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/MASTER_METADATA.csv")
random_codes <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/pilot/id_date_code.csv")
random_codes$plasma.barcode <- gsub("D1PN8A", "D1_PN8A", random_codes$plasma.barcode)
random_codes$plasma.barcode <- gsub("D1JLGS", "DIJLGS", random_codes$plasma.barcode)
random_codes$plasma.barcode <- gsub("D1KWT2", "D1_KWT2", random_codes$plasma.barcode)
random_codes$plasma.barcode <- gsub("D1EF4A", "DEF4A", random_codes$plasma.barcode)
random_codes$date <- as.Date(random_codes$date)

slim_musical_metadata <- musical_metadata %>%
  mutate(day_annotation=if_else(day_annotation==84, -1, day_annotation))%>%
  select(combined_id, combined_date, enrolltype, day_annotation, gender_categorical, ageyrs, qpcr, TEMP)%>%
  mutate(id=combined_id, date=combined_date, class=enrolltype, timepoint=paste("t", day_annotation, sep=""), temperature=TEMP)%>%
  select(-combined_id, -combined_date, -enrolltype, -day_annotation)%>%
  mutate("study"="MUSICAL", date=as.Date(date))


combo_frame2 <- merge(slim_musical_metadata, random_codes, by=c("id", "date"))

nulisa <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/pilot/nulisa_data.csv")

wide_nulisa <- nulisa %>%
  pivot_longer(cols = colnames(nulisa)[2:ncol(nulisa)], names_to = "plasma.barcode", values_to = "concentration")

wide_nulisa <- inner_join(wide_nulisa, combo_frame2, by="plasma.barcode")

pilot_nulisa <- wide_nulisa %>%
  mutate("time_class"=paste(class, timepoint, sep='_'),
         "age_class"=if_else(.$id %in% c(268, 324, 137, 176, 353, 161, 363, 571, 10766, 10794, 10842), "child", "adult"),
         id=factor(id))%>%
  filter(age_class=="child")%>%
  mutate(timepoint=case_when(timepoint=="t-1"~"baseline",
                             timepoint=="t0"~"day0",
                             timepoint=="t14"~"day14",
                             timepoint=="t7"~"day7"),
         timepoint=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")),
         barcode=plasma.barcode,
         sample_id=paste("X", id, time_class, barcode))%>%
  group_by(targetName) %>%
  mutate("z_conc"=scale(concentration, center = TRUE, scale = TRUE))%>%
  ungroup()%>%
  group_by(id)%>%
  mutate("mean_z_conc"=median(z_conc), study="pilot")%>%
  select(targetName, sample_id, concentration, barcode, id, gender_categorical, class, timepoint, z_conc, mean_z_conc, ageyrs, study, qpcr, temperature)


  ## reading in big batch data ####
metadata <- readxl::read_excel("~/postdoc/stanford/plasma_analytes/MUSICAL/big_data/Immunology List _MUSICAL.xlsx")

list_of_files <- list.files(path = "~/postdoc/stanford/plasma_analytes/MUSICAL/big_data/", pattern = "Report.csv", full.names = TRUE)

list_of_batches <- list(vector(length = length(list_of_files)))
# list_of_batches <- list()

system.time(for(i in 1:length(list_of_files)){
  
  tmp <- read.csv(list_of_files[i], header = TRUE)
  list_of_batches[[i]] <- tmp
  
})

list_of_long_batches <- lapply(list_of_batches, function(x) x %>%
                                 pivot_longer(cols= starts_with("X"),
                                              names_to = "sample_id",
                                              values_to = "concentration")%>%
                                 select(targetName, sample_id, concentration))

big_df <- do.call(rbind, list_of_long_batches)

big_df1 <- big_df %>%
  mutate(barcode= substr(sample_id, nchar(sample_id)-5, nchar(sample_id)),
         id=factor(metadata$id[match(barcode, metadata$plasma1_barcode)]),
         gender_categorical=factor(metadata$gender_categorical[match(barcode, metadata$plasma1_barcode)]),
         class=metadata$infectiontype[match(barcode, metadata$plasma1_barcode)],
         ageyrs=metadata$ageyrs[match(barcode, metadata$plasma1_barcode)],
         timepoint_imm=metadata$timepoint_imm[match(barcode, metadata$plasma1_barcode)],
         qpcr=metadata$qpcr[match(barcode, metadata$plasma1_barcode)],
         temperature=metadata$temperature[match(barcode, metadata$plasma1_barcode)]
  )%>%
  filter(barcode != "Repeat")%>%
  mutate("timepoint"=case_when(timepoint_imm==-2~"bad_baseline",
                               timepoint_imm==-1~"baseline",
                               timepoint_imm==0~"day0",
                               timepoint_imm==7~"day7",
                               timepoint_imm==14~"day14",
                               timepoint_imm==28~"day28")
         
  )%>%
  mutate(timepoint=factor(timepoint, levels=c("bad_baseline", "baseline", "day0", "day7", "day14", "day28")))%>%
  group_by(targetName) %>%
  mutate(z_conc=scale(concentration, center = TRUE, scale = TRUE))%>%
  ungroup()%>%
  group_by(id)%>%
  mutate("mean_z_conc"=median(z_conc), study="big_batch")%>%
  ungroup()%>%
  select(-timepoint_imm)



  ## putting it all together ####
# 
big_df2 <- rbind(big_df1, pilot_nulisa)
# big_df2$timepoint <- 
# samples_to_ditch <- big_df2%>%
#   filter(mean_z_conc <= -0.4)%>%
#   select(sample_id)%>%
#   filter(!duplicated(sample_id))

write.csv(big_df2, "~/postdoc/stanford/plasma_analytes/MUSICAL/unclean_musical_combo_with_metadata.csv", row.names = FALSE)
# read clean data ####

big_df <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/unclean_musical_combo_with_metadata.csv")
big_df2 <- big_df %>%
  mutate(id=paste(study, id, sep="_"),
         timepoint=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")), 
         qpcr=as.numeric(qpcr),
         qpcr_cat=case_when(qpcr==0~"0", 
                            qpcr>0 & qpcr<10~ ">1",
                            qpcr>10 & qpcr<100~ ">10",
                            qpcr>100 & qpcr<1000~ ">10e2",
                            qpcr>1000 & qpcr<10000~ ">10e3",
                            qpcr>10000 & qpcr<100000~ ">10e4",
                            qpcr>100000 & qpcr<1000000~ ">10e5",
                            qpcr>1000000 & qpcr<10000000~ ">10e6",
                            qpcr>10000000 ~ ">10e7"),
         qpcr_cat=factor(qpcr_cat, levels=c("0", ">1", ">10", ">10e2", ">10e3", ">10e4", ">10e5", ">10e6", ">10e7")),
         temperature = as.numeric(temperature),
         temperature_cat = case_when(temperature < 38 ~ "<38",
                                     temperature >38 & temperature <39 ~ ">38",
                                     temperature >39 & temperature <40 ~ ">39",
                                     temperature >40 & temperature <41 ~ ">40",
                                     temperature >41 & temperature <42 ~ ">41"))%>%
  filter(mean_z_conc > -0.4)%>%
  filter(sample_id %notin% c("X384_S.1_D1V93U",
                             "X744_A.1_D1KT5Z",
                             "X323_A14_D1DZ2D",
                             "X496_NM7_D1CAYS",
                             "X316_S.1_D12FNR",
                             "X667_NM7_D1GBYA",
                             "X 176 S_t14 D1FXRJ",
                             "X164_NM0_D1WA6Q",
                             "X 176 S_t7 D1VFF7",
                             "X219_S.1_D1C4ZM",
                             "X132_S0_D1AUZD"
                             ))%>%
  filter(barcode %notin% c("D19E2G", "D1FK67", "D1SSPJ"))
 
  
# big heatmap ####
  ## unclean ####
t_nulisa <-t(as.matrix(big_df %>%
                         # filter(mean_z_conc > -0.4)%>%
                         select(sample_id, targetName, z_conc)%>%
                         pivot_wider(names_from = targetName, values_from = z_conc)))


colnames(t_nulisa) <- t_nulisa[1,]

num_nulisa_matrix <- t_nulisa[-1,]
class(num_nulisa_matrix) <- "numeric"

short_num_nulisa_matrix <- num_nulisa_matrix

# z_col_fun <- circlize::colorRamp2(c(min(short_num_nulisa_matrix)/1.5, 0, abs(min(short_num_nulisa_matrix)/1.5)), c("#123499", "white", "#FFA500"))
z_col_fun <- circlize::colorRamp2(c(min(-10.90915)/1.5, 0, abs(min(-10.90915)/1.5)), c("#123499", "white", "#FFA500"))

# big_heatmap_split <- substr(colnames(t_nulisa), 5, nchar(colnames(t_nulisa)))
# 
# big_heatmap_split <- factor(big_heatmap_split, levels=sort(unique(big_heatmap_split))[c(1:5, 7, 6)])


big_heatmap <- Heatmap(matrix = short_num_nulisa_matrix,
                       cluster_rows = TRUE,
                       cluster_columns=TRUE,
                       show_row_dend = FALSE,
                       show_column_dend = FALSE,
                       show_heatmap_legend = TRUE,
                       # column_split = big_heatmap_split,
                       name = "z score",
                       column_names_gp = gpar(fontsize = 3),
                       row_names_gp = gpar(fontsize = 0),
                       row_names_side = "left",
                       col = z_col_fun,
                       column_names_rot = 90)



png("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/unclean_z_heat.png", width=16, height = 9, res = 444, units = "in")
draw(big_heatmap)
dev.off()

  ## clean data ####

t_nulisa <-t(as.matrix(big_df2 %>%
                         # filter(mean_z_conc > -0.4)%>%
                         select(sample_id, targetName, z_conc)%>%
                         pivot_wider(names_from = targetName, values_from = z_conc)))


colnames(t_nulisa) <- t_nulisa[1,]

num_nulisa_matrix <- t_nulisa[-1,]
class(num_nulisa_matrix) <- "numeric"

short_num_nulisa_matrix <- num_nulisa_matrix

# z_col_fun <- circlize::colorRamp2(c(min(short_num_nulisa_matrix)/1.5, 0, abs(min(short_num_nulisa_matrix)/1.5)), c("#123499", "white", "#FFA500"))
z_col_fun <- circlize::colorRamp2(c(min(-10.90915)/1.5, 0, abs(min(-10.90915)/1.5)), c("#123499", "white", "#FFA500"))

# big_heatmap_split <- substr(colnames(t_nulisa), 5, nchar(colnames(t_nulisa)))
# 
# big_heatmap_split <- factor(big_heatmap_split, levels=sort(unique(big_heatmap_split))[c(1:5, 7, 6)])


big_heatmap <- Heatmap(matrix = short_num_nulisa_matrix,
                       cluster_rows = TRUE,
                       cluster_columns=TRUE,
                       show_row_dend = FALSE,
                       show_column_dend = FALSE,
                       show_heatmap_legend = TRUE,
                       # column_split = big_heatmap_split,
                       name = "z score",
                       column_names_gp = gpar(fontsize = 3),
                       row_names_gp = gpar(fontsize = 0),
                       row_names_side = "left",
                       col = z_col_fun,
                       column_names_rot = 90)



png("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/clean_z_heat.png", width=16, height = 9, res = 444, units = "in")
draw(big_heatmap)
dev.off()

# linear regression ####
## symptomatic only ####
fdr_cutoff = 0.05

base_zero_contrast <- t(matrix(c(0,1,0,0, rep(0, 0))))
base_7_contrast <- t(matrix(c(0,0,1,0, rep(0, 0))))
base_14_contrast <- t(matrix(c(0,0,0,1, rep(0, 0))))
zero_14_contrast <- t(matrix(c(0,-1,0,1, rep(0, 0))))



symp_only_purff <- big_df2 %>%
  filter(class=="S", timepoint!="day28", timepoint!="bad_baseline")%>%
  group_by(targetName)%>%
  nest() %>%
  # mutate(model=map(data, ~lm(concentration~timepoint+ageyrs+id, data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~timepoint+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_zero=map(model, ~multcomp::glht(., base_zero_contrast)),
         base_zero_p=map_dbl(base_zero, ~summary(.)$test$pvalues)) %>%
  mutate(base_14=map(model, ~multcomp::glht(., base_14_contrast)),
         base_14_p=map_dbl(base_14, ~summary(.)$test$pvalues)) %>%
  mutate(base_7=map(model, ~multcomp::glht(., base_7_contrast)),
         base_7_p=map_dbl(base_7, ~summary(.)$test$pvalues))%>%
  ungroup()%>%
  mutate(base_zero_padj=p.adjust(base_zero_p, method="BH"),
         base_7_padj=p.adjust(base_7_p, method="BH"),
         base_14_padj=p.adjust(base_14_p, method="BH")
  )


symp_results_table <- symp_only_purff %>%
  dplyr::select(targetName,
                base_zero_padj,
                base_7_padj,
                base_14_padj)%>%
  ungroup()

#79 with lm(); 159 with lme4;
symp_sig_base_zero <- symp_results_table %>%
  filter(base_zero_padj<fdr_cutoff)%>%
  arrange(base_zero_padj)%>%
  select(targetName)

#0 with lm(); 36 with lme4
symp_sig_base_7 <- symp_results_table %>%
  filter(base_7_padj<fdr_cutoff)%>%
  arrange(base_7_padj)%>%
  select(targetName)

#2 with lme4
symp_sig_base_14 <- symp_results_table %>%
  filter(base_14_padj<fdr_cutoff)%>%
  select(targetName)

da_boxplot_theme <- theme(legend.position = "none",
                          axis.title = element_blank())


big_df2 %>%
  filter(class=="S", timepoint%notin%c("day28"), timepoint!="bad_baseline")%>%
  filter(targetName %in% c("IL21", "IFNG"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
  # geom_line(aes(group=id))+
  # geom_point()+#
  geom_boxplot(outliers = FALSE)+
  # geom_violin(draw_quantiles = c(0,0.25,0.5,0.75), color="white")+
  ggtitle("symptomatic malaria")+
  facet_wrap(~targetName, scales = "free")+
  scale_fill_manual(values=viridis::magma(5))+
  theme_minimal()


for(i in 1:ceiling(nrow(symp_sig_base_zero)/16)){
  
plt <- big_df2 %>%
  filter(class=="S", timepoint!="day28", timepoint!="bad_baseline")%>%
  filter(targetName %in% symp_sig_base_zero$targetName[seq((i-1)*16+1, i*16)])%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
  # geom_point()+#
  # geom_line(aes(group=id))+
  geom_boxplot(outliers = FALSE)+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  ggtitle("upregulated during malaria")+
  facet_wrap(~targetName, scales = "free")+
  scale_fill_manual(values=viridis::magma(5))+
  theme_minimal()+
  da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/symp_base_zero", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
  }


for(i in 1:ceiling(nrow(symp_sig_base_7)/16)){
  
  plt <- big_df2 %>%
    filter(class=="S", timepoint!="day28", timepoint!="bad_baseline")%>%
    filter(targetName %in% symp_sig_base_7$targetName[seq((i-1)*16+1, i*16)])%>%
    mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
    ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
    # geom_point()+#
    # geom_line(aes(group=id))+
    geom_boxplot(outliers = FALSE)+
    # geom_violin(draw_quantiles = seq(0,1,0.25))+
    ggtitle("upregulated on day 7 post malaria")+
    facet_wrap(~targetName, scales = "free")+
    scale_fill_manual(values=viridis::magma(5))+
    theme_minimal()+
    da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/symp_base_7", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}

analytes_da_7_not_0 <- symp_sig_base_7$targetName[(symp_sig_base_7$targetName %notin% symp_sig_base_zero$targetName)]

analytes_da_7_not_0_plot <- big_df2 %>%
  filter(class=="S", timepoint!="day28", timepoint!="bad_baseline")%>%
  filter(targetName %in% analytes_da_7_not_0)%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
  # geom_point()+#
  # geom_line(aes(group=id))+
  geom_boxplot(outliers = FALSE)+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  ggtitle("upregulated on day 7 post malaria")+
  facet_wrap(~targetName, scales = "free", ncol=4)+
  scale_fill_manual(values=viridis::magma(5))+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/symp_base_7_not_0.png", analytes_da_7_not_0_plot, height=4, width=8, dpi=444, bg="white")


## asymptomatic only ####
# base_zero_contrast <- t(matrix(c(0,1,0,0, rep(0,44))))
# base_14_contrast <- t(matrix(c(0,0,1,0, rep(0,44))))
# # base_28_contrast <- t(matrix(c(0,0,0)))
# zero_14_contrast <- t(matrix(c(0,-1,1,0, rep(0,44))))

base_zero_contrast <- t(matrix(c(0,1,0, rep(0,0))))
base_14_contrast <- t(matrix(c(0,0,1, rep(0,0))))
# base_28_contrast <- t(matrix(c(0,0,0)))
zero_14_contrast <- t(matrix(c(0,-1,1, rep(0,0))))

asymp_only_purff <- big_df2 %>%
  filter(class=="A", timepoint!="day28", timepoint!="bad_baseline")%>%
  group_by(targetName)%>%
  nest() %>%
  mutate(model=map(data, ~lme4::lmer(concentration~timepoint+(1|id), data=.))) %>%
  # mutate(model=map(data, ~lm(concentration~timepoint+ageyrs+id, data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_zero=map(model, ~multcomp::glht(., base_zero_contrast)),
         base_zero_p=map_dbl(base_zero, ~summary(.)$test$pvalues)) %>%
  mutate(base_14=map(model, ~multcomp::glht(., base_14_contrast)),
         base_14_p=map_dbl(base_14, ~summary(.)$test$pvalues)) %>%
  mutate(zero_14=map(model, ~multcomp::glht(., zero_14_contrast)),
         zero_14_p=map_dbl(zero_14, ~summary(.)$test$pvalues))%>%
  # mutate(base_28=map(model, ~multcomp::glht(., base_28_contrast)),
  #        base_28_p=map_dbl(base_28, ~summary(.)$test$pvalues))%>%
  ungroup()%>%
  mutate(base_zero_padj=p.adjust(base_zero_p, method="BH"),
         base_14_padj=p.adjust(base_14_p, method="BH"),
         zero_14_padj=p.adjust(zero_14_p, method="BH"),
         # base_28_padj=p.adjust(base_28_p, method="BH")
  )


asymp_results_table <- asymp_only_purff %>%
  dplyr::select(targetName,
                base_zero_padj,
                base_14_padj,
                zero_14_padj
  )%>%
  ungroup()

#0 with lm; 43 with lme4
asymp_sig_base_zero <- asymp_results_table %>%
  filter(base_zero_padj<fdr_cutoff)%>%
  arrange(base_zero_padj)%>%
  select(targetName)

#8 with lm; 52 with lme4
asymp_sig_base_14 <- asymp_results_table %>%
  filter(base_14_padj<fdr_cutoff)%>%
  arrange(base_14_padj)%>%
  select(targetName)

#0 with either
asymp_sig_zero_14 <- asymp_results_table %>%
  filter(zero_14_padj<fdr_cutoff)%>%
  arrange(zero_14_padj)%>%
  select(targetName)


for(i in 1:ceiling(nrow(asymp_sig_base_zero)/16)){
  
  plt <- big_df2 %>%
    filter(class=="A", timepoint!="day28", timepoint!="bad_baseline")%>%
    filter(targetName %in% asymp_sig_base_zero$targetName[seq((i-1)*16+1, i*16)])%>%
    mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
    ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
    # geom_point()+#
    # geom_line(aes(group=id))+
    geom_boxplot(outliers = FALSE)+
    # geom_violin(draw_quantiles = seq(0,1,0.25))+
    ggtitle("regulated during asymptomatic parasitemia")+
    facet_wrap(~targetName, scales = "free")+
    scale_fill_manual(values=viridis::magma(5))+
    theme_minimal()+
    da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/asymp_base_zero", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}


for(i in 1:ceiling(nrow(asymp_sig_base_14)/16)){
  
  plt <- big_df2 %>%
    filter(class=="A", timepoint!="day28", timepoint!="bad_baseline")%>%
    filter(targetName %in% asymp_sig_base_14$targetName[seq((i-1)*16+1, i*16)])%>%
    mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
    ggplot(aes(x=factor(timepoint), y=concentration, fill=timepoint))+
    # geom_point()+#
    # geom_line(aes(group=id))+
    geom_boxplot(outliers = FALSE)+
    # geom_violin(draw_quantiles = seq(0,1,0.25))+
    ggtitle("upregulated on day 14 of asymptomatic parasitemia")+
    facet_wrap(~targetName, scales = "free")+
    scale_fill_manual(values=viridis::magma(5))+
    theme_minimal()+
    da_boxplot_theme
  
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/asymp_base_14", (i-1)*16+1, i*16, ".png", sep="_"), plt, height=8, width=8, dpi=444, bg="white")
  
}

# 
base_zero_asymp_plot <- big_df2 %>%
  filter(class=="A", timepoint!="day28", timepoint!="bad_baseline")%>%
  filter(targetName %in% c("IL10", "LAG3", "GZMA", "PDCD1", "IL2RA", "CTLA4"))%>%
  ggplot(., aes(x=timepoint, y=concentration, fill=timepoint))+
  geom_boxplot()+
  ggtitle("regulated during asymptomatic parasitemia")+
  facet_wrap(~targetName, scales = "free", ncol = 3)+
  scale_fill_manual(values=c(viridis::magma(4)))+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/select_asymp_base_zero.png", base_zero_asymp_plot, width=8, height=6, bg="white", dpi=444)

# base_14_asymp_plot <- big_df2 %>%
#   filter(class=="A", timepoint!="day28", timepoint!="bad_baseline")%>%
#   filter(targetName %in% asymp_sig_base_14$targetName)%>%
#   ggplot(., aes(x=timepoint, y=concentration, fill=timepoint))+
#   geom_boxplot()+
#   ggtitle("upregulated at day 14 of asymptomatic parasitemia")+
#   facet_wrap(~targetName, scales = "free", ncol = 4)+
#   scale_fill_manual(values=c(viridis::magma(4)))+
#   theme_minimal()+
#   da_boxplot_theme
# 
# ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/asymp_base_14.png", base_14_asymp_plot, width=8, height=4, bg="white", dpi=444)

## baselines only ####
base_base_contrast <- t(matrix(c(0,1)))


base_only_purff <- big_df2 %>%
  filter(class %in% c("S", "A"), timepoint=="baseline")%>%
  group_by(targetName)%>%
  mutate(bino_class=if_else(class=="S", 1, 0))%>%
  nest() %>%
  # mutate(model=map(data, ~lm(concentration~class+id, data=.))) %>%
  # mutate(model=map(data, ~lme4::lmer(concentration~class+(1|id), data=.))) %>%
  # this one "works"
  mutate(model=map(data, ~glm(bino_class~concentration+id, data=., family="binomial"))) %>%
  # mutate(model=map(data, ~lme4::lmer(concentration~class+(1|id), data=.)))%>%
  # mutate(model=map(data, ~lme4::lmer(concentration~class+ageyrs+gender_categorical+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  # mutate(base_base=map(model, ~multcomp::glht(., base_base_contrast)),
  #        base_base_p=map_dbl(base_base, ~summary(.)$test$pvalues)) %>%
  # this one "works" too
  mutate(base_base_p=map_dbl(summary, ~coef(.)[140])) %>%
  # mutate(base_14=map(model, ~multcomp::glht(., base_14_contrast)),
  #        base_14_p=map_dbl(base_14, ~summary(.)$test$pvalues)) %>%
  # mutate(zero_14=map(model, ~multcomp::glht(., zero_14_contrast)),
  #        zero_14_p=map_dbl(zero_14, ~summary(.)$test$pvalues))%>%
  # mutate(zero_28=map(model, ~multcomp::glht(., zero_28_contrast)),
  #        zero_28_p=map_dbl(zero_28, ~summary(.)$test$pvalues))%>%
  ungroup()%>%
  mutate(base_base_padj=p.adjust(base_base_p, method="BH"),
         # base_14_padj=p.adjust(base_14_p, method="BH"),
         # zero_14_padj=p.adjust(zero_14_p, method="BH"),
         # zero_28_padj=p.adjust(zero_28_p, method="BH")
  )

base_only_results_table <- base_only_purff %>%
  dplyr::select(targetName,
                base_base_p,
                base_base_padj)%>%
  ungroup()

#18 with glm
base_only_sig <- base_only_results_table %>%
  filter(base_base_padj<fdr_cutoff)%>%
  select(targetName)


base_only_plot <- big_df2 %>%
  filter(timepoint=="baseline", class %in% c("S", "A"))%>%
  # filter(targetName %in% c("CXCL8", "MMP8", "MMP9", "MP9", "OSM", "PTX3", "TNFSF11"))%>%
  filter(targetName %in% base_only_sig$targetName)%>%
  ggplot(aes(x=class, y=concentration, fill=targetName))+
  # geom_violin(draw_quantiles = c(seq(0,1,0.25)), colour = "darkgrey")+
  geom_boxplot(outliers = FALSE)+
  ggtitle("baselines before parasitemia")+
  facet_wrap(~targetName, scales = "free", ncol=6)+
  scale_fill_manual(values=c(viridis::magma(18)))+
  theme_minimal()+
  theme(legend.position = "none")+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/base_only_bino_glm.png", base_only_plot, width=8, height=6, bg="white", dpi=444)


big_df2 %>%
  filter(timepoint=="baseline", class %in% c("S", "A"))%>%
  # filter(targetName %in% c("CXCL8", "MMP8", "MMP9", "MP9", "OSM", "PTX3", "TNFSF11"))%>%
  filter(targetName %in% base_only_sig$targetName)%>%
  ggplot(aes(x=class, y=concentration, fill=gender_categorical))+
  geom_boxplot()+
  # geom_violin(draw_quantiles = c(seq(0,1,0.25)), colour = "darkgrey")+
  # geom_point()+
  # geom_line(aes(group=id))+
  ggtitle("baselines before parasitemia")+
  facet_wrap(~targetName, scales = "free")+
  # scale_fill_manual(values=c(viridis::magma(10)))+
  theme_minimal()+
  theme()


## big plots ####
big_df2 %>%
  # filter(targetName %in% c("MPO", "MMP8", "MMP9", "OSM"), class %in% c("A", "S"), timepoint !="day7", timepoint!="day28")%>%
  filter(targetName %in% base_only_sig$targetName, class %in% c("A", "S"), timepoint !="day7", timepoint!="day28")%>%
  ggplot(aes(x=timepoint, y=concentration, fill=interaction(class, timepoint)))+
  # geom_line(aes(group=id), position = position_dodge(width = .75))+
  geom_boxplot()+
  geom_point(position = position_dodge(width=0.75))+# geom_boxplot(outliers = FALSE)+
  ggtitle("symptomatic malaria")+
  facet_wrap(~targetName, scales = "free")+
  # scale_fill_manual(values=c(viridis::magma(5)))+
  theme_minimal()

big_df2 %>%
  filter(targetName %in% c("CXCL10", "IL10", "IFNB", "IFNG"), class %in% c("A", "S"), timepoint%in%c("baseline", "day0"))%>%
  ggplot(aes(x=timepoint, y=concentration, fill=interaction(class, timepoint)))+
  geom_boxplot()+
  geom_point(position = position_dodge(width=0.75))+# geom_boxplot(outliers = FALSE)+
  ggtitle("symptomatic malaria")+
  facet_wrap(~targetName, scales = "free")+
  # scale_fill_manual(values=c(viridis::magma(5)))+
  theme_minimal()

## nmf ####
base_zero_contrast <- t(matrix(c(0,1,0,0)))
base_7_contrast <- t(matrix(c(0,0,1,0)))
base_14_contrast <- t(matrix(c(0,0,0,1)))
zero_14_contrast <- t(matrix(c(0,-1,0,1)))



nmf_only_purff <- big_df2 %>%
  filter(class=="NM", timepoint!="day28")%>%
  group_by(targetName)%>%
  nest() %>%
  mutate(model=map(data, ~lme4::lmer(concentration~timepoint+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_zero=map(model, ~multcomp::glht(., base_zero_contrast)),
         base_zero_p=map_dbl(base_zero, ~summary(.)$test$pvalues)) %>%
  mutate(base_14=map(model, ~multcomp::glht(., base_14_contrast)),
         base_14_p=map_dbl(base_14, ~summary(.)$test$pvalues)) %>%
  mutate(base_7=map(model, ~multcomp::glht(., base_7_contrast)),
         base_7_p=map_dbl(base_7, ~summary(.)$test$pvalues))%>%
  mutate(zero_14=map(model, ~multcomp::glht(., zero_14_contrast)),
         zero_14_p=map_dbl(zero_14, ~summary(.)$test$pvalues))%>%
  ungroup()%>%
  mutate(base_zero_padj=p.adjust(base_zero_p, method="BH"),
         base_7_padj=p.adjust(base_7_p, method="BH"),
         base_14_padj=p.adjust(base_14_p, method="BH"),
         zero_14_padj=p.adjust(zero_14_p, method="BH")
  )


nmf_results_table <- nmf_only_purff %>%
  dplyr::select(targetName,
                base_zero_padj,
                base_7_padj,
                base_14_padj,
                zero_14_padj)%>%
  ungroup()

nmf_sig_base_zero <- nmf_results_table %>%
  filter(base_zero_padj<fdr_cutoff)%>%
  arrange(desc(base_zero_padj))%>%
  select(targetName)

nmf_sig_base_7 <- nmf_results_table %>%
  filter(base_7_padj<fdr_cutoff)%>%
  select(targetName)

nmf_sig_base_14 <- nmf_results_table %>%
  filter(base_14_padj<fdr_cutoff)%>%
  select(targetName)

nmf_sig_zero_14 <- nmf_results_table %>%
  filter(zero_14_padj<fdr_cutoff)%>%
  select(targetName)




nmf_plot <- big_df2 %>%
  filter(targetName %in% c("CX3CL1", "CCL2", "CXCL10", "IL10"), class %in% c("NM"), timepoint!="day28")%>%
  ggplot(aes(x=timepoint, y=concentration, fill=interaction(class, timepoint)))+
  geom_boxplot()+
  geom_point(position = position_dodge(width=0.75))+# geom_boxplot(outliers = FALSE)+
  ggtitle("non-malarial fever")+
  facet_wrap(~targetName, scales = "free")+
  scale_fill_manual(values=c(viridis::magma(4)))+
  theme_minimal()+
  da_boxplot_theme

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/nmf_plot.png", nmf_plot, width=4, height=4, bg="white", dpi=444)

# PCA ####

id_columns <- c("sample_id", "id", "gender_categorical", "ageyrs", "timepoint", "class", "study", "qpcr", "qpcr_cat")
wide_df2 <- big_df2 %>%
  filter(targetName !="CTSS")%>%
  filter(sample_id %notin% c("X384_S.1_D1V93U", "X744_A.1_D1KT5Z", "X323_A14_D1DZ2D", "X496_NM7_D1CAYS", "X316_S.1_D12FNR", "X667_NM7_D1GBYA", "X 176 S_t14 D1FXRJ", "X164_NM0_D1WA6Q"))%>%
  pivot_wider(names_from = targetName, values_from = concentration, id_cols = all_of(id_columns))


rownames(wide_df2) <- wide_df2$sample_id
# each row = measurement; each column = feature
big_pca <-  prcomp(wide_df2[,(length(id_columns)+1):ncol(wide_df2)], center = T)
pca_plot_data <- as.data.frame(cbind(wide_df2, big_pca$x))

loadings_df <- data.frame(big_pca$rotation)
loadings_df$targetName <- rownames(loadings_df)
loadings_df$targetName <- factor(loadings_df$targetName, levels = loadings_df$targetName[order(loadings_df$PC1)])



pc1_cols <- colorspace::sequential_hcl(nrow(loadings_df), palette = "Purple Yellow")
names(pc1_cols) <- loadings_df$targetName[order(loadings_df$PC1)]

PC1_plot <- ggplot(loadings_df, aes(x=factor(targetName, levels = targetName[order(loadings_df$PC1)]), y=PC1, fill=targetName))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = pc1_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")

PC2_plot <- ggplot(loadings_df, aes(x=factor(targetName, levels = targetName[order(PC2)]), y=PC2, fill=targetName))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = pc1_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")

PC3_plot <- ggplot(loadings_df, aes(x=factor(targetName, levels = targetName[order(PC3)]), y=PC3, fill=targetName))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = pc1_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")


age_pca_plot <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  mutate("over_six"=if_else(ageyrs>6, "over 6y", "under 6y"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(., aes(x=PC1, y=PC2, color=over_six))+
  geom_point()+
  # stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=c("violet", "chartreuse"))+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/age_pca_plot.png", age_pca_plot, width=4, height=4, bg="white", dpi=444)


sex_pca_plot <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  mutate("over_six"=if_else(ageyrs>6, "over 6y", "under 6y"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(., aes(x=PC1, y=PC2, color=gender_categorical))+
  geom_point()+
  # stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=c("darkblue", "orange"))+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/sex_pca_plot.png", sex_pca_plot, width=4, height=4, bg="white", dpi=444)


time_pca_plot <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  mutate("over_six"=if_else(ageyrs>6, "over 6y", "under 6y"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(., aes(x=PC1, y=PC2, color=timepoint))+
  geom_point()+
  # stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=viridis::magma(5))+
  # facet_wrap(~timepoint)+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/time_pca_plot.png", time_pca_plot, width=4, height=4, bg="white", dpi=444)


study_pca_plot <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  mutate("over_six"=if_else(ageyrs>6, "over 6y", "under 6y"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(., aes(x=PC1, y=PC2, color=study))+
  geom_point()+
  # stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=c("black", "red"))+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/study_pca_plot.png", study_pca_plot, width=4, height=4, bg="white", dpi=444)


class_pca_plot <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  mutate("over_six"=if_else(ageyrs>6, "over 6y", "under 6y"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(., aes(x=PC1, y=PC2, color=class))+
  geom_point()+
  # stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=viridis::magma(6))+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/class_pca_plot.png", class_pca_plot, width=4, height=4, bg="white", dpi=444)



qpcr_pca_plot <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  mutate("over_six"=if_else(ageyrs>6, "over 6y", "under 6y"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  ggplot(., aes(x=PC1, y=PC2, color=qpcr_cat))+
  geom_point()+
  # stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=colorspace::sequential_hcl(n = 7, "LaJolla", rev = TRUE))+
  theme(legend.title = element_blank(),
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/age_pca_plot.png", age_pca_plot, width=4, height=4, bg="white", dpi=444)

## PC || covariate examination####

loop_data <- pca_plot_data %>%
  filter(timepoint %notin% c("day28", "bad_baseline"))%>%
  mutate("over_six"=if_else(ageyrs>6, "over 6y", "under 6y"))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))

plot_list <- list()

for(color_by in c("gender_categorical", "class", "over_six", "study")){
for(i in seq(1,9, by=2)){
  
  xvar <- paste("PC", i, sep="")
  yvar <- paste("PC", i+1, sep="")
  
  # print(head(loop_data[,xvar]))
  # print(head(loop_data[,yvar]))
  
  plt <- ggplot(loop_data, aes_string(x=xvar, y=yvar, color=color_by))+
    geom_point()+
    xlab(paste0("PC", i, " ", data.frame(summary(big_pca)[6])[2,i]*100, "%", sep = ""))+
    ylab(paste0("PC", i+1, " ", data.frame(summary(big_pca)[6])[2,i+1]*100, "%", sep = ""))+
    theme_minimal()+
    theme(legend.title = element_blank())+
    scale_color_manual(values=viridis::magma(length(unique(loop_data[,color_by]))+1))
  #ggrepel::geom_label_repel(aes_string(label = "sample_id"), show.legend = FALSE)
  
  plot_list[[i]] <- plt
  
}
  big_plot <- plot_list[[1]] | plot_list[[3]] | plot_list[[5]] | plot_list[[7]] | plot_list[[9]]
  ggsave(paste("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/many_pca_", color_by, ".png", sep=""), big_plot, height = 5, width=25, dpi=444, bg="white")
  
}



# model selection ####
  ## lm() covariate inclusion####
model_selection <- big_df2 %>%
  filter(timepoint!="day28", timepoint!="bad_baseline")%>%
  group_by(targetName)%>%
  nest() %>%
  mutate(simple_model=map(data,  ~lm(concentration~timepoint+class+id, data=.)), simple_AIC=map_dbl(simple_model, ~AIC(.)))%>%
  mutate(age_model=map(data,     ~lm(concentration~timepoint+class+ageyrs+id, data=.)), age_AIC=map_dbl(age_model, ~AIC(.)))%>%
  mutate(sex_model=map(data,     ~lm(concentration~timepoint+class+gender_categorical+id, data=.)), sex_AIC=map_dbl(sex_model, ~AIC(.)))%>%
  mutate(age_sex_model=map(data, ~lm(concentration~timepoint+class+ageyrs+gender_categorical+id, data=.)), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))%>%
  rowwise()%>%
  mutate(lowest=min(simple_AIC, age_AIC, sex_AIC, age_sex_AIC))%>%
  ungroup()%>%
  mutate(across(c(simple_AIC, age_AIC, sex_AIC, age_sex_AIC), ~ .x - lowest))

table(model_selection$simple_AIC==0 & model_selection$age_AIC>4)#13
table(model_selection$age_AIC==0 & model_selection$simple_AIC>4)#131
table(model_selection$sex_AIC==0 & model_selection$simple_AIC>4)#2
table(model_selection$age_sex_AIC==0 & model_selection$age_AIC>4)#80
# --> inlcude age as covariate

## lm() model structure####

model_selection <- big_df2 %>%
  filter(timepoint!="day28", timepoint!="bad_baseline")%>%
  group_by(targetName)%>%
  nest() %>%
  mutate(additive_model=map(data,  ~lm(concentration~timepoint+class+ageyrs+id, data=.)), additive_model_AIC=map_dbl(additive_model, ~AIC(.)))%>%
  mutate(age_interact=map(data,    ~lm(concentration~timepoint+class*ageyrs+id, data=.)), age_interact_AIC=map_dbl(age_interact, ~AIC(.)))%>%
  mutate(class_interact=map(data,     ~lm(concentration~timepoint*class+ageyrs+id, data=.)), class_interact_AIC=map_dbl(class_interact, ~AIC(.)))%>%
  # mutate(age_sex_model=map(data, ~lm(concentration~timepoint+class+ageyrs+gender_categorical+id, data=.)), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))%>%
  rowwise()%>%
  mutate(lowest=min(additive_model_AIC, age_interact_AIC, class_interact_AIC))%>%
  ungroup()%>%
  mutate(across(c(additive_model_AIC, age_interact_AIC, class_interact_AIC), ~ .x - lowest))

table(model_selection$additive_model_AIC==0)#129
table(model_selection$age_interact_AIC==0 & model_selection$additive_model_AIC>4)#4
table(model_selection$class_interact_AIC==0 & model_selection$additive_model_AIC>4)#66
table(model_selection$age_sex_AIC==0 & model_selection$age_AIC>4)#47

# --> simple additive model is best


## lmer () ####
mixed_model_selection <- big_df2 %>%
  filter(timepoint!="day28", timepoint!="bad_baseline")%>%
  group_by(targetName)%>%
  mutate(log_qpcr=log10(qpcr+0.1))%>%
  nest() %>%
  mutate(simple_model=map(data,  ~lme4::lmer(concentration~timepoint+class+(1|id), data=.)), simple_AIC=map_dbl(simple_model, ~AIC(.)))%>%
  mutate(qpcr_model=map(data,     ~lme4::lmer(concentration~timepoint+class+log_qpcr+(1|id), data=.)), qpcr_AIC=map_dbl(qpcr_model, ~AIC(.)))%>%
  mutate(age_model=map(data,     ~lme4::lmer(concentration~timepoint+class+ageyrs+(1|id), data=.)), age_AIC=map_dbl(age_model, ~AIC(.)))%>%
  mutate(sex_model=map(data,     ~lme4::lmer(concentration~timepoint+class+gender_categorical+(1|id), data=.)), sex_AIC=map_dbl(sex_model, ~AIC(.)))%>%
  mutate(age_sex_model=map(data, ~lme4::lmer(concentration~timepoint+class+ageyrs+gender_categorical+(1|id), data=.)), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))%>%
  rowwise()%>%
  mutate("lowest"=min(simple_AIC, age_AIC, sex_AIC, age_sex_AIC, qpcr_AIC))%>%
  ungroup()%>%
  mutate(across(c(simple_AIC, age_AIC, sex_AIC, age_sex_AIC, qpcr_AIC), ~ .x - lowest))
  
table(mixed_model_selection$qpcr_AIC==0) #101
table(mixed_model_selection$simple_AIC==0) #97
table(mixed_model_selection$qpcr_AIC==0 & mixed_model_selection$simple_AIC>4) #79
table(mixed_model_selection$simple_AIC==0 & mixed_model_selection$qpcr_AIC>4) #48

table(mixed_model_selection$age_AIC==0 & mixed_model_selection$simple_AIC>4)#22
table(mixed_model_selection$sex_AIC==0 & mixed_model_selection$simple_AIC>4)#3
table(mixed_model_selection$age_sex_AIC==0 & mixed_model_selection$age_AIC>4)#1

# --> use simple mixed model
# > table(mixed_model_selection$sex_AIC==0 & mixed_model_selection$simple_AIC>4)
# FALSE  TRUE 
# 237    13 
# [1] "ANGPT1"  "CD40LG"  "CRP"     "CTF1"    "GZMB"    "IL12p70" "IL17F"   "IL1RN"   "IL27"    "IRAK4"   "MIF"     "PTX3"    "TNFSF14"

# > table(mixed_model_selection$age_AIC==0 & mixed_model_selection$simple_AIC>4)
# FALSE  TRUE
# 213    37 
# [1] "AGRP"          "BST2"          "CCL17"         "CCL28"         "CD70"          "CD80"          "CSF2"          "CTLA4"        
# [9] "CXCL2"         "EGF"           "FGF2"          "GDF15"         "GFAP"          "IFNA1; IFNA13" "IFNG"          "IKBKG"        
# [17] "IL17RB"        "IL22"          "IL2RA"         "IL33"          "IL7"           "KLRK1"         "LGALS9"        "MUC16"        
# [25] "NCR1"          "NGF"           "S100A9"        "SCG2"          "SELE"          "SELP"          "TEK"           "TGFB1"        
# [33] "THBS2"         "TIMP1"         "TNFSF15"       "VEGFD"         "VSNL1" 


# sandbox ####

mean_base_concs <- big_df2 %>%
  filter(timepoint=="baseline")%>%
  group_by(targetName)%>%
  summarise("mean"=mean(concentration), "var"=var(concentration))%>%
  ggplot(aes(x=mean, y=var))+
  geom_point()+
  geom_text(aes(label=targetName))
model_selection <- big_df2 %>%
  filter(timepoint!="day28", timepoint!="bad_baseline")%>%
  group_by(targetName)%>%
  nest() %>%
  mutate(simple_model=map(data,  ~lm(concentration~timepoint*class+id, data=.)), simple_AIC=map_dbl(simple_model, ~AIC(.)))%>%
  mutate(age_model=map(data,     ~lm(concentration~timepoint*class+ageyrs+id, data=.)), age_AIC=map_dbl(age_model, ~AIC(.)))%>%
  mutate(sex_model=map(data,     ~lm(concentration~timepoint*class+gender_categorical+id, data=.)), sex_AIC=map_dbl(sex_model, ~AIC(.)))%>%
  mutate(age_sex_model=map(data, ~lm(concentration~timepoint*class+ageyrs+gender_categorical+id, data=.)), age_sex_AIC=map_dbl(age_sex_model, ~AIC(.)))%>%
  rowwise()%>%
  mutate(lowest=min(simple_AIC, age_AIC, sex_AIC, age_sex_AIC))%>%
  ungroup()%>%
  mutate(simplest_is_best=lowest==simple_AIC,
         age_is_best=lowest==age_AIC,
         sex_is_best=lowest==sex_AIC,
         age_sex_is_best=lowest==age_sex_AIC)

model_selection$targetName[model_selection$simplest_is_best]
model_selection$targetName[model_selection$age_is_best]
model_selection$targetName[model_selection$sex_is_best]
model_selection$targetName[model_selection$age_sex_is_best]






#most variable features

topVarProts <- head(order(rowVars(cpm(d_norm)), decreasing = TRUE))



# big_df2 <- big_df %>%
#   mutate(id=substr(sample_id, 2,4),
#          class=substr(sample_id, 6,6),
#          barcode= substr(sample_id, nchar(sample_id)-5, nchar(sample_id)))%>%
#   rowwise()%>%
#   mutate("sample_type"=substr(sample_id,
#                               stringi::stri_locate_all(pattern = "_", sample_id, fixed=TRUE)[[1]][1]+1,
#                               stringi::stri_locate_all(pattern = "_", sample_id, fixed=TRUE)[[1]][2]-1))%>%
#   ungroup()%>%
#   mutate("timepoint"=case_when(sample_type=="A.1"~"baseline",
#                                sample_type=="A.2"~"drop",
#                                sample_type=="A0"~"day0",
#                                sample_type=="A14"~"day14",
#                                sample_type=="A2.1"~"drop",
#                                sample_type=="A20"~"day0",
#                                sample_type=="A214"~"day14",
#                                sample_type=="A28"~"day28",
#                                sample_type=="NM.1"~"baseline",
#                                sample_type=="NM0"~"day0",
#                                sample_type=="NM14"~"day14",
#                                sample_type=="NM7"~"day7",
#                                sample_type=="PS0"~"PS0",
#                                sample_type=="S.1"~"baseline",
#                                sample_type=="S.2"~"drop",
#                                sample_type=="S0"~"day0",
#                                sample_type=="S14"~"day14",
#                                sample_type=="S28"~"day28",
#                                sample_type=="S7"~"day7")
#   )%>%
#   mutate(timepoint=factor(timepoint, levels=c("baseline", "day0", "day7", "day14", "day28")))%>%
#   group_by(targetName) %>%
#   mutate(z_conc=scale(concentration, center = TRUE, scale = TRUE))%>%
#   ungroup()%>%
#   group_by(id)%>%
#   mutate("mean_z_conc"=median(z_conc))%>%
#   ungroup()%>%
#   filter(mean_z_conc > -0.35)
# parasitemia inclussion?


base_zero_contrast <- t(matrix(c(0,1,0,0)))
base_14_contrast <- t(matrix(c(0,0,1,0)))
# base_28_contrast <- t(matrix(c(0,0,0)))
zero_14_contrast <- t(matrix(c(0,-1,1,0)))

base_zero_para_contrast <- t(matrix(c(0,1,0,1)))
base_14_para_contrast <- t(matrix(c(0,0,1,1)))

asymp_only_purff <- big_df2 %>%
  filter(class=="A", timepoint!="day28", timepoint!="bad_baseline")%>%
  group_by(targetName)%>%
  mutate(log_pcr=log10(qpcr+0.1),
         log_pcr_dich=if_else(log_pcr>=2.69, "over 500k/ml", "under 500k/ml"))%>%
  nest() %>%
  # mutate(model=map(data, ~lm(concentration~timepoint+ageyrs+id, data=.))) %>%
  mutate(model=map(data, ~lme4::lmer(concentration~timepoint+log_pcr+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_zero=map(model, ~multcomp::glht(., base_zero_contrast)),
         base_zero_p=map_dbl(base_zero, ~summary(.)$test$pvalues)) %>%
  mutate(base_zero_para=map(model, ~multcomp::glht(., base_zero_para_contrast)),
         base_zero_para_p=map_dbl(base_zero_para, ~summary(.)$test$pvalues),
         base_zero_para_coef=map_dbl(summary, ~coef(.)[5]),
         base_zero_para_coef_p=map_dbl(summary, ~coef(.)[5]))%>%
  mutate(base_14_para=map(model, ~multcomp::glht(., base_14_para_contrast)),
         base_14_para_p=map_dbl(base_14_para, ~summary(.)$test$pvalues),
         base_14_para_coef=map_dbl(summary, ~coef(.)[5]),
         base_14_para_coef_p=map_dbl(summary, ~coef(.)[5]))%>%
  ungroup()%>%
  mutate(base_zero_padj=p.adjust(base_zero_p, method="BH"),
         base_14_para_padj=p.adjust(base_14_para_p, method="BH"),
         base_zero_para_padj=p.adjust(base_zero_para_p, method="BH")
  )



asymp_results_table <- asymp_only_purff %>%
  dplyr::select(targetName,
                base_zero_padj,
                base_zero_para_padj,
                base_zero_para_coef,
                base_14_para_padj)%>%
  ungroup()

#79 with lm(); 159 with lme4;
asymp_sig_base_zero <- asymp_results_table %>%
  filter(base_zero_padj<fdr_cutoff)


#0 with lm(); 36 with lme4
asymp_sig_base_0_para <- asymp_only_purff %>%
  filter(base_zero_para_padj<fdr_cutoff)%>%
  arrange(desc(base_zero_para_coef))

asymp_sig_base_14_para <- asymp_only_purff %>%
  filter(base_14_para_padj<fdr_cutoff)%>%
  arrange(desc(base_14_para_padj))


#2 with lme4
symp_sig_base_14 <- symp_results_table %>%
  filter(base_14_padj<fdr_cutoff)%>%
  select(targetName)




big_df2 %>%
  filter(class=="A", timepoint!="day28", timepoint!="bad_baseline")%>%
  filter(targetName %in% head(asymp_sig_base_14_para$targetName, n=29))%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  mutate(log_pcr=log10(qpcr+0.1),
         log_pcr_dich=if_else(log_pcr>=2.69, "over 500k/ml", "under 500k/ml"))%>%
  ggplot(aes(x=log_pcr, y=concentration))+
  geom_point(aes(color=timepoint))+#
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman", color="red")+
  # geom_line(aes(group=id, color=qpcr_cat))+
  # geom_boxplot(aes(fill=log_pcr_dich), outliers = FALSE)+
  # geom_violin(draw_quantiles = seq(0,1,0.25))+
  # ggtitle("regulated during asymptomatic parasitemia")+
  facet_wrap(~targetName, scales = "free")+
  scale_color_manual(values=viridis::magma(6))+
  theme_minimal()

