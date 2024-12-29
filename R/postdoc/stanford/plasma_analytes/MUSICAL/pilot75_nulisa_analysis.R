library(purrr)
library(tidyr)
library(dplyr)
`%notin%` <- Negate(`%in%`)

musical_metadata <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/MASTER_METADATA.csv")
random_codes <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/pilot/id_date_code.csv")
random_codes$plasma.barcode <- gsub("D1PN8A", "D1_PN8A", random_codes$plasma.barcode)
random_codes$plasma.barcode <- gsub("D1JLGS", "DIJLGS", random_codes$plasma.barcode)
random_codes$plasma.barcode <- gsub("D1KWT2", "D1_KWT2", random_codes$plasma.barcode)
random_codes$plasma.barcode <- gsub("D1EF4A", "DEF4A", random_codes$plasma.barcode)


micdrop_codes <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/micdrop_nulisa_sample_codes.csv")
micdrop_codes$study <- "MICDROP"
micdrop_codes$plasma.barcode <- gsub("3YHJI", "X3YHJI", micdrop_codes$plasma.barcode)
micdrop_codes$plasma.barcode <- gsub("8UPHI", "X8UPHI", micdrop_codes$plasma.barcode)
#careful!!! not sure what the 426NI sample is
micdrop_codes$plasma.barcode <- gsub("QF9VI", "X426NI", micdrop_codes$plasma.barcode)

# subset(colnames(nulisa),  colnames(nulisa)%notin% c(micdrop_codes$plasma.barcode, random_codes$plasma.barcode))


slim_musical_metadata <- musical_metadata %>%
  mutate(day_annotation=if_else(day_annotation==84, -1, day_annotation))%>%
  # select(combined_id, combined_date, enrolltype, day_annotation)%>%
  mutate(id=combined_id, date=combined_date, class=enrolltype, timepoint=paste("t", day_annotation, sep=""))%>%
  select(-combined_id, -combined_date, -enrolltype, -day_annotation)%>%
  mutate("study"="MUSICAL")

combo_frame <- merge(slim_musical_metadata, random_codes, by=c("id", "date"))
combo_frame2 <- rbind(combo_frame, micdrop_codes)

nulisa <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/nulisa_data.csv")

wide_nulisa <- nulisa %>%
  pivot_longer(cols = colnames(nulisa)[2:ncol(nulisa)], names_to = "plasma.barcode", values_to = "concentration")
  

wide_nulisa <- inner_join(wide_nulisa, combo_frame2, by="plasma.barcode")

wide_nulisa <- wide_nulisa %>%
  mutate("time_class"=paste(class, timepoint, sep='_'),
         "age_class"=if_else(.$id %in% c(268, 324, 137, 176, 353, 161, 363, 571, 10766, 10794, 10842), "child", "adult"))%>%
  mutate(timepoint=factor(timepoint, levels=c("t-1", "t0", "t7", "t14")),
         age_class=factor(age_class, levels=c("child", "adult")),
         time_class=factor(time_class, levels=c("A_t-1",
                                                "A_t0",
                                                "A_t14",
                                                "S_t-1",
                                                "S_t0",
                                                "S_t7",
                                                "S_t14",
                                                "A_8", 
                                                "A_24", 
                                                "A_52",
                                                "S_24", 
                                                "A_pool")))%>%
  group_by(targetName) %>%
  mutate(z_conc=scale(concentration, center = TRUE, scale = TRUE))
  
         

# n=250, consistently, across 71 smaples :)
# wide_nulisa %>%
#   group_by(id, timepoint, class)%>%
#   summarise("n"=n())%>%
#   print(n=71)

# linear regression ####
  ## symp data through time ####
#child contrasts with day7
# base_zero_contrast <- t(matrix(c(0,1,0,0)))
# base_14_contrast <- t(matrix(c(0,0,1,0)))
# base_7_contrast <- t(matrix(c(0,0,0,1)))
# 
# zero_7_contrast <- t(matrix(c(0,-1,0,1)))
# zero_14_contrast <- t(matrix(c(0,-1,1,0)))
# 


#adult contrasts
base_zero_contrast <- t(matrix(c(0,1,0,0,1)))
base_14_contrast <- t(matrix(c(0,0,0,1,0)))
base_7_contrast <- t(matrix(c(0,0,1,0,0)))

zero_14_contrast <- t(matrix(c(0,-1,0,1,-1)))

symp_only_purff <- wide_nulisa %>%
  filter(plasma.barcode %notin% c("D19E2G", "D1FK67", "D1SSPJ"))%>%
  filter(study=="MUSICAL")%>%
  filter(class=="S")%>%
  group_by(targetName) %>%
  nest() %>%
  mutate(model=map(data, ~lme4::lmer(concentration~factor(timepoint)+age_class+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_zero=map(model, ~multcomp::glht(., base_zero_contrast)),
         base_zero_p=map_dbl(base_zero, ~summary(.)$test$pvalues)) %>%
  mutate(base_14=map(model, ~multcomp::glht(., base_14_contrast)),
         base_14_p=map_dbl(base_14, ~summary(.)$test$pvalues)) %>%
  mutate(base_7=map(model, ~multcomp::glht(., base_7_contrast)),
         base_7_p=map_dbl(base_7, ~summary(.)$test$pvalues))%>%
  mutate(zero_14=map(model, ~multcomp::glht(., zero_14_contrast)),
         zero_14_p=map_dbl(zero_14, ~summary(.)$test$pvalues)) %>%
  ungroup()%>%
  mutate(base_zero_padj=p.adjust(base_zero_p, method="BH"),
         base_7_padj=p.adjust(base_7_p, method="BH"),
         base_14_padj=p.adjust(base_14_p, method="BH"),
         zero_14_padj=p.adjust(zero_14_p, method="BH")
  )

symp_results_table <- symp_only_purff %>%
  mutate(coef=map_dbl(model, ~coef(.)$id[1,2]+coef(.)$id[1,3]))%>%
  dplyr::select(targetName, coef,
                base_zero_padj,
                base_7_padj,
                base_14_padj,
                zero_14_padj)%>%
  ungroup()

fdr_cutoff=0.05
# WHEN REFERENCE GROUP WAS ADULT BASELINE
#80 when 0; 38 when -1 (ageChild); 48 when 1
#80 across the board (ignore age class)
#249 when intercept =-1
#250 when intercept =1

#WHEN REFERENCE GROUP WAS CHILD BASELINE
#80 when 0; when 1 (ageAdult); 38 when 1; 48 when -1
symp_sig_base_zero <- symp_results_table %>%
  filter(base_zero_padj<fdr_cutoff)%>%
  select(targetName)


sig_base_7 <- symp_results_table %>%
  filter(base_7_padj<fdr_cutoff)%>%
  select(targetName)

symp_sig_base_14 <- symp_results_table %>%
  filter(base_14_padj<fdr_cutoff)%>%
  select(targetName)

symp_sig_zero_14 <- symp_results_table %>%
  filter(zero_14_padj<fdr_cutoff)%>%
  select(targetName)


sig_symp_lineplot <- wide_nulisa %>%
  filter(study=="MUSICAL")%>%
  # filter(class=="S", age_class=="child")%>%
  filter(class=="S")%>%
  filter(targetName %in% symp_sig_base_zero$targetName)%>%
  ggplot(., aes(x=timepoint, y=concentration, fill=age_class))+
  ggtitle("differentially abundant proteins baseline vs. day 0 in symptomatic children")+
  # ggtitle("differentially abundant proteins baseline vs. day 0 in symptomatic adults")+
  geom_boxplot()+
  # geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  facet_wrap(~targetName, scales = "free")+
  theme_minimal()

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/figures/sig_symp_child_base_zero.png", sig_symp_lineplot, width=20, height=20, bg="white", dpi=444)
# ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/figures/sig_symp_adult_base_zero.png", sig_symp_lineplot, width=20, height=20, bg="white", dpi=444)

adult_targets <- sig_base_zero$targetName
# [1] "BDNF"        "BMP7"        "CCL11"       "CCL17"       "CCL19"       "CCL28"      
# [7] "CCL5"        "CD40LG"      "CLEC4A"      "CRP"         "CTF1"        "CXCL12"     
# [13] "CXCL2"       "CXCL3"       "CXCL5"       "EPO"         "ICOSLG"      "IFNB1"      
# [19] "IL17A"       "IL17A|IL17F" "IL17RB"      "IL2"         "IL23"        "IL36G"      
# [25] "IRAK4"       "LTA|LTB"     "MIF"         "MMP1"        "NGF"         "PDGFA"      
# [31] "SCG2"        "TAFA5"       "THPO"        "TNFSF11"     "TNFSF12"     "TNFSF13"    
# [37] "TNFSF15"     "TREM2"       "VEGFA"       "VSNL1"


child_targets <- sig_base_zero$targetName
# [1] "AGRP"     "AREG"     "C1QA"     "CALCA"    "CCL17"    "CCL2"     "CCL20"   
# [8] "CCL24"    "CCL4"     "CCL8"     "CRP"      "CX3CL1"   "CXCL10"   "CXCL11"  
# [15] "CXCL9"    "FGF21"    "FTH1"     "GDF15"    "GDF2"     "GZMA"     "GZMB"    
# [22] "IFNB1"    "IFNG"     "IL10"     "IL18BP"   "IL1B"     "IL1RL1"   "IL1RN"   
# [29] "IL27"     "IL4"      "IL6"      "LAG3"     "LCN2"     "LIF"      "LILRB2"  
# [36] "MICB"     "MMP12"    "NGF"      "PGF"      "PTX3"     "S100A12"  "TGFB3"   
# [43] "TNF"      "TNFRSF1A" "TNFRSF1B" "TNFSF11"  "TNFSF14"

# table(child_targets %in% adult_targets)
# FALSE  TRUE 
# 42     5 

# table(adult_targets %in% child_targets)
# FALSE  TRUE 
# 35     5


# wide_nulisa %>%
#   filter(targetName %in% sig_base_14$targetName)%>%
#   ggplot(., aes(x=timepoint, y=concentration))+
#   geom_boxplot()+
#   facet_wrap(~targetName)+
#   theme_minimal()
  ## asymp data through time ####


#asymp contrasts
base_zero_contrast <- t(matrix(c(0,1,0)))
base_14_contrast <- t(matrix(c(0,0,1)))
zero_14_contrast <- t(matrix(c(0,-1,1)))

asymp_only_purff <- wide_nulisa %>%
  filter(study=="MUSICAL", timepoint!="t7")%>%
  filter(class=="A")%>%
  group_by(targetName) %>%
  nest() %>%
  mutate(model=map(data, ~lme4::lmer(concentration~factor(timepoint)+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_zero=map(model, ~multcomp::glht(., base_zero_contrast)),
         base_zero_p=map_dbl(base_zero, ~summary(.)$test$pvalues)) %>%
  mutate(base_14=map(model, ~multcomp::glht(., base_14_contrast)),
         base_14_p=map_dbl(base_14, ~summary(.)$test$pvalues)) %>%
  # mutate(base_7=map(model, ~multcomp::glht(., base_7_contrast)),
  #        base_7_p=map_dbl(base_7, ~summary(.)$test$pvalues))%>%
  ungroup()%>%
  mutate(base_zero_padj=p.adjust(base_zero_p, method="BH"),
         # base_7_padj=p.adjust(base_7_p, method="BH"),
         base_14_padj=p.adjust(base_14_p, method="BH")
  )

asymp_results_table <- asymp_only_purff %>%
  mutate(coef=map_dbl(model, ~coef(.)$id[1,2]+coef(.)$id[1,3]))%>%
  dplyr::select(targetName, coef,
                base_zero_padj,
                # base_7_padj,
                base_14_padj)%>%
  ungroup()

asymp_sig_base_zero <- asymp_results_table %>%
  filter(base_zero_padj<fdr_cutoff)%>%
  select(targetName)

# sig_base_7 <- symp_results_table %>%
#   filter(base_7_padj<fdr_cutoff)%>%
#   select(targetName)

# sig_base_14 <- asymp_results_table %>%
#   filter(base_14_padj<fdr_cutoff)%>%
#   select(targetName)
# 

sig_asymp_line <- wide_nulisa %>%
  filter(study=="MUSICAL", timepoint!="t7")%>%
  filter(class=="A")%>%
  filter(targetName %in% asymp_sig_base_zero$targetName)%>%
  ggplot(., aes(x=timepoint, y=concentration))+
  ggtitle("asypmtomatic")+
  geom_boxplot()+
  facet_wrap(~targetName, scales = "free")+
  theme_minimal()


  ## compare baseline between classes ####
base_class_contrast <- t(matrix(c(0,1,0)))

base_only_purff <- wide_nulisa %>%
  filter(plasma.barcode %notin% c("D19E2G", "D1FK67", "D1SSPJ"))%>%
  filter(study=="MUSICAL", timepoint=="t-1")%>%
  group_by(targetName) %>%
  nest() %>%
  mutate(model=map(data, ~lme4::lmer(concentration~class+age_class+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_class=map(model, ~multcomp::glht(., base_class_contrast)),
         base_class_p=map_dbl(base_class, ~summary(.)$test$pvalues)) %>%
  ungroup()%>%
  mutate(base_class_padj=p.adjust(base_class_p, method="BH"))



base_only_results_table <- base_only_purff %>%
  mutate(coef=map_dbl(model, ~coef(.)$id[1,2]+coef(.)$id[1,3]))%>%
  dplyr::select(targetName, coef,
                base_class_padj)%>%
  ungroup()

sig_base_only <- base_only_purff %>%
  filter(base_class_p<0.05)%>%
  select(targetName)

# sig_base_7 <- symp_results_table %>%
#   filter(base_7_padj<fdr_cutoff)%>%
#   select(targetName)

# sig_base_14 <- asymp_results_table %>%
#   filter(base_14_padj<fdr_cutoff)%>%
#   select(targetName)
# 
wide_nulisa %>% 
  filter(plasma.barcode %notin% c("D19E2G", "D1FK67", "D1SSPJ"))%>%
  filter(study=="MUSICAL", timepoint=="t-1")%>%
  filter(targetName %in% sig_base_only$targetName)%>%
  ggplot(., aes(x=timepoint, y=concentration, fill=class))+
  ggtitle("baseline differences")+
  geom_boxplot()+
  facet_wrap(~targetName, scales = "free")+
  theme_minimal()

# complex heatmap ####
  ## symptomatic children ####


# inferno <- colorspace::sequential_hcl("inferno", n=21)
# col_inferno <- circlize::colorRamp2(seq(0,20, by=1), inferno)
t_nulisa <-t(as.matrix(wide_nulisa %>%
                         filter(study=="MUSICAL", timepoint %in% c("t-1", "t0"), age_class=="child", class=="S")%>%
                         mutate(sample_id=paste(id, time_class))%>%
                         select(sample_id, targetName, z_conc)%>%
                         pivot_wider(names_from = targetName, values_from = z_conc)))


colnames(t_nulisa) <- t_nulisa[1,]

num_nulisa_matrix <- t_nulisa[-1,]
class(num_nulisa_matrix) <- "numeric"

short_num_nulisa_matrix <- subset(num_nulisa_matrix, rownames(num_nulisa_matrix)%in%symp_sig_base_zero$targetName)

z_col_fun <- circlize::colorRamp2(c(min(short_num_nulisa_matrix), 0, abs(min(short_num_nulisa_matrix))), c("#123499", "black", "#FFA500"))

big_heatmap <- Heatmap(matrix = short_num_nulisa_matrix,
                       cluster_rows = TRUE,
                       cluster_columns=TRUE,
                       show_row_dend = FALSE,
                       show_column_dend = TRUE,
                       show_heatmap_legend = TRUE,
                       column_split = substr(colnames(t_nulisa), 5, nchar(colnames(t_nulisa))),
                       # name = "Pearson r",
                       #cluster_columns = FALSE,
                       column_names_gp = gpar(fontsize = 6),
                       row_names_gp = gpar(fontsize = 6),
                       row_names_side = "left",
                       col = z_col_fun,
                       column_names_rot = 45)

png("~/postdoc/stanford/plasma_analytes/MUSICAL/figures/sig_symp_child_base_heatmap.png", height=12, width = 5, units = "in", res=444)
draw(big_heatmap, padding=unit(c(2,2,2,2), "mm"))
dev.off()

  ## asymptomatic children ####

t_nulisa <-t(as.matrix(wide_nulisa %>%
                         filter(study=="MUSICAL", timepoint %in% c("t-1", "t0"), age_class=="child", class=="A")%>%
                         mutate(sample_id=paste(id, time_class))%>%
                         select(sample_id, targetName, z_conc)%>%
                         pivot_wider(names_from = targetName, values_from = z_conc)))


colnames(t_nulisa) <- t_nulisa[1,]

num_nulisa_matrix <- t_nulisa[-1,]
class(num_nulisa_matrix) <- "numeric"

short_num_nulisa_matrix <- subset(num_nulisa_matrix, rownames(num_nulisa_matrix)%in%asymp_sig_base_zero$targetName)

z_col_fun <- circlize::colorRamp2(c(min(short_num_nulisa_matrix), 0, abs(min(short_num_nulisa_matrix))), c("#123499", "black", "#FFA500"))

big_heatmap <- Heatmap(matrix = short_num_nulisa_matrix,
                       cluster_rows = TRUE,
                       cluster_columns=TRUE,
                       show_row_dend = FALSE,
                       show_column_dend = TRUE,
                       show_heatmap_legend = TRUE,
                       column_split = substr(colnames(t_nulisa), 5, nchar(colnames(t_nulisa))),
                       # name = "Pearson r",
                       #cluster_columns = FALSE,
                       column_names_gp = gpar(fontsize = 6),
                       row_names_gp = gpar(fontsize = 6),
                       row_names_side = "left",
                       col = z_col_fun,
                       column_names_rot = 45)

png("~/postdoc/stanford/plasma_analytes/MUSICAL/figures/sig_asymp_child_base_heatmap.png", height=12, width = 5, units = "in", res=444)
draw(big_heatmap, padding=unit(c(2,2,2,2), "mm"))
dev.off()


  ## BIG HEATMAP ####
t_nulisa <-t(as.matrix(wide_nulisa %>%
                         filter(study=="MUSICAL")%>%
                         mutate(sample_id=interaction(id, time_class))%>%
                         select(sample_id, targetName, z_conc)%>%
                         pivot_wider(names_from = targetName, values_from = z_conc)))


colnames(t_nulisa) <- t_nulisa[1,]

num_nulisa_matrix <- t_nulisa[-1,]
class(num_nulisa_matrix) <- "numeric"

short_num_nulisa_matrix <- num_nulisa_matrix

z_col_fun <- circlize::colorRamp2(c(min(short_num_nulisa_matrix)/2, 0, abs(min(short_num_nulisa_matrix)/2)), c("#123499", "black", "#FFA500"))

big_heatmap_split <- substr(colnames(t_nulisa), 5, nchar(colnames(t_nulisa)))

big_heatmap_split <- factor(big_heatmap_split, levels=sort(unique(big_heatmap_split))[c(1:5, 7, 6)])

big_heatmap <- Heatmap(matrix = short_num_nulisa_matrix,
                       cluster_rows = TRUE,
                       cluster_columns=FALSE,
                       show_row_dend = FALSE,
                       show_column_dend = FALSE,
                       show_heatmap_legend = TRUE,
                       column_split = big_heatmap_split,
                       name = "z score",
                       #cluster_columns = FALSE,
                       
                       column_names_gp = gpar(fontsize = 7),
                       row_names_gp = gpar(fontsize = 0),
                       row_names_side = "left",
                       col = z_col_fun,
                       column_names_rot = 90)

png("~/postdoc/stanford/plasma_analytes/MUSICAL/figures/BIG_HEATMAP.png", height=8, width = 8, units = "in", res=444)
draw(big_heatmap, padding=unit(c(2,2,2,2), "mm"))
dev.off()

# combined model fitting ####

#big purrrf contrasts
base_zero_adult_contrast <- t(matrix(c(0,1,0,0,0,0)))
base_zero_child_contrast <- t(matrix(c(0,1,0,0,0,1)))
base_zero_class_contrast <- t(matrix(c(0,1,0,0,1,0)))


big_purff <- wide_nulisa %>%
  filter(study=="MUSICAL")%>%
  group_by(targetName) %>%
  nest() %>%
  mutate(model=map(data, ~lme4::lmer(concentration~factor(timepoint)+class+age_class+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_zero_adult=map(model, ~multcomp::glht(., base_zero_adult_contrast)),
         base_zero_adult_p=map_dbl(base_zero_adult, ~summary(.)$test$pvalues)) %>%
  mutate(base_zero_child=map(model, ~multcomp::glht(., base_zero_child_contrast)),
         base_zero_child_p=map_dbl(base_zero_child, ~summary(.)$test$pvalues)) %>%
  mutate(base_zero_class=map(model, ~multcomp::glht(., base_zero_class_contrast)),
         base_zero_class_p=map_dbl(base_zero_class, ~summary(.)$test$pvalues))%>%
  ungroup()%>%
  mutate(base_zero_adult_padj=p.adjust(base_zero_adult_p, method="BH"),
         base_zero_child_padj=p.adjust(base_zero_child_p, method="BH"),
         base_zero_class_padj=p.adjust(base_zero_class_p, method="BH")
  )

sig_base_zero_adult_padj <- filter(big_purff, base_zero_adult_padj <0.05)
sig_base_zero_child_padj <- filter(big_purff, base_zero_child_padj <0.05)
sig_base_zero_class_padj <- filter(big_purff, base_zero_class_padj <0.05)

wide_nulisa %>%
  filter(study=="MUSICAL")%>%
  filter(timepoint %in% c("t-1", "t0"))%>%
  filter(targetName %in% sig_base_zero_class_padj$targetName)%>%
  ggplot(., aes(x=timepoint, y=concentration))+
  # geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  ggtitle("differentially abundant protein in adults")+
  geom_boxplot(aes(fill=age_class))+
  geom_point(aes(color=id, group=age_class), position = position_dodge(width=0.75))+
  # geom_line(aes(group=id))+
  facet_wrap(class~targetName, scales = "free")+
  theme_minimal()



t_nulisa <-t(as.matrix(wide_nulisa %>%
                         filter(study=="MUSICAL", timepoint %in% c("t-1", "t0"), age_class=="child", class=="A")%>%
                         mutate(sample_id=paste(id, time_class))%>%
                         select(sample_id, targetName, z_conc)%>%
                         pivot_wider(names_from = targetName, values_from = z_conc)))


colnames(t_nulisa) <- t_nulisa[1,]

num_nulisa_matrix <- t_nulisa[-1,]
class(num_nulisa_matrix) <- "numeric"

short_num_nulisa_matrix <- subset(num_nulisa_matrix, rownames(num_nulisa_matrix)%in%asymp_sig_base_zero$targetName)

z_col_fun <- circlize::colorRamp2(c(min(short_num_nulisa_matrix), 0, abs(min(short_num_nulisa_matrix))), c("#123499", "black", "#FFA500"))

big_heatmap <- Heatmap(matrix = short_num_nulisa_matrix,
                       cluster_rows = TRUE,
                       cluster_columns=TRUE,
                       show_row_dend = FALSE,
                       show_column_dend = TRUE,
                       show_heatmap_legend = TRUE,
                       column_split = substr(colnames(t_nulisa), 5, nchar(colnames(t_nulisa))),
                       # name = "Pearson r",
                       #cluster_columns = FALSE,
                       column_names_gp = gpar(fontsize = 6),
                       row_names_gp = gpar(fontsize = 6),
                       row_names_side = "left",
                       col = z_col_fun,
                       column_names_rot = 45)

png("~/postdoc/stanford/plasma_analytes/MUSICAL/figures/sig_asymp_child_base_heatmap.png", height=12, width = 5, units = "in", res=444)
draw(big_heatmap, padding=unit(c(2,2,2,2), "mm"))
dev.off()


## compare baseline between classes ####
base_class_contrast <- t(matrix(c(0,1,0)))

base_only_purff <- wide_nulisa %>%
  filter(plasma.barcode %notin% c("D19E2G", "D1FK67", "D1SSPJ"))%>%
  filter(study=="MUSICAL", timepoint=="t-1")%>%
  group_by(targetName) %>%
  nest() %>%
  mutate(model=map(data, ~lme4::lmer(concentration~class+age_class+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_class=map(model, ~multcomp::glht(., base_class_contrast)),
         base_class_p=map_dbl(base_class, ~summary(.)$test$pvalues)) %>%
  ungroup()%>%
  mutate(base_class_padj=p.adjust(base_class_p, method="BH"))



base_only_results_table <- base_only_purff %>%
  mutate(coef=map_dbl(model, ~coef(.)$id[1,2]+coef(.)$id[1,3]))%>%
  dplyr::select(targetName, coef,
                base_class_padj)%>%
  ungroup()

sig_base_only <- base_only_purff %>%
  filter(base_class_p<0.05)%>%
  select(targetName)

# sig_base_7 <- symp_results_table %>%
#   filter(base_7_padj<fdr_cutoff)%>%
#   select(targetName)

# sig_base_14 <- asymp_results_table %>%
#   filter(base_14_padj<fdr_cutoff)%>%
#   select(targetName)
# 
wide_nulisa %>% 
  filter(plasma.barcode %notin% c("D19E2G", "D1FK67", "D1SSPJ"))%>%
  filter(study=="MUSICAL", timepoint=="t-1")%>%
  filter(targetName %in% sig_base_only$targetName)%>%
  ggplot(., aes(x=timepoint, y=concentration, fill=class))+
  ggtitle("baseline differences")+
  geom_boxplot()+
  facet_wrap(~targetName, scales = "free")+
  theme_minimal()

wide_nulisa %>% 
  filter(plasma.barcode %notin% c("D19E2G", "D1FK67", "D1SSPJ"))%>%
  filter(study=="MUSICAL")%>%
  # filter(targetName %in% sig_base_only$targetName)%>%
  ggplot(., aes(x=concentration))+
  ggtitle("baseline differences")+
  geom_histogram()+
  facet_wrap(~targetName, scales = "free")+
  theme_minimal()





# incorporating parasitemia ####

fdr_cutoff=0.05

base_zero_contrast <- t(matrix(c(0,1,0)))
base_zero_para_contrast <- t(matrix(c(0,1,1)))

symp_only_purff <- wide_nulisa %>%
  filter(plasma.barcode %notin% c("D19E2G", "D1FK67", "D1SSPJ"))%>%
  mutate("string_id"=paste("s", id, sep=''))%>%
  filter(study=="MUSICAL", class=="S", age_class=="child", timepoint %in% c("t-1", "t0"))%>%
  mutate(log_qpcr=log10(qpcr+0.1))%>%
  group_by(targetName) %>%
  nest() %>%
  mutate(model=map(data, ~lme4::lmer(concentration~factor(timepoint)+log_qpcr+(1|string_id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_zero=map(model, ~multcomp::glht(., base_zero_contrast)),
         base_zero_p=map_dbl(base_zero, ~summary(.)$test$pvalues)) %>%
  mutate(base_zero_para=map(model, ~multcomp::glht(., base_zero_para_contrast)),
         base_zero_para_p=map_dbl(base_zero_para, ~summary(.)$test$pvalues)) %>%
  ungroup()%>%
  mutate(base_zero_padj=p.adjust(base_zero_p, method="BH"),
         base_zero_para_padj=p.adjust(base_zero_para_p, method="BH")
  )

symp_results_table <- symp_only_purff %>%
  mutate(time_coef=map(model, ~stats::coef(.)$string_id[1,2]),
         qpcr_coef=map(model, ~coef(.)$string_id[1,3])
  )%>%
  dplyr::select(targetName, time_coef, qpcr_coef,
                base_zero_padj,
                base_zero_para_padj)%>%
  ungroup()



symp_sig_base_zero <- symp_results_table %>%
  filter(base_zero_padj<fdr_cutoff)%>%
  select(targetName)

symp_sig_base_para_zero <- symp_results_table %>%
  filter(base_zero_para_padj<fdr_cutoff)%>%
  select(targetName)

wide_nulisa %>%
  filter(study=="MUSICAL", timepoint %in% c("t-1", "t0", "t7"))%>%
  # filter(targetName %in% symp_sig_base_para_zero$targetName[1:49])%>%
  filter(targetName %in% c("IL10","CCL2"))%>%
  mutate(log_qpcr=log10(qpcr+1))%>%
  ggplot(., aes(x=timepoint, y=concentration))+
  # ggtitle("top 20 DA proteins baseline vs. day 0 of symptomatic malaria across all individuals")+
  # ggtitle("differentially abundant proteins baseline vs. day 0 in symptomatic adults")+
  geom_point(size=8, aes(fill=log_qpcr, shape=class))+
  # geom_smooth(method="lm")+
  ggpubr::stat_cor(label.x = 1)+
  # geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  facet_wrap(~targetName, scales = "free")+
  # scale_fill_manual(values = symp_time_palette)+
  scale_shape_manual(values = c(21, 24))+
  theme_minimal()+
  symp_time_theme


wide_nulisa %>%
  filter(study=="MUSICAL", timepoint %in% c("t-1", "t0"))%>%
  # filter(targetName %in% symp_sig_base_para_zero$targetName[1:49])%>%
  filter(targetName %in% c("KDR", "CCL20", "IL1R1", "IL6R"))%>%
  mutate(log_qpcr=log10(qpcr+1))%>%
  ggplot(., aes(x=log_qpcr, y=concentration))+
  # ggtitle("top 20 DA proteins baseline vs. day 0 of symptomatic malaria across all individuals")+
  # ggtitle("differentially abundant proteins baseline vs. day 0 in symptomatic adults")+
  geom_point(size=4, aes(fill=timepoint), shape=21)+
  # geom_smooth(method="lm")+
  ggpubr::stat_cor(label.x = 1)+
  # geom_boxplot()+
  # geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  facet_wrap(~targetName, scales = "free")+
  # scale_fill_manual(values = symp_time_palette)+
  # scale_shape_manual(values = c(21, 24))+
  theme_minimal()+
  symp_time_theme
