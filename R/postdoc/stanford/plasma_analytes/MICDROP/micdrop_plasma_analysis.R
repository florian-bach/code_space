library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(ComplexHeatmap)

`%notin%` <- Negate(`%in%`)

fdr_cutoff=0.01


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
micdrop_codes$timepoint=paste(micdrop_codes$timepoint, "w", sep="")


slim_musical_metadata <- musical_metadata %>%
  mutate(day_annotation=if_else(day_annotation==84, -1, day_annotation))%>%
  select(combined_id, combined_date, enrolltype, day_annotation)%>%
  mutate(id=combined_id, date=combined_date, class=enrolltype, timepoint=paste("t", day_annotation, sep=""))%>%
  select(-combined_id, -combined_date, -enrolltype, -day_annotation)%>%
  mutate("study"="MUSICAL")

combo_frame <- merge(slim_musical_metadata, random_codes, by=c("id", "date"))
combo_frame2 <- rbind(combo_frame, micdrop_codes)

nulisa <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/pilot/nulisa_data.csv")

wide_nulisa <- nulisa %>%
  pivot_longer(cols = colnames(nulisa)[2:ncol(nulisa)], names_to = "plasma.barcode", values_to = "concentration")


wide_nulisa <- inner_join(wide_nulisa, combo_frame2, by="plasma.barcode")
wide_nulisa <- wide_nulisa %>%
  mutate("time_class"=paste(class, timepoint, sep='_'),
         "age_class"=if_else(.$id %in% c(268, 324, 137, 176, 353, 161, 363, 571, 10766, 10794, 10842), "child", "adult"))%>%
  mutate(timepoint=factor(timepoint, levels=c("t-1", "t0", "t7", "t14", "8w", "24w", "52w", "poolw")),
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

# statistics ####

musical_for_stats <- wide_nulisa %>%
  filter(plasma.barcode %notin% c("D19E2G", "D1FK67", "D1SSPJ"))%>%
  filter(study=="MUSICAL", timepoint=="t-1", age_class=="child")

micdrop_for_stats <-  wide_nulisa %>%
  dplyr::filter(study=="MICDROP",  age_class=="child")

micdrop_vs_musical <- rbind(musical_for_stats, micdrop_for_stats)

study_contrast <- t(matrix(c(0,1)))

mic_vs_mus_purf <- micdrop_vs_musical %>%
  group_by(targetName) %>%
  nest() %>%
  mutate(model=map(data, ~lme4::lmer(concentration~study+(1|id), data=.)))%>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(study_effect=map(model, ~multcomp::glht(., study_contrast)),
         study_effect_p=map_dbl(study_effect, ~summary(.)$test$pvalues)) %>%
  ungroup()%>%
  mutate(study_effect_padj=p.adjust(study_effect_p, method="BH"))

study_results_table <- mic_vs_mus_purf %>%
  dplyr::select(targetName, study_effect_padj)

top20_age <- symp_results_table%>%
  top_n(-60, wt=study_effect_padj)%>%
  arrange(study_effect_padj)


# stats top20 ####
plottable_micdrop <- wide_nulisa %>%
  filter(study=="MICDROP", age_class=="child")%>%
  filter(targetName %in% top20_age$targetName[1:20])

plottable_musical <- wide_nulisa %>%
  filter(study=="MUSICAL", timepoint=="t-1")%>%
  filter(targetName %in% top20_age$targetName[1:20])%>%
  mutate(timepoint=if_else(age_class=="child", "MUSICAL\nchildren", "MUSICAL\nadults"),
         timepoint=factor(timepoint, levels=c("MUSICAL\nchildren", "MUSICAL\nadults")))


top20_targets <- ggplot()+
  geom_point(aes(x=timepoint, y=concentration, color=id), data = plottable_micdrop)+
  geom_line(aes(x=timepoint, y=concentration, color=id, group=id), data = plottable_micdrop)+
  geom_violin(aes(x=timepoint, y=concentration), draw_quantiles = seq(0,1,0.25), data = plottable_musical)+
  facet_wrap(~targetName, scales = "free")+
  scale_color_manual(values=c("darkred", "#708090", "#6082B6"))+
  theme_minimal()+
  symp_time_theme+
  theme(legend.position="none", 
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8),
        axis.text.y=element_text(size=8))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/top20_targets.png", top20_targets, bg="white", width=8, height=8, dpi=444)


# stats top20-40 ####
plottable_micdrop <- wide_nulisa %>%
  filter(study=="MICDROP", age_class=="child")%>%
  filter(targetName %in% top20_age$targetName[21:40])

plottable_musical <- wide_nulisa %>%
  filter(study=="MUSICAL", timepoint=="t-1")%>%
  filter(targetName %in% top20_age$targetName[21:40])%>%
  mutate(timepoint=if_else(age_class=="child", "MUSICAL\nchildren", "MUSICAL\nadults"),
         timepoint=factor(timepoint, levels=c("MUSICAL\nchildren", "MUSICAL\nadults")))


top20_40_targets <- ggplot()+
  geom_point(aes(x=timepoint, y=concentration, color=id), data = plottable_micdrop)+
  geom_line(aes(x=timepoint, y=concentration, color=id, group=id), data = plottable_micdrop)+
  geom_violin(aes(x=timepoint, y=concentration), draw_quantiles = seq(0,1,0.25), data = plottable_musical)+
  facet_wrap(~targetName, scales = "free")+
  scale_color_manual(values=c("darkred", "#708090", "#6082B6"))+
  theme_minimal()+
  symp_time_theme+
  theme(legend.position="none", 
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8),
        axis.text.y=element_text(size=8))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/top20_40_targets.png", top20_40_targets, bg="white", width=8, height=8, dpi=444)




# stats top40-60 ####
plottable_micdrop <- wide_nulisa %>%
  filter(study=="MICDROP", age_class=="child")%>%
  filter(targetName %in% top20_age$targetName[41:60])

plottable_musical <- wide_nulisa %>%
  filter(study=="MUSICAL", timepoint=="t-1")%>%
  filter(targetName %in% top20_age$targetName[41:60])%>%
  mutate(timepoint=if_else(age_class=="child", "MUSICAL\nchildren", "MUSICAL\nadults"),
         timepoint=factor(timepoint, levels=c("MUSICAL\nchildren", "MUSICAL\nadults")))


top40_60_targets <- ggplot()+
  geom_point(aes(x=timepoint, y=concentration, color=id), data = plottable_micdrop)+
  geom_line(aes(x=timepoint, y=concentration, color=id, group=id), data = plottable_micdrop)+
  geom_violin(aes(x=timepoint, y=concentration), draw_quantiles = seq(0,1,0.25), data = plottable_musical)+
  facet_wrap(~targetName, scales = "free")+
  scale_color_manual(values=c("darkred", "#708090", "#6082B6"))+
  theme_minimal()+
  symp_time_theme+
  theme(legend.position="none", 
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8),
        axis.text.y=element_text(size=8))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/top40_60_targets.png", top40_60_targets, bg="white", width=8, height=8, dpi=444)



# malaria targets ####

malaria_targets_to_plot <- c("CXCL10", "IL10", "IFNG", "TNF", "CRP", "CALCA")

plottable_micdrop <- wide_nulisa %>%
  filter(study=="MICDROP", age_class=="child")%>%
  filter(targetName %in% malaria_targets_to_plot)

plottable_musical <- wide_nulisa %>%
  filter(study=="MUSICAL", timepoint=="t-1")%>%
  filter(targetName %in% malaria_targets_to_plot)%>%
  mutate(timepoint=if_else(age_class=="child", "MUSICAL\nchildren", "MUSICAL\nadults"),
         timepoint=factor(timepoint, levels=c("MUSICAL\nchildren", "MUSICAL\nadults")))


malaria_targets <- ggplot()+
  geom_point(aes(x=timepoint, y=concentration, color=id), data = plottable_micdrop)+
  geom_line(aes(x=timepoint, y=concentration, color=id, group=id), data = plottable_micdrop)+
  geom_violin(aes(x=timepoint, y=concentration), draw_quantiles = seq(0,1,0.25), data = plottable_musical)+
  facet_wrap(~targetName, scales = "free")+
  scale_color_manual(values=c("darkred", "#708090", "#6082B6"))+
  theme_minimal()+
  symp_time_theme+
  theme(legend.position="none", 
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/malaria_targets.png", malaria_targets, bg="white", width=8, height=5, dpi=444)

# child vs. adult targets ####

age_targets_to_plot <- c("IL17B", "MMP3", "SPP1", "CTLA4", "IL12B", "IL20", "IL2RA", "IL34", "TNFSF11")

plottable_micdrop <- wide_nulisa %>%
  filter(study=="MICDROP", age_class=="child")%>%
  filter(targetName %in% age_targets_to_plot)

plottable_musical <- wide_nulisa %>%
  filter(study=="MUSICAL", timepoint=="t-1")%>%
  filter(targetName %in% age_targets_to_plot)%>%
  mutate(timepoint=if_else(age_class=="child", "MUSICAL\nchildren", "MUSICAL\nadults"),
         timepoint=factor(timepoint, levels=c("MUSICAL\nchildren", "MUSICAL\nadults")))


age_targets_to_plot <- ggplot()+
  geom_point(aes(x=timepoint, y=concentration, color=id), data = plottable_micdrop)+
  geom_line(aes(x=timepoint, y=concentration, color=id, group=id), data = plottable_micdrop)+
  geom_violin(aes(x=timepoint, y=concentration), draw_quantiles = seq(0,1,0.25), data = plottable_musical)+
  facet_wrap(~targetName, scales = "free")+
  scale_color_manual(values=c("darkred", "#708090", "#6082B6"))+
  theme_minimal()+
  symp_time_theme+
  theme(legend.position="none", 
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/age_targets_to_plot.png", age_targets_to_plot, bg="white", width=8, height=8, dpi=444)


# t cell polarisers

t_cell_polarisers <- c("IL12B", "IL12p70",
                       "IL4", "TNF", "TGFB1", "IL1B", 
                       "IL6", "IL23"
                       )

plottable_micdrop <- wide_nulisa %>%
  filter(study=="MICDROP", age_class=="child")%>%
  filter(targetName %in% t_cell_polarisers)

plottable_musical <- wide_nulisa %>%
  filter(study=="MUSICAL", timepoint=="t-1")%>%
  filter(targetName %in% t_cell_polarisers)%>%
  mutate(timepoint=if_else(age_class=="child", "MUSICAL\nchildren", "MUSICAL\nadults"),
         timepoint=factor(timepoint, levels=c("MUSICAL\nchildren", "MUSICAL\nadults")))


tcell_targets_to_plot <- ggplot()+
  geom_point(aes(x=timepoint, y=concentration, color=id), data = plottable_micdrop)+
  geom_line(aes(x=timepoint, y=concentration, color=id, group=id), data = plottable_micdrop)+
  geom_violin(aes(x=timepoint, y=concentration), draw_quantiles = seq(0,1,0.25), data = plottable_musical)+
  facet_wrap(~targetName, scales = "free", nrow=2)+
  scale_color_manual(values=c("darkred", "#708090", "#6082B6"))+
  theme_minimal()+
  symp_time_theme+
  theme(legend.position="none", 
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8),
        axis.text.y=element_text(size=8))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/tcell_targets_to_plot.png", tcell_targets_to_plot, bg="white", width=8, height=5, dpi=444)






# wide_nulisa %>%
#   filter(study=="MICDROP", timepoint!="poolw")%>%
#   filter(targetName %in% symp_plottable_targets_base_0$targetName)%>%
#   ggplot(., aes(x=factor(timepoint), y=concentration, color=id, group=id))+
#   # ggtitle("top 20 DA proteins baseline vs. day 0 of symptomatic malaria across all individuals")+
#   # ggtitle("differentially abundant proteins baseline vs. day 0 in symptomatic adults")+
#   # geom_boxplot()+
#   # geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
#   geom_point()+
#   geom_line()+
#   facet_wrap(~targetName, scales = "free")+
#   scale_color_manual(values=c("darkred", "#708090", "#6082B6"))+
#   theme_minimal()+
#   symp_time_theme


# DPSP pool analysis ####

musical_baselines <- wide_nulisa %>%
  filter(timepoint=="t-1", study=="MUSICAL")

pool_data <- filter(wide_nulisa, id=="DPSP_pool")

plot_data <- rbind(musical_baselines, pool_data)

plot_data <- plot_data%>%
  group_by(targetName)%>%
  mutate(z_conc = scale(concentration, center = TRUE, scale = TRUE))%>%
  ungroup()

dpsp_specific <- plot_data %>%
  filter(id=="DPSP_pool")%>%
  filter(abs(z_conc)>2)
  # slice_max(abs(z_conc), n=36)
  
dpsp_pool <- filter(plot_data, id=="DPSP_pool", targetName %in% dpsp_specific$targetName)
sans_pool <- filter(plot_data, id!="DPSP_pool", targetName %in% dpsp_specific$targetName)
# pool_data <- filter(wide_nulisa, id=="DPSP_pool", targetName %in% malaria_targets_to_plot)



pool_plot <- ggplot()+
  geom_violin(aes(x=age_class, y=z_conc), draw_quantiles = seq(0,1,0.25), data = sans_pool)+
  geom_point(aes(x=age_class, y=z_conc), data = dpsp_pool)+
  facet_wrap(~targetName, scales = "free")+
  scale_color_manual(values=c("darkred", "#708090", "#6082B6"))+
  theme_minimal()+
  # symp_time_theme+
  theme(legend.position="none",
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8),
        axis.text.y=element_text(size=8))

ggsave("~/postdoc/stanford/plasma_analytes/dpsp/dpsp_pool_plot.png", height = 40, width = )


# single target visualiser ####
# 
# symp_time_palette <- c("white", "darkred", "orange")
# names(symp_time_palette) <- c("t-1", "t0", "t14")
# 
# asymp_time_palette <- c("lightgrey", "#444444", "#222")
# names(asymp_time_palette) <- c("t-1", "t0", "t14")
# 
# time_class_palette <- c(unname(symp_time_palette[1:2]), asymp_time_palette[1:2])
# names(time_class_palette) <- c("t-1.S", "t0.S", "t-1.A", "t0.A")
# 
# wide_nulisa %>%
#   filter(study=="MUSICAL")%>%
#   # filter(age_class=="child")%>%
#   filter(timepoint %in% c("t-1", "t0"))%>%
#   filter(targetName=="LAG3")%>%
#   ggplot(., aes(x=age_class, y=concentration, fill=interaction(timepoint, class)))+
#   # ggtitle("top20 of 50 proteins DA across all children")+
#   # ggtitle("differentially abundant proteins baseline vs. day 0 in symptomatic adults")+
#   geom_boxplot()+
#   # geom_point(shape=21)+
#   # geom_path()+
#   # geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
#   facet_wrap(~targetName, scales = "free")+
#   scale_fill_manual(values = (time_class_palette))+
#   theme_minimal()+
#   symp_time_theme+
#   theme(legend.title = element_blank())

simple_metadata <- wide_nulisa%>%
  filter(study=="MICDROP")%>%
  ungroup()%>%
  distinct(plasma.barcode, id, timepoint)

write.csv("~/postdoc/stanford/plasma_analytes/MICDROP/pilot_metadata.csv")
