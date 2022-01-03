library(tidyr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

`%notin%` <- Negate(`%in%`)



vol_pal <- rev(c("#FFC800",
                 "#EE000C",
                 "#FF8000",
                 "#0080FF",
                 "#009B95",
                 "#E54787",
                 "#8000FF",
                 "#66CCFF",
                 "#4B0055"))

names(vol_pal) <- c("v313", "v315", "v320", "v306", "v301", "v308", "v305", "v304", "v310")


vivax_colours <- list("v02" = "#FB9A99",
                      "v03" = "#E31A1C",
                      "v05" = "#A6CEE3",
                      "v06" = "#1F78B4",
                      "v07" = "#F0E442",
                      "v09" = "#E69F00")


vivax_vol_pal <- unlist(unname(vivax_colours))
names(vivax_vol_pal) <- names(vivax_colours)



volunteer_colours <- list("v02" = "#FB9A99",
                          "v03" = "#E31A1C",
                          "v05" = "#A6CEE3",
                          "v06" = "#1F78B4",
                          "v07" = "#F0E442",
                          "v09" = "#E69F00")


volunteer_palette <- unlist(unname(volunteer_colours))
names(volunteer_palette) <- names(volunteer_colours)



vivax_falci_palette <- c(vol_pal, vivax_vol_pal)





vivax_legend <- read.csv("~/PhD/plasma/vac69a/big_plasma_table.csv")
vivax_legend$timepoint <- gsub("DoD", "Diagnosis", vivax_legend$timepoint)

#vivax_legend[,3:ncol(vivax_legend)] <- log2(vivax_legend[,3:ncol(vivax_legend)])



vivax_legend <- vivax_legend %>%
  mutate(Volunteer = gsub("00", "0", Volunteer)) %>%
  mutate(Timepoint = gsub("C-1", "Baseline", timepoint, fixed=TRUE)) %>%
  mutate(Timepoint = gsub("+", "", Timepoint, fixed = T)) %>%
  mutate(Species = "P. vivax")

vivax_legend$Timepoint <- factor(vivax_legend$Timepoint, levels=c("Baseline", "Diagnosis", "T6", "C45"))
vivax_legend$timepoint <- NULL

vivax_legend <- pivot_longer(vivax_legend, cols = colnames(vivax_legend)[2:40], names_to = "Analyte", values_to = "Concentration")

#'vivax_legend$Sample_ID <- paste(vivax_legend$Volunteer, vivax_legend$timepoint, sep="_")





proto_falci_legend <- read.csv("~/PhD/plasma/vac63/Plasma_all_data_logdata_withVEGF_and_TGFB.csv", header=T)


bad_falci_timepoints <- unique(proto_falci_legend$timepoint_with_DoD)

good_falci_timepoints <- c("Baseline", "C28", "C6", "Diagnosis", "C7", "C8", "C45", "C90",
                           "T6", "C9", "C10", "C11", "C12", "C14", "C60", "T11",
                           "T18", "T8")

falci_timepoint_replacement <- setNames(good_falci_timepoints, bad_falci_timepoints)


proto_falci_legend$timepoint <- stringr::str_replace_all(proto_falci_legend$timepoint_with_DoD, falci_timepoint_replacement)

proto_falci_legend$timepoint <- gsub("D+", "T", proto_falci_legend$timepoint, fixed = TRUE)
proto_falci_legend$timepoint <- gsub("+", "", proto_falci_legend$timepoint, fixed = TRUE)


#proto_falci_legend <- subset(proto_falci_legend, proto_falci_legend$N_infection=="First")


falci_legend <- data.frame("Volunteer"=gsub("V", "v", proto_falci_legend$Volunteer_code, fixed=TRUE),
                           "Timepoint"=proto_falci_legend$timepoint,
                           "Species"="P. falciparum",
                           "Analyte"=proto_falci_legend$Analyte,
                           "Concentration"=proto_falci_legend$Concentration)

falci_legend <- cbind(falci_legend, proto_falci_legend[,16:56])

falci_legend$Analyte <- gsub(" (g/L)", "", falci_legend$Analyte, fixed=T)
falci_legend$Analyte <- gsub(" (pg/ml)", "", falci_legend$Analyte, fixed = TRUE)
falci_legend$Analyte <- gsub("(pg/ml)", "", falci_legend$Analyte, fixed = TRUE)

falci_legend$Analyte <- gsub("-", "", falci_legend$Analyte, fixed=T)
falci_legend$Analyte <- gsub(".", "", falci_legend$Analyte, fixed=T)
falci_legend$Analyte <- gsub("α", "a", falci_legend$Analyte, fixed=T)
falci_legend$Analyte <- gsub("β", "b", falci_legend$Analyte, fixed=T)
falci_legend$Analyte <- gsub("γ", "y", falci_legend$Analyte, fixed=T)
falci_legend$Analyte <- gsub("IL1b", "IL1B", falci_legend$Analyte, fixed=T)
falci_legend$Analyte <- gsub("(free active)", "", falci_legend$Analyte, fixed=T)
falci_legend$Analyte <- gsub("(IU/L)", "", falci_legend$Analyte, fixed=T)
falci_legend$Analyte <- gsub("IL12(p70)", "IL12p70", falci_legend$Analyte, fixed=T)


combo_data <- rbind(falci_legend, select(vivax_legend, colnames(falci_legend)))

#long_combo_data <- pivot_longer(combo_data, cols = colnames(combo_data)[4:42], names_to = "Analyte", values_to = "Concentration")


interesting_analytes <- scan("~/PhD/plasma/vac69a/analytes_sorted_by_padj.txt", what="", skip = 1)[c(1:3, 5:9)]

long_combo_data <- combo_data %>%
  filter(Timepoint %in% c("Baseline", "Diagnosis", "T6"))#%>%
  #filter(Analyte %in% interesting_analytes)

#long_combo_data$Analyte <- gsub("IFNy", "IFNγ", long_combo_data$Analyte, fixed=T)

long_combo_data <- filter(long_combo_data, Analyte %in% names(subset(list_of_adj_pvalues_dod, list_of_adj_pvalues_dod<=0.05)
))

vivax_falci_plasma_comparison <- ggplot(long_combo_data, aes(x=Timepoint, y=Concentration))+
  geom_boxplot(aes(fill=Species))+
  facet_wrap(~Analyte, scales="free", ncol = 5)+
  geom_point(aes(color=Volunteer, group=Species, shape=Species), position = position_jitterdodge(jitter.width = 0, dodge.width = 0.75))+
  scale_fill_manual(values=rev(c("#fec200", "#db0085")))+
  scale_y_log10()+
  ylab("plasma concentration (pg / mL)")+
  guides(color=guide_legend(title = "volunteer", direction = "horizontal", ncol = 11),
         fill=guide_legend(title = "species", label.theme = element_text(face="italic"), ncol=1),
         shape=guide_legend(title = "species", label.theme = element_text(face="italic")))+
  #scale_colour_manual(values=vivax_falci_palette)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.spacing.x = unit(1.5, "mm"),
        axis.text.x = element_text(angle = 45, hjust=1))


#ggsave("~/postdoc/BSI2021/vivax_falci_plasma_comparison_big_with_reinfected.png", vivax_falci_plasma_comparison, height=25, width=10, dpi=444, bg="white")
#ggsave("~/postdoc/BSI2021/vivax_falci_plasma_comparison_with_T6.png", vivax_falci_plasma_comparison, height=6, width=10, dpi=444, bg="white")
ggsave("~/postdoc/BSI2021/vivax_falci_all_sig_dod_analytes.png", vivax_falci_plasma_comparison, height=6, width=9, dpi=444, bg="white")

ggsave("~/postdoc/BSI2021/big_vivax_falci_plasma_comparison.png", vivax_falci_plasma_comparison, height=14, width=9, dpi=444, bg="white")




all_analytes <- scan("~/PhD/plasma/vac69a/analytes_sorted_by_padj.txt", what="", skip = 1)
combo_data$Analyte <- gsub("TGFb1freeactive", "TGFb1", combo_data$Analyte)

glm_data <- subset(combo_data, combo_data$Analyte %in% all_analytes)
glm_data <- subset(glm_data, glm_data$Timepoint %in% c("Baseline", "Diagnosis", "T6"))

glm_data$Concentration <- log10(glm_data$Concentration)

list_of_dfs_for_glm <- split(glm_data, glm_data$Analyte)

#list_of_models <- lapply(list_of_dfs_for_glm, function(x) lm(Concentration~timepoint+Volunteer, data=x))
#list_of_models <- lapply(list_of_dfs_for_glm, function(x) nlme::lme(concentration~timepoint, random=~1|Volunteer, data=x))
list_of_models <- lapply(list_of_dfs_for_glm, function(x) lme4::lmer(Concentration~Timepoint+Species+(1|Volunteer), data=x))

# 
dod_contrast <- t(matrix(c(0,1,0,1)))
t6_contrast <- t(matrix(c(0,0,1,1)))

list_of_tests_t6 <- lapply(list_of_models, function(x) multcomp::glht(x, t6_contrast))
list_of_pvalues_t6 <- sapply(list_of_tests_t6, function(x) summary(x)$test$pvalues)

#
#
list_of_adj_pvalues_t6 <- sort(p.adjust(list_of_pvalues_t6, method = "fdr"))
subset(list_of_adj_pvalues_t6, list_of_adj_pvalues_t6<=0.05)
subset(list_of_pvalues_t6, list_of_pvalues_t6<=0.05)



list_of_tests_dod <- lapply(list_of_models, function(x) multcomp::glht(x, dod_contrast))
list_of_pvalues_dod <- sapply(list_of_tests_dod, function(x) summary(x)$test$pvalues)

list_of_adj_pvalues_dod <- sort(p.adjust(list_of_pvalues_dod, method = "fdr"))
subset(list_of_adj_pvalues_dod, list_of_adj_pvalues_dod<=0.05)
subset(list_of_pvalues_dod, list_of_pvalues_dod<=0.05)


adj_results_df <- data.frame("dod"=list_of_adj_pvalues_dod, "t6"=list_of_adj_pvalues_t6)
non_adj_results_df <- data.frame("dod"=list_of_pvalues_dod, "t6"=list_of_pvalues_t6)
species_estimate_df <- sapply(list_of_models, function(x)summary(x)$coef[4,1])

write.csv(adj_results_df, "~/PhD/manuscripts/vac69a/jci_corrections/plasma_stats_adj_pvalue.csv")
write.csv(non_adj_results_df, "~/PhD/manuscripts/vac69a/jci_corrections/plasma_stats_NON_adj_pvalue.csv")



# haematology vivax vs. falci ####
# vivax


library(ggplot2)
library(tidyr)
library(dplyr)

`%notin%` <- Negate(`%in%`)

fig1_theme <- theme(axis.title.x = element_blank(),
                    legend.title = element_text(size = 9), 
                    legend.text = element_text(size = 9),
                    axis.title=element_text(size=10))


volunteer_colours <- list("v02" = "#FB9A99",
                          "v03" = "#E31A1C",
                          "v05" = "#A6CEE3",
                          "v06" = "#1F78B4",
                          "v07" = "#F0E442",
                          "v09" = "#E69F00")


volunteer_palette <- unlist(unname(volunteer_colours))
names(volunteer_palette) <- names(volunteer_colours)




haem_data <- data.table::fread("~/PhD/clinical_data/vac69a/haem.csv")

haem_data$volunteer <- gsub("69010", "v", haem_data$trial_number)

long_haem_data <- haem_data %>%
  select(volunteer, timepoint, platelets, lymphocytes) %>%
  #select(trial_number, timepoint, lymphocytes) %>%
  filter(timepoint %in% c("_C_1", "_C1_7", "_C8_14", "_C21", "_EP", "_T6", "_C90")) %>%
  gather(Cell, Frequency, c(platelets,lymphocytes))





supp_haem_data <- haem_data %>%
  select(volunteer, timepoint, wbc, neutrophils,  monocytes, eosinophils, haemaglobin) %>%
  #select(trial_number, timepoint, lymphocytes) %>%
  filter(timepoint %in% c("_C_1", "_C1_7", "_C8_14", "_C21", "_EP", "_T6", "_C90")) %>%
  gather(Cell, Frequency, c(wbc, neutrophils,  monocytes, eosinophils, haemaglobin))


supp_haem_data$Cell <- paste(
  toupper(substr(supp_haem_data$Cell, 1,1)),
  substr(supp_haem_data$Cell, 2,nchar(supp_haem_data$Cell)),
  sep="")

supp_haem_data$Cell <- factor(supp_haem_data$Cell, levels=c("Wbc", "Neutrophils",  "Monocytes", "Eosinophils", "Haemaglobin"))


bad_timepoints <- c("_C_1", "_C1_7", "_C8_14", "_C21", "_EP", "_T6", "_C90")
great_timepoints <- c("Baseline", "C7 am", "C14 am", "Diagnosis", "T1", "T6", "C90")

time_dic <- setNames(great_timepoints, bad_timepoints)

long_haem_data$timepoint <- stringr::str_replace_all(long_haem_data$timepoint, time_dic)
supp_haem_data$timepoint <- stringr::str_replace_all(supp_haem_data$timepoint, time_dic)


supp_haem_data$Cell <- gsub("Haemaglobin", "Haemoglobin", supp_haem_data$Cell)

haemoglobin_data <-  filter(supp_haem_data, Cell == "Haemoglobin")
supp_haem_data <- filter(supp_haem_data, Cell != "Haemoglobin")



lymphsss <- filter(long_haem_data, Cell=="lymphocytes")
# 
# thrombos_lymphs <- ggplot(long_haem_data, aes(x=factor(timepoint, levels=c("Baseline", "C7 am", "C14 am", "Diagnosis", "T1", "T6", "C90")), y=Frequency*1000, color=volunteer, group=volunteer))+
#   scale_fill_manual(values=volunteer_palette)+
#   scale_color_manual(values=volunteer_palette)+
#   geom_line(aes(color=volunteer), size=0.9)+
#   geom_point(fill="white", stroke=1, shape=21, size=0.9)+
#   theme_minimal()+
#   facet_wrap(~Cell, scales="free")+
#   xlab("Timepoint")+
#   ylab(expression(Cells~"/"~mu*L~blood))+
#   guides(color=guide_legend(title="Volunteer", override.aes = list(size=1)))+
#   fig1_theme+
#   scale_y_continuous(label=scales::comma)+
#   theme(plot.title = element_text(hjust=0.5),
#         axis.text.x = element_text(hjust=1, angle=45, size=8), 
#         axis.title.x = element_blank(),
#         legend.position = "none",
#         strip.text = element_text(size=10))
# 
# 
# ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/thrombos_lymphs.pdf", thrombos_lymphs, width=6, height=2.2)



lymph_plot <- ggplot(lymphsss, aes(x=factor(timepoint, levels=c("Baseline", "C7 am", "C14 am", "Diagnosis", "T1", "T6", "C90")), y=Frequency*1000, color=volunteer, group=volunteer))+
  scale_fill_manual(values=volunteer_palette)+
  scale_color_manual(values=volunteer_palette)+
  geom_line(aes(color=volunteer), size=0.9)+
  geom_point(fill="white", stroke=1, shape=21, size=0.9)+
  theme_minimal()+
  facet_wrap(~Cell, scales="free")+
  xlab("Timepoint")+
  ylab(expression(Cells~"/"~mu*L~blood))+
  guides(color=guide_legend(title="Volunteer", override.aes = list(size=1)))+
  fig1_theme+
  scale_y_continuous(label=scales::comma)+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8), 
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.text = element_text(size=10))


#ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/lymphs.pdf", lymph_plot, width=3, height=2.2)

# falci #


#vol_pal <- colorspace::qualitative_hcl(n=9, palette = "dark3")[c(1,9,8,2:7)]
vol_pal <-c("#FFC800",
            "#EE000C",
            "#FF8000",
            "#0080FF",
            "#009B95",
            "#E54787",
            "#8000FF",
            "#66CCFF",
            "#4B0055")

names(vol_pal) <- c("v313", "v315", "v320", "v306", "v301", "v308", "v305", "v304", "v310")



para_vol_pal <- vol_pal
names(para_vol_pal ) <- c("313", "315", "320", "1039", "1040", "1061", "1068", "1075", "6032")

lymph_vol_pal <- para_vol_pal
names(lymph_vol_pal) <- paste("v", names(para_vol_pal), sep="")


time_col <- colorspace::sequential_hcl(5, palette = "Purple Yellow")



vac63c_lymph <- read.csv("~/PhD/clinical_data/vac63c/VAC063_haem_all_sequenced_WBC_real_percent.csv")

vac63c_lymph <- dplyr::filter(vac63c_lymph, Leukocytes=="Lymphocytes")

vac63c_lymph <- select(vac63c_lymph, Volunteer_code, trial_number, timepoint, N_infection, Leukocytes, cell_counts)
vac63c_lymph <- dplyr::filter(vac63c_lymph, timepoint %in% c("C-1", "Diagnosis", "D+6"))
#vac63c_lymph <- dplyr::filter(vac63c_lymph, N_infection %in% c("First", "Diagnosis", "Third"))
vac63c_lymph$timepoint <- gsub("D+6", "T6", vac63c_lymph$timepoint, fixed=T)
vac63c_lymph$timepoint <- gsub("C-1", "Baseline", vac63c_lymph$timepoint, fixed=T)

vac63c_lymph$Volunteer_code <- gsub("V", "v", vac63c_lymph$Volunteer_code, fixed=T)
#vac63c_lymph<- filter(vac63c_lymph, vac63c_lymph$Volunteer_code %in% names(lymph_vol_pal))

n_infection_palette <- c(rgb(5,50,80, maxColorValue = 255),
                         rgb(250, 100, 0, maxColorValue = 255),
                         rgb(40,210,250, maxColorValue = 255))


group_lymph_vac63c_box <- ggplot(vac63c_lymph, aes(x=factor(timepoint, levels=c("Baseline", "Diagnosis", "T6")), y=cell_counts*1000, fill=N_infection))+
  geom_boxplot()+
  theme_minimal()+
  ylab(expression(Lymphocytes~"/"~mu*L~blood))+
  labs(fill="Infection")+
  scale_y_continuous(label=scales::comma)+
  #scale_color_manual(values=vol_pal)+
  scale_fill_manual(values=n_infection_palette)+
  guides(colour=guide_legend(title="Volunteer", override.aes = list(alpha=1)))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=12, angle=45, hjust=1),
        plot.margin = unit(c(2,2,2,2), "mm"))


ggsave("~/postdoc/BSI2021/group_lymph_vac63c_box.png", group_lymph_vac63c_box, width=5, height=5, bg="white", dpi=444)
ggsave("~/postdoc/BSI2021/group_lymph_vac63c_box.pdf", group_lymph_vac63c_box, width=5, height=5, bg="white")




# lgd <- get_legend(group_lymph_vac63c_box)
# group_lymph_vac63c_box <- group_lymph_vac63c_box + theme(legend.position = "none")

# make common x axis title
timepoint.grob <- textGrob("Timepoint", 
                           gp=gpar(fontsize=14))

#add to plot
vac63c_lymphocytes_figure <- plot_grid(indie_lymph_vac63c, group_lymph_vac63c_box, align="h", axis="tblr")
#vac63c_lymphocytes_figure <- grid.arrange(arrangeGrob(vac63c_lymphocytes_figure, bottom = timepoint.grob))

lymphsss <- filter(lymphsss, timepoint %in% c("Baseline", "Diagnosis", "T6"))

vac63c_lymph <- filter(vac63c_lymph, N_infection == "First")

vivax_falci_lymphs <- data.frame("volunteer"=c(vac63c_lymph$Volunteer_code, lymphsss$volunteer),
                                 "timepoint"=c(vac63c_lymph$timepoint, lymphsss$timepoint),
                                 "lymphocytes"=c(vac63c_lymph$cell_counts, lymphsss$Frequency),
                                 "species"=c(rep("P. falciparum", nrow(vac63c_lymph)), rep("P. vivax", nrow(lymphsss))))


vivax_falci_lymphocytes <- ggplot(vivax_falci_lymphs, aes(x=factor(timepoint, levels=c("Baseline", "Diagnosis", "T6")), y=lymphocytes*1000, ))+
  geom_boxplot(aes(fill=species))+
  #geom_line(aes(group=interaction(species, volunteer)), alpha=0.1)+
  #geom_point(aes(colour=volunteer, group=species), position = position_dodge(width = 0.75))+
  theme_minimal()+
  xlab("Timepoint")+
  scale_y_continuous(label=scales::comma)+
  ylab(expression(Lymphocytes~"/"~mu*L~blood))+
  scale_fill_manual(values=rev(c("#fec200", "#db0085")))+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"),
        axis.text.x = element_text(size=12, angle=45, hjust=1),
        plot.margin = unit(c(2,2,2,2), "mm"))

ggsave("~/postdoc/BSI2021/vivax_falci_lymphocytes.png", vivax_falci_lymphocytes, width=5, height=5, bg="white", dpi=444)
ggsave("~/postdoc/BSI2021/vivax_falci_lymphocytes.pdf", vivax_falci_lymphocytes, width=5, height=5, bg="white")





interesting_analytes <- scan("~/PhD/plasma/vac69a/analytes_sorted_by_padj.txt", what="", skip = 1)[1:16]

pc_data <- combo_data %>%
  filter(Timepoint %in% c("Baseline", "Diagnosis"))%>%
  filter(Analyte %in% interesting_analytes)

wide_pc_data <- pivot_wider(pc_data, names_from = Analyte, values_from = Concentration)

big_pca <-  prcomp(log10(wide_pc_data[,4:ncol(wide_pc_data)]), center = T)
big_pca2 <- cbind(wide_pc_data[, 1:3], big_pca$x)

time_col <- colorspace::sequential_hcl(5, palette = "Purple Yellow")
time_col[5] <- "peachpuff1"
Volunteer <- c("v02" = "#FB9A99","v03" = "#E31A1C","v05" = "#A6CEE3", "v06" = "#1F78B4", "v07" = "#F0E442")
Timepoint <- c("Baseline"=time_col[4], "Diagnosis"=time_col[2], "T6"=time_col[1], "C45"=time_col[5])


PCA_timepoint_col <- ggplot(big_pca2, aes(x=PC1, y=PC2))+
  geom_text(aes(label=Species, color=Timepoint), size=2.5, fontface="bold")+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values = Timepoint[1:2])+
  theme(legend.position = "right",
        axis.text = element_text(size=10),
        panel.border = element_rect(color="black", fill=NA),
        axis.title = element_text(size=12),
        #axis.text = element_blank(),
        
        plot.title = element_text(size=14, hjust=0.5)
  )


PCA_species_col <- ggplot(big_pca2, aes(x=PC1, y=PC2))+
  geom_text(aes(label=Timepoint, color=Species), size=2.5, fontface="bold")+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=rev(c("#fec200", "#db0085")))+
  theme(legend.position = "right",
        axis.text = element_text(size=10),
        panel.border = element_rect(color="black", fill=NA),
        axis.title = element_text(size=12),
        #axis.text = element_blank(),
        
        plot.title = element_text(size=14, hjust=0.5)
  )


cowplot::plot_grid(PCA_species_col, PCA_timepoint_col)




ggsave("~/PhD/figures_for_thesis/chapter_1/1_sig_plasma_pca.png", all_vols_together_vol_color, height=3, width=4.5)
