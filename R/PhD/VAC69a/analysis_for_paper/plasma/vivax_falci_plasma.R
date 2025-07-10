
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


proto_falci_legend <- subset(proto_falci_legend, proto_falci_legend$N_infection=="First")


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


#interesting_analytes <- scan("~/PhD/plasma/vac69a/analytes_sorted_by_padj.txt", what="", skip = 1)[c(1:3, 5:9)]
interesting_analytes <- scan("~/PhD/plasma/vac69a/analytes_sorted_by_padj.txt", what="", skip = 1)[1:12]

long_combo_data <- combo_data %>%
  filter(Timepoint %in% c("Baseline", "Diagnosis", "T6"))%>%
  filter(Analyte %in% interesting_analytes)

long_combo_data$Analyte <- gsub("IFNy", "IFNγ", long_combo_data$Analyte, fixed=T)

# long_combo_data <- filter(long_combo_data, Analyte %in% names(subset(list_of_adj_pvalues_dod, list_of_adj_pvalues_dod<=0.05)
# ))

long_combo_data <- filter(long_combo_data, n_infection == "First")

vivax_falci_plasma_comparison <- ggplot(long_combo_data, aes(x=Timepoint, y=Concentration))+
  geom_boxplot(aes(fill=Species), outlier.alpha = 1)+
  facet_wrap(~Analyte, scales="free", ncol = 4)+
  #geom_point(aes(color=Volunteer, group=Species, shape=Species), position = position_jitterdodge(jitter.width = 0, dodge.width = 0.75))+
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
        legend.title = element_blank(),
        legend.spacing.x = unit(1.5, "mm"),
        axis.text.x = element_text(angle = 45, hjust=1))


#ggsave("~/postdoc/BSI2021/vivax_falci_plasma_comparison_big_with_reinfected.png", vivax_falci_plasma_comparison, height=25, width=10, dpi=444, bg="white")
#ggsave("~/postdoc/BSI2021/vivax_falci_plasma_comparison_with_T6.png", vivax_falci_plasma_comparison, height=6, width=10, dpi=444, bg="white")
ggsave("~/postdoc/presentations/BSI2021/vivax_falci_all_sig_dod_analytes.png", vivax_falci_plasma_comparison, height=6, width=9, dpi=444, bg="white")
ggsave("~/postdoc/presentations/BSI2021/vivax_falci_all_sig_dod_analytes.pdf", vivax_falci_plasma_comparison, height=6, width=9, bg="white", device = cairo_pdf)
ggsave("~/postdoc/presentations/BSI2021/vivax_falci_all_sig_dod_analytes_tall.pdf", vivax_falci_plasma_comparison, height=9, width=9, bg="white", device = cairo_pdf)

ggsave("~/postdoc/presentations/BSI2021/first_second_third_vivax_falci_plasma_comparison.png", vivax_falci_plasma_comparison, height=17, width=11.5, dpi=444, bg="white")
ggsave("~/postdoc/presentations/BSI2021/first_vivax_falci_plasma_comparison.png", vivax_falci_plasma_comparison, height=17, width=11.5, dpi=444, bg="white")




all_analytes <- scan("~/PhD/plasma/vac69a/analytes_sorted_by_padj.txt", what="", skip = 1)
combo_data$Analyte <- gsub("TGFb1freeactive", "TGFb1", combo_data$Analyte)

glm_data <- subset(combo_data, combo_data$Analyte %in% all_analytes)
glm_data <- subset(glm_data, glm_data$Timepoint %in% c("Baseline", "Diagnosis", "T6"))

glm_data$Concentration <- log10(glm_data$Concentration)

list_of_dfs_for_glm <- split(glm_data, glm_data$Analyte)

#list_of_models <- lapply(list_of_dfs_for_glm, function(x) lm(Concentration~timepoint+Volunteer, data=x))
#list_of_models <- lapply(list_of_dfs_for_glm, function(x) nlme::lme(Concentration~Timepoint+Species, random=~1|Volunteer, data=x))
list_of_models <- lapply(list_of_dfs_for_glm, function(x) lme4::lmer(Concentration~Timepoint+Species+(1|Volunteer), REML=TRUE, data=x))

# 
dod_contrast <- t(matrix(c(0,1,0,1)))
t6_contrast <- t(matrix(c(0,0,1,1)))
# 
list_of_tests_t6 <- lapply(list_of_models, function(x) multcomp::glht(model = x, linfct = t6_contrast))
t6_results <- lapply(list_of_tests_t6, function(x) cbind("coef"=summary(x)$test$coefficients, 
                                                         "pvalue"=summary(x)$test$pvalues))

# 
# #
# #
# list_of_adj_pvalues_t6 <- sort(p.adjust(list_of_pvalues_t6, method = "fdr"))
# subset(list_of_adj_pvalues_t6, list_of_adj_pvalues_t6<=0.05)
# subset(list_of_pvalues_t6, list_of_pvalues_t6<=0.05)
# 


list_of_tests_dod <- lapply(list_of_models, function(x) multcomp::glht(x, dod_contrast))
dod_results <- lapply(list_of_tests_dod, function(x) cbind("coef"=summary(x)$test$coefficients, 
                                                           "pvalue"=summary(x)$test$pvalues))

dod_results_df <- data.frame("analyte"=names(dod_results),
                             "contrast"="base_dod_species",
                             do.call(rbind, dod_results))

#coef>1 means more in vivax, coef<1 means more in falciparum
dod_results_df$padj <- p.adjust(dod_results_df$pvalue)
subset(dod_results_df, dod_results_df$padj<=0.05)





t6_results_df <- data.frame("analyte"=names(t6_results),
                            "contrast"="base_t6_species",
                            do.call(rbind, t6_results))

#coef>1 means more in vivax, coef<1 means more in falciparum
t6_results_df$padj <- p.adjust(t6_results_df$pvalue)
subset(t6_results_df, t6_results_df$padj<=0.05)


adj_results_df <- data.frame("dod"=list_of_adj_pvalues_dod, "t6"=list_of_adj_pvalues_t6)
non_adj_results_df <- data.frame("dod"=list_of_pvalues_dod, "t6"=list_of_pvalues_t6)
species_estimate_df <- sapply(list_of_models, function(x)summary(x)$coef[4,1])

write.csv(adj_results_df, "~/PhD/manuscripts/vac69a/jci_corrections/plasma_stats_adj_pvalue.csv")
write.csv(non_adj_results_df, "~/PhD/manuscripts/vac69a/jci_corrections/plasma_stats_NON_adj_pvalue.csv")


