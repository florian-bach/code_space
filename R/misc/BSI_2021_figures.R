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

vivax_legend[,3:ncol(vivax_legend)] <- log2(vivax_legend[,3:ncol(vivax_legend)])



vivax_legend <- vivax_legend %>%
  mutate(Volunteer = gsub("00", "0", Volunteer)) %>%
  mutate(timepoint = gsub("C-1", "Baseline", timepoint)) %>%
  mutate(timepoint = gsub("+", "", timepoint, fixed = T)) %>%
  mutate(Species = "P. vivax")

vivax_legend$timepoint <- factor(vivax_legend$timepoint, levels=c("Baseline", "Diagnosis", "T6", "C45"))
vivax_legend <- arrange(vivax_legend, timepoint, Volunteer)
#'vivax_legend$Sample_ID <- paste(vivax_legend$Volunteer, vivax_legend$timepoint, sep="_")


proto_falci_legend <- read.csv("~/PhD/plasma/vac63/Paired_data_log2transf_proteome.csv", header=T)

bad_falci_timepoints <- unique(proto_falci_legend$timepoint_with_DoD)
good_falci_timepoints <- c("Baseline", "C6", "C7", "C8", "Diagnosis", "C45", "T6", "C9", "C10", "C28", "C11", "C12", "C14")

falci_timepoint_replacement <- setNames(good_falci_timepoints, bad_falci_timepoints)

proto_falci_legend$timepoint_with_DoD <- stringr::str_replace_all(proto_falci_legend$timepoint_with_DoD, falci_timepoint_replacement)

proto_falci_legend <- subset(proto_falci_legend, proto_falci_legend$N_infection=="First")

falci_legend <- data.frame("Volunteer"=gsub("V", "v", proto_falci_legend$Volunteer_code, fixed=TRUE),
                           "timepoint"=proto_falci_legend$timepoint_with_DoD,
                           "Species"="P. falciparum")

falci_legend <- cbind(falci_legend, 2^proto_falci_legend[,6:44])

colnames(falci_legend) <- gsub(".", "", colnames(falci_legend), fixed=T)
colnames(falci_legend) <- gsub("α", "a", colnames(falci_legend), fixed=T)
colnames(falci_legend) <- gsub("β", "b", colnames(falci_legend), fixed=T)
colnames(falci_legend) <- gsub("γ", "y", colnames(falci_legend), fixed=T)
colnames(falci_legend) <- gsub("IL1b", "IL1B", colnames(falci_legend), fixed=T)

combo_data <- rbind(falci_legend, vivax_legend)

long_combo_data <- pivot_longer(combo_data, cols = colnames(combo_data)[4:42], names_to = "Analyte", values_to = "Concentration")


interesting_analytes <- scan("~/PhD/plasma/vac69a/analytes_sorted_by_padj.txt", what="", skip = 1)[1:12]

long_combo_data <- long_combo_data %>%
  filter(timepoint %in% c("Baseline", "Diagnosis", "T6"))%>%
  filter(Analyte %in% interesting_analytes)

vivax_falci_plasma_comparison <- ggplot(long_combo_data, aes(x=timepoint, y=2^Concentration))+
  geom_boxplot(aes(fill=Species))+
  facet_wrap(~Analyte, scales="free")+
  geom_point(aes(color=Volunteer, group=Species, shape=Species), position = position_jitterdodge(jitter.width = 0, dodge.width = 0.75))+
  scale_fill_manual(values=rev(c("#fec200", "#db0085")))+
  scale_y_log10()+
  #scale_colour_manual(values=vivax_falci_palette)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1))


ggsave("~/postdoc/BSI2021/vivax_falci_plasma_comparison.png", vivax_falci_plasma_comparison, height=7, width=9, dpi=444, bg="white")




list_of_dfs_for_glm <- split(long_combo_data, long_combo_data$Analyte)

list_of_models <- lapply(list_of_dfs_for_glm, function(x) lm(Concentration~timepoint+Volunteer, data=x))
#list_of_models <- lapply(list_of_dfs_for_glm, function(x) nlme::lme(concentration~timepoint, random=~1|Volunteer, data=x))
#list_of_models <- lapply(list_of_dfs_for_glm, function(x) lme4::lmer(concentration~timepoint+(1|Volunteer), data=x))

# 
dod_contrast <- t(matrix(c(0,0,1,0)))
t6_contrast <- t(matrix(c(0,0,0,1)))

list_of_tests <- lapply(list_of_models, function(x) multcomp::glht(x, dod_contrast))
list_of_pvalues <- sapply(list_of_tests, function(x) summary(x)$test$pvalues)
#
#
list_of_adj_pvalues <- sort(p.adjust(list_of_pvalues, method = "fdr"))
subset(list_of_adj_pvalues, list_of_adj_pvalues<=0.05)

#
