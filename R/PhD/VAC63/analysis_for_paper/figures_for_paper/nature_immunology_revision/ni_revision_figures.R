library(dplyr)
library(tidyr)
library(ggplot2)

n_infection_palette <- c(rgb(5,50,80, maxColorValue = 255),
                         rgb(250, 100, 0, maxColorValue = 255),
                         rgb(40,210,250, maxColorValue = 255))

# T6 data only, missing vivax reinfection
combo_alt_lineage_acti_data <- read.csv("~/PhD/manuscripts/vac69a/jci_corrections/vac69a_b_vac63c_lineage_acti_data_.csv")

combo_alt_lineage_acti_data$n_infection <- ifelse(combo_alt_lineage_acti_data$volunteer %in% c("v301", "v302", "v304", "v305", "v306", "v307", "v308", "v310"), "Third", "First")

#T6 data only, WITH vivax reinfection, but no falciparum data
vac69b_data <-  read.csv("~/PhD/manuscripts/vac69a/jci_corrections/vac69ab_t_lineage_acti_alt.csv")

combo_frame <- rbind(combo_alt_lineage_acti_data, vac69b_data)
combo_frame <- combo_frame[!duplicated(combo_frame),]

combo_frame$species <- ifelse(nchar(combo_frame$volunteer)==3, "P. vivax", "P.falciparum")
combo_frame$first <- ifelse(combo_frame$n_infection=="First", "First", "Reinfection")

combo_frame <- subset(combo_frame, combo_frame$lineage!="NKT")
#
# 
# combo_alt_lineage_acti_data <- combo_alt_lineage_acti_data%>%
#   filter(volunteer %notin% c("v301", "v302", "v304", "v305", "v306", "v307", "v308", "v310"), lineage %in% c("CD4", "Treg"))
# 
# #combo_alt_lineage_acti_data <- subset(combo_alt_lineage_acti_data, combo_alt_lineage_acti_data$volunteer %notin% c("v313", "v315", "v320"))
# 
pearson_results <- combo_frame %>%
   group_by(lineage) %>%
   do(broom::tidy(cor.test(.$activated, .$alt, method="pearson")))

spearman_results <- combo_frame %>%
  group_by(lineage) %>%
  do(broom::tidy(cor.test(.$activated, .$alt, method="spearman")))

write.table(data.frame(rbind(pearson_results, spearman_results)), "~/postdoc/grant_stuff/vivax_falci_all_tcell_alt_correlations.csv", row.names = FALSE, quote = FALSE, sep=",")

# 
# lineage estimate statistic p.value parameter conf.low conf.high method          alternative
# <chr>      <dbl>     <dbl>   <dbl>     <int>    <dbl>     <dbl> <chr>           <chr>      
#   1 CD4       0.572      2.95  0.00848        18    0.173     0.809 Pearson's prod… two.sided  
# 2 CD8       0.535      2.68  0.0152         18    0.121     0.790 Pearson's prod… two.sided  
# 3 DN        0.0376     0.160 0.875          18   -0.412     0.472 Pearson's prod… two.sided  
# 4 gd        0.627      3.41  0.00311        18    0.255     0.837 Pearson's prod… two.sided  
# 5 MAIT      0.677      3.90  0.00104        18    0.335     0.861 Pearson's prod… two.sided  
# 6 Treg      0.586      3.06  0.00668        18    0.193     0.816 Pearson's prod… two.sided 
# 
# combo_alt_lineage_acti_data$species <- ifelse(nchar(combo_alt_lineage_acti_data$volunteer)==3, "P. vivax", "P.falciparum")
# 
# 
# combo_alt_lineage_acti_data$lineage <- gsub("CD4", "CD4 (r=0.636, p=0.0404)", combo_alt_lineage_acti_data$lineage)
# combo_alt_lineage_acti_data$lineage <- gsub("Treg", "Treg (r=0.618, p=0.0478)", combo_alt_lineage_acti_data$lineage)

(all_lineages_alt_plot <- ggplot(combo_frame, aes(x=alt, y=activated/100))+
  geom_point(aes(colour=first, shape=species))+
  geom_smooth(method="lm", colour="black")+
  ylab("Fraction of Lineage Activated")+
  xlab("ALT (IU / L)")+
  # guides(shape = guide_legend(title="Species"),
  #        color = guide_legend(title="Volunteer",
  #                             override.aes = list(shape= c(rep(16,6),
  #                                                          rep(17,3)))),
  # )+
  facet_wrap(~lineage, scales = "free")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 2))+
  scale_colour_manual(values=n_infection_palette[c(1,3)])+
  theme_minimal()+
  theme(strip.text=element_text(size=12),
        legend.title = element_blank()))

ggsave("~/PhD/manuscripts/vac63c/nature_immunology/pdf/alt_tcell_activation_correlation_vac69ab_vac63c.pdf", all_lineages_alt_plot, height=4, width=6.5, bg="white")
ggsave("~/PhD/manuscripts/vac63c/nature_immunology/png/alt_tcell_activation_correlation_vac69ab_vac63c.png", all_lineages_alt_plot, height=4, width=6.5, dpi=444, bg="white")


# more granular stuff ####

# vivax data

all_data <- read.table("~/PhD/manuscripts/vac63c/nature_immunology/vac69a_b_vac63c_cluster_frequencies.csv")
all_cd4 <- subset(all_data, lineage=="CD4")
all_cd4[grep("cytotoxic*", all_cd4$cluster_id),]
