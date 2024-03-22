library(broom)
library(ggplot2)
library(tidyr)
library(dplyr)


positive_hits <- read.csv("~/postdoc/stanford/misc/basset_peptides/positive_hits.csv")
shortnames <- read.csv("~/postdoc/stanford/misc/basset_peptides/shortnames_to_HITS_meanZscore_Round2_Pfonly_250kfil_techclean_3zscorefil_5patientsfil (2).csv")
ugandan_metadata <- read.csv("~/postdoc/stanford/misc/basset_peptides/Ugandan_samples_metadata.csv")

# time since last infection, Incidence of symptomatic malaria per year, Household annual EIR (infective bites / person),
# age, positive at the time of infection

hits_of_interest <- c("t127813",
                      "t127814",
                      "t127815",
                      "t121177",
                      "t121180",
                      "t123838",
                      "t91292")

zscores_of_interest <- filter(shortnames, shortname %in% hits_of_interest)

long_zscores_of_interest <- zscores_of_interest%>%
  pivot_longer(cols = colnames(zscores_of_interest)[-c(1, 288)], values_to = "z_score", names_to = "id")%>%
  pivot_wider(names_from = shortname, values_from = z_score, id_cols = id)%>%
  mutate(Barcode=gsub(".", "-", id, fixed = TRUE))

merged_data <- inner_join(long_zscores_of_interest, ugandan_metadata, by = "Barcode")

long_merge <- merged_data %>%
  pivot_longer(cols=all_of(hits_of_interest), names_to = "peptide", values_to = "z_score")

last_inf_plot <- ggplot(long_merge, aes(x=time.since.last.infection, y=z_score))+
  geom_point(aes(color=Barcode), alpha=0.5)+
  geom_smooth(method="lm")+
  facet_wrap(~peptide, scales = "free")+
  theme_minimal()+
  xlab("Time Since Last Infection")+
  theme(legend.position="none")
  
ggsave("~/postdoc/stanford/misc/basset_peptides/figures/time.since.last.infection.png", last_inf_plot, height=6, width=8, bg="white", dpi=444)
  

incidence_plot <- ggplot(long_merge, aes(x=Incidence.of.symptomatic.malaria.per.year, y=z_score))+
    geom_point(aes(color=Barcode), alpha=0.5)+
    geom_smooth(method="lm")+
    facet_wrap(~peptide, scales = "free")+
    theme_minimal()+
    xlab("Incidence of Symptomatic Malaria per Year")+
    theme(legend.position="none")

ggsave("~/postdoc/stanford/misc/basset_peptides/figures/incidence_plot.png", incidence_plot, height=6, width=8, bg="white", dpi=444)

  
eir_plot <- ggplot(long_merge, aes(x=Household.annual.EIR..infective.bites...person., y=z_score))+
    geom_point(aes(color=Barcode), alpha=0.5)+
    geom_smooth(method="lm")+
    facet_wrap(~peptide, scales = "free")+
    theme_minimal()+
    xlab("Household annual EIR")+
    theme(legend.position="none")
  
ggsave("~/postdoc/stanford/misc/basset_peptides/figures/eir_plot.png", eir_plot, height=6, width=8, bg="white", dpi=444)

age_plot <- long_merge %>%
  # filter(age<10)%>%
  ggplot(., aes(x=age, y=z_score))+
    geom_point(aes(color=Barcode), alpha=0.5)+
    # scale_y_binned(breaks = c(seq(0, 70, by=10)))+
    geom_smooth(method="lm")+
    facet_wrap(~peptide, scales = "free")+
    theme_minimal()+
    xlab("Age")+
    theme(legend.position="none")
  
ggsave("~/postdoc/stanford/misc/basset_peptides/figures/age_plot.png", age_plot, height=6, width=8, bg="white", dpi=444)

pos_plot <- ggplot(long_merge, aes(x=positive.at.the.time.of.infection, y=z_score))+
    geom_point(aes(color=Barcode), alpha=0.5)+
    geom_boxplot()+
    geom_smooth(method="lm")+
    facet_wrap(~peptide, scales = "free")+
    theme_minimal()+
    xlab("Positive at The Time of Infection")+
    theme(legend.position="none")

ggsave("~/postdoc/stanford/misc/basset_peptides/figures/positive_plot.png", pos_plot, height=6, width=8, bg="white", dpi=444)



spearman_last_infection <- long_merge %>%
  group_by(peptide) %>%
  do(broom::tidy(cor.test(.$z_score, .$time.since.last.infection, method="spearman")))%>%
  write.csv(., "~/postdoc/stanford/misc/basset_peptides/spearman_last_infection.csv")

spearman_incidence <- long_merge %>%
  group_by(peptide) %>%
  do(broom::tidy(cor.test(.$z_score, .$Incidence.of.symptomatic.malaria.per.year, method="spearman")))%>%
  write.csv(., "~/postdoc/stanford/misc/basset_peptides/spearman_incidence.csv")


spearman_eir <- long_merge %>%
  group_by(peptide) %>%
  do(broom::tidy(cor.test(.$z_score, .$Household.annual.EIR..infective.bites...person., method="spearman")))%>%
  write.csv(., "~/postdoc/stanford/misc/basset_peptides/spearman_eir.csv")

spearman_age <- long_merge %>%
  group_by(peptide) %>%
  do(broom::tidy(cor.test(.$z_score, .$age, method="spearman")))%>%
  write.csv(., "~/postdoc/stanford/misc/basset_peptides/spearman_age.csv")


spearman_positive <- long_merge %>%
  group_by(peptide) %>%
  do(broom::tidy(cor.test(.$z_score, .$positive.at.the.time.of.infection, method="spearman")))%>%
  write.csv(., "~/postdoc/stanford/misc/basset_peptides/spearman_positive.csv")






# uncomment to run pearson correlations
# 
# pearson_last_infection <- long_merge %>%
#   group_by(peptide) %>%
#   do(broom::tidy(cor.test(.$z_score, .$time.since.last.infection, method="pearson")))%>%
#   write.csv(., "~/postdoc/stanford/misc/basset_peptides/pearson_last_infection.csv")
# 
# pearson_incidence <- long_merge %>%
#   group_by(peptide) %>%
#   do(broom::tidy(cor.test(.$z_score, .$Incidence.of.symptomatic.malaria.per.year, method="pearson")))%>%
#   write.csv(., "~/postdoc/stanford/misc/basset_peptides/pearson_incidence.csv")
# 
# 
# pearson_eir <- long_merge %>%
#   group_by(peptide) %>%
#   do(broom::tidy(cor.test(.$z_score, .$Household.annual.EIR..infective.bites...person., method="pearson")))%>%
#   write.csv(., "~/postdoc/stanford/misc/basset_peptides/pearson_eir.csv")
# 
# pearson_age <- long_merge %>%
#   group_by(peptide) %>%
#   do(broom::tidy(cor.test(.$z_score, .$age, method="pearson")))%>%
#   write.csv(., "~/postdoc/stanford/misc/basset_peptides/pearson_age.csv")

