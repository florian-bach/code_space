library(magrittr)
library(ggplot2)
library(ggrepel)


# DoD sheet
dod_data <- data.table::fread("~/PhD/RNAseq/vac69a/cytoscape/vivax_falciparum_dod_all.csv", header = T, stringsAsFactors = F)
dod_data <- subset(dod_data, grepl("5", dod_data$GOLevels))




# positive means enriched in vivax, negative means enriched in falciparum
dod_data$Cluster_Difference <- dod_data$`%Genes Cluster #2`-dod_data$`%Genes Cluster #1`



threshold <- 30

#DoD
dod_falci_rich_data <- subset(dod_data, dod_data$Cluster_Difference< -threshold)
dod_vivax_rich_data <- subset(dod_data, dod_data$Cluster_Difference> threshold)

vivax_x_limits <- c(80, NA)
vivax_y_limits <- c(20, NA)

dod_dot_plot <- ggplot(dod_data, aes(x=`%Genes Cluster #2`, y=`%Genes Cluster #1`))+
  theme_minimal()+
  scale_color_gradient2(high="#fec200", low="#db0085", midpoint = 0)+
  scale_y_continuous(breaks=seq(0,100,by=10), labels = scales::label_number(suffix="%"), limits = c(-5, 105), expand=c(0,0))+
  scale_x_continuous(breaks=seq(0,100,by=10), labels = scales::label_number(suffix="%"), limits = c(-5, 105), expand=c(0,0))+
  ylab(expression('GO term enrichement'~italic("P. falciparum")~'at Diagnosis'))+
  xlab(expression('GO term enrichement'~italic("P. vivax")~'at Diagnosis'))+
  ggitle("Diagnosis")+
  ggforce::geom_circle(aes(x0=50, y0=50, r=sqrt((0.5*threshold)^2+(0.5*threshold)^2)), fill="grey", color=NA, alpha=0.2, inherit.aes = F)+
  # geom_rect(aes(ymin=50-(threshold*0.5), ymax=50+(threshold*0.5),
  #               xmin=50-(threshold*0.5), xmax=50+(threshold*0.5)),fill="grey", color=NA, alpha=0.2, inherit.aes = F )+
  geom_point(aes(color=Cluster_Difference))+
  # geom_label_repel(data=dod_vivax_rich_data, aes(label = stringr::str_wrap(GOTerm, 25)), size=1.8,
  #                  box.padding   = 0.35, nudge_x=-20, nudge_y=-30,
  #                  point.padding = 0.5, segment.alpha = 0.2, ylim  = vivax_y_limits, xlim  = vivax_x_limits)+
  theme(legend.position = "none",
        axis.text = element_text(size=6, ),
        axis.title = element_text(size=8),
        plot.margin=unit(c(0.5,0,0.5,0.5),"cm"))

dod_hist_plot <- ggplot(dod_data, aes(y=`%Genes Cluster #1`))+
  theme_minimal()+
  xlab("# GO Terms")+
  scale_y_continuous(breaks = seq(0,100, by=20), limits=c(-5, 105))+
  scale_x_continuous(breaks = seq(10,60, by=20))+
  geom_histogram(aes(fill = ..y..), orientation = "y", binwidth = 1, color="darkgrey", size=0.15)+
  scale_fill_gradient2(low="#fec200", high="#db0085", midpoint = 50)+
  theme(axis.title.y =  element_blank(),
        axis.text.y =  element_blank(),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=6),
        legend.position = "none",
        plot.margin=unit(c(0.5,0.5,0.5,0.3),"cm"))


vivax_falciparum_dod <- cowplot::plot_grid(dod_dot_plot, dod_hist_plot, ncol=2, rel_widths = c(3,1), align="h", axis="b")
ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/vivax_falciparum_dod_GO.png", vivax_falciparum_dod, height = 3, width=3.9)





#### How to name Nodes ####

library(magrittr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readxl)


dod_data <- read_xls("~/PhD/RNAseq/vac69a/cytoscape/VIVAX_DOD_ALL/VIVAX_DOD_ALL_RESULTS_TABLE.xls")
t6_data <-  read_xls("~/PhD/RNAseq/vac69a/cytoscape/VIVAX_T6_ALL/VIVAX_T6_ALL_RESULTS_TABLE")

base_dod_t6_data <- read_xls("~/PhD/RNAseq/vac69a/cytoscape/VIVAX_BASE_DOD_T6/VIVAX_BASE_DOD_T6_Results_table.xls")

vivax_falciparum_dod <- read_xls("~/PhD/RNAseq/vac69a/cytoscape/VIVAX_FALCIPARUM_DOD/vivax_falciparum_dod_results_table.xls")
# colnames(vivax_falciparum_dod) <- gsub(".", " ", colnames(vivax_falciparum_dod), fixed=T)
# colnames(vivax_falciparum_dod) <- gsub("  ", " ", colnames(vivax_falciparum_dod), fixed=T)
# colnames(vivax_falciparum_dod) <- gsub("X", "", colnames(vivax_falciparum_dod), fixed=T)
# colnames(vivax_falciparum_dod) <- gsub(" Associated Genes","% Associated Genes",colnames(vivax_falciparum_dod), fixed=T)


vivax_falciparum_t6 <-read_xls("~/PhD/RNAseq/vac69a/cytoscape/VIVAX_FALCIPARUM_T6/VIVAX_FALCIPARUM_T6_Results_Table")

data <- vivax_falciparum_t6
# data <- t6_data

#filter GOTerms to be level 5 only
dod_level_siz <- subset(data, grepl("5", data$GOLevels))




names_dod <- dod_level_siz %>%
  filter(`% Associated Genes`>10) %>%
  group_by(GOGroups) %>%
  top_n(-3, `Term PValue`) %>%
  ungroup() %>%
  #top_n(-60, `Group PValue`) %>%
  arrange(`Group PValue`, `Term PValue`) %>%
  select(GOTerm, `Term PValue`, GOGroups, `Group PValue`) 

names_dod <- names_dod[!duplicated(names_dod$GOTerm),]

names_dod <- names_dod %>%
  group_by(GOGroups) %>%
  top_n(-1, `Term PValue`) 



subset(names_dod, names_dod$GOGroups==dod_level_siz$GOGroups[grep("*kinetochore*", dod_level_siz$GOTerm)])

