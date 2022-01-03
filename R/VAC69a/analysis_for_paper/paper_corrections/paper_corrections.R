# preamble ####

library(ggplot2)
library(dplyr)
library(tidyr)

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



vac63c_haem <- read.csv("~/PhD/clinical_data/vac63c/VAC063_haem_all_sequenced_WBC_real_percent.csv")

vac63c_haem$Volunteer_code <- gsub("V", "v", vac63c_haem$Volunteer_code)

vac63c_haem$timepoint <- gsub("D+6", "T6", vac63c_haem$timepoint, fixed=T)
vac63c_haem$timepoint <- gsub("C-1", "Baseline", vac63c_haem$timepoint, fixed=T)
vac63c_haem$timepoint <- gsub("C+90", "C90", vac63c_haem$timepoint, fixed=T)

vac63c_lymph <- vac63c_haem %>%
  filter(N_infection=="First", Leukocytes=="Lymphocytes", timepoint %in% c("Baseline", "Diagnosis", "T6", "C90"), Volunteer_code %in% c("v1039", "v1040", "v1061", "v1065", "v1067", "v1068", "v1075", "v313", "v315", "v320", "v6032", "v806", "v818")) %>%
  select(Volunteer_code, trial_number, timepoint, N_infection, Leukocytes, cell_counts) 


vac63c_lymph_premerge <- data.frame("volunteer" = vac63c_lymph$Volunteer_code,
                                    "timepoint" = vac63c_lymph$timepoint,
                                    "Cell" = vac63c_lymph$Leukocytes,
                                    "Frequency" = vac63c_lymph$cell_counts,
                                    "Species" = "P. falciparum")

redo_dod_data <- readxl::read_xls("~/PhD/RNAseq/vac69a/cytoscape/vivax_falci_dod_all_redo/vivax_falci_dod_all_redo/ClueGOResultTable-0.xls")



#whole blood RNAseq normalisation ####

# raw_data <- data.table::fread("/home/flobuntu/PhD/RNAseq/vac69a/All_vivax_falcip_wholeblood_RawTranscriptCounts_19Feb2020.xls")
# 
# vivax_raw <- subset(raw_data, select=c(TRUE, TRUE, grepl("*vivax*", colnames(raw_data))[-c(1,2)]))
# 


all_unique_genes <- data.table::fread("/home/flobuntu/PhD/RNAseq/vac69a/all/xls/all_unique_genes_cleaned.csv", header = T, stringsAsFactors = F)

list_of_all_unique <- split(all_unique_genes, all_unique_genes$file_name)
list_of_all_sig_unique <- lapply(list_of_all_unique, function(x)filter(x, padj<0.05))

DoD_Baseline_sig <- list_of_all_sig_unique[[3]]

super_sig_dod <- subset(DoD_Baseline_sig, abs(DoD_Baseline_sig$log2FoldChange)>log2(3))

write.table(super_sig_dod$Symbol, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, file = "/home/flobuntu/PhD/RNAseq/vac69a/all/xls/gene_lists/vivax_super_sig_dod_fc3.txt")


super_sig_cluego <- readxl::read_xls("/home/flobuntu/PhD/RNAseq/vac69a/cytoscape/super_sig_dod/super_sig_dod_xls")

top_10_go_dod_super_sig <- super_sig_cluego %>%
  slice_min(n=30, order_by=`Term PValue`) %>%
  select(GOTerm)

top_10_go_dod_sig <- unique(redo_dod_data %>%
  slice_min(n=43, order_by=`Term PValue`) %>%
  summarise(GOTerm)
)

table(top_10_go_dod_super_sig$GOTerm %in% top_10_go_dod_sig$GOTerm)
# FALSE  TRUE 
# 21     9


table(super_sig_cluego$GOTerm %in% redo_dod_data$GOTerm)

# FALSE  TRUE 
# 8    22 
subset(super_sig_cluego, !super_sig_cluego$GOTerm %in% redo_dod_data$GOTerm, select = GOTerm)
                                            
# 1 negative regulation of cell cycle G2/M phase transition
# 2 regulation of defense response to virus                
# 3 2'-5'-oligoadenylate synthetase activity               
# 4 regulation of ribonuclease activity                    
# 5 mononuclear cell migration                             
# 6 chemokine receptor binding                             
# 7 granulocyte migration                                  
# 8 chemokine activity  


# so you do get slightly different GO Terms pop out- I reckon this has to do with the disparity of list length- also, slicing out
# all DEGs with low-ish fold changes leaves you with not necessarily a representative list of pathways
# longer lists may be penalised because you need a larger percentage of associated genes found to rise above the noise?
# may need to think about that a littl more..

# vivax falci transcriptomics correction ####

was_used <- scan("~/PhD/vivax_falci_dod_was_used.txt", what="")

dod_unpaired <- read.csv("~/PhD/RNAseq/vac63c/First_DoD_unpaired_all_VAC063data.csv")


#this is the one that was used, only contains A & B
dod_genelist <- read.csv("~/PhD/RNAseq/vac63c/FirstDoD_genelist_padj_0.05.csv")

# this is the one that should have been used
dod_unpaired <- read.csv("~/PhD/RNAseq/vac63c/First_DoD_unpaired_all_VAC063data.csv")

write.table(dod_unpaired$Symbol, file = "~/PhD/RNAseq/vac63c/USE_THIS_VAC63ABC_DOD.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


redo_dod_data <- subset(redo_dod_data, grepl("6", redo_dod_data$GOLevels))

# redo_dod_data$Cluster_Difference <- redo_dod_data$`%Genes Cluster #1.adjusted.to.input.list.size` -redo_dod_data$`%Genes Cluster #2.adjusted.to.input.list.size`
redo_dod_data$Cluster_Difference <- redo_dod_data$`%Genes Cluster #1` -redo_dod_data$`%Genes Cluster #2`

threshold <- 30 # forcing at least a 35/65 split

dod_falci_rich_data <- subset(redo_dod_data, redo_dod_data$Cluster_Difference< -threshold)
dod_vivax_rich_data <- subset(redo_dod_data, redo_dod_data$Cluster_Difference> threshold)

dod_data <- redo_dod_data

vivax_x_limits <- c(80, NA)
vivax_y_limits <- c(20, NA)

dod_dot_plot <- ggplot(dod_data, aes(x=`%Genes Cluster #1`, y=`%Genes Cluster #2`))+
  theme_minimal()+
  scale_color_gradient2(high="#fec200", low="#db0085", midpoint = 0)+
  scale_y_continuous(breaks=seq(0,100,by=10), labels = scales::label_number(suffix="%"), limits = c(-5, 105), expand=c(0,0))+
  scale_x_continuous(breaks=seq(0,100,by=10), labels = scales::label_number(suffix="%"), limits = c(-5, 105), expand=c(0,0))+
  ylab(expression('GO term enrichement'~italic("P. falciparum")~'at Diagnosis'))+
  xlab(expression('GO term enrichement'~italic("P. vivax")~'at Diagnosis'))+
  ggtitle("Diagnosis")+
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

dod_hist_plot <- ggplot(dod_data, aes(y=`%Genes Cluster #2`))+
  theme_minimal()+
  xlab("# GO Terms")+
  scale_y_continuous(breaks = seq(0, 100, by=20), limits=c(-5, 105))+
  scale_x_continuous(breaks = seq(0,80, by=20))+
  geom_histogram(aes(fill = ..y..), orientation = "y", binwidth = 1, color="darkgrey", size=0.15)+
  #stat_bin(binwidth=1, geom="text", aes(label=..count..)) +
  scale_fill_gradient2(low="#fec200", high="#db0085", midpoint = 50)+
  theme(axis.title.y =  element_blank(),
        axis.text.y =  element_blank(),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=6),
        legend.position = "none",
        plot.margin=unit(c(0.5,0.5,0.5,0.3),"cm"))

vivax_falciparum_dod <- cowplot::plot_grid(dod_dot_plot, dod_hist_plot, ncol=2, rel_widths = c(3,1), align="h", axis="bt")

ggsave("~/PhD/manuscripts/vac69a/jci_corrections/vivax_falci_cluego_circle.png", vivax_falciparum_dod,  height = 3, width=3.63)




# redo cluego analysis with very high fold change cutoff#
# show key inflammatory markers with fold cahnge cutoff of at least 4 (or something)



# dod and t6 seperately for ClueGO ####

solo_dod_data <- readxl::read_xls("~/PhD/RNAseq/vac69a/cytoscape/VIVAX_DOD_ALL/old_files/VIVAX_DOD_ALL_RESULTS_TABLE.xls")

solo_t6_data <- readxl::read_xls("~/PhD/RNAseq/vac69a/cytoscape/VIVAX_T6_ALL/VIVAX_T6_ALL_RESULTS_TABLE")

combo_dod_t6_data <- readxl::read_xls("~/PhD/RNAseq/vac69a/cytoscape/VIVAX_BASE_DOD_T6/VIVAX_BASE_DOD_T6_Results_table.xls")

#235
solo_dod_terms <- solo_dod_data$GOTerm

#149
solo_t6_terms <- solo_t6_data$GOTerm

#376 unique
sep_dod_t6_terms <- c(solo_dod_terms, solo_t6_terms)

#289 unique; 88 less than when doing things seperately
combo_dod_t6_terms <- combo_dod_t6_data$GOTerm

# there are 40 specific GO terms that only appear when combining the lists whiwch are a mixture of immune activation, proliferation,
# and general cell biology stuff;
combo_specific <- subset(combo_dod_t6_terms, !combo_dod_t6_terms %in% sep_dod_t6_terms)

#there are a lot (170) of GO terms here that are specific to T6 that don't pop out when combining lists. they are all related to cell cycle
#regulation though, so we're not missing out on biology
sep_specific <- subset(sep_dod_t6_terms, !sep_dod_t6_terms %in% combo_dod_t6_terms)

# i think there's a point to be made that if you wanted to most deeply understand these transcriptomic signatures you should really
# analyse dod and t6 seperately, but the way we're doing it doesn't change the main message that's being told by the data

# it's worth remembering though, that having gene lists of unequal length means that you're drowning out signal possibly- not an issue
# in this dataset, but worth keeping in mind



# full blood count variation ####



fig1_theme <- theme(axis.title.x = element_blank(),
                    legend.title = element_text(size = 9), 
                    legend.text = element_text(size = 9),
                    axis.title=element_text(size=10))





haem_data <- data.table::fread("~/PhD/clinical_data/vac69a/haem.csv")

haem_data$volunteer <- gsub("69010", "v", haem_data$trial_number)


supp_haem_data <- haem_data %>%
  select(volunteer, timepoint, wbc, neutrophils,  monocytes, eosinophils, lymphocytes, platelets) %>%
  #select(trial_number, timepoint, lymphocytes) %>%
  filter(timepoint %in% c("_C_1", "_C1_7", "_C21", "_EP", "_T6", "_C90")) %>%
  gather(Cell, Frequency, c(wbc, neutrophils,  monocytes, eosinophils, lymphocytes, platelets))


supp_haem_data$Cell <- paste(
  toupper(substr(supp_haem_data$Cell, 1,1)),
  substr(supp_haem_data$Cell, 2,nchar(supp_haem_data$Cell)),
  sep="")

supp_haem_data$Cell <- gsub("Wbc", "Total White Cells", supp_haem_data$Cell)

supp_haem_data$Cell <- factor(supp_haem_data$Cell, levels=c("Total White Cells", "Lymphocytes", "Monocytes", "Neutrophils", "Eosinophils", "Platelets"))


bad_timepoints <- c("_C_1", "_C1_7", "_C8_14", "_C21", "_EP", "_T6", "_C90")
great_timepoints <- c("Baseline", "C7", "C14 am", "Diagnosis", "T1", "T6", "C90")

time_dic <- setNames(great_timepoints, bad_timepoints)

supp_haem_data$timepoint <- stringr::str_replace_all(supp_haem_data$timepoint, time_dic)


supp_haem_data <- filter(supp_haem_data, Cell != "Haemoglobin")

wide_haem_data <- pivot_wider(supp_haem_data, names_from = Cell, values_from = Frequency)

wide_haem_perc <- wide_haem_data %>%
  select(Neutrophils, Monocytes, Eosinophils, Lymphocytes)
  
wide_haem_perc <- wide_haem_perc/wide_haem_data$`Total White Cells`
wide_haem_perc <- cbind("volunteer"=wide_haem_data$volunteer, "timepoint"=wide_haem_data$timepoint, wide_haem_perc)

long_haem_perc <- pivot_longer(wide_haem_perc, cols = c(Neutrophils, Monocytes, Eosinophils, Lymphocytes), names_to = "Cell", values_to = "% of CD45+")

long_haem_perc <- filter(long_haem_perc, timepoint%in%c("Baseline", "Diagnosis", "T6", "C90"))

mean_haem <- long_haem_perc %>%
  group_by(Cell, timepoint) %>%
  summarise(Mean=mean(`% of CD45+`))


#long_haem_perc$`% of CD45+` <- ifelse(long_haem_perc$timepoint=="C90", long_haem_perc$`% of CD45+`*6/4, long_haem_perc$`% of CD45+`)

haem_palette <- c("#201F41", "#c90e46", "#eba715", "#eb4b15")
names(haem_palette) <- c("Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils")



summary_long_haem_perc <- long_haem_perc %>%
  group_by(Cell, timepoint) %>%
  summarise("mean"=mean(`% of CD45+`))



(blood_comp_plot <- ggplot(summary_long_haem_perc, aes(x=factor(timepoint, levels=c("Baseline", "C7", "C14 am", "Diagnosis", "T1", "T6", "C90")),
                                              y=mean,
                                              fill=factor(Cell)))+
  geom_bar(position = "stack", stat="identity", color=NA, size=0)+
  scale_y_continuous(label=scales::percent)+
  geom_text(aes(color=Cell, x=timepoint, y=Mean, label=paste(round(Mean*100, digits=2), "%", sep="")), data = mean_haem, position = position_stack(vjust = .5))+
  theme_minimal()+
  ylab("% of White Blood Cells")+
  scale_fill_manual(values=haem_palette[rev(c(1,3,2,4))])+
  scale_color_manual(values = c(rep("black", 3), "white"))+
  guides(fill=guide_legend(title="Cell Type"),
         color=guide_none())+
  ggtitle(expression(paste("Whole Blood Cellular Composition First ", italic("P. vivax"), " Infection")))+
  theme(plot.title = element_text(hjust=0.5),
          axis.text.x = element_text(hjust=1, angle=45, size=8),
          axis.title.x = element_blank(),
          legend.position = "right",
          strip.text = element_text(size=10))
)


ggsave("~/PhD/manuscripts/vac69a/jci_corrections/final_figures/vivax_whole_blood_cell_comp.pdf", blood_comp_plot,  height = 5, width=5.5, bg="white")

long_haem_falci <- vac63c_haem %>%
  filter()

# 
# supp_haem_data_plots <- ggplot(supp_haem_data, aes(x=factor(timepoint, levels=c("Baseline", "C7 am", "C14 am", "Diagnosis", "T1", "T6", "C90")), y=Frequency*1000, color=volunteer, group=volunteer))+
#   scale_fill_manual(values=volunteer_palette)+
#   scale_color_manual(values=volunteer_palette)+
#   geom_line(aes(color=volunteer), size=0.9)+
#   geom_point(fill="white", stroke=1, shape=21, size=0.9)+
#   theme_minimal()+
#   facet_wrap(~Cell, scales="free", nrow = 2)+
#   xlab("Timepoint")+
#   ylab(expression(Cells~"/"~mu*L~blood))+
#   guides(color=guide_legend(title="Volunteer", override.aes = list(size=1)))+
#   fig1_theme+
#   scale_y_continuous(label=scales::comma)+
#   theme(plot.title = element_text(hjust=0.5),
#         axis.text.x = element_text(hjust=1, angle=45, size=8),
#         axis.title.x = element_blank(),
#         legend.position = "right",
#         strip.text = element_text(size=10))
# 

# ggsave("~/PhD/figures_for_thesis/chapter_1/1_haem_counts.png", supp_haem_data_plots, height=4, width=7)


vac63a_haem <- read.csv("~/PhD/clinical_data/vac63a/vac063a_haem.csv")
vac63a_haem$study_name <- "vac63a"

vac63b_haem <- read.csv("~/PhD/clinical_data/vac63b/vac063b_haem.csv")
vac63b_haem$study_name <- "vac63b"

vac63c_haem <- read.csv("~/PhD/clinical_data/vac63c/vac063c_haem.csv")
vac63c_haem$study_name <- "vac63c"

all_colnames <- unique(c(colnames(vac63a_haem),
                         colnames(vac63b_haem),
                         colnames(vac63c_haem)))

harmonise_colnames <- subset(all_colnames, all_colnames %in% colnames(vac63a_haem) & all_colnames %in% colnames(vac63b_haem) & all_colnames %in% colnames(vac63c_haem))

falci_haem <- data.frame(rbind(vac63a_haem[,harmonise_colnames],
                               vac63b_haem[,harmonise_colnames],
                               vac63c_haem[,harmonise_colnames]))


lousy_timepoints <- unique(falci_haem$timepoint)
# [1] "C1_"        "C121_"      "TREATMENT_" "C28_"       "C1_12"      "D7"         "DIAG"       "C_1"        "C90"        "BC28"      
# [11] "CE"         "SCR"        "C1_13"      "DoD"        "C45"        "EV"         "T6" 

# [1] "Baseline"  "C6"        "Diagnosis" "C28"       "C6"        "Screening" "Diagnosis" "Baseline"  "C90"       "C28"      
# [11] "Extra"     "Screening" "C6"        "Diagnosis"       "T6"        "Extra"     "C45"  

great_timepoints <- c("Baseline", "C6", "Diagnosis", "C28", "C6", "Screening", "Diagnosis", "Baseline", "C90", "C28",
                      "Extra", "Screening", "C6", "Diagnosis", "C45", "Extra", "T6")

haem_timepoint_replacement <- setNames(great_timepoints, lousy_timepoints)
falci_haem$timepoint <- stringr::str_replace_all(falci_haem$timepoint, haem_timepoint_replacement)


long_falci_haem <- falci_haem %>%
  select(c(timepoint, trial_number, wbc, neutrophils, monocytes, eosinophils, lymphocytes))%>%
  filter(timepoint %in% c("Baseline", "C6", "Diagnosis", "C28", "T6", "C45")) %>%
  filter(as.character(trial_number) %in% c("6301039", "6301040", "6301061", "6301068", "6301075", "6306032",
                                                       "6301806", "6301818", "6301065", "6301067",
                                                       "6301313", "6301315", "6301320"))

long_falci_perc <- data.frame("timepoint"=long_falci_haem$timepoint, "volunteer"=long_falci_haem$trial_number, long_falci_haem[,4:7]/long_falci_haem$wbc)

long_falci_perc <-  pivot_longer(long_falci_perc, cols = c(neutrophils, monocytes, eosinophils, lymphocytes), names_to = "Cell", values_to = "% of CD45+")

long_falci_perc$timepoint <- ifelse(long_falci_perc$timepoint=="C28", "Memory", long_falci_perc$timepoint)
long_falci_perc$timepoint <- ifelse(long_falci_perc$timepoint=="C45", "Memory", long_falci_perc$timepoint)

long_falci_perc <- filter(long_falci_perc, timepoint %in% c("Baseline", "Diagnosis", "T6", "Memory"))


long_falci_perc$Cell <- paste(
  toupper(substr(long_falci_perc$Cell, 1,1)),
  substr(long_falci_perc$Cell, 2,nchar(long_falci_perc$Cell)),
  sep="")

falci_mean_haem <- long_falci_perc %>%
  group_by(Cell, timepoint) %>%
  summarise(Mean=mean(`% of CD45+`))


long_falci_perc$`% of CD45+` <- ifelse(long_falci_perc$timepoint=="T6", long_falci_perc$`% of CD45+`*13/3, long_falci_perc$`% of CD45+`)

(falci_blood_comp_plot <- ggplot(long_falci_perc, aes(x=factor(timepoint, levels=c("Baseline", "Diagnosis", "T6", "Memory")),
                                               y=`% of CD45+`/13,
                                               fill=factor(Cell)))+
    geom_bar(position = "stack", stat="identity")+
    scale_y_continuous(label=scales::percent)+
    geom_text(aes(x=timepoint, y=Mean, label=paste(round(Mean*100, digits=2), "%", sep="")), data = falci_mean_haem, position = position_stack(vjust = .5))+
    theme_minimal()+
    ylab("% of CD45+")+
    scale_fill_brewer(type="qual", palette = "Dark2", direction=-1)+
    guides(fill=guide_legend(title="Cell Type"))+
    ggtitle("Whole Blood Cellular Composition First Pf Infection")+
    theme(plot.title = element_text(hjust=0.5),
          axis.text.x = element_text(hjust=1, angle=45, size=8),
          axis.title.x = element_blank(),
          legend.position = "right",
          strip.text = element_text(size=10))
)



#memory means C28 for everybody except v313, v315 & 320
ggsave("~/PhD/manuscripts/vac69a/jci_corrections/final_figures/falciparum_whole_blood_cell_comp.pdf", falci_blood_comp_plot,  height = 5, width=5, bg="white")

#individual haem through time plots

longer_falci_haem <- pivot_longer(long_falci_haem, cols = c("wbc", "neutrophils", "monocytes", "eosinophils", "lymphocytes"), names_to = "Cell", values_to = "Cell_Count")

vac63c_prim_indie_haem_plot <- ggplot(subset(longer_falci_haem, longer_falci_haem$trial_number %in% c(6301313, 6301315, 6301320)), aes(x=factor(timepoint, levels=c("Baseline", "C6", "Diagnosis", "T6", "C28", "C45", "C90")), y=Cell_Count, group=factor(trial_number), colour=factor(trial_number)))+
  geom_point()+
  geom_line()+
  facet_wrap(~Cell, scales="free")+
  ylab("10^6 Cells / mL")+
  theme_minimal()+
  guides(colour=guide_legend(title="Volunteer"))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust=0.5))

ggsave("~/PhD/manuscripts/vac69a/jci_corrections/final_figures/vac63c_prim_indie_haem_plot.pdf", vac63c_prim_indie_haem_plot,  height = 5, width=8, bg="white")



# whole blood RNAseq heatmaps ####


# run quick_heatmaps.R lines 1:24 before this

purvesh_list <- c("TGFBI", "DEFA4", "LY86", "BATF",
                  "HK3", "HLA-DPB1")

flo_Inflammation_markers <- list(`Inflammation Markers` = c("STAT1", "STAT2",  "IRF1", "IRF2", "IRF7", "IRF9", "MYD88",
                                                            "TICAM1", "TICAM2", "IDO1", "IDO2", "ACOD1", "GBP1",
                                                            "GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "SOD2", "SOD3",
                                                            "HIF1A", "SECTM1", "ICAM1",
                                                            "CD40", "CD274", "PDCD1LG2"))

flo_More_Inflammation_Markers <- list(`Cytokines & Chemokines` = c("CXCL11", "CXCL10", "CCL2", "CCL25", "IL27", "TNFSF13B", "IL1RN",
                                                                   "TNF", "IL15", "IL1B", "CSF1", "TNFSF13",
                                                                   "IFNA1"
)
)



inflam1 <- quick_gene_heatmaps(flo_Inflammation_markers, sort_by = "DoD_Baseline")$data
#inflam1a <- quick_gene_heatmaps(flo_Inflammation_markers, sort_by = "T6_Baseline")$data

combo_inflam1 <-inflam1

combo_inflam1$file_name <- gsub("_", " vs. ", combo_inflam1$file_name)

combo_inflam1$log2FoldChange <- ifelse(combo_inflam1$padj<=0.05, combo_inflam1$log2FoldChange, 0)

inflam1_plot <- ggplot(combo_inflam1, aes(x=file_name, y=factor(Symbol, levels=inflam1$Symbol[order(inflam1$log2FoldChange, decreasing = FALSE)])))+
  geom_tile(aes(fill=log2FoldChange))+
  viridis::scale_fill_viridis(option = "B", breaks=c(2,4,6), limits=c(0,6.5))+
  theme_void()+
  ggtitle("Acute Phase Response Genes \n")+
  guides(fill=guide_colorbar(title="log2FC", nbin=30))+
  theme(
    axis.text.y = element_text(size=7, angle = 0, hjust=1),
    axis.text.x = element_text(size=7, angle = 0, hjust=0.5),
    plot.title = element_text(hjust=0.5),
    legend.title = element_text(),
    plot.margin=margin(0,0,1,1),
    legend.position = "right"
    #legend.box.margin=margin(0,0,0,0)
  )
ggsave("~/PhD/manuscripts/vac69a/jci_corrections/final_figures/acute_phase_heat.pdf", inflam1_plot, height=5, width = 2.5, bg="white")


  

inflam2 <- quick_gene_heatmaps(flo_More_Inflammation_Markers, sort_by = "DoD_Baseline")$data
#inflam2a <- quick_gene_heatmaps(flo_More_Inflammation_Markers, sort_by = "T6_Baseline")$data
  

combo_inflam2 <- inflam2

combo_inflam2$file_name <- gsub("_", " vs. ", combo_inflam2$file_name)

combo_inflam2$log2FoldChange <- ifelse(combo_inflam2$padj<=0.05, combo_inflam2$log2FoldChange, 0)



inflam2_plot <- ggplot(combo_inflam2, aes(x=file_name, y=factor(Symbol, levels=inflam2$Symbol[order(inflam2$log2FoldChange, decreasing = FALSE)])))+
  geom_tile(aes(fill=log2FoldChange))+
  viridis::scale_fill_viridis(option = "B", breaks=c(2,4,6), limits=c(0,6.5))+
  theme_void()+
  ggtitle("Cytokines & Chemokines \n")+
  guides(fill=guide_colorbar(title="log2FC", nbin=30))+
  theme(
    axis.text.y = element_text(size=7, angle = 0, hjust=1),
    axis.text.x = element_text(size=7, angle = 0, hjust=0.5),
    plot.title = element_text(hjust=0.5),
    legend.title = element_text(),
    plot.margin=margin(0,0,1,1),
    legend.position = "right"
    #legend.box.margin=margin(0,0,0,0)
  )


ggsave("~/PhD/manuscripts/vac69a/jci_corrections/final_figures/cyto_chemo_heat.pdf", inflam2_plot, height=5, width = 2.5, bg="white")




# vivax falciparum parasitaemia comparison ####

vac63c_parasitaemia <- read.csv("~/PhD/clinical_data/vac63c/VAC063_parasitaemias_all.csv")

vac63c_parasitaemia$Timepoint <- paste("C", vac63c_parasitaemia$Timepoint, sep="")

long_vac63c_parasitaemia <- gather(vac63c_parasitaemia, vol_id, parasitaemia, colnames(vac63c_parasitaemia)[2:ncol(vac63c_parasitaemia)])

long_vac63c_parasitaemia$Volunteer <- substr(long_vac63c_parasitaemia$vol_id, 1, 5)
long_vac63c_parasitaemia$N_infection <- ifelse(grepl("First", long_vac63c_parasitaemia$vol_id)==T, "First",
                                               ifelse(grepl("Second", long_vac63c_parasitaemia$vol_id)==T, "Second", "Third"))


long_vac63c_parasitaemia$Volunteer <- gsub("X", "", long_vac63c_parasitaemia$Volunteer)
long_vac63c_parasitaemia$Volunteer <- gsub("_", "", long_vac63c_parasitaemia$Volunteer)

long_vac63c_parasitaemia <- long_vac63c_parasitaemia[!is.na(long_vac63c_parasitaemia$parasitaemia),]


long_vac63c_parasitaemia$Volunteerf <- paste("v", long_vac63c_parasitaemia$Volunteer, sep="")


long_vac63c_parasitaemia_all <- long_vac63c_parasitaemia %>%
  #filter(Volunteer %in% c("313", "315", "320")) %>%
  filter(N_infection %in% c("First"))

long_vac63c_parasitaemia_prim <- long_vac63c_parasitaemia %>%
  filter(Volunteer %in% c("313", "315", "320")) %>%
  filter(N_infection %in% c("First"))



latest_falci_parasitaemia_all <- long_vac63c_parasitaemia_all %>%
  group_by(Volunteerf) %>%
  slice_max(parasitaemia)
  

latest_falci_parasitaemia_prim <- long_vac63c_parasitaemia_prim %>%
  group_by(Volunteerf) %>%
  slice_max(parasitaemia)


(vac63c_indie_paras <- ggplot(data=long_vac63c_parasitaemia, aes(x=Timepoint, y=parasitaemia, group=vol_id))+
    geom_line(aes(color=Volunteerf))+
    geom_point(aes(fill=Volunteerf), shape=21, size=2, stroke=0.1, color="black")+
    scale_y_log10()+
    theme_minimal()+
    #scale_colour_manual(values=vol_pal)+
    #scale_fill_manual(values=vol_pal)+
    ylab("Parasites / mL")+
    xlab("Days Post Infection")+
    guides(color=guide_legend(title = "Volunteer"),
           fill=guide_legend(title = "Volunteer"))+
    theme(#legend.position = "none",
      axis.text = element_text(size=12),
      axis.title = element_text(size=14),
      axis.title.x = element_blank(),
      legend.title = element_text(size=12),
      legend.text = element_text(size=11)))


#vac69 parasitaemia

vivax_para_data <- read.csv("~/PhD/clinical_data/vac69a/parasitaemia/better_vac69a_parasitaemia_ultimate_qc.csv", header=T)
parasitaemias <- gather(vivax_para_data, Timepoint, Genomes, colnames(vivax_para_data)[5:ncol(vivax_para_data)])
parasitaemias$Genomes <- as.numeric(parasitaemias$Genomes)
# get rid of garbage timepoints that mess up graph

parasitaemias$Timepoint <- ifelse(grepl(".5", parasitaemias$Timepoint, fixed = T), paste(parasitaemias$Timepoint, "pm"),  paste(parasitaemias$Timepoint, "am"))
parasitaemias$Timepoint <- gsub(".5", "", parasitaemias$Timepoint, fixed=T)
parasitaemias$Timepoint <- gsub("D0", "Baseline", parasitaemias$Timepoint, fixed=T)

parasitaemias$Timepoint <- gsub("D", "C", parasitaemias$Timepoint, fixed=T)


highest_vivax_parasitaemias <- parasitaemias %>%
  filter(Treatment=="before Treatment") %>%
  group_by(Volunteer) %>%
  slice_max(Genomes)


# change 3 to 13 if including all falciparum volunteers
combo_dod_para_prim <- data.frame("volunteer" = c(highest_vivax_parasitaemias$Volunteer,latest_falci_parasitaemia_prim$Volunteerf),
                             "parasitaemia" = c(highest_vivax_parasitaemias$Genomes, latest_falci_parasitaemia_prim$parasitaemia),
                              "species" = c(rep("P. vivax", 6), rep("P. falciparum", 3))
                             #"species" = c(rep("P. vivax", 6), rep("P. falciparum", 3))
                             
                             )

combo_dod_para_all <- data.frame("volunteer" = c(highest_vivax_parasitaemias$Volunteer,latest_falci_parasitaemia_all$Volunteerf),
                             "parasitaemia" = c(highest_vivax_parasitaemias$Genomes, latest_falci_parasitaemia_all$parasitaemia),
                             "species" = c(rep("P. vivax", 6), rep("P. falciparum", 13))
                             #"species" = c(rep("P. vivax", 6), rep("P. falciparum", 3))
                             
)


wilcox.test(subset(combo_dod_para_all, combo_dod_para_all$species=="P. vivax")$parasitaemia,
            subset(combo_dod_para_all, combo_dod_para_all$species=="P. falciparum")$parasitaemia)
#p-value = 0.6388

combo_dod_para_plot_all <- ggplot(combo_dod_para_all, aes(x=species, y=parasitaemia))+
  geom_boxplot(aes(fill=species))+
  scale_y_log10(limits=c(100,200000), breaks=c(1000,10000,100000))+
  #scale_y_continuous(breaks = c(5000, 10000, 20000, 50000, 100000, 200000))+
  ylab("Peak Genome Copies / mL")+
  #scale_colour_manual(values=subset(vivax_falci_palette, names(vivax_falci_palette) %in% unique(combo_dod_para$volunteer)))+
  scale_shape_manual(values = list("P. vivax"= 16,
                                   "P. falciparum" = 17), guide="none")+
  scale_fill_manual(values=rev(c("#fec200", "#db0085")))+
  # guides(color = guide_legend(title="Volunteer",
  #                             override.aes = list(shape= c(rep(17,3),
  #                                                          rep(16,6)))),
  #        fill = guide_legend(title="Species"))+ 
  theme_minimal()+
  theme(axis.text = element_text(size=12),
    axis.title = element_text(size=14),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=11))

ggsave("~/PhD/manuscripts/vac69a/jci_corrections/final_figures/vivax_falci_dod_para_13v6_log.pdf", combo_dod_para_plot_all, height=4, width=5, bg="white")




combo_dod_para_plot_prim <- ggplot(combo_dod_para_prim, aes(x=species, y=parasitaemia))+
  geom_boxplot(aes(fill=species))+
  ylab("Peak Genome Copies / mL")+
  scale_fill_manual(values=rev(c("#fec200", "#db0085")))+
  scale_y_continuous(limits=c(0,20000))+
  theme_minimal()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=11))


ggsave("~/PhD/manuscripts/vac69a/jci_corrections/final_figures/vivax_falci_dod_para_3v6.pdf", combo_dod_para_plot_prim, height=4, width=5, bg="white")



# p-value = 0.6388
wilcox.test(x=subset(combo_dod_para_all, combo_dod_para_all$species=="P. vivax")$parasitaemia,
            y=subset(combo_dod_para_all, combo_dod_para_all$species=="P. falciparum")$parasitaemia)



# vac69a vac63c first infection lymphopenia ####


vac69a_lymph_premerge <- supp_haem_data %>%
  filter(Cell=="Lymphocytes", timepoint %in% c("Baseline", "Diagnosis", "T6", "C90")) %>%
  mutate("Species"="P. vivax")


vivax_falci_lymph_merge <- rbind(vac69a_lymph_premerge, vac63c_lymph_premerge)

indie_lymph_plot <- ggplot(vivax_falci_lymph_merge, aes(x=factor(timepoint, levels=c("Baseline", "Diagnosis", "T6")), y=Frequency, group=volunteer))+
  geom_point(aes(color=Species))+
  geom_line(aes(color=Species))+
  scale_color_manual(values=rev(c("#fec200", "#db0085")))+
  theme_minimal()+
  xlab("Timepoint")+
  ylab(bquote('Lymphocytes ('*10^6~cells~'/ mL)'))+
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=22),
        axis.title.x = element_blank())


box_lymph_plot_all <- ggplot(vivax_falci_lymph_merge, aes(x=factor(timepoint, levels=c("Baseline", "Diagnosis", "T6", "C90")), y=Frequency))+
    geom_boxplot(aes(fill=Species))+
    #geom_point(aes(color=volunteer, group=Species, shape=Species), position = position_jitterdodge(jitter.width = 0, dodge.width = 0.75), )+
    scale_fill_manual(values=rev(c("#fec200", "#db0085")))+
    scale_colour_manual(na.value = "darkgrey", values=subset(vivax_falci_palette, names(vivax_falci_palette) %in% unique(vivax_falci_lymph_merge$volunteer)))+
    theme_minimal()+
    # scale_shape_manual(values = list("P. vivax"= 16,
    #                                  "P. falciparum" = 17), guide="none")+
    #guides(color = guide_legend(title="Volunteer", override.aes = list(shape=c(rep(17,3), rep(16,6)))))+
    xlab("Timepoint")+
    ylab(bquote('Lymphocytes ('*10^6~cells~'/ mL)'))+
    theme(axis.text = element_text(size=16),
          axis.text.x = element_text(angle = 90, vjust = 0.5),
          axis.title = element_text(size=22),
          legend.title = element_blank(),
          axis.title.x = element_blank())

ggsave("~/PhD/manuscripts/vac69a/jci_corrections/final_figures/vivax_falci_lymphocytes_6v13.pdf", box_lymph_plot_all, height=5.5, width=5.5, bg="white")


vivax_falci_lymph_merge_prim <- filter(vivax_falci_lymph_merge, volunteer %notin% c("v1039", "v1040", "v1061", "v1065", "v1067", "v1068", "v1075", "v6032", "v806", "v818"))

box_lymph_plot_prim <- ggplot(vivax_falci_lymph_merge_prim, aes(x=factor(timepoint, levels=c("Baseline", "Diagnosis", "T6", "C90")), y=Frequency))+
  geom_boxplot(aes(fill=Species))+
  #geom_point(aes(color=volunteer, group=Species, shape=Species), position = position_jitterdodge(jitter.width = 0, dodge.width = 0.75), )+
  scale_fill_manual(values=rev(c("#fec200", "#db0085")))+
  scale_colour_manual(na.value = "darkgrey", values=subset(vivax_falci_palette, names(vivax_falci_palette) %in% unique(vivax_falci_lymph_merge$volunteer)))+
  theme_minimal()+
  # scale_shape_manual(values = list("P. vivax"= 16,
  #                                  "P. falciparum" = 17), guide="none")+
  #guides(color = guide_legend(title="Volunteer", override.aes = list(shape=c(rep(17,3), rep(16,6)))))+
  xlab("Timepoint")+
  ylab(bquote('Lymphocytes ('*10^6~cells~'/ mL)'))+
  theme(axis.text = element_text(size=16),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.title = element_text(size=22),
        axis.title.x = element_blank(),
        legend.title = element_blank())

ggsave("~/PhD/manuscripts/vac69a/jci_corrections/final_figures/vivax_falci_lymphocytes_6v3.pdf", box_lymph_plot_prim, height=5.5, width=5.5, bg="white")



vivax_falci_lymph_merge_prim <- subset(vivax_falci_lymph_merge,
                                       ifelse(vivax_falci_lymph_merge$Species=="P. falciparum",
                                              ifelse(vivax_falci_lymph_merge$volunteer %in% c("v313", "v315", "v320"), TRUE, FALSE), TRUE))

                                       
box_lymph_plot_prim <- ggplot(vivax_falci_lymph_merge_prim, aes(x=factor(timepoint, levels=c("Baseline", "Diagnosis", "T6")), y=Frequency))+
  geom_boxplot(aes(fill=Species))+
  scale_fill_manual(values=rev(c("#fec200", "#db0085")))+
  scale_colour_manual(na.value = "darkgrey", values=subset(vivax_falci_palette, names(vivax_falci_palette) %in% unique(vivax_falci_lymph_merge$volunteer)))+
  theme_minimal()+
  xlab("Timepoint")+
  ylab(bquote('Lymphocytes ('*10^6~cells~'/ mL)'))+
  theme(axis.text = element_text(size=16),
        axis.text.x = element_text(angle = 90, vjust=0.5),
        axis.title = element_text(size=22),
        axis.title.x = element_blank())


ggsave("~/PhD/manuscripts/vac69a/jci_corrections/vivax_falci_lymphocytes_prim.pdf", box_lymph_plot_prim, height=5.5, width=5.5, bg="white")

# alt correlation ####


# 
# #chapter_3_figures.R line 1782; T6 data only
# 
# falci_data <- alt_corr_data
# # correlation_figure.R line 708; T6 data only
# vivax_data <- lineage_freqs
# 
# combo_alt_lineage_acti_data <- data.frame("volunteer"=c(vivax_data$volunteer, as.character(falci_data$volunteer)),
#                                           "lineage"=c(vivax_data$lineage, as.character(falci_data$lineage)),
#                                           "activated"=c(vivax_data$fraction_of_lineage_activated, falci_data$perc_acti*100),
#                                           "alt"=c(vivax_data$alt, falci_data$alt))
#                                           
# combo_alt_lineage_acti_data$lineage <- gsub("gamma delta", "gd", combo_alt_lineage_acti_data$lineage)
# 
#write.csv(combo_alt_lineage_acti_data, "~/PhD/manuscripts/vac69a/jci_corrections/combo_alt_lineage_acti_data.csv", row.names = FALSE)



# here we made a df called vac69b_lineage_activation & vac69b_lin_acti_alt on line 284 in "vac69b_activation.R" to calculate the lineage-wise
# activation at T6 and obtined alt data as well to buff out our first infection dataset

combo_alt_lineage_acti_data <- rbind(combo_alt_lineage_acti_data, vac69b_lin_acti_alt[,1:4])
write.csv(combo_alt_lineage_acti_data, "~/PhD/manuscripts/vac69a/jci_corrections/vac69a_b_vac63c_lineage_acti_data_.csv", row.names = FALSE)


combo_alt_lineage_acti_data <- read.csv("~/PhD/manuscripts/vac69a/jci_corrections/vac69a_b_vac63c_lineage_acti_data_.csv")

combo_alt_lineage_acti_data <- combo_alt_lineage_acti_data%>%
  filter(volunteer %notin% c("v301", "v302", "v304", "v305", "v306", "v307", "v308", "v310"), lineage %in% c("CD4", "Treg"))

#combo_alt_lineage_acti_data <- subset(combo_alt_lineage_acti_data, combo_alt_lineage_acti_data$volunteer %notin% c("v313", "v315", "v320"))

combo_alt_lineage_acti_data %>%
  group_by(lineage) %>%
  do(broom::tidy(cor.test(.$activated, .$alt, method="pearson")))

combo_alt_lineage_acti_data$species <- ifelse(nchar(combo_alt_lineage_acti_data$volunteer)==3, "P. vivax", "P.falciparum")


combo_alt_lineage_acti_data$lineage <- gsub("CD4", "CD4 (r=0.636, p=0.0404)", combo_alt_lineage_acti_data$lineage)
combo_alt_lineage_acti_data$lineage <- gsub("Treg", "Treg (r=0.618, p=0.0478)", combo_alt_lineage_acti_data$lineage)

alt_cd4_treg_plot <- ggplot(combo_alt_lineage_acti_data, aes(x=alt, y=activated/100))+
  geom_point(aes(colour=volunteer, shape=species))+
  geom_smooth(method="lm", colour="black")+
  ylab("Fraction of Lineage Activated")+
  xlab("ALT (IU / L)")+
  # guides(shape = guide_legend(title="Species"),
  #        color = guide_legend(title="Volunteer",
  #                             override.aes = list(shape= c(rep(16,6),
  #                                                          rep(17,3)))),
  # )+
  scale_y_continuous(labels = scales::percent_format(accuracy = 2), limits=c(0,NA))+
  #scale_colour_manual(values=vivax_falci_palette[c(10:15, 1:3)])+
  facet_wrap(~lineage, scales = "free_y")+
  theme_minimal()+
  theme(strip.text=element_text(size=12))

ggsave("~/PhD/manuscripts/vac69a/jci_corrections/final_figures/alt_tcell_activation_correlation_vac69ab.pdf", height=4, width=6.5, bg="white")



l# ALL VOLUNTEERS INCLUDING THIRD INFECTION PEARSON#
# 
# A tibble: 7 x 9
# # Groups:   lineage [7]
# lineage estimate statistic p.value parameter  conf.low conf.high method                              alternative
# <chr>      <dbl>     <dbl>   <dbl>     <int>     <dbl>     <dbl> <chr>                               <chr>      
#   1 CD4       0.555     2.40   0.0319         13  0.0593       0.831 Pearson's product-moment correlati… two.sided  
# 2 CD8       0.512     2.15   0.0512         13 -0.000805     0.811 Pearson's product-moment correlati… two.sided  
# 3 DN       -0.0227   -0.0820 0.936          13 -0.529        0.495 Pearson's product-moment correlati… two.sided  
# 4 gd        0.632     2.94   0.0115         13  0.177        0.864 Pearson's product-moment correlati… two.sided  
# 5 MAIT      0.704     3.57   0.00339        13  0.300        0.894 Pearson's product-moment correlati… two.sided  
# 6 NKT       0.447     1.32   0.227           7 -0.309        0.857 Pearson's product-moment correlati… two.sided  
# 7 Treg      0.685     3.39   0.00487        13  0.265        0.886 Pearson's product-moment correlati… two.sided  
# 

# ALL VOLUNTEERS INCLUDING THIRD INFECTION SPEARMAN#

# # A tibble: 7 x 6
# # Groups:   lineage [7]
# lineage estimate statistic  p.value method                          alternative
# <chr>      <dbl>     <dbl>    <dbl> <chr>                           <chr>      
# 1 CD4       0.821       100. 0.000258 Spearman's rank correlation rho two.sided  
# 2 CD8       0.643       200  0.0117   Spearman's rank correlation rho two.sided  
# 3 DN       -0.0786      604  0.783    Spearman's rank correlation rho two.sided  
# 4 gd        0.479       292  0.0735   Spearman's rank correlation rho two.sided  
# 5 MAIT      0.439       314  0.103    Spearman's rank correlation rho two.sided  
# 6 NKT       0.15        102  0.708    Spearman's rank correlation rho two.sided  
# 7 Treg      0.625       210  0.0149   Spearman's rank correlation rho two.sided  
# 





#ONLY FIRST FALCRIPARUM INFECTION PEARSON#

# A tibble: 7 x 9
# Groups:   lineage [7]
# lineage estimate statistic p.value parameter conf.low conf.high method                               alternative
# <chr>      <dbl>     <dbl>   <dbl>     <int>    <dbl>     <dbl> <chr>                                <chr>      
# 1 CD4        0.411      1.19  0.272          7  -0.348      0.845 Pearson's product-moment correlation two.sided  
# 2 CD8        0.380      1.09  0.313          7  -0.380      0.834 Pearson's product-moment correlation two.sided  
# 3 DN         0.750      3.00  0.0199         7   0.171      0.944 Pearson's product-moment correlation two.sided  
# 4 gd         0.577      1.87  0.104          7  -0.141      0.897 Pearson's product-moment correlation two.sided  
# 5 MAIT       0.707      2.65  0.0332         7   0.0810     0.933 Pearson's product-moment correlation two.sided  
# 6 NKT        0.739      1.10  0.470          1  NA         NA     Pearson's product-moment correlation two.sided  
# 7 Treg       0.621      2.10  0.0743         7  -0.0734     0.910 Pearson's product-moment correlation two.sided  



#ONLY FIRST FALCRIPARUM INFECTION SPEARMAN#

# A tibble: 7 x 6
# Groups:   lineage [7]
# lineage estimate statistic p.value method                          alternative
# <chr>      <dbl>     <dbl>   <dbl> <chr>                           <chr>      
# 1 CD4        0.717        34  0.0369 Spearman's rank correlation rho two.sided  
# 2 CD8        0.483        62  0.194  Spearman's rank correlation rho two.sided  
# 3 DN         0.717        34  0.0369 Spearman's rank correlation rho two.sided  
# 4 gd         0.533        56  0.148  Spearman's rank correlation rho two.sided  
# 5 MAIT       0.467        64  0.213  Spearman's rank correlation rho two.sided  
# 6 NKT        0.5           2  1      Spearman's rank correlation rho two.sided  
# 7 Treg       0.683        38  0.0503 Spearman's rank correlation rho two.sided  



# parasite multiplication rate vac68 vac69 ####
vivax_para_model_data <- filter(parasitaemias, Treatment=="before Treatment", Genomes>100)


vivax_para_model_data$Num_Timepoint <- gsub(" pm", ".3", vivax_para_model_data$Timepoint)
vivax_para_model_data$Num_Timepoint <- gsub(" am", "", vivax_para_model_data$Num_Timepoint)
vivax_para_model_data$Num_Timepoint <- as.numeric(gsub("C", "", vivax_para_model_data$Num_Timepoint))

vivax_para_model_data$Genomes <- log10(vivax_para_model_data$Genomes)




ggplot(vivax_para_model_data, aes(x=Num_Timepoint, y=Genomes, colour=Volunteer))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept = 3)+
  theme_minimal()+
  facet_wrap(~Volunteer, scales="free_x")



vivax_para_model_data_split <- split(vivax_para_model_data, vivax_para_model_data$Volunteer)

vivax_para_models <- lapply(vivax_para_model_data_split, function(x) lm(Genomes~Num_Timepoint, data=x))

results <- data.frame(t(sapply(vivax_para_models, function(x) x$coefficients)))



A <- matrix(c(1, 2, -1unkn, 2), 1, 1)
b <- c(2,1)
matlib::showEqn(A, b)


growth_fun_v2 <- function(x) 10^(results$Num_Timepoint[1]*x)-10
growth_fun_v3 <- function(x) 10^(results$Num_Timepoint[2]*x)-10
growth_fun_v5 <- function(x) 10^(results$Num_Timepoint[3]*x)-10
growth_fun_v6 <- function(x) 10^(results$Num_Timepoint[4]*x)-10
growth_fun_v7 <- function(x) 10^(results$Num_Timepoint[5]*x)-10
growth_fun_v9 <- function(x) 10^(results$Num_Timepoint[6]*x)-10


model_prediction_v2 <- function(x) 10^(-0.5732747+0.2660029*x)
model_prediction_v3 <- function(x) 10^(-0.8611585+0.3905966*x)
model_prediction_v5 <- function(x) 10^(-0.7830211+0.3159702*x)
model_prediction_v6 <- function(x) 10^(-0.8030917+0.2974601*x)
model_prediction_v7 <- function(x) 10^(-1.1575100+0.3294735*x) 
model_prediction_v9 <- function(x) 10^(-1.9519750+0.3445428*x) 

list_of_predictions <- c(model_prediction_v2,
                         model_prediction_v3,
                         model_prediction_v5, 
                         model_prediction_v6, 
                         model_prediction_v7, 
                         model_prediction_v9)

para_plot_list <- vector(mode="list", length=6)



for(i in 1:6){
  data <- vivax_para_model_data_split[[i]]
  data$Genomes <- 10^data$Genomes
  plot <- ggplot(data, aes_(x=quote(Num_Timepoint), y=quote(Genomes), colour=quote(Volunteer)))+
    geom_point()+
    geom_line()+
    scale_color_manual(values=volunteer_palette)+
    scale_y_log10(limits=c(100, 20000), labels=scales::label_scientific())+
    #scale_y_continuous(limits=c(100, 20000), labels=scales::label_scientific())+
    geom_function(fun=list_of_predictions[[i]], linetype = "dashed")+
    theme_minimal()+
    ggtitle(paste("48h PMR ~", round((10^results$Num_Timepoint[i])^2, digits = 2)))+
    theme(legend.position = "none")
  
  para_plot_list[[i]] <- plot
  
}



cowplot::plot_grid(plotlist = para_plot_list, ncol = 3)    



vac68_data <- read.csv("~/PhD/oxford/vac68/parasitaemias.csv")

long_vac68_data <- pivot_longer(vac68_data, cols=2:ncol(vac68_data), names_to = "Num_Timepoint", values_to = "Genomes")

long_vac68_data$Volunteer <- gsub("01-004", "v04", long_vac68_data$Volunteer)
long_vac68_data$Volunteer <- gsub("01-008", "v08", long_vac68_data$Volunteer)

long_vac68_data$Num_Timepoint <- gsub("D", "", long_vac68_data$Num_Timepoint)
long_vac68_data$Num_Timepoint <- gsub(".5", ".3", long_vac68_data$Num_Timepoint)
long_vac68_data$Num_Timepoint <- as.numeric(long_vac68_data$Num_Timepoint)

long_vac68_data$Genomes <- log10(long_vac68_data$Genomes)


long_vac68_data <- na.omit(long_vac68_data)
long_vac68_data <- subset(long_vac68_data, long_vac68_data$Genomes>=2)


vac68_model_data_split <- split(long_vac68_data, long_vac68_data$Volunteer)

vac68_para_models <- lapply(vac68_model_data_split, function(x) lm(Genomes~Num_Timepoint, data=x))

vac_68_results <- data.frame(t(sapply(vac68_para_models, function(x) x$coefficients)))

model_prediction_v4 <- function(x) 10^(-2.812654+0.5103351*x)
model_prediction_v8 <- function(x) 10^(-2.368950+0.4809315*x)


vac68_vac69_data <- c(vivax_para_model_data_split, vac68_model_data_split)

vac68_vac69_results <- rbind(results[,1:2], vac_68_results)



growth_fun_v4 <- function(x) 10^(vac68_vac69_results$Num_Timepoint[7]*x)-10
growth_fun_v8 <- function(x) 10^(vac68_vac69_results$Num_Timepoint[8]*x)-10


vac68_vac69_results$time_to_10x <- rep(0,8)

vac68_vac69_results$time_to_100x[1] <- uniroot(growth_fun_v2, lower = 0, upper = 10)$root
vac68_vac69_results$time_to_100x[2] <- uniroot(growth_fun_v3, lower = 0, upper = 10)$root
vac68_vac69_results$time_to_100x[3] <- uniroot(growth_fun_v5, lower = 0, upper = 10)$root
vac68_vac69_results$time_to_100x[4] <- uniroot(growth_fun_v6, lower = 0, upper = 10)$root
vac68_vac69_results$time_to_100x[5] <- uniroot(growth_fun_v7, lower = 0, upper = 10)$root
vac68_vac69_results$time_to_100x[6] <- uniroot(growth_fun_v9, lower = 0, upper = 10)$root

vac68_vac69_results$time_to_100x[7] <- uniroot(growth_fun_v4, lower = 0, upper = 10)$root
vac68_vac69_results$time_to_100x[8] <- uniroot(growth_fun_v8, lower = 0, upper = 10)$root


list_of_predictions <- c(model_prediction_v2,
                         model_prediction_v3,
                         model_prediction_v5, 
                         model_prediction_v6, 
                         model_prediction_v7, 
                         model_prediction_v9,
                         model_prediction_v4,
                         model_prediction_v8)



para_plot_list <- vector(mode="list", length=8)


vac68_vac69_colours <- list("v02" = "#FB9A99",
                            "v03" = "#E31A1C",
                            "v05" = "#A6CEE3",
                            "v06" = "#1F78B4",
                            "v07" = "#F0E442",
                            "v09" = "#E69F00",
                            "v04" = "deeppink",
                            "v08" = "purple")

vac68_vac69_palette <- unlist(unname(vac68_vac69_colours))
names(vac68_vac69_palette) <- names(vac68_vac69_colours)



for(i in 1:8){
  data <- vac68_vac69_data[[i]]
  data$Genomes <- 10^data$Genomes
  plot <- ggplot(data, aes_(x=quote(Num_Timepoint), y=quote(Genomes), colour=quote(Volunteer)))+
    geom_point()+
    geom_line()+
    ylab("Genomes / mL")+
    xlab("DPI")+
    scale_color_manual(values=vac68_vac69_palette)+
    scale_y_log10(limits=c(100, 30000), labels=scales::label_scientific())+
    scale_x_continuous(breaks = seq(from=8, to=17,by=1))+
    #scale_y_continuous(limits=c(100, 20000), labels=scales::label_scientific())+
    geom_function(fun=list_of_predictions[[i]], linetype = "dashed")+
    theme_minimal()+
    ggtitle(paste(unique(data$Volunteer), "48h PMR ~", round((10^vac68_vac69_results$Num_Timepoint[i])^2, digits = 2)))+
    theme(legend.position = "none",
          axis.title.x = element_text(size=8))
  
  para_plot_list[[i]] <- plot
  
}



all_para_plot <- cowplot::plot_grid(plotlist = para_plot_list, ncol = 3)    

ggsave("~/PhD/manuscripts/vac69a/jci_corrections/vivax_pmr_plot.pdf", width=8, height=6, bg="white")

# reworked figure 2 heatmaps ####
vivax_t6_data <- read.csv("~/PhD/RNAseq/vac69a/all/xls/all_unique_genes_cleaned.csv")

vivax_t6_data$Infection <- "P. vivax"

# 
diana_list <- scan("~/PhD/RNAseq/vac69a/all/xls/gene_lists/diana_cell_cycle_genes.txt", what="NULL", sep = ",")
# 
# 

flo_Inflammation_markers <- c("STAT1", "STAT2",  "IRF1", "IRF2", "IRF7", "IRF9", "MYD88",
                              "TICAM1", "TICAM2", "IDO1", "IDO2", "ACOD1", "GBP1",
                              "GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "SOD2", "SOD3",
                              "HIF1A", "SECTM1", "ICAM1", "CD40", "CD274", "PDCD1LG2", "CXCL11",
                              "CXCL10", "CCL2", "CCL25", "IL27", "TNFSF13B", "IL1RN", "TNF",
                              "IL15", "IL1B", "CSF1", "TNFSF13", "IFNA1")
)


# inflammatory_list <- c("STAT1", "STAT2",  "IRF1", "IRF2", "IRF7", "IRF9", "MYD88",
#                        "TICAM1", "TICAM2", "TLR4", "IDO1", "IDO2", "ACOD1", "GBP1",
#                        "GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "SOD1", "SOD2", "SOD3",
#                        "S100A8", "S100A9", "HIF1A", "HMOX1", "HMOX2", "SECTM1", "ICAM1",
#                        "CD40", "PDCD1", "CD274", "PDCD1LG2", "CXCL11", "CXCL10", "CCL2",
#                        "CCL25", "IL27", "CCL23", "TNFSF13B", "IL1RN", "TNF", "IL15",
#                        "IL1B", "CSF1", "TNFSF13", "TGFB1", "IL1A", "IL18", "IL18", "IL7", 
#                        "IL12B", "CSF2", "LTA", "IFNB1", "IFNA1", "IL10", "IL6", "IL12A",
#                        "CXCL8", "AIM2", "OASL", "MX1", "MX2", "DDX58", "IL18", "CXCL9", "IL6", "IL21")
# 
# #shorter
# inflammatory_list <- c("STAT1", "STAT2",  "IRF1", "IRF2", "IRF7", "IRF9", "MYD88",
#                        "TLR4", "IDO1", "IDO2", "ACOD1", "GBP1",
#                        "GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "SOD1", "SOD2", "SOD3",
#                        "S100A8", "S100A9", "HIF1A", "HMOX1", "HMOX2", "ICAM1",
#                        "CD40", "PDCD1", "CD274", "PDCD1LG2", "CXCL11", "CXCL10", "CCL2",
#                        "CCL25", "IL27", "CCL23", "TNFSF13B", "IL1RN", "TNF", "IL15",
#                        "IL1B", "CSF1", "TNFSF13", "TGFB1", "IL1A", "IL18", "IL18", "IL7", 
#                        "IL12B", "CSF2", "LTA", "IFNB1", "IFNA1", "IL10", "IL6", "IL12A",
#                        "CXCL8", "AIM2", "OASL", "MX1", "MX2", "DDX58", "IL18", "CXCL9", "IL6", "IL21")



inflammatory_list <- c("IL1RN","IL1B",  "TNF", "IL15", "IL27", "STAT1", "STAT2", "GBP1", "GBP2", "OASL", "MX1", "MX2",
                       "IRF1", "IRF9", "MYD88", "DDX58", "SOD2", "HIF1A", "CD274", "PDCD1LG2", "CASP1", "CASP3","CASP8",
                       "BAK1", "BCL2L1", "BCL2L13", "ADAM17", "CCR1", "CCR2", "CCL2", "CCL3", #"CCL4", "CCL5",
                       "CSF1", "CXCL10", "IDO1", "TNFa", "TLR", "NLRP", "MDA5", "MAVS", "AIM2")

diana_list <- unlist(strsplit(c("POLA1, MKI67, POLA2, POLD1, POLD3, POLE2, PCNA, FANCE, FANCG, FANCI, FANCL, RAD51, RAD51AP1, RAD54L, CNE2, CDC25A, CDC6, CCNA2, CCNB1, CCNB2, CDC25C, CENPA, CENPE, CENPF, CENPI"),
                              split = ", ",))


tp_keep <- "T6_Baseline";chosen_list <- diana_list
tp_keep <- "DoD_Baseline";chosen_list <- inflammatory_list



faves_in_vivax <- vivax_t6_data %>%
  #filter(file_name %in% tp_keep) %>%
  #filter(file_name == "T6_Baseline") %>%
  #filter(padj<=0.05)%>%
  #filter(abs(log2FoldChange)>=log2(1.5))%>%
  select(Symbol, log2FoldChange, padj, file_name) %>%
  arrange(log2FoldChange)


phil_cycle_faves <- c("CDC25A", "CCNB2", "CDK1", "CDC20", "RAD51", "BIRC5", "CEP55", "NEK2", "CENPI", "CDC6", "MKI67", "EXO1", "AURKB", "CENPF", "EZH2", "HELLS")

chosen_list <- unique(c(diana_list, phil_cycle_faves, "STAT1", "STAT2",  "IRF1", "IRF7", "MYD88",
                        "IDO1", "ACOD1", "GBP1",
                        "GBP2", "GBP4", "GBP5", "GBP6", "SOD2",
                        "HIF1A", "SECTM1", "ICAM1", "CD40", "CD274", "PDCD1LG2",
                        "CXCL10", "CCL2", "CCL25", "IL27", "IL1RN", "TNF",
                        "IL15", "IL1B", "CSF1"))




#combo_data <- subset(faves_in_vivax, faves_in_vivax$Symbol %in% c(chosen_list))
combo_data <- subset(faves_in_vivax, faves_in_vivax$Symbol %in% purvesh_list)

combo_data$log2FoldChange <- ifelse(combo_data$padj>0.05, 0, combo_data$log2FoldChange)
combo_data$log2FoldChange <- ifelse(abs(combo_data$log2FoldChange)<log2(1.5), 0, combo_data$log2FoldChange)


gene_matrix <- data.frame(combo_data %>%
                            arrange(desc(log2FoldChange)) %>%                          
                            select(Symbol, file_name, log2FoldChange) %>%
                            filter(log2FoldChange>log2(1.5))%>%
                            pivot_wider(names_from = file_name, values_from = log2FoldChange))

gene_matrix <- gene_matrix[gene_matrix[,2]!=0,]

#gene_matrix <- gene_matrix[ifelse(gene_matrix$T6_Baseline==0&gene_matrix$DoD_Baseline==0, FALSE, TRUE), ]

rownames(gene_matrix) <- gene_matrix$Symbol

gene_matrix$Symbol <- NULL

gene_matrix <- as.matrix(gene_matrix)


inferno <- c(colorspace::sequential_hcl("inferno", n=10))


col_fun_rna_t6 <- circlize::colorRamp2(breaks=seq(from = -0.5,
                                                  to = max(gene_matrix),
                                                  by = max(gene_matrix)/8),
                                       colors=inferno)



# rownames(gene_matrix) <- gsub("_", " relative to\n ", rownames(gene_matrix))
# rownames(gene_matrix) <- gsub("DoD", "Diagnosis", rownames(gene_matrix))

colnames(gene_matrix) <- gsub("_", " relative to\n ", colnames(gene_matrix))

gene_heatmap <- Heatmap(matrix = gene_matrix,
                        cluster_rows = FALSE,
                        show_heatmap_legend = TRUE,
                        column_title ="Selected Cell\nCycle Genes",
                        column_title_gp = gpar(just="center", fontsize=10),
                        cluster_columns = FALSE,
                        row_names_gp = gpar(fontsize = 10, fontface="italic"),
                        heatmap_legend_param = list(legend_position = "bottom",
                                                    col=col_fun_rna_dod,
                                                    title = "log2FC",
                                                    legend_direction = "horizontal",
                                                    title_position = "topcenter",
                                                    #legend_width = unit(6.2, "cm"),
                                                    border = FALSE),
                        column_names_gp = gpar(fontsize = 8, just="center"),
                        show_row_names = TRUE,
                        row_names_side = "left",
                        col = col_fun_rna_dod,
                        column_names_rot = 45)



pdf("/home/flobuntu/PhD/manuscripts/vac69a/jci_corrections/vivax_t6_heat.pdf", width = 2.5, height=8.2)
draw(gene_heatmap,  heatmap_legend_side = "bottom",
     padding=unit(c(2,8,2,8), "mm")
)
dev.off()





faves_in_vivax_dod <- vivax_t6_data %>%
  #filter(file_name %in% tp_keep) %>%
  filter(file_name == "DoD_Baseline") %>%
  #filter(padj<=0.05)%>%
  #filter(abs(log2FoldChange)>=log2(1.5))%>%
  select(Symbol, log2FoldChange, padj, file_name) %>%
  arrange(log2FoldChange)



combo_data2 <- subset(faves_in_vivax_dod, faves_in_vivax_dod$Symbol %in% purvesh_list)

#combo_data2 <- subset(faves_in_vivax_dod, faves_in_vivax_dod$Symbol %in% chosen_list)

combo_data2$log2FoldChange <- ifelse(combo_data2$padj>0.05, 0, combo_data2$log2FoldChange)
combo_data2$log2FoldChange <- ifelse(abs(combo_data2$log2FoldChange)<log2(1.5), 0, combo_data2$log2FoldChange)



gene_matrix2 <- data.frame(combo_data2 %>%
                            arrange(desc(log2FoldChange)) %>%                          
                            select(Symbol, file_name, log2FoldChange) %>%
                            filter(log2FoldChange>log2(1.5))%>%
                            pivot_wider(names_from = file_name, values_from = log2FoldChange))

gene_matrix2 <- gene_matrix2[gene_matrix2[,2]!=0,]

#gene_matrix <- gene_matrix[ifelse(gene_matrix$T6_Baseline==0&gene_matrix$DoD_Baseline==0, FALSE, TRUE), ]

rownames(gene_matrix2) <- gene_matrix2$Symbol

gene_matrix2$Symbol <- NULL

gene_matrix2 <- as.matrix(gene_matrix2)



colnames(gene_matrix2) <- gsub("_", " relative to\n ", colnames(gene_matrix2))



col_fun_rna_dod <- circlize::colorRamp2(breaks=seq(from = -0.5,
                                                   to = max(gene_matrix2),
                                                   by = max(gene_matrix2)/9),
                                        colors=inferno)




gene_heatmap2 <- Heatmap(matrix = gene_matrix2,
                        cluster_rows = FALSE,
                        show_heatmap_legend = TRUE,
                        column_title ="Selected \nInflammatory Genes",
                        column_title_gp = gpar(just="center", fontsize=10),
                        cluster_columns = FALSE,
                        row_names_gp = gpar(fontsize = 10, fontface="italic"),
                        heatmap_legend_param = list(legend_position = "bottom",
                                                    col=col_fun_rna_dod,
                                                    title = "log2FC",
                                                    legend_direction = "horizontal",
                                                    title_position = "topcenter",
                                                    #legend_width = unit(6.2, "cm"),
                                                    border = FALSE),
                        column_names_gp = gpar(fontsize = 8, just="center"),
                        show_row_names = TRUE,
                        row_names_side = "left",
                        col=col_fun_rna_dod,
                        column_names_rot = 45)



pdf("/home/flobuntu/PhD/manuscripts/vac69a/jci_corrections/vivax_dod_heat.pdf", width = 2.5, height=8.2)
draw(gene_heatmap2,  heatmap_legend_side = "bottom",
     padding=unit(c(2,8,2,8), "mm")
)
dev.off()


# all genese heatmap checkerboard

combo_data2 <- combo_data2[order(combo_data2$log2FoldChange),]

all_genes_df <- as.matrix(data.frame("T6_Baseline"=as.numeric(combo_data$log2FoldChange)[match(combo_data2$Symbol, combo_data$Symbol)],
                                     "Diagnosis_Baseline"=as.numeric(combo_data2$log2FoldChange)
                                     
))
rownames(all_genes_df) <- combo_data2$Symbol

all_genes_df <- subset(all_genes_df, !(all_genes_df[,1]==0&all_genes_df[,2]==0))

all_genes_matrix <- t(all_genes_df)

rownames(all_genes_matrix) <- gsub("_", " relative to\n ", rownames(all_genes_matrix))



# col_fun_all <- circlize::colorRamp2(c(-max(all_genes_df), 0, max(all_genes_df)), c("#0859C6", "black", "#FFA500"))
# 
# col_fun_rna_all<- circlize::colorRamp2(breaks=seq(from = -0.5,
#                                                    to = max(gene_matrix2),
#                                                    by = max(gene_matrix2)/9),
#                                         colors=c("#0859C6", inferno[1:9]))
# 


col_fun_rna_all <- circlize::colorRamp2(breaks=seq(from = 0,
                                                   to = max(all_genes_matrix),
                                                   by = max(all_genes_matrix)/9),
                                        colors=inferno)



gene_heatmap3 <- Heatmap(matrix = all_genes_matrix,
                         cluster_columns = FALSE,
                         cluster_rows=FALSE,
                         show_heatmap_legend = TRUE,
                         show_column_dend = FALSE,
                         #column_title ="Selected \nInflammatory Genes",
                         #column_title_gp = gpar(just="center", fontsize=10),
                         column_names_gp = gpar(fontsize = 10, fontface="italic"),
                         heatmap_legend_param = list(legend_position = "bottom",
                                                     #col=col_fun_all,
                                                     title = "log2FC",
                                                     legend_direction = "horizontal",
                                                     title_position = "topcenter",
                                                     #legend_width = unit(6.2, "cm"),
                                                     border = FALSE),
                         show_row_names = TRUE,
                         row_names_side = "left",
                         col=col_fun_rna_all,
                         column_names_rot = 90)



pdf("/home/flobuntu/PhD/manuscripts/vac69a/jci_corrections/vivax_dod_t6_heat3.pdf", width = 8, height=3)
draw(gene_heatmap3,  heatmap_legend_side = "bottom"
     #padding=unit(c(2,8,2,8), "mm")
)
dev.off()





