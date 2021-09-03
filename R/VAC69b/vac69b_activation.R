library(CATALYST)
#library(diffcyt)
# library(dplyr)
# library(tidyr)
#library(ggplot2)

`%notin%` <- Negate(`%in%`)


setwd("~/PhD/cytof/vac69b/T_cells_only/comped/recomped")
# 
# fcs <- list.files(path = "~/PhD/cytof/vac69b/T_cells_only", pattern = "fcs")
# 
# #vac69b_flowset <- flowCore::read.flowSet(fcs)
# 
# 
# 
# md <- data.frame("file_name"=fcs,
#                  "volunteer"= substr(fcs, 1,3),
#                  "timepoint" = substr(fcs, nchar(fcs)-6, nchar(fcs)-4))
# 
# md$timepoint <- gsub("ine", "baseline", md$timepoint)
# md$sample_id <- paste(md$volunteer, md$timepoint, sep="_")
# md$batch <- ifelse(md$volunteer %in% c("v11", "v21"), "batch_1", "batch_2")
# 
# write.csv(md, "metadata.csv", row.names = F)

#split clustering to accomodate your laptops lousy 8gb of ram

md <- read.csv("~/PhD/cytof/vac69b/T_cells_only/metadata.csv")
md <- subset(md, md$timepoint %in% c("baseline", "dod", "ep6", "c56"))


# the internal panel has lost the Gaussian channels and acquired the cluster_id channel-
# we need to make the internal panel match with the one we pass to CATALYST so let's quickly fix that
#internal_panel <- colnames(flowCore::read.FCS(md$file_name[1]))

panel <- read.csv("~/PhD/cytof/vac69b/T_cells_only/vac69b_panel.csv", header=T)
panel[63, 1]="icd"
# # this neatly graps all the channels associated with metals i.e. everything but gaussian & time
# panel <- subset(panel, grepl("Di", panel$fcs_colname))
# panel <- rbind(panel, c("cluster_id", "cluster_id", "none"))
# 
# write.csv(panel, "vac69b_panel.csv", row.names = F)

sce <- prepData(x=md$file_name, panel, md, md_cols =
                  list(file = "file_name", id = "sample_id", factors = c("batch", "volunteer", "timepoint")))


non_t_cell_markers <- c("CD45", "CD3",  "CD14", "CD16")

refined_markers <- c("CD4",
                     "CD8",
                     "Vd2",
                     "Va72",
                     "CXCR5",
                     "CD38",
                     #"CD69",
                     "HLADR",
                     "ICOS",
                     "CD28",
                     "PD1",
                     #"TIM3",
                     "CD95",
                     "BCL2",
                     "CD27",
                     "Perforin",
                     "GZB",
                     #"TCRgd",
                     "Tbet",
                     "Eomes",
                     #"RORgt",
                     #"GATA3",
                     "CTLA4",
                     "Ki67",
                     "CD127",
                     "CD56",
                     #"CD16",
                     "CD161",
                     "CD49d",
                     "CD25",
                     "FoxP3",
                     "CD39",
                     "CX3CR1",
                     "CD57",
                     "CD45RA",
                     "CD45RO",
                     "CCR7")


all_markers <- c(non_t_cell_markers, refined_markers)



set.seed(123);sce <- CATALYST::cluster(sce, features = refined_markers, xdim = 10, ydim = 10, maxK = 45)



#t_cell_phenoh <- plotExprHeatmap(x = sce, by = "cluster", row_clust = TRUE, col_clust = FALSE, k = "meta50", bars = TRUE, features = all_markers)
t_cell_pheno <- plotExprHeatmap(x = sce, by = "cluster", row_clust = FALSE, col_clust = FALSE, k = "meta45", bars = TRUE, features = all_markers)

#ctrl_whole_blod_cluster_phenotype <- plotExprHeatmap(x = sce, by = "cluster", row_clust = TRUE, col_clust = FALSE, k = "meta40", bars = TRUE, features = all_markers)

# pdf("./figures/t_cell_phenoh.pdf", height = 8, width = 9)
# t_cell_phenoh
# dev.off()


pdf("./figures/t_cell_pheno2.pdf", height = 8, width = 9)
t_cell_pheno
dev.off()



meta45_table <- read.csv("/home/flobuntu/PhD/cytof/vac69b/T_cells_only/rough_vac69b_meta45_merging_table.csv", header=T)
merged_daf <- CATALYST::mergeClusters(sce, k = "meta45", table = meta45_table, id = "flo_merge", overwrite = TRUE)

# playing with gating ####
t6 <- filterSCE(sce, timepoint=="baseline")

t6_fcs <- sce2fcs(sce, split_by = "sample_id")

ggcyto(t6_fcs, aes(x="PD1", y="CXCR5"))+
  geom_hex(bins = 60)+
  #geom_density2d(colour = "black", bins=13)+
  scale_y_flowCore_fasinh(a=5)+
  scale_x_flowCore_fasinh(a=5)



vac69b_gs <- GatingSet(t6_fcs)


gs_add_gating_method(vac69b_gs,
                     alias = "CD4_CD45RO",
                     pop = "+", parent = "root",
                     dims = "Nd145Di,Ho165Di",
                     gating_method = "boundary",
                     gating_args = "min = c(30, 30), max=c(400,400)")


gs_add_gating_method(vac69b_gs,
                     alias = "Tfh",
                     pop = "+", parent = "CD4_CD45RO",
                     dims = "Gd155Di, Sm149Di",
                     gating_method = "boundary",
                     gating_args = "min = c(30, 30), max=c(500,500)")


# gs_add_gating_method(vac69b_gs,
#                      alias = "Tfh",
#                      pop = "+", parent = "CD4_memory",
#                      dims = "Sm149Di",
#                      gating_method = "mindensity")
#
#
# df <- gs_pop_get_stats(vac69b_gs,
#                        type = "percent",
#                        nodes = c("CD4_memory", "Tfh"))
#
# df
#

ggcyto(vac69b_gs, aes(x="CXCR5", y="PD1"), subset="CD4_CD45RO")+
  #geom_hex(bins=60)+
  facet_null()+
  facet_wrap(~name, nrow=5)+
  geom_point(alpha=0.1, size=0.8)+
  geom_gate("Tfh")+
  geom_stats()+
  scale_y_flowCore_fasinh(a=5)+
  scale_x_flowCore_fasinh(a=5)






# differential abundance for plotting ####

library(diffcyt)

ei <- metadata(sce)$experiment_info
vac69b_design <- createDesignMatrix(ei, c("timepoint", "volunteer"))
colnames(vac69b_design)
  
  
t6_contrast <- createContrast(c(c(0,0,0,1), rep(0, 4)))
  
  
  
da_t6_all <- diffcyt(merged_daf,
                         design = vac69b_design,
                         contrast = t6_contrast,
                         analysis_type = "DA",
                         method_DA = "diffcyt-DA-edgeR",
                         clustering_to_use = "flo_merge",
                         verbose = T)
  

diffy_data <- vac69a.cytof::diffcyt_boxplot(da_t6_all, merged_daf, FDR=1, logFC=0)$data 
diffy_data_count <- vac69a.cytof::diffcyt_boxplot(da_t6_all, merged_daf, FDR=1, logFC=0, counts = TRUE)$data 

write.csv(diffy_data, "/home/flobuntu/PhD/cytof/vac69b/T_cells_only/comped/recomped/vac69b_meta45_all_freqs.csv", row.names = FALSE)
write.csv(diffy_data_count, "/home/flobuntu/PhD/cytof/vac69b/T_cells_only/comped/recomped/vac69b_meta45_all_counts.csv", row.names = FALSE)


# figures ####

library(dplyr)
library(ggplot2)
library(tidyr)

vac69b_data <- read.csv("/home/flobuntu/PhD/cytof/vac69b/T_cells_only/comped/recomped/vac69b_meta45_all_freqs.csv", header = T)
vac69a_data <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_cluster_freqs.csv", header = T)



lousy_timepoints <- unique(vac69b_data$timepoint)

good_timepoints <- c("Baseline", "C56", "Diagnosis", "T6")
timepoint_replacement <- setNames(good_timepoints, lousy_timepoints)

vac69b_data$timepoint <- stringr::str_replace_all(vac69b_data$timepoint, timepoint_replacement)
vac69b_data$lineage <- ifelse(grepl("CD4", vac69b_data$cluster_id), "CD4", NA)
vac69b_data$lineage <- ifelse(grepl("CD8", vac69b_data$cluster_id), "CD8", vac69b_data$lineage)
vac69b_data$lineage <- ifelse(grepl("Treg", vac69b_data$cluster_id), "Treg", vac69b_data$lineage)
vac69b_data$lineage <- ifelse(grepl("MAIT", vac69b_data$cluster_id), "MAIT", vac69b_data$lineage)
vac69b_data$lineage <- ifelse(grepl("gd", vac69b_data$cluster_id), "gd", vac69b_data$lineage)
vac69b_data$lineage <- ifelse(grepl("DN", vac69b_data$cluster_id), "DN", vac69b_data$lineage)
vac69b_data$lineage <- ifelse(grepl("DP", vac69b_data$cluster_id), "DP", vac69b_data$lineage)

vac69a_data$lineage <- ifelse(grepl("CD4", vac69a_data$cluster_id), "CD4", NA)
vac69a_data$lineage <- ifelse(grepl("CD8", vac69a_data$cluster_id), "CD8", vac69a_data$lineage)
vac69a_data$lineage <- ifelse(grepl("Treg", vac69a_data$cluster_id), "Treg", vac69a_data$lineage)
vac69a_data$lineage <- ifelse(grepl("MAIT", vac69a_data$cluster_id), "MAIT", vac69a_data$lineage)
vac69a_data$lineage <- ifelse(grepl("gamma delta", vac69a_data$cluster_id), "gd", vac69a_data$lineage)
vac69a_data$lineage <- ifelse(grepl("DN", vac69a_data$cluster_id), "DN", vac69a_data$lineage)
vac69a_data$lineage <- ifelse(grepl("DP", vac69a_data$cluster_id), "DP", vac69a_data$lineage)

vac69a_data$volunteer <- gsub("V", "v", vac69a_data$volunteer, fixed=T)
vac69a_data$timepoint <- gsub("DoD", "Diagnosis", vac69a_data$timepoint, fixed=T)

vac69a_data$n_infection <- "First"
vac69b_data$n_infection <- ifelse(vac69b_data$volunteer %in% c("v11", "v21"), "First", "Second")



combo_data <- rbind(select(vac69a_data, cluster_id, sample_id, volunteer, timepoint, frequency, lineage, n_infection),
                    select(vac69b_data, cluster_id, sample_id, volunteer, timepoint, frequency, lineage, n_infection)
                    )

plottable_data <- filter(combo_data, timepoint %in% c("Baseline, Diagnosis", "T6"))

non_naive_perc <- plottable_data %>%
  group_by(volunteer, timepoint, n_infection) %>%
  filter(grepl("non-naive", cluster_id, fixed = F))%>%
  summarise(sum(frequency))
  

plottable_data <- subset(plottable_data, grepl("activated", plottable_data$cluster_id))  
plottable_data$volunteer <- factor(plottable_data$volunteer)

lineage_palette <- c("#AA3377", "#EE6677", "#4477AA", "#66CCEE", "#228833", "#FFA500", "#BBBBBB")
names(lineage_palette) <-c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "Resting")



volunteer_colours <- list("v02" = "#FB9A99",
                          "v03" = "#E31A1C",
                          "v05" = "#A6CEE3",
                          "v06" = "#1F78B4",
                          "v07" = "#F0E442",
                          "v09" = "#E69F00",
                          "v11" = "#c3015c",
                          "v21" = "#6d0133")



volunteer_palette <- unlist(unname(volunteer_colours))
names(volunteer_palette) <- names(volunteer_colours)





activation_stacked_barchart <- ggplot(plottable_data, aes(x=volunteer, y=frequency/100, fill=lineage))+
  geom_bar(stat="identity", position="stack")+
  #geom_text(aes(y=(total_cd3/100)+0.01, label=paste0(round(total_cd3, digits = 1), "%", sep='')))+
  theme_minimal()+
  facet_wrap(timepoint~n_infection)+
  ggtitle("Overall T cell activation")+
  scale_fill_manual(values=lineage_palette, labels=c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "Resting"))+
  scale_y_continuous(name = "Fraction of CD3+ T cells activated", labels=scales::percent_format(accuracy = 1))+
  #ylim(0,25)+
  #geom_text(aes(label=cluster_id), position = position_stack(vjust = .5))+
  theme(#legend.position = "none",
    plot.title = element_text(hjust=0.5, size=8),
    strip.text = element_text(),
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(hjust=0.5, angle=45),
    panel.spacing.x = unit(0.8,"lines"),
    axis.title.y = element_text(size=7),
    panel.grid.minor.y = element_blank(),
    legend.position = "none",)



#### sandbox ####

plottable_data <- filter(plottable_data, lineage=="CD4")  

prim_summary_t6 <-  filter(plottable_data, n_infection=="First")
  
prim_sig_activation_lolli <- ggplot()+
    geom_bar(data=prim_summary_t6,aes(y=factor(volunteer, levels=c("v02", "v03", "v05", "v06", "v07", "v09", "v11", "v21")), x=frequency/100, fill=lineage), stat="identity", position="stack", width=1)+
    theme_minimal()+
    ggtitle("T Cell Activation at T6 First Vivax Infection")+
    scale_fill_manual(values=lineage_palette)+
    #scale_x_reverse(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1), limits=c(0.25, 0), breaks=c(25, 20, 15, 10, 5, 0)/100)+
    scale_x_reverse(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1), na.value=0, limits=c(0.15, 0),  breaks=c(15, 10, 5, 0)/100)+
    guides(fill=guide_legend(nrow = 6,keyheight = unit(4, "mm"), keywidth = unit(8, "mm")))+
    theme(plot.title = element_text(hjust=0.5, size=11),
          strip.text = element_text(hjust=0.5, size=10, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(color = "white"),
          strip.placement = "outside",
          #plot.margin = unit(c(1,0,1,5), "mm"),
          legend.position = "bottom",
          legend.text = element_text(size=6),
          legend.title = element_blank())
  
  
  cluster_leg <- cowplot::get_legend(prim_sig_activation_lolli)
  
  prim_sig_activation_lolli <- prim_sig_activation_lolli+theme(legend.position = "none")
  
  
  ter_summary_t6 <-  filter(plottable_data, n_infection=="Second")
  
  
 ter_sig_activation_lolli <- ggplot()+
    geom_bar(data=ter_summary_t6,aes(y=factor(volunteer, levels=c("v02", "v03", "v05", "v06", "v07", "v09", "v11", "v21")), x=frequency/100, fill=lineage), stat="identity", position="stack", width=1)+
    theme_minimal()+
    ggtitle("T Cell Activation at T6 Second Vivax Infection")+
    scale_fill_manual(values=lineage_palette)+
    scale_y_discrete(drop=FALSE)+
    #scale_x_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1), na.value=0, limits=c(0, 0.25),  breaks=c(25, 20, 15, 10, 5, 0)/100)+
    scale_x_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1), na.value=0, limits=c(0, 0.15),  breaks=c(15, 10, 5, 0)/100)+
    guides(fill=guide_legend(reverse = TRUE, nrow = 5,keyheight = unit(2, "mm"), keywidth = unit(4, "mm")))+
    theme(plot.title = element_text(hjust=0.5, size=11),
          strip.text = element_text(hjust=0.5, size=10, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.placement = "outside",
          #plot.margin = unit(c(1,5,1,0), "mm"),
          legend.position = "none",
          legend.text = element_text(size=6),
          legend.title = element_blank())
  
  
  

combo_lolli <- cowplot::plot_grid(prim_sig_activation_lolli, ter_sig_activation_lolli, cluster_leg, nrow=1, rel_widths = c(5,5,1),align="v", axis = "tbrl")
ggsave("/home/flobuntu/PhD/cytof/vac69b/T_cells_only/comped/recomped/figures/vac69a_b_activation_lollipop2.pdf", combo_lolli, height=4, width=7.5)
 
  
