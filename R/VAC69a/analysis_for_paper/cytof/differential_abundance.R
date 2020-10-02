## extra script, just for stats ###

library(CATALYST)
library(diffcyt)
library(vac69a.cytof)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(dplyr)
library(scales)
library(ggplot2)

#daf <- read_small("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/", proportional = T, event_number = 1000)
merged_daf <- read_full("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")

#merged_daf <- filterSCE(merged_daf, timepoint != "C10")


ei <- metadata(merged_daf)$experiment_info
#ei$timepoint <- factor(ei$timepoint, levels=c("C10", "Baseline", "DoD", "T6"))

# alt <- ifelse(ei$timepoint=="T6"&ei$volunteer!="V05", TRUE, FALSE)
# ei$alt <- alt
#ei$alt <- ei$alt>40
design <- createDesignMatrix(ei, c("timepoint", "volunteer"))
# 

#design <- model.matrix(~ei$time+ei$volunteer:ei$alt)
# batch_design <- createDesignMatrix(ei, c("timepoint", "t"))

FDR_cutoff <- 0.05

pairwise_contrast_t6 <- createContrast(c(c(0, 0, 0, 1), rep(0,5)))
pairwise_contrast_dod <- createContrast(c(c(0, 0, 1, 0), rep(0,5)))
pairwise_contrast_c10 <- createContrast(c(c(0, 1, 0, 0), rep(0,5)))

# pairwise_contrast_dod_t6 <- createContrast(c(c(0, 0, -1, 1), rep(0,5)))
# pairwise_contrast_t6 <- createContrast(c(c(0, 0, 0, 1), rep(0,5), 1))
# pairwise_contrast_dod <- createContrast(c(c(0, 0, 1, 0), rep(0,5), 1))

#pairwise_contrast_t6 <- createContrast(c(0, 0, 1, 1))



da_c10 <- diffcyt(merged_daf,
                  design = design,
                  contrast = pairwise_contrast_c10,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "flo_merge",
                  verbose = T)


da_dod <- diffcyt(merged_daf,
                  design = design,
                  contrast = pairwise_contrast_dod,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "flo_merge",
                  verbose = T)

# vd2s <- grep("gamma", levels(cluster_ids(merged_daf, k="flo_merge")), fixed=T, value = T)
# merged_daf2 <- filterSCE(merged_daf,  cluster_id %in% vd2s, k="flo_merge")


da_t6 <- diffcyt(merged_daf,
                 design = design,
                 #contrast = contrast_t6,
                 contrast = pairwise_contrast_t6,
                 analysis_type = "DA",
                 method_DA = "diffcyt-DA-edgeR",
                 clustering_to_use = "flo_merge",
                 verbose = T)

table(rowData(da_c10$res)$p_adj < FDR_cutoff)
# FALSE
#    42
table(rowData(da_dod$res)$p_adj < FDR_cutoff)
# # FALSEFALSE  TRUE
# 34     8
table(rowData(da_t6$res)$p_adj < FDR_cutoff)
# # FALSEFALSE  TRUE
# 28     14


# with C+10 being the intercept term, c(-1, 1) contrasts can be made between baseline and dod/t6
# the significant clusters at dod are 10, 12, 14, 37, 3, 6, 43, 45, (n=8) p_adjust range from 2.5e-6 to 2.3e-2
# the significant clusters at t6 are 37, 15, 43, 10, 36, 12, 42, 14, 6, 8, 3, 29, 18, 19 (n=14), p+adjust range from 5.2e-52 to 4.4e-2

# these are exactly the same clusters, when baseline is the intercept term and we set it to 0 for the contrasts
# this is good news!!!

# using a design matrix with timepoint and batch reduces our dod clusters to 0 and the t6 clusters to 6 (37, 43, 36, 15, 42, 12)
# removing the cluster term returns 0 significant clusters at dod and 5 at t6 (37, 32, 36, 15, 42) so there's a small effect;
# all the mismatched clusters between those models are detected in the full design matrix
plotDiffHeatmap(merged_daf, da_c10, th = FDR_cutoff, normalize = TRUE, hm1 = F, top_n = 25)
plotDiffHeatmap(merged_daf, da_dod, th = FDR_cutoff, normalize = TRUE, hm1 = F, top_n = 10)
plotDiffHeatmap(merged_daf, da_t6, th = FDR_cutoff, normalize = TRUE, hm1 = F, top_n = 30)




# glmms #####

#this one works, don't change
da_formula1 <- createFormula(ei, cols_fixed = "timepoint",
                             cols_random = "volunteer")

# this one you're allowe to play with
da_formula2 <- createFormula(ei, cols_fixed = c("timepoint"),
                             cols_random = c("batch"))



glm_contrast_c10 <- createContrast(c(0, 1, 0, 0))
glm_contrast_dod <- createContrast(c(0, 0, 1, 0))
glm_contrast_t6 <- createContrast(c(0, 0, 0, 1))

glm_c10 <- diffcyt(merged_daf,
                      formula = da_formula2,
                      contrast = glm_contrast_c10,
                      analysis_type = "DA",
                      method_DA = "diffcyt-DA-GLMM",
                      clustering_to_use = "meta45",
                      verbose = T)

glm_dod <- diffcyt(merged_daf,
                      formula = da_formula2,
                      contrast = glm_contrast_dod,
                      analysis_type = "DA",
                      method_DA = "diffcyt-DA-GLMM",
                      clustering_to_use = "meta45",
                      verbose = T)

glm_t6 <- diffcyt(merged_daf,
                     contrast = glm_contrast_t6,
                     formula = da_formula2,
                     analysis_type = "DA",
                     method_DA = "diffcyt-DA-GLMM",
                     clustering_to_use = "meta45",
                     verbose = T)



table(rowData(glm_c10$res)$p_adj < FDR_cutoff)
# FALSE  TRUE 
# 18    27
table(rowData(glm_dod$res)$p_adj < FDR_cutoff)
# FALSE  TRUE 
# 9    36  
table(rowData(glm_t6$res)$p_adj < FDR_cutoff)
# FALSE  TRUE 
# 3    42


plotDiffHeatmap(merged_daf, glm_c10, th = FDR_cutoff, normalize = TRUE, hm1 = F, top_n = 25)
plotDiffHeatmap(merged_daf, glm_dod, th = FDR_cutoff, normalize = TRUE, hm1 = F, top_n = 25)
plotDiffHeatmap(merged_daf, glm_t6, th = FDR_cutoff, normalize = TRUE, hm1 = F, top_n = 25)



glm_c10_sig <- topTable(glm_c10, top_n = 5)$cluster_id
glm_dod_sig <- topTable(glm_dod, top_n = 18)$cluster_id
glm_t6_sig <- topTable(glm_t6, top_n = 24)$cluster_id


edger_dod_sig <- topTable(da_dod, top_n = 8)$cluster_id
edger_t6_sig <- topTable(da_t6, all = T, show_logFC = T)


#write.csv(edger_t6_sig, "/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/edger_t6_FC_padj.csv")

all(edger_dod_sig %in% glm_dod_sig) # TRUE
all(edger_t6_sig %in% glm_t6_sig) # TRUE


# BOXPLOTS OF CLUSTER COUNTS AND FREQUENCIES ####

all_cluster_counts <- diffcyt_boxplot(da_t6, merged_daf, counts=T, FDR=1, logFC = 0)
all_cluster_log_counts <- diffcyt_boxplot(da_t6, merged_daf, counts=T, FDR=1, logFC = 0)+scale_y_log10(name = "Number of Cells in Sample")
all_cluster_freqs <- diffcyt_boxplot(da_t6, merged_daf, counts=F, FDR=1, logFC = 0)+theme(axis.title = element_text(size=12))

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/all_clusters_counts.png",all_cluster_counts , height = 12, width=18)# works
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/all_clusters_freqs.png",all_cluster_freqs, height = 12, width=18) #works
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/all_clusters_log_counts.png",all_cluster_log_counts, height = 12, width=18)# works

write.csv(all_cluster_freqs$data,"~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_cluster_freqs.csv")

da_dod_freq_box <- diffcyt_boxplot(da_dod, merged_daf, FDR=0.05, logFC=1)
da_t6_freq_box <- diffcyt_boxplot(da_t6, merged_daf, FDR=0.05, logFC=1)

da_dod_count_box <- diffcyt_boxplot(da_dod, merged_daf, FDR=0.05, counts=T, logFC=1)
da_t6_count_box <- diffcyt_boxplot(da_t6, merged_daf, FDR=0.05, counts=T, logFC=1)

da_dod_log_count_box <- diffcyt_boxplot(da_dod, merged_daf, FDR=0.05, counts=T, logFC=1)+scale_y_log10()
da_t6_log_count_box <- diffcyt_boxplot(da_t6, merged_daf, FDR=0.05, counts=T, logFC=1)+scale_y_log10()


ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/diffcyt/edgeR/da_dod_freq_box005_log2fc.png", da_dod_freq_box, height = 7, width = 7)# works
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/diffcyt/edgeR/da_t6_freq_box005_log2fc.png", da_t6_freq_box, height = 7, width = 14)# works

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/diffcyt/edgeR/da_dod_count_box005_log2fc.png", da_dod_count_box, height = 7, width = 7)# works
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/diffcyt/edgeR/da_t6_count_box005_log2fc.png", da_t6_count_box, height = 7, width = 14)# works

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/diffcyt/edgeR/da_dod_log_count_box005_log2fc.png", da_dod_log_count_box, height = 7, width = 7)# works
ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/diffcyt/edgeR/da_t6_log_count_box005_log2fc.png", da_t6_log_count_box, height = 7, width = 14)# works



# stacked barcharts ####


t6_map_data <- da_t6_log_count_box$data

t6_barchart_data <- subset(t6_map_data, t6_map_data$timepoint=="T6")

sig_t6_clusters <- as.character(unique(t6_barchart_data$cluster_id))

write.csv(sig_t6_clusters, "/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/sig_t6_clusters.csv", row.names = F)

ggplot(t6_barchart_data, aes(x=cluster_id, y=frequency))+
  geom_boxplot(aes(fill=cluster_id))+
  geom_point(aes(shape=volunteer))+
  theme_minimal()



t6_barchart_data <- t6_barchart_data %>%
  group_by(volunteer) %>%
  mutate(., total_cd3=sum(frequency))





flo_merge_cd3_stacked_barchart <- ggplot(t6_barchart_data, aes(x=timepoint, y=frequency/100, group=volunteer))+
  geom_bar(stat="identity", position="stack", aes(fill=cluster_id))+
  geom_text(aes(y=(total_cd3/100)+0.01, label=paste0(round(total_cd3, digits = 1), "%", sep='')))+
  theme_minimal()+
  facet_wrap(~volunteer)+
  ggtitle("Significant Clusters at T6")+
  scale_y_continuous(name = "Percentage of CD3+ T cells\n\n", labels=scales::percent_format(accuracy = 1))+
  #ylim(0,25)+
  #geom_text(aes(label=cluster_id), position = position_stack(vjust = .5))+
  theme(#legend.position = "none",
    plot.title = element_text(hjust=0.5, size=15),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_text(hjust=0.5, size=12, face = "bold"),
    panel.grid.minor.y = element_blank())

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/diffcyt/edgeR/flo_merge_cd3_stacked_barchart.png", flo_merge_cd3_stacked_barchart, height=7, width=11)


# garbage code you'll likely never need ####
#heatmap of normalised significant cluster frequencies #

# t6_map_data$scaled_trans_norm <- scales::rescale(t6_map_data$trans_norm_freq,
#                                                  from=min(t6_map_data$trans_norm_freq, 0), to=c(-2,0)
#                                                  )


asinTransform <- function(x){asin(sqrt(x))}

t6_map_data <- da_t6_log_count_box$data

t6_map_data <- t6_map_data%>%
  group_by(cluster_id) %>%
  #group_by(volunteer) %>%
  mutate(trans_freq=asin(sqrt(frequency/100))) %>%
  mutate(max_freq=max(frequency)) %>%
  mutate(trans_norm_freq=scale(trans_freq, center = TRUE, scale = TRUE)) %>%
  ungroup()

#write.csv(t6_map_data, "/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/t6_map_data.csv")
all_t6_data <- all_cluster_counts$data
write.csv(all_t6_data, "/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_t6_data.csv")



t6_col_levels <- unique(t6_map_data[order(t6_map_data$timepoint, t6_map_data$volunteer),]$sample_id)

t6_map_data$sample_id <- factor(t6_map_data$sample_id, levels=t6_col_levels) 


#color bullshit (flo_berlin best) 


berlin <- colorspace::diverging_hcl(palette = "Berlin",15)
flo_berlin <- berlin[c(4:15)]


cute <- colorspace::diverging_hcl(palette="Blue-Yellow 2", 12)
yellow_half <- colorspace::diverging_hcl(palette="Blue-Yellow 2", 12)[6:12]
blue_half <- colorspace::diverging_hcl(palette="Blue-Yellow 3", 12)[1:5]

flo_blue_yellow <- c(blue_half, yellow_half)

cute <- rev(colorspace::sequential_hcl(palette="YlGnBu", 12))


flo_blue_red <-  c("#061B29", "#FFFFFF", "#411C19")




berlin <- colorspace::diverging_hcl(palette = "Berlin",16)

red_half <- rev(berlin[12:16])
blue_half <- rev(berlin[1:5])

yellow_half <- rev(c("#FFD400", "#FFDD3C", "#FFEA61", "#FFF192", "#FFFFB7"))


yellow_half <- rev(c("#FFC30B", "#FFD400", "#FFFFB7"))


#show_col(c(blue_half[3:5], yellow_half))

# flo diff heatmap #


  (da_edger_t6_heatmap <- ggplot(t6_map_data, aes(x=factor(sample_id, levels=t6_col_levels), y=cluster_id))+
  #geom_tile(aes(fill=viridis_rescaler(t6_map_data$trans_norm_freq)))+
  geom_tile(aes(fill=trans_norm_freq))+
  scale_fill_gradientn(#colours = colorspace::diverging_hcl(palette = "Berlin",8), 
                       #colours=rev(c("#FFC125", "#FFFFFF", "#4169E1")),
                       colours=c(blue_half[3:5], yellow_half),
                       #values = rescale(c(min(t6_map_data$trans_norm_freq), 0, max(t6_map_data$trans_norm_freq)), to=c(0,1)), 
                       values = rescale(c(min(t6_map_data$trans_norm_freq), 0, max(t6_map_data$trans_norm_freq)), to=c(0,1)), 
                       guide = guide_colourbar(nbin = 1000))+
  theme_minimal()+
  ggtitle("T6, FDR=0.05, logFC=1")+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.text = element_text(size=14),
        legend.position="right",
        legend.title = element_blank())
)

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/diffcyt/edgeR/da_edger_t6_heatmap.png", da_edger_t6_heatmap, height=7, width=11)





dod_map_data <- da_dod_log_count_box$data

dod_map_data <- dod_map_data%>%
  group_by(cluster_id) %>%
  #group_by(volunteer) %>%
  mutate(trans_freq=asin(sqrt(frequency/100))) %>%
  mutate(max_freq=max(frequency)) %>%
  mutate(trans_norm_freq=scale(trans_freq, center = TRUE, scale = TRUE))



dod_col_levels <- unique(dod_map_data[order(dod_map_data$timepoint, dod_map_data$volunteer),]$sample_id)

dod_map_data$sample_id <- factor(dod_map_data$sample_id, levels=dod_col_levels) 






(da_edger_dod_heatmap <- ggplot(dod_map_data, aes(x=factor(sample_id, levels=col_levels), y=cluster_id))+
    #geom_tile(aes(fill=viridis_rescaler(dod_map_data$trans_norm_freq)))+
    geom_tile(aes(fill=trans_norm_freq))+
    scale_fill_gradientn(#colours = colorspace::diverging_hcl(palette = "Berlin",8), 
      #colours=rev(c("#FFC125", "#FFFFFF", "#4169E1")),
      colours=c(blue_half[3:5], yellow_half),
      values = rescale(c(min(dod_map_data$trans_norm_freq), 0, max(dod_map_data$trans_norm_freq)), to=c(0,1)), 
      guide = guide_colourbar(nbin = 1000))+
    ggtitle("DoD, FDR=0.05, logFC=1")+
    theme_minimal()+
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1),
          axis.text = element_text(size=14),
          legend.position="right",
          legend.title = element_blank())
)

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/diffcyt/edgeR/da_edger_dod_heatmap.png", da_edger_dod_heatmap, height=7, width=11)




heat <-  grid::grid.grabExpr(
  plotClusterHeatmap(merged_daf, hm2 = "state_markers",
                   k = "meta45",
                   m = "flo_merge",
                   cluster_anno = FALSE,
                   draw_freqs = FALSE,
                   scale=TRUE, draw_dend = FALSE)
)


ggsave("heat.png", heat, width=20, height=15)


















viridis_rescaler <- function(x){
  ifelse(x<0,
         scales::rescale(x,
                         to=c(-max(x),0),
                         from=c(min(x, na.rm = TRUE), 0)
                         ), 
         x
  )
}
