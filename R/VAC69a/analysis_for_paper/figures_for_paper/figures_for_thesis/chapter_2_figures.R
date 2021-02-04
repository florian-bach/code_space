#boxplots of all cluster frequencies ####

merged_daf <- read_full("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")

ei <- metadata(merged_daf)$experiment_info

design <- createDesignMatrix(ei, c("timepoint", "volunteer"))

pairwise_contrast_t6 <- createContrast(c(c(0, 0, 0, 1), rep(0,5)))

da_t6 <- diffcyt(merged_daf,
                 design = design,
                 contrast = pairwise_contrast_t6,
                 analysis_type = "DA",
                 method_DA = "diffcyt-DA-edgeR",
                 clustering_to_use = "flo_merge",
                 verbose = T)

all_cluster_freqs <- diffcyt_boxplot(da_t6, merged_daf, counts=F, FDR=1, logFC = 0)+
  #scale_y_continuous(breaks=waiver(), n.breaks = 5)+
  # ylab("Cell Counts")+
  theme(axis.text.x = element_text(hjust=1, angle=30))

all_cluster_freqs$data$cluster_id <- gsub("gamma delta", "γδ", all_cluster_freqs$data$cluster_id, fixed=T)

# ggsave("/home/flobuntu/PhD/cytof/vac69a/final_figures_for_paper/supp_all_clusters_freqs.pdf", all_cluster_freqs, height = 12, width=18)# works
# ggsave("/home/flobuntu/PhD/figures_for_thesis/chapter_2/2_all_clusters_log_counts.png", all_cluster_freqs, height = 12 , width = 11)# works
ggsave("/home/flobuntu/PhD/figures_for_thesis/chapter_2/2_all_clusters_freqs.png", all_cluster_freqs, height = 12 , width = 11)# works





# lymph counts ####

library(ggplot2)
library(dplyr)

setwd("~/PhD/clinical_data/vac69a/")



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




data <- read.csv("lymph_counts.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(data)[2:3] <- c("B Cells", "T Cells")
data <- tidyr::gather(data, cell_type, count, colnames(data)[2:3])  
data$count <- data$count*1000  
#data$cell_typef <- factor(data$cell_type, levels=c("B cells", "T cells"))
data <- subset(data, cell_type=="T Cells")

data$Volunteer <- gsub("v", "v0", data$Volunteer, fixed = T)

data$Timepoint <- data$Timepoint %>%
  gsub("C-1", "Baseline", .) %>%
  gsub("T+6", "T6", ., fixed = T) %>%
  gsub("DoD", "Diagnosis", .)
  

t_cell_counts <- ggplot(data, aes(x=factor(Timepoint, levels=c("Baseline", "Diagnosis", "T6")), y=count, group=Volunteer))+
  geom_line(aes(color=Volunteer), size=0.9)+
  geom_point(aes(color=Volunteer), fill="white", stroke=1, shape=21, size=0.9)+
  scale_colour_manual(values=volunteer_palette)+
  ylab(expression(paste("T Cells / ", mu*"L")))+
  fig1_theme+
  theme_minimal()+
  guides(color=guide_legend(nrow = 1))+
  ylim(0,1750)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "white", color = "white"))

ggsave("~/PhD/figures_for_thesis/chapter_2/2_t_cell_counts.png", t_cell_counts, height = 3, width = 4.2, units = "in")  



#ds limma ####



daf <- read_full("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")

coarse_table <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/merging_tables/most_coarse_merge.csv", header=T, stringsAsFactors = F)



#get rid of spaces at beginning of string

coarse_table$new_cluster <- factor(coarse_table$new_cluster)

merged_daf<- mergeClusters(daf, k = "mod_meta45", table = coarse_table, id = "coarse_merge")


plotClusterHeatmap(merged_daf, hm2=NULL,
                   k = "mod_meta45",
                   m = "coarse_merge",
                   cluster_anno = FALSE,
                   draw_freqs = TRUE,
                   scale = TRUE
)



ei <- metadata(merged_daf)$experiment_info

design <- createDesignMatrix(ei, c("timepoint", "batch"))

levels(ei$timepoint)

###If a design matrix has been used, the entries of contrast correspond to the columns of the design
#matrix and the length of contrast equals the number of columns in the design matrix. If a model formula
#has been used, the entries correspond to the levels of the fixed effect terms;
#and the length equals the number of levels of the fixed effect terms.

FDR_cutoff <- 0.05

# edgeR models with all coefficient for volunteer####
# pairwise_contrast_t6 <- createContrast(c(c(0, 0, 0, 1), rep(0,5)))
# pairwise_contrast_dod <- createContrast(c(c(0, 0, 1, 0), rep(0,5)))
# pairwise_contrast_c10 <- createContrast(c(c(0, 1, 0, 0), rep(0,5)))

# copntrasts for models excluding volunteer

pairwise_contrast_t6 <- createContrast(c(c(0, 0, 0, 1), rep(0,1)))
pairwise_contrast_dod <- createContrast(c(c(0, 0, 1, 0), rep(0,1)))
pairwise_contrast_c10 <- createContrast(c(c(0, 1, 0, 0), rep(0,1)))

# how to test only activation markers
states <- names(marker_classes(merged_daf))
states <- states[-c(1:21,
                    match(c("CCR7", "CD127", "CD20", "CXCR5", "CD103", "TCRgd", "TIM3", "IntegrinB7", "CD56", "CD3", "CD49d", "CD45RA", "CD4", "Vd2", "Va72", "CD161", "FoxP3", "CD45RO"), states),
                    53, 54, 56:length(states)
)
]

# how to test ALLL markers
states <- names(marker_classes(merged_daf))
states <- states[-c(1:21,53, 54, 56:length(states),
                    match(c("CD20", "TIM3", "KG", "CD3"), states))]



logic <- names(marker_classes(merged_daf)) %in% states

metadata(merged_daf)$id_state_markers <- states

# ds_c10 <- diffcyt(merged_daf,                                            
#                   design = design, contrast = pairwise_contrast_c10,                    
#                   analysis_type = "DS", method_DS = "diffcyt-DS-limma",       
#                   clustering_to_use = "coarse_merge", verbose = TRUE)               


ds_dod <- diffcyt(merged_daf,                                            
                  design = design, contrast = pairwise_contrast_dod,                    
                  analysis_type = "DS", method_DS = "diffcyt-DS-limma",       
                  clustering_to_use = "coarse_merge", verbose = TRUE)               


ds_t6 <- diffcyt(merged_daf, markers_to_test = rep(TRUE, 28),                         
                 design = design, contrast = pairwise_contrast_t6,                    
                 analysis_type = "DS",  method_DS =  "diffcyt-DS-limma",         
                 clustering_to_use = "coarse_merge", verbose = TRUE)



topTable(ds_dod, top_n = 5, order_by = "cluster_id",
         show_meds = TRUE, format_vals = TRUE, digits = 3)
topTable(ds_t6, top_n = 55, order_by = "cluster_id",
         show_meds = TRUE, format_vals = TRUE, digits = 3)


# res_DA_c10 <- topTable(ds_c10, all = TRUE, show_logFC = T)
# table(res_DA_c10$p_adj <= 0.05)

res_DA_dod <- topTable(ds_dod, all = TRUE, show_logFC = T)
table(res_DA_dod$p_adj <= 0.05)

res_DA_t6 <- topTable(ds_t6, all = TRUE, show_logFC = T)
table(res_DA_t6$p_adj <= 0.05)

sigs_t6 <- subset(res_DA_t6, p_adj<=0.05)
sigs_dod <- subset(res_DA_dod, p_adj<=0.05)

spider_sigs_t6 <- data.frame(sigs_t6$cluster_id, sigs_t6$marker_id, sigs_t6$logFC)
#spider_sigs_t6 <- data.frame(sigs_dod$cluster_id, sigs_dod$marker_id, sigs_dod$logFC)

#spider_sigs_t6$direction <- ifelse(sigs_t6$logFC>=0, "up", "down")

#trimmed_spider_sigs_t6 <- subset(spider_sigs_t6, 2^spider_sigs_t6$sigs_t6.logFC>= 1.1 |  2^spider_sigs_t6$sigs_t6.logFC<=0.9)

#trimmed_spider_sigs_t6$sigs_t6.cluster_id <- gsub("gamma delta", expression(alpha) , trimmed_spider_sigs_t6$sigs_t6.cluster_id)

(ds_limma_t6_heatmap <- ggplot(data=spider_sigs_t6, aes(x= sigs_t6.cluster_id, y=sigs_t6.marker_id, label=round(2^sigs_t6.logFC, digits = 2)))+
    #geom_tile(aes(fill=rescale(sigs_t6.logFC, to=c(-5, 5))))+
    geom_tile(aes(fill=round(2^sigs_t6.logFC, digits = 2)))+
    geom_text(parse = TRUE, aes(color=2^spider_sigs_t6$sigs_t6.logFC>0.8))+
    scale_color_manual(values = c("TRUE"="black", "FALSE"="white"), guide=FALSE)+
    scale_fill_gradientn("Fold Change", values = scales::rescale(c(2^min(spider_sigs_t6$sigs_t6.logFC), 1, 2^max(spider_sigs_t6$sigs_t6.logFC)), to=c(0,1)), colors = c("#0859C6","white","#FFA500"))+
    scale_x_discrete(labels=c("gamma delta"=expression(paste(gamma, delta))))+
    theme_minimal()+
    ggtitle("Fold Change in Intensity of Marker Expression\nper Cluster at T6 relative to Baseline")+
    theme(axis.title = element_blank(),
          axis.text = element_text(size=14, color = "black"),
          axis.text.x = element_text(hjust=1, angle = 30),
          legend.position = "right",
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
)

ggsave("/home/flobuntu/PhD/figures_for_thesis/chapter_2/ds_limma_01.png", ds_limma_t6_heatmap, width = 7, height = 7)







shplit_cells <- function(x, by) {
  stopifnot(is.character(by), by %in% colnames(colData(x)))
  cd <- data.frame(colData(x))
  cd$cluster_id <- cluster_ids(x, "coarse_merge")
  dt <- data.table::data.table(cd, i = seq_len(ncol(x)))
  dt_split <- split(dt, by = by, sorted = TRUE, flatten = FALSE)
  map_depth(dt_split, length(by), "i")
}


ahgg <- function(x, by, fun = c("median", "mean", "sum")) {
  fun <- switch(match.arg(fun), 
                median = rowMedians, mean = rowMeans, sum = rowSums)
  cs <- shplit_cells(x, by)
  pb <- map_depth(cs, -1, function(i) {
    if (length(i) == 0) return(numeric(nrow(x)))
    fun(assay(x, "exprs")[, i, drop = FALSE])
  })
  map_depth(pb, -2, function(u) as.matrix(data.frame(
    u, row.names = rownames(x), check.names = FALSE)))
}



ms <- ahgg(merged_daf[refined_markers$refined_markers, ], by = c("cluster_id", "sample_id"))
ms <- lapply(ms, reshape2::melt, varnames = c("antigen", "sample_id"))
ms <- bind_rows(ms, .id = "cluster_id")

ms$timepoint <- substr(ms$sample_id, 5, nchar(as.character(ms$sample_id)))

scaled_ms <- ms %>%
  group_by(cluster_id, antigen) %>%
  mutate("scaled_value"=scale(value)) %>%
  ungroup()


#scaled_ms <- subset(scaled_ms, scaled_ms$antigen %in% states[-(c(3, 11:14))])

sig_ds_contrasts <- subset(rowData(ds_t6$res), rowData(ds_t6$res)$p_adj<0.05)

sig_ds_contrasts <- paste(sig_ds_contrasts$cluster_id, " (", sig_ds_contrasts$marker_id, ")", sep="")



scaled_ms$volunteer <- substr(scaled_ms$sample_id, 0, 3)


wide_scaled_ms <- scaled_ms %>%
  #group_by(cluster_id, antigen) %>%
  #mutate("mean_scaled"=mean(scaled_value)) %>%
  #ungroup()%>%
  mutate("contrast"=paste(cluster_id, " (", antigen, ")", sep="")) %>%
  filter(contrast %in% sig_ds_contrasts) %>%
  select(scaled_value, contrast, sample_id, antigen) %>%
  pivot_wider(names_from = sample_id, values_from = scaled_value)




wide_scaled_ms <- data.frame(wide_scaled_ms)
rownames(wide_scaled_ms) <- wide_scaled_ms$contrast



wide_scaled_ms$contrast <- NULL
num_wide_scaled_ms <- as.matrix(wide_scaled_ms[, 2:ncol(wide_scaled_ms)])

rownames(num_wide_scaled_ms) <- gsub("gamma delta", "γδ", rownames(num_wide_scaled_ms), fixed=T)


library(ComplexHeatmap)

inferno <- colorspace::sequential_hcl("inferno", n=8)
col_fun_ds_limma <- circlize::colorRamp2(c(min(num_wide_scaled_ms), 0, max(num_wide_scaled_ms)), c("#0859C6", "black", "#FFA500"))

limma_col_split <- factor(substr(colnames(wide_scaled_ms)[2:ncol(wide_scaled_ms)], 5, 20))
limma_row_split_lineage <- factor(substr(rownames(wide_scaled_ms),0, 4))

limma_row_split_lineage <- limma_row_split_lineage %>%
  gsub("Doub", "DN", .) %>%
  gsub(" ", "", ., fixed = TRUE) %>%
  gsub("gamm", "γδ", .)

limma_row_split_antigen <- factor(wide_scaled_ms$antigen)

double_whammy <- cbind(limma_row_split_lineage, limma_row_split_antigen)


median_cluster_heat <- Heatmap(matrix = num_wide_scaled_ms,
                               cluster_rows = FALSE,
                               name = "Scaled Marker Expression",
                               cluster_columns = FALSE,
                               row_names_side = "left",
                               col = col_fun_ds_limma,
                               column_names_gp = gpar(fontsize = 11),
                               #split = factor(wide_scaled_ms$antigen),
                               column_split = limma_col_split,
                               split = limma_row_split_lineage,
                               rect_gp = gpar(col = "white"),
                               show_heatmap_legend = TRUE,
                               column_names_rot = 45,
                               heatmap_legend_param = list(legend_position = "top",
                                                           col=col_fun_ds_limma,
                                                           title = "Normalised Marker Expression",
                                                           legend_direction = "horizontal",
                                                           title_position = "topcenter",
                                                           legend_width = unit(6.2, "cm"),
                                                           border = FALSE)
                               #width = unit(16, "cm"),
                               #height = unit(16*9/28, "cm")
)

png("~/PhD/figures_for_thesis/chapter_2/all_markers_ds_limma.png", width = 8, height=10, units="in", res=400)
draw(median_cluster_heat, heatmap_legend_side = "bottom")
dev.off()


# activation stacked barhcart plus pie ####


library(vac69a.cytof)
library(dplyr)
summary <- data.table::fread("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/activation_barchart_data")
summary$volunteer <- gsub("V", "v", summary$volunteer, fixed=T)

summary$timepoint <- gsub("DoD", "Diagnosis", summary$timepoint, fixed=T)


lineage_palette <- c("#AA3377", "#EE6677", "#4477AA", "#66CCEE", "#228833", "#FFA500", "#BBBBBB")
names(lineage_palette) <-c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "Resting")



activation_stacked_barchart <- ggplot(summary, aes(x=volunteer, y=frequency/100, fill=lineage))+
  geom_bar(stat="identity", position="stack")+
  #geom_text(aes(y=(total_cd3/100)+0.01, label=paste0(round(total_cd3, digits = 1), "%", sep='')))+
  theme_minimal()+
  facet_wrap(~timepoint, ncol=4)+
  #ggtitle("Overall T cell activation")+
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


# lineage pies #

#get cluster freqs and cluster lineage freqs
barchart_data <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/cluster_and_lineage_stats_for_barcharts.csv", header=T, stringsAsFactors = FALSE)

#keep only activated clusters at T6
activated_clusters <- grep("activated", unique(barchart_data$cluster_id), value = T)
barchart_data <- subset(barchart_data, barchart_data$cluster_id %in% activated_clusters)
barchart_data <- subset(barchart_data, barchart_data$timepoint=="T6")

# get average lineage activation
barchart_data <- barchart_data %>%
  group_by(timepoint, lineage) %>%
  mutate("lin_activation" = sum(lin_freq)/6) %>%
  ungroup()


# construct ggplotable dataframe; now this table contains each average for all volunteers, so get rid of duplicates, then append
# the df with 1 -the same values so that the activated and nonactivated percentages are in the same column

barchart_plot_data <- data.frame("lineage" = barchart_data$lineage, "lin_activation" = barchart_data$lin_activation)
barchart_plot_data <- barchart_plot_data[!duplicated(barchart_plot_data), ]
barchart_plot_data <- rbind(barchart_plot_data, data.frame("lineage"=barchart_plot_data$lineage, "lin_activation"=1-barchart_plot_data$lin_activation))

# add activated/resting column, make it so that the greek letters render in ggplot
barchart_plot_data$state <- c(as.character(barchart_plot_data$lineage[1:6]), rep("Resting", 6))
barchart_plot_data$state <- gsub("gd", expression(paste(gamma, delta)), barchart_plot_data$state)
barchart_plot_data$lineage <- gsub("gd", expression(paste(gamma, delta)), barchart_plot_data$lineage)


# we have to rename the colour palette because otherwise the strip.text won't parse properly
lineage_palette <- c("#AA3377", "#EE6677", "#4477AA", "#66CCEE", "#228833", "#FFA500", "#BBBBBB")
names(lineage_palette) <-c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "Resting")

pie_lineage_palette <- lineage_palette
names(pie_lineage_palette)<- c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "Resting")

#order lineage and state columns for plotting
barchart_plot_data$lineage <- factor(barchart_plot_data$lineage, levels=c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "Resting"))
barchart_plot_data$state <- factor(barchart_plot_data$state, levels=rev(c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "Resting")))

lineage_activation_pies <- ggplot(barchart_plot_data)+
  geom_bar(stat="identity", position = "stack", width=1, aes(x="", y=lin_activation, fill=state))+
  theme_minimal()+
  coord_polar(theta = "y", start = 0)+
  facet_wrap(~lineage, labeller = label_parsed, ncol=2, strip.position = "left")+
  scale_fill_manual(values=pie_lineage_palette)+
  guides()+
  ggtitle("Average Activation in Major\nT cell Lineages at T6")+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"),
        plot.margin = unit(c(0,0,0,0), "mm"),
        plot.title = element_text(hjust = 0.5, size=7),
        strip.text.y.left = element_text(hjust=0.5, size=7, face = "bold", angle = 0),
        legend.position = "none")



#put the bar and pie charts together with cowlplot

panel_f <- cowplot::plot_grid(activation_stacked_barchart, lineage_activation_pies, ncol = 2, rel_widths = c(4,1), rel_heights = c(1,1), axis = "bt", align = "hv")

ggsave("~/PhD/figures_for_thesis/chapter_2/activation_stack_and_pie.png", panel_f, height=2, width=8)
``



# individual variation figures ####


library(ggplot2)
library(dplyr)
library(tidyr)

volunteer_colours <- list("v02" = "#FB9A99",
                          "v03" = "#E31A1C",
                          "v05" = "#A6CEE3",
                          "v06" = "#1F78B4",
                          "v07" = "#F0E442",
                          "v09" = "#E69F00")


volunteer_palette <- unlist(unname(volunteer_colours))
names(volunteer_palette) <- names(volunteer_colours)



time_col=colorspace::sequential_hcl(5, palette = "Purple Yellow")

time_palette <- c("Baseline"=time_col[4], "C10"=time_col[3], "Diagnosis"=time_col[2], "T6"=time_col[1])

# CyTOF MDS ####

# cluster frequency

cytof_data <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_t6_data.csv", header = T, stringsAsFactors = F)
#cytof_data <- subset(cytof_data, cytof_data$timepoint!="C10")

#cytof_data$trans_freq=scale(asin(sqrt(cytof_data$frequency/100)), center = TRUE, scale = TRUE)

cytof_data$trans_freq=asin(sqrt(cytof_data$frequency/100))


# wide_cytof <- data.frame(cytof_data %>%
#                            select(cluster_id, sample_id, frequency) %>%
#                            pivot_wider(names_from = sample_id, values_from = frequency))

wide_cytof <- data.frame(cytof_data %>%
                           select(cluster_id, sample_id, frequency) %>%
                           spread(sample_id, frequency)
)

rownames(wide_cytof) <- gsub("ï", "i", wide_cytof$cluster_id)
wide_cytof <- as.matrix(select(wide_cytof, -cluster_id))
class(wide_cytof) <- "double"

activated_cytof <- subset(wide_cytof, grepl("activated", rownames(wide_cytof)))

cytof_mds <- data.frame(cmdscale(robCompositions::aDist(t(wide_cytof))))

# only activated clusters- individual variation is reduced overall, but still doesn't converge much thorugh time,
# v05 & v09 cluster sperate from the other volutneers

activated_cytof_mds <- data.frame(cmdscale(robCompositions::aDist(t(activated_cytof))))


colnames(cytof_mds) <- c("MDS1", "MDS2")

cytof_mds$Sample_ID <- rownames(cytof_mds)
cytof_mds$Timepoint <- substr(cytof_mds$Sample_ID, 5, nchar(cytof_mds$Sample_ID))
cytof_mds$Volunteer <- substr(cytof_mds$Sample_ID, 1, 3)


colnames(activated_cytof_mds) <- c("MDS1", "MDS2")

activated_cytof_mds$Sample_ID <- rownames(activated_cytof_mds)
activated_cytof_mds$Timepoint <- substr(activated_cytof_mds$Sample_ID, 5, nchar(activated_cytof_mds$Sample_ID))
activated_cytof_mds$Volunteer <- substr(activated_cytof_mds$Sample_ID, 1, 3)




aitchison_cytof <- ggplot(cytof_mds, aes(x=MDS1, y=MDS2))+
  geom_point(aes(shape=Timepoint, color=Volunteer, fill=Volunteer))+
  #ggrepel::geom_label_repel(aes_string(label = "sample_id"), show.legend = FALSE)+ 
  theme_minimal()+
  ylab("")+
  ggtitle("Cluster Frequency Composition")+
  #scale_shape_manual(values = c("Baseline"=21, "C10"=24, "Diagnosis"=22, "T6"=3))+
  scale_color_manual(values = volunteer_palette)+
  scale_fill_manual(values = volunteer_palette)+
  guides(color=guide_legend(title="Volunteer",override.aes = list(size = 1)),
         shape=guide_legend(title="Timepoint", override.aes = list(size = 1)),
         fill=guide_none())+
  theme(legend.position = "none", plot.title = element_text(size=10, hjust=0.5))
        



#median marker expression

merged_daf <- vac69a.cytof::read_full("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")


limma_markers <- c("CD38','ICOS','CD27','Tbet','Ki67','FoxP3','CD127','PD1','BCL2','GZB','CTLA4','CD25','CD95','HLADR','CD28")

cs_by_s <- split(seq_len(ncol(merged_daf)), merged_daf$sample_id)
es <- as.matrix(SummarizedExperiment::assay(merged_daf, "exprs"))
ms <- vapply(cs_by_s, function(cs) Biobase::rowMedians(es[, cs, drop = FALSE]), 
             numeric(nrow(merged_daf)))
rownames(ms) <- rownames(merged_daf)
# state markers, type markers or both
ms <- subset(ms, rownames(ms)%in%climma_markers)

mds <- limma::plotMDS(ms, plot = FALSE)
df <- data.frame(MDS1 = mds$x, MDS2 = mds$y)
md <- S4Vectors::metadata(merged_daf)$experiment_info
m <- match(rownames(df), md$sample_id)
df <- data.frame(df, md[m, ])

#df2 <- filter(df, timepoint %in% c("Baseline", "DoD"))

state_markers_mds <- ggplot(df, aes(x=MDS1, y=MDS2))+
  geom_point(aes(shape=timepoint, color=volunteer, fill=volunteer))+
  #ggrepel::geom_label_repel(aes_string(label = "sample_id"), show.legend = FALSE)+ 
  theme_minimal()+
  xlab("")+
  ggtitle("Cell State Marker Expression")+
  #scale_shape_manual(values = c("Baseline"=21, "C10"=24, "Diagnosis"=22, "T6"=3))+
  scale_color_manual(values = volunteer_palette)+
  scale_fill_manual(values = volunteer_palette)+
  guides(color=guide_legend(title="Volunteer",override.aes = list(size = 1)),
         shape=guide_legend(title="Timepoint", override.aes = list(size = 1)),
         fill=guide_none())+
  theme(legend.title = element_text(size=8),
        legend.text = element_text(size=6),
        plot.title = element_text(size=10, hjust=0.5))

indie_var_lgd <- cowplot::get_legend(state_markers_mds)

state_markers_mds <- state_markers_mds+theme(legend.position = "none")

setwd("~/PhD/plasma/vac69a/")


volunteer_colours <- list("v02" = "#FB9A99",
                          "v03" = "#E31A1C",
                          "v05" = "#A6CEE3",
                          "v06" = "#1F78B4",
                          "v07" = "#F0E442",
                          "v09" = "#E69F00")


volunteer_palette <- unlist(unname(volunteer_colours))
names(volunteer_palette) <- names(volunteer_colours)

long_data <- read.csv("fancy_greek_long_data_sans9.csv", header = TRUE, stringsAsFactors = FALSE)
sig_analytes <- scan("~/PhD/plasma/vac69a/analytes_sorted_by_padj.txt", what="", skip = 1)

sig_analytes[1] <- "IFNγ"


#filter for top 12 analytes
long_data <- subset(long_data, long_data$analyte %in% sig_analytes[1:24])

data2 <- spread(long_data, analyte, concentration)
# split_data <- split(data2, data2$Volunteer)


data3 <- data2
data3[,3:ncol(data3)] <- log10(data3[,3:ncol(data3)])


big_pca <-  prcomp(data3[,3:ncol(data3)], center = T)
big_pca2 <- cbind(data3[, 1:2], big_pca$x)


plasma_pca <- ggplot(big_pca2, aes(x=PC1, y=PC2))+
  geom_point(aes(shape=timepoint, color=Volunteer, fill=Volunteer))+
  #xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  #ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  #ggrepel::geom_label_repel(aes_string(label = "sample_id"), show.legend = FALSE)+ 
  theme_minimal()+
  xlab("")+
  ylab("")+
  ggtitle("Plasma Marker PCA")+
  #scale_shape_manual(values = c("Baseline"=21, "C10"=24, "Diagnosis"=22, "T6"=3))+
  scale_color_manual(values = volunteer_palette)+
  scale_fill_manual(values = volunteer_palette)+
  # guides(color=guide_legend(title="Volunteer",override.aes = list(size = 1)),
  #        shape=guide_legend(title="Timepoint", override.aes = list(size = 1)),
  #        fill=guide_none())+
  theme(legend.position = "none", plot.title = element_text(size=10, hjust=0.5))


indie_var_plot <- cowplot::plot_grid(state_markers_mds, aitchison_cytof, plasma_pca, indie_var_lgd, nrow = 1, rel_widths = c(3,3,3,1))

ggsave("~/PhD/figures_for_thesis/chapter_2/indie_var_plot.pdf", indie_var_plot, height = 3.5, width=8)

