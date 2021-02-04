library(dplyr)
library(tidyr)
library(ggplot2)


`%notin%` <- Negate(`%in%`)
#making the cytof data for correlations ####
cytof_data <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_t6_data.csv", header = T, stringsAsFactors = F)
cytof_data <- subset(cytof_data, cytof_data$timepoint!="C10")

cytof_data$trans_freq=scale(asin(sqrt(cytof_data$frequency/100)), center = TRUE, scale = TRUE)

wide_cytof <- data.frame(cytof_data %>%
                           select(cluster_id, sample_id, frequency) %>%
                           pivot_wider(names_from = sample_id, values_from = frequency))

rownames(wide_cytof) <- wide_cytof$cluster_id
wide_cytof$cluster_id <- NULL
wide_cytof_sans_v09 <- subset(wide_cytof, select=!grepl("v09", colnames(wide_cytof)))
#wide_cytof_sans_v09 <- subset(wide_cytof_sans_v09, select=!grepl("C10", colnames(wide_cytof)))



#making clinical chemistry data ####
biochem_data <- read.csv("~/PhD/clinical_data/vac69a/biochem.csv")

colnames(biochem_data)[1] <- "Volunteer"

biochem_data$Volunteer <- paste("v", substr(biochem_data$Volunteer, 6, 7), sep='')

short <- filter(biochem_data, timepoint %in% c("_C_1", "_T6", "_EP"))
short2 <- select(short,  Volunteer, timepoint, sodium, potassium, creatinine, bilirubin, alt, alkphos, albumin, ast, ggt)

short2$timepoint <- as.character(short2$timepoint)

short2$timepoint[short2$timepoint=="_C_1"] <- "Baseline"
short2$timepoint[short2$timepoint=="_T6"] <- "T6"
short2$timepoint[short2$timepoint=="_EP"] <- "Diagnosis"

short2 <- filter(short2, timepoint %in% c("Baseline", "T6", "Diagnosis"))


mofa_biochem <- t(short2)
colnames(mofa_biochem) <- paste(mofa_biochem[1,], mofa_biochem[2,], sep="_") 
mofa_biochem <- mofa_biochem[3:nrow(mofa_biochem),]

#mofa_biochem <- as.matrix(mofa_biochem[,colnames(mofa_biochem)%in%colnames(wide_cytof)]) # <3
class(mofa_biochem) <- "double"

biochem_data <- log2(mofa_biochem[1:7,])
#biochem_data <- ifelse(is.infinite(biochem_data), NA, biochem_data)
biochem_data[5,1] <- 3.7
biochem_data <- biochem_data[c(3,5:7), ]

alt_timecourse <- mofa_biochem[5,order(colnames(mofa_biochem))]
alt_timecourse[7] <- 13

# making plasma anlayte data for correlations ####


dataplasma <- read.csv("~/PhD/plasma/vac69a/Vivax_plasma_analytes2_no_inequalities.csv")
dataplasma <- subset(dataplasma, dataplasma$Volunteer!="v009")
t_dataplasma <- t(dataplasma)#

colnames(t_dataplasma) <- paste(t_dataplasma[1,], t_dataplasma[2,], sep='_')

dataplasma <- data.frame(t_dataplasma[3:nrow(t_dataplasma),])

colnames(dataplasma) <- gsub("C.1", "Baseline", colnames(dataplasma))
colnames(dataplasma) <- gsub("v00", "v0", colnames(dataplasma))
colnames(dataplasma) <- gsub(".", "", colnames(dataplasma), fixed=T)
colnames(dataplasma) <- gsub("DoD", "Diagnosis", colnames(dataplasma))

#convert to numerical matri
plasma_data <- as.matrix(dataplasma[, colnames(dataplasma) %in% colnames(wide_cytof)]) # <3

changing_analytes <- sig_analytes <- scan("~/PhD/plasma/vac69a/analytes_sorted_by_padj.txt", what="", skip = 1)
changing_analytes <- changing_analytes[1:12]

rownames(plasma_data) <- gsub("pg.ml.", "", rownames(plasma_data), fixed = T)
rownames(plasma_data) <- gsub(".", "",rownames (plasma_data), fixed=T)

plasma_data <- subset(plasma_data, rownames(plasma_data) %in% changing_analytes)

class(plasma_data) <- "double" # <3

#put it together, names() makes view names
#log transform plasma data

log_plasma_data <- log2(plasma_data)

#this "alt_timecourse" variable was made ~50 lines below- don't ask...
log_plasma_data <- rbind(log_plasma_data, alt=log2(alt_timecourse[!grepl("v09", names(alt_timecourse), fixed=T)]))

biochem_data_sans_v9 <- subset(biochem_data, select=!grepl("v09", colnames(biochem_data)))


#plasma distance traveled####
plasma_mds <- limma::plotMDS(log_plasma_data, plot = F)

plasma_df <- data.frame(MDS1 = plasma_mds$x, MDS2 = plasma_mds$y)

plasma_df$Sample_ID <- rownames(plasma_df)
plasma_df$Timepoint <- substr(plasma_df$Sample_ID, 5, nchar(plasma_df$Sample_ID))
plasma_df$Volunteer <- substr(plasma_df$Sample_ID, 1, 3)

volunteer_colours <- list("v02" = "#FB9A99",
                          "v03" = "#E31A1C",
                          "v05" = "#A6CEE3",
                          "v06" = "#1F78B4",
                          "v07" = "#F0E442",
                          "v09" = "#E69F00")


volunteer_palette <- unlist(unname(volunteer_colours))
names(volunteer_palette) <- names(volunteer_colours)


plasma_mds_plot <- ggplot(plasma_df, aes(x=MDS1, y=MDS2, color=Volunteer))+
  geom_point(aes(shape=Timepoint))+
  #ggrepel::geom_label_repel(aes_string(label = "sample_id"), show.legend = FALSE)+ 
  theme_minimal()+
  scale_color_manual(values = volunteer_palette)+
  theme()


ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/plasma_mds.pdf", plasma_mds_plot, height=4, width=4)

arrow_pca <- subset(plasma_df, plasma_df$Timepoint %in% c("Baseline", "Diagnosis"))

wide_arrow_data <- arrow_pca[, c(1,2,4,5)]
wide_arrow_data <- pivot_wider(wide_arrow_data, names_from = Timepoint, values_from = c(MDS1, MDS2))



arrow_pca_plot <- ggplot(arrow_pca, aes(x=MDS1, y=MDS2))+
  #geom_point(aes(color=Volunteer))+
  geom_segment(data=wide_arrow_data,
               aes(x= MDS1_Baseline, xend=MDS1_Diagnosis,
                   y= MDS2_Baseline, yend=MDS2_Diagnosis,
                   color=Volunteer), arrow =arrow(length = unit(0.2, "cm")))+
  scale_color_manual(values=volunteer_palette)+
  xlab("PC1")+
  ylab("PC2")+
  #xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  #ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  theme(axis.text = element_text(size=10),
    #panel.border = element_rect(color="black", fill=NA),
    axis.title = element_text(size=12),
    legend.position = "none",
    plot.title = element_text(size=14, hjust=0.5)
  )

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/arrow_pca_plot.pdf", arrow_pca_plot, height=3, width =3)


plasma_df <- data.frame(MDS1 = plasma_mds$x, MDS2 = plasma_mds$y)
plasma_df$Sample_ID <- rownames(plasma_df)
plasma_df$Timepoint <- substr(plasma_df$Sample_ID, 5, nchar(plasma_df$Sample_ID))
plasma_df$Volunteer <- substr(plasma_df$Sample_ID, 1, 3)



distance_traveled <- filter(plasma_df,Timepoint%in%c("Baseline", "Diagnosis"))
distance_traveled_dfs <- split(distance_traveled, distance_traveled$Volunteer)

distance_traveled_mds1 <- lapply(distance_traveled_dfs, function(x) x$sum_sq_MDS2=(x$MDS2[2]-x$MDS2[1])^2)

distance_traveled_mds2 <- lapply(distance_traveled_dfs, function(x) x$sum_sq_MDS1=(x$MDS1[2]-x$MDS1[1])^2)

distance_traveled <- lapply(names(distance_traveled_mds2), function(x)
  data.frame(Volunteer=x, "MDS1"=as.numeric(distance_traveled_mds1[[x]]), "MDS2"=as.numeric(distance_traveled_mds2[[x]])))

distance_frame <- do.call(rbind, distance_traveled)

distance_frame$distance <- sqrt(apply(distance_frame[,2:3], 1, sum))

pca_distance <- distance_frame$distance
names(pca_distance) <- c("v02", "v03", "v05", "v06", "v07")


# put it all together, make big correlation matrix ####
big_table <- rbind(wide_cytof_sans_v09, log_plasma_data, biochem_data_sans_v9[c(1,3,4),])

# DoD plasma change; log2FC
plasma_dod_fc <- 2^subset(log_plasma_data, select=grepl("Diagnosis", colnames(log_plasma_data)))/2^subset(log_plasma_data, select=grepl("Baseline", colnames(log_plasma_data)))
colnames(plasma_dod_fc) <- substr(colnames(plasma_dod_fc), 1,3)
plasma_dod_fc <- subset(log2(plasma_dod_fc), rownames(plasma_dod_fc) %notin% c("alt", "IL18", "Ang2"))

#t6 plasma change
plasma_t6_fc <- 2^subset(log_plasma_data, select=grepl("T6", colnames(log_plasma_data)))/2^subset(log_plasma_data, select=grepl("Baseline", colnames(log_plasma_data)))
colnames(plasma_t6_fc) <- substr(colnames(plasma_dod_fc), 1,3)
plasma_t6_fc <- subset(log2(plasma_t6_fc), rownames(plasma_t6_fc) %in% c("alt", "IL18", "Ang2"))

#T6 cytof change

wide_cytof <- data.frame(cytof_data %>%
                           select(cluster_id, sample_id, frequency) %>%
                           pivot_wider(names_from = sample_id, values_from = frequency))

# wide_cytof <- data.frame(cytof_data %>%
#                            select(cluster_id, sample_id, trans_freq) %>%
#                            pivot_wider(names_from = sample_id, values_from = trans_freq))
# 

rownames(wide_cytof) <- gsub("ï", "i", wide_cytof$cluster_id)
wide_cytof <- as.matrix(select(wide_cytof, -cluster_id))
class(wide_cytof) <- "double"



cytof_t6_fc <- subset(wide_cytof, select=grepl("T6", colnames(wide_cytof)))/subset(wide_cytof, select=grepl("Baseline", colnames(wide_cytof)))
colnames(cytof_t6_fc) <- substr(colnames(cytof_t6_fc), 1,3)
cytof_t6_fc <- subset(log2(cytof_t6_fc), select = colnames(cytof_t6_fc)!="v09")
cytof_t6_fc <- rbind(cytof_t6_fc, plasma_t6_fc)

#skip to line ~220 to make pca_distance metric, line ~125 to make changing analytes
sig_cytof_t6_fc <- subset(cytof_t6_fc, grepl("activated", rownames(cytof_t6_fc)))

sig_plasma_dod_fc <- subset(plasma_dod_fc, rownames(plasma_dod_fc)%in%c(changing_analytes, "alt", "Ang2"))

# cell counts 


haem_count_data <- read.csv("~/PhD/clinical_data/vac69a/haem.csv", header=T, stringsAsFactors = F)
haem_count_data <- haem_count_data[1:(nrow(haem_count_data)-5),]

colnames(haem_count_data)[1] <- "volunteer"

haem_count_data$volunteer <- paste("v", substr(haem_count_data$volunteer, 6, 7), sep='')

lousy_haem_timepoints <- unique(haem_count_data$timepoint)
good_haem_timepoints <- c("C90", "T6", "C7 am", "Baseline", "C28", "EV", "C14 am", "Diagnosis", "Screening", "T1")

haem_time_dic <- setNames(good_haem_timepoints, lousy_haem_timepoints)

haem_count_data$timepoint <- stringr::str_replace_all(haem_count_data$timepoint, haem_time_dic)

haem_count_data <- haem_count_data %>%
  select(volunteer, timepoint, !grep("_ae", colnames(.)[13:ncol(.)], fixed = TRUE, value = TRUE))

haem_count_data <- haem_count_data[,c(1:2, 13:ncol(haem_count_data))]

long_haem_count_data <- haem_count_data %>%
  gather(cell, count, colnames(haem_count_data)[4:ncol(haem_count_data)]) %>%
  group_by(volunteer, cell) %>%
  summarise(T1_FC = count[timepoint == "T1"] / count[timepoint == "Baseline"]) %>%
  #summarise(DoD_FC = count[timepoint == "Diagnosis"] / count[timepoint == "Baseline"]) %>%
  spread(volunteer,T1_FC)
  #spread(volunteer,DoD_FC)#<3

haem_count_fc <- as.matrix(long_haem_count_data[,2:6])
rownames(haem_count_fc) <- long_haem_count_data$cell


haem_count_plot_data <- haem_count_data %>%
  gather(cell, count, colnames(haem_count_data)[4:ncol(haem_count_data)]) %>%
  filter(cell %in% c("eosinophils", "monocytes", "neutrophils", "wbc"))

haem_count_plot_data$cell <- paste(
  toupper(substr(haem_count_plot_data$cell, 1,1)),
  substr(haem_count_plot_data$cell, 2,nchar(haem_count_plot_data$cell)),
  sep="")

haem_count_plot_data$cell <- gsub("Wbc", "Total White Cells", haem_count_plot_data$cell)

fig1_theme <- theme(axis.title.x = element_blank(),
                    legend.title = element_text(size = 9), 
                    legend.text = element_text(size = 9),
                    axis.title=element_text(size=10))


haem_count_plot <- ggplot(filter(haem_count_plot_data, timepoint %in% c("Baseline", "C7 am", "C14 am", "Diagnosis", "T1", "T6", "C90")), aes(x=factor(timepoint, levels=c("Baseline", "C7 am", "C14 am", "Diagnosis", "T1", "T6", "C90")), y=count, color=volunteer, group=volunteer))+
  geom_line(aes(color=volunteer), size=1.1)+
  geom_point(fill="white", stroke=1, shape=21)+
  facet_wrap(~ cell, scales = "free", ncol=4)+
  #scale_y_continuous(trans = "log2", labels=scales::comma)+
  ylab(expression(paste("Cell Counts (", 10^6, "/mL)")))+
  theme_minimal()+
  guides(color=guide_legend(override.aes = list("size"=0.1)))+
  scale_color_manual(values=volunteer_palette)+
  fig1_theme+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(hjust=1, angle=45, size=6),
        strip.background = element_rect(fill = "white", color = "white"))


ggsave(filename = "~/PhD/cytof/vac69a/final_figures_for_paper/supp_haem_count_plot.pdf", haem_count_plot, width = 8, height=2.5)








# all_max_parasitaemias <- data.frame(#max_parasiatemia=c(4907, 8054, 16733, 7464, 21870, 15051),
#                    #`Day of Diagnosis`=c(15.5, 12.5, 15.5, 15.5, 16, 16.5),
#                    #max_ae=c(3,7,4,6,2,2),
#                    #"Lymphopenia" = c(0.3891892, 0.3989637, 0.3055556, 0.1685393, 0.3450479, 0.3073394),
#                    #"Highest Body Temperature" = c(38.5, 39.7, 37.1, 37.9, 38.7, 37.7)
#                    #dose=c(1,1,0.2,0.2,0.05,0.05)
#                    )
# 
# matrix_names <- colnames(all_max_parasitaemias)
# all_max_parasitaemias <- t(all_max_parasitaemias)
# rownames(all_max_parasitaemias) <- matrix_names
# 
# colnames(all_max_parasitaemias) <- c('v02', 'v03', 'v05', 'v06', 'v07', 'v09')

#big_fc_table <- rbind(sig_plasma_dod_fc, sig_cytof_t6_fc, plasma_t6_fc, all_max_parasitaemias[,1:5])

big_fc_table <- rbind(sig_plasma_dod_fc, sig_cytof_t6_fc, plasma_t6_fc)

# rownames(big_fc_table)[nrow(big_fc_table)] <- "Highest Body Temperature"
# 
# big_fc_table <- rbind(big_fc_table, cd4_memory_t6_summary$CD4_memory_percentage[1:5])
# rownames(big_fc_table)[nrow(big_fc_table)] <- "CD4 Memory % Baseline"

# baseline_table <- subset(big_table, select=grepl("Baseline", colnames(big_table)))
# dod_table <-  subset(big_table, select=grepl("Diagnosis", colnames(big_table)))
# t6_table <-  subset(big_table, select=grepl("T6", colnames(big_table)))

timepoint <- big_fc_table

distance <- "euclidean" # manhattan/euclidean maybe

spearman <- cor(t(timepoint), method = "pearson")

baseline_dist <- dist(spearman, method = distance, diag = FALSE, upper = FALSE, p = 2)
baseline_hclust <- hclust(baseline_dist)

#check.names=FALSE here makes sure that the +/- symbols parse and spaces aren't dots
colnames(spearman) <- gsub(".", " ", colnames(spearman), fixed=T)
rownames(spearman) <- gsub(".", " ", rownames(spearman), fixed=T)

baseline_spearman_df  <- data.frame(spearman, check.names = FALSE)
colnames(baseline_spearman_df) <- gsub(".", " ", colnames(baseline_spearman_df), fixed=T)


pearson_matrix <- as.matrix(baseline_spearman_df[rownames(baseline_spearman_df)[rev(baseline_hclust$order)],colnames(baseline_spearman_df)[rev(baseline_hclust$order)]])



col_fun_pearson <- circlize::colorRamp2(c(min(pearson_matrix), 0, max(pearson_matrix)), c("#0859C6", "black", "#FFA500"))

#doesnt work for pdf..
# colnames(pearson_matrix) <- gsub("IFNy", "IFNγ", colnames(pearson_matrix))
# rownames(pearson_matrix) <-  gsub("IFNy", "IFNγ",rownames(pearson_matrix))

library(ComplexHeatmap)

pearson_heatmap <- Heatmap(matrix = pearson_matrix,
                           cluster_rows = TRUE,
                           cluster_columns=TRUE,
                           show_row_dend = FALSE,
                           show_column_dend = TRUE,
                           show_heatmap_legend = TRUE,
                           name = "Pearson r",
                           #cluster_columns = FALSE,
                           column_names_gp = gpar(fontsize = 6),
                           row_names_gp = gpar(fontsize = 6),
                           row_names_side = "left",
                           col = col_fun_pearson,
                           column_names_rot = 45)

pdf("/home/flobuntu/PhD/cytof/vac69a/final_figures_for_paper/pearson_heatmap.pdf", height=4.5, width = 5)
draw(pearson_heatmap, padding=unit(c(2,2,2,2), "mm"))
dev.off()

# png("/home/flobuntu/PhD/cytof/vac69a/final_figures_for_paper/pearson_heatmap.png", height=4.5, width = 5, res = 900, units = "in")
# draw(pearson_heatmap, padding=unit(c(2,2,2,2), "mm"))
# dev.off()



# here's the plot using ggplot2:


# baseline_spearman_df$cluster_id_x <- rownames(baseline_spearman_df)
# 
# 
# long_baseline_spearman <- gather(baseline_spearman_df, cluster_id_y, ro, colnames(baseline_spearman_df)[1:ncol(baseline_spearman_df)-1])
# 
# 
# corr_matrix_plot <- ggplot(long_baseline_spearman, aes(x=factor(cluster_id_x, levels = colnames(spearman)[baseline_hclust$order]),
#                                                        y=factor(cluster_id_y, levels=colnames(spearman)[baseline_hclust$order])
#                                                        )
#                            )+
#   geom_tile(aes(fill=ro))+
#   #geom_text(aes(label=round(ro, digits=2)), size=1)+
#   ggplot2::scale_fill_gradient2(low = "#0859C6", mid="black", high="#FFA500", midpoint = 0, breaks=c(-1,0,1), limits=c(-1, 1))+
#   labs(fill = "r")+
#   theme_minimal()+
#   guides(fill=guide_colorbar(ticks = FALSE))+
#   theme(axis.title = element_blank(),
#         plot.background = element_blank(),
#         panel.grid = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         plot.title = element_text(hjust=0.5),
#         axis.text = element_text(size=6),
#         legend.title = element_text(size=10),
#         legend.text = element_text(size=8),
#         legend.title.align = 0.1,
#         legend.key.size=unit(2, "mm"),
#         legend.key.height = unit(5, "mm"))
# 
# #ggsave("~/PhD/multi_omics/pearson_euclidean_corr_matrix_fc.pdf", corr_matrix_plot, width=10, height=8)
# ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/sig_only_pearson_euclidean_corr_matrix_fc.pdf", corr_matrix_plot, height=4.2, width = 5, dpi=1200)
# #ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/sig_only_spearman_euclidean_corr_matrix_fc.pdf", corr_matrix_plot, height=4.2, width = 5, dpi=1200)


#log2(absolute) corr matrix ####
# 
# sig_cytof <- subset(wide_cytof_sans_v09, grepl("activated", rownames(wide_cytof_sans_v09)))
# sig_cytof <- subset(sig_cytof, select=grepl("T6", colnames(sig_cytof)))
# sig_cytof <- log2(sig_cytof)
# 
# sig_plasma_dod <- subset(log_plasma_data, rownames(log_plasma_data)%in%c(changing_analytes, "alt", "Ang2"))
# sig_plasma_dod <- subset(sig_plasma_dod, select=grepl("T6", colnames(sig_plasma_dod)))
# sig_plasma_dod <- subset(sig_plasma_dod, rownames(sig_plasma_dod) %notin% c("alt", "IL18"))
# 
# 
# sig_plasma_t6 <- subset(log_plasma_data, select=grepl("T6", colnames(log_plasma_data)))
# sig_plasma_t6 <- subset(sig_plasma_t6, rownames(sig_plasma_t6) %in% c("alt", "IL18", "Ang2"))
# 

# 
# all_max_parasitaemias <- data.frame(vol=c('v02', 'v03', 'v05', 'v06', 'v07', 'v09'),
#                                     max_parasiatemia=c(4907, 8054, 16733, 7464, 21870, 15051),
#                                     DoD=c(15.5, 12.5, 15.5, 15.5, 16, 16.5),
#                                     max_ae=c(3,7,4,6,2,2))
# clin_data <- t(all_max_parasitaemias[1:5,2:4])
# colnames(clin_data) <- paste(c('v02', 'v03', 'v05', 'v06', 'v07'), "_T6", sep='')
# 
# big_table <- rbind(sig_plasma_dod, sig_cytof, sig_plasma_t6, clin_data)
# 
# timepoint <- big_table
# 
# distance <- "euclidean" # manhattan/euclidean maybe
# 
# spearman <- cor(t(timepoint), method = "pearson")
# 
# baseline_dist <- dist(spearman, method = distance, diag = FALSE, upper = FALSE, p = 2)
# baseline_hclust <- hclust(baseline_dist)
# 
# #check.names=FALSE here makes sure that the +/- symbols parse and spaces aren't dots
# baseline_spearman_df  <- data.frame(spearman, check.names = FALSE)
# baseline_spearman_df$cluster_id_x <- rownames(baseline_spearman_df)
# 
# long_baseline_spearman <- gather(baseline_spearman_df, cluster_id_y, ro, colnames(baseline_spearman_df)[1:ncol(baseline_spearman_df)-1])
# 
# 
# 
# absolute_corr_matrix_plot <- ggplot(long_baseline_spearman, aes(x=factor(cluster_id_x, levels = colnames(spearman)[baseline_hclust$order]), y=factor(cluster_id_y, levels=colnames(spearman)[baseline_hclust$order])))+
#   geom_tile(aes(fill=ro))+
#   ggplot2::scale_fill_gradient2(low = "#0859C6", mid="black", high="#FFA500", midpoint = 0, breaks=c(-1,0,1), limits=c(-1, 1))+
#   labs(fill = expression(rho))+
#   theme(axis.title = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         plot.title = element_text(hjust=0.5),
#         axis.text = element_text(size=4),
#         legend.title = element_text(size=10),
#         legend.text = element_text(size=8),
#         legend.title.align = 0.1,
#         legend.key.size=unit(2, "mm"),
#         legend.key.height = unit(5, "mm"))
# 
# #ggsave("~/PhD/multi_omics/pearson_euclidean_corr_matrix_fc.pdf", corr_matrix_plot, width=10, height=8)
# ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/sig_only_pearson_euclidean_corr_matrix_absolute.pdf", absolute_corr_matrix_plot, height=4.2, width = 5)



#individual correlation plots ####

# T cell activation
# cd3 activation data
summary <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_t6_data.csv", stringsAsFactors = F)

summary$volunteer <- gsub("V", "v", summary$volunteer, fixed=T)


t_cell_summary <- summary %>%
  filter(timepoint=="T6") %>%
  filter(., grepl("activated", cluster_id)) %>%
  group_by(volunteer) %>%
  summarise(sum_cd3=sum(frequency))
  


combo <- data.frame("distance"=distance_frame$distance,
                    "alt"=t(2^t(log_plasma_data["alt",grepl("T6", colnames(log_plasma_data))])), 
                    t_cell_summary[1:5,]
                    #max_temp = all_max_parasitaemias$Highest.Body.Temperature[1:5],
                    #lymphopenia=all_max_parasitaemias$Lymphopenia[1:5]
                    )


volunteer_colours <- list("v02" = "#FB9A99",
                          "v03" = "#E31A1C",
                          "v05" = "#A6CEE3",
                          "v06" = "#1F78B4",
                          "v07" = "#F0E442",
                          "v09" = "#E69F00")
volunteer_palette <- unlist(unname(volunteer_colours))
names(volunteer_palette) <- names(volunteer_colours)

alt_plasma_corr_plot <- ggplot(combo, aes(x=distance, y=alt))+
  geom_point(aes(colour=volunteer))+
  xlab("Distance Traveled Plasma PCA")+
  ylab("ALT at T6")+
  geom_smooth(method="lm", se=T, fill="lightgrey")+
  
  scale_color_manual(name="Volunteer", values=volunteer_palette)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title = element_text(size=7),
        axis.text = element_text(size=6))

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/alt_plasma_corr_plot.pdf", alt_plasma_corr_plot, height=2, width=2)

cd3_plasma_corr_plot <- ggplot(combo, aes(x=distance, y=sum_cd3))+
  geom_point(aes(colour=volunteer))+
  xlab("Distance Traveled Plasma PCA")+
  ylab("CD3+ T cell activation")+
  geom_smooth(method="lm", se=T, fill="lightgrey")+
  
  scale_color_manual(name="Volunteer", values=volunteer_palette)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title = element_text(size=7),
        axis.text = element_text(size=6))


ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/cd3_plasma_corr_plot.pdf", cd3_plasma_corr_plot, height=2, width=2)




fever_plasma_corr_plot <- ggplot(combo, aes(x=distance, y=max_temp))+
  geom_point(aes(colour=volunteer))+
  xlab("Distance Traveled Plasma PCA")+
  ylab("Maximum Temperature")+
  geom_smooth(method="lm", se=T, fill="lightgrey")+
  
  scale_color_manual(name="Volunteer", values=volunteer_palette)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title = element_text(size=7),
        axis.text = element_text(size=6))


ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/fever_plasma_corr_plot.pdf", fever_plasma_corr_plot, height=2, width=2)

# cor.test(combo$distance, combo$max_temp)
# t = -1.4309, df = 3, p-value = 0.2479
# cor 
# -0.6368913 






# combo2 <- t_cell_summary
# combo2$alt <- complete_alt_timecourse[match(combo2$volunteer, names(complete_alt_timecourse))]
combo2 <- data.frame("alt"=unname(alt_timecourse[grepl("T6", names(alt_timecourse), fixed = T)]), 
                    t_cell_summary, "max_ae"=all_max_parasitaemias$max_ae)

cd3_alt_corr_plot <- ggplot(combo2, aes(y=alt, x=sum_cd3))+
  geom_point(aes(colour=volunteer))+
  ylab("ALT at T6")+
  xlab("CD3+ T cell activation")+
  scale_color_manual(values = volunteer_palette)+
  theme_minimal()+
  geom_smooth(method="lm", se=T, fill="lightgrey")+
  theme(legend.title=element_blank(),
        axis.title = element_text(size=7),
        axis.text = element_text(size=6))

library(cowplot)

vol_lgd <- get_legend(cd3_alt_corr_plot)
cd3_alt_corr_plot <- cd3_alt_corr_plot+theme(legend.position = "none")

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/cd3_alt_corr_plot.pdf", cd3_alt_corr_plot, height=2, width=2)



cd3_ae_corr_plot <- ggplot(combo2, aes(y=max_ae, x=sum_cd3))+
  geom_point(aes(colour=volunteer))+
  ylab("# of AEs around Diagnosis")+
  xlab("CD3+ T cell activation")+
  scale_color_manual(values = volunteer_palette)+
  theme_minimal()+
  scale_y_continuous(breaks = seq(0,10, by=2), limits = c(0,10))+
  geom_smooth(method="lm", se=T, fill="lightgrey")+
  theme(legend.title=element_blank(),
        axis.title = element_text(size=7),
        axis.text = element_text(size=6),
        legend.position = "none")

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/cd3_ae_corr_plot.pdf", cd3_ae_corr_plot, height=2, width=2)

# > cor.test(combo2$sum_cd3, combo2$max_ae)
# t = 1.9593, df = 4, p-value = 0.1217
#      cor 
# 0.699801


# lineage specific activation ####
summary <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_t6_data.csv", stringsAsFactors = F)

summary$lineage <- ifelse(grepl("CD4", summary$cluster_id), "CD4", "dunno")
summary$lineage <- ifelse(grepl("CD8", summary$cluster_id), "CD8", summary$lineage)
summary$lineage <- ifelse(grepl("Treg", summary$cluster_id), "Treg", summary$lineage)
summary$lineage <- ifelse(grepl("MAIT", summary$cluster_id), "MAIT", summary$lineage)
summary$lineage <- ifelse(grepl("gamma", summary$cluster_id), "gamma delta", summary$lineage)
summary$lineage <- ifelse(grepl("DN", summary$cluster_id), "DN", summary$lineage)


try <- summary %>%
  filter(timepoint=="T6") %>%
  group_by(volunteer, lineage) %>%
  filter(., grepl("activated", cluster_id)) %>%
  summarise(lin_act_freq=sum(frequency))

try2 <- summary %>%
  filter(timepoint=="T6") %>%
  group_by(volunteer, lineage) %>%
  filter(., !grepl("activated", cluster_id)) %>%
  summarise(lin_rest_freq=sum(frequency)) %>%
  select(lin_rest_freq)



lineage_freqs <- cbind(try, "lin_rest_freq"=try2$lin_rest_freq)

lineage_freqs <- lineage_freqs %>%
  group_by(volunteer) %>%
  mutate("lin_tot_freq"=lin_act_freq+lin_rest_freq) %>%
  mutate("fraction_of_lineage_activated"=lin_act_freq/lin_tot_freq*100) 

# 
# lineage_freqs$lin_tot_freq <- sum(lineage_freqs$lin_act_freq, lineage_freqs$lin_rest_freq)
# 
# lineage_freqs$fraction_of_lineage_activated <- (lineage_freqs$lin_act_freq/lineage_freqs$lin_tot_freq)*100

  lineage_spcific_activation <- ggplot(lineage_freqs, aes(x=volunteer, y=fraction_of_lineage_activated, fill=volunteer))+
  geom_bar(stat="identity")+
  facet_wrap(~lineage, scales="free")+
  theme_minimal()+
  ylab("Percentage of Lineage activated")+
  scale_fill_manual(values=volunteer_palette)+
  theme(legend.title = element_blank(),
        axis.title.x = element_blank())


lineage_freqs_sansv09 <- subset(lineage_freqs, lineage_freqs$volunteer!="v09")

lineage_freqs_sansv09 <- lineage_freqs_sansv09 %>%
  group_by(lineage) %>%
  mutate("distance"=distance_frame$distance)

lineage_freqs_sansv09 %>%
  do(broom::tidy(cor.test(.$fraction_of_lineage_activated, .$distance, method="spearman"))) #%>%


# cytof_conf_data <- filter(lineage_freqs_sansv09, lineage=="CD4")
# 
# ggplot(cytof_conf_data, aes(y=fraction_of_lineage_activated, x=distance))+
#   geom_point(aes(colour=volunteer))+
#   #xlab("Distance Traveled Plasma PCA")+
#   xlab("Systemic Inflammation")+
#   geom_smooth(method="lm", se=T, fill="lightgrey")+
#   ylab("CD4 T cell Lineage Activation")+
#   scale_color_manual(name="Volunteer", values=volunteer_palette)+
#   theme_minimal()+
#   scale_y_continuous(limits = c(0,NA), labels = scales::label_percent(accuracy = 1))+
#   theme(legend.title=element_blank(),
#         plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"),
#         axis.title = element_text(size=7),
#         axis.text = element_text(size=6))

# lineage     estimate statistic p.value parameter conf.low conf.high method                               alternative
# <chr>          <dbl>     <dbl>   <dbl>     <int>    <dbl>     <dbl> <chr>                                <chr>      
#   1 CD4          -0.423    -0.809    0.477         3   -0.951     0.732 Pearson's product-moment correlation two.sided  
# 2 CD8          -0.0279   -0.0483   0.965         3   -0.888     0.876 Pearson's product-moment correlation two.sided  
# 3 DN           -0.769    -2.09     0.128         3   -0.984     0.351 Pearson's product-moment correlation two.sided  
# 4 gamma delta   0.496     0.989    0.396         3   -0.687     0.959 Pearson's product-moment correlation two.sided  
# 5 MAIT          0.443     0.855    0.455         3   -0.721     0.953 Pearson's product-moment correlation two.sided  
# 6 Treg         -0.599    -1.30     0.286         3   -0.969     0.601 Pearson's product-moment correlation two.sided  


biochem2 <- data.frame(t(as.matrix(biochem_data)))
biochem2$Sample_ID <- rownames(biochem2)
biochem2$Timepoint <- substr(biochem2$Sample_ID, 5, nchar(biochem2$Sample_ID))
biochem2$Volunteer <- substr(biochem2$Sample_ID, 1, 3)

complete_alt_timecourse <- 2^biochem2[biochem2$Timepoint=="T6" , "alt"]
names(complete_alt_timecourse) <- biochem2[biochem2$Timepoint=="T6" , "Volunteer"]

lineage_freqs$alt <- complete_alt_timecourse[match(lineage_freqs$volunteer, names(complete_alt_timecourse))]
lineage_freqs$lymphopenia <- all_max_parasitaemias$Lymphopenia[match(lineage_freqs$volunteer, all_max_parasitaemias$vol)]

lineage_freqs$max_ae <- all_max_parasitaemias$max_ae[match(lineage_freqs$volunteer, all_max_parasitaemias$vol)]

lineage_freqs %>%
  group_by(lineage) %>%
  do(broom::tidy(cor.test(.$fraction_of_lineage_activated, .$alt, method="spearman"))) %>%


# lineage     estimate statistic p.value parameter conf.low conf.high method                               alternative
# <chr>          <dbl>     <dbl>   <dbl>     <int>    <dbl>     <dbl> <chr>                                <chr>      
#   1 CD4            0.791     2.59   0.0610         4  -0.0576     0.976 Pearson's product-moment correlation two.sided  
# 2 CD8            0.455     1.02   0.365          4  -0.566      0.925 Pearson's product-moment correlation two.sided  
# 3 DN             0.668     1.80   0.147          4  -0.313      0.959 Pearson's product-moment correlation two.sided  
# 4 gamma delta    0.107     0.214  0.841          4  -0.772      0.845 Pearson's product-moment correlation two.sided  
# 5 MAIT           0.147     0.297  0.781          4  -0.755      0.856 Pearson's product-moment correlation two.sided  
# 6 Treg           0.816     2.82   0.0478         4   0.0122     0.979 Pearson's product-moment correlation two.sided  

lineage_freqs %>%
  group_by(lineage) %>%
  do(broom::tidy(cor.test(.$fraction_of_lineage_activated, .$max_ae, method="pearson")))

# lineage     estimate statistic p.value parameter conf.low conf.high method                               alternative
# <chr>          <dbl>     <dbl>   <dbl>     <int>    <dbl>     <dbl> <chr>                                <chr>      
#   1 CD4            0.700      1.96  0.122          4 -0.259       0.964 Pearson's product-moment correlation two.sided  
# 2 CD8            0.809      2.75  0.0514         4 -0.00858     0.978 Pearson's product-moment correlation two.sided  
# 3 gamma delta    0.886      3.83  0.0187         4  0.266       0.988 Pearson's product-moment correlation two.sided  
# 4 MAIT           0.836      3.04  0.0382         4  0.0754      0.982 Pearson's product-moment correlation two.sided  
# 5 Treg           0.562      1.36  0.245          4 -0.458       0.943 Pearson's product-moment correlation two.sided  


lineage_freqs %>%
  group_by(lineage) %>%
  do(broom::tidy(cor.test(.$fraction_of_lineage_activated, .$lymphopenia, method="pearson")))
  
# lineage     estimate statistic p.value parameter conf.low conf.high method                               alternative
# <chr>          <dbl>     <dbl>   <dbl>     <int>    <dbl>     <dbl> <chr>                                <chr>      
# 1 CD4           0.130      0.262   0.807         4   -0.762     0.852 Pearson's product-moment correlation two.sided  
# 2 CD8          -0.464     -1.05    0.354         4   -0.927     0.557 Pearson's product-moment correlation two.sided  
# 3 DN            0.0779     0.156   0.883         4   -0.783     0.837 Pearson's product-moment correlation two.sided  
# 4 gamma delta  -0.175     -0.355   0.741         4   -0.864     0.742 Pearson's product-moment correlation two.sided  
# 5 MAIT          0.197      0.401   0.709         4   -0.732     0.869 Pearson's product-moment correlation two.sided  
# 6 Treg          0.121      0.244   0.819         4   -0.766     0.849 Pearson's product-moment correlation two.sided  


activation_ae_corr_plot_data <- filter(lineage_freqs, lineage %notin% c("DN")) 
#lineage_freqs$lineage <- factor(lineage_freqs$lineage, levels=c("Treg", "CD4", "gamma delta", "MAIT"))

lineage_freqs_sansv09$lineage <- gsub("gamma delta", "γδ", lineage_freqs_sansv09$lineage)

lineage_activation_ae_corr_plot <- ggplot(lineage_freqs_sansv09, aes(x=distance, y=fraction_of_lineage_activated/100))+
  geom_point(aes(colour=volunteer))+
  #xlab("Distance Traveled Plasma PCA")+
  xlab("Plasma PCA Distance Traveled")+
  geom_smooth(method="lm", se=T, fill="lightgrey")+
  ylab("T cell Lineage Activation")+
  scale_color_manual(name="Volunteer", values=volunteer_palette)+
  facet_wrap(~lineage, scales="free", ncol=6)+
  theme_minimal()+
  scale_y_continuous(limits = c(0,NA), labels = scales::label_percent(accuracy = 1))+
  theme(legend.title=element_blank(),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"),
        axis.title = element_text(size=7),
        axis.text = element_text(size=6))

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/lineage_activation_plasma_corr_plot.png", lineage_activation_ae_corr_plot, height=2, width=7)
ggsave("~/PhD/figures_for_thesis/chapter_2/lineage_activation_plasma_corr_plot.png", lineage_activation_ae_corr_plot, height=2, width=7)




activation_alt_corr_plot_data <- filter(lineage_freqs, lineage %in% c("CD4", "gamma delta", "MAIT", "Treg")) 
activation_alt_corr_plot_data$lineage <- factor(activation_alt_corr_plot_data$lineage, levels=c("Treg", "CD4", "gamma delta", "MAIT"))

lineage_activation_alt_corr_plot <- ggplot(activation_alt_corr_plot_data, aes(x=alt, y=fraction_of_lineage_activated/100))+
  geom_point(aes(colour=volunteer))+
  #xlab("Distance Traveled Plasma PCA")+
  xlab("ALT at T6")+
  geom_smooth(method="lm", se=T, fill="lightgrey")+
  ylab("T cell Lineage Activation")+
  scale_color_manual(name="Volunteer", values=volunteer_palette)+
  facet_wrap(~lineage, scales="free", ncol=4)+
  theme_minimal()+
  scale_y_continuous(limits = c(0,NA), labels = scales::label_percent(accuracy = 1))+
  theme(legend.title=element_blank(),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"),
        axis.title = element_text(size=7),
        axis.text = element_text(size=6))

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/lineage_activation_alt_corr_plot.pdf", lineage_activation_alt_corr_plot, height=2, width=6)



top_row <- cowplot::plot_grid(arrow_pca_plot, cd3_plasma_corr_plot,  alt_plasma_corr_plot, ncol=3,
                              rel_widths = c(1,1,1,1), align = "v", axis="b")

#ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/all_correlations.pdf", both_rows, height=5, width=8)
ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/arrow_pca_plasma_cd3_alt_correlations.pdf", top_row, height=2, width=6)



#more correlations, parasitaemia, day of diagnosis

combo <- cbind(combo, all_max_parasitaemias[1:5, 2:3])
combo <- combo[,!duplicated(rownames(combo))]

ggplot(combo, aes(x=distance, y=DoD))+
  geom_point(aes(colour=volunteer))+
  xlab("Distance Traveled Plasma PCA")+
  ylab("Day of Diagnosis")+
  geom_smooth(method="lm", se=T, fill="lightgrey")+
  scale_color_manual(name="Volunteer", values=volunteer_palette)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title = element_text(size=7),
        axis.text = element_text(size=6))

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/distance_vs_dod.pdf")

ggplot(combo, aes(x=distance, y=max_parasiatemia))+
  geom_point(aes(colour=volunteer))+
  xlab("Distance Traveled Plasma PCA")+
  ylab("Max. Parasitaemia")+
  geom_smooth(method="lm", se=T, fill="lightgrey")+
  scale_color_manual(name="Volunteer", values=volunteer_palette)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title = element_text(size=7),
        axis.text = element_text(size=6))

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/distance_vs_max_parasites.pdf")


combo2 <- cbind(combo2, all_max_parasitaemias[, 2:3])
combo2 <- combo2[,!duplicated(rownames(combo2))]

ggplot(combo2, aes(x=sum_cd3, y=DoD))+
  geom_point(aes(colour=volunteer))+
  xlab("CD3+ T cell activation")+
  ylab("Day of Diagnosis")+
  geom_smooth(method="lm", se=T, fill="lightgrey")+
  scale_color_manual(name="Volunteer", values=volunteer_palette)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title = element_text(size=7),
        axis.text = element_text(size=6))

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/tcell_activation_vs_dod.pdf")


ggplot(combo2, aes(x=sum_cd3, y=max_parasiatemia))+
  geom_point(aes(colour=volunteer))+
  xlab("CD3+ T cell activation")+
  ylab("Max. Parasitaemia")+
  geom_smooth(method="lm", se=T, fill="lightgrey")+
  scale_color_manual(name="Volunteer", values=volunteer_palette)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title = element_text(size=7),
        axis.text = element_text(size=6))

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/tcell_activation_vs_max_parasiatemia.pdf")

