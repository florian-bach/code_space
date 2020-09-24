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
changing_analytes <- changing_analytes[1:18]

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


plasma_mds_plot <- ggplot(plasma_df, aes(x=MDS1, y=MDS2, color=Volunteer))+
  geom_point(aes(shape=Timepoint))+
  #ggrepel::geom_label_repel(aes_string(label = "sample_id"), show.legend = FALSE)+ 
  theme_minimal()+
  scale_color_manual(values = my_paired_palette)+
  theme()


ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/plasma_mds.png", plasma_mds, height=4, width=4)

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

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/arrow_pca_plot.png", arrow_pca_plot, height=3, width =3)


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

df <- do.call(rbind, distance_traveled)

df$distance <- sqrt(apply(df[,2:3], 1, sum))

pca_distance <- df$distance
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

big_fc_table <- rbind(sig_plasma_dod_fc, sig_cytof_t6_fc, plasma_t6_fc, "pca_distance"=pca_distance)

#big_fc_table <- rbind(plasma_dod_fc, cytof_t6_fc, "pca_distance"=pca_distance)

baseline_table <- subset(big_table, select=grepl("Baseline", colnames(big_table)))
dod_table <-  subset(big_table, select=grepl("Diagnosis", colnames(big_table)))
t6_table <-  subset(big_table, select=grepl("T6", colnames(big_table)))

timepoint <- big_fc_table

distance <- "euclidean" # manhattan/euclidean maybe

spearman <- cor(t(timepoint), method = "pearson")

baseline_dist <- dist(spearman, method = distance, diag = FALSE, upper = FALSE, p = 2)
baseline_hclust <- hclust(baseline_dist)

#check.names=FALSE here makes sure that the +/- symbols parse and spaces aren't dots
baseline_spearman_df  <- data.frame(spearman, check.names = FALSE)
baseline_spearman_df$cluster_id_x <- rownames(baseline_spearman_df)

long_baseline_spearman <- gather(baseline_spearman_df, cluster_id_y, ro, colnames(baseline_spearman_df)[1:ncol(baseline_spearman_df)-1])



corr_matrix_plot <- ggplot(long_baseline_spearman, aes(x=factor(cluster_id_x, levels = colnames(spearman)[baseline_hclust$order]), y=factor(cluster_id_y, levels=colnames(spearman)[baseline_hclust$order])))+
  geom_tile(aes(fill=ro))+
  #viridis::scale_fill_viridis(option="A")+
  ggplot2::scale_fill_gradient2(low = "#0859C6", mid="black", high="#FFA500", midpoint = 0, breaks=c(-1,0,1), limits=c(-1, 1))+
  labs(fill = expression(rho))+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust=0.5),
        axis.text = element_text(size=4),
        legend.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title.align = 0.1,
        legend.key.size=unit(2, "mm"),
        legend.key.height = unit(5, "mm"))

#ggsave("~/PhD/multi_omics/pearson_euclidean_corr_matrix_fc.png", corr_matrix_plot, width=10, height=8)
ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/sig_only_pearson_euclidean_corr_matrix_fc.png", corr_matrix_plot, height=4.2, width = 5)

#individual correlation plots

# T cell activation
# cd3 activation data
summary <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_t6_data.csv", stringsAsFactors = F)

summary$volunteer <- gsub("V", "v", t_cell_summary$volunteer, fixed=T)


t_cell_summary <- summary %>%
  group_by(volunteer, timepoint) %>%
  filter(timepoint=="T6") %>%
  filter(., grepl("activated", cluster_id)) %>%
  summarise(sum_cd3=sum(frequency))
  


combo <- data.frame("distance"=df$distance,
                    "alt"=t(2^t(log_plasma_data["alt",grepl("T6", colnames(log_plasma_data))])), 
                    t_cell_summary[1:5,])


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

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/alt_plasma_corr_plot.png", alt_plasma_corr_plot, height=2, width=2)

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


ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/cd3_plasma_corr_plot.png", cd3_plasma_corr_plot, height=2, width=2)

# combo2 <- t_cell_summary
# combo2$alt <- complete_alt_timecourse[match(combo2$volunteer, names(complete_alt_timecourse))]

lm_eqn = function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}

cd3_alt_corr_plot <- ggplot(combo, aes(y=alt, x=sum_cd3))+
  geom_point(aes(colour=volunteer))+
  ylab("ALT at T6")+
  xlab("CD3+ T cell activation")+
  scale_color_manual(name="Volunteer", values=volunteer_palette)+
  theme_minimal()+
  geom_smooth(method="lm", se=T, fill="lightgrey")+
  
  #geom_text(aes(x = 12, y = 125, label = lm_eqn(lm(alt ~ sum_cd3, combo))), parse = TRUE, size=1)+
  #geom_smooth(method="lm")+
  theme(legend.position = "none",
        axis.title = element_text(size=7),
        axis.text = element_text(size=6))

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/cd3_alt_corr_plot.png", cd3_alt_corr_plot, height=2, width=2)

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

lineage_freqs_sansv09 %>%
  group_by(lineage) %>%
  mutate("distance"=df$distance)%>%
  do(broom::tidy(cor.test(.$fraction_of_lineage_activated, .$distance, method="pearson")))
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

lineage_freqs %>%
  group_by(lineage) %>%
  do(broom::tidy(cor.test(.$fraction_of_lineage_activated, .$alt, method="pearson")))
 

lineage_activation_alt_corr_plot <- ggplot(lineage_freqs, aes(x=alt, y=fraction_of_lineage_activated))+
  geom_point(aes(colour=volunteer))+
  #xlab("Distance Traveled Plasma PCA")+
  xlab("ALT at T6")+
  geom_smooth(method="lm", se=T, fill="lightgrey")+
  ylab("T cell Lineage Activation")+
  scale_color_manual(name="Volunteer", values=volunteer_palette)+
  facet_wrap(~lineage, scales="free")+
  theme_minimal()+
  theme(plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"),
        axis.title = element_text(size=7),
        axis.text = element_text(size=6))

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/lineage_activation_alt_corr_plot.png", lineage_activation_alt_corr_plot, height=3, width=6)

top_row <- cowplot::plot_grid(cd3_plasma_corr_plot,  alt_plasma_corr_plot, cd3_alt_corr_plot, ncol=3)
both_rows <- cowplot::plot_grid(top_row, lineage_activation_alt_corr_plot, nrow=2, rel_heights = c(1,2), align = "v", axis="b")

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/all_correlations.png", both_rows, height=6, width=8)
