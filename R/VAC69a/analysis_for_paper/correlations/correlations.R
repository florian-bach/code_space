# big correlation matrix containing plasma, biochem and cytof data ####

library(dplyr)
library(ggplot2)


my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")
time_col=colorspace::sequential_hcl(5, palette = "Purple Yellow")

time_palette <- c("Baseline"=time_col[4], "C10"=time_col[3], "DoD"=time_col[2], "T6"=time_col[1])


`%notin%` <- Negate(`%in%`)

corr_matrix_theme <-
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1),
        plot.title = element_text(hjust=0.5),
        axis.text = element_text(size=7),
        legend.title.align = 0.2)


# 
# log_plasma_data <- read.csv("~/PhD/multi_omics/all_log_plasma_data_with_alt.csv", header=T, stringsAsFactors = F, row.names = 1)
# 
# wide_cytof <- read.csv("~/PhD/multi_omics/cytof_cluster_freqs.csv", header=T, stringsAsFactors = F, row.names = 1)
# 
# biochem_data <- read.csv("~/PhD/multi_omics/mofa_log2_biochem_data.csv", header=T, stringsAsFactors = F, row.names = 1)
# 


### PLASMA DATA generation ####


# data <- read.csv("C:/Users/Florian/PhD/oxford/vac69/plasma_analytes_vivax_no_inequalitites.csv")
# dataplasma <- read.csv("/Users/s1249052/PhD/plasma/vac69a/Vivax_plasma_analytes2_no_inequalities.csv")
dataplasma <- read.csv("~/PhD/plasma/vac69a/Vivax_plasma_analytes2_no_inequalities.csv")
dataplasma <- subset(dataplasma, dataplasma$Volunteer!="v009")
t_dataplasma <- t(dataplasma)#

colnames(t_dataplasma) <- paste(t_dataplasma[1,], t_dataplasma[2,], sep='_')

dataplasma <- data.frame(t_dataplasma[3:nrow(t_dataplasma),])

colnames(dataplasma) <- gsub("C.1", "Baseline", colnames(dataplasma))
colnames(dataplasma) <- gsub("v00", "V0", colnames(dataplasma))
colnames(dataplasma) <- gsub("T", "DoD", colnames(dataplasma))
colnames(dataplasma) <- gsub("DoD.6", "T6", colnames(dataplasma))

#convert to numerical matri
plasma_data <- as.matrix(dataplasma[,colnames(dataplasma)%in%colnames(wide_cytof)]) # <3

changing_analytes <- sig_analytes <- scan("~/PhD/plasma/vac69a/analytes_sorted_by_padj.txt", what="", skip = 1)
changing_analytes <- changing_analytes[1:18]

rownames(plasma_data) <- gsub("pg.ml.", "", rownames(plasma_data), fixed = T)
rownames(plasma_data) <- gsub(".", "",rownames (plasma_data), fixed=T)

plasma_data <- subset(plasma_data, rownames(plasma_data)%in%changing_analytes)

class(plasma_data) <- "double" # <3

#put it together, names() makes view names
#log transform plasma data

log_plasma_data <- log2(plasma_data)

#this "alt_timecourse" variable was made ~50 lines below- don't ask...
log_plasma_data <- rbind(log_plasma_data, alt=log2(alt_timecourse[!grepl("V09", names(alt_timecourse), fixed=T)]))

#write.csv(log_plasma_data, "~/PhD/multi_omics/all_log_plasma_data_with_alt.csv", row.names = T)
# log_plasma_data <- read.csv("~/PhD/multi_omics/log_plasma_data_with_alt.csv", header=T, stringsAsFactors = F, row.names = 1)

wide_cytof <- read.csv("~/PhD/multi_omics/cytof_cluster_freqs.csv", header=T, stringsAsFactors = F, row.names = 1)

# log_plasma_data <- rbind(log_plasma_data, substr(colnames(log_plasma_data), 1,3), substr(colnames(log_plasma_data), 5,nchar(colnames(log_plasma_data))))
# rownames(log_plasma_data)[(nrow(log_plasma_data)-1):nrow(log_plasma_data)] <- c("Volunteer", "Timepoint")


###### clinical & biochem data


## biochem data generation ####

#data <- read.csv("/Users/s1249052/PhD/clinical data/vac69a/biochem.csv")
data <- read.csv("~/PhD/clinical_data/vac69a/biochem.csv")

colnames(data)[1] <- "Volunteer"

data$Volunteer <- paste("V", substr(data$Volunteer, 6, 7), sep='')

short <- filter(data, timepoint %in% c("_C_1", "_T6", "_EP"))
short2 <- select(short,  Volunteer, timepoint, sodium, potassium, creatinine, bilirubin, alt, alkphos, albumin, ast, ggt)

short2$timepoint <- as.character(short2$timepoint)

short2$timepoint[short2$timepoint=="_C_1"] <- "Baseline"
short2$timepoint[short2$timepoint=="_T6"] <- "T6"
short2$timepoint[short2$timepoint=="_EP"] <- "DoD"

short2 <- filter(short2, timepoint %in% c("Baseline", "T6", "DoD"))


mofa_biochem <- t(short2)
colnames(mofa_biochem) <- paste(mofa_biochem[1,], mofa_biochem[2,], sep="_") 
mofa_biochem <- mofa_biochem[3:nrow(mofa_biochem),]

mofa_biochem <- as.matrix(mofa_biochem[,colnames(mofa_biochem)%in%colnames(wide_cytof)]) # <3
class(mofa_biochem) <- "double"

biochem_data <- log2(mofa_biochem[1:7,])
#biochem_data <- ifelse(is.infinite(biochem_data), NA, biochem_data)
biochem_data[5,1] <- 3.7
biochem_data <- biochem_data[c(3,5:7), ]

#write.csv(biochem_data, "~/PhD/multi_omics/mofa_log2_biochem_data.csv", row.names = T)

wide_cytof <- read.csv("~/PhD/multi_omics/cytof_cluster_freqs.csv", header=T, stringsAsFactors = F, row.names = 1)

alt_timecourse <- mofa_biochem[5,order(colnames(mofa_biochem))]
alt_timecourse[7] <- 13
cor(t(wide_cytof), alt_timecourse)


cytof_data$alt <- alt_timecourse[match(cytof_data$sample_id, names(alt_timecourse))]

correlation_dfs <- split(cytof_data, cytof_data$cluster_id)
names(correlation_dfs) <- lapply(correlation_dfs, function(x)unique(x$cluster_id))
correlation_results <- lapply(correlation_dfs, function(x)cor.test(x$frequency, x$alt, method="pearson"))

corr_res <- lapply(correlation_results, function(x) cbind("p_value"=x$p.value, "rho"=x$estimate))

corr_res <- data.frame(do.call(rbind, corr_res))
corr_res <- corr_res[order(corr_res$p_value),]
rownames(corr_res) <- names(correlation_results)
corr_res$p_adj <- p.adjust(corr_res$p_value, "BH")

subset(corr_res, corr_res$p_adj<=0.05)

#                               p_value        rho        p_adj
# activated  CD4 CM         1.851518e-06 0.8763470 5.924857e-05
# activated  CD8 EM         1.139985e-04 0.7849951 1.823976e-03
# activated  MAIT           5.362225e-04 0.7331917 5.719707e-03
# activated  Treg EM        1.134533e-03 0.7030937 8.858416e-03 # notably this one isn't significant when not accounting for ALT bins
# activated  Vd2+           1.384127e-03 0.6944552 8.858416e-03
# activated HLADR+ DN TEMRA 7.048489e-03 0.6111414 3.759194e-02 # notably this one isn't significant when not accounting for ALT bins


# correlations between PCA distance traveled and T cell activation ####

#PCA distance traveled 
plasma_mds <- limma::plotMDS(log_plasma_data, plot = FALSE)

plasma_df <- data.frame(MDS1 = plasma_mds$x, MDS2 = plasma_mds$y)
plasma_df$Sample_ID <- rownames(plasma_df)
plasma_df$Timepoint <- substr(plasma_df$Sample_ID, 5, nchar(plasma_df$Sample_ID))
plasma_df$Volunteer <- substr(plasma_df$Sample_ID, 1, 3)



distance_traveled <- filter(plasma_df,Timepoint%in%c("Baseline", "DoD"))
distance_traveled_dfs <- split(distance_traveled, distance_traveled$Volunteer)

distance_traveled_mds1 <- lapply(distance_traveled_dfs, function(x) x$sum_sq_MDS2=(x$MDS2[2]-x$MDS2[1])^2)

distance_traveled_mds2 <- lapply(distance_traveled_dfs, function(x) x$sum_sq_MDS1=(x$MDS1[2]-x$MDS1[1])^2)

distance_traveled <- lapply(names(distance_traveled_mds2), function(x)
  data.frame(Volunteer=x, "MDS1"=as.numeric(distance_traveled_mds1[[x]]), "MDS2"=as.numeric(distance_traveled_mds2[[x]])))

df <- do.call(rbind, distance_traveled)

df$distance <- sqrt(apply(df[,2:3], 1, sum))

pca_distance <- df$distance
names(pca_distance) <- c("V02", "V03", "V05", "V06", "V07")
# T cell activation
# cd3 activation data
summary <- data.table::fread("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/activation_barchart_data")

t_cell_summary <- summary %>%
  group_by(volunteer, timepoint) %>%
  summarise(sum_cd3=sum(frequency)) %>%
  filter(timepoint=="T6")

combo <- data.frame("distance"=df$distance, "alt"=2^t(log_plasma_data["alt",grepl("T6", colnames(log_plasma_data))]),  t_cell_summary[1:5,])


alt_plasma_corr_plot <- ggplot(combo, aes(x=distance, y=alt))+
  geom_point(aes(colour=volunteer))+
  xlab("Distance Traveled Plasma PCA")+
  ylab("alt")+
  scale_color_manual(name="Volunteer", values=my_paired_palette)+
  theme_minimal()

# ggsave("~/PhD/multi_omics/alt_plasma_corr_plot.png", alt_plasma_corr_plot, height=3, width=3)

cd3_plasma_corr_plot <- ggplot(combo, aes(x=distance, y=sum_cd3))+
  geom_point(aes(colour=volunteer))+
  xlab("Distance Traveled Plasma PCA")+
  ylab("CD3+ T cell activation")+
  scale_color_manual(name="Volunteer", values=my_paired_palette)+
  theme_minimal()

# ggsave("~/PhD/multi_omics/cd3_alt_corr_plot.png", cd3_plasma_corr_plot, height=3, width=3)

combo2 <- t_cell_summary
combo2$alt <- complete_alt_timecourse[match(combo2$volunteer, names(complete_alt_timecourse))]




cd3_alt_corr_plot <- ggplot(combo2, aes(x=alt, y=sum_cd3))+
  geom_point(aes(colour=volunteer))+
  #xlab("Distance Traveled Plasma PCA")+
  xlab("alt")+
  ylab("CD3+ T cell activation")+
  scale_color_manual(name="Volunteer", values=my_paired_palette)+
  theme_minimal()


# ggsave("~/PhD/multi_omics/cd3_alt_corr_plot.png", cd3_alt_corr_plot, height=3, width=3)

# cor.test(combo$distance, combo$sum_cd3, method = "pearson")
# data:  combo$distance and combo$sum_cd3
#t = -0.2647, df = 3, p-value = 0.8084
#   cor
#-0.1510704


# > cor.test(combo2$alt, combo2$sum_cd3, method = "pearson")
# t = 2.36, df = 4, p-value = 0.07766
# 95 percent confidence interval:
#  -0.1277660  0.9724059
# sample estimates:
#       cor 
# 0.7629004 
# 
# > cor.test(combo$alt, combo$distance, method = "pearson")
# data:  combo$alt and combo$distance
# t = -1.4777, df = 3, p-value = 0.236
#   cor 
# -0.6490362 



# lineage specific activation
wide_cytof <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/cluster_and_lineage_stats_for_barcharts.csv", header = T, stringsAsFactors = F)

#wide_cytof$lineage <- cluster_dic$lineage[match(wide_cytof$cluster_id, cluster_dic$new_cluster)]

try <- wide_cytof %>%
  filter(timepoint=="T6") %>%
  group_by(volunteer, lineage) %>%
  filter(., grepl("activated", cluster_id)) %>%
  summarise(lin_act_freq=sum(frequency))

try2 <- wide_cytof %>%
  filter(timepoint=="T6") %>%
  group_by(volunteer, lineage) %>%
  filter(., !grepl("activated", cluster_id)) %>%
  summarise(lin_rest_freq=sum(frequency)) %>%
  select(lin_rest_freq)

lineage_freqs <- cbind(try, try2)

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
  scale_fill_manual(values=my_paired_palette)+
  theme(legend.title = element_blank(),
        axis.title.x = element_blank())


# how to make faceted plot where each combination of lineages is plotted against each other???

# lineage_freqs$lienage2 <- lineage_freqs$lineage
# lineage_spcific_activation2 <- ggplot(lineage_freqs, aes(x=fraction_of_lineage_activated,y=fraction_of_lineage_activated, fill=volunteer, color=volunteer))+
#   geom_point()+
#   facet_grid(lienage2~lineage, scales="free")+
#   theme_minimal()+
#   ylab("Percentage of Lineage activated")+
#   scale_fill_manual(values=my_paired_palette)+
#   scale_color_manual(values=my_paired_palette)+
#   
#   theme(legend.title = element_blank(),
#         axis.title.x = element_blank())

# ggsave("~/PhD/multi_omics/lineage_specific_activation2.png", lineage_spcific_activation2, height=5, width=7)


lineage_freqs_sansv09 <- subset(lineage_freqs, lineage_freqs$volunteer!="V09")

lineage_freqs_sansv09 %>%
  group_by(lineage) %>%
  mutate("distance"=df$distance)%>%
  do(broom::tidy(cor.test(.$fraction_of_lineage_activated, .$distance, method="pearson")))

# A tibble: 6 x 9
# Groups:   lineage [6]
# lineage estimate statistic p.value parameter conf.low conf.high method                               alternative
# <fct>      <dbl>     <dbl>   <dbl>     <int>    <dbl>     <dbl> <chr>                                <chr>      
#   1 CD4      -0.424    -0.811   0.477          3   -0.951    0.732  Pearson's product-moment correlation two.sided  
# 2 CD8      -0.0281   -0.0487  0.964          3   -0.888    0.876  Pearson's product-moment correlation two.sided  
# 3 DN       -0.899    -3.55    0.0380         3   -0.993   -0.0805 Pearson's product-moment correlation two.sided  
# 4 gd        0.323     0.590   0.597          3   -0.782    0.938  Pearson's product-moment correlation two.sided  
# 5 MAIT      0.442     0.854   0.456          3   -0.722    0.953  Pearson's product-moment correlation two.sided  
# 6 Treg     -0.599    -1.30    0.285          3   -0.969    0.600  Pearson's product-moment correlation two.sided  


lineage_freqs_sansv09_distance <- lineage_freqs_sansv09 %>%
  group_by(lineage) %>%
  mutate("distance"=df$distance)


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
# 
  # lineage estimate statistic p.value parameter conf.low conf.high method                               alternative
  # <chr>      <dbl>     <dbl>   <dbl>     <int>    <dbl>     <dbl> <chr>                                <chr>      
  #   1 CD4        0.791     2.59   0.0610         4  -0.0576     0.976 Pearson's product-moment correlation two.sided  
  # 2 CD8        0.455     1.02   0.365          4  -0.566      0.925 Pearson's product-moment correlation two.sided  
  # 3 DN         0.692     1.92   0.128          4  -0.273      0.963 Pearson's product-moment correlation two.sided  
  # 4 gd         0.281     0.587  0.589          4  -0.687      0.890 Pearson's product-moment correlation two.sided  
  # 5 MAIT       0.147     0.297  0.781          4  -0.755      0.856 Pearson's product-moment correlation two.sided  
  # 6 Treg       0.816     2.82   0.0478         4   0.0122     0.979 Pearson's product-moment correlation two.sided  


lineage_activation_alt_corr_plot <- ggplot(lineage_freqs, aes(x=alt, y=fraction_of_lineage_activated))+
  geom_point(aes(colour=volunteer))+
  #xlab("Distance Traveled Plasma PCA")+
  xlab("alt")+
  ylab("T cell Lineage Activation")+
  scale_color_manual(name="Volunteer", values=my_paired_palette)+
  facet_wrap(~lineage, scales="free")+
  theme_minimal()

# ggsave("~/PhD/multi_omics/lineage_activation_alt_corr_plot.png", lineage_activation_alt_corr_plot, height=5, width=7)


lineage_activation_plasma_corr_plot <- ggplot(lineage_freqs_sansv09_distance, aes(x=distance, y=fraction_of_lineage_activated))+
  geom_point(aes(colour=volunteer))+
  xlab("Distance Traveled Plasma PCA")+
  ylab("T cell Lineage Activation")+
  scale_color_manual(name="Volunteer", values=my_paired_palette)+
  facet_wrap(~lineage, scales="free")+
  theme_minimal()

ggsave("~/PhD/multi_omics/lineage_activation_alt_plot.png", lineage_activation_alt_corr_plot, height=5, width=7)








wide_cytof_sans_v09 <- subset(wide_cytof, select=!grepl("V09", colnames(wide_cytof)))
#wide_cytof_sans_v09 <- subset(wide_cytof_sans_v09, select=!grepl("C10", colnames(wide_cytof)))

biochem_data_sans_v9 <- subset(biochem_data, select=!grepl("V09", colnames(biochem_data)))

big_table <- rbind(wide_cytof_sans_v09, log_plasma_data, biochem_data_sans_v9[c(1,3,4),])

# DoD plasma change; log2FC
plasma_dod_fc <- 2^subset(log_plasma_data, select=grepl("DoD", colnames(log_plasma_data)))/2^subset(log_plasma_data, select=grepl("Baseline", colnames(log_plasma_data)))
colnames(plasma_dod_fc) <- substr(colnames(plasma_dod_fc), 1,3)
plasma_dod_fc <- subset(log2(plasma_dod_fc), rownames(plasma_dod_fc) %notin% c("alt", "IL18", "Ang2"))

#t6 plasma change
plasma_t6_fc <- 2^subset(log_plasma_data, select=grepl("T6", colnames(log_plasma_data)))/2^subset(log_plasma_data, select=grepl("Baseline", colnames(log_plasma_data)))
colnames(plasma_t6_fc) <- substr(colnames(plasma_dod_fc), 1,3)
plasma_t6_fc <- subset(log2(plasma_t6_fc), rownames(plasma_t6_fc) %in% c("alt", "IL18", "Ang2"))

#T6 cytof change
cytof_t6_fc <- subset(wide_cytof, select=grepl("T6", colnames(wide_cytof)))/subset(wide_cytof, select=grepl("Baseline", colnames(wide_cytof)))
colnames(cytof_t6_fc) <- substr(colnames(cytof_t6_fc), 1,3)
cytof_t6_fc <- subset(log2(cytof_t6_fc), select = colnames(cytof_t6_fc)!="V09")
cytof_t6_fc <- rbind(cytof_t6_fc, plasma_t6_fc)

#skip to line ~220 to make pca_distance metric, line ~125 to make changing analytes
sig_cytof_t6_fc <- subset(cytof_t6_fc, grepl("activated", rownames(cytof_t6_fc)))

sig_plasma_dod_fc <- subset(plasma_dod_fc, rownames(plasma_dod_fc)%in%c(changing_analytes, "alt", "Ang2"))

big_fc_table <- rbind(sig_plasma_dod_fc, sig_cytof_t6_fc, plasma_t6_fc, "pca_distance"=pca_distance)
  
#big_fc_table <- rbind(plasma_dod_fc, cytof_t6_fc, "pca_distance"=pca_distance)

baseline_table <- subset(big_table, select=grepl("Baseline", colnames(big_table)))
dod_table <-  subset(big_table, select=grepl("DoD", colnames(big_table)))
t6_table <-  subset(big_table, select=grepl("T6", colnames(big_table)))

timepoint <- big_fc_table



#euclidean + spearman/pearson: alt smack in the middle of the activated T cell clusters

distance <- "euclidean" # manhattan/euclidean maybe

spearman <- cor(t(timepoint), method = "pearson")

baseline_dist <- dist(spearman, method = distance, diag = FALSE, upper = FALSE, p = 2)
baseline_hclust <- hclust(baseline_dist)

#check.names=FALSE here makes sure that the +/- symbols parse and spaces aren't dots
baseline_spearman_df  <- data.frame(spearman, check.names = FALSE)
baseline_spearman_df$cluster_id_x <- rownames(baseline_spearman_df)

long_baseline_spearman <- gather(baseline_spearman_df, cluster_id_y, ro, colnames(baseline_spearman_df)[1:ncol(baseline_spearman_df)-1])


# ggplot(long_baseline_spearman, aes(x=factor(long_baseline_spearman$cluster_id_x, levels = specific_order), y=factor(long_baseline_spearman$cluster_id_y, levels=specific_order)))+
# corr_matrix_plot <- ggplot(long_baseline_spearman, aes(x=factor(cluster_id_x, levels = colnames(spearman)[baseline_hclust$order]), y=factor(cluster_id_y, levels=colnames(spearman)[baseline_hclust$order])))+
#   geom_tile(aes(fill=ro))+
#   #viridis::scale_fill_viridis(option="A")+
#   ggplot2::scale_fill_gradient2(low = "#0859C6", mid="black", high="#FFA500", midpoint = 0, breaks=seq(-1,1, by=0.5), limits=c(-1,1))+
#   labs(fill = expression(rho))+
#   corr_matrix_theme
# 
# #ggsave("~/PhD/multi_omics/pearson_euclidean_corr_matrix_fc.png", corr_matrix_plot, width=10, height=8)
# ggsave("~/PhD/multi_omics/sig_only_pearson_euclidean_corr_matrix_fc.png", corr_matrix_plot, height=6, width = 7)

  # Cytof data generation ####
cytof_data <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/cluster_and_lineage_stats_for_barcharts.csv", header = T, stringsAsFactors = F)
cytof_data <- subset(cytof_data, cytof_data$timepoint!="C10")

cytof_data$trans_freq=scale(asin(sqrt(cytof_data$frequency/100)), center = TRUE, scale = TRUE)

wide_cytof <- data.frame(cytof_data %>%
                           select(cluster_id, sample_id, frequency) %>%
                           pivot_wider(names_from = sample_id, values_from = frequency))

# wide_cytof <- data.frame(cytof_data %>%
#                            select(cluster_id, sample_id, trans_freq) %>%
#                            pivot_wider(names_from = sample_id, values_from = trans_freq))

rownames(wide_cytof) <- gsub("ï", "i", wide_cytof$cluster_id)
wide_cytof <- as.matrix(select(wide_cytof, -cluster_id))
class(wide_cytof) <- "double"

#write.csv(wide_cytof, "~/PhD/multi_omics/cytof_cluster_freqs.csv", row.names = T)
