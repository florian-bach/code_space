log_plasma_data <- read.csv("~/PhD/multi_omics/all_log_plasma_data_with_alt.csv", header=T, stringsAsFactors = F, row.names = 1)

wide_cytof <- read.csv("~/PhD/multi_omics/cytof_cluster_freqs.csv", header=T, stringsAsFactors = F, row.names = 1)

biochem_data <- read.csv("~/PhD/multi_omics/mofa_log2_biochem_data.csv", header=T, stringsAsFactors = F, row.names = 1)

`%notin%` <- Negate(`%in%`)

corr_matrix_theme <-
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1),
        plot.title = element_text(hjust=0.5),
        axis.text = element_text(size=7),
        legend.title.align = 0.2)



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


big_fc_table <- rbind(plasma_dod_fc, cytof_t6_fc)


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
corr_matrix_plot <- ggplot(long_baseline_spearman, aes(x=factor(cluster_id_x, levels = colnames(spearman)[baseline_hclust$order]), y=factor(cluster_id_y, levels=colnames(spearman)[baseline_hclust$order])))+
  geom_tile(aes(fill=ro))+
  #viridis::scale_fill_viridis(option="A")+
  ggplot2::scale_fill_gradient2(low = "#0859C6", mid="black", high="#FFA500", midpoint = 0, breaks=seq(-1,1, by=0.5), limits=c(-1,1))+
  labs(fill = expression(rho))+
  corr_matrix_theme

ggsave("~/PhD/multi_omics/pearson_euclidean_corr_matrix_fc.png", corr_matrix_plot, width=10, height=8)

  # Cytof data generation ####
cytof_data <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/cluster_counts_and_freqs.csv", header = T, stringsAsFactors = F)
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

write.csv(wide_cytof, "~/PhD/multi_omics/cytof_cluster_freqs.csv", row.names = T)

### PLASMA DATA generation ####


data <- read.csv("C:/Users/Florian/PhD/oxford/vac69/plasma_analytes_vivax_no_inequalitites.csv")
dataplasma <- read.csv("/Users/s1249052/PhD/plasma/vac69a/Vivax_plasma_analytes2_no_inequalities.csv")
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

#plasma_data <- subset(plasma_data, rownames(plasma_data)%in%changing_analytes)

class(plasma_data) <- "double" # <3

#put it together, names() makes view names
#log transform plasma data

log_plasma_data <- log2(plasma_data)

#this "alt_timecourse" variable was made ~50 lines below- don't ask...
log_plasma_data <- rbind(log_plasma_data, alt=log2(alt_timecourse[!grepl("V09", names(alt_timecourse), fixed=T)]))

write.csv(log_plasma_data, "~/PhD/multi_omics/all_log_plasma_data_with_alt.csv", row.names = T)
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

write.csv(biochem_data, "~/PhD/multi_omics/mofa_log2_biochem_data.csv", row.names = T)

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
