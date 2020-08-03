cytof_data <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/cluster_counts_and_freqs.csv", header = T, stringsAsFactors = F)
cytof_data <- subset(cytof_data, cytof_data$timepoint!="C10")

cytof_data$trans_freq=scale(asin(sqrt(cytof_data$frequency/100)), center = TRUE, scale = TRUE)

wide_cytof <- data.frame(cytof_data %>%
                           select(cluster_id, sample_id, trans_freq) %>%
                           pivot_wider(names_from = sample_id, values_from = trans_freq))

rownames(wide_cytof) <- gsub("ï", "i", wide_cytof$cluster_id)
wide_cytof <- as.matrix(select(wide_cytof, -cluster_id))
class(wide_cytof) <- "double"

### PLASMA DATA ####


# data <- read.csv("C:/Users/Florian/PhD/oxford/vac69/plasma_analytes_vivax_no_inequalitites.csv")
#dataplasma <- read.csv("/Users/s1249052/PhD/plasma/vac69a/Vivax_plasma_analytes2_no_inequalities.csv")
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

log_plasma_data <- rbind(log_plasma_data, "alt"=biochem_data["alt",!grepl("V05", colnames(biochem_data))])


# log_plasma_data <- rbind(log_plasma_data, substr(colnames(log_plasma_data), 1,3), substr(colnames(log_plasma_data), 5,nchar(colnames(log_plasma_data))))
# rownames(log_plasma_data)[(nrow(log_plasma_data)-1):nrow(log_plasma_data)] <- c("Volunteer", "Timepoint")


###### clinical & biochem data


## biochem ####

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
