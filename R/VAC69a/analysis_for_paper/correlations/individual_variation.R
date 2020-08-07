library(ggplot2)
library(dplyr)
library(tidyr)
my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")
time_col=colorspace::sequential_hcl(5, palette = "Purple Yellow")

time_palette <- c("Baseline"=time_col[4], "C10"=time_col[3], "DoD"=time_col[2], "T6"=time_col[1])

# CyTOF MDS ####

# cluster frequency

cytof_data <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/cluster_counts_and_freqs.csv", header = T, stringsAsFactors = F)
#cytof_data <- subset(cytof_data, cytof_data$timepoint!="C10")

#cytof_data$trans_freq=scale(asin(sqrt(cytof_data$frequency/100)), center = TRUE, scale = TRUE)

cytof_data$trans_freq=asin(sqrt(cytof_data$frequency/100))


# wide_cytof <- data.frame(cytof_data %>%
#                            select(cluster_id, sample_id, frequency) %>%
#                            pivot_wider(names_from = sample_id, values_from = frequency))

wide_cytof <- data.frame(cytof_data %>%
                           select(cluster_id, sample_id, trans_freq) %>%
                           pivot_wider(names_from = sample_id, values_from = trans_freq))


rownames(wide_cytof) <- gsub("ï", "i", wide_cytof$cluster_id)
wide_cytof <- as.matrix(select(wide_cytof, -cluster_id))
class(wide_cytof) <- "double"



cytof_mds <- data.frame(cmdscale(robCompositions::aDist(t(wide_cytof))))

colnames(cytof_mds) <- c("MDS1", "MDS2")

cytof_mds$Sample_ID <- rownames(cytof_mds)
cytof_mds$Timepoint <- substr(cytof_mds$Sample_ID, 5, nchar(cytof_mds$Sample_ID))
cytof_mds$Volunteer <- substr(cytof_mds$Sample_ID, 1, 3)



aitchison_cytof <- ggplot(cytof_mds, aes(x=MDS1, y=MDS2))+
  geom_point(aes(shape=Timepoint, colour=Volunteer))+
  theme_minimal()+
  scale_color_manual(values = my_paired_palette)
                      
ggsave("~/PhD/multi_omics/Aitchison's_CyTOF.png", aitchison_cytof)     


#median marker expression

merged_daf <- vac69a.cytof::read_full("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")

cs_by_s <- split(seq_len(ncol(merged_daf)), merged_daf$sample_id)
es <- as.matrix(SummarizedExperiment::assay(merged_daf, "exprs"))
  ms <- vapply(cs_by_s, function(cs) Biobase::rowMedians(es[, cs, drop = FALSE]), 
               numeric(nrow(merged_daf)))
  rownames(ms) <- rownames(merged_daf)
  # state markers, type markers or both
  ms <- subset(ms, rownames(ms)%in%c(c(CATALYST::state_markers(merged_daf)[-5], "Perforin")))
  
  mds <- limma::plotMDS(ms, plot = FALSE)
  df <- data.frame(MDS1 = mds$x, MDS2 = mds$y)
  md <- S4Vectors::metadata(merged_daf)$experiment_info
  m <- match(rownames(df), md$sample_id)
  df <- data.frame(df, md[m, ])
  
  #df2 <- filter(df, timepoint %in% c("Baseline", "DoD"))
  
  state_markers_mds <- ggplot(df, aes(x=MDS1, y=MDS2, color=volunteer))+
    geom_point(aes(shape=timepoint))+
    #ggrepel::geom_label_repel(aes_string(label = "sample_id"), show.legend = FALSE)+ 
    theme_minimal()+
    scale_color_manual(values = my_paired_palette)
  
  ggsave("~/PhD/multi_omics/state_markers_mds.png", state_markers_mds)
    # 
  # ggplot()+
  # geom_point(df, aes_(x=df$MDS1, y=df$MDS2, color=df$volunteer))+  
  # geom_segment(aes_(x=df2$MDS1[grepl("Baseline", df2$timepoint)],
  #                  xend=df2$MDS1[grepl("DoD", df2$timepoint)],
  #                  y=df2$MDS2[grepl("Baseline", df2$timepoint)],
  #                  yend=df2$MDS2[grepl("DoD", df2$timepoint)]),
  #                  arrow =arrow(length = unit(0.2, "cm")))+





# plasma mds ####
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
  
  log_plasma_data <- log2(plasma_data)

# sig_analytes <- scan("~/PhD/plasma/vac69a/analytes_sorted_by_padj.txt", what="", skip = 1)[1:12]
#   
# sig_log_plasma_data <- subset(log_plasma_data, rownames(log_plasma_data)%in%sig_analytes)
#   
plasma_mds <- limma::plotMDS(log_plasma_data, plot = FALSE)

plasma_df <- data.frame(MDS1 = plasma_mds$x, MDS2 = plasma_mds$y)

plasma_df$Sample_ID <- rownames(plasma_df)
plasma_df$Timepoint <- substr(plasma_df$Sample_ID, 5, nchar(plasma_df$Sample_ID))
plasma_df$Volunteer <- substr(plasma_df$Sample_ID, 1, 3)


plasma_mds <- ggplot(plasma_df, aes(x=MDS1, y=MDS2, color=Volunteer))+
  geom_point(aes(shape=Timepoint))+
  #ggrepel::geom_label_repel(aes_string(label = "sample_id"), show.legend = FALSE)+ 
  theme_minimal()+
  scale_color_manual(values = my_paired_palette)


ggsave("~/PhD/multi_omics/plasma_mds.png", plasma_mds)


sig_analytes <- scan("~/PhD/plasma/vac69a/analytes_sorted_by_padj.txt", what="", skip = 1)[1:12]

