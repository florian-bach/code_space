# CATALYST ####
library(CATALYST)

`%notin%` <- Negate(`%in%`)


hard_downsample <- function(fs, event_number){
  flowCore::fsApply(fs, function(ff){
    idx <- sample.int(nrow(ff), min(event_number, nrow(ff)))
    ff[idx,]
  })
}

musical_panel <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/musical_panel_edit.csv")

musical_metadata <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/peacoQC_metadata.csv")
musical_metadata$batch <- as.character(musical_metadata$batch)

# long_file_list <- list.files("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/peacoQC/PeacoQC_results/fcs_files/", full.names = TRUE)
# long_file_list <- list.files("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/big_fcs/", full.names = TRUE)
# files_to_read <- subset(musical_metadata, subject_id==1217, select=file_path)

metadata_to_read <- subset(musical_metadata, class ==  %in% c(268, 324, 137, 176, 353, 161, 363, 571))

musical_flowset <- read.ncdfFlowSet(metadata_to_read$file_path)

set.seed(1234)
downsampled_musical_flowset <- musical_flowset
# downsampled_musical_flowset <- musical_flowset

musical_flowset <- NULL
gc(full = TRUE)

sce <- prepData(downsampled_musical_flowset,
                musical_panel,
                metadata_to_read,
                transform = TRUE,
                md_cols =
                  list(file = "file_name",
                       id = "sample_id",
                       factors = c("subject_id", "timepoint", "class", "batch")))

downsampled_musical_flowset <- NULL
gc(full = TRUE)

# p <- plotExprs(sce, color_by = "batch")
# p$facet$params$ncol <- 6
# p

cluster_markers <- musical_panel$antigen[-c(1,2,5,40:42, 44:48)]
type_markers <- musical_panel$antigen[musical_panel$class=="type"&!is.na(musical_panel$class)]

# sce <- filterSCE(sce, batch!=41824)

sce <- CATALYST::cluster(sce, features = cluster_markers, xdim = 12, ydim = 12, maxK = 50)

# plotCounts(sce, group_by = "sample_id", color_by = "class")



big_heat <- plotMultiHeatmap(live_sce,
                             k = "meta30",
                             m="meta30",
                             hm1=type_markers,
                             # hm2 = NA,
                             row_clust = FALSE,
                             col_clust = FALSE,
                             # col_dend = TRUE,
                             scale = "first")

png("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/figures/big_heat.png", width=20, height=12, units = "in", res=400)
ComplexHeatmap::draw(big_heat)
dev.off()

# 
# system.time(sce <- runDR(sce,
#                          dr="UMAP",
#                          cells=1000,
#                          features=cluster_markers)
# )
# 
# plotDR(sce,  color_by = "CD3", scale = TRUE)
# plotDR(sce,  color_by = "Dead", scale = TRUE)
# 
# ggplot2::ggsave("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/marker_umap.png", marker_umap, height = 10, width=10, bg="white")

# filter dead ####

# live_sce <- filterSCE(sce, cluster_id %notin% c(27, 25, 19, 11:13, 2), k="meta30")

# differential abundance ####
library(diffcyt)
  ## kids only ####


ei <- metadata(kids_only)$experiment_info

colnames(design)
# [1] "(Intercept)"      "timepointday 0"   "timepointday 14"  "timepointday 7"   "classsymptomatic" "subject_id161"   
# [7] "subject_id176"    "subject_id268"    "subject_id324"    "subject_id353"    "subject_id363"    "subject_id571"  
#timepoint=baseline is dummy; 
design <- createDesignMatrix(ei, c("timepoint", "class", "subject_id"))

FDR_cutoff <- 0.05


baseline_w_o_symp <- createContrast(c(c(0, 0, 0, 0, 1), rep(0,7)))
d14_w_o_symp <- createContrast(c(c(0, 0, 1, 0, 1), rep(0,7)))

# pairwise_contrast_dod <- createContrast(c(c(0, 0, 1, 0), rep(0,5)))
# pairwise_contrast_c10 <- createContrast(c(c(0, 1, 0, 0), rep(0,5)))

# pairwise_contrast_dod_t6 <- createContrast(c(c(0, 0, -1, 1), rep(0,5)))
# pairwise_contrast_t6 <- createContrast(c(c(0, 0, 0, 1), rep(0,5), 1))
# pairwise_contrast_dod <- createContrast(c(c(0, 0, 1, 0), rep(0,5), 1))

#pairwise_contrast_t6 <- createContrast(c(0, 0, 1, 1))



da_baseline_class <- diffcyt(kids_only,
                  design = design,
                  contrast = baseline_w_o_symp,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "meta35",
                  verbose = T)


d14_w_o_symp <- createContrast(c(c(0, 0, 1, 0, 1), rep(0,7)))

base_d14_class <- diffcyt(kids_only,
                  design = design,
                  contrast = d14_w_o_symp,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "meta35",
                  verbose = T)


table(rowData(da_baseline_class$res)$p_adj < FDR_cutoff)
# FALSE
#    42
table(rowData(base_d14_class$res)$p_adj < FDR_cutoff)
# # FALSEFALSE  TRUE
# 34     8



# with C+10 being the intercept term, c(-1, 1) contrasts can be made between baseline and dod/t6
# the significant clusters at dod are 10, 12, 14, 37, 3, 6, 43, 45, (n=8) p_adjust range from 2.5e-6 to 2.3e-2
# the significant clusters at t6 are 37, 15, 43, 10, 36, 12, 42, 14, 6, 8, 3, 29, 18, 19 (n=14), p+adjust range from 5.2e-52 to 4.4e-2

# these are exactly the same clusters, when baseline is the intercept term and we set it to 0 for the contrasts
# this is good news!!!

# using a design matrix with timepoint and batch reduces our dod clusters to 0 and the t6 clusters to 6 (37, 43, 36, 15, 42, 12)
# removing the cluster term returns 0 significant clusters at dod and 5 at t6 (37, 32, 36, 15, 42) so there's a small effect;
# all the mismatched clusters between those models are detected in the full design matrix
da <- rowData(da_baseline_class$res)
plotDiffHeatmap(kids_only, da, all=TRUE)


da <- rowData(base_d14_class$res)
plotDiffHeatmap(kids_only, da, all=TRUE)

plotDiffHeatmap(merged_daf, da_dod, th = FDR_cutoff, normalize = TRUE, hm1 = F, top_n = 10)
plotDiffHeatmap(merged_daf, da_t6, th = FDR_cutoff, normalize = TRUE, hm1 = F, top_n = 30)




# how accesory files were made #### 

  ## make new metadata ####
file_list <- list.files("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_big_fcs/")
long_file_list <- list.files("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_big_fcs/", full.names = TRUE)

batch <- paste("batch_", substr(file_list, 0, 6), sep="")


short_file_names <- substr(file_list, 39, nchar(file_list))
short_file_names <- gsub("final-", "", short_file_names)

subject_id <- substr(short_file_names, 1, 4)
subject_id <- gsub("_", "", subject_id)

timepoint <- substr(short_file_names, 7, 11)
timepoint <- dplyr::case_when(timepoint=="14_01"~"day 14",
                              timepoint=="7_01_"~"day 7",
                              timepoint=="0_01_"~"day 0",
                              timepoint=="basel"~"baseline",
                              timepoint=="mune_"~"immune",
                              timepoint=="IMCct"~"control")

class <- substr(short_file_names, 4, 6)

class <- dplyr::case_when(class=="7_H"~"control",
                          class=="_A_"~"asympmtomatic",
                          class=="_S_"~"symptomatic",
                          class=="_im"~"immune")

cytof_metadata <- data.frame("batch"=as.character(batch),
                             "subject_id"=subject_id,
                             "timepoint"=timepoint,
                             "class"=class,
                             "file_name"=file_list,
                             "file_path"=long_file_list,
                             "sample_id"=paste(batch, subject_id, class, timepoint, sep="_"))

write.csv(cytof_metadata, "~/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_raw_metadata.csv", row.names = FALSE)

  ## make old metadata ####
file_list <- list.files("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/sandbox_out/peacoQC/PeacoQC_results/fcs_files/", pattern = ".fcs$")
long_file_list <- list.files("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/sandbox_out/peacoQC/PeacoQC_results/fcs_files/", pattern = ".fcs$", full.names = TRUE)

batch <- substr(file_list, 0, 6)


short_file_names <- substr(file_list, 39, nchar(file_list))
short_file_names <- gsub("final-", "", short_file_names)

subject_id <- substr(short_file_names, 1, 4)
subject_id <- gsub("_", "", subject_id)

timepoint <- substr(short_file_names, 7, 11)
timepoint <- dplyr::case_when(timepoint=="14_01"~"day 14",
                              timepoint=="7_01_"~"day 7",
                              timepoint=="0_01_"~"day 0",
                              timepoint=="basel"~"baseline",
                              timepoint=="mune_"~"immune",
                              timepoint=="IMCct"~"control")

class <- substr(short_file_names, 4, 6)

class <- dplyr::case_when(class=="7_H"~"control",
                          class=="_A_"~"asympmtomatic",
                          class=="_S_"~"symptomatic",
                          class=="_im"~"immune")

cytof_metadata <- data.frame("batch"=as.character(batch),
                             "subject_id"=subject_id,
                             "timepoint"=timepoint,
                             "class"=class,
                             "file_name"=file_list,
                             "file_path"=long_file_list,
                             "sample_id"=paste(batch, subject_id, timepoint, sep="_"))

write.csv(cytof_metadata, "~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/peacoqc_metadata.csv", row.names = FALSE)


  ## get panel ####
raw_musical_panel <- data.frame(premessa::read_parameters(long_file_list[1]))

musical_panel <- data.frame("fcs_colname"=rownames(raw_musical_panel),
                            "antigen"= substr(raw_musical_panel[,1], 7, 12))

musical_panel$antigen <- ifelse(nchar(musical_panel$antigen)<3, musical_panel$fcs_colname, musical_panel$antigen)
musical_panel$antigen <- gsub("length", "Event_length", musical_panel$antigen)

musical_panel$class <- "type"

write.csv(musical_panel, "~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/musical_panel.csv", row.names = FALSE)

