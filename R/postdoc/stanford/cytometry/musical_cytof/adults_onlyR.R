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

musical_metadata <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_raw_metadata.csv")
musical_metadata$batch <- as.character(musical_metadata$batch)

# long_file_list <- list.files("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/peacoQC/PeacoQC_results/fcs_files/", full.names = TRUE)
# long_file_list <- list.files("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/big_fcs/", full.names = TRUE)
# files_to_read <- subset(musical_metadata, subject_id==1217, select=file_path)

metadata_to_read <- subset(musical_metadata, subject_id %in% c(174, 487, 575, 627, 643) & timepoint != 7)

musical_flowset <- flowCore::read.flowSet(metadata_to_read$file_path)

musical_flowset <- hard_downsample(musical_flowset, event_number = 111000)



sce <- prepData(musical_flowset,
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
state_markers <- musical_panel$antigen[musical_panel$class=="state"&!is.na(musical_panel$class)]

# sce <- filterSCE(sce, batch!=41824)
set.seed(1234)
sce <- CATALYST::cluster(sce, features = cluster_markers, xdim = 12, ydim = 12, maxK = 50)

# plotCounts(sce, group_by = "sample_id", color_by = "class")



big_heat <- plotMultiHeatmap(sce,
                             k = "meta45",
                             m="meta30",
                             hm1=c(type_markers,state_markers),
                             hm2 = "abundances",
                             row_clust = FALSE,
                             col_clust = TRUE)

png("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/figures/adults_only_big_heat.png", width=20, height=12, units = "in", res=400)
ComplexHeatmap::draw(big_heat)
dev.off()


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


ei <- metadata(sce)$experiment_info

colnames(design)
# [1] "(Intercept)"      "timepointday 0"   "timepointday 14"  "timepointday 7"   "classsymptomatic" "subject_id161"   
# [7] "subject_id176"    "subject_id268"    "subject_id324"    "subject_id353"    "subject_id363"    "subject_id571"  
#timepoint=baseline is dummy; 
design <- createDesignMatrix(ei, c("timepoint", "subject_id"))

FDR_cutoff <- 0.05


baseline_asymp <-  createContrast(c(1, 0, 0, 0, 0))


base_0 <-    createContrast(c(c(0, 1, 0, 0, 0)))
base_14 <-    createContrast(c(c(0, 0, 1, 0, 0)))
zero_14 <-    createContrast(c(c(0, 1, -1, 0, 0)))

da_base_asymp <- diffcyt(sce,
                         design = design,
                         contrast = baseline_asymp,
                         analysis_type = "DA",
                         method_DA = "diffcyt-DA-edgeR",
                         clustering_to_use = "meta45",
                         verbose = T)

# 
da_base_0 <- diffcyt(sce,
                          design = design,
                          contrast = base_0,
                          analysis_type = "DA",
                          method_DA = "diffcyt-DA-edgeR",
                          clustering_to_use = "meta45",
                          verbose = T)
# 
da_base_14 <- diffcyt(sce,
                             design = design,
                             contrast = base_14,
                             analysis_type = "DA",
                             method_DA = "diffcyt-DA-edgeR",
                             clustering_to_use = "meta45",
                             verbose = T)

da_0_14 <- diffcyt(sce,
                      design = design,
                      contrast = zero_14,
                      analysis_type = "DA",
                      method_DA = "diffcyt-DA-edgeR",
                      clustering_to_use = "meta45",
                      verbose = T)

FDR_cutoff=0.05

table(rowData(da_base_0$res)$p_adj < FDR_cutoff)
# FALSE  TRUE 
# 30    13
table(rowData(da_base_14$res)$p_adj < FDR_cutoff)
# FALSE  TRUE 
# 42     1 
table(rowData(da_0_14$res)$p_adj < FDR_cutoff)
# FALSE  TRUE 
# 27    16




plotDiffHeatmap(sce, rowData(da_base_0$res))

plotDiffHeatmap(sce, rowData(da_base_14$res), all=TRUE)

plotDiffHeatmap(sce, rowData(da_0_14$res), all=TRUE)


## visualise diffcyt results ####
library(tidyr)
library(dplyr)
library(ggplot2)

results_table <- data.frame(diffcyt::topTable(da_0_14, all=T, show_counts = TRUE, show_props = TRUE))

long_results <- results_table %>%
  pivot_longer(cols = starts_with("props"), names_to = "sample_prop", values_to = "freq")%>%
  pivot_longer(cols = starts_with("counts"), names_to = "sample_counts", values_to = "count")%>%
  mutate("subject_id"=substr(sample_prop, 20, 22),
         "timepoint"=substr(sample_prop, nchar(sample_prop)-5, nchar(sample_prop)),
         "class"=substr(sample_prop, 24, 30))%>%
  mutate(timepoint=case_when(timepoint=="seline"~"baseline",
                             timepoint=="_day.0"~ "day 0",
                             timepoint=="day.14"~ "day 14",
                             timepoint=="_day.7"~ "day 7"),
         class=case_when(class=="asympmt"~"asymptomatic",
                         class=="symptom"~"symptomatic"))



time_vs_symp <- long_results %>%
  filter(p_adj<0.05)%>%
  ggplot(aes(x=timepoint, y=freq))+
  geom_point(aes(color=subject_id))+
  geom_line(aes(color=subject_id, group=subject_id))+
  # geom_boxplot(aes(fill=timepoint), outlier.shape = NA)+
  facet_wrap(~cluster_id, scales="free")+
  theme_minimal()+
  scale_fill_manual(values = colorspace::sequential_hcl(5, palette = "Purple Yellow")[1:3])

ggsave("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/figures/kids_only_cluster_freqs_time_vs_symp.png", time_vs_symp, width=15, height=15, dpi=400, bg="white")

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
# file_list <- list.files("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/big_fcs/")
# long_file_list <- list.files("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/big_fcs/", full.names = TRUE)
# 
# batch <- substr(file_list, 0, 6)
# 
# 
# short_file_names <- substr(file_list, 39, nchar(file_list))
# short_file_names <- gsub("final-", "", short_file_names)
# 
# subject_id <- substr(short_file_names, 1, 4)
# subject_id <- gsub("_", "", subject_id)
# 
# timepoint <- substr(short_file_names, 7, 11)
# timepoint <- dplyr::case_when(timepoint=="14_01"~"day 14",
#                               timepoint=="7_01_"~"day 7",
#                               timepoint=="0_01_"~"day 0",
#                               timepoint=="basel"~"baseline",
#                               timepoint=="mune_"~"immune",
#                               timepoint=="IMCct"~"control")
# 
# class <- substr(short_file_names, 4, 6)
# 
# class <- dplyr::case_when(class=="7_H"~"control",
#                           class=="_A_"~"asympmtomatic",
#                           class=="_S_"~"symptomatic",
#                           class=="_im"~"immune")
# 
# cytof_metadata <- data.frame("batch"=as.character(batch),
#                              "subject_id"=subject_id,
#                              "timepoint"=timepoint,
#                              "class"=class,
#                              "file_name"=file_list,
#                              "file_path"=long_file_list,
#                              "sample_id"=paste(batch, subject_id, timepoint, sep="_"))
# 
# write.csv(cytof_metadata, "~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/raw_metadata.csv", row.names = FALSE)


## get panel ####
raw_musical_panel <- data.frame(premessa::read_parameters(long_file_list[1]))

musical_panel <- data.frame("fcs_colname"=rownames(raw_musical_panel),
                            "antigen"= substr(raw_musical_panel[,1], 7, 12))

musical_panel$antigen <- ifelse(nchar(musical_panel$antigen)<3, musical_panel$fcs_colname, musical_panel$antigen)
musical_panel$antigen <- gsub("length", "Event_length", musical_panel$antigen)

musical_panel$class <- "type"

write.csv(musical_panel, "~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/musical_panel.csv", row.names = FALSE)

