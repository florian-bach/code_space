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

metadata_to_read <- subset(musical_metadata, subject_id %in% c(268, 324, 137, 176, 353, 161, 363, 571))

musical_flowset <- ncdfFlow::read.ncdfFlowSet(metadata_to_read$file_path)

musical_flowset <- hard_downsample(musical_flowset, event_number = 100000)



sce <- prepData(musical_flowset,
                musical_panel,
                metadata_to_read,
                transform = TRUE,
                md_cols =
                  list(file = "file_name",
                       id = "sample_id",
                       factors = c("subject_id", "timepoint", "class", "batch")))

musical_flowset <- NULL
gc(full = TRUE)

# p <- plotExprs(sce, color_by = "batch", features = "Dead")
# p$facet$params$ncol <- 6
# p

cluster_markers <- musical_panel$antigen[-c(1,2,4:5,40:42, 44:48)]
type_markers <- musical_panel$antigen[musical_panel$class=="type"&!is.na(musical_panel$class)]

# sce <- filterSCE(sce, batch!=41824)
set.seed(1234)
sce <- CATALYST::cluster(sce, features = cluster_markers, xdim = 12, ydim = 12, maxK = 50)

# plotCounts(sce, group_by = "sample_id", color_by = "class")



big_heat <- plotMultiHeatmap(sce,
                             k = "meta45",
                             # m="meta30",
                             hm1=cluster_markers,
                             hm2 = "abundances",
                             row_clust = TRUE,
                             col_clust = TRUE)

png("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/sandbox_out/peacoQC/peacoq_big_heat.png", width=20, height=12, units = "in", res=400)
ComplexHeatmap::draw(big_heat)
dev.off()


system.time(sce <- runDR(sce,
                         dr="UMAP",
                         cells=5000,
                         features=cluster_markers)
)

# plotDR(sce,  color_by = "CD3", scale = TRUE)
# plotDR(sce,  color_by = "Dead", scale = TRUE)

marker_umap <- plotDR(sce,  color_by = type_markers, scale = TRUE)

ggplot2::ggsave("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/marker_umap.png", marker_umap, height = 20, width=20, bg="white")

# filter dead ####
sce <- filterSCE(sce, cluster_id %notin% c(1, 2, 4, 9, 13, 22, 28, 31, 37), k="meta45")
sce <- CATALYST::cluster(sce, features = cluster_markers, xdim = 12, ydim = 12, maxK = 50)

# differential abundance ####
library(diffcyt)
## kids only ####


ei <- metadata(sce)$experiment_info

colnames(design)
# [1] "(Intercept)"      "timepointday 0"   "timepointday 14"  "timepointday 7"   "classsymptomatic" "subject_id161"   
# [7] "subject_id176"    "subject_id268"    "subject_id324"    "subject_id353"    "subject_id363"    "subject_id571"  
#timepoint=baseline is dummy; 
design <- createDesignMatrix(ei, c("timepoint", "class", "subject_id"))

FDR_cutoff <- 0.05


baseline_symp <-  createContrast(c(c(0, 0, 0, 1), rep(0,8)))

base_0_symp <-    createContrast(c(c(0, 1, 0, 0, 1), rep(0,7)))
zero_14_o_symp <- createContrast(c(c(0, 1, -1, 0, 1), rep(0,7)))
zero_7_o_symp <-  createContrast(c(c(0, 1, 0, -1, 1), rep(0,7)))

base_0_asymp <-    createContrast(c(c(-1, 1, 0, 0, 0), rep(0,7)))
zero_14_o_asymp <- createContrast(c(c(0, 1, -1, 0, 0), rep(0,7)))
zero_7_o_asymp <-  createContrast(c(c(0, 1, 0, -1, 0), rep(0,7)))


# baseline_asymp <- createContrast(c(c(1, 0, 0, 0, 0), rep(0,7)))

# da_base_asymp <- diffcyt(sce,
#                           design = design,
#                           contrast = baseline_asymp,
#                           analysis_type = "DA",
#                           method_DA = "diffcyt-DA-edgeR",
#                           clustering_to_use = "meta45",
#                           verbose = T)

da_base_0_symp <- diffcyt(sce,
                             design = design,
                             contrast = base_0_symp,
                             analysis_type = "DA",
                             method_DA = "diffcyt-DA-edgeR",
                             clustering_to_use = "meta45",
                             verbose = T)

da_zero_14_o_symp <- diffcyt(sce,
                             design = design,
                             contrast = zero_14_o_symp,
                             analysis_type = "DA",
                             method_DA = "diffcyt-DA-edgeR",
                             clustering_to_use = "meta45",
                             verbose = T)




da_zero_7_o_symp <- diffcyt(sce,
                             design = design,
                             contrast = zero_7_o_symp,
                             analysis_type = "DA",
                             method_DA = "diffcyt-DA-edgeR",
                             clustering_to_use = "meta45",
                             verbose = T)




da_base_0_asymp <- diffcyt(sce,
                          design = design,
                          contrast = base_0_asymp,
                          analysis_type = "DA",
                          method_DA = "diffcyt-DA-edgeR",
                          clustering_to_use = "meta45",
                          verbose = T)

da_zero_14_o_asymp <- diffcyt(sce,
                             design = design,
                             contrast = zero_14_o_asymp,
                             analysis_type = "DA",
                             method_DA = "diffcyt-DA-edgeR",
                             clustering_to_use = "meta45",
                             verbose = T)




da_zero_7_o_asymp <- diffcyt(sce,
                            design = design,
                            contrast = zero_7_o_asymp,
                            analysis_type = "DA",
                            method_DA = "diffcyt-DA-edgeR",
                            clustering_to_use = "meta45",
                            verbose = T)




FDR_cutoff=0.05

table(rowData(da_base_0_symp$res)$p_adj < FDR_cutoff)
# FALSE
#    45
table(rowData(da_zero_14_o_symp$res)$p_adj < FDR_cutoff)
# FALSE
#    45
table(rowData(da_zero_7_o_symp$res)$p_adj < FDR_cutoff)
# FALSE
#    45
table(rowData(da_base_0_asymp$res)$p_adj < FDR_cutoff)
# FALSE
#    45
table(rowData(da_zero_14_o_asymp$res)$p_adj < FDR_cutoff)
# FALSE
#    45
table(rowData(da_zero_7_o_asymp$res)$p_adj < FDR_cutoff)
# FALSE
#    45


da <- rowData(da_base_0_asymp$res)
plotDiffHeatmap(sce, da, all=TRUE)

da <- rowData(da_base_0_symp$res)
plotDiffHeatmap(sce, da, all=TRUE)



  ## visualise diffcyt results ####
library(tidyr)
library(dplyr)
library(ggplot2)

results_table <- data.frame(diffcyt::topTable(da_zero_14_o_symp, all=T, show_counts = TRUE, show_props = TRUE))

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


symp_vs_time <- long_results %>%
  filter(p_adj<0.1)%>%
  ggplot(aes(x=factor(timepoint, levels=c("baseline", "day 0", "day 7", "day 14")), y=freq))+
  geom_point(aes(color=subject_id))+
  geom_line(aes(color=subject_id, group=subject_id))+
  # geom_boxplot(aes(fill=class), outlier.shape = NA)+
  facet_wrap(class~cluster_id, scales="free")+
  theme_minimal()+
  scale_fill_manual(values = list("asymptomatic"="darkgrey", "symptomatic"="darkred"))

ggsave("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/figures/kids_only_cluster_freqs_symp_vs_time.png", symp_vs_time, width=15, height=15, dpi=400, bg="white")


time_vs_symp <- long_results %>%
  filter(timepoint!="day 7")%>%
  ggplot(aes(x=class, y=freq))+
  # geom_point(aes(color=subject_id))+
  # geom_line(aes(color=subject_id, group=subject_id))+
  geom_boxplot(aes(fill=timepoint), outlier.shape = NA)+
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
file_list <- list.files("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_big_fcs/single_cells/", pattern = ".fcs$")
long_file_list <- list.files("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_big_fcs/single_cells/", pattern = ".fcs$", full.names = TRUE)

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

write.csv(cytof_metadata, "~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/single_cell_metadata.csv", row.names = FALSE)


## get panel ####
raw_musical_panel <- data.frame(premessa::read_parameters(long_file_list[1]))

musical_panel <- data.frame("fcs_colname"=rownames(raw_musical_panel),
                            "antigen"= substr(raw_musical_panel[,1], 7, 12))

musical_panel$antigen <- ifelse(nchar(musical_panel$antigen)<3, musical_panel$fcs_colname, musical_panel$antigen)
musical_panel$antigen <- gsub("length", "Event_length", musical_panel$antigen)

musical_panel$class <- "type"

write.csv(musical_panel, "~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/musical_panel.csv", row.names = FALSE)
