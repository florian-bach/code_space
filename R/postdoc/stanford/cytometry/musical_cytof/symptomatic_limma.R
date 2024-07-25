# CATALYST ####
library(CATALYST)

`%notin%` <- Negate(`%in%`)

kids_palette <- colorspace::sequential_hcl(8, palette = "Emrld")
names(kids_palette) <- as.character(c(268, 324, 137, 176, 353, 161, 363, 571))

adults_palette <- colorspace::sequential_hcl(5, palette = "OrRd")[1:4]
names(adults_palette) <- as.character(c(575, 627, 563, 643))


hard_downsample <- function(fs, event_number){
  flowCore::fsApply(fs, function(ff){
    idx <- sample.int(nrow(ff), min(event_number, nrow(ff)))
    ff[idx,]
  })
}

musical_panel <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/musical_panel_edit.csv")[,2:4]

musical_metadata <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/single_cell_metadata.csv")
musical_metadata$batch <- as.character(musical_metadata$batch)

# long_file_list <- list.files("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/peacoQC/PeacoQC_results/fcs_files/", full.names = TRUE)
# long_file_list <- list.files("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/big_fcs/", full.names = TRUE)
# files_to_read <- subset(musical_metadata, subject_id==1217, select=file_path)

metadata_to_read <- subset(musical_metadata, class == "symptomatic")

musical_flowset <- ncdfFlow::read.ncdfFlowSet(metadata_to_read$file_path)
set.seed(1234)
musical_flowset <- hard_downsample(musical_flowset, event_number = 100000)



sce <- prepData(musical_flowset,
                musical_panel,
                metadata_to_read,
                transform = FALSE,
                md_cols =
                  list(file = "file_name",
                       id = "sample_id",
                       factors = c("subject_id", "timepoint", "class", "batch")),
                panel_cols =
                  list(channel = "fcs_colname",
                       antigen = "antigen",
                       class = "class"))


musical_flowset <- NULL
gc(full = TRUE)

assay(sce, "exprs") <- assay(sce, "counts")
# p <- plotExprs(sce, color_by = "batch", features = "Dead")
# p$facet$params$ncol <- 6
# p

cluster_markers <- musical_panel$antigen[-c(1,2,4:5,40:42, 44:48)]
type_markers <- musical_panel$antigen[musical_panel$class=="type"&!is.na(musical_panel$class)]


cluster_markers <- c("CD3",
                     "CD4",
                     "CD8",
                     "TCRgd",
                     "Vd2TCR",
                     "CD161",
                     "CD56",
                     "CD94",
                     "CD19",
                     "CD20",
                     "CD21",
                     "IgD",
                     "CD38",
                     "CD27",
                     "CD25",
                     "CD127",
                     "CD39",
                     "PD-1",
                     "CXCR5",
                     "ICOS",
                     "CCR6",
                     "CXCR3",
                     "CD49b",
                     "HLADR",
                     "CD14",
                     "CD16",
                     "CD85j",
                     "CD11b",
                     "CD11c",
                     "CD123",
                     "CD33",
                     "CD86",
                     "LAG3",
                     "CD57",
                     "CD45RA",
                     "CCR7")

# sce <- filterSCE(sce, batch!=41824)

sce <- CATALYST::cluster(sce,
                         features = cluster_markers,
                         xdim = 12, ydim = 12, maxK = 50, seed = 1234)

# phenotypic heatmaps & UMAPs####
##remove putative doublets ####
# sce <- filterSCE(sce, cluster_id %notin% c(25, 22), k="meta45")
# 

complex_merge <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/complex_merge.csv", header = TRUE)
simple_merge <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/simple_merge.csv", header = TRUE)

sce <- mergeClusters(sce, k="meta50", complex_merge, "complex_merge", overwrite = TRUE)
sce <- mergeClusters(sce, k="meta50", simple_merge, "simple_merge", overwrite = TRUE)

library(diffcyt)

ei <- metadata(sce)$experiment_info
ei$age_class <- ifelse(ei$subject_id %in% c(268, 324, 137, 176, 353, 161, 363, 571), "child", "adult")
FDR_cutoff <- 0.1

design <- createDesignMatrix(ei, c("timepoint", "batch"))


colnames(design)
# [1] "(Intercept)"     "timepointday 0"  "timepointday 14" "timepointday 7"  "subject_id161"   "subject_id176"  
# [7] "subject_id268"   "subject_id324"   "subject_id353"   "subject_id363"   "subject_id563"   "subject_id571"  
# [13] "subject_id575"   "subject_id627"   "subject_id643" 
#timepoint=baseline is dummy; 

base_0 <-    createContrast(c(c(0, 1, 0, 0), rep(0,4)))
base_14 <- createContrast(c(c(0, 0, 1, 0), rep(0,4)))
base_7 <-  createContrast(c(c(0, 0, 0, 1), rep(0,4)))

zero_7 <-  createContrast(c(c(0, -1, 0, 1), rep(0,4)))
zero_14 <-  createContrast(c(c(0, 1, 1, 0), rep(0,4)))

cluster_markers <- musical_panel$antigen[-c(1,2,4:5,40:42, 44:48)]
type_markers <- musical_panel$antigen[musical_panel$class=="type"&!is.na(musical_panel$class)]
state_markers <- subset(cluster_markers, cluster_markers %notin%type_markers)
state_logic <- musical_panel$antigen %in% state_markers

# state_markers(sce) <- sate_markers
# marker_classes(sce) <- c("state",
#                          "none",
#                          "none")
ds_base_0 <- diffcyt(sce,
                     design = design,
                     contrast = base_0,
                     # markers_to_test = rep(TRUE, 17),
                     analysis_type = "DS",
                     method_DS = "diffcyt-DS-limma",
                     clustering_to_use = "simple_merge",
                     verbose = T)

ds <- rowData(ds_base_0$res)

sce$sample_id <- factor(sce$sample_id, levels=metadata_to_read$sample_id[order(metadata_to_read$timepoint)])

map <- plotDiffHeatmap(sce, ds)
map@matrix <- map@matrix[,order(substr(colnames(map@matrix), nchar(colnames(map@matrix))-15, nchar(colnames(map@matrix))))]

png("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/symptomatic_only/symptomatic_only_time_batch_limma.png", width=15, height=10, units = "in", res=400)
ComplexHeatmap::draw(map)
dev.off()

