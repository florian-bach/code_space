library(CATALYST)
library(tidyr)
library(dplyr)
set.seed(1234)

`%notin%` <- Negate(`%in%`)


hard_downsample <- function(fs, event_number){
  flowCore::fsApply(fs, function(ff){
    idx <- sample.int(nrow(ff), min(event_number, nrow(ff)))
    ff[idx,]
  })
}

# fcs_files <- list.files("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/Aurora Spectral Flow/For Arefin/fcs/", pattern = ".fcs")
# full_paths <-  list.files("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/Aurora Spectral Flow/For Arefin/fcs", pattern = ".fcs", full.names = TRUE)

fcs_files <- list.files("~/postdoc/stanford/cytometry/spectral/MUSICAL/lrs_all_leukocytes/", pattern = ".fcs")
full_paths <-  list.files("~/postdoc/stanford/cytometry/spectral/MUSICAL/lrs_all_leukocytes/", pattern = ".fcs", full.names = TRUE)

musical_panel <- readxl::read_excel("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/Aurora Spectral Flow/For Arefin/MUSICALSpectralFlowPanel.xlsx")
musical_panel$class <- "type"
musical_panel$Fluorophore <- paste(musical_panel$Fluorophore, "A", sep="-")

kylie_metadata <- readxl::read_excel("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/Aurora Spectral Flow/For Arefin/MUSICALSpectralFlowMetadata.xlsx")
colnames(kylie_metadata)[1] <- "file_name" 

kylie_metadata$live_name <- paste(kylie_metadata$experiment, kylie_metadata$file_name, sep="_")
# kylie_metadata$live_name <- gsub("mus", "livecells_mus", kylie_metadata$live_name, fixed = T)
kylie_metadata$live_name <- gsub(".fcs", "_Live Cells.fcs", kylie_metadata$live_name, fixed = T)

# malaria_metadata <- kylie_metadata%>%
#   filter(live_name %in% fcs_files)%>%
#   filter(!grepl("lrs", id), !is.na(file_name), )

# weird_files <- c("mus6_D4 Well_040_Live Cells.fcs",
#                  "mus6_D5 Well_041_Live Cells.fcs",
#                  "mus7_D8 Well_044_Live Cells.fcs",
#                  "mus7_D9 Well_045_Live Cells.fcs",
#                  "mus7_D10 Well_046_Live Cells.fcs")

metadata_to_read <- kylie_metadata%>%
  # filter(live_name %notin% weird_files)%>%
  filter(grepl("lrs", id), !is.na(file_name), experiment %notin% c("mus5"))

metadata_to_read$file_path <- full_paths[match(metadata_to_read$live_name, fcs_files)]
# set.seed(1234)


musical_flowset <- flowCore::read.flowSet(files = metadata_to_read$file_path, truncate_max_range = FALSE)
# # set.seed(1234)
# musical_flowset <- hard_downsample(musical_flowset, event_number = 10000)
# 
# co_factors <- read.csv("~/postdoc/stanford/cytometry/spectral/MUSICAL/trans_coefs.csv")
# 
# trans_musical_flowset <- flowSpecs::arcTrans(musical_flowset,  transCoFacs = co_factors$coef, transNames = colnames(musical_flowset)[8:38])

sce <- prepData(musical_flowset,
                FACS = T,
                musical_panel,
                metadata_to_read,
                transform = FALSE,
                md_cols =
                  list(file = "live_name",
                       id = "live_name",
                       factors = c("id", "timepoint", "experiment")),
                panel_cols =
                  list(channel = "Fluorophore",
                       antigen = "Marker",
                       class = "class"))

assay(sce, "exprs") <- assay(sce, "counts")

cluster_markers <- musical_panel$Marker[-match(c("IgD", "TCR VD1", "LiveDead", "CD45"), musical_panel$Marker)]

# plotExprs(sce, features = cluster_markers, color_by = "id")


sce <- CATALYST::cluster(sce,
                         features = cluster_markers,
                         xdim = 12, ydim = 12, maxK = 50, seed = 1234)


plotting_markers <- c("CD19", "CD20", "IgG", "IgM", "IgD", "CD21", "CD27", "HLA DR",  "CD11c", "CD1d", "CD123", "CD33", "CD14", "CD16", "CD56", "CD7", "CD85j", "CD24", "CD2", "TCR VD1", "TCR VD2", "CD57", "CD127", "CD3", "CD8", "CCR7CD197", "CD45RA")

(unclean_heatmap <- plotExprHeatmap(sce,
                                    by="cluster_id",
                                    assay = "exprs",
                                    fun="median",
                                    k="meta35",
                                    row_anno = TRUE,
                                    bars = T,
                                    perc = TRUE,
                                    scale = "first",
                                    q=0.05,
                                    features=plotting_markers[-c(20)],
                                    row_clust = TRUE,
                                    hm_pal = viridis::inferno(n=5),
                                    col_clust = FALSE))

png("~/postdoc/stanford/cytometry/spectral/MUSICAL/omiq_normed/no_correction_clustered_heatmap_meta35.png", width=12, height=8, units="in", res=444)
ComplexHeatmap::draw(unclean_heatmap)
dev.off()

first_merge <- readxl::read_excel("~/postdoc/stanford/cytometry/spectral/MUSICAL/lrs_all_leukocytes/first_merge.xlsx", sheet = 2)
first_merge$new <- gsub("__", "", first_merge$new)
first_merge$new <- gsub("^_", "", first_merge$new)
first_merge$new <- gsub("___", "_", first_merge$new)

sce <- mergeClusters(sce, k="meta40", first_merge, "first_merge")


(first_merge_heatmap <- plotExprHeatmap(sce,
                                    by="cluster_id",
                                    assay = "exprs",
                                    fun="median",
                                    k="first_merge",
                                    row_anno = TRUE,
                                    bars = T,
                                    perc = TRUE,
                                    scale = "first",
                                    q=0.05,
                                    features=plotting_markers[-c(20)],
                                    row_clust = TRUE,
                                    hm_pal = viridis::inferno(n=5),
                                    col_clust = FALSE))

png("~/postdoc/stanford/cytometry/spectral/MUSICAL/omiq_normed/first_merge_heatmap.png", width=12, height=7, units="in", res=444)
ComplexHeatmap::draw(first_merge_heatmap)
dev.off()
