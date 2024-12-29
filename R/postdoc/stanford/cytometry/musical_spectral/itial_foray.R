library(CATALYST)
library(tidyr)
library(dplyr)


hard_downsample <- function(fs, event_number){
  flowCore::fsApply(fs, function(ff){
    idx <- sample.int(nrow(ff), min(event_number, nrow(ff)))
    ff[idx,]
  })
}

fcs_files <- list.files("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/Aurora Spectral Flow/For Arefin/FCSFilesLiveSinglets", pattern = ".fcs")
full_paths <-  list.files("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/Aurora Spectral Flow/For Arefin/FCSFilesLiveSinglets", pattern = ".fcs", full.names = TRUE)

musical_panel <- readxl::read_excel("~/Downloads/MUSICALSpectralFlowPanel.xlsx")
musical_panel$class <- "type"
musical_panel$Fluorophore <- paste(musical_panel$Fluorophore, "A", sep="-")

kylie_metadata <- readxl::read_excel("~/Downloads/MUSICALSpectralFlowMetadata.xlsx")
colnames(kylie_metadata)[1] <- "file_name" 

kylie_metadata$live_name <- paste(kylie_metadata$experiment, kylie_metadata$file_name, sep="_")
# kylie_metadata$live_name <- gsub("mus", "livecells_mus", kylie_metadata$live_name, fixed = T)
kylie_metadata$live_name <- gsub(".fcs", "_Live Cells.fcs", kylie_metadata$live_name, fixed = T)

malaria_metadata <- kylie_metadata%>%
  filter(live_name %in% fcs_files)%>%
  filter(!grepl("lrs", id), !is.na(file_name))

set.seed(1234)
metadata_to_read <- malaria_metadata[sample(1:159, 85, replace = F),]

metadata_to_read$file_path <- full_paths[match(metadata_to_read$live_name, fcs_files)]
# set.seed(1234)


musical_flowset <- flowCore::read.flowSet(files = metadata_to_read$file_path, truncate_max_range = FALSE)
set.seed(1234)
musical_flowset <- hard_downsample(musical_flowset, event_number = 10000)

co_factors <- read.csv("~/postdoc/stanford/cytometry/spectral/MUSICAL/trans_coefs.csv")

trans_musical_flowset <- flowSpecs::arcTrans(musical_flowset,  transCoFacs = co_factors$coef, transNames = colnames(musical_flowset)[8:38])

sce <- prepData(trans_musical_flowset,
                FACS = T,
                musical_panel,
                metadata_to_read,
                transform = FALSE,
                md_cols =
                  list(file = "live_name",
                       id = "file_name",
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
                                    k="meta40",
                                    row_anno = TRUE,
                                    bars = T,
                                    perc = TRUE,
                                    scale = "first",
                                    q=0.05,
                                    features=plotting_markers[-c(5, 20)],
                                    row_clust = TRUE,
                                    hm_pal = viridis::inferno(n=5),
                                    col_clust = FALSE))



slim_sce <- filterSCE(sce, cluster_id %notin% c(14), k="meta40")


