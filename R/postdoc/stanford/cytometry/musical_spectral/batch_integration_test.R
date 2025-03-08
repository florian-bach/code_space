library(CATALYST)
library(tidyr)
library(dplyr)
set.seed(1234)

`%notin%` <- Negate(`%in%`)


# changes for cyCombine
# musical_panel$Fluorophore <- paste(musical_panel$Fluorophore, "-A | ", musical_panel$Marker, "-norm", sep="")

# changes for fdaNorm
# musical_panel$Fluorophore <- paste(musical_panel$Marker, "-norm", sep="")

# changes for cytoNorm
# musical_panel$Fluorophore <- paste(musical_panel$Fluorophore, "-A-Norm", sep="")

path_locations <- c("/Users/fbach/postdoc/stanford/cytometry/spectral/MUSICAL/lrs_all_leukocytes/batch_integration_tests/fdaNorm/omiq_exported_data/",
                    "/Users/fbach/postdoc/stanford/cytometry/spectral/MUSICAL/lrs_all_leukocytes/batch_integration_tests/cyCombine/omiq_exported_data/",
                    "/Users/fbach/postdoc/stanford/cytometry/spectral/MUSICAL/lrs_all_leukocytes/batch_integration_tests/cytoNorm/omiq_exported_data/")

musical_panel <- readxl::read_excel("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/Aurora Spectral Flow/For Arefin/MUSICALSpectralFlowPanel.xlsx")

fluorophore_names <- list("fdaNorm"=paste(musical_panel$Marker, "-norm", sep=""),
                          "cycombine"=paste(musical_panel$Fluorophore, "-A | ", musical_panel$Marker, "-norm", sep=""),
                          "cytoNorm"=paste(musical_panel$Fluorophore, "-A-Norm", sep=""))

for(i in 1:3){

path_location <- path_locations[i]
fcs_files <- list.files(path_location, pattern = ".fcs")
full_paths <-  list.files(path_location, pattern = ".fcs", full.names = TRUE)


musical_panel <- mutate(musical_panel, Fluorophore =  fluorophore_names[[i]] )
musical_panel$class <- "type"


kylie_metadata <- readxl::read_excel("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/Aurora Spectral Flow/For Arefin/MUSICALSpectralFlowMetadata.xlsx")
colnames(kylie_metadata)[1] <- "file_name" 

kylie_metadata$live_name <- paste(kylie_metadata$experiment, kylie_metadata$file_name, sep="_")
# kylie_metadata$live_name <- gsub("mus", "livecells_mus", kylie_metadata$live_name, fixed = T)
kylie_metadata$live_name <- gsub(".fcs", "_Live Cells.fcs", kylie_metadata$live_name, fixed = T)



metadata_to_read <- kylie_metadata%>%
  filter(grepl("lrs", id), !is.na(file_name), experiment %notin% c("mus5"))

metadata_to_read$file_path <- full_paths[match(metadata_to_read$live_name, fcs_files)]

musical_flowset <- flowCore::read.flowSet(files = metadata_to_read$file_path, truncate_max_range = FALSE)



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

cluster_markers <- musical_panel$Marker[-match(c("IgD", "TCR VD1"), musical_panel$Marker)]

# plotExprs(sce, features = cluster_markers, color_by = "id")


sce <- CATALYST::cluster(sce,
                         features = cluster_markers,
                         xdim = 12, ydim = 12, maxK = 50, seed = 1234)


plotting_markers <- c("CD19", "CD20", "IgG", "IgM", "IgD", "CD21", "CD27", "HLA DR",  "CD11c", "CD1d", "CD123", "CD33", "CD14", "CD16", "CD56", "CD7", "CD85j", "CD24", "CD2", "TCR VD1", "TCR VD2", "CD57", "CD127", "CD3", "CD8", "CCR7CD197", "CD45RA")

unclean_heatmap <- plotExprHeatmap(sce,
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
                                    col_clust = FALSE)

png(paste(path_location, names(fluorophore_names)[i], "_meta35_heatmap.png", sep=""), width=12, height=8, units="in", res=444)
ComplexHeatmap::draw(unclean_heatmap)
dev.off()

print(paste(names(fluorophore_names)[i], "done"))
}






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





kylie_sheet <- read.csv("~/Downloads/raw spectral single stain panel - Sheet1.csv")

bad_names <- read.csv("~/Downloads/bad_names.csv", check.names = F)

long_bad_names <- bad_names%>%
  pivot_longer(cols=2:72, names_to = "column_name", values_to = "laser_name")%>%
  mutate(marker_name=kylie_sheet$Marker[match(laser_name, kylie_sheet$Laser.name)],
         channel_name=kylie_sheet$Fluorophore[match(laser_name, kylie_sheet$Laser.name)])%>%
  mutate(marker_name=if_else(is.na(marker_name), laser_name, marker_name))%>%
  mutate(channel_name=if_else(is.na(channel_name), laser_name, channel_name))


good_primary_names <- pivot_wider(long_bad_names, names_from = column_name, values_from = channel_name, id_cols = "Filename")
# colnames(good_primary_names)[2:72] <- good_primary_names[2,2:72]

good_secondary_names <- pivot_wider(long_bad_names, names_from = column_name, values_from = marker_name, id_cols = "Filename")
# colnames(good_secondary_names)[2:72] <- good_primary_names[2,2:72]

write.csv(good_primary_names, "~/Downloads/good_primary_names.csv", row.names = F)
write.csv(good_secondary_names, "~/Downloads/good_secondary_names.csv", row.names = F)

