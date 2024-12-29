library(ggplot2)
library(flowCore)
library(tidyr)
library(dplyr)
library(flowSpecs)

# read in data ####

# this is a function for downsampling flowsets to a specified number of events
hard_downsample <- function(fs, event_number){
  flowCore::fsApply(fs, function(ff){
    idx <- sample.int(nrow(ff), min(event_number, nrow(ff)))
    ff[idx,]
  })
}

fcs_files <- list.files("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/Aurora Spectral Flow/For Arefin/FCSFilesLiveSinglets", pattern = ".fcs")
full_paths <-  list.files("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/Aurora Spectral Flow/For Arefin/FCSFilesLiveSinglets", pattern = ".fcs", full.names = TRUE)

# read in panel and adapt for CATALYST to match the information stored in the FCS files
musical_panel <- readxl::read_excel("~/Downloads/MUSICALSpectralFlowPanel.xlsx")
musical_panel$class <- "type"
musical_panel$Fluorophore <- paste(musical_panel$Fluorophore, "A", sep="-")

#read in metadata & add a filename column to reflects the live singlet name 
kylie_metadata <- readxl::read_excel("~/Downloads/MUSICALSpectralFlowMetadata.xlsx")
colnames(kylie_metadata)[1] <- "file_name" 
kylie_metadata$live_name <- paste(kylie_metadata$experiment, kylie_metadata$file_name, sep="_")
kylie_metadata$live_name <- gsub(".fcs", "_Live Cells.fcs", kylie_metadata$live_name, fixed = T)

#subset metadata to only inlcude cells of interest, here, all US healthy control samples, LRS
metadata_to_read <- kylie_metadata%>%
  dplyr::filter(live_name %in% fcs_files)%>%
  dplyr::filter(grepl("^lrs", id))

#add a column that contains the path to the live singlet fcs file
metadata_to_read$file_path <-  full_paths[match(metadata_to_read$live_name, fcs_files)]

# read in flowset; truncate_max_range=FALSE is important for spectral flow
big_musical_flowset <- flowCore::read.flowSet(files = metadata_to_read$file_path, truncate_max_range = FALSE)

# downsample, setting seed so that it's reproducible
set.seed(1234)
musical_flowset <- hard_downsample(big_musical_flowset, event_number = 50000)


# transform data and inspect quality ####
#read in cofactors for asinh transformation
trans_coefs <- read.csv("~/postdoc/stanford/cytometry/spectral/MUSICAL/trans_coefs.csv", row.names = FALSE)

#transform data with bespoke cofactors;
# this version seems OK, but I think we should try coming up with experiment-specific transformations,
# certainly for a couple of markers like IgD and TCR Vd1 that looks like it might be a good idea.
trans_musical_flowset <- flowSpecs::arcTrans(musical_flowset,  transCoFacs = co_factors$coef, transNames = colnames(musical_flowset)[8:38])

# alternative transformation with same cofactor for every channel
equal_trans_musical_flowset <- flowSpecs::arcTrans(musical_flowset,  transCoFacs = 6000, transNames = colnames(musical_flowset)[8:38])

# turn into long format for plotting (i.e. each cell is one row), fix channel names, add marker column
long_musical_flowset <- flowSpecs::flowSet2LongDf(musical_flowset)%>%
  mutate(live_name=names)%>%
  #add metadata information
  left_join(., metadata_to_read, by = "live_name")%>%
  # make even longer, so that every cell is split into n rows, where n is the number of channels
  pivot_longer(cols=colnames(.)[8:38], names_to = "channel", values_to = "expression")%>%
  mutate(channel=gsub(".A", "-A", channel, fixed = TRUE),
         channel=gsub(".Cy", "-Cy", channel, fixed=TRUE),
         channel=gsub(".Fire", "-Fire", channel, fixed=TRUE),
         channel=gsub(".Dazzle", "-Dazzle", channel, fixed=TRUE),
         channel=gsub(".", " ", channel, fixed = TRUE),
         channel=gsub("5 5", "5.5", channel, fixed = TRUE))%>%
  # add marker information for every channel
  mutate(marker=musical_panel$Marker[match(channel, musical_panel$Fluorophore)])


long_trans_musical_flowset <- flowSpecs::flowSet2LongDf(trans_musical_flowset) %>%
  mutate(live_name=names)%>%
  left_join(., metadata_to_read, by = "live_name")%>%
  pivot_longer(cols=colnames(.)[8:38], names_to = "channel", values_to = "expression")%>%
  mutate(channel=gsub(".A", "-A", channel, fixed = TRUE),
         channel=gsub(".Cy", "-Cy", channel, fixed=TRUE),
         channel=gsub(".Fire", "-Fire", channel, fixed=TRUE),
         channel=gsub(".Dazzle", "-Dazzle", channel, fixed=TRUE),
         channel=gsub(".", " ", channel, fixed = TRUE),
         channel=gsub("5 5", "5.5", channel, fixed = TRUE))%>%
  mutate(marker=musical_panel$Marker[match(channel, musical_panel$Fluorophore)])


long_equal_trans_musical_flowset <- flowSpecs::flowSet2LongDf(equal_trans_musical_flowset)%>%
  mutate(live_name=names)%>%
  left_join(., metadata_to_read, by = "live_name")%>%
  pivot_longer(cols=colnames(.)[8:38], names_to = "channel", values_to = "expression")%>%
  mutate(channel=gsub(".A", "-A", channel, fixed = TRUE),
         channel=gsub(".Cy", "-Cy", channel, fixed=TRUE),
         channel=gsub(".Fire", "-Fire", channel, fixed=TRUE),
         channel=gsub(".Dazzle", "-Dazzle", channel, fixed=TRUE),
         channel=gsub(".", " ", channel, fixed = TRUE),
         channel=gsub("5 5", "5.5", channel, fixed = TRUE))%>%
  mutate(marker=musical_panel$Marker[match(channel, musical_panel$Fluorophore)])


# set common plot modification for histograms
histogram_theme <- theme(legend.position = "bottom", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=7))


# make a bunch of hisograms for each channel
no_trans <- ggplot(long_musical_flowset, aes(x=expression, color=experiment))+
  geom_density()+
  facet_wrap(~marker, scales = "free")+
  ggtitle("no transformation")+
  theme_minimal()+
  scale_x_continuous(labels = scales::label_scientific())+
  histogram_theme

ggsave("~/postdoc/stanford/cytometry/spectral/MUSICAL/no_trans_histogram.png", no_trans, width=12, height=12, bg="white", dpi=444)

uni_trans <- ggplot(long_equal_trans_musical_flowset, aes(x=expression, color=experiment))+
  geom_density()+
  facet_wrap(~marker, scales = "free")+
  ggtitle("cofactor = 6000")+
  theme_minimal()+
  scale_x_continuous(labels = scales::label_scientific(digits = 0, ))+
  histogram_theme

ggsave("~/postdoc/stanford/cytometry/spectral/MUSICAL/average_uni_trans_histogram.png", uni_trans, width=12, height=12, bg="white", dpi=444)

multi_trans <- ggplot(long_trans_musical_flowset, aes(x=expression, color=experiment))+
  geom_density()+
  facet_wrap(~marker, scales = "free")+
  ggtitle("bespoke cofactors")+
  theme_minimal()+
  scale_x_continuous(labels = scales::label_scientific())+
  histogram_theme

ggsave("~/postdoc/stanford/cytometry/spectral/MUSICAL/average_multi_trans_histogram.png", multi_trans, width=12, height=12, bg="white", dpi=444)


# load data into CATALYST ####
sce <- prepData(trans_musical_flowset,
                FACS = T,ass
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

# add expression assay
assay(sce, "exprs") <- assay(sce, "counts")

# define markers for clustering, everything except CD45 & LiveDead (because they were used in pregating);
# and no IgD because those stains look a little problematic for now
cluster_markers <- musical_panel$Marker[-match(c("IgD", "TCR VD1", "LiveDead", "CD45"), musical_panel$Marker)]


# flowSOM clustering ####
# run Flowsome with (almost) default values, including seed for reproducibility
sce <- CATALYST::cluster(sce,
                         features = cluster_markers,
                         xdim = 12, ydim = 12, maxK = 50, seed = 1234)

# plotting markers include IgD and TCR Vd1 for inspection
plotting_markers <- c("CD19", "CD20", "IgG", "IgM", "IgD", "CD21", "CD27", "HLA DR",  "CD11c", "CD1d", "CD123", "CD33", "CD14", "CD16", "CD56", "CD7", "CD85j", "CD24", "CD2", "TCR VD1", "TCR VD2", "CD57", "CD127", "CD3", "CD8", "CCR7CD197", "CD45RA")

# make a heatmap
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

# generally it's a good idea to overcluster (say to 50) and then merge clusters of similar phenotypes;
# this can be tricky and is often an iterative process
# NB this is just an example file path, I haven't done any merging yet
rough_merge <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/rough_activated_50_merge.csv", header = TRUE)
sce <- mergeClusters(sce, k="meta50", rough_merge, "rough_merge", overwrite = T)


# dimensionality reduction ####
sce <- runDR(sce,
             seed=1234,
             dr="UMAP",
             features=cluster_markers)

subset_umap <- plotDR(sce,  color_by = "rough_merge", scale = TRUE, k_pal = rough_merge_palette)
