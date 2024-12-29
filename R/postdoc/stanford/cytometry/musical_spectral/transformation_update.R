library(ggplot2)
library(flowCore)
library(tidyr)
library(dplyr)


hard_downsample <- function(fs, event_number){
  flowCore::fsApply(fs, function(ff){
    idx <- sample.int(nrow(ff), min(event_number, nrow(ff)))
    ff[idx,]
  })
}

# fcs_files <- list.files("~/Downloads/fcs-selected", pattern = ".fcs")
# full_paths <-  list.files("~/Downloads/fcs-selected", pattern = ".fcs", full.names = TRUE)
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

metadata_to_read <- kylie_metadata%>%
  dplyr::filter(live_name %in% fcs_files)%>%
  dplyr::filter(grepl("^lrs", id))

metadata_to_read$file_path <-  full_paths[match(metadata_to_read$live_name, fcs_files)]
# set.seed(1234)
# musical_flowset <- hard_downsample(musical_flowset, event_number = 10000)

# musical_flowset <- ncdfFlow::read.ncdfFlowSet(metadata_to_read$file_path)
big_musical_flowset <- flowCore::read.flowSet(files = metadata_to_read$file_path, truncate_max_range = FALSE)

set.seed(1234)
musical_flowset <- hard_downsample(big_musical_flowset, event_number = 50000)


# old
# co_factors = c(80000, 40000, 8000, 100000, 6000, 30000,
#                12000, 8000,8000, 500, 3000, 2000, 10000,
#                200, 8000, 8000, 8000, 5000, 5000, 5000,
#                8000, 8000, 4000, 8000, 6000, 100000, 
#                3000, 100000, 7000, 7000, 70000)

co_factors = data.frame("coef"=c(8000, 6000, 8000, 6000, 6000, 30000,
               12000, 8000,8000, 5000, 3000, 6000, 6000,
               200, 8000, 8000, 8000, 5000, 5000, 5000,
               8000, 8000, 4000, 8000, 12000, 6000, 
               3000, 100000, 7000, 7000, 6000),
               "Fluorophore"=colnames(musical_flowset[[1]])[8:38],
               "Marker"=musical_panel$Marker[match(colnames(musical_flowset)[8:38], musical_panel$Fluorophore)])

co_factors$coef[co_factors$Marker=="CCR7CD197"] <- 6000
co_factors$coef[co_factors$Marker=="CD127"] <- 12000
co_factors$coef[co_factors$Marker=="CD2"] <- 12000
co_factors$coef[co_factors$Marker=="IgM"] <- 6000
co_factors$coef[co_factors$Marker=="CD24"] <- 6000
co_factors$coef[co_factors$Marker=="CD56"] <- 3000
co_factors$coef[co_factors$Marker=="CD16"] <- 1000
co_factors$coef[co_factors$Marker=="TCR VD1"] <- 12000
trans_coefs <- write.csv("~/postdoc/stanford/cytometry/spectral/MUSICAL/trans_coefs.csv", row.names = FALSE)

trans_musical_flowset <- flowSpecs::arcTrans(musical_flowset,  transCoFacs = co_factors$coef, transNames = colnames(musical_flowset)[8:38])

equal_trans_musical_flowset <- flowSpecs::arcTrans(musical_flowset,  transCoFacs = 6000, transNames = colnames(musical_flowset)[8:38])

long_musical_flowset <- flowSpecs::flowSet2LongDf(musical_flowset)%>%
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


histogram_theme <- theme(legend.position = "bottom", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=7))


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





interaction_colors <- c(colorspace::sequential_hcl("Greens", n=5)[1:4], colorspace::sequential_hcl("Blues", n=5)[1:4], colorspace::sequential_hcl("Magenta", n=5)[1:4])

uni_trans <- ggplot(long_equal_trans_musical_flowset, aes(x=expression, color=interaction(experiment, id)))+
  geom_density()+
  facet_wrap(~marker+id, scales = "free")+
  ggtitle("cofactor = 6000")+
  theme_minimal()+
  scale_color_manual(values = interaction_colors)+
  scale_x_continuous(labels = scales::label_scientific())+
  histogram_theme

ggsave("~/postdoc/stanford/cytometry/spectral/MUSICAL/indie_uni_trans_histogram.png", uni_trans, width=24, height=40, bg="white", dpi=444)

multi_trans <- ggplot(long_trans_musical_flowset, aes(x=expression, color=interaction(experiment, id)))+
  geom_density()+
  facet_wrap(~marker+id, scales = "free")+
  ggtitle("bespoke cofactors")+
  theme_minimal()+
  scale_color_manual(values = interaction_colors)+
  scale_x_continuous(labels = scales::label_scientific())+
  histogram_theme

ggsave("~/postdoc/stanford/cytometry/spectral/MUSICAL/indie_multi_trans_histogram.png", multi_trans, width=24, height=40, bg="white", dpi=444)






i=match("C", musical_panel$Marker);hist(exprs(musical_flowset[[1]])[,i+8], main = "no transformation",
                                           xlab = musical_panel$Marker[i], breaks = 200)

i=match("CD3", musical_panel$Marker);hist(exprs(trans_musical_flowset[[1]])[,i+8], main = "multi transformation",
         xlab = musical_panel$Marker[i], breaks = 200)

i=match("CD56", musical_panel$Marker);hist(exprs(equal_trans_musical_flowset[[1]])[,i+8], main = "uni transformation",
         xlab = musical_panel$Marker[i], breaks = 200)
