library(flowCore)
library(CATALYST)

## code to prepare `DATASET` dataset goes here
setwd("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")
working_directory <- getwd()

# read in the metadata file. they MUST have the columns "file_name" and "sample_id" the rest is optional and depends on your experiment. in our case we'll add "timepoint", "batch" and "volunteer".
# then we order the table- I don't think this is necessary, but for some reason I wrote it that way at the time soooo
md <- read.csv("meta_data.csv", header=TRUE, stringsAsFactors = FALSE)
md$timepoint <- factor(md$timepoint, levels = c("Baseline", "C8", "C10", "C12", "Diagnosis", "T6"))
md <- md[with(md, order(md$volunteer, md$timepoint)),]

#read in panel; must have colnames of "fcs_colname" "antigen" "marker_class"
panel <- read.csv("VAC69_PANEL.CSV", header = TRUE, stringsAsFactors = FALSE)

### read in fcs files using flowCore; the way i wrote this means that you can filter the md table (if you only want to work on one sample type, volunteer or whatever, and thus only read in the files you need)
vac69a <- flowCore::read.flowSet(md$file_name)

hard_downsample <- function(fs, event_number){
  flowCore::fsApply(fs, function(ff){
    idx <- sample.int(nrow(ff), min(event_number, nrow(ff)))
    ff[idx,]
  })
}

smaller_vac69a <- hard_downsample(vac69a, 5000)

vac69a_single_cell_experiment <- CATALYST::prepData(smaller_vac69a,
                          panel=panel,
                          md=md,
                          md_cols = list(file = "file_name", id = "sample_id", factors = c("timepoint", "batch", "volunteer")),
                          panel_cols = list(channel = "fcs_colname", antigen = "antigen", class = "marker_class")
)

setwd("~/code_space/my_packages/wiebkeTOF/")
usethis::use_data(vac69a_single_cell_experiment, overwrite = TRUE)
