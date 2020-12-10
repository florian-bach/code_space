pkgname <- "vac69a.cytof"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "vac69a.cytof-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('vac69a.cytof')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("flo_umap")
### * flo_umap

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: flo_umap
### Title: This function makes pretty UMAP plots. It takes a dataframe as
###   input and uses the UMAP1 and UMAP2 columns as X and Y coordinates and
###   any other named column as the coloring variable. It can facet by one
###   variable.
### Aliases: flo_umap

### ** Examples

sce <- read_small("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/", 500)
df <- prep_sce_for_ggplot(sce)
cd25_plot <- flo_umap(df, "CD25", "timepoint")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("flo_umap", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("prep_sce_for_ggplot")
### * prep_sce_for_ggplot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: prep_sce_for_ggplot
### Title: This function makes data.frames when passed
###   SingleCellExperiments. This allows making ggplot objects from scratch
###   using sce data.
### Aliases: prep_sce_for_ggplot

### ** Examples

big_table <- prep_sce_for_ggplot(daf)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("prep_sce_for_ggplot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read_full")
### * read_full

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read_full
### Title: This function simplifies and streamlines reading the full vac69a
###   cytof dataset into R It relies on the packages flowCore, CATALYST and
###   deplyr and their dependcies (e.g. SingleCellExperiment)
### Aliases: read_full

### ** Examples


# you can call the resulting variable anything other than sce (sce is the internal name) #
daf <- read_full("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read_full", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read_small")
### * read_small

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read_small
### Title: This function simplifies and streamlines reading a subsampled
###   version of the vac69a cytof dataset into R It relies on the packages
###   flowCore, CATALYST and deplyr and their dependcies (e.g.
###   SingleCellExperiment)
### Aliases: read_small

### ** Examples


# you can call the resulting variable anything other than sce (sce is the internal name) #
daf <- read_small("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/", proportional=TRUE, event_number=3000)






base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read_small", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
