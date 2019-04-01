library(FlowSOM) 
library(readxl) 
library(flowCore) 
library(premessa)

setwd("/Users/s1249052/PhD/flow_data/vac69a/t_cells_only/better_gating/")

files_list <- list.files(path=".", pattern="*.fcs")

concatenate_fcs_files(files_list, "concat02.fcs")

flo_set <- read.flowSet(files_list, transformation = FALSE, truncate_max_range = FALSE)

fcs_raw <- list()
for (i in files_list){
  print(i)
  j <- read.FCS(i, transformation = FALSE, truncate_max_range = FALSE)
  assign(paste(substr(i, nchar(i)-10, nchar(i)-4)), j)
  fcs_raw[[paste(substr(i, nchar(i)-10, nchar(i)-4))]] <- j}

head(fcs_raw)

df_lol <- as.matrix(fcs_raw[[1]]@parameters@data[["desc"]])
channel_names <- df_lol[1:69]

md <- read_excel("metadata.xlsx")
sample_ids <- rep(md$volunteer, fsApply(fcs_raw, nrow))





no_of_dfs <- seq(1,30)

for (i in no_of_dfs){fsApply(fcs_raw[[i]], )
  
}

## arcsinh transformation and column subsetting
fcs <- fsApply(flo_set, function(x, cofactor = 5){
  colnames(x) <- channel_names
  expr <- exprs(x)
  expr <- asinh(expr[, colnames(x)] / cofactor)
  exprs(x) <- expr
  x
})

fcs

expr <- fsApply(fcs, exprs)
dim(expr)


fsom <- ReadInput(fcs, transform = FALSE, scale = FALSE)
set.seed(1234)
som <- BuildSOM(fsom, colsToUse = lineage_markers)
