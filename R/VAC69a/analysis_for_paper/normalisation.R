library(Biobase)
library(data.table)
library(flowCore)
library(CATALYST)
library(premessa)


# concatenate
setwd("/media/flobuntu/data/cytof_imd/12022019/IMD/non_normalised/")
files <- list.files(pattern="*.fcs")

concatenate_fcs_files(files, output.file = "/media/flobuntu/data/cytof_imd/12022019/IMD/non_normalised/concat.fcs")
non_normalised_concat <- read.FCS("/media/flobuntu/data/cytof_imd/12022019/IMD/non_normalised/concat.fcs")

#normalise and write out
normalised_ff <- normCytof(x=non_normalised_concat, y="dvs", k=80, plot=FALSE)
write.flowSet(normalised_ff[1], "/media/flobuntu/data/cytof_imd/12022019/IMD/normalised/")


## delete 110cd channel ####
normalised_concat <- read.FCS("/media/flobuntu/data/cytof_imd/12022019/IMD/normalised/normalised_concat.fcs")
data <- exprs(normalised_concat)
data.table::fwrite(data, "normalised_concat.csv")

csv <- data.table::fread("normalised_concat.csv")
csv <- dplyr::select(csv, -Cd110Di, -V1)

data.table::fwrite(csv, "/media/flobuntu/data/cytof_imd/12022019/IMD/normalised/normalised_cd100_removed_concat.csv")

## convert to fcs ####


## Use to list the .csv files in the working directory -- important, the only CSV files in the directory should be the one desired for analysis. If more than one are found, only the first file will be used
PrimaryDirectory <- setwd("/media/flobuntu/data/cytof_imd/12022019/IMD/normalised")
FileNames <- "normalised_cd100_removed_concat.csv"

#FileNames <- list.files(path=PrimaryDirectory, pattern = ".csv")     # see a list of CSV files
#as.matrix(FileNames) # See file names in a list

## Read data from Files into list of data frames
DataList=list() # Creates and empty list to start 

for (File in FileNames) { # Loop to read files into the list
  tempdata <- fread(File, check.names = FALSE)
  File <- gsub(".csv", "", File)
  DataList[[File]] <- tempdata
}

rm(tempdata)
AllSampleNames <- names(DataList)

## Chech data quality
head(DataList)

##### END USER INPUT ##

x <- Sys.time()
x <- gsub(":", "-", x)
x <- gsub(" ", "_", x)

# newdir <- paste0("Output_CSV-to-FCS", "_", x)

setwd(PrimaryDirectory)
# dir.create(paste0(newdir), showWarnings = FALSE)
# setwd(newdir)


for(i in c(1:length(AllSampleNames))){
  data_subset <- DataList[i]
  data_subset <- rbindlist(as.list(data_subset))
  dim(data_subset)
  a <- names(DataList)[i]
  
  metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
  
  ## Create FCS file metadata - ranges, min, and max settings
  #metadata$range <- apply(apply(data_subset,2,range),2,diff)
  metadata$minRange <- apply(data_subset,2,min)
  metadata$maxRange <- apply(data_subset,2,max)
  
  data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
  head(data_subset.ff)
  write.FCS(data_subset.ff, paste0(a, ".fcs"))
}

# debarcoding #####
# custom_isotope_list <- c(CATALYST::isotope_list,list(BCKG=190))
# custom_isotope_list$Cd
# # [1] 106 108 110 111 112 113 114 116
# custom_isotope_list$Cd <- c(106, 108, 111, 112, 113, 114, 116)


normalised_concat <- read.FCS("/media/flobuntu/data/cytof_imd/12022019/IMD/normalised/normalised_cd100_removed_concat.fcs")

data(sample_key)
sample_key <- sample_key[1:15,]
#15-25 min
re0 <- assignPrelim(x=normalised_concat, y=sample_key, verbose=FALSE)
re0

re <- estCutoffs(x=re0)
plotYields(x=re, which=0, plotly=FALSE)

re <- applyCutoffs(x = re)


# plotMahal(x = re, which = "B3")
# 

names <- read.csv("/home/flobuntu/Documents/barcode_key_20190212.csv", stringsAsFactors = F, header=F)

outFCS(x=re, y=normalised_concat, out_nms=names[,2], out_path = "/media/flobuntu/data/cytof_imd/12022019/IMD/debarcoded", verbose=T)


setwd("/media/flobuntu/Backups/IMD1402/IMD/debarcoded")
debarcoded <- read.flowSet(files = list.files(path = ".", pattern = "*.fcs")[2:16])
uno <- read.FCS("V02_DoD.fcs")

ggcyto(debarcoded, aes(x=Nd142Di, y=Tm169Di))+
  geom_hex(bins=100)+
  scale_x_flowJo_biexp(equal.space = T)+
  scale_y_flowJo_biexp(equal.space = T)


#fix description in fcs file (this was lost during the fcs to csv to fcs transfromation)
setwd("~/PhD/cytof/vac69a/reprocessed/")
debarcoded_files <- list.files("~/PhD/cytof/vac69a/reprocessed/", pattern = "*fcs")

batch1 <-read.flowSet(debarcoded_files[16:30])

batch2 <- read.flowSet(debarcoded_files[1:15])

params1 <- parameters(batch1[[1]])
pData(params1)

desc <- description(batch1[[1]])

c('114Cd_CD14', '115In_CD57', '141Pr_HLA-DR', '142Nd_BCL2', '143Nd_CD45RA', '144Nd_GZB', '145Nd_CD4', '146Nd_Vd2'
'147Sm_CD20', '148Nd_ICOS', '149Sm_CXCR5', '150Nd_CD95', '151Eu_CD103', '153Eu_Va7.2', '154Sm_TIM-3', '155Gd_PD1'
'156Gd_CD161', '158Gd_CD27', '159Tb_FoxP3', '160Gd_CTLA4', '161Dy_Tbet', '162Dy_IntegrinB7', '163Dy_CD28', '164Dy_Ki-67'
'165Ho_CD45RO', '166Er_CD56', '167Er_CCR7', '168Er_CD127', '169Tm_CD38', '171Yb_CD49d', '172Yb_CD25', '173Yb_CD39'
'174Yb_CLA', '175Lu_Perforin', '176Yb_CX3CR1', '198Pt_CD8', '209Bi_CD16')

panel <- read.csv("/home/flobuntu/PhD/cytof/vac69a/T_cells_only/fcs/VAC69_PANEL.CSV")

#fsApply(batch1, function(x){pData(parameters(x))$desc <- panel$antigen[1:64]})

ggcyto(batch1, aes(x=Pt194Di, y=Ce140Di))+
  geom_hex(bins=90000)+
  # scale_x_flowCore_fasinh(b=5)+
  # scale_y_flowCore_fasinh(b=5)+
  scale_x_flowJo_biexp(maxValue = 1000)+
  scale_y_flowJo_biexp(maxValue = 1000)+
  #theme_minimal()+
  scale_fill_viridis(option="B")

# files <- list.files(pattern="*.fcs")
# 
# raw_data <- read.flowSet(files)
# 
# ff <- concatFCS(raw_data, fn="raw_data_concat.fcs")
# 
# write.FCS(ff, "raw_data_concat.fcs", what = "numeric")
# 
# non_normalised_concat <- read.FCS("raw_data_concat.fcs")







# concatenate_fcs_files(files, output.file = "./concat/raw_data_concat.fcs")
# non_normalised_concat <- read.FCS("./concat/raw_data_concat.fcs")
# 
# normalised_ff <- normCytof(x=non_normalised_concat, y="dvs", k=80, plot=FALSE)
# 
# write.flowSet(normalised_ff[1], "./concat/normalised_concat.fcs")
# 

