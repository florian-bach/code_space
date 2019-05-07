

##### USER INPUT #####

# Install packages if required

# Load packages
library('flowCore')

# Set working directory
getwd()
setwd("/Users/s1249052/PhD/cytof/vac69a/T_cells_only/") # set to your desired working folder (directory)
PrimaryDirectory <- getwd()

# Find file names of .fcs files in the current working directory
FileNames <- list.files(path=PrimaryDirectory, pattern = ".fcs")
FileNames


# Chose which .fcs file to read into 'data' -- rename 'sample_data.fcs' to whatever file you want to read
data <- exprs(read.FCS("14c10_sample1_01_0_0_T_cells.fcs", transformation = FALSE))
data

# Give a name for your output .csv file (don't foret to add '.csv' at the end of the name)
csvfilename <- "vac69a_07_t+6.csv"


##### END USER INPUT #####

# Create an 'output' folder
setwd(PrimaryDirectory)
dir.create("Output", showWarnings = FALSE)
setwd("Output")

# Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
write.csv(data, csvfilename)

setwd(PrimaryDirectory)


column_names <- colnames(read.csv("vac69a_07_t+6.csv", header=T))

      
column_names <- c("Cd114Di", "In115Di","Pr141Di", "Nd142Di", "Nd143Di", "Nd144Di", "Nd145Di", "Nd146Di", "Sm147Di", "Nd148Di", "Sm149Di", "Nd150Di", "Eu151Di", "Eu153Di", "Sm154Di", "Gd155Di", "Gd156Di", "Gd158Di",
"Tb159Di", "Gd160Di", "Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di", "Ho165Di", "Er166Di", "Er167Di", "Er168Di",
"Tm169Di", "Yb171Di", "Yb172Di", "Yb173Di", "Yb174Di", "Lu175Di", "Yb176Di", "Pt198Di", "Bi209Di")

column_targets <- c("114Cd_CD14 <Cd114Di>", "115In_CD57 <In115Di>", "141Pr_HLA-DR <Pr141Di>", "142Nd_BCL-2 <Nd142Di>", "143Nd_CD45RA <Nd143Di>", "144Nd_GZB <Nd144Di>",
"145Nd_CD4 <Nd145Di>", "146Nd_Vd2 <Nd146Di>", "147Sm_CD20 <Sm147Di>", "148Nd_ICOS <Nd148Di>", "149Sm_CXCR5 <Sm149Di>", "150Nd_CD95 <Nd150Di>", "151Eu_CD103 <Eu151Di>",
"153Eu_Va7.2 <Eu153Di>",	"154Sm_TIM-3 <Sm154Di>",	"155Gd_PD1 <Gd155Di>",	"156Gd_CD161 <Gd156Di>",	"158Gd_CD27 <Gd158Di>",
"159Tb_FoxP3 <Tb159Di>", "160Gd_CTLA4 <Gd160Di>", "161Dy_Tbet <Dy161Di>", "162Dy_IntegrinB7 <Dy162Di>", "163Dy_CD28 <Dy163Di>", "164Dy_Ki-67 <Dy164Di>", "165Ho_CD45RO <Ho165Di>",
"166Er_CD56 <Er166Di>", "167Er_CCR7 <Er167Di>", "168Er_CD127 <Er168Di>", "169Tm_CD38 <Tm169Di>", "171Yb_CD49d <Yb171Di>",
"172Yb_CD25 <Yb172Di>", "173Yb_CD39 <Yb173Di>", "174Yb_CLA <Yb174Di>", "175Lu_Perforin <Lu175Di>", "176Yb_CX3CR1 <Yb176Di>",
"198Pt_CD8 <Pt198Di>", "209Bi_CD16 <Bi209Di>")

library(stringr)
str_sub(column_targets, start=1, end=-10)

'114Cd_CD14', '115In_CD57', '141Pr_HLA-DR', '142Nd_BCL-2', '143Nd_CD45RA', '144Nd_GZB', '145Nd_CD4', '146Nd_Vd2'
'147Sm_CD20', '148Nd_ICOS', '149Sm_CXCR5', '150Nd_CD95', '151Eu_CD103', '153Eu_Va7.2', '154Sm_TIM-3', '155Gd_PD1'
'156Gd_CD161', '158Gd_CD27', '159Tb_FoxP3', '160Gd_CTLA4', '161Dy_Tbet', '162Dy_IntegrinB7', '163Dy_CD28', '164Dy_Ki-67'
'165Ho_CD45RO', '166Er_CD56', '167Er_CCR7', '168Er_CD127', '169Tm_CD38', '171Yb_CD49d', '172Yb_CD25', '173Yb_CD39'
'174Yb_CLA', '175Lu_Perforin', '176Yb_CX3CR1', '198Pt_CD8', '209Bi_CD16'



######

######

# CSV to FCS
# Coverting .csv file data into an .fcs file
# Thomas Ashhurst
# 2017-09-13
# github.com/sydneycytometry
# .fcs file reading and writing adapted from https://gist.github.com/yannabraham/c1f9de9b23fb94105ca5


##### USER INPUT #####


library('Biobase')
library('data.table')
setwd("/Users/s1249052/PhD/flow data/vac69a/t cells only/experiment_210618_files/csv_by_person/all/big_one")
## Use to list the .csv files in the working directory -- important, the only CSV files in the directory should be the one desired for analysis. If more than one are found, only the first file will be used
PrimaryDirectory <- getwd()
FileNames <- list.files(path=PrimaryDirectory, pattern = ".csv")     # see a list of CSV files
as.matrix(FileNames) # See file names in a list

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



##### END USER INPUT #####

x <- Sys.time()
x <- gsub(":", "-", x)
x <- gsub(" ", "_", x)

newdir <- paste0("Output_CSV-to-FCS", "_", x)

setwd(PrimaryDirectory)
dir.create(paste0(newdir), showWarnings = FALSE)
setwd(newdir)


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


########     LOOP TO CONVERT MULTIPLE FCS TO CSV

FileNames <- list.files(path=PrimaryDirectory, pattern = ".fcs")     # see a list of CSV files
as.matrix(FileNames)

for (File in FileNames) { # Loop to read files into the list
  data <- exprs(read.FCS(File, transformation = FALSE))
  write.csv(data, gsub(".fcs", ".csv", File))
}


