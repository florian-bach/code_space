

##### USER INPUT #####

# Install packages if required

# Load packages
library('flowCore')

# Set working directory
getwd()
setwd("/Users/s1249052/Downloads") # set to your desired working folder (directory)
PrimaryDirectory <- getwd()

# Find file names of .fcs files in the current working directory
FileNames <- list.files(path=PrimaryDirectory, pattern = ".fcs")
FileNames

# Chose which .fcs file to read into 'data' -- rename 'sample_data.fcs' to whatever file you want to read
data <- exprs(read.FCS("818-1_cd4_viSNE.fcs", transformation = FALSE))
data

# Give a name for your output .csv file (don't foret to add '.csv' at the end of the name)
csvfilename <- "818cd4onlybaseline.csv"


##### END USER INPUT #####

# Create an 'output' folder
setwd(PrimaryDirectory)
dir.create("Output", showWarnings = FALSE)
setwd("Output")

# Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
write.csv(data, csvfilename)

setwd(PrimaryDirectory)
