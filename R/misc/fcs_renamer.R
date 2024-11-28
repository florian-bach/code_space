library(premessa)

# reorder / rename the channel names ####
trial_file <- "~/Library/CloudStorage/Box-Box/Jagannathan_Lab_Folder/PROJECTS/Tr1/MUSICAL_Tr1_Flow/All_Data/20241115_malaria_tr1_MUS5/Plate 1/iRBC_G1_G01.fcs"
fcs_dir <- "~/Library/CloudStorage/Box-Box/Jagannathan_Lab_Folder/PROJECTS/Tr1/MUSICAL_Tr1_Flow/All_Data/20241115_malaria_tr1_MUS5/Plate 1"
# the true positive is: CXCR3 on BV421; CXCR6 on FITC
# the wrong label is CXCR3 BB515; CXCR6 on BV421
file_list <- list.files(fcs_dir, pattern = ".fcs", full.names = T)
# parameters <- read_parameters(file_list[18:length(file_list)])
parameters <- read_parameters(file_list[18:length(file_list)])

new_row_names <- rownames(parameters)
new_row_names[7] <- "BB515-A"
new_row_names[12] <- "BV421-A" 

# channels <- paste(rownames(parameters), parameters[,1])
channels <- new_row_names
# new_row_names[7] <- "BB515-A"
# new_row_names[12] <- "BV421-A" 

parameters$Remove <- FALSE
parameters$Parameter <- channels
parameters$Parameter <- ifelse(parameters$Parameter=="", rownames(parameters), parameters$Parameter)
# parameters$Parameter <- ifelse(is.na(parameters$Parameter), parameters[,ncol(parameters)], parameters$Parameter)

# apply changes and write out
rename_parameters_in_files(fcs_dir, "renamed2", parameters)


# didn't work as well####


# rename descriptions
trial_file <- "~/Library/CloudStorage/Box-Box/Jagannathan_Lab_Folder/PROJECTS/Tr1/MUSICAL_Tr1_Flow/All_Data/20241115_malaria_tr1_MUS5/Plate 1/iRBC_G1_G01.fcs"
fcs_dir <- "~/Library/CloudStorage/Box-Box/Jagannathan_Lab_Folder/PROJECTS/Tr1/MUSICAL_Tr1_Flow/All_Data/20241115_malaria_tr1_MUS5/Plate 1/"
# the true positive is: CXCR3 on BV421; CXCR6 on FITC
# the wrong label is CXCR3 BB515; CXCR6 on BV421

bad_parameters <- read_parameters(trial_file)

bad_parameters$Remove <- FALSE
bad_parameters$Channel <- rownames(bad_parameters)
bad_parameters$Parameter <-bad_parameters$iRBC_G1_G01.fcs
bad_parameters$Parameter <- ifelse(bad_parameters$Parameter=="", bad_parameters$Channel, bad_parameters$Parameter)

bad_parameters$Parameter[grep("BV421-A", rownames(bad_parameters))] <- "CXCR3"
bad_parameters$Parameter[grep("BB515-A", rownames(bad_parameters))] <-"CXCR6"

good_paramers <- bad_parameters
good_paramers$Channel <- NULL
good_paramers$Parameter <- paste(rownames(good_paramers), good_paramers$Parameter)

good_paramers$iRBC_G1_G01.fcs <- good_paramers$Parameter
# apply changes and write out
rename_parameters_in_files(fcs_dir, "renamed", good_paramers)








# try 2 ####
trial_file <- "~/Library/CloudStorage/Box-Box/Jagannathan_Lab_Folder/PROJECTS/Tr1/MUSICAL_Tr1_Flow/All_Data/20241115_malaria_tr1_MUS5/Plate 1/iRBC_G1_G01.fcs"
fcs_dir <- "~/Library/CloudStorage/Box-Box/Jagannathan_Lab_Folder/PROJECTS/Tr1/MUSICAL_Tr1_Flow/All_Data/20241115_malaria_tr1_MUS5/Plate 1"
# the true positive is: CXCR3 on BV421; CXCR6 on FITC
# the wrong label is CXCR3 BB515; CXCR6 on BV421
file_list <- list.files(fcs_dir, pattern = ".fcs", full.names = T)
parameters <- read_parameters(file_list[18])
# parameters <- read_parameters(file_list[18])

new_row_names <- rownames(parameters)
new_row_names[7] <- "BB515-A"
new_row_names[12] <- "BV421-A" 

# channels <- paste(rownames(parameters), parameters[,1])


parameters$Remove <- FALSE

parameters$Parameter <- new_row_names

rownames(parameters) <- new_row_names
# rownames(parameters) <- rownames(parameters)[c(1:6, 12, 8:11, 7, 13:24)]
# parameters$Parameter <- ifelse(parameters$Parameter=="", rownames(parameters), parameters$Parameter)
# parameters$Parameter <- ifelse(is.na(parameters$Parameter), parameters[,ncol(parameters)], parameters$Parameter)

# apply changes and write out
rename_parameters_in_files(fcs_dir, "renamed2", parameters)

