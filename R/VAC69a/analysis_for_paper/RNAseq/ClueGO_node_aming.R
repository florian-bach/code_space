library(magrittr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readxl)


dod_data <- read_xls("~/PhD/RNAseq/vac69a/cytoscape/VIVAX_DOD_ALL/VIVAX_DOD_ALL_RESULTS_TABLE.xls")
t6_data <-  read_xls("~/PhD/RNAseq/vac69a/cytoscape/VIVAX_T6_ALL/VIVAX_T6_ALL_RESULTS_TABLE")

base_dod_t6_data <- read_xls("~/PhD/RNAseq/vac69a/cytoscape/VIVAX_BASE_DOD_T6/VIVAX_BASE_DOD_T6_Results_table.xls")

vivax_falciparum_dod <- read_xls("~/PhD/RNAseq/vac69a/cytoscape/VIVAX_FALCIPARUM_DOD/vivax_falciparum_dod_results_table.xls")
# colnames(vivax_falciparum_dod) <- gsub(".", " ", colnames(vivax_falciparum_dod), fixed=T)
# colnames(vivax_falciparum_dod) <- gsub("  ", " ", colnames(vivax_falciparum_dod), fixed=T)
# colnames(vivax_falciparum_dod) <- gsub("X", "", colnames(vivax_falciparum_dod), fixed=T)
# colnames(vivax_falciparum_dod) <- gsub(" Associated Genes","% Associated Genes",colnames(vivax_falciparum_dod), fixed=T)


vivax_falciparum_t6 <-read_xls("~/PhD/RNAseq/vac69a/cytoscape/VIVAX_FALCIPARUM_T6/VIVAX_FALCIPARUM_T6_Results_Table")

data <- vivax_falciparum_dod
# data <- t6_data


#filter GOTerms to be level 5 only
data <- subset(data, grepl("5", data$GOLevels))


names_dod <- data %>%
  filter(`% Associated Genes`>10) %>%
  arrange(`Term PValue`) %>%
  group_by(GOGroups) %>%
  top_n(30, GOGroups) %>%
  #group_by(GOGroups) %>%
  #top_n(-1, GOTerm) %>%
  select(GOTerm, `Term PValue`, GOGroups, `Group PValue`) 

names_dod <- names_dod[!duplicated(names_dod$GOGroups),]
names_dod <- names_dod[!duplicated(names_dod$GOTerm),]

View(subset(names_dod, names_dod$GOGroups %in% dod_level_siz$GOGroups[grep("*leukocyte chemotaxis*", dod_level_siz$GOTerm)]))




# names_dod <- dod_level_siz %>%
#   filter(`% Associated Genes`>10) %>%
#   group_by(GOGroups) %>%
#   top_n(-3, `Term PValue`) %>%
#   ungroup() %>%
#   #top_n(-60, `Group PValue`) %>%
#   arrange(`Group PValue`, `Term PValue`) %>%
#   select(GOTerm, `Term PValue`, GOGroups, `Group PValue`) 

names_dod <- names_dod[!duplicated(names_dod$GOTerm),]

names_dod <- names_dod %>%
  group_by(GOGroups) %>%
  top_n(-1, `Term PValue`) 




# write.csv(names_dod, "~/PhD/RNAseq/vac69a/cytoscape/VIVAX_FALCIPARUM_DOD/vivax_falciparum_dod_flo_group_names.csv")
# 
# write.table(names_dod$GOTerm, "~/PhD/RNAseq/vac69a/cytoscape/VIVAX_DOD_ALL/Top20_GOTerms_Vivax_DoD.txt", sep="\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)



names_t6 <- dod_level_siz %>%
  group_by(GOGroups) %>%
  filter(`% Associated Genes`>20) %>%
  #'group_by(gg1) %>%
  top_n(-3, `Term PValue`) %>%
  ungroup() %>%
  top_n(-60, `Group PValue`) %>%
  arrange(`Group PValue`, `Term PValue`) %>%
  select(GOTerm, `Term PValue`, GOGroups, `Group PValue`) 

names_t6 <- names_t6[!duplicated(names_t6$GOTerm),]

names_t6 <- names_t6 %>%
  group_by(GOGroups) %>%
  top_n(-1, `Term PValue`)

write.table(names_t6$GOTerm, "~/PhD/RNAseq/vac69a/cytoscape/VIVAX_T6_ALL/Top20_GOTerms_Vivax_T6.txt", sep="\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# [1] "regulation of programmed cell death"                     "pattern recognition receptor signaling pathway"         
# [3] "regulation of type I interferon production"              "response to interferon-alpha"                           
# [5] "negative regulation of cell cycle G2/M phase transition" "response to interferon-beta"                            
# [7] "response to lipopolysaccharide"                          "macrophage activation"                                  
# [9] "rRNA binding"                                            "protein kinase binding"                                 
# [11] "small GTPase mediated signal transduction"               "purine nucleotide binding"    



# garbage code (works tho) that splits GOGroup string into multiples and creates a column for each ####
#get rid of spaces and annoying brackets
# dod_level_siz$GOGroups <- gsub("[", "", dod_level_siz$GOGroups, fixed=T)
# dod_level_siz$GOGroups <- gsub("]", "", dod_level_siz$GOGroups, fixed=T)
# dod_level_siz$GOGroups <- gsub(" ", "", dod_level_siz$GOGroups, fixed=T)
# 
# #split the strings of GOGroups, the output is a nested list
# try <- strsplit(dod_level_siz$GOGroups, ",", fixed=T)

# use a for loop to make a df that contains 3 columns (max number of GOgroups) and populate with the different GOGroups
# associated with the GOterms

# try_df <- data.frame()
# 
# for (i in 1:nrow(dod_level_siz)){
# 
#   if(length(unlist(try[i]))==1)
#     new_row=data.frame("gg1"=unlist(try[i])[1], "gg2"=NA, "gg3"=NA)
# 
#   else if(length(unlist(try[i]))==2)
#     new_row=data.frame("gg1"=unlist(try[i])[1], "gg2"=unlist(try[i])[2], "gg3"=NA)
# 
#   else if(length(unlist(try[i]))==3)
#     new_row=data.frame("gg1"=unlist(try[i])[1], "gg2"=unlist(try[i])[2], "gg3"=unlist(try[i])[3])
# 
#   try_df <- rbind(try_df, new_row)
# 
# }
# 
# dod_level_siz <- cbind(dod_level_siz, try_df)
