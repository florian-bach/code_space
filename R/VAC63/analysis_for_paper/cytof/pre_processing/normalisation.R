# library(Biobase)
# library(data.table)
  library(flowCore)
  library(CATALYST)
#library(premessa)


#    CONCATENATION   ####

#setwd("/media/flobuntu/Backups/cytof_imd/vac63c/unprocessed/09082019/09082019_unprocessed_imd/not_normalised_fcs")
#setwd("/media/flobuntu/Backups/cytof_imd/vac63c/unprocessed/13082019/13082019_unprocessed/IMD/not_normalised")
setwd("/media/flobuntu/Backups/cytof_imd/vac63c/unprocessed/15082019/IMD/not_normalised")




files <- list.files(pattern="*.fcs")

dir.create("concat")
system.time(concatenate_fcs_files(files, output.file = "./concat/non_normalised_concat.fcs"))





#    NORMALISATION    ####

#setwd("/media/flobuntu/Backups/cytof_imd/vac63c/unprocessed/09082019/09082019_unprocessed_imd/not_normalised_fcs/concat")
#setwd("/media/flobuntu/Backups/cytof_imd/vac63c/unprocessed/13082019/13082019_unprocessed/IMD/not_normalised/concat")
#setwd("/media/flobuntu/Backups/cytof_imd/vac63c/unprocessed/15082019/IMD/not_normalised/concat")



# loop didn't run well, becaause the RAM allocation sucked(?) and it failed after the first successful iteration
dirs <- c("/media/flobuntu/Backups/cytof_imd/vac63c/unprocessed/09082019/09082019_unprocessed_imd/not_normalised_fcs/concat",
          "/media/flobuntu/Backups/cytof_imd/vac63c/unprocessed/13082019/13082019_unprocessed/IMD/not_normalised/concat",
          "/media/flobuntu/Backups/cytof_imd/vac63c/unprocessed/15082019/IMD/not_normalised/concat",
          "/media/flobuntu/Backups/cytof_imd/vac63c/beads/")

 
  
setwd(dirs[4])  

non_normalised_concat <- read.FCS("comp_beads_06_cells_found.fcs")

normalised_ff <- normCytof(x=non_normalised_concat, y="dvs", k=80, plot=FALSE, verbose=T);
#dir.create("normalised2")

write.flowSet(normalised_ff, outdir = "normalised")



  
##    DELETE  110Cd CHANNEL  ####
  
  #setwd("/media/flobuntu/Backups/cytof_imd/vac63c/unprocessed/09082019/09082019_unprocessed_imd/not_normalised_fcs/concat")
  #setwd("/media/flobuntu/Backups/cytof_imd/vac63c/unprocessed/15082019/IMD/not_normalised/concat")
  setwd("/media/flobuntu/Backups/cytof_imd/vac63c/processed/normalised/")
  
  
  #make premessa table to rename/remove channels
  vac63c_premessa_table <- read_parameters(list.files(".", pattern="*.fcs"))
  vac63c_premessa_table$Remove <- ifelse(vac63c_premessa_table[,1]=="110Cd", TRUE, FALSE)
  
  #get rid of metal in channel name
  
  # prm <- ifelse(grepl("_", vac63c_premessa_table$non_normalised_concat.fcs),
  #               substr(vac63c_premessa_table$non_normalised_concat.fcs, 7, 12),
  #               vac63c_premessa_table$non_normalised_concat.fcs)
  # 
  # prm[2:3] <- c("Event_Length", "CD45")
  # 
  
  
  prm <-rownames(vac63c_premessa_table)

  vac63c_premessa_table$Parameter <- prm

  write.csv(vac63c_premessa_table, "vac63c_premessa_table.csv")
  
  
  
  # apply changes and write out
  system.time(rename_parameters_in_files(".", "renamed", vac63c_premessa_table))
  
  
  
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

