# library(Biobase)
# library(data.table)
library(flowCore)
library(CATALYST)



#    CONCATENATION   ####

#setwd("/media/flobuntu/Backups/cytof_imd/vac63c/unprocessed/09082019/09082019_unprocessed_imd/not_normalised_fcs")
#setwd("/media/flobuntu/Backups/cytof_imd/vac63c/unprocessed/13082019/13082019_unprocessed/IMD/not_normalised")
setwd("~/flo_r_biz/")




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

normalised_ff <- normCytof(x=non_normalised_concat, y="dvs", k=80, plot=FALSE, verbose=T, remove_beads=T);
dir.create("normalised2")

write.flowSet(normalised_ff, outdir = "normalised2")



  
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

  #write.csv(vac63c_premessa_table, "vac63c_premessa_table.csv")
  
  
  
  # apply changes and write out
  system.time(rename_parameters_in_files(".", "renamed", vac63c_premessa_table))
  

  
  
# CALCULATE COMPENSATION ####
  
  # get single-stained control samples
  # completely irrelevant whther the normalisation beads are still in the file btw
  ss_exp <- read.FCS("~/PhD/cytof/beads_and_tuning_reports/aug2019/normalised_renamed/vac63c_comp_beads.fcs")
  
    # specify mass channels stained for
  
  custom_isotope_list <- c(CATALYST::isotope_list,list("BCKG"=190))
  custom_isotope_list$Cd <- custom_isotope_list$Cd[-3]
  
  
  #bc_ms <-  c(89, 102, 104:106, 108, 110:116, 120, 127, 131, 133, 138, 140:156, 158:176, 190:195, 198, 208, 209)
  
  # vac63c_premessa_table <- read.csv("~/PhD/cytof/vac63c/vac63c_premessa_table.csv", header = T, stringsAsFactors = F, row.names = 1)
  # masses <- as.numeric(gsub("[a-zA-Z ]", "", rownames(vac63c_premessa_table)))
  # masses <- na.omit(masses)
  # 
  #including 190 here is curcial otherwise the compCytof function won't accept the custom isotope list
  bc_ms <-  c(89, 114, 115, 141:156, 158:176, 198, 209)
  # debarcode
  re <- assignPrelim(x=ss_exp, y=bc_ms, verbose=TRUE)
  re <- estCutoffs(x=re)
  re <- applyCutoffs(x=re)
  # compute spillover matrix
  spillMat <- computeSpillmat(x=re)
  # 
  write.csv(spillMat, "~/PhD/cytof/beads_and_tuning_reports/aug2019/aug2019_spillmat_final.csv")
  # 
# APPLY COMPENSATION TO FLOWSET####
# setwd("~/PhD/cytof/vac63c/normalised_renamed/debarcoded_not_comped/")
# vac63c <- read.flowSet(path = "~/PhD/cytof/vac63c/normalised_renamed/debarcoded_not_comped/", pattern = "*fcs")
# 
# test <- read.FCS("~/PhD/cytof/vac69a/reprocessed/relabeled/V02_DoD.fcs")
# comped_nnls <- compCytof(x=test, y=spillMat, method="nnls", isotope_list=custom_isotope_list, out_path="~/PhD/cytof/vac69a/reprocessed/reprocessed_comped/")
  
spillMat <- read.csv("~/PhD/cytof/beads_and_tuning_reports/aug2019/aug2019_spillmat_final.csv", header = T, row.names = 1)
test <- flowCore::read.FCS("~/PhD/cytof/vac63c/normalised_renamed/batch2_normalised_renamed_concat.fcs")
system.time(CATALYST::compCytof(x=test, y=spillMat[,-55], method="nnls", out_path="batch2_normalised_renamed_comped_concat.fcs"))
# 
# comp <- function(ff){compCytof(x=ff,
#                                y=spillMat,
#                                out_path="~/PhD/cytof/vac63c/normalised_renamed_comped/",
#                                method="nnls",
#                                isotope_list=custom_isotope_list)}
# 
#   compCytof(x=batch1,
#             isotope_list=custom_isotope_list,
#             y=spillMat,
#             out_path="~/PhD/cytof/vac63c/normalised_renamed_comped/",
#             method="nnls")  
#   

  # 
  # 
  # filez <- list.files(pattern = "*fcs")
  # 
  # try <- flowCore::read.flowSet(filez[1:3])
      
  
# debarcoding #####


# custom_isotope_list <- c(CATALYST::isotope_list,list(BCKG=190))
# custom_isotope_list$Cd
# # [1] 106 108 110 111 112 113 114 116
# custom_isotope_list$Cd <- c(106, 108, 111, 112, 113, 114, 116)
 

batch1 <- "/home/irina/flo_r_biz/comped/batch1_normalised_renamed_concat_comped.fcs"
batch2 <- "/home/irina/flo_r_biz/comped//batch2_normalised_renamed_concat_comped.fcs"
batch3 <- "/home/irina/flo_r_biz/comped/batch3_normalised_renamed_concat_comped.fcs"


batch1_names <- "/home/irina/flo_r_biz/batch1_debarcoder.csv"
batch2_names <- "/home/irina/flo_r_biz/batch2_debarcoder.csv"
batch3_names <- "/home/irina/flo_r_biz/batch3_debarcoder.csv"

batch1_names <- read.csv(batch1_names, stringsAsFactors = F, header=F)
batch2_names <- read.csv(batch2_names, stringsAsFactors = F, header=F)
batch3_names <- read.csv(batch3_names, stringsAsFactors = F, header=F)

data(sample_key)

# this step is important because the algorithm will look for every barcode present in the sample_key object, even unused ones
batch1_sample_key <- sample_key[rownames(sample_key)%in%batch1_names[,1],]
batch2_sample_key <- sample_key[rownames(sample_key)%in%batch2_names[,1],]
batch3_sample_key <- sample_key[rownames(sample_key)%in%batch3_names[,1],]


flz <- batch3
smpl_key <- batch3_sample_key
nmz <- batch3_names[,2]


normalised_concat <- read.FCS(flz)

#sample_key <- sample_key[1:15,]
#15-25 min
system.time(re0 <- assignPrelim(x=normalised_concat, y=smpl_key, verbose=TRUE))
re0

re <- estCutoffs(x=re0)
plotYields(x=re, which=0, plotly=FALSE)

re <- applyCutoffs(x = re)
re


system.time(outFCS(x=re,
                   y=normalised_concat,
                   out_nms=nmz,
                   out_path = "flo_r_biz/comped/debarcoded/",
                   verbose=T)
            )




# other stuff ####

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

