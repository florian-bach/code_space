library(immunarch)

sce <- readRDS("/Volumes/lab_prasj/BIG_Flo/kassie_bd_rhapsody/harmonised_scaled_cleaned_sce.rds")


example_contig_annotation <- read.csv("/Volumes/fbach/vac63c/first_all_contig_annotations.csv")
immdata_10x <- repLoad("/Volumes/fbach/vac63c/first_all_contig_annotations.csv")
# [1] "barcode"               "is_cell"               "contig_id"             "high_confidence"       "length"               
# [6] "chain"                 "v_gene"                "d_gene"                "j_gene"                "c_gene"               
# [11] "full_length"           "productive"            "fwr1"                  "fwr1_nt"               "cdr1"                 
# [16] "cdr1_nt"               "fwr2"                  "fwr2_nt"               "cdr2"                  "cdr2_nt"              
# [21] "fwr3"                  "fwr3_nt"               "cdr3"                  "cdr3_nt"               "fwr4"                 
# [26] "fwr4_nt"               "reads"                 "umis"                  "raw_clonotype_id"      "raw_consensus_id"     
# [31] "exact_subclonotype_id"
# 
# data.frame("barcode", "is_cell", "contig_id", "high_confidence", "length",
#  "chain","v_gene", "d_gene", "j_gene", "c_gene",
# "full_length","productive", "fwr1","fwr1_nt","cdr1",
# "cdr1_nt","fwr2","fwr2_nt","cdr2","cdr2_nt",
# "fwr3","fwr3_nt","cdr3","cdr3_nt","fwr4",
# "fwr4_nt","reads","umis","raw_clonotype_id","raw_consensus_id",
# "exact_subclonotype_id")

list_repertoire <- readRDS("/Volumes/lab_prasj/BIG_Flo/kassie_bd_rhapsody/list_repertoire.RDS")
long_repertoire <- do.call(rbind, list_repertoire)
colnames(long_repertoire)[grep("cdr3_aa1", colnames(long_repertoire))]="CDR3.aa"

shorter_repertoire <- long_repertoire %>%
  filter(CTaa!="_")%>%
  group_by(sample_id, CTaa)%>%
  add_count(name = "Clones")%>%
  ungroup()%>%
  group_by(sample_id)%>%
  filter(!duplicated(CTaa))%>%
  mutate("donor"=substr(sample_id, 0,2),
         "stim"=substr(sample_id, nchar(sample_id)-3,nchar(sample_id)))%>%
  ungroup()
  
  
list_repertoire <- split(shorter_repertoire, shorter_repertoire$sample_id)

metadata <- shorter_repertoire%>%
  reframe(sample_id, donor, stim)%>%
  filter(!duplicated(sample_id))%>%
  mutate(stim=ifelse(stim=="n_d0", "d0", stim))

div_chao <- repDiversity(list_repertoire, "chao1")

# Hill numbers
div_hill <- repDiversity(list_repertoire, "hill")

# D50
div_d50 <- repDiversity(list_repertoire, "d50")

# Ecological diversity measure
div_div <- repDiversity(list_repertoire, "div")

#gini
div_gini <- repDiversity(list_repertoire, "gini")

p1 <- vis(div_chao, .by = "stim", meta = metadata)
p2 <- vis(div_chao, ,.)

p4 <- vis(div_d50)
p5 <- vis(div_d50, .by = "Status", .meta = first_clonos$meta)
p6 <- vis(div_div)

kmers <- getKmers(list_repertoire, 5)
p1 <- vis(kmers, .head = 5)
p2 <- vis(kmers, .head = 10)
p3 <- vis(kmers, .head = 30)

(p1 + p2) / p3

kp <- kmer_profile(kmers[[1]])
p1 <- vis(kp)
p2 <- vis(kp, .plot = "seq")

p1 + p2
