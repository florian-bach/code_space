setwd("/home/flobuntu/PhD/RNAseq/vac69a/all/xls/gene_lists/")

up_dod <- scan(file = "DoD_Baseline_UP_significant_symbol_only.txt", sep="\t", what = "")
up_t6 <- scan(file = "T6_Baseline_UP_significant_symbol_only.txt", sep="\t", what = "")

#genes up at T6 that aren't up on DOD; 172/280 ~ 60%
t6_specific <- up_t6[!up_t6 %in% up_dod]
t6_specific <- t6_specific[order(t6_specific)]
t6_specific

vivax_t6_faves <- c("CD38", "CTLA4", "CXCR6", "GZMA", "ICOS", "IL21", "JCHAIN", "MPO", "MKI67")

vivax_t6_faves[!vivax_t6_faves %in% t6_specific]
# "CD38"  "MKI67"
vivax_t6_faves[vivax_t6_faves %in% t6_specific]
#"CTLA4" "CXCR6" "GZMA"  "ICOS"  "IL21"  "MPO"

write.table(t6_specific, "t6_specific_up.txt", sep = "\t", quote = F, row.names = F, col.names = F)


falci_up_dod <- scan("~/PhD/RNAseq/vac63c/vac63a_b_sig_dod_baseline_up.txt", sep="\t", what = "")
