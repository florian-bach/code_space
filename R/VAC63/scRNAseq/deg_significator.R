first_deg <- read.csv("/Volumes/fbach/vac63c/final_revision/edger_first_DE_lrt.csv", row.names = 1)
third_deg <- read.csv("/Volumes/fbach/vac63c/final_revision/edger_third_DE_lrt.csv", row.names = 1)

T6_deg <- read.csv("/Volumes/fbach/vac63c/final_revision/T6_DE_lrt.csv")


first_deg %>%
  filter(p_val_adj < 0.05)%>%
  group_by(cell_type) -> sig_first

third_deg %>%
  filter(p_val_adj < 0.05)%>%
  group_by(cell_type) -> sig_third

T6_deg %>%
  filter(p_val_adj < 0.05)%>%
  group_by(cell_type) -> sig_t6d

write.csv(sig_first, "~/postdoc/edinburgh/scRNAseq/scg/sig_edger_first_DE_lrt.csv")
write.csv(sig_third, "~/postdoc/edinburgh/scRNAseq/scg/sig_edger_third_DE_lrt.csv")
write.csv(sig_t6d, "~/postdoc/edinburgh/scRNAseq/scg/sig_T6_DE_lrt.csv")





complicated_data %>%
  group_by(first_inf_before_six)
