pras_metadata2 <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/big_data/fixed_musical_metadata_final.csv")
musical_rnaseq <- read.csv("~/postdoc/stanford/clinical_data/MUSICAL/MUSICSL_RNA-Seq_Pilot.csv")

musical_rnaseq <- musical_rnaseq %>%
  mutate(id=combined_id, date=combined_date)

combo <- inner_join(musical_rnaseq, pras_metadata2, by = c("id", "date"))

# missing <- anti_join(musical_rnaseq, combo, by = c("combined_id", "combined_date"))

write.csv(combo, "~/postdoc/stanford/clinical_data/MUSICAL/new_rna-seq_metadata.csv", row.names = F)
