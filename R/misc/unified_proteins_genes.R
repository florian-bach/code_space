ihm_genes <- read.csv("~/postdoc/stanford/misc/unified_immune_metric/sig_genes.csv", header=TRUE)

long_ihm_genes <- ihm_genes%>%
  pivot_longer(cols=paste("Feature", 1:98, sep=""), values_to = "gene")%>%
  select(-name)
  
long_ihm_genes[c(grep("^RPS", long_ihm_genes$gene),grep("^RPL", long_ihm_genes$gene)),]

ihm_proteins <- read.csv("~/postdoc/stanford/misc/unified_immune_metric/sig_proteins.csv")
table(ihm_proteins$EntrezGeneSymbol %in% nulisa$targetName)
overlap <- data.frame(nulisa$targetName[nulisa$targetName %in% ihm_proteins$EntrezGeneSymbol])
colnames(overlap)="symbol"
write.table(overlap, "~/postdoc/stanford/misc/unified_immune_metric/nulisa_sig_somalogic_overlap.csv", row.names = FALSE)
