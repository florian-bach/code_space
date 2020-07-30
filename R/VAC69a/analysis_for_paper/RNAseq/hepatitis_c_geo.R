#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE7123", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))
# group names for all samples
gsms <- paste0("102X102XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXX102X102X102XXX102XXX102XXX102XXX102XXX102XXX10",
               "2XXX102XXX102XXX102XXX102XXX102XXX102XXX102XXX102X",
               "X102XXX102XXX102XXX102XXX102XXX102XXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX102XX",
               "X102XXX102XX102XXX102XXXXXXXXXXXXXXXXXXXXXXXXXX")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G2-G0, G1-G0, G2-G1, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=2500)
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","G1...G0", "G2...G0", "Gene.symbol","Gene.title"), adj.P.Val<0.05)
#write.table(tT, file=stdout(), row.names=F, sep="\t")








vivax_dod <- scan("gene_lists/DoD_Baseline_ALL_significant_symbol_only.txt", what = "")

vivax_hepc_overlap <- subset(vivax_dod, vivax_dod %in% tT$Gene.symbol)




falci_dod <- scan("~/PhD/RNAseq/vac63c/vac63a_b_sig_dod_baseline_all.txt", what = "")

falci_hepc_overlap <- subset(falci_dod, falci_dod %in% tT$Gene.symbol)




