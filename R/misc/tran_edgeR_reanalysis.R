library(edgeR)


urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE52166", "file=GSE52166_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")
ref_tbl <- tbl 

# load gene annotations 
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID


gsms <- paste0("XXXXXXXXXXXXXXXX00X2X2X0X0X2X2X01X2X1X0X2X2X1XXXX0",
               "XXX0X1X0X0X0X2X0XXX2X0X2X1X1X2X2X1X1X2X0X1X2XXX12X",
               "1X1XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
tbl <- tbl[ ,sel]


# group membership for samples
gs <- factor(sml)
groups <- make.names(c("immune","early fever","delayed fever"))
levels(gs) <- groups
sample_info <- data.frame(Group = gs, row.names = colnames(subset_tbl))

early_and_immune <- tbl[,gs%in%c("immune", "early.fever")]

# d <- DGEList(counts=tbl, group=factor(gs), genes = annot)
d <- DGEList(counts=early_and_immune, group=subset(gs, gs%in%c("immune", "early.fever")), genes = annot)

#filter
# apply(d$counts, 2, sum) # total gene counts per sample

keep <- rowSums(cpm(d)>4) >= 10
d <- d[keep,]
dim(d)

#normalize
d_norm <- calcNormFactors(d)

topVarGenes <- head(order(rowVars(cpm(d_norm)), decreasing = TRUE), 6743)
filtered_d <- d_norm[topVarGenes,]

# filtered_d <- estimateCommonDisp(filtered_d, verbose=TRUE)
# et <- exactTest(filtered_d)
# sig_tags <- topTags(et, p.value = 0.05)

metadata <- data.frame("id"=1:ncol(early_and_immune),
                       "group"=as.character(subset(gs, gs%in%c("immune", "early.fever"))))

design <- model.matrix(~metadata$id+metadata$group)

filtered_d <- estimateGLMCommonDisp(filtered_d, design, verbose=TRUE)
filtered_d <- estimateGLMTrendedDisp(filtered_d, design)
filtered_d <- estimateGLMTagwiseDisp(filtered_d, design)

fit <- glmFit(filtered_d, design)

lrt <- glmLRT(fit)
tags <- topTags(lrt,p.value = 0.05)
