library(ggplot2)
library(affy)
library(limma)
library(annotate)
library(statmod)
library(hgu133a.db)

setwd("/Users/s1249052/PhD/transcriptomics/children/")

#### read in all CEL files in directory

files_list <-list.files(pattern=".CEL")

MyData <- ReadAffy(filenames=files_list[c(1:5, 11:15)])

#### make metadata datagrame
pd <- data.frame(population = factor(rep(c("Control", "Severe", "Uncomplicated"), each=5)), replicate = rep(seq(1,5), 3))
rownames(pd) <- sampleNames(MyData)

#### feature info
vl <- list(population = "control is healthy, severe is severe anaemia, uncomplicated is non-severe malaria", replicate = "1, 2, 3, 4, 5 arbitrary numbering, each chip contains 4 kids' RNA")

#### store as phenotpyic data in affy object

phenoData(MyData) <- AnnotatedDataFrame(pd, varMetadata=as.data.frame(c(vl[[1]], vl[[2]])), dimLabels=c("rowNames", "columnNames"))

#### convert probe level data to express measure
eset <- rma(MyData)

# #### convert to easy, tasty matrix
e <- exprs(eset)
# 
# #### make log fold change
# control <- which(eset$population == "Control")
# mild <- which(eset$population == "Uncomplicated")
# severe <- which(eset$population == "Severe")
# 
# d <- rowMeans(e[, control]) - rowMeans(e[, severe])
# 
# #### use genefilter to apply t test to each row in dataset
# tt <- genefilter::rowttests(subset(eset, eset$population %in% c("Severe", "Control")), "population")
# tt <- genefilter::rowttests(eset, "population")
# 
# 
# #### make volcano plot
# lod <- -log10(tt$p.value)
# 
# plot(d, lod, cex = 0.25, main = "Volcano plot for $t$-test", xlab="Fold Change", ylab="log p value")+
#   abline(h = 2)
# 
# 
#### test for differential expression with limma

##  
design <- model.matrix(~0 + eset$population)
colnames(design) <- c("Control", "Severe", "Uncomplicated")

### more stuff to help limma...

blocks=rep(1:3,each=5)
corfit <- duplicateCorrelation(eset, design=design, block = blocks)

fit <- lmFit(eset, design, block = blocks, correlation = )

healthy_mild <- makeContrasts(Control - Uncomplicated, levels=design)
mild_severe <- makeContrasts(Uncomplicated - Severe, levels=design)
healthy_severe <- makeContrasts(Control - Severe, levels=design)

fit2 <- contrasts.fit(fit, contrasts=healthy_mild)
fit3  <- eBayes(fit2)

fit4 <- contrasts.fit(fit, contrasts=mild_severe)
fit5  <- eBayes(fit34)
                
fit6 <- contrasts.fit(fit, contrasts=healthy_severe)
fit7  <- eBayes(fit6)



## extract significant genes

mild_results <- decideTests(fit2, adjust.method = "BH")
severe_mild_results <- decideTests(fit3, adjust.method = "BH")
severe_results <- decideTests(fit4, adjust.method = "BH")

mild_degs <- as.data.frame(subset(mild_results, mild_results@.Data !=0))
severe_mild_degs <- as.data.frame(subset(severe_mild_results, severe_mild_results@.Data !=0))
severe_degs <- as.data.frame(subset(severe_results, severe_results@.Data !=0))

#####   1255 DEGs

## make top n list
deg <- topTable(fit2, coef=1, n=1255, adjust = "BH")

table(rownames(degs) %in% rownames(deg))
### TRUE

### get real gene names
deg$symbol <- getSYMBOL(rownames(deg), "hgu133a.db")

### get rid of duplicate genes

non_redundant <- as.data.frame(deg %>%
                                 group_by(symbol) %>%
                                 arrange(abs(logFC)) %>%
                                 slice(1)
)

###  volcano plot

(test <- ggplot(non_redundant, aes(x=logFC, y=-log10(P.Value), text=symbol))+
    geom_point(aes(color=logFC))+
    scale_color_gradientn(colors=colorspace::diverge_hcl(5, "Blue-Red 2"))+
    # ggrepel::geom_text_repel(data=subset(non_redundant, non_redundant$Symbol %in% tab$symbol), aes(label=Symbol), color="red")+
    #geom_text(data=subset(non_redundant, non_redundant$Symbol %in% tab$symbol), aes(label=Symbol), color="red")+
    theme_bw()+
    geom_hline(yintercept = 2))+
  theme(legend.title.align = 0.5)

plotly::ggplotly(test, tooltip="text")






gene_names <- degs$symbol
#### 
tab <- topTable(ebayes, coef = 2, adjust.method = "BH", n = 100)
genenames <- rownames(tab)

#### convert probe IDs to gene names

annotation(eset) #  "hgu133a"

ll <- getEG(genenames, "hgu133a.db")
sym <- getSYMBOL(genenames, "hgu133a.db")

tt$symbol <- getSYMBOL(rownames(tt), "hgu133a.db")

tab <- data.frame(sym, signif(tab[, -1], 3))

#### ggplot biz

df <- as.data.frame(d)
df$symbol <- getSYMBOL(rownames(df), "hgu133a.db")

colnames(df) <- c("log2_FC", "Symbol")
df$p_value <- tt[match(tt$symbol, df$Symbol), "p.value"]

tab <- topTable(ebayes, coef = 2, adjust.method = "BH", n = 50)
tab$symbol <- getSYMBOL(rownames(tab), "hgu133a.db")

non_redundant <- as.data.frame(df %>%
  group_by(Symbol) %>%
  arrange(abs(p_value)) %>%
  slice(1)
)


(test <- ggplot(non_redundant, aes(x=log2_FC, y=-log10(p_value), text=Symbol))+
  geom_point(data=subset(non_redundant, -log10(p_value)>2), aes(color=log2_FC, alpha=-log10(p_value)))+
  geom_point(aes(alpha=-1))+
  scale_color_gradientn(colors=colorspace::diverge_hcl(5, "Blue-Red 2"))+
  # ggrepel::geom_text_repel(data=subset(non_redundant, non_redundant$Symbol %in% tab$symbol), aes(label=Symbol), color="red")+
  #geom_text(data=subset(non_redundant, non_redundant$Symbol %in% tab$symbol), aes(label=Symbol), color="red")+
  theme_bw()+
  geom_hline(yintercept = 2))+
  theme(legend.title.align = 0.5)

plotly::ggplotly(test, tooltip="text")



htmlpage(ll, filename = "GeneList1.html", title = "HTML report", othernames = tab, table.head = c("Locus ID", colnames(tab)),table.center = TRUE)

         