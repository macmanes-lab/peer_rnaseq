```

source("http://bioconductor.org/biocLite.R")
biocLite(c("Rgraphviz", "topGO", "edgeR", "tximport", "readr", "locfit", "statmod", "gplots", "AnnotationDbi", "impute", "GO.db", "preprocessCore", "KEGG.db", "topGO", "packageNames", "hgu133a.db", "hgu95av2.db", "annotate", "hgu133plus2.db", "SNPlocs.Hsapiens.dbSNP.20100427", "minet", "OrderedList", "mygene"))
install.packages(c("matrixStats", "Hmisc", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival", "WGCNA", "mygene", "lsmeans", "doBy"))
orgCodes = c("Hs", "Mm","Gg")
orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6))
packageNames = paste("org.", orgCodes, orgExtensions, ".db", sep="")
library(edgeR)
library(tximport)
library(readr)
library(locfit)
library(statmod)
library(gplots)
library(topGO)
library(Rgraphviz)
library(doBy)

dir <- setwd('/Users/macmanes/Dropbox/peer_rnaseq/')
tx2gene <- read.csv(file.path(dir, "peer.tximport"), header=T)

all_samples <- read.table(file.path(dir, "counts/samples.list"), header=F)
all_files <- file.path(dir, "counts/", all_samples$V1, "quant.sf")
names(all_files) <- paste0(all_samples$V1)
txi.stress <- tximport(all_files, type = "salmon", tx2gene = tx2gene, reader = read_tsv)

cts <- txi.stress$counts
normMat <- txi.stress$length
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))

phenotype <- factor(c(rep("dry",17),rep("wet",18), rep("wild",10)))

write.table(cts, file="peerkidney_tximport.counts", sep = "," , row.names = TRUE)

peer_for_analysis <- read.csv("/Users/macmanes/Dropbox/peer_rnaseq/peerkidney_tximport.counts", row.names="Name")


data <- read.csv("/Users/macmanes/Dropbox/peer_rnaseq/peer.tximport.kidney.counts1", row.names="name")
phenotype <- factor(c(rep("dry",17),rep("wet",18), rep("wild",10)))
design <- model.matrix( ~0+phenotype )


object <- DGEList(counts=data[,4:48], gene=data[,1:3], group=phenotype)
object <- estimateGLMCommonDisp( object, design )
object <- estimateGLMTrendedDisp( object, design )
object <- estimateGLMTagwiseDisp( object, design )
object <- calcNormFactors(object)
object <- estimateCommonDisp(object)
object <- estimateTagwiseDisp(object)
fit <- glmFit( object, design )
contrast_W_v_D <- glmLRT( fit, contrast=makeContrasts( phenotypewet-phenotypedry, levels=design ) )
summary(decideTestsDGE(contrast_W_v_D, adjust.method="fdr", p.value=0.01))

#   [,1] 
#-1   252
#0  13231
#1    213

detags <- rownames(contrast_W_v_D)[as.logical(decideTestsDGE(contrast_W_v_D, adjust.method="fdr", p.value=0.01))]
topTags(contrast_W_v_D, n=30)
```
