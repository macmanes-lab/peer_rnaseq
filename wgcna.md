```
#WGCNA

###Setup and MOdule 1

```r

source("http://bioconductor.org/biocLite.R")
biocLite(c("gplots", "edgeR", "AnnotationDbi", "impute", "GO.db", "preprocessCore", "KEGG.db", "topGO", "packageNames", "hgu133a.db", "hgu95av2.db", "annotate", "hgu133plus2.db", "SNPlocs.Hsapiens.dbSNP.20100427", "minet", "OrderedList"))
install.packages(c("matrixStats", "Hmisc", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival", "WGCNA")) 
library("WGCNA")
library(edgeR)
library(gplots)
allowWGCNAThreads()
orgCodes = c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg");
orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6));
packageNames = paste("org.", orgCodes, orgExtensions, ".db", sep="");
data <- read.csv("~/Dropbox/peer_rnaseq/peer.kidney.counts", row.names="Name")
phenotype <- factor(c(rep("dry",19),rep("wet",19)))
design <- model.matrix( ~0+phenotype )
keep <- rowSums(cpm(data) > 5) >= 19
data <- data[keep,]
object <- DGEList(counts=data, group=phenotype)
object <- estimateGLMCommonDisp( object, design )
object <- estimateGLMTrendedDisp( object, design )
object <- estimateGLMTagwiseDisp( object, design )
object <- calcNormFactors(object)
object <- estimateCommonDisp(object)
object <- estimateTagwiseDisp(object)
fit <- glmFit( object, design )





datExpr1 = as.data.frame(t(object$pseudo.counts))
dim(datExpr1)
sampleTree = hclust(dist(datExpr1), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
mouseSamples <- rownames(datExpr1)
phys <- read.csv("~/Dropbox/peer_rnaseq/phys4.csv")
traitRows = match(mouseSamples, phys$ID)
datTraits <- phys[traitRows, -1]
rownames(datTraits) = phys[traitRows, 1]

sampleTree2 = hclust(dist(datExpr1), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)


plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

save(datExpr1, datTraits, file = "mouse-01-dataInput.RData")



```


### Module 2

```r
lnames = load(file = "mouse-01-dataInput.RData");
powers = c(c(1:10), seq(from = 2, to=30, by=2))
sft = pickSoftThreshold(datExpr1, powerVector = powers, verbose = 5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net = blockwiseModules(datExpr1, power = 5,
                       TOMType = "unsigned", minModuleSize = 40,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "mouseSampleTOM",
                       verbose = 3)


mergedColors = labels2colors(net$colors)

plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.003,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
file = "mouse-02-networkConstruction-auto.RData")



```



## Make  Cool plot of trait/module relationship

```r
lnames = load(file = "mouse-01-dataInput.RData");
lnames = load(file = "mouse-02-networkConstruction-auto.RData")
nGenes = ncol(datExpr1);
nSamples = nrow(datExpr1);
MEs0 = moduleEigengenes(datExpr1, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mfcol = c(1,1))

png("heatmap.png", width = 6.5, height = 4.5, units = 'in', res = 300)
par(mai = c(.8, 1.4, .5, .3), mfrow = c(1,1), bty = 'n')

labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
dev.off()

```




###GS v MM Correlation

```r

weight = as.data.frame(datTraits$Sodium);
names(weight) = "weight"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr1, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr1, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

module = "tan"
column = match(module, modNames);
moduleGenes = moduleColors==module;
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for body weight",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

```


###GO Term Enrichment

```r

GO term enrichnment stuff:

blast to mus to get entrex ids
cat entresids | awk '{print $2}' | paste -d, - list2 > import.txt


annot <- read.csv("~/Dropbox/peer_rnaseq/import.txt")
probes = names(datExpr1)
probes2annot = match(probes, object$genes$entrez_id)
allLLIDs = object$genes$entrez_id[probes2annot]
intModules = c("red", "yellow")

for (module in intModules)
{
modGenes = (moduleColors==module)
modLLIDs = allLLIDs[modGenes];
fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
write.table(as.data.frame(modLLIDs), file = fileName,
row.names = FALSE, col.names = FALSE)
}

fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
row.names = FALSE, col.names = FALSE)

GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "mouse", nBestP = 20)

tab = GOenr$bestPTerms[[4]]$enrichment
names(tab)
write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)
keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
options(width=95)
screenTab

```


### GO for diff expression
###POS
```bash

annot <- read.csv("~/Dropbox/peer_rnaseq/import.txt")
posdiffexp <- read.table("~/Dropbox/peer_rnaseq/positive.txt", quote="\"", comment.char="")
posdiffexp2annot = match(posdiffexp$V1, annot$transid)
pos = annot$entrez[posdiffexp2annot]
GOpos = GOenrichmentAnalysis(posdiffexp$V1, pos, organism = "mouse", nBestP = 20)
tab = GOpos$bestPTerms[[4]]$enrichment
names(tab)
write.table(tab, file = "GOpos.csv", sep = ",", quote = TRUE, row.names = FALSE)
```

###NEG

```bash

annot <- read.csv("~/Dropbox/peer_rnaseq/import.txt")
negdiffexp <- read.table("~/Dropbox/peer_rnaseq/negative.txt", quote="\"", comment.char="")
negdiffexp2annot = match(negdiffexp$V1, annot$transid)
neg = annot$entrez[negdiffexp2annot]
GOneg = GOenrichmentAnalysis(negdiffexp$V1, pos, organism = "mouse", nBestP = 20)
tab = GOneg$bestPTerms[[4]]$enrichment
names(tab)
write.table(tab, file = "GOneg.csv", sep = ",", quote = TRUE, row.names = FALSE)
```
