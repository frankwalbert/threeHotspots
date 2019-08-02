# this code creates the epistasis plot between the OAF1 & OLE1 eQTLs from the Albert, Bloom, et al eLife data

library(tidyverse)

# get expression data; these and all other objects below can be obtained from Albert, Bloom, et al eLife 2018
load("gbatch_fact.RData")
load("ODcov.RData")
load("log2_t.tpm.matrix.RData")

# genotypes
load("gdata_42k.RData")

# gene annotations
geneAnnotation = read.table("ensemblGenes_ensembl83_160307_MOD.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE)
rownames(geneAnnotation) = geneAnnotation[,1]
allNames <- geneAnnotation[, "geneName"]
names(allNames) <- geneAnnotation[,1]
allNamesInv <- names(allNames)
names(allNamesInv) <- allNames

# cleaned-up expression data
phenoBatchOD = apply(t.tpm.matrix, 2, function(y){scale(lm(y ~ gbatch.fact + OD.cov)$res)})
rownames(phenoBatchOD) <- sapply(rownames(t.tpm.matrix), function(x){strsplit(x, "-")[[1]][1]})

# plot the expression levels of a given pair
# (function is from Albert/Bloom eLife paper)
plotOnePair <- function(gene, marker1, marker2){
    plotList <- list(
    phenoBatchOD[rownames(gdata)[gdata[,marker1] == -1 & gdata[,marker2] == -1], gene],
    phenoBatchOD[rownames(gdata)[gdata[,marker1] == -1 & gdata[,marker2] == 1], gene],
    phenoBatchOD[rownames(gdata)[gdata[,marker1] == 1 & gdata[,marker2] == -1], gene],
    phenoBatchOD[rownames(gdata)[gdata[,marker1] == 1 & gdata[,marker2] == 1], gene]
    )
    ggDat <- data.frame(expression = unlist(plotList), genotype=rep(c("BB", "BR", "RB", "RR"), times = sapply(plotList, length)))
    p <- ggplot(ggDat, aes(x = genotype, y = expression)) + geom_boxplot(color="darkblue", outlier.shape=NA) + theme_light() + labs(x="Genotype", y="Normalized expression level", title=paste(gene, allNames[gene], sep=" - "), subtitle=paste(marker1, marker2, sep=" and ")) +
    geom_jitter(width=0.2, alpha=0.4)
    print(p)
    #boxplot(plotList, names=c("BB", "BR", "RB", "RR"), main=paste(noAdditiveNoHotspot[i, "trait"], allNames[noAdditiveNoHotspot[i, "trait"]], sep=" - "))
}

# FAA4 trait, OAF1 & OLE1 markers
pdf("epistasis_FAA4_OAF1_OLE1.pdf")
    plotOnePair("YMR246W", "chrI:48890_C/T", "chrVII:397449_G/A")
dev.off()

# tests:
testOnePair <- function(gene, marker1, marker2){
    m1 <- gdata[,marker1]
    m2 <- gdata[,marker2]
    trait <- phenoBatchOD[, gene]
    h1 <- lm(trait ~ m1 + m2 + m1:m2)
    print(summary(h1))
    h0 <- lm(trait ~ m1 + m2)
    anova(h0, h1)
}

# FAA4
testOnePair("YMR246W", "chrVII:397449_G/A", "chrI:48890_C/T")
#Res.Df    RSS Df Sum of Sq      F    Pr(>F)
#1   1009 593.25
#2   1008 586.48  1    6.7675 11.631 0.0006742 ***
