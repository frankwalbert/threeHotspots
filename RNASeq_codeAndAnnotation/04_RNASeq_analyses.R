# run in R/3.6.0
#library(ggplot2)
#library(gplots)
library(tidyverse)
library(DESeq2)
library(readxl)
library(readr)
library(qvalue)
#library(pheatmap)
#library(RColorBrewer)
library(sva)


readFolder <- "YOUR_OUTPUT_FOLDER/reads/"
samples <- dir(readFolder)[!str_detect(dir(readFolder), "trimmomatic")]
samples <- samples[!str_detect(samples, "QC")]
wells <- sapply(samples, function(x){strsplit(x, "_")[[1]][3]})
genotype <- sapply(samples, function(x){paste(strsplit(x, "_")[[1]][5:length(strsplit(x, "_")[[1]])], collapse="_")})

geneAnnotation = read.table("ensemblGenes_ensembl83_160307_MOD.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE) # provided
rownames(geneAnnotation) = geneAnnotation[,1]
allNames <- geneAnnotation[, "geneName"]
names(allNames) <- geneAnnotation[,1]
allNames[which(allNames == "")] <- names(allNames)[which(allNames == "")]


load("B.Forward.RData") # provided; from Albert/Bloom eLife 2018
load("B.cis.Forward.RData") # provided; from Albert/Bloom eLife 2018
hotspots <- c("chrI:48890_C/T", "chrVII:397449_G/A", "chrIV:215523_T/C")
phenotypers <- c("YMR246W", "YMR246W", "YHR094C")
localGenes <- list("YAL051W", c("YGL056C", "YGL055W"), "YDL138W")
specialGenes <- list(c(), c(), c("YMR300C", "YGL234W", "YLR359W", "YAR015W", "YOR128C", "YDR408C")) # for highlighting additional genes, e.g. the ADE genes in RGT2

# add local effects for OLE1 & SDS23
B["chrVII:397449_G/A", "YGL056C"] <- B.cis[which(abs(B.cis[,"YGL056C"]) == max(abs((B.cis[,"YGL056C"])))), "YGL056C"]
B["chrVII:397449_G/A", "YGL055W"] <- B.cis[which(abs(B.cis[,"YGL055W"]) == max(abs((B.cis[,"YGL055W"])))), "YGL055W"]


# CAREFUL: sampleInfo is initially in a different order than the samples vector; be sure to always address via the vector
# or better, reorder sampleInfo as we do below
sampleInfo <- data.frame(read_excel("RNASeq_ALLinfo_RPFA0002_190513_mod.xlsx", col_names=TRUE)) # provided
rownames(sampleInfo) <- sampleInfo$fileName
# reorder to match "samples"
sampleInfo <- sampleInfo[samples,]
# make some batches factors (not numerics)
sampleInfo$Harvest_batch <- as.factor(sampleInfo$Harvest_batch)
sampleInfo$RNA_isolation_batch <- as.factor(sampleInfo$RNA_isolation_batch)
#OD <- sampleInfo[samples,"OD600"]


geneCounts <- sapply(samples, function(x){
    thisDat <- read.table(paste(readFolder, x, "/abundance.tsv", sep=""), sep="\t", head=TRUE)
    thisDat[,"est_counts"]
})
rownames(geneCounts) <- read.table(paste(readFolder, samples[1], "/abundance.tsv", sep=""), sep="\t", head=TRUE, stringsAsFactors=FALSE)$target_id
# write out for paper/GEO
# write.table(geneCounts, file="geneCounts4paper.txt", quote=FALSE, sep="\t")


TINPerGene <- sapply(wells, function(x){
    theseTINs <- read.table(paste(readFolder, "QC/RSeqQC/pseudoalignments_", x, ".tin.xls", sep=""), head=TRUE, stringsAsFactors=FALSE)
    rownames(theseTINs) <- theseTINs[,1]
    theseTINs <- theseTINs[rownames(geneCounts),]
    theseTINs[,"TIN"]
})
rownames(TINPerGene) <- rownames(geneCounts)


meanTINPerGene <- rowMeans(TINPerGene)
# effective length is annotation based, and therefore can be pulled from any one sample
effectiveLengthPerGene <- read.table(paste(readFolder, samples[1], "/abundance.tsv", sep=""), sep="\t", head=TRUE)[,"eff_length"]
names(effectiveLengthPerGene) <- read.table(paste(readFolder, samples[1], "/abundance.tsv", sep=""), sep="\t", head=TRUE)[,"target_id"]


geneCountsFiltered <- geneCounts[rowMin(round(geneCounts)) > 0 & apply(TINPerGene, 1, min) > 0 & effectiveLengthPerGene > 0,]
# down to 5400 from 6713
# the "round" is important, otherwise get Inf norm factors for three genes with "counts" very close to, but not exactly, zero

# write.table(geneCountsFiltered, file="geneCountsFiltered4paper.txt", quote=FALSE, sep="\t")


############################
# SVA

# need a DESeq model first to cull expression values from
# can be any in terms of factors, just be sure to below include the ones we want
thisColDat <- data.frame(genotype, sampleInfo[,c("OD600", "Harvest_batch", "RNA_isolation_batch", "RIN_score", "RIN_Conc", "Qubit_Conc", "library_Prep_batch")])
dds <- DESeqDataSetFromMatrix(countData = round(geneCountsFiltered), colData = thisColDat, design = ~ Harvest_batch + RNA_isolation_batch + library_Prep_batch + genotype)
dds <- DESeq(dds, betaPrior=TRUE)

# SVA
dat4SVA  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat4SVA) > 1
dat4SVA  <- dat4SVA[idx, ]
mod  <- model.matrix(~ Harvest_batch + RNA_isolation_batch + library_Prep_batch + genotype, colData(dds))
mod0 <- model.matrix(~ Harvest_batch + RNA_isolation_batch + library_Prep_batch, colData(dds))
svseq <- svaseq(dat4SVA, mod, mod0, n.sv = 2)

ddssva <- dds
thisColDat_SVA <- thisColDat
thisColDat_SVA$SV1 <- svseq$sv[,1]
thisColDat_SVA$SV2 <- svseq$sv[,2]

ddssva <- DESeqDataSetFromMatrix(countData = round(geneCountsFiltered), colData = thisColDat_SVA, design = ~ Harvest_batch + RNA_isolation_batch + library_Prep_batch + SV1 + SV2 + genotype)
ddssva <- DESeq(ddssva, betaPrior=TRUE)

# what does this look like in PCA?
pdf("PCA_SVA.pdf")
vsd <- vst(ddssva, blind=FALSE)
plotPCA(vsd, intgroup=c("SV1"))
plotPCA(vsd, intgroup="SV2")
plotPCA(vsd, intgroup="genotype")
plotPCA(vsd, intgroup="OD600")
plotPCA(vsd, intgroup="Harvest_batch") # <= this looks fairly important
plotPCA(vsd, intgroup="RNA_isolation_batch")
plotPCA(vsd, intgroup="RIN_score")
plotPCA(vsd, intgroup="RIN_Conc")
plotPCA(vsd, intgroup="Qubit_Conc")
plotPCA(vsd, intgroup="library_Prep_batch")
dev.off()



##################
# DESeq2

# initially just do each vs the FAA4-derived WT
otherGTs <- unique(genotype[genotype != "FAA4-GFP_to_BY"])
# if test against the RGT2-derived WT:
# otherGTs <- unique(genotype[genotype != "rgt2KanMX_to_BY"])


DESeq_results_pairwise <- lapply(otherGTs, function(thatGT){
    #res <- results(ddssva, contrast=c("genotype", thatGT, "FAA4-GFP_to_BY"))
    # if test against the RGT2-derived WT:
    res <- results(ddssva, contrast=c("genotype", thatGT, "rgt2KanMX_to_BY"))
        res <- res[order(res$pvalue),]
        res <- res[!is.na(res$pvalue),]
        res
})
names(DESeq_results_pairwise) <- otherGTs
# sort the DESeq results in same order as the hotspots!!!
DESeq_results_pairwise <- DESeq_results_pairwise[c("OAF1_L63S", "OLE1_FAR_RM", "RGT2_V539I", "rgt2KanMX_to_BY")]

#save(DESeq_results_pairwise, file="R_DESeq_results_pairwise_allKnownCovariates_2SVAs_no_RINScore_RINConc_QbitConc_OD_190524.RData")
# if test against the RGT2-derived WT:
#save(DESeq_results_pairwise, file="R_DESeq_results_pairwise_allKnownCovariates_2SVAs_no_RINScore_RINConc_QbitConc_OD_RGT2WT_190524.RData")

# glue on hotspot effects; write out for paper data
load("R_DESeq_results_pairwise_allKnownCovariates_2SVAs_no_RINScore_RINConc_QbitConc_OD_190524.RData")
resultTableForPaperOAF1 <- cbind(DESeq_results_pairwise[["OAF1_L63S"]][intersect(rownames(DESeq_results_pairwise[["OAF1_L63S"]]), colnames(B)),], B["chrI:48890_C/T", intersect(rownames(DESeq_results_pairwise[["OAF1_L63S"]]), colnames(B))])
colnames(resultTableForPaperOAF1)[7] <- "hotspotEffect"
#write.table(resultTableForPaperOAF1, file="resultTableForPaperOAF1.txt", quote=FALSE, sep="\t")

resultTableForPaperOLE1 <- cbind(DESeq_results_pairwise[["OLE1_FAR_RM"]][intersect(rownames(DESeq_results_pairwise[["OLE1_FAR_RM"]]), colnames(B)),], B["chrVII:397449_G/A", intersect(rownames(DESeq_results_pairwise[["OLE1_FAR_RM"]]), colnames(B))])
colnames(resultTableForPaperOLE1)[7] <- "hotspotEffect"
#write.table(resultTableForPaperOLE1, file="resultTableForPaperOLE1.txt", quote=FALSE, sep="\t")

load("R_DESeq_results_pairwise_allKnownCovariates_2SVAs_no_RINScore_RINConc_QbitConc_OD_RGT2WT_190524.RData")
resultTableForPaperRGT2 <- cbind(DESeq_results_pairwise[["RGT2_V539I"]][intersect(rownames(DESeq_results_pairwise[["RGT2_V539I"]]), colnames(B)),], B["chrIV:215523_T/C", intersect(rownames(DESeq_results_pairwise[["RGT2_V539I"]]), colnames(B))])
colnames(resultTableForPaperRGT2)[7] <- "hotspotEffect"
#write.table(resultTableForPaperRGT2, file="resultTableForPaperRGT2.txt", quote=FALSE, sep="\t")




##############################
# compare to hotspot effects

# rename this plot for whether or not horizontal line of genes without a hotspot effect is included
pdf("DE_vs_hotspotEffects_pairwise.pdf", width=15, height=5, useDingbats = FALSE)
par(mfrow=c(1,3))
# or, if combined wts:
    for(hs in 1:3){
        
        #theseQvals <- qvalue(DESeq_results_pairwise[[hs]]$pvalue[DESeq_results_pairwise[[hs]]$pvalue < 1])$q
        
        commonGenes <- intersect(rownames(DESeq_results_pairwise[[hs]]), colnames(B)[B[hotspots[hs],] != 0])
        commonSigGenes <- intersect(rownames(DESeq_results_pairwise[[hs]])[DESeq_results_pairwise[[hs]]$pvalue < 0.05], colnames(B)[B[hotspots[hs],] != 0])
        commonVerySigGenes <- intersect(rownames(DESeq_results_pairwise[[hs]])[DESeq_results_pairwise[[hs]]$padj < 0.1], colnames(B)[B[hotspots[hs],] != 0])
        ##commonQValueGenes <- intersect(rownames(DESeq_results_pairwise[[hs]])[theseQvals <= 0.2], colnames(B)[B[hotspots[hs],] != 0])
        
        # comment out next block if want to *not* show the horizontal line:
        #commonGenes <- intersect(rownames(DESeq_results_pairwise[[hs]]), colnames(B))
        #commonSigGenes <- intersect(rownames(DESeq_results_pairwise[[hs]])[DESeq_results_pairwise[[hs]]$pvalue < 0.05], colnames(B))
        #commonVerySigGenes <- intersect(rownames(DESeq_results_pairwise[[hs]])[DESeq_results_pairwise[[hs]]$padj < 0.1], colnames(B))
        ##commonQValueGenes <- intersect(rownames(DESeq_results_pairwise[[hs]])[theseQvals <= 0.2], colnames(B))

        plot(DESeq_results_pairwise[[hs]][commonGenes,]$log2FoldChange, B[hotspots[hs], commonGenes], col="#00000022",
        xlab="log2 differential expression", ylab="Hotspot effect", main=paste(hotspots[hs], " vs ", "BY", " effects", sep=""))
        abline(h=0, lty=2, col="grey")
        abline(v=0, lty=2, col="grey")
        abline(0, 1, lty=2, col="grey")
        points(DESeq_results_pairwise[[hs]][commonSigGenes,]$log2FoldChange, B[hotspots[hs], commonSigGenes], pch=19, col="#FF000044")
        
        cors <- list(
        cor.test(DESeq_results_pairwise[[hs]][commonGenes,]$log2FoldChange, B[hotspots[hs], commonGenes], method="s"),
        cor.test(DESeq_results_pairwise[[hs]][commonSigGenes,]$log2FoldChange, B[hotspots[hs], commonSigGenes], method="s"),
        try({cor.test(DESeq_results_pairwise[[hs]][commonVerySigGenes,]$log2FoldChange, B[hotspots[hs], commonVerySigGenes], method="s")})
        #try({cor.test(DESeq_results_pairwise[[hs]][commonQValueGenes,]$log2FoldChange, B[hotspots[hs], commonQValueGenes], method="s")})
        )
        #names(cors) <- c("all", "sig", "verySig", "qVal")
        names(cors) <- c("all", "sig", "verySig")
        
        for(i in 1:length(cors)){
            if(class(cors[[i]])=="try-error"){cors[[i]] <- data.frame(est=NA, p.value=NA)}
        }
        
        formattedPs <- sapply(cors, function(x){
            if(is.na(x)){return(NA)}
            if(x$p.value < 0.05){pValFormatted <- formatC(x$p.value, format = "e", digits = 1)}
            else{pValFormatted <- round(x$p.value, 3)}
            pValFormatted
        })
        
        legPos <- c("topleft", "bottomleft")
        names(legPos) <- c("BY", "RM")
        legend(legPos["BY"], box.lty=0, legend=c(
        paste("rho all: ", round(cors[[1]]$est, 4), sep=""),
        paste("p all: ", formattedPs[1], sep=""),
        paste("rho p < 0.05: ", round(cors[[2]]$est, 4), sep=""),
        paste("p p < 0.05: ", formattedPs[2], sep=""),
        paste("rho 10% FDR: ", round(cors[[3]]$est, 4), sep=""),
        paste("p 10% FDR: ", formattedPs[3], sep="")
        #paste("rho 10% qValue: ", round(cors[[4]]$est, 2), sep=""),
        #paste("p 10% qValue: ", formattedPs[4], sep="")
        ))
        
        #points(DESeq_results_pairwise[[hs]][commonQValueGenes,]$log2FoldChange, B[hotspots[hs], commonQValueGenes], col="orange", cex=1.5)
        points(DESeq_results_pairwise[[hs]][commonVerySigGenes,]$log2FoldChange, B[hotspots[hs], commonVerySigGenes], col="purple")
        
        points(DESeq_results_pairwise[[hs]][phenotypers[hs],]$log2FoldChange, B[hotspots[hs], phenotypers[hs]], col="red", cex=2, lwd=2)
        for(j in 1:length(localGenes[[hs]])){
            points(DESeq_results_pairwise[[hs]][localGenes[[hs]][j],]$log2FoldChange, B[hotspots[hs], localGenes[[hs]][j]], col="blue", cex=2, lwd=2)
        }
        if(length(specialGenes[[hs]] > 0)){
            for(j in 1:length(specialGenes[[hs]])){
                points(DESeq_results_pairwise[[hs]][specialGenes[[hs]][j],]$log2FoldChange, B[hotspots[hs], specialGenes[[hs]][j]], col="orange", cex=2, lwd=2)
                text(jitter(DESeq_results_pairwise[[hs]][specialGenes[[hs]][j],]$log2FoldChange), B[hotspots[hs], specialGenes[[hs]][j]], labels=allNames[specialGenes[[hs]][j]], cex=0.7, color="orange")
            }
        }
        text(jitter(DESeq_results_pairwise[[hs]][commonSigGenes,]$log2FoldChange), B[hotspots[hs], commonSigGenes], labels=allNames[commonSigGenes], cex=0.5)
    }
dev.off()


# result tables:
for(i in names(DESeq_results_pairwise)){
    tempTable <- cbind(rownames(DESeq_results_pairwise[[i]]), cbind(allNames[rownames(DESeq_results_pairwise[[i]])], DESeq_results_pairwise[[i]]))
    colnames(tempTable)[1] <- "gene"
    colnames(tempTable)[2] <- "geneName"
    write.table(tempTable, file=paste("DESeq_results_", i, ".txt", sep=""), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
}


####################
# are there any strong B effects that have no DE at all?
hs <- 3
thres=0.4
names(allNames[intersect(rownames(DESeq_results_pairwise[[hs]])[DESeq_results_pairwise[[hs]]$pvalue > thres], colnames(B)[abs(B[hotspots[hs],]) > thres])])
#write.table(intersect(rownames(DESeq_results_pairwise[[hs]]), colnames(B)), file="testedGenes_RGT2.txt", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
# for OAF1, none at B>0.4, ALR1 only at >.3, 20 at >0.2
# putting the 20 through the SGD GO finder gives nothing

# OLE1; 12 at > 0.2; weak enrichment in transposition (!); but then we also have those in the DE only?
# SDS23 and YDR210W-B at >0.4 and 0.3

# RGT2 (with correct WT):
# 46 at 0.2: best GO: "'de novo' IMP biosynthetic process" 6/46 expect 9/7166, 1.5e-9
# 10 at 0.3: ARF2, MRK1, SNF3, MMP1, FCY21, ADE2, ADE5,7, SRF1, YPL113C, ADE13; the three ADE genes again drive "'de novo' IMP biosynthetic process": 3/10, 9/7166, 1e-5
# 2 at 0.4:  "SRF1"

# the ADE thing is interesing. what ARE their results?
DESeq_results_pairwise[[3]]["YLR359W",] # ADE13 p = 0.95, logFC = -0.002
DESeq_results_pairwise[[3]]["YOR128C",] # ADE2 p = 0.78, logFC = -0.01
DESeq_results_pairwise[[3]]["YGL234W",] # ADE5,7 p = 0.82, logFC = -0.01
# actually, it looks like NONE of the ADE genes show up?
# the 46 contain:  "ADE1"   "ADE13"    "ADE2"    "ADE4"  "ADE5,7"    "ADE8", and the final ADE17 shows up in 0.1
# but the RGT2 region does not contain an obvious ADE gene that could be a secondary gene. Maybe another variant in RGT2?

# ARF2 local effect:
B["chrIV:215523_T/C", "YDL137W"]
# 0.4094716

DESeq_results_pairwise[[3]]["YDL137W",]
# logFC=0.04, p=0.32

#####
# how many of the strong hotspot targets does the DE catch?
thres=0.4
for(hs in 1:3){
    print(c(length(intersect(rownames(DESeq_results_pairwise[[hs]])[DESeq_results_pairwise[[hs]]$pvalue < 0.05], colnames(B)[abs(B[hotspots[hs],]) > thres])),
    length(intersect(rownames(DESeq_results_pairwise[[hs]]), colnames(B)[abs(B[hotspots[hs],]) > thres]))))
}
# FAA4-wt:
# 0.4: 15/15, 9/11, 5/8
# 0.3: 26/28, 12/15, 10/25
# 0.2: 38/70, 19/35, 26/86

# RGT2-wt:
# 0.4: 12/15, 8/11, 6/8
# 0.3: 20/28, 12/15, 11/25
# 0.2: 36/70, 20/35, 27/86


# and what about their direction?
for(hs in 1:3){
    thisDat <- cbind(DESeq_results_pairwise[[hs]][intersect(rownames(DESeq_results_pairwise[[hs]]), colnames(B)),], B[hotspots[hs], intersect(rownames(DESeq_results_pairwise[[hs]]), colnames(B))])
    colnames(thisDat)[length(colnames(thisDat))] <- "BValue"
    thisDat <- thisDat[complete.cases(thisDat),]
    
    print(c(
        nrow(thisDat),
        length(which(thisDat$padj < 0.1)),
        length(which(thisDat$padj < 0.1 & thisDat$BValue != 0)),
        length(which(thisDat$padj < 0.1 & sign(thisDat$BValue) == sign(thisDat$log2FoldChange) & thisDat$BValue != 0))
    ))
}

# FDR10%
# RGT wt
# OAF1 5323   82   39   33
# OLE1 5323   60   30   29
# RGT2 5249   44   38   38

# FAA4 wt
# OAF1 5323   42   30   29
# OLE1 5323   18    9    9
# RGT2 5249   39   25   25


# p< 0.05
#RGT2 wt
# OAF1 5323  316   84   61
# OLE1 5323  205   64   59
# RGT2 5249  182   91   83

# FAA4 wt
# OAF1 5323  143   61   56
# OLE1 5323   93   41   38
# RGT2 5249  212  105   94



# how many significant genes?
sapply(1:3, function(hs){
    length(which(DESeq_results_pairwise[[hs]]$pvalue < 0.05))
})
# RGT2 (correct WT): 183
# OAF1: 147, OLE1: 95

sapply(1:3, function(hs){
    length(which(DESeq_results_pairwise[[hs]]$padj < 0.1))
})
# RGT2: 45
# OAF1: 45, OLE1: 18
