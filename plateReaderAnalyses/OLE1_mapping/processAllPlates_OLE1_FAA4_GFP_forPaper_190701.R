# assume these are all in the same folder
plateReaderExports <- c(
    "180103_RMandBY_FAA4-GFP_OLE1allelesEXPORT.xlsx",
    "180202_RMandBY_FAA4-GFP_OLE1allelesEXPORT.xlsx",
    "180218_RMandBY_FAA4-GFP_OLE1allelesEXPORT.xlsx",
    "180304_BY_FAA4-GFP_OLE1allelesEXPORT.xlsx",
    "180319_FAA4-GFP_OLE1allelesEXPORT.xlsx",
    "181211_BYRM_FAA4-GFP_OLE1_EXPORT.xlsx"
)

# these have to match the order of the plate reader exports above!
plateReaderDescriptionFiles <- c(
    "180103_OLE1_PlateData_Repaired.xlsx",
    "180202_OLE1_PlateData.xlsx",
    "180218_OLE1_PlateData.xlsx",
    "180304_OLE1_PlateData.xlsx",
    "180319_OLE1_PlateData.xlsx",
    "181211_OLE1_PlateData.xlsx"
)
inputData <- cbind(plateReaderExports, plateReaderDescriptionFiles)

# assuming this is run from the folder that contains the two folders below:
inputPath <- "."
outputPath <- "."

# functions
source("Maggie_MidLogRatio3_asFunction_withGrowthcurver.R")
source("functionsAcrossPlates.R")

MaggieOutList <- lapply(1:nrow(inputData), function(i){
    #MaggieOutList <- lapply(7, function(i){
    print(inputData[i, 1])
    runMaggieOnPlates(
    dataFile = inputData[i, 1],
    descriptionFile = inputData[i, 2],
    outTextFile=paste(strsplit(inputData[i, 1], "_")[[1]][1], "output.txt", sep="_"),
    outPlotFile=paste(strsplit(inputData[i, 1], "_")[[1]][1], "output.pdf", sep="_"),
    userOD=0.25, outputFolder=outputPath, inputFolder=inputPath
    )
})
names(MaggieOutList) <- sapply(plateReaderDescriptionFiles, function(x){strsplit(x, "_")[[1]][1]})
#save(MaggieOutList, file="R_MaggieOutList_190701.RData")


# put them all into one table
dat <- c()
for (plate in names(MaggieOutList)){
 dat <- rbind(dat, cbind(MaggieOutList[[plate]], plate))
}

# remove rows with blanks
dat <- dat[is.na(dat$Blank),]
# remove excluded rows with missing data (one for the case of FAA4 GFP)
dat <- dat[-188,]

# exclude strains based on Sheila's info from 7/1/2019
dat <- dat[!dat$catDesc %in%
    c("BY_FAA4-GFP_OLE1(RM)_1_", "BY_FAA4-GFP_OLE1(RM)_13_", "BY_FAA4-GFP_OLE1(RM)_14_",
    "BY_FAA4-GFP_OLE1(FARRM)_4_",
    "BY_FAA4-GFP_OLE1(RMFARBY)_2_"
),]

colnames(dat)[colnames(dat) == "catDesc"] <- "clone"
dat <- dat[dat$clone != "____",] # remove rows with no data; i.e. wells that were empty
dat$Desc_3 <- factor(dat$Desc_3, levels=c("wt", "sds23::KanMX", "SDS23(RM)", "OLE1(DAmP)", "OLE1(5Hyg3Kan)", "OLE1(RM)", "OLE1(ORFRM)", "OLE1(RM_ORF_SNP)", "OLE1(proRM)", "OLE1(FARRM)", "OLE1(RMFARBY)", "OLE1(BY)", "OLE1(ORFBY)", "OLE1(proBY)", "OLE1(FARBY)"))


dat$Desc_conc <- factor(paste(dat$Desc_1, dat$Desc_2, dat$Desc_3))
dat$Desc_conc <- factor(dat$Desc_conc, levels=c(
"BY noGFP wt", "BY noGFP OLE1(5Hyg3Kan)", "BY FAA4-GFP wt", "BY FAA4-GFP sds23::KanMX", "BY FAA4-GFP SDS23(RM)", "BY FAA4-GFP OLE1(DAmP)", "BY FAA4-GFP OLE1(5Hyg3Kan)", "BY FAA4-GFP OLE1(RM)", "BY FAA4-GFP OLE1(ORFRM)", "BY FAA4-GFP OLE1(RM_ORF_SNP)", "BY FAA4-GFP OLE1(proRM)", "BY FAA4-GFP OLE1(FARRM)", "BY FAA4-GFP OLE1(BY)", "BY FAA4-GFP OLE1(ORFBY)", "BY FAA4-GFP OLE1(RMFARBY)",
"RM noGFP wt", "RM noGFP OLE1(5Hyg3Kan)", "RM FAA4-GFP wt", "RM FAA4-GFP OLE1(DAmP)", "RM FAA4-GFP OLE1(5Hyg3Kan)", "RM FAA4-GFP OLE1(BY)", "RM FAA4-GFP OLE1(proBY)", "RM FAA4-GFP OLE1(FARBY)"))

# log-transform the GFP readings
# do these out here, not within the functions
dat[,c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- log(dat[,c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")])

# remove background fluorescence
grandMeansBYnoGFP <- apply(dat[dat$Desc_conc == "BY noGFP wt", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")], 2, function(x){mean(x, na.rm=TRUE)})
grandMeansRMnoGFP <- apply(dat[dat$Desc_conc == "RM noGFP wt", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")], 2, function(x){mean(x, na.rm=TRUE)})
# subtract them (have confirmed that the t(t()) structure works
dat[dat$Desc_1 == "BY", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- t(t(dat[dat$Desc_1 == "BY", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")]) - grandMeansBYnoGFP)
dat[dat$Desc_1 == "RM", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- t(t(dat[dat$Desc_1 == "RM", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")]) - grandMeansRMnoGFP)


# stats across all plates
resultsAcrossPlates_FAA4 <- statsAcrossPlates(dat[dat$Desc_2 %in% c("FAA4-GFP"),])
#save(resultsAcrossPlates_FAA4, file="R_resultsAcrossPlates_FAA4_190701.RData")

resultsAcrossPlates_FAA4_allPairwise <- statsAcrossPlates(dat[dat$Desc_2 %in% c("FAA4-GFP"),], compareThese="all")
#save(resultsAcrossPlates_FAA4_allPairwise, file="R_resultsAcrossPlates_FAA4_190701_allPairwise.RData")

# write out p-values & log2FoldChanges for supplement table
pValTable <- tableFromPairwise(resultsAcrossPlates_FAA4_allPairwise[[1]])
for(i in names(pValTable)){
    write.table(pValTable[[i]], file=paste("pVals_", i,".txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
}

# one plot per plate
plotPlatesSeparately(dat[dat$Desc_2 %in% c("noGFP", "FAA4-GFP"),])

datForPlot <- dat[dat$Desc_conc %in% c("BY FAA4-GFP wt", "BY FAA4-GFP OLE1(DAmP)", "BY FAA4-GFP OLE1(RM)", "BY FAA4-GFP OLE1(ORFRM)", "BY FAA4-GFP OLE1(proRM)", "BY FAA4-GFP OLE1(FARRM)", "BY FAA4-GFP OLE1(RMFARBY)", "RM FAA4-GFP wt", "RM FAA4-GFP OLE1(FARBY)"),]

# plot plates together
plotPlatesTogether(datForPlot[datForPlot$Desc_2 %in% c("noGFP", "FAA4-GFP"),])


# epistasis
# promoter based
epistasisAcrossPlates(dat[dat$Desc_2 %in% c("FAA4-GFP"),], allele1Genos=c("BY FAA4-GFP wt", "RM FAA4-GFP OLE1(proBY)"), allele2Genos=c("RM FAA4-GFP wt", "BY FAA4-GFP OLE1(proRM)"))
# Ratios_Midlog Ratios_UserDef     Ratios_Sat    Ratios_TMid     growthRate       capacity
# 0.31806374     0.32950494     0.50395668     0.07613819     0.06856859     0.17094579

# FAR variant based
epistasisAcrossPlates(dat[dat$Desc_2 %in% c("FAA4-GFP"),], allele1Genos=c("BY FAA4-GFP wt", "RM FAA4-GFP OLE1(FARBY)"), allele2Genos=c("RM FAA4-GFP wt", "BY FAA4-GFP OLE1(FARRM)"))
#Ratios_Midlog Ratios_UserDef     Ratios_Sat    Ratios_TMid     growthRate       capacity
#0.763779607    0.377068714    0.003029041    0.805084422    0.012771381    0.104409956 
