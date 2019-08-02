# assume these are all in the same folder
plateReaderExports <- c(
    "171203_RMBY_FAA4GFPoaf1_EXPORT.xlsx",
    "180404_FAA4-GFP_OAF1allelesEXPORT.xlsx"
)

# these have to match the order of the plate reader exports above!
plateReaderDescriptionFiles <- c(
    "171203_OAF1_PlateData_YFA.xlsx",
    "180404_OAF1_PlateData.xlsx"
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
#save(MaggieOutList, file="R_MaggieOutList_180906.RData")



# put them all into one table for stats
dat <- c()
for (plate in names(MaggieOutList)){
 dat <- rbind(dat, cbind(MaggieOutList[[plate]], plate))
}

# remove rows with blanks
dat <- dat[is.na(dat$Blank),]

colnames(dat)[colnames(dat) == "catDesc"] <- "clone"
dat <- dat[dat$clone != "____",] # remove rows with no data; i.e. wells that were empty
#dat <- dat[dat$Desc_2 != "noGFP",] # remove strains without GFP
dat <- dat[dat$Desc_3 != "OAF1(T184I)",] # remove strains without GFP# remove un-needed genotype

dat$Desc_3 <- factor(dat$Desc_3, levels=c("wt", "oaf1::KanMX", "OAF1(RM)", "OAF1(ChA)", "OAF1(ChB)", "OAF1(L63S)", "OAF1(I184T)", "OAF1(BY)", "OAF1(S63L)", "OAF1(ChA1RM)"))


dat$Desc_conc <- factor(paste(dat$Desc_1, dat$Desc_2, dat$Desc_3))
dat$Desc_conc <- factor(dat$Desc_conc, levels=c(
"BY noGFP wt", "BY FAA4-GFP wt", "BY FAA4-GFP oaf1::KanMX", "BY FAA4-GFP OAF1(RM)", "BY FAA4-GFP OAF1(ChA)", "BY FAA4-GFP OAF1(ChB)", "BY FAA4-GFP OAF1(L63S)", "BY FAA4-GFP OAF1(I184T)",
"RM noGFP wt", "RM FAA4-GFP wt", "RM FAA4-GFP oaf1::KanMX", "RM FAA4-GFP OAF1(BY)", "RM FAA4-GFP OAF1(S63L)", "RM FAA4-GFP OAF1(ChA1RM)"))

# log-transform the GFP readings
# do these out here, not within the functions
dat[,c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- log2(dat[,c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")])

# remove background fluorescence
grandMeansBYnoGFP <- apply(dat[dat$Desc_1 == "BY" & dat$Desc_2 == "noGFP" & dat$Desc_3 == "wt", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")], 2, function(x){mean(x, na.rm=TRUE)})
grandMeansRMnoGFP <- apply(dat[dat$Desc_1 == "RM" & dat$Desc_2 == "noGFP" & dat$Desc_3 == "wt", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")], 2, function(x){mean(x, na.rm=TRUE)})
# subtract them (have confirmed that the t(t()) structure works
dat[dat$Desc_1 == "BY", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- t(t(dat[dat$Desc_1 == "BY", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")]) - grandMeansBYnoGFP)
dat[dat$Desc_1 == "RM", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- t(t(dat[dat$Desc_1 == "RM", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")]) - grandMeansRMnoGFP)

dat <- dat[dat$Desc_2 != "noGFP",] # remove strains without GFP

# stats across all plates
resultsAcrossPlates <- statsAcrossPlates(dat)
#save(resultsAcrossPlates, file="R_resultsAcrossPlates_190704.RData")

resultsAcrossPlatesPairwise <- statsAcrossPlates(dat[dat$Desc_2 %in% c("FAA4-GFP"),], compareThese="all")
#save(resultsAcrossPlatesPairwise, file="R_resultsAcrossPlatesPairwise_190704.RData")

# write out p-values & log2FoldChanges for supplement table
pValTable <- tableFromPairwise(resultsAcrossPlatesPairwise[[1]])
for(i in names(pValTable)){
    write.table(pValTable[[i]], file=paste("pVals_", i,".txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
}


# one plot per plate
plotPlatesSeparately(dat)

# plot plates together
plotPlatesTogether(dat, redAlleles = c("BY FAA4-GFP OAF1(RM)", "BY FAA4-GFP OAF1(ChA)", "BY FAA4-GFP OAF1(L63S)", "RM noGFP wt", "RM FAA4-GFP wt", "RM FAA4-GFP OAF1(ChA1RM)", "RM FAA4-GFP OAF1(T184I)"), blueAlleles = c("BY noGFP wt", "BY FAA4-GFP wt", "BY FAA4-GFP OAF1(ChB)", "BY FAA4-GFP OAF1(I184T)", "RM FAA4-GFP OAF1(BY)", "RM FAA4-GFP OAF1(S63L)"))

plotPlatesTogether(dat)

# version with only the two WTs and the causal edit:
plotPlatesTogether(dat[dat$Desc_3 %in% c("wt", "OAF1(L63S)", "OAF1(S63L)"),])

# epistasis
epistasisAcrossPlates(dat[dat$Desc_2 %in% c("FAA4-GFP"),], allele1Genos=c("BY FAA4-GFP wt", "RM FAA4-GFP OAF1(S63L)"), allele2Genos=c("RM FAA4-GFP wt", "BY FAA4-GFP OAF1(L63S)"))
# linear
# Ratios_Midlog Ratios_UserDef     Ratios_Sat    Ratios_TMid     growthRate       capacity
# 0.29111451     0.00363691     0.22589216     0.01410285     0.52393450     0.18427699

# log2
#Ratios_Midlog Ratios_UserDef     Ratios_Sat    Ratios_TMid     growthRate       capacity
#0.7548679      0.1680101      0.7068599      0.4379090      0.5158809      0.1276365
