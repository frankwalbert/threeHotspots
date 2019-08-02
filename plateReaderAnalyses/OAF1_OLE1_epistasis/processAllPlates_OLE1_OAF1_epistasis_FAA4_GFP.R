# assume these are all in the same folder
plateReaderExports <- c(
    "190313_BYRM_OLE1OAF1alleles_FAA4-GFP_EXPORT.xlsx",
    "190402_BYRM_OLE1OAF1alleles_FAA4-GFP_EXPORT.xlsx",
    "190430_BYRM_OLE1OAF1alleles_FAA4-GFP_EXPORT.xlsx"
)

# these have to match the order of the plate reader exports above!
plateReaderDescriptionFiles <- c(
    "190313_OLE1OAF1_PlateData.xlsx",
    "190402_OLE1OAF1_PlateData.xlsx",
    "190430_OLE1OAF1_PlateData.xlsx"
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
#save(MaggieOutList, file="R_MaggieOutList_190524.RData")


# put them all into one table
dat <- c()
for (plate in names(MaggieOutList)){
 dat <- rbind(dat, cbind(MaggieOutList[[plate]], plate))
}

# remove rows with blanks
dat <- dat[is.na(dat$Blank),]

colnames(dat)[colnames(dat) == "catDesc"] <- "clone"
dat <- dat[dat$clone != "____",] # remove rows with no data; i.e. wells that were empty
dat$Desc_3 <- factor(dat$Desc_3, levels=c("wt", "OLE1(FAR_RM)", "OAF1(L63S)", "OLE1(FAR_RM)OAF1(L63S)", "OLE1(FAR_BY)", "OAF1(S63L)", "OLE1(FAR_BY)OAF1(S63L)"))


dat$Desc_conc <- factor(paste(dat$Desc_1, dat$Desc_2, dat$Desc_3))
dat$Desc_conc <- factor(dat$Desc_conc, levels=c(
    "BY noGFP wt", "BY FAA4-GFP wt", "BY FAA4-GFP OLE1(FAR_RM)", "BY FAA4-GFP OAF1(L63S)", "BY FAA4-GFP OLE1(FAR_RM)OAF1(L63S)",
    "RM noGFP wt", "RM FAA4-GFP wt", "RM FAA4-GFP OLE1(FAR_BY)", "RM FAA4-GFP OAF1(S63L)", "RM FAA4-GFP OLE1(FAR_BY)OAF1(S63L)", "RM noGFP OLE1(FAR_BY)OAF1(S63L)"
))

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
#resultsAcrossPlates_FAA4 <- statsAcrossPlates(dat[dat$Desc_2 %in% c("FAA4-GFP"),])
resultsAcrossPlates_FAA4 <- statsAcrossPlates(dat[dat$Desc_2 %in% c("FAA4-GFP"),], compareThese="all")
#save(resultsAcrossPlates_FAA4, file="R_resultsAcrossPlates_epistasis_1900524_log.RData")


# one plot per plate
plotPlatesSeparately(dat[dat$Desc_2 %in% c("noGFP", "FAA4-GFP"),])

# remove genotypes that won't be used in paper
dat <- dat[dat$Desc_conc %in% c(
    "BY FAA4-GFP wt", "BY FAA4-GFP OAF1(L63S)", "BY FAA4-GFP OLE1(FAR_RM)", "BY FAA4-GFP OLE1(FAR_RM)OAF1(L63S)",
    "RM FAA4-GFP wt", "RM FAA4-GFP OAF1(S63L)", "RM FAA4-GFP OLE1(FAR_BY)", "RM FAA4-GFP OLE1(FAR_BY)OAF1(S63L)"),]

# FOR POSTER & PAPER: Remove RM completey
#dat <- dat[dat$Desc_conc %in% c("BY FAA4-GFP wt", "BY FAA4-GFP OAF1(L63S)", "BY FAA4-GFP OLE1(FAR_RM)", "BY FAA4-GFP OLE1(FAR_RM)OAF1(L63S)"),]


# plot plates together
plotPlatesTogether(dat[dat$Desc_2 %in% c("FAA4-GFP"),],
    redAlleles=c("BY FAA4-GFP OLE1(FAR_RM)", "BY FAA4-GFP OLE1(FAR_RM)OAF1(L63S)"),
    blueAlleles=c("BY FAA4-GFP wt", "BY FAA4-GFP OAF1(L63S)"))


# epistasis
# here, we don't care about the background, we want to test between alleles
# to still use the same function, we need to hack the data a bit
# make one of the alleles part of the background
datForEpistasis <- dat
datForEpistasis$Desc_1[datForEpistasis$Desc_3 %in% c("OAF1(L63S)", "OLE1(FAR_RM)OAF1(L63S)")] <- "BY_OAF1(L63S)"
datForEpistasis$Desc_1[datForEpistasis$Desc_3 == c("OAF1(S63L)", "OLE1(FAR_BY)OAF1(S63L)")] <- "RM_OAF1(S63L)"

# do in BY, for which we have replicate clones:
epistasisAcrossPlates(datForEpistasis[datForEpistasis$Desc_2 %in% c("FAA4-GFP") & datForEpistasis$Desc_1 %in% c("BY", "BY_OAF1(L63S)"),], allele1Genos=c("BY FAA4-GFP wt", "BY FAA4-GFP OAF1(L63S)"), allele2Genos=c("BY FAA4-GFP OLE1(FAR_RM)", "BY FAA4-GFP OLE1(FAR_RM)OAF1(L63S)"), background1="BY", background2="BY_OAF1(L63S)")
# results as of May 2019, three plates of the newly constructed double-edit strain:
#Ratios_Midlog Ratios_UserDef     Ratios_Sat    Ratios_TMid     growthRate       capacity
#0.030739006    0.008935163    0.006215084    0.002177304    0.093801988    0.414319585
