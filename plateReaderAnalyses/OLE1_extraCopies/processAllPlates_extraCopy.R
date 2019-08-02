# assume these are all in the same folder
plateReaderExports <- c(
    "190502_BY_FAA4-GFP_extracopiesEXPORT.xlsx"
)

# these have to match the order of the plate reader exports above!
plateReaderDescriptionFiles <- c(
    "190502_BY_FAA4-GFP_extracopies_PlateData.xlsx"
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
dat <- dat[dat$Desc_2 != "noGFP",] # remove strains without GFP

dat$Desc_3 <- factor(dat$Desc_3, levels=c("pRS415", "pRS415-OLE1", "pRS415-SDS23", "pRS415-OLE1-SDS23p"))


dat$Desc_conc <- factor(paste(dat$Desc_1, dat$Desc_2, dat$Desc_3))
dat$Desc_conc <- factor(dat$Desc_conc, levels=c(
"BY FAA4-GFP pRS415", "BY FAA4-GFP pRS415-OLE1", "BY FAA4-GFP pRS415-SDS23", "BY FAA4-GFP pRS415-OLE1-SDS23p"
))

# log-transform the GFP readings
# do these out here, not within the functions
dat[,c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- log2(dat[,c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")])


# stats across all plates
resultsAcrossPlatesPairwise <- statsAcrossPlates(dat[dat$Desc_2 %in% c("FAA4-GFP"),], compareThese="all")
#save(resultsAcrossPlatesPairwise, file="R_resultsAcrossPlates_190530_log.RData")

# one plot per plate
plotPlatesSeparately(dat)

# plot plates together
plotPlatesTogether(dat)

