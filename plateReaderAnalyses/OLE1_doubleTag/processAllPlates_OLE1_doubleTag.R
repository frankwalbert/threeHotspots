# this experiment has both FAA4-GFP and OLE1-mCherry
# but Maggie's code only reads data for one color
# so, we could build a version of the source that reads two tables from one data file
# but that is work
# instead, how about we make two versions of the input data
# one just as is, such that GFP gets read
# one where we copy the mCherry data to the position of the GFP data
# that way, can process the same way
# BE CAREFUL WITH FILENAMES
# make separate plots for the two colors and then combine in Illustrator


# assume these are all in the same folder
plateReaderExports <- c(
    "190227_BY_FAA4-GFP_OLE1-mCherryEXPORT.xlsx",
    "190314_FAA4-GFP_OLE1mCherryEXPORT.xlsx",
    "190227_BY_FAA4-GFP_OLE1-mCherryEXPORT_mCherry.xlsx",
    "190314_FAA4-GFP_OLE1mCherryEXPORT_mCherry.xlsx"
)

# these have to match the order of the plate reader exports above!
plateReaderDescriptionFiles <- c(
    "190227_PlateData.xlsx",
    "190314_PlateData.xlsx",
    "190227_PlateData_mCherry.xlsx",
    "190314_PlateData_mCherry.xlsx"
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
    outTextFile=paste(inputData[i, 1], "output.txt", sep="_"),
    outPlotFile=paste(inputData[i, 1], "output.pdf", sep="_"),
    userOD=0.25, outputFolder=outputPath, inputFolder=inputPath
    )
})
#names(MaggieOutList) <- sapply(plateReaderDescriptionFiles, function(x){strsplit(x, "_")[[1]][1]})
names(MaggieOutList) <- plateReaderDescriptionFiles
#save(MaggieOutList, file="R_MaggieOutList_190616.RData")



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

dat$Desc_3 <- factor(dat$Desc_3, levels=c("wt", "OLE1(RM)", "OLE1(FAR_RM)", "OLE1-mCherry", "OLE1(RM)-mCherry", "OLE1(FAR_RM)-mCherry"))


dat$Desc_conc <- factor(paste(dat$Desc_1, dat$Desc_2, dat$Desc_3))
dat$Desc_conc <- factor(dat$Desc_conc, levels=c(
"BY FAA4-GFP wt", "BY FAA4-GFP OLE1(RM)", "BY FAA4-GFP OLE1(FAR_RM)", "BY FAA4-GFP OLE1-mCherry", "BY FAA4-GFP OLE1(RM)-mCherry", "BY FAA4-GFP OLE1(FAR_RM)-mCherry", "BY OLE1-mCherry OLE1-mCherry", "BY OLE1-mCherry OLE1(RM)-mCherry", "BY OLE1-mCherry OLE1(FAR_RM)-mCherry"
))

# log-transform the GFP readings
# do these out here, not within the functions
dat[,c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- log2(dat[,c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")])


# stats across all plates
resultsAcrossPlatesPairwise_GFP <- statsAcrossPlates(dat[dat$Desc_2 %in% c("FAA4-GFP"),], compareThese="all")
#save(resultsAcrossPlatesPairwise_GFP, file="R_resultsAcrossPlates_190530_log_GFP.RData")
resultsAcrossPlatesPairwise_mCherry <- statsAcrossPlates(dat[dat$Desc_2 %in% c("OLE1-mCherry"),], compareThese="all")
#save(resultsAcrossPlatesPairwise_mCherry, file="R_resultsAcrossPlates_190530_log_mCherry.RData")

# one plot per plate
plotPlatesSeparately(dat)

# plot plates together
# RENAME!
plotPlatesTogether(dat)
plotPlatesTogether(dat[dat$Desc_2 == "FAA4-GFP",])
plotPlatesTogether(dat[dat$Desc_2 == "OLE1-mCherry",])

# for paper main figure: colors separate, just double-tag, just variant
plotPlatesTogether(dat[dat$Desc_2 == "FAA4-GFP" & dat$Desc_conc %in% c("BY FAA4-GFP OLE1-mCherry", "BY FAA4-GFP OLE1(FAR_RM)-mCherry"),])
plotPlatesTogether(dat[dat$Desc_2 == "OLE1-mCherry" & dat$Desc_conc %in% c("BY OLE1-mCherry OLE1-mCherry", "BY OLE1-mCherry OLE1(FAR_RM)-mCherry"),])
