# assume these are all in the same folder
plateReaderExports <- c(
    "170630_HXT1GFPrgt2swapsEXPORT.xlsx",
    "170808_RGT2_EXPORT.xlsx",
    "170902_HXT1GFP_rgt2EXPORT.xlsx",
    "171115_RMBY_HXT1GFPrgt2EXPORT.xlsx"
)

# these have to match the order of the plate reader exports above!
plateReaderDescriptionFiles <- c(
    "rgt2datafileswithindelgenotypes/170630_RGT2_2percent_Glu_PlateData_YFA_INDEL.xlsx",
    "rgt2datafileswithindelgenotypes/170808_RGT2_PlateData_YFA_INDEL.xlsx",
    "rgt2datafileswithindelgenotypes/170902_RGT2_PlateData_YFA_INDEL.xlsx",
    "rgt2datafileswithindelgenotypes/171115_RGT2_PlateData_YFA_INDEL.xlsx"
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
#save(MaggieOutList, file="R_MaggieOutList_190113.RData")



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
dat$Desc_3 <- factor(dat$Desc_3, levels=c("wt", "rgt2::KanMX", "RGT2(RM)", "RGT2(725BY)", "RGT2(639BY)", "RGT2(503RM)", "RGT2(RM596BY)", "RGT2(RM539BY)", "RGT2(639RM)", "RGT2(503BY)", "RGT2(V539I)", "RGT2_5'UTR(RM)", "RGT2(BY)", "RGT2(K596N)", "RGT2(I539V)"))


dat$Desc_conc <- factor(paste(dat$Desc_1, dat$Desc_2, dat$Desc_3))
dat$Desc_conc <- factor(dat$Desc_conc, levels=c(
"BY noGFP wt", "BY HXT1-GFP wt", "BY HXT1-GFP rgt2::KanMX", "BY HXT1-GFP RGT2(RM)", "BY HXT1-GFP RGT2_5'UTR(RM)", "BY HXT1-GFP RGT2(725BY)", "BY HXT1-GFP RGT2(639BY)", "BY HXT1-GFP RGT2(639RM)", "BY HXT1-GFP RGT2(503BY)", "BY HXT1-GFP RGT2(503RM)", "BY HXT1-GFP RGT2(RM596BY)", "BY HXT1-GFP RGT2(V539I)", "BY HXT1-GFP RGT2(RM539BY)",
"RM noGFP wt", "RM HXT1-GFP wt", "RM HXT1-GFP rgt2::KanMX", "RM HXT1-GFP RGT2(BY)", "RM HXT1-GFP RGT2(725BY)", "RM HXT1-GFP RGT2(K596N)", "RM HXT1-GFP RGT2(I539V)"))

# log-transform the GFP readings
# do these out here, not within the functions
dat[,c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- log2(dat[,c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")])

# remove background fluorescence
grandMeansBYnoGFP <- apply(dat[dat$Desc_1 == "BY" & dat$Desc_2 == "noGFP" & dat$Desc_3 == "wt", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")], 2, function(x){mean(x, na.rm=TRUE)})
grandMeansRMnoGFP <- apply(dat[dat$Desc_1 == "RM" & dat$Desc_2 == "noGFP" & dat$Desc_3 == "wt", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")], 2, function(x){mean(x, na.rm=TRUE)})
# subtract them (have confirmed that the t(t()) structure works
dat[dat$Desc_1 == "BY", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- t(t(dat[dat$Desc_1 == "BY", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")]) - grandMeansBYnoGFP)
dat[dat$Desc_1 == "RM", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- t(t(dat[dat$Desc_1 == "RM", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")]) - grandMeansRMnoGFP)


# throw out alleles we don't need to plot:
dat <- dat[dat$Desc_conc %in% c("BY HXT1-GFP wt", "BY HXT1-GFP rgt2::KanMX", "BY HXT1-GFP RGT2(RM)", "BY HXT1-GFP RGT2(639BY)", "BY HXT1-GFP RGT2(503RM)", "BY HXT1-GFP RGT2(RM539BY)", "BY HXT1-GFP RGT2(639RM)", "BY HXT1-GFP RGT2(503BY)", "BY HXT1-GFP RGT2(V539I)", "BY HXT1-GFP RGT2_5'UTR(RM)", "RM HXT1-GFP wt", "RM HXT1-GFP rgt2::KanMX", "RM HXT1-GFP RGT2(BY)", "RM HXT1-GFP RGT2(I539V)"),]
# drop the kanMX deletes?
#dat <- dat[!dat$Desc_conc %in% c("BY HXT1-GFP rgt2::KanMX", "RM HXT1-GFP rgt2::KanMX"),]

# drop one 5'UTR swap with inconsistent indel
dat <- dat[dat$Desc_5 != "YFA0351",]
# drop BY RGT2(RM) strains in which the indel genotype is not confirmed
dat <- dat[!dat$Desc_5 %in% c("YFA0315", "YFA0316", "YFA0317", "YFA0318"),]
# drop additional strains with unknown indel genotype
dat <- dat[!dat$Desc_5 %in% "not filed - remove",]


# stats across all plates
resultsAcrossPlates <- statsAcrossPlates(dat)
#save(resultsAcrossPlates, file="R_resultsAcrossPlates_190704.RData")

resultsAcrossPlates_allPairwise <- statsAcrossPlates(dat[dat$Desc_2 %in% c("HXT1-GFP"),], compareThese="all")
#save(resultsAcrossPlates_allPairwise, file="R_resultsAcrossPlates_allPairwise_190704.RData")

# write out p-values & log2FoldChanges for supplement table
pValTable <- tableFromPairwise(resultsAcrossPlates_allPairwise[[1]])
for(i in names(pValTable)){
    write.table(pValTable[[i]], file=paste("pVals_", i,".txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
}


# one plot per plate
plotPlatesSeparately(dat)

# plot plates together
plotPlatesTogether(dat)

# epistasis
epistasisAcrossPlates(dat[dat$Desc_2 %in% c("HXT1-GFP"),], allele1Genos=c("BY HXT1-GFP wt", "RM HXT1-GFP RGT2(I539V)"), allele2Genos=c("RM HXT1-GFP wt", "BY HXT1-GFP RGT2(V539I)"))
# log
#Ratios_Midlog Ratios_UserDef     Ratios_Sat    Ratios_TMid     growthRate       capacity
#0.096511887    0.006462214    0.111325320    0.099530524    0.038813189    0.020271060


# all plates combined, but with that one indel genotype colored in (in the paper figure, we changed this to dots indicating indel alleles using Illustrator)
# pull out the code from the function and adapt manually for this edge case
pdf("allPlates_together_RGT2_indel.pdf", width=14, height=8, useDingbats=FALSE)
    datPlateCorrected <- dat
    alleleColors <- rep("#000000", nrow(datPlateCorrected))
    alleleColors[datPlateCorrected$"chrIV.213315" == "RM"] <- "#FF0000"
    alleleColors[datPlateCorrected$"chrIV.213315" == "BY"] <- "#0000FF"
    datPlateCorrected$alleleColors <- alleleColors
    datPlateCorrected$alleleColorsFaint <- paste(alleleColors, "44", sep="")
    # in case we want none of these fancy faint colors:
    #datPlateCorrected$alleleColorsFaint <- "black"
    
    datMeans <- datPlateCorrected %>%
    group_by(Desc_conc, clone) %>% summarise(
    meanMid = mean(Ratios_Midlog, na.rm=TRUE),
    meanUser = mean(Ratios_UserDef, na.rm=TRUE),
    meanSat = mean(Ratios_Sat, na.rm=TRUE),
    meanTMid = mean(Ratios_TMid, na.rm=TRUE),
    meanR = mean(growthRate, na.rm=TRUE),
    meanCap = mean(capacity, na.rm=TRUE),
    background = factor(Desc_1[1]),
    alleleColors = factor(alleleColors[1]),
    )
    # or, a non-collapsed version, keeping technical replicates separate
    # this is for convencience in plotting
    datMeansNonCollapsed <- datPlateCorrected %>% mutate(
    meanMid = Ratios_Midlog,
    meanUser = Ratios_UserDef,
    meanSat = Ratios_Sat,
    meanTMid = Ratios_TMid,
    meanR = growthRate,
    meanCap = capacity,
    plate = as.factor(plate)
    )
    
    xLabs <- c("log2(GFP/OD) at half maximum OD point", "log2(GFP/OD) at OD=0.25", "log2(GFP/OD) in saturation", "log2(GFP/OD) at inflection point", "log2(Growth rate)", "log2(Capacity)")
    colID <- c("meanMid", "meanUser", "meanSat", "meanTMid", "meanR", "meanCap")
    
    for(i in 1:length(colID)){
        p <- ggplot(datMeans, aes(Desc_conc, eval(parse(text=colID[i])))) +
        geom_boxplot(aes(fill=background, color="grey"), outlier.color=NA) + # deactivate if want all black & white
        scale_fill_manual(values=c("#0000FF44", "#FF000044")) + # this will sort the backgrounds alphabetically, which works for BY/RM but will become ugly for more/other backgrounds, deactivate if want all black & white
        #geom_boxplot(outlier.color=NA, color="grey") + # activate if want all BW
        scale_colour_identity() +
        #geom_boxplot(outlier.color=NA) + # outlier suppression seems OK since we do show all the measurements anyway
        geom_point(data = datMeansNonCollapsed, aes(shape=plate, group=clone, color=alleleColors), position=position_dodge(width=0.6)) + # , color="#00000044"
        geom_line(data = datMeansNonCollapsed, aes(group=clone, color=alleleColorsFaint), position=position_dodge(width=0.6)) +
        theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "") + labs(y = xLabs[i]) + theme(legend.position="none") #+ #scale_shape_manual(values=0:(nlevels(datMeansNonCollapsed$plate)*2)) +
        print(p)
    }
dev.off()
