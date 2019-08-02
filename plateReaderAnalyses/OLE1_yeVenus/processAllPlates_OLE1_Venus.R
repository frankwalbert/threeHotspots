# assume these are all in the same folder
plateReaderExports <- c(
    "181012_OLE1pVenus_EXPORT.xlsx",
    "181102_OLE1andSDS23pVenus_EXPORT.xlsx", # single colonies
    "181103_OLE1andSDS23pVenus_pools_EXPORT.xlsx", # pools
    "181206_OAF1_OLE1pVenus_EXPORT.xlsx" # OAF1 alleles in BY background
)

# these have to match the order of the plate reader exports above!
plateReaderDescriptionFiles <- c(
    "181012_OLE1pVenus_PlateData.xlsx",
    "181102_OLE1andSDS23pVenus_PlateData.xlsx",
    "181103_OLE1andSDS23pVenus_pools_PlateData.xlsx",
    "181206_OLE1pVenus_OAF1_PlateData_modified.xlsx"
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
#save(MaggieOutList, file="R_MaggieOutList_181217_OLE1_Venus.RData")


# put them all into one table
dat <- c()
for (plate in names(MaggieOutList)){
 dat <- rbind(dat, cbind(MaggieOutList[[plate]], plate))
}

# remove rows with blanks
dat <- dat[is.na(dat$Blank),]

colnames(dat)[colnames(dat) == "catDesc"] <- "clone"
dat <- dat[dat$clone != "____",] # remove rows with no data; i.e. wells that were empty

# remove the plate 181012, which only has entire BY and RM, no FAR etc
dat <- dat[dat$plate != "181012",]

dat$Desc_3 <- factor(dat$Desc_3, levels=c("noVenus", "OLE1(BY)", "OLE1(RM)", "OLE1(BY_FAR_RM)", "OLE1(RM_FAR_BY)", "OLE1(BY_Astretch_RM)", "OLE1(RM_Astretch_BY)", "SDS23(BY)", "SDS23(RM)", "SDS23(BY_FAR_RM)", "SDS23(RM_FAR_BY)"))


dat$Desc_conc <- factor(paste(dat$Desc_1, dat$Desc_2, dat$Desc_3))
dat$Desc_conc <- factor(dat$Desc_conc, levels=c(
"BY noGFP noVenus", "BY Venus OLE1(BY)", "BY Venus OLE1(RM)", "BY Venus OLE1(BY_FAR_RM)", "BY Venus OLE1(RM_FAR_BY)", "BY Venus OLE1(BY_Astretch_RM)", "BY Venus OLE1(RM_Astretch_BY)", "BY Venus SDS23(BY)", "BY Venus SDS23(RM)", "BY Venus SDS23(BY_FAR_RM)", "BY Venus SDS23(RM_FAR_BY)", "RM noGFP noVenus", "RM Venus OLE1(BY)", "RM Venus OLE1(RM)", "RM Venus OLE1(BY_FAR_RM)", "RM Venus OLE1(RM_FAR_BY)", "RM Venus OLE1(BY_Astretch_RM)", "RM Venus OLE1(RM_Astretch_BY)", "RM Venus SDS23(BY)", "RM Venus SDS23(RM)", "RM Venus SDS23(BY_FAR_RM)", "RM Venus SDS23(RM_FAR_BY)",
    "BY-OAF1(L63S) noGFP noVenus", "BY-OAF1(L63S) Venus OLE1(BY)", "BY-OAF1(L63S) Venus OLE1(BY_FAR_RM)"))


# log-transform the GFP readings
# do these out here, not within the functions
dat[,c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- log(dat[,c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")])

# background-correct
# remove background fluorescence
grandMeansBYnoGFP <- apply(dat[dat$Desc_1 == "BY" & dat$Desc_2 == "noGFP" & dat$Desc_3 == "noVenus", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")], 2, function(x){mean(x, na.rm=TRUE)})
grandMeansRMnoGFP <- apply(dat[dat$Desc_1 == "RM" & dat$Desc_2 == "noGFP" & dat$Desc_3 == "noVenus", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")], 2, function(x){mean(x, na.rm=TRUE)})
#grandMeansBY_OAF_noGFP <- apply(dat[dat$Desc_1 == "BY-OAF1(L63S)" & dat$Desc_2 == "noGFP" & dat$Desc_3 == "noVenus", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")], 2, function(x){mean(x, na.rm=TRUE)})

# subtract them (have confirmed that the t(t()) structure works
dat[dat$Desc_1 == "BY", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- t(t(dat[dat$Desc_1 == "BY", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")]) - grandMeansBYnoGFP)
dat[dat$Desc_1 == "RM", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- t(t(dat[dat$Desc_1 == "RM", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")]) - grandMeansRMnoGFP)
#dat[dat$Desc_1 == "BY-OAF1(L63S)", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- t(t(dat[dat$Desc_1 == "BY-OAF1(L63S)", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")]) - grandMeansBY_OAF_noGFP)



# stats across all plates
# WT comparison does not apply here
resultsAcrossPlates_Venus <- statsAcrossPlates(dat[dat$Desc_2 %in% c("Venus"),], compareThese="all")
#save(resultsAcrossPlates_Venus, file="R_resultsAcrossPlates_OLE1_181217_log.RData")

resultsAcrossPlates_Venus[[1]][["RM"]][["OLE1(DAmP)"]]


# one plot per plate
plotPlatesSeparately(dat[dat$Desc_2 %in% c("noGFP", "Venus"),])

# plot plates together
# focus on the initial experiment:
dat <- dat[dat$Desc_conc %in% c("BY Venus OLE1(BY)", "BY Venus OLE1(RM)", "BY Venus OLE1(BY_FAR_RM)", "BY Venus OLE1(RM_FAR_BY)", "BY Venus SDS23(BY)", "BY Venus SDS23(RM)", "BY Venus SDS23(BY_FAR_RM)", "BY Venus SDS23(RM_FAR_BY)", "RM Venus OLE1(BY)", "RM Venus OLE1(RM)", "RM Venus OLE1(BY_FAR_RM)", "RM Venus OLE1(RM_FAR_BY)", "RM Venus SDS23(BY)", "RM Venus SDS23(RM)", "RM Venus SDS23(BY_FAR_RM)", "RM Venus SDS23(RM_FAR_BY)"),]
dat <- dat[dat$Desc_1 %in% c("BY", "RM"),]

resultsAcrossPlates_Venus <- statsAcrossPlates(dat[dat$Desc_2 %in% c("Venus"),], compareThese="all")
#save(resultsAcrossPlates_Venus, file="R_resultsAcrossPlates_Venus_190617_BYRM.RData")

# write out p-values & log2FoldChanges for supplement table
pValTable <- tableFromPairwise(resultsAcrossPlates_Venus[[1]])
for(i in names(pValTable)){
    write.table(pValTable[[i]], file=paste("pVals_", i,".txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
}



plotPlatesTogether(dat[dat$Desc_2 %in% c("Venus"),],
    redAlleles=c("BY Venus OLE1(RM)", "BY Venus OLE1(BY_FAR_RM)", "BY Venus SDS23(RM)", "BY Venus SDS23(BY_FAR_RM)", "RM Venus OLE1(RM)", "RM Venus OLE1(BY_FAR_RM)", "RM Venus SDS23(RM)", "RM Venus SDS23(BY_FAR_RM)"),
    blueAlleles=c("BY Venus OLE1(BY)", "BY Venus OLE1(RM_FAR_BY)", "BY Venus SDS23(BY)", "BY Venus SDS23(RM_FAR_BY)", "RM Venus OLE1(BY)", "RM Venus OLE1(RM_FAR_BY)", "RM Venus SDS23(BY)", "RM Venus SDS23(RM_FAR_BY)"))


# OR, comparing the OAF1 alleles in the background:
dat <- dat[dat$Desc_conc %in% c("BY Venus OLE1(BY)", "BY Venus OLE1(BY_FAR_RM)", "BY-OAF1(L63S) Venus OLE1(BY)", "BY-OAF1(L63S) Venus OLE1(BY_FAR_RM)"),]
dat <- dat[dat$Desc_1 %in% c("BY", "BY-OAF1(L63S)"),]
# we need to remove the earlier plates from this plot: one BY-BY plate has lots of data but only carries this genotype
# all values from this plate are quite high, making it hard to see what the epistasis test below (which ignores this plate) picks up
dat <- dat[dat$plate == "181206",]
resultsAcrossPlates_Venus_OAF1 <- statsAcrossPlates(dat[dat$Desc_2 %in% c("Venus"),], compareThese="all")
#save(resultsAcrossPlates_Venus_OAF1, file="R_resultsAcrossPlates_Venus_190617_OAF1.RData")

plotPlatesTogether(dat[dat$Desc_2 %in% c("Venus"),],
redAlleles=c("BY Venus OLE1(BY_FAR_RM)", "BY-OAF1(L63S) Venus OLE1(BY_FAR_RM)"),
blueAlleles=c("BY Venus OLE1(BY)", "BY-OAF1(L63S) Venus OLE1(BY)"))


# epistasis
# here, we have several ways to test the causal variant!
# BY vs BY->FAR_RM
epistasisAcrossPlates(dat[dat$Desc_2 %in% c("Venus"),], allele1Genos=c("BY Venus OLE1(BY)", "RM Venus OLE1(BY)"), allele2Genos=c("BY Venus OLE1(BY_FAR_RM)", "RM Venus OLE1(BY_FAR_RM)"))
# combined, log
#Ratios_Midlog Ratios_UserDef     Ratios_Sat    Ratios_TMid     growthRate       capacity
#0.55763149     0.51810114     0.23567174     0.53469396     0.02902314     0.04899912


# RM vs RM->FAR_BY
epistasisAcrossPlates(dat[dat$Desc_2 %in% c("Venus"),], allele1Genos=c("BY Venus OLE1(RM_FAR_BY)", "RM Venus OLE1(RM_FAR_BY)"), allele2Genos=c("BY Venus OLE1(RM)", "RM Venus OLE1(RM)"))
# combined, log
#Ratios_Midlog Ratios_UserDef     Ratios_Sat    Ratios_TMid     growthRate       capacity
#0.2450018      0.2494878      0.1721821      0.2905904      0.6063299      0.3501570


# the two swaps against each other
epistasisAcrossPlates(dat[dat$Desc_2 %in% c("Venus"),], allele1Genos=c("BY Venus OLE1(RM_FAR_BY)", "RM Venus OLE1(RM_FAR_BY)"), allele2Genos=c("BY Venus OLE1(BY_FAR_RM)", "RM Venus OLE1(BY_FAR_RM)"))
# combined, log
#Ratios_Midlog Ratios_UserDef     Ratios_Sat    Ratios_TMid     growthRate       capacity
#0.8376355      0.9233405      0.8636665      0.8870997      0.7457185      0.1347977


# BY vs BY-OAF1
# "combined" (i.e. in this case this works out to just one plate with all genotypes)
epistasisAcrossPlates(dat[dat$Desc_2 %in% c("Venus"),], background1="BY", background2="BY-OAF1(L63S)", allele1Genos=c("BY Venus OLE1(BY)", "BY-OAF1(L63S) Venus OLE1(BY)"), allele2Genos=c("BY Venus OLE1(BY_FAR_RM)", "BY-OAF1(L63S) Venus OLE1(BY_FAR_RM)"))
#Ratios_Midlog Ratios_UserDef     Ratios_Sat    Ratios_TMid     growthRate       capacity
#8.397853e-06   5.531897e-06   7.099996e-04   1.498759e-07   8.511231e-01   9.948203e-01


# does the OAF1 allele affect Venus?
# current code doesn't deal, since the pairwise above runs separately on each background...
# the background correction I use in the plots could be problematic here, since we are directly testing for a difference in expression between backgrounds!
# i.e., the exact thing we're otherwise trying to avoid
thisDat <- dat[dat$plate == "181206" & dat$Desc_3 != "noVenus",]

retPerComparison <- sapply(c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity"), function(trait){
    subDat = thisDat[, c(trait, "Desc_1", "Desc_3", "clone", "plate", "Desc_conc")]
    colnames(subDat) <- c("thisTrait", "background", "genotype", "clone", "plate", "Desc_conc")
    subDat$allele <- subDat$genotype
    
    if(length(unique(subDat$plate)) == 1 & length(unique(subDat$clone)) < nrow(subDat)){
        h0 = lmer(thisTrait ~ allele + (1|clone), data = subDat, REML=FALSE)
        h1 = lmer(thisTrait ~ allele + background + (1|clone), data = subDat, REML=FALSE)
        fixedEffects <- fixef(h1)
        print(fixedEffects)
    }
    #print(summary(h1))
    #print(anova(h1, h0))
    
    retInner <- c(anova(h1, h0)$P[2], sum(fixedEffects[1:2]), sum(fixedEffects), log2(sum(fixedEffects) / sum(fixedEffects[1:2])))
    names(retInner) <- c("pValue", "baselineTraitValue", "alternativeTraitValue", "log2FC_AltVsBase")
    #retInner <- c(anova(h1, h0)$P[2])
    return(retInner)
})
retPerComparison

#                       Ratios_Midlog Ratios_UserDef   Ratios_Sat  Ratios_TMid    growthRate    capacity
#pValue                 8.761040e-15   2.435077e-10 9.854272e-17 4.965808e-06  7.260807e-11  0.78781343
#baselineTraitValue     3.061511e+00   3.124583e+00 3.923070e+00 3.175085e+00 -9.968915e-02 -0.01473894
#alternativeTraitValue  3.459430e+00   3.378628e+00 4.467996e+00 3.340797e+00 -2.778944e-02 -0.01620565
#log2FC_AltVsBase       1.762903e-01   1.127738e-01 1.876451e-01 7.339686e-02 -1.842900e+00  0.13686443
