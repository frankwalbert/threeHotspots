library(stringr)
library(ggplot2)

# assume these are all in the same folder
plateReaderExports <- c(
    "171117_RGT2_sugars_EXPORT.xlsx",
    "171120_RGT2_sugars_EXPORT.xlsx",
    "171129_RGT2_glucose_EXPORT.xlsx",
    "190129_HXT-GFP_RGT2glucoseEXPORT.xlsx"
)

# these have to match the order of the plate reader exports above!
plateReaderDescriptionFiles <- c(
    "171117_RGT2_sugars_PlateData.xlsx",
    "171120_RGT2_sugars_PlateData.xlsx",
    "171129_RGT2_glucose_PlateData.xlsx",
    "190129_RGT2_glucose_PlateData.xlsx"
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
#save(MaggieOutList, file="R_MaggieOutList_190131.RData")



# put them all into one table for stats
dat <- c()
for (plate in names(MaggieOutList)){
 dat <- rbind(dat, cbind(MaggieOutList[[plate]], plate))
}

# this is unique to this experiment:
# we need it because on the plate with 0.5 glucose, all the other numbers x become "x.0"
dat$catDesc <- str_replace(dat$catDesc, "\\.0", "")

# remove rows with blanks
dat <- dat[is.na(dat$Blank),]

colnames(dat)[colnames(dat) == "catDesc"] <- "clone"
dat <- dat[dat$clone != "____",] # remove rows with no data; i.e. wells that were empty

# remove galactose and sucrose
dat <- dat[!dat$Desc_5 %in% c("gal2", "suc2", "suc4"),]
dat$Desc_5 <- as.numeric(dat$Desc_5)

# remove a strain (YFA0275) that on the 190129 plate had very high HXT1 expression at higher glucose
# remove it by the wells it's in
dat <- dat[!dat$Wells %in% c("B1", "B3", "B5", "B7", "B9", "B11"),]


dat$Desc_3 <- factor(dat$Desc_3, levels=c("wt", "RGT2(V539I)", "RGT2(I539V)"))

dat$Desc_conc <- factor(paste(dat$Desc_1, dat$Desc_2, dat$Desc_3, dat$Desc_5))
dat$Desc_conc <- factor(dat$Desc_conc, levels=c(
"BY noGFP wt 0.5", "BY HXT1-GFP wt 0.5", "BY HXT1-GFP RGT2(V539I) 0.5", "RM noGFP wt 0.5", "RM HXT1-GFP RGT2(I539V) 0.5", "RM HXT1-GFP wt 0.5",
"BY noGFP wt 1", "BY HXT1-GFP wt 1", "BY HXT1-GFP RGT2(V539I) 1", "RM noGFP wt 1", "RM HXT1-GFP RGT2(I539V) 1", "RM HXT1-GFP wt 1",
"BY noGFP wt 2", "BY HXT1-GFP wt 2", "BY HXT1-GFP RGT2(V539I) 2", "RM noGFP wt 2", "RM HXT1-GFP RGT2(I539V) 2", "RM HXT1-GFP wt 2",
"BY noGFP wt 4", "BY HXT1-GFP wt 4", "BY HXT1-GFP RGT2(V539I) 4", "RM noGFP wt 4", "RM HXT1-GFP RGT2(I539V) 4", "RM HXT1-GFP wt 4",
"BY noGFP wt 8", "BY HXT1-GFP wt 8", "BY HXT1-GFP RGT2(V539I) 8", "RM noGFP wt 8", "RM HXT1-GFP RGT2(I539V) 8", "RM HXT1-GFP wt 8",
"BY noGFP wt 12", "BY HXT1-GFP wt 12", "BY HXT1-GFP RGT2(V539I) 12", "RM noGFP wt 12", "RM HXT1-GFP RGT2(I539V) 12", "RM HXT1-GFP wt 12",
"BY noGFP wt 16", "BY HXT1-GFP wt 16", "BY HXT1-GFP RGT2(V539I) 16", "RM noGFP wt 16", "RM HXT1-GFP RGT2(I539V) 16", "RM HXT1-GFP wt 16"))

# log it all
dat[,c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- log2(dat[,c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")])

# remove background fluorescence
grandMeansBYnoGFP <- apply(dat[dat$Desc_1 == "BY" & dat$Desc_2 == "noGFP" & dat$Desc_3 == "wt", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")], 2, function(x){mean(x, na.rm=TRUE)})
grandMeansRMnoGFP <- apply(dat[dat$Desc_1 == "RM" & dat$Desc_2 == "noGFP" & dat$Desc_3 == "wt", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")], 2, function(x){mean(x, na.rm=TRUE)})
# subtract them (have confirmed that the t(t()) structure works
dat[dat$Desc_1 == "BY", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- t(t(dat[dat$Desc_1 == "BY", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")]) - grandMeansBYnoGFP)
dat[dat$Desc_1 == "RM", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- t(t(dat[dat$Desc_1 == "RM", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")]) - grandMeansRMnoGFP)

# one plot per plate
plotPlatesSeparately(dat)

# plot plates together
# for this, since we already know the causal allele, we can color the plot a bit; see modifications in the source function

# can throw out the noGFPs for higher doses
# throw out all the untagged strains
dat <- dat[!dat$Desc_conc %in% c("BY noGFP wt 0.5", "BY noGFP wt 1", "BY noGFP wt 2", "BY noGFP wt 4", "BY noGFP wt 8", "BY noGFP wt 12", "BY noGFP wt 16", "RM noGFP wt 0.5", "RM noGFP wt 1", "RM noGFP wt 2", "RM noGFP wt 4", "RM noGFP wt 8", "RM noGFP wt 12", "RM noGFP wt 16"),]


#plotPlatesTogether(dat)
plotPlatesTogether(dat,
    redAlleles = dat$Desc_conc[which((str_detect(dat$Desc_conc, "BY") & str_detect(dat$Desc_conc, "RGT2")) | (str_detect(dat$Desc_conc, "RM") & str_detect(dat$Desc_conc, "HXT1-GFP wt")))],
    blueAlleles = dat$Desc_conc[which((str_detect(dat$Desc_conc, "RM") & str_detect(dat$Desc_conc, "RGT2")) | (str_detect(dat$Desc_conc, "BY") & str_detect(dat$Desc_conc, "HXT1-GFP wt")))]
)

# throw out 16% & 0.5% (growth was very poor here; interestingly at 0.5% the effect is dropping again)
plotPlatesTogether(dat[dat$Desc_5 < 16 & dat$Desc_5 > 0.5,],
    redAlleles = dat$Desc_conc[which((str_detect(dat$Desc_conc, "BY") & str_detect(dat$Desc_conc, "RGT2")) | (str_detect(dat$Desc_conc, "RM") & str_detect(dat$Desc_conc, "HXT1-GFP wt")))],
    blueAlleles = dat$Desc_conc[which((str_detect(dat$Desc_conc, "RM") & str_detect(dat$Desc_conc, "RGT2")) | (str_detect(dat$Desc_conc, "BY") & str_detect(dat$Desc_conc, "HXT1-GFP wt")))]
)


# looks like RM grows to much higher capacity than BY, and a bit faster too
# both strains grow worse at higher glucose

# model this with a model that controls for plate & clone as random effects
# background, allele (BY vs RM), effect of glucose level (maybe model as log2?) as fixed effects


thisDat <- dat[dat$Desc_2 == "HXT1-GFP", c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity", "Desc_3", "clone", "plate", "Desc_conc", "Desc_1", "Desc_5")]
colnames(thisDat) <- c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity", "genotype", "clone", "plate", "Desc_conc", "background", "glucose")
thisDat$plate <- as.factor(thisDat$plate)
thisDat$allele <- rep("BY", nrow(thisDat))
thisDat$allele[thisDat$background == "RM" & thisDat$genotype == "wt"] <- "RM"
thisDat$allele[thisDat$background == "BY" & thisDat$genotype == "RGT2(V539I)"] <- "RM"
# log2 the glucose concentration
thisDat$glucose <- log2(thisDat$glucose)

# can remove the 0.5% concentration:
# thisDat <- thisDat[thisDat$glucose != log2(0.5),]

# can drop the 16mM concentration
#thisDat <- thisDat[thisDat$glucose < log2(16),]

# see below for results when these concentrations are included vs excluded (as we report in the paper)

# version with DROPPING the factor rather than ADDING it
retPerComparison <- lapply(c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity"), function(trait){
    subDat = thisDat[, c(trait, "background", "genotype", "allele", "clone", "plate", "Desc_conc", "glucose")]
    colnames(subDat)[1] <- "thisTrait"
    
    # base model with all fixed effects and without interactions
    h1 = lmer(thisTrait ~ background + allele + glucose + (1|plate) + (1|clone), data = subDat, REML=FALSE)
    
    # fixed effect tests
    h0_noBackground = lmer(thisTrait ~ allele + glucose + (1|plate) + (1|clone), data = subDat, REML=FALSE)
    h0_noAllele = lmer(thisTrait ~ background + glucose + (1|plate) + (1|clone), data = subDat, REML=FALSE)
    h0_noGlucose = lmer(thisTrait ~ background + allele + (1|plate) + (1|clone), data = subDat, REML=FALSE)
    
    fixedEffectResults <- list(anova(h1, h0_noBackground), anova(h1, h0_noAllele), anova(h1, h0_noGlucose))
    fixedEffectPValues <- c(anova(h1, h0_noBackground)$P[2], anova(h1, h0_noAllele)$P[2], anova(h1, h0_noGlucose)$P[2])
    names(fixedEffectResults) <- c("background", "allele", "glucose")
    names(fixedEffectPValues) <- c("background", "allele", "glucose")
    
    # interactions
    # all three pairwise interactions, then DROP one at a time
    h1_3pairwise = lmer(thisTrait ~ background + allele + glucose + background:allele + allele:glucose + background:glucose + (1|plate) + (1|clone), data = subDat, REML=FALSE)
    
    # DROPPING the given term
    # interaction of allele effect with background
    h1_BA = lmer(thisTrait ~ background + allele + glucose + allele:glucose + background:glucose + (1|plate) + (1|clone), data = subDat, REML=FALSE)
    # interaction of allele effect with glucose
    h1_AG = lmer(thisTrait ~ background + allele + glucose + background:allele + background:glucose + (1|plate) + (1|clone), data = subDat, REML=FALSE)
    # interaction of background effect with glucose
    h1_BG = lmer(thisTrait ~ background + allele + glucose + background:allele + allele:glucose + (1|plate) + (1|clone), data = subDat, REML=FALSE)
    
    pairwiseResults <- list(anova(h1_3pairwise, h1_BA), anova(h1_3pairwise, h1_AG), anova(h1_3pairwise, h1_BG))
    pairwisePValues <- list(anova(h1_3pairwise, h1_BA)$P[2], anova(h1_3pairwise, h1_AG)$P[2], anova(h1_3pairwise, h1_BG)$P[2])
    names(pairwiseResults) <- c("BA", "AG", "BG")
    names(pairwisePValues) <- c("BA", "AG", "BG")
    
    # three-way: does the allele effect vary with glucose in a strain dependent manner?
    # compare to model with all three pairwise above
    h1_3way = lmer(thisTrait ~ background + allele + glucose + background:allele + allele:glucose + background:glucose + background:glucose:allele + (1|plate) + (1|clone), data = subDat, REML=FALSE)
    threewayResult <- anova(h1_3pairwise, h1_3way)
    threewayPValue <- anova(h1_3pairwise, h1_3way)$P[2]
    
    #return(list(fixedEffectResults, fixedEffectPValues, pairwiseResults, pairwisePValues, threewayResult, threewayPValue))
    return(list(fixedEffectPValues, pairwisePValues, threewayPValue))
})
names(retPerComparison) <- c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")
retPerComparison

# at 0.5% there is a smaller effect, which breaks the linear trend from 1% down
# interesting!

# all concentrations
#$Ratios_TMid[[1]]
#background       allele      glucose
#6.082826e-05 1.680300e-11 4.641803e-58
#$Ratios_TMid[[2]]
#$Ratios_TMid[[2]]$BA
#[1] 0.9642069
#$Ratios_TMid[[2]]$AG
#[1] 2.568653e-07
#$Ratios_TMid[[2]]$BG
#[1] 4.972464e-10
#$Ratios_TMid[[3]]
#[1] 0.3982904

# 0.5% and 16% removed
#$Ratios_TMid[[1]]
#background       allele      glucose
#2.160069e-04 2.780207e-12 6.194472e-38
#$Ratios_TMid[[2]]
#$Ratios_TMid[[2]]$BA
#[1] 0.2531641
#$Ratios_TMid[[2]]$AG
#[1] 6.679547e-13
#$Ratios_TMid[[2]]$BG
#[1] 1.060415e-10
#$Ratios_TMid[[3]]
#[1] 0.02949309

# AG effect is there in both cases




# plotting the effect size
# effects per glucose across plates, extract effect estimate and SD (ACTUALLY, they are standard errors, NOT sd!
effectPerGlucoseAcrossPlates <- lapply(c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity"), function(thisTrait){
    sapply(sort(unique(thisDat$glucose)), function(x){
        subDat <- thisDat[thisDat$glucose == x, c(thisTrait, "clone", "plate", "background", "glucose", "allele", "genotype")]
        colnames(subDat)[1] <- "trait"
        subDat <- subDat[!is.na(subDat$trait),]
        #print(head(subDat))
        BYMeanWT <- mean(subDat[subDat$background == "BY" & subDat$genotype == "wt", "trait"], na.rm=TRUE)
        BYMeanAlt <- mean(subDat[subDat$background == "BY" & subDat$genotype != "wt", "trait"], na.rm=TRUE)
        RMMeanWT <- mean(subDat[subDat$background == "RM" & subDat$genotype == "wt", "trait"], na.rm=TRUE)
        RMMeanAlt <- mean(subDat[subDat$background == "RM" & subDat$genotype != "wt", "trait"], na.rm=TRUE)
        ## NOTE: we invert the RM effects to make them same direction as BY effects and ensure consistent polarization in the plots
        ret <- c(BYMeanWT - BYMeanAlt, RMMeanAlt - RMMeanWT, NA, NA)
        names(ret) <- c("BY_grandMean", "RM_grandMean", "BY_SD", "RM_SD")
        # base model with all fixed effects and without interactions
        # for this experiment, each "clone" (a real clone combined with a glucose level) was only run once per plate, but sometimes there were multiple plates
        # so there is no grouping by clone within plate
        # separately per background:
        for(bg in c("BY", "RM")){
            thisSubDat <- subDat[subDat$background == bg,]
            if (length(unique(thisSubDat$plate)) < length(thisSubDat$plate) & length(unique(thisSubDat$plate)) > 1 & length(unique(thisSubDat$allele)) == 2){
                h1 = lmer(trait ~ allele + (1|plate) + (1|clone), data = thisSubDat, REML=FALSE)
                ret[paste(bg, "SD", sep="_")] <- sqrt(diag(vcov(h1)))[2] # see https://stackoverflow.com/questions/51996252/r-code-for-pulling-out-fixed-effect-standard-errors-in-lme4-package
            }
            if (length(unique(thisSubDat$plate)) == 1 & length(unique(thisSubDat$allele)) == 2){
                h1 = lm(trait ~ allele, data = thisSubDat)
                ret[paste(bg, "SD", sep="_")] <- sqrt(diag(vcov(h1)))[2]
            }
        }
        ret
    })
})
names(effectPerGlucoseAcrossPlates) <- c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")
#effectPerGlucoseAcrossPlates

# SWITCH TO GGPLOT
pdf("effectPerGlucose_withSD_ggplot.pdf", width=9, height=8, useDingbats=FALSE)
par(mfrow=c(3, 2))
ggPlots <- list()
for(i in 1:length(effectPerGlucoseAcrossPlates)){
    ggDat <- data.frame(t(effectPerGlucoseAcrossPlates[[i]]))
    ggDat$glucose <- sort(unique(thisDat$glucose))
    ggDat <- melt(ggDat, id="glucose")
    ggDat$strain <- sapply(as.character(ggDat$variable), function(x){strsplit(x, "_")[[1]][1]})
    ggDat$MeanOrSD <- sapply(as.character(ggDat$variable), function(x){strsplit(x, "_")[[1]][2]})
    ggDat <- ggDat[,names(ggDat) != "variable"]
    ggDat <- dcast(ggDat, glucose + strain ~ MeanOrSD)
    # if wish to exclude the 16%
    #ggDat <- ggDat[ggDat$glucose < 4,]
    # if wish to exclude the 0.5%
    #ggDat <- ggDat[ggDat$glucose >= 0,]
    p <- ggplot(ggDat, aes(x=glucose, y=grandMean, group=strain, color=strain)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = grandMean-SD, ymax = grandMean+SD), width=.2, position=position_dodge(0.05)) +
    labs(title=names(effectPerGlucoseAcrossPlates)[i], x="Glucose concentration (log2(%))", y = "Allele effect (BY - RM allele)") +
    theme_classic() +
    scale_color_manual(values=c('blue','red'))
    ggPlots[[i]] <- p
}
grid.arrange(ggPlots[[1]], ggPlots[[2]], ggPlots[[3]], ggPlots[[4]], ggPlots[[5]], ggPlots[[6]], nrow=3)
dev.off()


