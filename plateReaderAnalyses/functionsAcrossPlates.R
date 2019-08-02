library(ggplot2)
library(gridExtra)
library(lme4)
library(dplyr)
library(lmerTest) # extends lmer to report p-values for fixed effects in mixed models; used in epistasis tests


############# stats ###############
# run separately for each background (BY or RM) and tag (none, each given gene)
# control for replicate ID and plate

# when this is run in "all" mode, this will really do all vs all, in both directions
# this is wasteful (every test is computed twice), but nice for looking up individual results from the resulting object without having to guess which genotype may have come first

# IMPORTANT NOTE: when input values have been log-transformed, the log2(ratio) the function reports as effect size is not meaningful. Instead, in this case, substract the two trait values from each other to get the correct log2(ratio)

statsAcrossPlates <- function(dat, compareThese="wt"){
    # I want to make sure there is no way that the clone IDs can be confused across genotypes
    # i.e., "1" is a different "1" for BY vs RM, the model must not aggregate those by accident
    # therefore, use the "catDesc" column to batch by replicate (note that we renamed it to "clone" above)

    if (!compareThese %in% c("wt", "all")){print("invalid value for compareThese. Needs to be wt or all"); exit()}

    if(compareThese == "wt"){pdfFileName <- "allComparisons_againstWT.pdf"}
    if(compareThese == "all"){pdfFileName <- "allComparisons_allPairwise.pdf"}
    
    pdf(pdfFileName, width=12, height=10)
    par(mfrow=c(2, 3))
    retAll <- lapply(sort(unique(as.character(dat$Desc_2)[as.character(dat$Desc_2) != "noGFP"]))[1], function(taggedGene){ # function will process just one tagged gene
        retBG <- lapply(sort(unique(as.character(dat$Desc_1))), function(background){ # do separately for BY and RM backgrounds
            if(compareThese == "wt"){
                groups1 <- "wt"
                # all tagged genotypes that are not the wildtype, line is so ugly to maintain level ordering
                groups2 <- levels(dat$Desc_3)[levels(dat$Desc_3) %in% unique(as.character(dat$Desc_3[dat$Desc_1 == background & dat$Desc_2 == taggedGene & dat$Desc_2 != "noGFP" & dat$Desc_3 != "wt"]))]
            }
            
            if(compareThese == "all"){
                groups1 <- levels(dat$Desc_3)[levels(dat$Desc_3) %in% unique(as.character(dat$Desc_3[dat$Desc_1 == background & dat$Desc_2 == taggedGene & dat$Desc_2 != "noGFP"]))]
                groups2 <- groups1
            }
            
            retGroup1 <- lapply(groups1, function(group1) {
                retGroup2 <- lapply(groups2[groups2 != group1], function(group2){
                    thisDat <- dat[dat$Desc_1 == background & dat$Desc_3 %in% c(group1, group2) & dat$Desc_2 == taggedGene, c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity", "Desc_3", "clone", "plate", "Desc_conc")]
                    
                    # one test per phenotype
                    ggPlots <- list()
                    retPerComparison <- sapply(c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity"), function(trait){
                        print(c(background, group1, group2, trait))
                        subDat = thisDat[, c(trait, "Desc_3", "clone", "plate", "Desc_conc")]
                        colnames(subDat) <- c("thisTrait", "genotype", "clone", "plate", "Desc_conc")
                        subDat$plate <- as.factor(subDat$plate)
                        subDat$genotype <- factor(subDat$genotype)
                        print(levels(subDat$genotype))
                        subDat$genotype <- relevel(subDat$genotype, ref=group1)
                        print(levels(subDat$genotype))
                        print(trait)
                        
                        # remove rows with no trait data (e.g. growth rates can fail)
                        subDat <- subDat[!is.na(subDat$thisTrait),]
                        
                        # filter out plates that only carry one of the given genotypes
                        genotypesPerPlate <- sapply(as.character(unique(subDat$plate)), function(x){length(unique(subDat$genotype[subDat$plate == x]))})
                        subDat <- subDat[genotypesPerPlate[as.character(subDat$plate)] > 1,]
                        if(nrow(subDat) == 0){return(NA)}
                        
                        if(length(unique(subDat$plate)) > 1 & length(unique(subDat$clone)) < length(subDat$clone)){
                            h0 = lmer(thisTrait ~ (1|plate) + (1|clone), data = subDat, REML=FALSE)
                            h1 = lmer(thisTrait ~ genotype + (1|plate) + (1|clone), data = subDat, REML=FALSE)
                        }
                        if(length(unique(subDat$plate)) > 1 & length(unique(subDat$clone)) == length(subDat$clone)){
                            h0 = lmer(thisTrait ~ (1|plate), data = subDat, REML=FALSE)
                            h1 = lmer(thisTrait ~ genotype + (1|plate), data = subDat, REML=FALSE)
                        }
                        if(length(unique(subDat$plate)) == 1 & length(unique(subDat$clone)) < nrow(subDat)){
                            h0 = lmer(thisTrait ~ (1|clone), data = subDat, REML=FALSE)
                            h1 = lmer(thisTrait ~ genotype + (1|clone), data = subDat, REML=FALSE)
                        }
                        if(length(unique(subDat$plate)) == 1 & length(unique(subDat$clone)) == nrow(subDat)){
                            h0 = lm(thisTrait ~ 1, data = subDat, REML=FALSE)
                            h1 = lm(thisTrait ~ genotype, data = subDat, REML=FALSE)
                        }
                        
                        # make plots along the way; one for each comparison so we can see the adjusted data
                        # NOTE we do NOT actually "correct" the data for plate effects (plates get incroporated in the model when needed) - so let's plot it raw
                        subDatPlateCorrected <- subDat
                        #subDatPlateCorrected$thisTrait <- residuals(lmer(thisTrait ~ (1|plate)))
                        datMeans <- subDatPlateCorrected %>%
                        group_by(Desc_conc, clone) %>% summarise(
                            meanTrait = mean(thisTrait, na.rm=TRUE),
                        )
                        # or, a non-collapsed version, keeping technical replicates separate
                        # this is for convencience in plotting
                        datMeansNonCollapsed <- subDatPlateCorrected %>% mutate(
                            meanTrait = thisTrait
                        )
                        pValFormatted <- NA
                        if(!is.na(anova(h1, h0)$P[2])){
                            if(anova(h1, h0)$P[2] < 0.001){pValFormatted <- formatC(anova(h1, h0)$P[2], format = "e", digits = 1)}
                            else{pValFormatted <- round(anova(h1, h0)$P[2], 3)}
                        }
                        
                        yLabs <- c("Minmax ratio", "OD=0.25 ratio", "Saturation ratio", "Inflection ratio", "Growth rate", "Capacity")
                        names(yLabs) <- c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")
                        p <- ggplot(datMeans, aes(Desc_conc, meanTrait)) +
                        geom_boxplot() + geom_point(data = datMeansNonCollapsed, aes(colour=clone, shape=plate, group=clone), position=position_dodge(width=0.6)) + geom_line(data = datMeansNonCollapsed, aes(colour=clone), position=position_dodge(width=0.6)) +
                        theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                        labs(x = "") + labs(y = yLabs[trait]) + theme(legend.position="none") +
                        annotate("text", x=-Inf, y=Inf, hjust = -0.2, vjust = 2, label=paste("p-value: ", pValFormatted, sep=""))
                        # push the plot object to ggPlots
                        ggPlots[[trait]] <<- p
                        
                        fixedEffects <- c(NA, NA)
                        if(class(h1) != "lm"){fixedEffects <- fixef(h1)}
                        if(class(h1) == "lm"){fixedEffects <- coef(h1)}
                        
                        retInner <- c(anova(h1, h0)$P[2], fixedEffects[1], sum(fixedEffects), log2(sum(fixedEffects) / fixedEffects[1]))
                        names(retInner) <- c("pValue", "baselineTraitValue", "alternativeTraitValue", "log2FC_AltVsBase")
                        return(retInner)
                    })
                    if(length(ggPlots) == 6){
                        grid.arrange(ggPlots[[1]], ggPlots[[2]], ggPlots[[3]], ggPlots[[4]], ggPlots[[5]], ggPlots[[6]], nrow=2)
                    }
                    
                    return(retPerComparison)
                })
                names(retGroup2) <- groups2[groups2 != group1]
                retGroup2
            })
            names(retGroup1) <- groups1
            retGroup1
            })
        names(retBG) <- sort(unique(as.character(dat$Desc_1)))
        retBG
    })
    dev.off()
    retAll
}

    ########## plots ################

    # plot each plate separately
plotPlatesSeparately <- function(dat){
    pdf("allPlates_separate.pdf", width=18, height=12)
    #par(mfrow=c(1, 3))
    for(thisPlate in unique(dat$plate)){
        datMeans <- dat %>% filter(plate == thisPlate) %>%
            group_by(Desc_conc, clone) %>% summarise(
                meanMid = mean(Ratios_Midlog, na.rm=TRUE),
                meanUser = mean(Ratios_UserDef, na.rm=TRUE),
                meanSat = mean(Ratios_Sat, na.rm=TRUE),
                meanTMid = mean(Ratios_TMid, na.rm=TRUE),
                meanR = mean(growthRate, na.rm=TRUE),
                meanCap = mean(capacity, na.rm=TRUE)
        )
        # or, a non-collapsed version, keeping technical replicates separate
        # this is for convencience in plotting
        datMeansNonCollapsed <- dat %>% filter(plate == thisPlate) %>% mutate(
            meanMid = Ratios_Midlog,
            meanUser = Ratios_UserDef,
            meanSat = Ratios_Sat,
            meanTMid = Ratios_TMid,
            meanR = growthRate,
            meanCap = capacity
        )

        xLabs <- c("Minmax ratio", "OD=0.25 ratio", "Saturation ratio", "Inflection ratio", "Growth rate", "Capacity")
        colID <- c("meanMid", "meanUser", "meanSat", "meanTMid", "meanR", "meanCap")

        ggPlots <- lapply(1:length(colID), function(i){
            p <- ggplot(datMeans, aes(Desc_conc, eval(parse(text=colID[i])))) +
            geom_boxplot() + geom_point(data = datMeansNonCollapsed, aes(colour=clone), position=position_dodge(width=0.6)) + geom_line(data = datMeansNonCollapsed, aes(colour=clone), position=position_dodge(width=0.6)) +
                theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                labs(x = "") + labs(y = xLabs[i]) + theme(legend.position="none") + ggtitle(thisPlate)
        })
        grid.arrange(ggPlots[[1]], ggPlots[[2]], ggPlots[[3]], ggPlots[[4]], ggPlots[[5]], ggPlots[[6]], nrow=2)
    }
    dev.off()
}

# DO NOT USE:
    #datPlateCorrected$Ratios_Midlog[which(!is.na(datPlateCorrected$Ratios_Midlog))] <- residuals(lmer(Ratios_Midlog ~ (1|plate), data=dat))
    #datPlateCorrected$Ratios_UserDef[which(!is.na(datPlateCorrected$Ratios_UserDef))] <- residuals(lmer(Ratios_UserDef ~ (1|plate), data=dat))
    #datPlateCorrected$Ratios_Sat[which(!is.na(datPlateCorrected$Ratios_Sat))] <- residuals(lmer(Ratios_Sat ~ (1|plate), data=dat))

# the input variable name is called 'datPlateCorrected' for historical reasons. The data are not, in fact, "corrected"
# the red and blue alleles allow painting the outsides of the boxplots
# e.g. according to alleles we know to be causal after mapping - we decided this leads the eye too much for presentation in the paper, but is still useful to have as an option


plotPlatesTogether <- function(datPlateCorrected, redAlleles=c(), blueAlleles=c()){
    pdf("allPlates_together.pdf", width=14, height=8, useDingbats=FALSE)
    #par(mfrow=c(1, 3))
        alleleColors <- rep("#000000", nrow(datPlateCorrected))
        if(length(redAlleles) > 0){
            alleleColors[datPlateCorrected$Desc_conc %in% redAlleles] <- "#FF0000"
        }
        if(length(blueAlleles) > 0){
            alleleColors[datPlateCorrected$Desc_conc %in% blueAlleles] <- "#0000FF"
        }
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
    #grid.arrange(ggPlots[[1]], ggPlots[[2]], ggPlots[[3]], nrow=1)
    dev.off()
}



#############
# epistasis
# we'll do this in a focused manner
# need a function that takes a vector of background1/2 and allele1/2 genotypes

epistasisAcrossPlates <- function(dat, background1="BY", background2="RM", allele1Genos, allele2Genos){
        thisDat <- dat[dat$Desc_1 %in% c(background1, background2) & dat$Desc_conc %in% c(allele1Genos, allele2Genos), c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity", "Desc_3", "clone", "plate", "Desc_conc", "Desc_1")]
        
        # filter out plates that do not carry all of the four genotypes
        genotypesPerPlate <- sapply(as.character(unique(thisDat$plate)), function(x){length(unique(thisDat$Desc_conc[thisDat$plate == x]))})
        thisDat <- thisDat[genotypesPerPlate[as.character(thisDat$plate)] == 4,]
        if(nrow(thisDat) == 0){return(NA)}
        
        retPerComparison <- sapply(c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity"), function(trait){
            subDat = thisDat[, c(trait, "Desc_1", "Desc_3", "clone", "plate", "Desc_conc")]
            colnames(subDat) <- c("thisTrait", "background", "genotype", "clone", "plate", "Desc_conc")
            subDat$plate <- as.factor(subDat$plate)
            subDat$allele <- rep("1", nrow(subDat))
            subDat$allele[subDat$Desc_conc %in% allele2Genos] <- "2"
            
            if(length(unique(subDat$plate)) > 1 & length(unique(subDat$clone)) < nrow(subDat)){
                h0 = lmer(thisTrait ~ background + allele + (1|plate) + (1|clone), data = subDat, REML=FALSE)
                h1 = lmer(thisTrait ~ background + allele + background:allele + (1|plate) + (1|clone), data = subDat, REML=FALSE)
            }
            if(length(unique(subDat$plate)) > 1 & length(unique(subDat$clone)) == nrow(subDat)){
                h0 = lmer(thisTrait ~ background + allele + (1|plate), data = subDat, REML=FALSE)
                h1 = lmer(thisTrait ~ background + allele + background:allele + (1|plate), data = subDat, REML=FALSE)
            }
            if(length(unique(subDat$plate)) == 1 & length(unique(subDat$clone)) < nrow(subDat)){
                h0 = lmer(thisTrait ~ background + allele + (1|clone), data = subDat, REML=FALSE)
                h1 = lmer(thisTrait ~ background + allele + background:allele + (1|clone), data = subDat, REML=FALSE)
            }
            if(length(unique(subDat$plate)) == 1 & length(unique(subDat$clone)) == nrow(subDat)){
                h0 = lm(thisTrait ~ background + allele, data = subDat, REML=FALSE)
                h1 = lm(thisTrait ~ background + allele + background:allele, data = subDat, REML=FALSE)
            }
            #print(anova(h1, h0))
            
            #retInner <- c(anova(h1, h0)$P[2], fixedEffects[1], sum(fixedEffects), log2(sum(fixedEffects) / fixedEffects[1]))
            #names(retInner) <- c("pValue", "baselineTraitValue", "alternativeTraitValue", "log2FC_AltVsBase")
            # default is to report just the interaction p-value
            retInner <- c(anova(h1, h0)$P[2])
            # can also print the full ANOVA:
            print(summary((h1)))
            #names(retInner) <- c("pValue", "baselineTraitValue", "alternativeTraitValue", "log2FC_AltVsBase")
            return(retInner)
        })
        retPerComparison
}



# take a pariwise result list and make a table that can be printed for paper supplement
# fill the upper triangle with pValues, the lower triangle with fold changes
tableFromPairwise <- function(pairwiseResults){
    ret2 <- lapply(names(pairwiseResults), function(thisStrain){
        thisRet <- pairwiseResults[[thisStrain]]
        theseNames <- names(thisRet)
        ret1 <- sapply(1:length(theseNames), function(thisGenotype){
            #sapply(names(thisRet[[thisGenotype]]), function(thatGenotype){
            sapply(1:length(theseNames), function(thatGenotype){
                ret <- NA
                if(class(thisRet[[theseNames[thisGenotype]]][[theseNames[thatGenotype]]]) == "matrix"){
                    if(thisGenotype > thatGenotype){
                        ret <- thisRet[[theseNames[thisGenotype]]][[theseNames[thatGenotype]]]["pValue","Ratios_TMid"]
                    }
                    if(thisGenotype < thatGenotype){
                        ret <- thisRet[[theseNames[thisGenotype]]][[theseNames[thatGenotype]]]["log2FC_AltVsBase", "Ratios_TMid"]
                    }
                }
                ret
            })
        })
        colnames(ret1) <- theseNames
        rownames(ret1) <- theseNames
        ret1
    })
    names(ret2) <- names(pairwiseResults)
    ret2
}



