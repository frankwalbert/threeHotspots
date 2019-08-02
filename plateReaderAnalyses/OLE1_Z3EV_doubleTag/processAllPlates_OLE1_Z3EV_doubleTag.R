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
    "190531_Z3EVOLE1-mCherry_FAA4-GFP_EXPORT_GFP.xlsx",
    "190531_Z3EVOLE1-mCherry_FAA4-GFP_EXPORT_mCherry.xlsx"
)

# these have to match the order of the plate reader exports above!
plateReaderDescriptionFiles <- c(
    "190531_Z3EVpOLE1_FAA4-GFP_PlateData_mod_GFP.xlsx",
    "190531_Z3EVpOLE1_FAA4-GFP_PlateData_mod_mCherry.xlsx"
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
#save(MaggieOutList, file="R_MaggieOutList_190611.RData")



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

dat$Desc_3 <- factor(dat$Desc_3, levels=c("wt", "OLE1(FAR_RM)", "Z3EVpOLE1"))


dat$Desc_conc <- factor(paste(dat$Desc_1, dat$Desc_2, dat$Desc_3, dat$Desc_5))
dat$Desc_conc <- factor(dat$Desc_conc, levels=c(
"BY OLE1-mCherry wt 0", "BY OLE1-mCherry wt 1", "BY OLE1-mCherry wt 2", "BY OLE1-mCherry wt 2.5", "BY OLE1-mCherry wt 3", "BY OLE1-mCherry wt 3.5", "BY OLE1-mCherry wt 4", "BY OLE1-mCherry wt 5", "BY OLE1-mCherry wt 6", "BY OLE1-mCherry wt 7", "BY OLE1-mCherry wt 8", "BY OLE1-mCherry wt 9", "BY OLE1-mCherry wt 10", "BY OLE1-mCherry wt 15", "BY OLE1-mCherry wt 30", "BY OLE1-mCherry wt 40",
"BY OLE1-mCherry OLE1(FAR_RM) 0", "BY OLE1-mCherry OLE1(FAR_RM) 1", "BY OLE1-mCherry OLE1(FAR_RM) 2", "BY OLE1-mCherry OLE1(FAR_RM) 2.5", "BY OLE1-mCherry OLE1(FAR_RM) 3", "BY OLE1-mCherry OLE1(FAR_RM) 3.5", "BY OLE1-mCherry OLE1(FAR_RM) 4", "BY OLE1-mCherry OLE1(FAR_RM) 5", "BY OLE1-mCherry OLE1(FAR_RM) 6", "BY OLE1-mCherry OLE1(FAR_RM) 7", "BY OLE1-mCherry OLE1(FAR_RM) 8", "BY OLE1-mCherry OLE1(FAR_RM) 9", "BY OLE1-mCherry OLE1(FAR_RM) 10", "BY OLE1-mCherry OLE1(FAR_RM) 15", "BY OLE1-mCherry OLE1(FAR_RM) 30", "BY OLE1-mCherry OLE1(FAR_RM) 40",
"BY OLE1-mCherry Z3EVpOLE1 0", "BY OLE1-mCherry Z3EVpOLE1 1", "BY OLE1-mCherry Z3EVpOLE1 2", "BY OLE1-mCherry Z3EVpOLE1 2.5", "BY OLE1-mCherry Z3EVpOLE1 3", "BY OLE1-mCherry Z3EVpOLE1 3.5", "BY OLE1-mCherry Z3EVpOLE1 4", "BY OLE1-mCherry Z3EVpOLE1 5", "BY OLE1-mCherry Z3EVpOLE1 6", "BY OLE1-mCherry Z3EVpOLE1 7", "BY OLE1-mCherry Z3EVpOLE1 8", "BY OLE1-mCherry Z3EVpOLE1 9", "BY OLE1-mCherry Z3EVpOLE1 10", "BY OLE1-mCherry Z3EVpOLE1 15", "BY OLE1-mCherry Z3EVpOLE1 30", "BY OLE1-mCherry Z3EVpOLE1 40",
"BY FAA4-GFP wt 0", "BY FAA4-GFP wt 1", "BY FAA4-GFP wt 2", "BY FAA4-GFP wt 2.5", "BY FAA4-GFP wt 3", "BY FAA4-GFP wt 3.5", "BY FAA4-GFP wt 4", "BY FAA4-GFP wt 5", "BY FAA4-GFP wt 6", "BY FAA4-GFP wt 7", "BY FAA4-GFP wt 8", "BY FAA4-GFP wt 9", "BY FAA4-GFP wt 10", "BY FAA4-GFP wt 15", "BY FAA4-GFP wt 30", "BY FAA4-GFP wt 40",
"BY FAA4-GFP OLE1(FAR_RM) 0", "BY FAA4-GFP OLE1(FAR_RM) 1", "BY FAA4-GFP OLE1(FAR_RM) 2", "BY FAA4-GFP OLE1(FAR_RM) 2.5", "BY FAA4-GFP OLE1(FAR_RM) 3", "BY FAA4-GFP OLE1(FAR_RM) 3.5", "BY FAA4-GFP OLE1(FAR_RM) 4", "BY FAA4-GFP OLE1(FAR_RM) 5", "BY FAA4-GFP OLE1(FAR_RM) 6", "BY FAA4-GFP OLE1(FAR_RM) 7", "BY FAA4-GFP OLE1(FAR_RM) 8", "BY FAA4-GFP OLE1(FAR_RM) 9", "BY FAA4-GFP OLE1(FAR_RM) 10", "BY FAA4-GFP OLE1(FAR_RM) 15", "BY FAA4-GFP OLE1(FAR_RM) 30", "BY FAA4-GFP OLE1(FAR_RM) 40",
"BY FAA4-GFP Z3EVpOLE1 0", "BY FAA4-GFP Z3EVpOLE1 1", "BY FAA4-GFP Z3EVpOLE1 2", "BY FAA4-GFP Z3EVpOLE1 2.5", "BY FAA4-GFP Z3EVpOLE1 3", "BY FAA4-GFP Z3EVpOLE1 3.5", "BY FAA4-GFP Z3EVpOLE1 4", "BY FAA4-GFP Z3EVpOLE1 5", "BY FAA4-GFP Z3EVpOLE1 6", "BY FAA4-GFP Z3EVpOLE1 7", "BY FAA4-GFP Z3EVpOLE1 8", "BY FAA4-GFP Z3EVpOLE1 9", "BY FAA4-GFP Z3EVpOLE1 10", "BY FAA4-GFP Z3EVpOLE1 15", "BY FAA4-GFP Z3EVpOLE1 30", "BY FAA4-GFP Z3EVpOLE1 40"
))

# log-transform the GFP readings
# do these out here, not within the functions
dat[,c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")] <- log2(dat[,c("Ratios_Midlog", "Ratios_UserDef", "Ratios_Sat", "Ratios_TMid", "growthRate", "capacity")])

#all are BY, so we don't need background subtraction/correction

# stats across all plates
# p-values don't really make sense for this analysis, certainly not all pairwise

# one plot per plate
plotPlatesSeparately(dat)

# plot plates together
# RENAME!
plotPlatesTogether(dat)
plotPlatesTogether(dat[dat$Desc_2 == "FAA4-GFP",])
plotPlatesTogether(dat[dat$Desc_2 == "OLE1-mCherry",])

# mean per group and the range (there are only two measures per strain & concentration)
datMeans <- dat %>% group_by(Desc_conc) %>% summarise(
    genotype = first(Desc_3),
    estradiol = mean(Desc_5),
    tag = first(Desc_2),
    expression=mean(Ratios_UserDef),
    SD=sd(Ratios_UserDef),
    minimum=min(Ratios_UserDef),
    maximum=max(Ratios_UserDef)
)

datMeansGrowth <- dat[dat$Desc_2 == "FAA4-GFP",] %>% group_by(Desc_conc) %>% summarise(
    genotype = first(Desc_3),
    estradiol = mean(Desc_5),
    growth=mean(growthRate),
    SDGrowth=sd(growthRate),
    minGrowth=min(growthRate),
    maxGrowth=max(growthRate)
)

pdf("estradiol.pdf", width=9, height=5)

print(ggplot(datMeans[datMeans$estradiol >= 4 & datMeans$estradiol < 30,], aes(x=estradiol, y=expression, color=tag)) +
geom_line(aes(linetype=genotype, size=genotype)) +
geom_pointrange(aes(ymin=expression-SD, ymax=expression+SD), position=position_jitter(width=0.01)) +
scale_color_manual(values=c("green", "red")) +
scale_linetype_manual(values=c("dotted", "twodash", "solid")) +
scale_size_manual(values=c(0.5, 0.5, 1.5)) +
xlab("Estradiol concentration (nM)") + ylab("Fluorescence") + theme_light() + scale_x_continuous(trans="log2", breaks = c(1, 2, 4, 5, 6, 8, 10, 15, 20, 30, 40)))

print(ggplot(datMeans[datMeans$estradiol >= 4,], aes(x=estradiol, y=expression, color=tag)) +
geom_line(aes(linetype=genotype, size=genotype)) +
geom_pointrange(aes(ymin=expression-SD, ymax=expression+SD), position=position_jitter(width=0.01)) +
scale_color_manual(values=c("green", "red")) +
scale_linetype_manual(values=c("dotted", "twodash", "solid")) +
scale_size_manual(values=c(0.5, 0.5, 1.5)) +
xlab("Estradiol concentration (nM)") + ylab("Fluorescence") + theme_light() + scale_x_continuous(trans="log2", breaks = c(1, 2, 4, 5, 6, 8, 10, 15, 20, 30, 40)))

print(ggplot(datMeans, aes(x=estradiol, y=expression, color=tag)) +
geom_line(aes(linetype=genotype, size=genotype)) +
geom_pointrange(aes(ymin=expression-SD, ymax=expression+SD), position=position_jitter(width=0.01)) +
scale_color_manual(values=c("green", "red")) +
scale_linetype_manual(values=c("dotted", "twodash", "solid")) +
scale_size_manual(values=c(0.5, 0.5, 1.5)) +
xlab("Estradiol concentration (nM)") + ylab("Fluorescence") + theme_light() + scale_x_continuous(trans="log2", breaks = c(0.5, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40)))

print(ggplot(datMeans, aes(x=estradiol, y=expression, color=tag)) +
geom_line(aes(linetype=genotype, size=genotype)) +
geom_pointrange(aes(ymin=expression-SD, ymax=expression+SD), position=position_jitter(width=0.01)) +
scale_color_manual(values=c("green", "red")) +
scale_linetype_manual(values=c("dotted", "twodash", "solid")) +
scale_size_manual(values=c(0.5, 0.5, 1.5)) +
xlab("Estradiol concentration (nM)") + ylab("Fluorescence") + theme_light())

print(ggplot(datMeansGrowth[datMeansGrowth$estradiol >= 4 & datMeansGrowth$estradiol < 30,], aes(x=estradiol, y=growth)) +
geom_line(aes(linetype=genotype, size=genotype)) +
geom_pointrange(aes(ymin=growth-SDGrowth, ymax=growth+SDGrowth), position=position_jitter(width=0.01)) +
scale_linetype_manual(values=c("dotted", "twodash", "solid")) +
scale_size_manual(values=c(0.5, 0.5, 1.5)) +
xlab("Estradiol concentration (nM)") + ylab("Growth rate") + theme_light() + scale_x_continuous(trans="log2", breaks = c(0.5, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40)))

print(ggplot(datMeansGrowth[datMeansGrowth$estradiol >= 3,], aes(x=estradiol, y=growth)) +
geom_line(aes(linetype=genotype, size=genotype)) +
geom_pointrange(aes(ymin=growth-SDGrowth, ymax=growth+SDGrowth), position=position_jitter(width=0.01)) +
scale_linetype_manual(values=c("dotted", "twodash", "solid")) +
scale_size_manual(values=c(0.5, 0.5, 1.5)) +
xlab("Estradiol concentration (nM)") + ylab("Growth rate") + theme_light() + scale_x_continuous(trans="log2", breaks = c(0.5, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40)))

print(ggplot(datMeansGrowth, aes(x=estradiol, y=growth)) +
geom_line(aes(linetype=genotype, size=genotype)) +
geom_pointrange(aes(ymin=growth-SDGrowth, ymax=growth+SDGrowth), position=position_jitter(width=0.01)) +
scale_linetype_manual(values=c("dotted", "twodash", "solid")) +
scale_size_manual(values=c(0.5, 0.5, 1.5)) +
xlab("Estradiol concentration (nM)") + ylab("Growth rate") + theme_light() + scale_x_continuous(trans="log2", breaks = c(0.5, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40)))

dev.off()
