# Frank 07/06/18:
# wrapping Maggie's code into a function
# other changes:
# catch cases where cells don't grow at all
# catch cases where cells don't get within 0.05 of the user defined threshold => don't report
# remove white space within catDesc

# load packages upon source
library(readxl)
library(stringr)
library(growthcurver)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)


runMaggieOnPlates <- function(dataFile, descriptionFile, outTextFile, outPlotFile="tempOUT_PLOT_MAGGIE.pdf", userOD=0.25, outputFolder=".", inputFolder="."){

outTextFile <- paste(outputFolder, outTextFile, sep="/")
outPlotFile <- paste(outputFolder, outPlotFile, sep="/")

# Midlog Ratio 3
# Midlog Analysis - with Desc Excel Doc
# New way to find MidLog 
# Au: Maggie Kliebhan


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read in Plate Reader Data


# Plate Reader Output Excel File

# Set this file path to the location of the excel file
#setwd("/Users/lutz0006/Documents/AlbertLab/CausalVariants/IRA2/PlateReaderResults_IRA2")


# Set this text to the name of the excel file 
file <- paste(inputFolder, dataFile, sep="/")

# Read in the Description File
# Set to location of desc. excel sheet
#setwd("/Users/lutz0006/Documents/AlbertLab/CausalVariants/IRA2/PlateReaderResults_IRA2")

# Specify File Name (and Sheet if applicable)

file_desc <- paste(inputFolder, descriptionFile, sep="/")

sheet_name <- 1 
# set to name or number
# defaults to sheet 1

# Create Description Dataframe
desc.df <- read_excel(file_desc, sheet = sheet_name, col_names = TRUE)



#**********

# Install the following package - ONLY RUN THIS ONCE, comment it out after first use
# install.packages("tidyverse")



#**********

# Read in the Data 


# Location of OD table
OD_loc <- which(read_excel(file, col_names = FALSE, range= "A1:A700") == "600")

# Create tables

# OD
OD <- read_excel(file, col_names=TRUE, range = paste("B", OD_loc +2,":CU", OD_loc +99, sep=""))    #assuming plate size to be standard size

# GFP
GFP <- read_excel(file, col_names=TRUE, range = paste("B", OD_loc +103,":CU", OD_loc +200, sep=""))

#Remove Extra Data
rm(OD_loc)




# Store OD and GFP Time for Grofit
OD_time <- OD$Time
GFP_time <- GFP$Time



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




# Define Blanks

# Blanks have contain either the word "blank" or a "Y" in their title

if("Blank" %in% colnames(desc.df)) { 
  # Only run if Blank Column Exists
  
  contains_blank <- grep("blank", desc.df$Blank, ignore.case = TRUE)
  contains_Y <- grep("Y", desc.df$Blank, ignore.case = TRUE)
  
  blank_num <- unique(c(contains_blank,contains_Y))
  blanks <- desc.df$Wells[blank_num]
  
  if(identical(blanks, character(0))) {warning("Blanks not defined - Use manual entry or change desc.df")}
  
  rm(contains_blank, contains_Y, blank_num)
  
}  else {
  warning("Blank column not defined - please edit excel sheet")
}



# Define Wells 
wells <- setdiff(desc.df$Wells, blanks) # does not include blanks in wells



# # Option to define wells
# wells <- c("A1", "A2", ...)

# # Option to remove specific wells
# removewells <- c("A1", "A2")

removewells <- desc.df$Wells[is.na(desc.df$Desc_1)]

wells <- wells[!is.element(wells, removewells)]




#############
# FRANK: Growthcurver inserted here, before OD gets modified below

# make a data.frame for growthcurver
OD.df <- data.frame(OD[,wells])
# time in hours
OD.df$time <- as.numeric((OD$Time - OD$Time[1])/60/60)


# aggregate all blank columns into one and add them to the df
OD.df$blank <- rowMeans(OD[,blanks])

# sometimes the plate did not finish the whole 24h; need to remove these time points
OD.df <- OD.df[sign(OD.df$time) >= 0,]
#print(OD.df)

gc_out <- SummarizeGrowthByPlate(OD.df, bg_correct = "blank", plot_fit="TRUE", plot_file=paste0(outPlotFile, "_growthCurver.pdf", sep=""))
rownames(gc_out) <- gc_out$sample
gc_out <- gc_out[gc_out$sample != "",]
gc_out[gc_out$note != "", 2:9] <- NA

# Determine Time Points at growthcurver's t_mid

rowIndicesTMid <- function(Data, GCOut, numberOfTimePoints = 3){
    Data <- Data[,rownames(GCOut)]
    flank <- floor((numberOfTimePoints-1)/2)   # Translate to the amount flanking on each side (prevents even numbers)
    indices <- sapply(colnames(Data), function(x){
        # catch cases where t_mid has a "note" and returns nothing
        if(GCOut[x, "note"] != ""){return(rep(NA, numberOfTimePoints))}
        # which timepoint is closest to t_mid?
        timeFromMid <- abs(OD.df$time - GCOut[x, "t_mid"])
        # set points after the max to NA
        max_index <- min(which(Data[,x]==max(Data[,x][3:length(Data[,x])], na.rm = TRUE)))
        timeFromMid[max_index:length(Data[,x])] <- NA
        
        k <- which(timeFromMid == min(timeFromMid, na.rm = TRUE))[1]
        if(is.na(k)){return(rep(NA, numberOfTimePoints))} # FRANK catch if well didn't grow
        ret <- (k-flank):(k+flank)
        if(ret[1] == 0){ret <- ret + 1} # necessary because this can sometimes pick the first point, then first index is 0, which cannot be
        ret
    })
    return(indices)
}

rowIndicesIwantTMid <- rowIndicesTMid(OD.df[,rownames(gc_out)], gc_out)

# Maggie's code resumes
#############






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Baseline Correct with Blank Data


# Running Well 
#   = a mean of all blank wells for each time point
#     can be changed to median

OD_avg_blankwell <- apply(OD, 1, function(x){
  mean(as.numeric(x[blanks]))
  # median(as.numeric(x[blanks]))
})


GFP_avg_blankwell <- apply(GFP, 1, function(x){
  mean(as.numeric(x[blanks]))
  # median(as.numeric(x[blanks]))
})




# Grand = avg all time points for the running well
grdCorr_OD <- mean(OD_avg_blankwell, na.rm = T)

grdCorr_GFP <- mean(GFP_avg_blankwell, na.rm = T)




# Blank Correct the Data

# Grand Corr
OD <- OD - grdCorr_OD
GFP <- GFP - grdCorr_GFP



# # Running Correction
# for (i in 1:nrow(OD)){
#   OD[i, ] <- OD[i, ] - OD_avg_blankwell[i]
# }
# for (i in 1:nrow(GFP)){
#   GFP[i, ] <- GFP[i, ] - GFP_avg_blankwell[i]
# }



#********

# Filter tables to contain only relevant wells
OD <- OD[,colnames(OD) %in% wells]
GFP <- GFP[,colnames(GFP) %in% wells]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# MidLog Ratios



# Determine Time Points at Mid Log

rowIndicesMidlog <- function(Data, numberOfTimePoints = 3){
  flank <- floor((numberOfTimePoints-1)/2)   # Translate to the amount flanking on each side (prevents even numbers)
  indices<- apply(Data, 2, function(x){
    log_OD <- log(x)
    log_OD[1:2] <- NA # Remove first two point for jumpy data
    log_OD[log_OD < -5] <- NA # Want to remove very small data (inaccurate)
    # FRANK: when cells didn't grow, (as in F8 on 18/04/09), this errors
    # we should throw such wells out
    if (length(which(!is.na(log_OD))) == 0){return(rep(NA, numberOfTimePoints))}
    
    max_log <- max(log_OD, na.rm = TRUE)
    min_log <- min(log_OD, na.rm = TRUE)
    halfmax <- mean(c(max_log, min_log))
    
    abs_dist <- abs(log_OD - halfmax)
    max_point <- which(log_OD == max_log)[1]
    abs_dist[max_point:length(x)] <- NA
    # FRANK: this can become all NAs, throw out this wekk
    if (length(which(!is.na(abs_dist))) == 0){return(rep(NA, numberOfTimePoints))}
    
    k <- which(abs_dist == min(abs_dist, na.rm = TRUE))[1]
    (k-flank):(k+flank)
 })
 return(indices)
}

# Set Time Points
rowIndicesIwantMidlog <- rowIndicesMidlog(OD)
# rowIndicesIwantMidLog <- rowIndicesMidLog(OD, numberOfTimePoints = )



#**********

# Calculate the Ratios at Mid Log

OD_values_Midlog <- sapply(colnames(OD), function(x){OD[rowIndicesIwantMidlog[,x],x]})

# Run the following
GFP_values_Midlog <- sapply(colnames(OD), function(x){GFP[rowIndicesIwantMidlog[,x],x]})



# Store Output
Ratios_Midlog <- apply(GFP_values_Midlog/OD_values_Midlog, 2, mean)

#Can also take a median!


#**********

# Calculate the Ratios at t_mid (i.e., the end of log, the "inflection point" according to the growthcurver manual)
# pragmatically, these are typically very close to 0.25, sometimes a bit higher

OD_values_TMid <- sapply(colnames(OD), function(x){OD[rowIndicesIwantTMid[,x],x]})
GFP_values_TMid <- sapply(colnames(OD), function(x){GFP[rowIndicesIwantTMid[,x],x]})

# Store Output
Ratios_TMid <- apply(GFP_values_TMid/OD_values_TMid, 2, mean)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Determine Time Points at a User-Defined OD

rowIndicesUser <- function(Data, userDefined, numberOfTimePoints = 3){
  flank <- floor((numberOfTimePoints-1)/2)   # Translate to the amount flanking on each side (prevents even numbers)
  indices <- apply(Data, 2, function(x){
    abs_dist <- abs(x-userDefined)
    # FRANK: catch cases where cells don't get to user defined (too low or too high from the start)
    if(min(abs_dist, na.rm = TRUE) > 0.05){return(rep(NA, numberOfTimePoints))}
    abs_dist[1:2] <- NA
    max_index <- min(which(x==max(x[3:length(x)], na.rm = TRUE)))
    abs_dist[max_index:length(x)] <- NA
    k <- which(abs_dist == min(abs_dist, na.rm = TRUE))[1]
    if(is.na(k)){return(rep(NA, numberOfTimePoints))} # FRANK catch if well didn't grow
    (k-flank):(k+flank)
  })
  return(indices)
}

rowIndicesIwantUser <- rowIndicesUser(OD, userDefined = userOD)

#**********

# Calculate the Ratios at Mid Log

OD_values_User <- sapply(colnames(OD), function(x){OD[rowIndicesIwantUser[,x],x]})

# Run the following
GFP_values_User <- sapply(colnames(OD), function(x){GFP[rowIndicesIwantUser[,x],x]})



# Store Output
Ratios_User <- apply(GFP_values_User/OD_values_User, 2, mean)

#Can also take a median!




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Determine Time Points at Saturation

# Function
rowIndicesSat <- function(Data, numberOfTimePoints = 3){
  flank <- floor((numberOfTimePoints-1)/2)   # Translate to the amount flanking on each side (prevents even numbers)
  indices <- apply(OD, 2, function(x){
    abs_dist <- x
    abs_dist[1:2] <- NA # remove to avoid jumps
    k <- min(which(x==max(x[3:length(x)], na.rm = TRUE)))[1] # k is timepoint at max OD
    # FRANK: If k is in the first half, something is fishy => do not return
    if (k < (length(x)/2)){return(rep(NA, numberOfTimePoints))}
    if(k+flank > length(x)) {k <- length(x) - flank} # Make sure that range remains within points
    floor((k-flank):(k+flank))
  })
  return(indices)
}



# Set Row Indices at Sat
rowIndicesIwantSat <- rowIndicesSat(OD)

# rowIndicesIwantSat <- rowIndicesSat(OD, numberOfTimePoints = )

#**********

# Calculate the Ratios at Saturation

OD_values_Sat <- sapply(colnames(OD), function(x){OD[rowIndicesIwantSat[,x],x]})


GFP_values_Sat <- sapply(colnames(OD), function(x){GFP[rowIndicesIwantSat[,x],x]})


Ratios_Sat <- apply(GFP_values_Sat/OD_values_Sat, 2, mean)
#Can also take a median!




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Output data


# Ratios in order of desc.df
#Ratios_Midlog[desc.df$Wells]


# Combined Desc Row
catDesc <- apply(desc.df[3:7], 1, function(x){
  x[is.na(x)] <- "" # store NA as nothing
  x <- str_trim(x)
  cat <- paste(x, collapse = "_") # paste together Desc columns
  #cat <- gsub(" *", "", cat)
  return(cat)
})
# IMPORTANT: Stored in order of wells in desc.df




#*********

# Store Output Data
ratios.df <- data.frame(
  "Wells" = desc.df$Wells, 
  catDesc, 
  
  "Ratios_Midlog" = Ratios_Midlog[desc.df$Wells], 
  "Ratios_UserDef" = Ratios_User[desc.df$Wells], 
  "Ratios_Sat" = Ratios_Sat[desc.df$Wells],
  "Ratios_TMid" = Ratios_TMid[desc.df$Wells],
  
  "OD_Midlog" = apply(OD_values_Midlog, 2, mean)[desc.df$Wells],
  "GFP_Midlog" = apply(GFP_values_Midlog, 2, mean)[desc.df$Wells],
  "OD_User" = apply(OD_values_User, 2, mean)[desc.df$Wells],
  "GFP_User" = apply(GFP_values_User, 2, mean)[desc.df$Wells],
  "OD_Sat" = apply(OD_values_Sat, 2, mean)[desc.df$Wells],
  "GFP_Sat" = apply(GFP_values_Sat, 2, mean)[desc.df$Wells],
  "OD_TMid" = apply(OD_values_TMid, 2, mean)[desc.df$Wells],
  "GFP_TMid" = apply(GFP_values_TMid, 2, mean)[desc.df$Wells],
  
  "growthRate" = gc_out[desc.df$Wells, "r"],
  "capacity" = gc_out[desc.df$Wells, "k"],
  "generationTime" = gc_out[desc.df$Wells, "t_gen"],

  desc.df)




#********

# # Reorder - Optional
# #   otherwise, will print in order of desc.df
# ratios.df <- ratios.df[order(substring(ratios.df$Well,2), ratios.df$Well),]





#********

# Option to output at this point
# Use only if in a rush 

# Data can be viewed using
#View(ratios.df)

# Otherwise, please continue to process growth data



# Change the name to where you want it to output

# Write a Tab Separated Table
write.table(ratios.df, outTextFile, sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Log OD
logOD <- log(OD)
# Ratio Table
Ratio <- GFP/OD


pdf(outPlotFile, width=8, height=8)
par(mfrow = c(2,2))
# Plot OD
for (i in colnames(OD)) {
  
  # OD 
  plot(OD[,i], main = i, ylab="OD600")
  points(rowIndicesIwantMidlog[,i], OD[rowIndicesIwantMidlog[,i],i], pch = 19, col = "blue")
  points(rowIndicesIwantUser[,i], OD[rowIndicesIwantUser[,i],i], pch = 19, col = "red")
  points(rowIndicesIwantSat[,i], OD[rowIndicesIwantSat[,i],i], pch = 19, col = "forest green")
  points(rowIndicesIwantTMid[,i], OD[rowIndicesIwantTMid[,i],i], col = "orange", cex=1.5)
  legend("topleft", box.lty=0, legend=c("1/2 (max-min)", "user defined", "max", "t_mid"), fill=c("blue", "red", "forest green", "orange"))
  
  # GFP
  plot(GFP[,i], main = "GFP", ylab="GFP")
  points(rowIndicesIwantMidlog[,i], GFP[rowIndicesIwantMidlog[,i],i], pch = 19, col = "blue")
  points(rowIndicesIwantUser[,i], GFP[rowIndicesIwantUser[,i],i], pch = 19, col = "red")
  points(rowIndicesIwantSat[,i], GFP[rowIndicesIwantSat[,i],i], pch = 19, col = "forest green")
  points(rowIndicesIwantTMid[,i], GFP[rowIndicesIwantTMid[,i],i], col = "orange", cex=1.5)
  
  # logOD
  if (length(which(!is.na(logOD[,i]))) == 0){plot.new()}
  else{
      plot(logOD[,i], main = "logOD", ylab="log(OD600)")
    points(rowIndicesIwantMidlog[,i], logOD[rowIndicesIwantMidlog[,i],i], pch = 19, col = "blue")
    points(rowIndicesIwantUser[,i], logOD[rowIndicesIwantUser[,i],i], pch = 19, col = "red")
    points(rowIndicesIwantSat[,i], logOD[rowIndicesIwantSat[,i],i], pch = 19, col = "forest green")
    points(rowIndicesIwantTMid[,i], logOD[rowIndicesIwantTMid[,i],i], col = "orange", cex=1.5)
  }
  
  # Ratios
  plot(Ratio[,i], main = "Ratios", ylab = "GFP / OD ratio", ylim=c(0, max(Ratio[,i][9:length(Ratio[,i])], na.rm=TRUE)))
  points(rowIndicesIwantMidlog[,i], Ratio[rowIndicesIwantMidlog[,i],i], pch = 19, col = "blue")
  points(rowIndicesIwantUser[,i], Ratio[rowIndicesIwantUser[,i],i], pch = 19, col = "red")
  points(rowIndicesIwantSat[,i], Ratio[rowIndicesIwantSat[,i],i], pch = 19, col = "forest green")
  points(rowIndicesIwantTMid[,i], Ratio[rowIndicesIwantTMid[,i],i], col = "orange", cex=1.5)
}
dev.off()

# plot that combines all curves for this plate; colored by genotype
combinedPlotFile <- paste0(strsplit(outPlotFile, ".pdf")[[1]][1], "_allGrowthCurvesCombined.pdf")
plotDatList <- list(OD, GFP, logOD, Ratio)
names(plotDatList) <- c("OD600", "GFP", "log(OD600)", "GFP / OD ratio")
yLowerLimits <- c(-0.01, 0, NA, 0)
names(yLowerLimits) <- names(plotDatList)

pdf(combinedPlotFile, width=20, height=12)
#par(mfrow = c(2,2))
ggPlots <- lapply(names(plotDatList), function(i){
    datForGG <- data.frame(t(plotDatList[[i]]))
    datForGG$background <- desc.df[match(colnames(plotDatList[[i]]), desc.df$Wells),]$Desc_1
    datForGG$genotype <- do.call(paste, desc.df[match(colnames(plotDatList[[i]]), desc.df$Wells),c("Desc_2", "Desc_3", "Desc_5")])
    datForGG$clone <- desc.df[match(colnames(plotDatList[[i]]), desc.df$Wells),]$Wells
    datForGG$WT <- desc.df[match(colnames(plotDatList[[i]]), desc.df$Wells),]$Desc_3
    datForGG$Tag <- desc.df[match(colnames(plotDatList[[i]]), desc.df$Wells),]$Desc_2
    datForGG_long <- melt(datForGG)
    
    plotTime <- OD.df$time
    names(plotTime) <- unique(datForGG_long$variable)[1:length(names(plotTime))]
    datForGG_long$time <- plotTime[datForGG_long$variable]
    
    cols <- colorRampPalette(brewer.pal(9, "Blues"))
    myPal <- rev(cols(length(unique(datForGG_long$genotype))))
    names(myPal) <- sort(unique(datForGG_long$genotype))
    myPal[names(myPal) == "noGFP wt"] <- "dark grey"
    myPal[str_detect(names(myPal), "-GFP wt")] <- "red"
    
    #p <- ggplot(data=datForGG_long[datForGG_long$WT != "wt",],
    p <- ggplot(data=datForGG_long,
        aes(x=time, y=value, colour=genotype, linetype=background, group=clone)) +
        geom_line() + ylim(yLowerLimits[i], max(datForGG_long$value[datForGG_long$time >= 2 & datForGG_long$value < abs(median(datForGG_long$value, na.rm=TRUE))*10], na.rm=TRUE)) +
        xlab("Time (hours)") + ylab(i) + scale_color_manual(values=myPal)
        #geom_line(data=datForGG_long[datForGG_long$WT == "wt" & datForGG_long$Tag == "noGFP",], aes(x=time, y=value), colour="black") +
        #geom_line(data=datForGG_long[datForGG_long$WT == "wt" & datForGG_long$Tag != "noGFP",], aes(x=time, y=value), colour="red")
})
grid.arrange(ggPlots[[1]], ggPlots[[2]], ggPlots[[3]], ggPlots[[4]], nrow=2)

dev.off()

# final return
return(ratios.df)

# final } to close the function:
}


