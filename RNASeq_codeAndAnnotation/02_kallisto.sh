#!/bin/bash

# using info from bioanalyzer chip on average library size; guessing at SD
# run one file at a time

inputFiles="Albert007SampleIDs_Indexes.txt" # provided in this code directory
inputFolder="YOUR_OUTPUT_FOLDER/trimmomatic/"
outputFolder="YOUR_OUTPUT_FOLDER/reads/"

while IFS=$'\t' read -r -a lineArray
do
    echo "well: ${lineArray[0]}"
    echo "genotype: ${lineArray[4]}"

    mkdir ${outputFolder}${lineArray[0]}"_"${lineArray[4]}

#ls ${inputFolder}${lineArray[0]}*.fastq.gz

# note that while bioanalyzer library has mean length of 336, this includes adapters – go with the standard 200 instead
# genomebam: the chromosomes file had a line ending issue!

    ~/src/kallisto_linux-v0.44.0/kallisto quant \
        -i Saccharomyces_cerevisiae.R64-1-1.cdna.all.kallistoIndex \ # build this first from the provided file
        -t 28 --single -l 200 -s 30 --rf-stranded -b 100 \
        -o ${outputFolder}${lineArray[0]}"_"${lineArray[4]} \
        --genomebam --gtf Saccharomyces_cerevisiae.R64-1-1.93_mod.gtf.gz \ # provided
        --chromosomes chrLengths_ensemblFormat.txt \ # provided
        ${inputFolder}${lineArray[0]}*.fastq.gz

done < "$inputFiles"
