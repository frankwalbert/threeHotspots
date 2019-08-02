#!/bin/bash

# using info from bioanalyzer chip on average library size; guessing at SD
# need to run one file at a time

inputFiles="Albert007SampleIDs_Indexes.txt" # provided in this code directory
inputFolder="READ_FOLDER"
outputFolder="YOUR_OUTPUT_FOLDER/trimmomatic/"

while IFS=$'\t' read -r -a lineArray
do
    echo "well: ${lineArray[0]}"
    echo "genotype: ${lineArray[4]}"

    java -jar PATH_TO_TRIMMOMATIC/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 14 \
    ${inputFolder}${lineArray[0]}*.fastq.gz \
    ${outputFolder}${lineArray[0]}"_"${lineArray[4]}.fastq.gz \
ILLUMINACLIP:/PATH_TO_TRIMMOMATIC/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

done < "$inputFiles"

# worked well; worst were 2C & 1A with 5.9% of reads dropped
