#!/usr/bin/env bash

#BSUB -J tb_download # job name
#BSUB -n 4          # number of tasksi/CPUs in job
#BSUB -q standard   # queue
#BSUB -e error.%J   # error file name in which %J is replaced by the job ID
#BSUB -o output.%J  # output file name in which %J is replaced by the job ID

#################################################################################

# Bash script to download sequence reads (fastq format)

# Make text file with the names of sample accession IDs that you want to download
# Here it's column 40 of a data file called CRyPTIC_reuse_table_20221019.csv

awk -F "\"*,\"*" '{print $40}' /nfs/research/jlees/jacqueline/gwas_datasets/tuberculosis/CRyPTIC_reuse_table_20221019.csv > CRyPTIC_accession_ids.txt

# This line just removes the first line since there is a header on the .csv I extracted the IDs from
sed -i '1d' CRyPTIC_accession_ids.txt

# Using ffq, read through the sample IDs and make a file with their ftp URLs

cat CRyPTIC_accession_ids.txt | while read line || [[ -n $line ]]; do ffq --ftp $line | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' >> complete_file.txt; done

# Remove quotes from the URLs in the file
sed -i 's/"//g' complete_file.txt

# Use wget to download the fastq files from the URLs in the previously created file

cat complete_file.txt | while read line || [[ -n $line ]]; do wget -nc $line ; done

# Completion notice

echo "Job successfully completed. That's pretty cool!"

