#!/bin/bash

# Directory containing the sequence reads
READS_DIR="/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/genotype/sequence_reads/"

# Output file for SAMPLEIDs
SAMPLEIDS_FILE="sample_ids.txt"

# Find unique SAMPLEIDs
ls ${READS_DIR}/*_1.fastq.gz | sed 's/_1.fastq.gz//' | xargs -n 1 basename > ${SAMPLEIDS_FILE}

