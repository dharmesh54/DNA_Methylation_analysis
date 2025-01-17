#!/bin/bash

# Set variables
REFERENCE="/mnt/c/Users/dharm/Downloads/Pupil_bio/human_genome/GCA_000001405.29_GRCh38.p14_genomic.fna"
OUTPUT_DIR="/mnt/c/Users/dharm/Downloads/Pupil_bio/sam_bam_files"
THREADS=2
SAMPLE_DIR="/mnt/c/Users/dharm/Downloads/Pupil_bio/Trimming/paired_end/"
mkdir -p $OUTPUT_DIR

# Iterate over Tumor and Normal samples
for SAMPLE_PATH in "$SAMPLE_DIR"/*; do
    SAMPLE_TYPE=$(basename "$SAMPLE_PATH")  # Extract sample type (e.g., Tumor, Normal)
    R1="${SAMPLE_PATH}_R1_paired.fq"
    R2="${SAMPLE_PATH}_R2_paired.fq"

    echo "Processing $SAMPLE_TYPE sample..."

    # Step 1: Align reads to reference genome
    bwa mem -t $THREADS $REFERENCE $R1 $R2 > $OUTPUT_DIR/${SAMPLE_TYPE}.sam
    echo "Alignment complete for $SAMPLE_TYPE."

    echo "$SAMPLE_TYPE sample processing finished."
done

echo "All samples processed. Results are in the $OUTPUT_DIR directory."
