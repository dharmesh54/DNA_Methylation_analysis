#!/bin/bash

# Define file paths
NORMAL_FILTERED_VCF="/mnt/c/Users/dharm/Downloads/Pupil_bio/filtered/PA221MH-lib09-P19-Norm_S1_L001.sorted_filtered.vcf.gz"
TUMOR_FILTERED_VCF="/mnt/c/Users/dharm/Downloads/Pupil_bio/filtered/PA220KH-lib09-P19-Tumor_S2_L001.sorted_filtered.vcf.gz"
OUTPUT_DIR="/mnt/c/Users/dharm/Downloads/Pupil_bio/filtered/output/"


# Create output directory
mkdir -p $OUTPUT_DIR

# Find variants present in tumor but absent in normal
bcftools isec \
    -C \
    $TUMOR_FILTERED_VCF \
    $NORMAL_FILTERED_VCF \
    -Ob \
    -o $OUTPUT_DIR/tumor_only.vcf

# Compress the resulting VCF in BGZF format
bgzip -c $OUTPUT_DIR/tumor_only.vcf > $OUTPUT_DIR/tumor_only.vcf.gz


# Index the resulting VCF
bcftools index $OUTPUT_DIR/tumor_only.vcf.gz
