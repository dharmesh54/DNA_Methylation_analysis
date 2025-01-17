#!/bin/bash

# Paths
REFERENCE="/mnt/c/Users/dharm/Downloads/Pupil_bio/human_genome/GCA_000001405.29_GRCh38.p14_genomic.fna"
BAM_DIR="/mnt/c/Users/dharm/Downloads/Pupil_bio/sam_bam_files/bam_file"
RESULTS_DIR="/mnt/c/Users/dharm/Downloads/Pupil_bio/results"
NORMALIZED_DIR="/mnt/c/Users/dharm/Downloads/Pupil_bio/normalized"
FILTERED_DIR="/mnt/c/Users/dharm/Downloads/Pupil_bio/filtered"
THREADS=4

# Create necessary directories
mkdir -p "$RESULTS_DIR" "$NORMALIZED_DIR" "$FILTERED_DIR"

# Step 1: Process BAM files and generate VCFs
echo "Starting BAM processing..."
for BAM_FILE in "$BAM_DIR"/*.bam; do
    SAMPLE=$(basename "$BAM_FILE" .bam)
    echo "Processing $SAMPLE..."

    # Step 1.1: Mark duplicates using Picard
    java -jar /mnt/c/Users/dharm/Downloads/Pupil_bio/picard.jar MarkDuplicates \
        I="$BAM_FILE" \
        O="$RESULTS_DIR/${SAMPLE}.dedup.bam" \
        M="$RESULTS_DIR/${SAMPLE}_dedup_metrics.txt" \
        REMOVE_DUPLICATES=true

    # Step 1.2: Index the deduplicated BAM file
    samtools index "$RESULTS_DIR/${SAMPLE}.dedup.bam"

    # Step 1.3: Variant calling with bcftools
    bcftools mpileup \
        --threads "$THREADS" \
        --max-depth 250 \
        -f "$REFERENCE" \
        "$RESULTS_DIR/${SAMPLE}.dedup.bam" | \
    bcftools call \
        --multiallelic-caller \
        -A \
        --ploidy GRCh38 \
        -Oz \
        -o "$RESULTS_DIR/${SAMPLE}.vcf.gz"

    # Index the VCF file
    bcftools index "$RESULTS_DIR/${SAMPLE}.vcf.gz"

    echo "$SAMPLE processing completed."
done

# Step 2: Normalize VCF files
echo "Starting normalization of VCF files..."
for VCF_FILE in "$RESULTS_DIR"/*.vcf.gz; do
    SAMPLE=$(basename "$VCF_FILE" .vcf.gz)
    echo "Normalizing $SAMPLE..."

    # Normalize VCF
    bcftools norm \
        -f "$REFERENCE" \
        "$VCF_FILE" \
        -Oz \
        -o "$NORMALIZED_DIR/${SAMPLE}_norm.vcf.gz"

    # Index normalized VCF
    bcftools index "$NORMALIZED_DIR/${SAMPLE}_norm.vcf.gz"

    echo "$SAMPLE normalization completed."
done

# Step 3: Filter VCF files
echo "Starting filtering of VCF files..."
for NORM_VCF in "$NORMALIZED_DIR"/*.vcf.gz; do
    SAMPLE=$(basename "$NORM_VCF" _norm.vcf.gz)
    echo "Filtering $SAMPLE..."

    # Filter VCF
    bcftools filter \
        -g3 -G10 \
        -e '%QUAL<20 || (AC<2 && %QUAL<25) || INFO/DP>99 || INFO/DP<4' \
        "$NORM_VCF" \
        -Oz \
        -o "$FILTERED_DIR/${SAMPLE}_filtered.vcf.gz"

    # Index filtered VCF
    bcftools index "$FILTERED_DIR/${SAMPLE}_filtered.vcf.gz"

    echo "$SAMPLE filtering completed."
done

echo "Pipeline completed! Normalized files are in $NORMALIZED_DIR, and filtered files are in $FILTERED_DIR."
