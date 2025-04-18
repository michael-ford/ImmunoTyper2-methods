#!/bin/bash

# Hard-coded paths
REFERENCE_FASTA="concatenated.fasta"
READS_PATH="HG00438_assigned_reads_IGLV2-14_04.fa"

# Define output filenames
OUTPUT_PREFIX="HG00438_assigned_reads_IGLV2-14_comparison_bowtie2"
BOWTIE2_INDEX="${OUTPUT_PREFIX}_index"
OUTPUT_SAM="${OUTPUT_PREFIX}_alignment.sam"
OUTPUT_BAM="${OUTPUT_PREFIX}_alignment.bam"
OUTPUT_SORTED_BAM="${OUTPUT_PREFIX}_alignment.sorted.bam"

# Check if input files exist
if [ ! -f "$REFERENCE_FASTA" ]; then
    echo "Error: Reference FASTA file not found at $REFERENCE_FASTA"
    exit 1
fi

if [ ! -e "$READS_PATH" ]; then
    echo "Error: Reads file or directory not found at $READS_PATH"
    exit 1
fi

echo "=== Bowtie2 Alignment Pipeline ==="
echo "Reference FASTA: $REFERENCE_FASTA"
echo "Reads Path: $READS_PATH"
echo "Output BAM: $OUTPUT_SORTED_BAM"
echo "==============================="

# Step 1: Build Bowtie2 index
echo "[1/4] Building Bowtie2 index..."
bowtie2-build "$REFERENCE_FASTA" "$BOWTIE2_INDEX"

# Step 2: Determine read type (single or paired-end)
# This is a simple check - modify as needed for your specific read naming convention
if ls "$READS_PATH"/*_R2* 1> /dev/null 2>&1 || ls "$READS_PATH"/*_2.* 1> /dev/null 2>&1; then
    echo "[2/4] Detected paired-end reads..."
    # Assuming standard Illumina naming convention (_R1/_R2 or _1/_2)
    READ1=$(ls "$READS_PATH"/*_R1* 2> /dev/null || ls "$READS_PATH"/*_1.* 2> /dev/null)
    READ2=$(ls "$READS_PATH"/*_R2* 2> /dev/null || ls "$READS_PATH"/*_2.* 2> /dev/null)
    
    # Run Bowtie2 with paired-end reads
    echo "Running Bowtie2 with paired-end reads: $READ1 and $READ2"
    # Add the -f flag to specify FASTA format input
    bowtie2 -a --end-to-end --very-sensitive --n-ceil C,100,0 --np 0 --ignore-quals --mp 2,2 --score-min C,-50,0 -L 10 \
        -x "$BOWTIE2_INDEX" -f -1 "$READ1" -2 "$READ2" -S "$OUTPUT_SAM"
else
    echo "[2/4] Assuming single-end reads..."
    
    # If READS_PATH is a directory, use all files; if it's a file, use it directly
    if [ -d "$READS_PATH" ]; then
        READS="$READS_PATH"/*
    else
        READS="$READS_PATH"
    fi
    
    # Run Bowtie2 with single-end reads
    echo "Running Bowtie2 with single-end reads: $READS"
    # Add the -f flag to specify FASTA format input
    bowtie2 -a --end-to-end --very-sensitive --n-ceil C,100,0 --np 0 --ignore-quals --mp 2,2 --score-min C,-50,0 -L 10 \
        -x "$BOWTIE2_INDEX" -f -U "$READS" -S "$OUTPUT_SAM"
fi

# Step 3: Convert SAM to BAM
echo "[3/4] Converting SAM to BAM..."
# Check if SAM file exists and has content before proceeding
if [ -s "$OUTPUT_SAM" ]; then
    samtools view -bS "$OUTPUT_SAM" > "$OUTPUT_BAM"
else
    echo "Error: SAM file is empty or does not exist. Bowtie2 alignment may have failed."
    exit 1
fi

# Step 4: Sort BAM file
echo "[4/4] Sorting BAM file..."
samtools sort "$OUTPUT_BAM" -o "$OUTPUT_SORTED_BAM"

# Optional: Index BAM file
echo "[+] Indexing BAM file..."
samtools index "$OUTPUT_SORTED_BAM"

# Optional: Remove intermediate files
echo "[+] Cleaning up intermediate files..."
rm -f "$OUTPUT_SAM" "$OUTPUT_BAM" "$BOWTIE2_INDEX".*

echo "=== Pipeline Complete ==="
echo "Final sorted BAM file: $OUTPUT_SORTED_BAM"
echo "BAM index file: ${OUTPUT_SORTED_BAM}.bai"
