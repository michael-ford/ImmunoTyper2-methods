#!/bin/bash

# Create output directory if it doesn't exist
mkdir -p ./prefix-consistency

# Get gene types dynamically
GENE_TYPES=$(ls -1 ../../HPRC-assembly-benchmarking/immunotyper-output/ | grep "v")

# Alternatively, use the hardcoded list if the dynamic approach fails
if [ -z "$GENE_TYPES" ]; then
    echo "Using hardcoded gene type list"
    GENE_TYPES=(
        "ighv"
        "igkv"
        "iglv"
        "trav"
        "trbv"
        "trdv"
        "trgv"
    )
fi

# Loop through each gene type
for GENE_TYPE in $GENE_TYPES; do
    echo "Processing $GENE_TYPE..."
    
    # Define paths
    # Convert gene type to uppercase for allele database and remove extra dash
    GENE_TYPE_UPPER=$(echo "$GENE_TYPE" | tr '[:lower:]' '[:upper:]')
    ALLELE_DB="../../ImmunoTyper2/immunotyper/data/allele_databases/${GENE_TYPE_UPPER}/${GENE_TYPE_UPPER}-IMGT-allele-db.fa"
    OUTPUT_DIR="./prefix-consistency/${GENE_TYPE}"
    BASE_DIR="../../HPRC-assembly-benchmarking/immunotyper-output"
    
    # Create output directory for this gene type
    mkdir -p "$OUTPUT_DIR"
    
    # Run the Python script with the appropriate parameters and redirect output
    python prefix_consistency_processing.py \
        --allele-db "$ALLELE_DB" \
        --output "$OUTPUT_DIR" \
        "$BASE_DIR" \
        "$GENE_TYPE" \
        >> make_prefix_consistency_results.stdout \
        2>> make_prefix_consistency_results.stderr
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "Successfully processed $GENE_TYPE"
    else
        echo "Error processing $GENE_TYPE"
    fi
    
    echo "------------------------"
done

echo "All gene types processed."