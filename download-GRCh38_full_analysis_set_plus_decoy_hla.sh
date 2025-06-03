#!/bin/bash

# Download GRCh38 full analysis set plus decoy HLA reference genome
# This script downloads the FASTA file and its index from the 1000 Genomes FTP site

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

# Base URL for 1000 Genomes GRCh38 reference
BASE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome"

# File names
FASTA_FILE="GRCh38_full_analysis_set_plus_decoy_hla.fa"
INDEX_FILE="GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"

# Create directory for reference files
mkdir -p reference_genome
cd reference_genome

echo "Downloading GRCh38 reference genome files..."
echo "This may take a while as the FASTA file is approximately 3.2GB"
echo

# Download the FASTA file
echo "Downloading ${FASTA_FILE}..."
if ! wget --progress=bar:force:noscroll "${BASE_URL}/${FASTA_FILE}"; then
    echo "Error: Failed to download ${FASTA_FILE}"
    exit 1
fi

# Download the FASTA index file
echo "Downloading ${INDEX_FILE}..."
if ! wget --progress=bar:force:noscroll "${BASE_URL}/${INDEX_FILE}"; then
    echo "Error: Failed to download ${INDEX_FILE}"
    exit 1
fi

echo
echo "Download completed successfully!"
echo "Files downloaded to $(pwd):"
echo "  - ${FASTA_FILE} ($(du -h ${FASTA_FILE} | cut -f1))"
echo "  - ${INDEX_FILE} ($(du -h ${INDEX_FILE} | cut -f1))"
echo

# Verify files exist and are not empty
if [[ -s "${FASTA_FILE}" && -s "${INDEX_FILE}" ]]; then
    echo "✓ Both files downloaded successfully and are not empty"
else
    echo "✗ Error: One or both files are missing or empty"
    exit 1
fi

echo
echo "Reference genome files are ready for use!"
echo "FASTA file: $(pwd)/${FASTA_FILE}"
echo "Index file: $(pwd)/${INDEX_FILE}"
