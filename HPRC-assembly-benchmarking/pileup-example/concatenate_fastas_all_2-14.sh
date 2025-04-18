#!/bin/bash

# Define output files
OUTPUT_FASTA="concatenated-all.fasta"
OUTPUT_BED="annotations-all.bed"

# Find all matching FASTA files and sort them
FASTA_FILES=($(ls iglv2-14_*.fasta | sort))

# Check if any matching files exist
if [ ${#FASTA_FILES[@]} -eq 0 ]; then
    echo "Error: No files matching iglv2-14_*.fasta found!"
    exit 1
fi

# Initialize variables
FINAL_SEQ=""
CURRENT_POS=100  # Start after initial N padding
BED_CONTENT=""
HEADERS=()
IDS=()

# Add initial padding
START_N=$(printf 'N%.0s' $(seq 1 100))
FINAL_SEQ="${START_N}"

# Process each file
for ((i=0; i<${#FASTA_FILES[@]}; i++)); do
    FASTA="${FASTA_FILES[$i]}"
    
    # Extract sequence header
    HEADER=$(grep ">" "$FASTA" | sed 's/>//')
    HEADERS+=("$HEADER")
    
    # Extract ID from header
    ID=$(echo "$HEADER" | cut -d'|' -f2)
    IDS+=("$ID")
    
    # Extract sequence (removing headers and newlines) and convert to uppercase
    SEQ=$(grep -v ">" "$FASTA" | tr -d '\n' | tr '[:lower:]' '[:upper:]')
    
    # Calculate length
    LEN=${#SEQ}
    
    # Add sequence to final sequence
    FINAL_SEQ="${FINAL_SEQ}${SEQ}"
    
    # Calculate positions for BED file
    SEQ_START=$CURRENT_POS
    SEQ_END=$((SEQ_START + LEN))
    
    # Add entry to BED content
    BED_CONTENT="${BED_CONTENT}concatenated_${IDS[0]}\t$SEQ_START\t$SEQ_END\t$HEADER\t1000\t+\t$SEQ_START\t$SEQ_END\t0,0,255\n"
    
    # Update current position
    CURRENT_POS=$SEQ_END
    
    # Add padding between sequences (except after the last one)
    if [ $i -lt $((${#FASTA_FILES[@]} - 1)) ]; then
        MIDDLE_N=$(printf 'N%.0s' $(seq 1 300))
        FINAL_SEQ="${FINAL_SEQ}${MIDDLE_N}"
        CURRENT_POS=$((CURRENT_POS + 300))
    fi
done

# Add final padding
END_N=$(printf 'N%.0s' $(seq 1 100))
FINAL_SEQ="${FINAL_SEQ}${END_N}"

# Create concatenated IDs string
CONCAT_IDS=$(IFS="_"; echo "${IDS[*]}")

# Create output FASTA file
echo ">concatenated_${CONCAT_IDS}" > "$OUTPUT_FASTA"
# Format sequence with 60 characters per line
echo "$FINAL_SEQ" | fold -w 60 >> "$OUTPUT_FASTA"

# Create BED file for IGV
echo -e "${BED_CONTENT}" > "$OUTPUT_BED"

echo "Files created:"
echo "  - $OUTPUT_FASTA (concatenated FASTA)"
echo "  - $OUTPUT_BED (BED annotation file for IGV)"
echo ""
echo "FASTA structure:"
echo "  - First 100 bases: N padding"

# Describe each sequence in the structure
CURRENT_POS=100
for ((i=0; i<${#FASTA_FILES[@]}; i++)); do
    FASTA="${FASTA_FILES[$i]}"
    SEQ=$(grep -v ">" "$FASTA" | tr -d '\n')
    LEN=${#SEQ}
    
    SEQ_START=$CURRENT_POS
    SEQ_END=$((SEQ_START + LEN))
    
    echo "  - Next $LEN bases: ${HEADERS[$i]} sequence (position $SEQ_START-$SEQ_END)"
    
    CURRENT_POS=$SEQ_END
    
    if [ $i -lt $((${#FASTA_FILES[@]} - 1)) ]; then
        echo "  - Next 300 bases: N padding"
        CURRENT_POS=$((CURRENT_POS + 300))
    fi
done

echo "  - Final 100 bases: N padding"

# Index the FASTA file
samtools faidx "$OUTPUT_FASTA"