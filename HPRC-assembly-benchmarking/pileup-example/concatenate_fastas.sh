#!/bin/bash

# Define input files
FASTA1="iglv2-14_03.fasta"
FASTA2="iglv2-14_04.fasta"
OUTPUT_FASTA="concatenated.fasta"
OUTPUT_BED="annotations.bed"

# Check if input files exist
if [ ! -f "$FASTA1" ] || [ ! -f "$FASTA2" ]; then
    echo "Error: Input FASTA files not found!"
    exit 1
fi

# Extract sequence headers
HEADER1=$(grep ">" "$FASTA1" | sed 's/>//')
HEADER2=$(grep ">" "$FASTA2" | sed 's/>//')

# Extract sequences (removing headers and newlines) and convert to uppercase
SEQ1=$(grep -v ">" "$FASTA1" | tr -d '\n' | tr '[:lower:]' '[:upper:]')
SEQ2=$(grep -v ">" "$FASTA2" | tr -d '\n' | tr '[:lower:]' '[:upper:]')

# Calculate lengths
LEN1=${#SEQ1}
LEN2=${#SEQ2}

# Generate padding sequences
START_N=$(printf 'N%.0s' $(seq 1 100))
MIDDLE_N=$(printf 'N%.0s' $(seq 1 300))
END_N=$(printf 'N%.0s' $(seq 1 100))

# Concatenate with padding
FINAL_SEQ="${START_N}${SEQ1}${MIDDLE_N}${SEQ2}${END_N}"

# Calculate positions for BED file
POS1_START=100
POS1_END=$((POS1_START + LEN1))
POS2_START=$((POS1_END + 300))
POS2_END=$((POS2_START + LEN2))

id1=$(echo "$HEADER1" | cut -d'|' -f2)
id2=$(echo "$HEADER2" | cut -d'|' -f2)

# Create output FASTA file
echo ">concatenated_${id1}_${id2}" > "$OUTPUT_FASTA"
# Format sequence with 60 characters per line
echo "$FINAL_SEQ" | fold -w 60 >> "$OUTPUT_FASTA"

# Create BED file for IGV
echo -e "concatenated_${id1}\t$POS1_START\t$POS1_END\t$HEADER1\t1000\t+\t$POS1_START\t$POS1_END\t0,0,255" > "$OUTPUT_BED"
echo -e "concatenated_${id2}\t$POS2_START\t$POS2_END\t$HEADER2\t1000\t+\t$POS2_START\t$POS2_END\t255,0,0" >> "$OUTPUT_BED"

echo "Files created:"
echo "  - $OUTPUT_FASTA (concatenated FASTA)"
echo "  - $OUTPUT_BED (BED annotation file for IGV)"
echo ""
echo "FASTA structure:"
echo "  - First 100 bases: N padding"
echo "  - Next $LEN1 bases: $HEADER1 sequence (position $POS1_START-$POS1_END)"
echo "  - Next 300 bases: N padding"
echo "  - Next $LEN2 bases: $HEADER2 sequence (position $POS2_START-$POS2_END)"
echo "  - Final 100 bases: N padding"

samtools faidx "$OUTPUT_FASTA"