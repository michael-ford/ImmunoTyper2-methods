#!/usr/bin/env python3
import pandas as pd
import sys
import argparse
import os
import glob
import re
from pathlib import Path
from Bio import SeqIO

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Extract sequences from DIGGER CSV files and IGLV2-14 FASTA files.')
    parser.add_argument('--basedir', default='../digger-results/', help='Base directory containing DIGGER results')
    parser.add_argument('--sample_id', default='HG00438', help='Sample ID to process')
    args = parser.parse_args()

    # Find directories containing the sample_id
    basedir = Path(args.basedir)
    matching_dirs = [d for d in basedir.glob('*') if d.is_dir() and args.sample_id in d.name]

    if not matching_dirs:
        print(f"No directories containing '{args.sample_id}' found in '{basedir}'")
        sys.exit(1)

    print(f"Found {len(matching_dirs)} directory(ies) matching '{args.sample_id}'")

    # Create a single FASTA file in the current directory
    output_fasta = Path(f"{args.sample_id}_2-14_sequences.fasta")

    # Collect all sequences
    sequence_count = 0

    with open(output_fasta, 'w') as fasta_out:
        for dir_path in matching_dirs:
            print(f"Processing directory: {dir_path}")

            # Find the IGL subdirectory
            igl_dir = dir_path / 'IGL'
            if not igl_dir.exists() or not igl_dir.is_dir():
                print(f"IGL directory not found in {dir_path}")
                continue

            # Process CSV files
            csv_files = list(igl_dir.glob('*.csv'))
            if csv_files:
                print(f"Found {len(csv_files)} CSV file(s) in {igl_dir}")
                for csv_file in csv_files:
                    print(f"Processing CSV file: {csv_file}")

                    # Extract the sample prefix from the directory name
                    dir_name = dir_path.name
                    match = re.match(r'([^.]+\.[^.]+)\.', dir_name)
                    if match:
                        sample_prefix = match.group(1)
                    else:
                        sample_prefix = dir_name

                    try:
                        df = pd.read_csv(csv_file)
                    except Exception as e:
                        print(f"Error reading {csv_file}: {e}")
                        continue

                    if 'seq' not in df.columns:
                        print(f"Required column 'seq' not found in {csv_file}")
                        continue

                    if 'imgt_match' not in df.columns:
                        print(f"Required column 'imgt_match' not found in {csv_file}")
                        print(f"Available columns: {', '.join(df.columns)}")
                        continue

                    mask = df.astype(str).apply(lambda x: x.str.contains('2-14', na=False)).any(axis=1)
                    filtered_df = df[mask]

                    if not filtered_df.empty:
                        file_seq_count = len(filtered_df)
                        print(f"Found {file_seq_count} row(s) containing '2-14'")
                        for _, row in filtered_df.iterrows():
                            clean_sample_prefix = sample_prefix.replace('.', '_').replace('*', '_')
                            clean_imgt_match = str(row['imgt_match']).replace('.', '_').replace('*', '_')
                            header = f">{clean_sample_prefix.replace('HG00438_', '')}-{clean_imgt_match}_c={row['contig'].replace('hg004382jahbc', '')}_s={row['start']}_e={row['end']}_{row['functional']}"
                            sequence = row['seq']
                            fasta_out.write(f"{header}\n{sequence}\n")
                            sequence_count += 1
                    else:
                        print(f"No rows containing '2-14' found in {csv_file}")
            else:
                print(f"No CSV files found in {igl_dir}")

        # Process FASTA files with pattern iglv2-14_*.fasta in the current directory
        current_dir = Path('.')
        fasta_files = list(current_dir.glob('iglv2-14_*.fasta'))
        if fasta_files:
            print(f"\nFound {len(fasta_files)} FASTA file(s) matching 'iglv2-14_*.fasta' in the current directory.")
            for fasta_file in fasta_files:
                print(f"Processing FASTA file: {fasta_file}")
                try:
                    for record in SeqIO.parse(fasta_file, "fasta"):
                        if "|" in record.id:
                            new_header = record.id.split('|')[1]
                            fasta_out.write(f">{new_header}\n{str(record.seq)}\n")
                            sequence_count += 1
                        else:
                            print(f"Warning: FASTA header '{record.id}' in file '{fasta_file}' does not contain '|'. Skipping.")
                except Exception as e:
                    print(f"Error reading or parsing FASTA file {fasta_file}: {e}")
        else:
            print(f"\nNo FASTA files matching 'iglv2-14_*.fasta' found in the current directory.")

    if sequence_count > 0:
        print(f"\nCreated FASTA file '{output_fasta}' with {sequence_count} sequences")
    else:
        print(f"\nNo matching sequences found. FASTA file '{output_fasta}' is empty.")

if __name__ == "__main__":
    main()