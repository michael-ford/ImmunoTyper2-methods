#!/bin/bash
# generate_sweep_commands.sh
# Generates commands for parameter sweep of immunotyper-SR

# Default parameters
wgs_dir="/data/fordmk/ImmunoTyper-expansion-methods/HPRC-assembly-benchmarking/wgs-samples"
bam_reference_type="hg38"
reference_fasta="/data/fordmk/ImmunoTyper-expansion-methods/GRCh38_full_analysis_set_plus_decoy_hla.fa"
output_base_dir="./"
runner_script="./run_immunotyper_sweep.sh"

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --wgs_dir)
      wgs_dir="$2"
      shift 2
      ;;
    --bam_reference_type)
      bam_reference_type="$2"
      shift 2
      ;;
    --reference_fasta)
      reference_fasta="$2"
      shift 2
      ;;
    --output_base_dir)
      output_base_dir="$2"
      shift 2
      ;;
    --runner_script)
      runner_script="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

# Print header as a comment
echo "# Generated parameter sweep commands $(date)" >&2
echo "# wgs_dir: ${wgs_dir}" >&2
echo "# output_base_dir: ${output_base_dir}" >&2
echo "# Using runner script: ${runner_script}" >&2

# Gene types to process
gene_types=("ighv" "igkv" "iglv" "trav" "trbv" "trdv" "trgv")

# Parameter sweep values (coarse grid)
landmark_groups=(4 8)
landmarks_per_group=(4 8)
stdev_coeffs=(1.0 2.0)

# Check if directory exists
if [ ! -d "$wgs_dir" ]; then
  echo "Error: Directory $wgs_dir does not exist" >&2
  exit 1
fi

# Find all CRAM files
cram_files=("${wgs_dir}"/*.cram)

# Check if any CRAM files were found
if [ ! -f "${cram_files[0]}" ]; then
  echo "No CRAM files found in $wgs_dir" >&2
  exit 1
fi

# Get full path of runner script
runner_script_path=$(readlink -f "${runner_script}")

# Process each CRAM file
for cram_file in "${cram_files[@]}"; do
  # Check if index exists
  if [ ! -f "${cram_file}.crai" ]; then
    echo "Warning: Index file not found for CRAM: ${cram_file}" >&2
    continue
  fi
  
  # Generate command for each gene type and parameter combination
  for gene_type in "${gene_types[@]}"; do
    for lg in "${landmark_groups[@]}"; do
      for lpg in "${landmarks_per_group[@]}"; do
        for stdev in "${stdev_coeffs[@]}"; do
          # Generate the command to run the immunotyper script with these parameters
          echo "bash \"${runner_script_path}\" --input \"${cram_file}\" --gene_type \"${gene_type}\" --landmark_groups ${lg} --landmarks_per_group ${lpg} --stdev_coeff ${stdev} --bam_reference_type \"${bam_reference_type}\" --reference_fasta \"${reference_fasta}\" --output_base_dir \"${output_base_dir}\""
        done
      done
    done
  done
done

echo "# Total commands: $(echo "$((${#cram_files[@]} * ${#gene_types[@]} * ${#landmark_groups[@]} * ${#landmarks_per_group[@]} * ${#stdev_coeffs[@]}))")" >&2