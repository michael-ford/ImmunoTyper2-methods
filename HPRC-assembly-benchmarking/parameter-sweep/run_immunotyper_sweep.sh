#!/bin/bash
# run_immunotyper_sweep.sh
# Script to run immunotyper-SR with specific parameters for parameter sweep

# Parse command line arguments
input_file=""
gene_type=""
landmark_groups=""
landmarks_per_group=""
stdev_coeff=""
bam_reference_type="hg38"
reference_fasta=""
output_base_dir="./parameter_sweep"

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input)
      input_file="$2"
      shift 2
      ;;
    --gene_type)
      gene_type="$2"
      shift 2
      ;;
    --landmark_groups)
      landmark_groups="$2"
      shift 2
      ;;
    --landmarks_per_group)
      landmarks_per_group="$2"
      shift 2
      ;;
    --stdev_coeff)
      stdev_coeff="$2"
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
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Check required arguments
if [ -z "$input_file" ] || [ -z "$gene_type" ] || [ -z "$landmark_groups" ] || [ -z "$landmarks_per_group" ] || [ -z "$stdev_coeff" ]; then
  echo "Error: Required arguments missing"
  echo "Usage: $0 --input <bam/cram file> --gene_type <gene type> --landmark_groups <value> --landmarks_per_group <value> --stdev_coeff <value> [--bam_reference_type <hg38/hg37>] [--reference_fasta <fasta file>] [--output_base_dir <base output directory>]"
  exit 1
fi

# Determine file type
file_extension="${input_file##*.}"
is_cram=false
if [ "$file_extension" == "cram" ]; then
  is_cram=true
  # Check if reference FASTA is provided for CRAM
  if [ -z "$reference_fasta" ]; then
    echo "Error: Reference FASTA is required for CRAM files. Please provide --reference_fasta parameter."
    exit 1
  fi
fi

# Get absolute path of input file
input_path=$(readlink -f "${input_file}")

# Get sample name from file path
sample_name=$(basename "${input_file}" ".${file_extension}")

# Create parameter-specific directory name
param_dir="lg${landmark_groups}_lpg${landmarks_per_group}_stdev${stdev_coeff}"
full_output_dir="${output_base_dir}/${gene_type}/${param_dir}"

# Create output directory
mkdir -p "${full_output_dir}"

# Determine the reference type argument
ref_arg=""
if [ "${bam_reference_type}" == "hg37" ]; then
    ref_arg="--hg37"
fi

# Add reference FASTA argument if input is CRAM
ref_fasta_arg=""
if [ "${is_cram}" == "true" ]; then
    ref_fasta_arg="--ref ${reference_fasta}"
fi

# Print information about the run
echo "Running immunotyper-SR with parameters:"
echo "  Sample: ${sample_name}"
echo "  Gene type: ${gene_type}"
echo "  Landmark groups: ${landmark_groups}"
echo "  Landmarks per group: ${landmarks_per_group}"
echo "  Stdev coefficient: ${stdev_coeff}"
echo "  Output directory: ${full_output_dir}"

# Run immunotyper-SR
cd "${full_output_dir}"
immunotyper-SR --gene_type "${gene_type}" \
  --solver gurobi \
  --multi_band_solutions \
  --no_vcf \
  --no_read_assignment \
  --landmark_groups "${landmark_groups}" \
  --landmarks_per_group "${landmarks_per_group}" \
  --stdev_coeff "${stdev_coeff}" \
  ${ref_arg} \
  ${ref_fasta_arg} \
  --output_dir ./ \
  --debug_log_path ./ \
  "${input_path}"

# Check if immunotyper-SR ran successfully
if [ $? -eq 0 ]; then
  echo "Parameter sweep run completed successfully"
  echo "Parameters: lg${landmark_groups}_lpg${landmarks_per_group}_stdev${stdev_coeff}"
  echo "Output saved to: ${full_output_dir}"
else
  echo "Error: immunotyper-SR failed"
  echo "Parameters: lg${landmark_groups}_lpg${landmarks_per_group}_stdev${stdev_coeff}"
fi