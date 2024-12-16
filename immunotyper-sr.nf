// immunotyper-sr.nf
// Expects parameters:
// - params.gene_type: Gene type for the Immunotyper process
// - params.bam_reference: Either hg38 or hg37

nextflow.enable.dsl=2

process Immunotyper {
    tag "${input_file.baseName}"

    publishDir "./", mode: 'rellink'

    input:
    tuple path(input_file), path(index)

    output:
    path "${input_file.baseName}/${input_file.baseName}-${params.gene_type.toUpperCase()}_allele_calls.txt"
    path "${input_file.baseName}/${input_file.baseName}-${params.gene_type.toUpperCase()}_functional_allele_calls.txt"

    script:
    def is_cram = input_file.extension == 'cram'
    if (is_cram && !params.reference_fasta) {
        error "Reference FASTA is required for CRAM files. Please provide --reference_fasta parameter."
    }
    """
    input_path=\$(readlink -f "${input_file}")

    mkdir -p "${input_file.baseName}"
    cd "${input_file.baseName}"

    # Determine the reference type argument
    ref_arg=""
    if [[ "${params.bam_reference_type}" == "hg37" ]]; then
        ref_arg="--hg37"
    fi

    # Add reference FASTA argument if input is CRAM
    ref_fasta_arg=""
    if [[ "${is_cram}" == "true" ]]; then
        ref_fasta_arg="--ref ${params.reference_fasta}"
    fi

    immunotyper-SR --gene_type "${params.gene_type}" \$ref_arg \$ref_fasta_arg --output_dir ./ --debug_log_path ./ "../${input_file}"
    """
}
