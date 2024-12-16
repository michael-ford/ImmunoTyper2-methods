// File: 1KGP_trios/workflow.nf

nextflow.enable.dsl=2

// Define default parameters
params.bam_reference_type = "hg38"  // Default reference type
params.wgs_dir = "${projectDir}/wgs-samples"
params.reference_fasta = "${projectDir}/GRCh38_full_analysis_set_plus_decoy_hla.fa"
params.output_dir = "immunotyper-output"

// List of gene types to process
gene_types = ['ighv', 'igkv', 'iglv', 'trav', 'trbv', 'trdv', 'trgv']

include { Immunotyper } from "./immunotyper.nf"

workflow {
    // Create channel of input BAM files and their indices
    bam_files = Channel
        .fromPath("${params.wgs_dir}/*.cram")
        .map { bam ->
            def bai = file("${bam}.crai")
            if (!bai.exists()) {
                log.error "BAI file not found for BAM: ${bam}"
                return null
            }
            return tuple(bam, bai)
        }
        .filter { it != null }  // Remove any null entries from failed samples

    // Create a cross product of BAM files and gene types
    process_inputs = bam_files
        .combine(Channel.from(gene_types))
        .map { bam, bai, gene_type -> 
            tuple(bam, bai, gene_type)
        }

    // Process each combination
    Immunotyper(process_inputs)

}