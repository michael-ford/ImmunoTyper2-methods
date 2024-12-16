// File: workflow.nf

nextflow.enable.dsl=2

params.gene_type = 'iglv'   // Default parameter; should be provided at runtime
params.bam_reference_type = "hg38" // Default reference type; can be overridden at runtime
params.iglv_sample_list_dir = "${projectDir}/../assembly-annotations/IGLV-reference-annotations/"
params.trav_sample_list_dir = "${projectDir}/../assembly-annotations/TRAV-reference-annotations/"
params.iglv_bam_dir = "${projectDir}/../1kgp-wgs/IGLV"
params.trav_bam_dir = "${projectDir}/../1kgp-wgs/TRAV"
params.reference_fasta = "${projectDir}/../../../1KGP_Trios/GRCh38_full_analysis_set_plus_decoy_hla.fa"

include { Immunotyper } from "${projectDir}/../../../immunotyper-sr.nf"

workflow {
    // Get sample IDs based on gene type from the appropriate directory
    def sample_list_dir = params.gene_type == 'iglv' ? 
        params.iglv_sample_list_dir : 
        params.trav_sample_list_dir
    
    // Create channel of sample IDs from directory names
    sample_ids = Channel
        .fromPath("${sample_list_dir}/*", type: 'dir')
        .map { dir -> dir.getName() }
    // sample_ids.view()

    // Create channel of input files (BAM/CRAM) and their indices matching the sample IDs
    bam_files = sample_ids.map { sample_id ->
        if (params.gene_type == 'iglv') {
            def bam_file = file("${params.iglv_bam_dir}/${sample_id}-vs-GRCh38.bam")
            def bai_file = file("${params.iglv_bam_dir}/${sample_id}-vs-GRCh38.bam.bai")
            
            if (!bam_file.exists()) {
                log.error "BAM file not found for sample ${sample_id}: ${bam_file}"
                return null
            }
            if (!bai_file.exists()) {
                log.error "BAI file not found for sample ${sample_id}: ${bai_file}"
                return null
            }
            
            tuple(bam_file, bai_file)
        } else {
            def cram_file = file("${params.trav_bam_dir}/${sample_id}.final.cram")
            def crai_file = file("${params.trav_bam_dir}/${sample_id}.final.cram.crai")
            
            if (!cram_file.exists()) {
                log.error "CRAM file not found for sample ${sample_id}: ${cram_file}"
                return null
            }
            if (!crai_file.exists()) {
                log.error "CRAI file not found for sample ${sample_id}: ${crai_file}"
                return null
            }
            
            tuple(cram_file, crai_file)
        }
    }
    .filter { it != null }  // Remove any null entries from failed samples
    // .take(1)  // Add this line to take only the first BAM file
    // bam_files.view()
    
    Immunotyper(bam_files)
}
