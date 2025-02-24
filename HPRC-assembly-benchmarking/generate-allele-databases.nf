#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define the input pattern for unaligned sequence files
params.input = "../ImmunoTyper2/immunotyper/data/allele_databases/*/*-IMGT-allele-db-no_duplicates.fa"


process run_vquest {
    input:
      tuple val(gene), path(unaligned_seq)
    output:
      // Emit tuple: [gene, vquest_airr.tsv, original unaligned_seq] for downstream processes
      tuple val(gene), path("vquest_airr.tsv"), path(unaligned_seq), emit: vquest_output
    script:
      """
      # Compute LOCUS as the first 3 characters of gene (e.g. IGHV -> IGH)
      GENE="${gene}"
      LOCUS=\${GENE:0:3}
      vquest --species human --receptorOrLocusType \${LOCUS} --fileSequences ${unaligned_seq}
      """
}

process parse_vquest {
    // Publish output FASTA to a gene-specific directory under allele-databases using model 'rellink'
    publishDir "./allele-databases", mode: 'rellink'
    input:
      tuple val(gene), path(vquest_airr)
    output:
      // The output file is written to a subdirectory named after the gene type.
      tuple val(gene), path("${gene}/${gene}-aligned.fasta"), emit: aligned_fasta
    script:
      """
      mkdir -p ${gene}
      # Parse the TSV: skip header, extract the 2nd '|' delimited field from the first column as header,
      # and use column 14 (sequence_alignment) as the sequence, converting it to uppercase.
      awk 'BEGIN{FS="\\t"} NR>1 { split(\$1, a, "|"); header=">" a[2]; seq=toupper(\$14); print header; print seq }' ${vquest_airr} > ${gene}/${gene}-aligned.fasta
      """
}

process parse_unaligned {
    // Publish output FASTA to a gene-specific directory under allele-databases using model 'rellink'
    publishDir "./allele-databases", mode: 'rellink'
    input:
      tuple val(gene), path(unaligned_seq)
    output:
      tuple val(gene), path("${gene}/${gene}.fasta"), emit: parsed_fasta
    script:
      """
      mkdir -p ${gene}
      # Parse the original FASTA: for each header, extract the 2nd '|' delimited field,
      # and convert all sequence lines to uppercase.
      awk '/^>/ { split(\$0, a, "|"); print ">" a[2] } !/^>/ { print toupper(\$0) }' ${unaligned_seq} > ${gene}/${gene}.fasta
      """
}


workflow {

    // Create a channel that emits a tuple: [gene_type, unaligned_seq file]
    Channel
        .fromPath(params.input)
        .map { file ->
            def gene = file.parent.getName()
            tuple(gene, file)
        }
        .set { unaligned_files_ch }

    // Run vquest for each input file
    run_vquest_out = unaligned_files_ch | run_vquest

    // Split the output tuple so that we can feed the vquest_airr.tsv and the original fasta separately.
    parse_vquest_in = run_vquest_out.map { gene, vquest_airr, unaligned_seq -> tuple(gene, vquest_airr) }
    parse_unaligned_in = run_vquest_out.map { gene, vquest_airr, unaligned_seq -> tuple(gene, unaligned_seq) }

    // Run parsing processes
    aligned_results = parse_vquest_in | parse_vquest
    parsed_results = parse_unaligned_in | parse_unaligned
}
