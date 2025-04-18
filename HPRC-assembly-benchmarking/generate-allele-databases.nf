#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define the input pattern for unaligned sequence files
params.input = "../ImmunoTyper2/immunotyper/data/allele_databases/*/*-IMGT-allele-db-no_duplicates.fa"

process parse_input_fasta {
    input:
      tuple val(gene), path(unaligned_seq)
    output:
      tuple val(gene), path("input_parsed.fasta"), emit: parsed_input
    script:
      """
      awk '/^>/ { split(\$0, a, "|"); print ">" a[2] } !/^>/ { print toupper(\$0) }' ${unaligned_seq} > input_parsed.fasta
      """
}

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
    input:
      tuple val(gene), path(vquest_airr), path(unaligned_seq)
    output:
      tuple val(gene), path("success_aligned.fasta"), path("success_ids.txt"), path("vquest_failed_ids.txt"), path(unaligned_seq), emit: parsed_output
    script:
      """
      touch vquest_failed_ids.txt
      awk 'BEGIN { FS="\\t" }
           NR > 1 {
             id = \$1;
             seq = \$14;
             if (seq != "") {
               print ">" id >> "success_aligned.fasta";
               print seq >> "success_aligned.fasta";
               print id >> "success_ids.txt";
             } else {
               print id >> "vquest_failed_ids.txt";
             }
           }' ${vquest_airr}
      """
}

process filter_unaligned {
    publishDir "./allele-databases", mode: 'rellink'
    input:
      tuple val(gene), path(success_aligned), path(success_ids), path(unaligned_seq)
    output:
      tuple val(gene), path("*/${gene}.fasta"), path(success_aligned), emit: filtered_output
    script:
      """
      GENE="${gene}"
      LOCUS=\${GENE:0:3}
      mkdir -p \${LOCUS}
            
      # Use grep with -A1 and --no-group-separator to include the header and following sequence line,
      # without printing the delimiter lines between groups.
      grep -A1 --no-group-separator -F -f ${success_ids} ${unaligned_seq} > \$LOCUS/${gene}.fasta
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

    // Parse the input fasta files to remove the IMGT headers
    parsed_input = parse_input_fasta(unaligned_files_ch)

    // Run vquest for each input file
    run_vquest_out = run_vquest(parsed_input)
    
    // Parse VQUEST output to extract only successful alignments
    parsed_vquest = parse_vquest(run_vquest_out)
    parsed_vquest_failed_removed = parsed_vquest
        .map { gene, success_aligned, success_ids, failed_ids, unaligned_seq ->
        tuple(gene, success_aligned, success_ids, unaligned_seq)}

    // Filter the unaligned sequences to keep only those successfully aligned by vquest
    filtered_unaligned = filter_unaligned(parsed_vquest_failed_removed)
}