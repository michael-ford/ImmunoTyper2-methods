nextflow.enable.dsl=2

params.sra_run_table1 = "SraRunTable-ERP120144.tsv"
params.sra_run_table2 = "SraRunTable-ERP114329.tsv"
params.out_dir = "wgs-samples"

process downloadReference {
    publishDir "./", mode: 'rellink'

    output:
    path "GRCh38_full_analysis_set_plus_decoy_hla.fa", emit: reference

    script:
    """
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
    """
}

process fetchTrioList {
    publishDir "./", mode: 'rellink'

    output:
    path "trio_samples.txt", emit: trio_samples

    script:
    """
    curl -L -o trio_samples.txt 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt'
    """
}

process filterTrioSamples {
    publishDir "./", mode: 'rellink'

    input:
    path trio_samples

    output:
    path "filtered_trio_samples.txt", emit: filtered_trio_samples
    path "trios_metadata.tsv", emit: trios_metadata

    script:
    """
    # Extract lines with more than one HG ID and write them to trios_metadata.tsv
    awk '{hg_count=0; for (i=1; i<=NF; i++) if (\$i ~ /^(HG|NA)/) hg_count++; if (hg_count > 1) print \$0}' ${trio_samples} > trios_metadata.tsv

    # Convert each HG ID in these lines to a new line in filtered_trio_samples.txt
    awk '{for (i=1; i<=NF; i++) if (\$i ~ /^(HG|NA)/) print \$i}' trios_metadata.tsv | sort -u > filtered_trio_samples.txt
    """
}


process mapSraAccession {
    publishDir "./", mode: 'rellink'

    input:
    tuple path(filtered_trio_samples), path(sra_run_table1), path(sra_run_table2)

    output:
    path "1kgp_sample_to_sra_accession.txt", emit: sample_to_sra

    script:
    """
    # Extract 1KGP sample IDs and corresponding SRA run IDs
    while read sample; do
        # Search in sra_run_table1
        sra_entry=\$(grep "\$sample" ${sra_run_table1} | sort -k5,5nr | head -n 1)
        
        # If not found, search in sra_run_table2
        if [ -z "\$sra_entry" ]; then
            sra_entry=\$(grep "\$sample" ${sra_run_table2} | sort -k34,34nr | head -n 1)
        fi

        # If still not found, exit with an error
        if [ -z "\$sra_entry" ]; then
            echo "Error: SRA accession ID for sample \$sample not found in either table." >&2
            exit 1
        fi

        # Extract SRA run ID from the selected entry (assuming it's in the first column)
        sra=\$(echo "\$sra_entry" | cut -f1)

        echo "\$sample,\$sra"
    done < "${filtered_trio_samples}" > 1kgp_sample_to_sra_accession.txt
    """
}


process mapSraAccessionTest {
    input:
    tuple path(filtered_trio_samples), path(sra_run_table1), path(sra_run_table2)

    output:
    path "1kgp_sample_to_sra_accession.txt", emit: sample_to_sra

    script:
    """
    echo "HG00403,ERR3241665" > 1kgp_sample_to_sra_accession.txt
    """
}

process generateURLs {
    publishDir "./", mode: 'rellink'

    input:
    path sample_to_sra_file

    output:
    path "download_urls.txt", emit: urls

    script:
    """
    while IFS=',' read -r sample sra_accession; do
        # Extract the first 6 characters of the SRA accession for subfolder structure
        sra_prefix=\${sra_accession:0:6}

        # Print the two required URLs for each sample
        echo "vol1/run/\${sra_prefix}/\${sra_accession}/\${sample}.final.cram"
        echo "vol1/run/\${sra_prefix}/\${sra_accession}/\${sample}.final.cram.crai"
    done < "${sample_to_sra_file}" > download_urls.txt
    """
}

process generateGlobusURLs {
    publishDir "./", mode: 'rellink'

    input:
    path urls_file

    output:
    path "globus_urls.txt", emit: globus_urls

    script:
    """
    while IFS= read -r line; do
        x=\${line##*/}
        echo "\$line ${params.out_dir}/\$x"
    done < ${urls_file} > globus_urls.txt
    """
}


workflow {
    // Download reference genome first
    reference = downloadReference()

    // Fetch the list of 1KGP samples (not used in this test, but kept for context)
    trio_list = fetchTrioList()

    // Filter for trios and split the outputs
    filter_trio_results = filterTrioSamples(trio_list)

    // Extract the two outputs from the filterTrioSamples process
    filtered_trio_samples_out = filter_trio_results.filtered_trio_samples

    // Create a channel of tuples combining filtered_trio_samples, sra_run_table1, and sra_run_table2
    sample_mapping_channel = filtered_trio_samples_out
        .map { filtered_trio_samples -> tuple(filtered_trio_samples, file(params.sra_run_table1), file(params.sra_run_table2)) }

    // Map 1KGP sample IDs to SRA accession IDs using the provided SRA CSV files
    sample_to_sra_file = mapSraAccession(sample_mapping_channel)
    // sample_to_sra_file = mapSraAccessionTest(sample_mapping_channel)

    // Generate URLs for downloading
    urls = generateURLs(sample_to_sra_file)

    // Generate Globus URLs
    globus_urls = generateGlobusURLs(urls.urls)
}
