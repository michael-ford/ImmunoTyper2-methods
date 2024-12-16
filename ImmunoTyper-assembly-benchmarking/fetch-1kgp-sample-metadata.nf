nextflow.enable.dsl=2

params.sra_run_table1 = "../SraRunTable-ERP120144.tsv"
params.sra_run_table2 = "../SraRunTable-ERP114329.tsv"
params.iglv_samples = "IGLV-1kgp-samples.txt"
params.trav_samples = "TRAV-1kgp-samples.txt"

process mapSraAccession {
    publishDir "./", mode: 'rellink'

    input:
    tuple path(sample_list), path(sra_run_table1), path(sra_run_table2), val(sample_type)

    output:
    tuple path("${sample_type}_sample_to_sra_accession.txt"), val(sample_type), emit: sample_to_sra

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
    done < "${sample_list}" > "${sample_type}_sample_to_sra_accession.txt"
    """
}

process generateURLs {
    publishDir "./", mode: 'rellink'

    input:
    tuple path(sample_to_sra_file), val(sample_type)

    output:
    tuple path("${sample_type}_download_urls.txt"), val(sample_type), emit: urls

    script:
    """
    while IFS=',' read -r sample sra_accession; do
        # Extract the first 6 characters of the SRA accession for subfolder structure
        sra_prefix=\${sra_accession:0:6}

        # Print the two required URLs for each sample
        echo "vol1/run/\${sra_prefix}/\${sra_accession}/\${sample}.final.cram"
        echo "vol1/run/\${sra_prefix}/\${sra_accession}/\${sample}.final.cram.crai"
    done < "${sample_to_sra_file}" > "${sample_type}_download_urls.txt"
    """
}

process generateGlobusURLs {
    publishDir "./", mode: 'rellink'

    input:
    tuple path(trav_urls), path(iglv_urls)

    output:
    path "globus_urls.txt", emit: globus_urls

    script:
    """
    # Initialize empty output file
    > globus_urls.txt
    
    # Process TRAV URLs
    while IFS= read -r line; do
        x=\${line##*/}
        echo "\$line data/1kgp-wgs/TRAV/\$x" >> globus_urls.txt
    done < ${trav_urls}
    
    # Process IGLV URLs
    while IFS= read -r line; do
        x=\${line##*/}
        echo "\$line data/1kgp-wgs/IGLV/\$x" >> globus_urls.txt
    done < ${iglv_urls}
    """
}

workflow {
    // Create channels for both sample types
    iglv_channel = Channel.fromPath(params.iglv_samples)
        .map { sample_list -> tuple(sample_list, file(params.sra_run_table1), file(params.sra_run_table2), 'IGLV') }
    
    trav_channel = Channel.fromPath(params.trav_samples)
        .map { sample_list -> tuple(sample_list, file(params.sra_run_table1), file(params.sra_run_table2), 'TRAV') }
    
    // Combine both channels
    sample_mapping_channel = iglv_channel.mix(trav_channel)
    
    // Map 1KGP sample IDs to SRA accession IDs
    sample_to_sra_files = mapSraAccession(sample_mapping_channel)
    
    // Generate URLs for downloading
    urls = generateURLs(sample_to_sra_files)
    
    // Split and group URLs by type
    trav_urls = urls.urls
        .filter { it[1] == 'TRAV' }
        .map { it[0] }
    
    iglv_urls = urls.urls
        .filter { it[1] == 'IGLV' }
        .map { it[0] }
    
    // Combine URL files into a tuple
    url_files = trav_urls
        .combine(iglv_urls)
    
    // Generate Globus URLs from all URL files
    globus_urls = generateGlobusURLs(url_files)
} 