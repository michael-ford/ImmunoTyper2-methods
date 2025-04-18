nextflow.enable.dsl=2

params.index_url             = 'https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/refs/heads/main/assembly_index/Year1_assemblies_v2_genbank.index'
params.index_publish_dir     = './'
params.assemblies_publish_dir = './assembly-sequences'

process DownloadIndex {
    publishDir params.index_publish_dir, mode: 'copy'
    output:
        path "Year1_assemblies_v2_genbank.index", emit: index_file
    script:
        """
        curl -L ${params.index_url} -o Year1_assemblies_v2_genbank.index
        """
}

process DownloadAssembly {
    publishDir params.assemblies_publish_dir, mode: 'rellink'
    input:
        val s3_url
    output:
        file("*")
    script:
        """
        fileName=\$(basename ${s3_url})
        aws s3 cp ${s3_url} \${fileName} --no-sign-request
        zcat \${fileName} | bgzip > \${fileName}.bgz
        mv \${fileName}.bgz \${fileName}
        """
}

workflow {
    // Download the index file
    index_file = DownloadIndex()

    // Parse the index file to extract S3 URLs (hap1 and hap2 FASTA)
    index_file
        .splitCsv(header: true, sep: '\t')
        .flatMap { row -> [ row.hap1_aws_fasta, row.hap2_aws_fasta ] }
        .filter { it && it.trim() }
        .set { s3_urls }

    // Download each assembly file from S3
    DownloadAssembly(s3_urls)
}
