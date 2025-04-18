#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Global parameters
params.cram_dir = "../wgs-samples/"  // Required parameter: directory containing CRAM files
params.output_dir = './'  // Base directory for outputs
params.threads = 4  // Number of threads to use for samtools processes
params.target_chr = "chr1"  // Chromosome to use for depth calculations
params.reference_genome = "../../GRCh38_full_analysis_set_plus_decoy_hla.fa"  // Reference genome (required for CRAM processing)

// Process 1: Downsample CRAM files by 1/3 (for 10x coverage)
process downsampleOneThird {
    publishDir "${params.output_dir}/10x-samples", mode: 'symlink'
    tag "${sample_id}"
    
    input:
        tuple val(sample_id), path(cram), path(crai)
        path reference_genome
    
    output:
        tuple val(sample_id), path("${sample_id}.10x.cram"), path("${sample_id}.10x.cram.crai"), emit: downsampled_cram
        
    script:
    """
    # Downsample the CRAM file to 1/3 of original
    samtools view -@ ${params.threads} -h -s 0.33 -o ${sample_id}.10x.cram --output-fmt CRAM --reference ${reference_genome} ${cram}
    
    # Index the downsampled CRAM file
    samtools index -@ ${params.threads} ${sample_id}.10x.cram
    """
}

// Process 2: Downsample CRAM files by 2/3 (for 20x coverage)
process downsampleTwoThirds {
    publishDir "${params.output_dir}/20x-samples", mode: 'symlink'
    tag "${sample_id}"
    
    input:
        tuple val(sample_id), path(cram), path(crai)
        path reference_genome
    
    output:
        tuple val(sample_id), path("${sample_id}.20x.cram"), path("${sample_id}.20x.cram.crai"), emit: downsampled_cram
        
    script:
    """
    # Downsample the CRAM file to 2/3 of original
    samtools view -@ ${params.threads} -h -s 0.67 -o ${sample_id}.20x.cram --output-fmt CRAM --reference ${reference_genome} ${cram}
    
    # Index the downsampled CRAM file
    samtools index -@ ${params.threads} ${sample_id}.20x.cram
    """
}

// Process 3: Calculate depth for 10x downsampled samples on a single chromosome
process calculateDepth10x {
    tag "${sample_id}"
    
    input:
        tuple val(sample_id), path(cram), path(crai)
        path reference_genome
    
    output:
        tuple val(sample_id), path("${sample_id}.10x.depth.txt"), emit: depth_file
        
    script:
    """
    # Calculate depth for a single chromosome only
    samtools depth -a -r ${params.target_chr} --reference ${reference_genome} ${cram} > ${sample_id}.10x.depth.txt
    """
}

// Process 4: Calculate depth for 20x downsampled samples on a single chromosome
process calculateDepth20x {
    tag "${sample_id}"
    
    input:
        tuple val(sample_id), path(cram), path(crai)
        path reference_genome
    
    output:
        tuple val(sample_id), path("${sample_id}.20x.depth.txt"), emit: depth_file
        
    script:
    """
    # Calculate depth for a single chromosome only
    samtools depth -a -r ${params.target_chr} --reference ${reference_genome} ${cram} > ${sample_id}.20x.depth.txt
    """
}

// Process 5: Generate depth statistics summary for 10x samples
process generateDepthStats10x {
    publishDir "${params.output_dir}/stats", mode: 'symlink'
    
    input:
        path(depth_files)
    
    output:
        path "10x_depth_stats.txt"
        
    script:
    """
    # Initialize the statistics file
    echo "Sample_ID\tMean_Depth\tMedian_Depth\t${params.target_chr}" > 10x_depth_stats.txt
    
    # Process each depth file
    for depth_file in ${depth_files}; do
        # Extract sample ID from filename
        sample_id=\$(basename \$depth_file .10x.depth.txt)
        
        # Calculate mean depth
        mean_depth=\$(awk '{ sum += \$3 } END { print sum/NR }' \$depth_file)
        
        # Calculate median depth
        median_depth=\$(awk '{ depths[NR] = \$3 } END { 
            asort(depths); 
            if (NR % 2) { print depths[int(NR/2) + 1] } 
            else { print (depths[NR/2] + depths[NR/2 + 1]) / 2.0 } 
        }' \$depth_file)
        
        # Append to the statistics file
        echo -e "\$sample_id\t\$mean_depth\t\$median_depth" >> 10x_depth_stats.txt
    done
    
    # Calculate group statistics
    mean_group=\$(awk 'NR>1 {sum+=\$2} END {print sum/(NR-1)}' 10x_depth_stats.txt)
    median_group=\$(awk 'NR>1 {a[NR-1]=\$2} END {
        asort(a);
        n=length(a);
        if (n % 2) { print a[int(n/2) + 1] }
        else { print (a[n/2] + a[n/2 + 1]) / 2.0 }
    }' 10x_depth_stats.txt)
    
    # Append group statistics
    echo -e "\nGroup Statistics:" >> 10x_depth_stats.txt
    echo -e "Mean Depth: \$mean_group" >> 10x_depth_stats.txt
    echo -e "Median Depth: \$median_group" >> 10x_depth_stats.txt
    """
}

// Process 6: Generate depth statistics summary for 20x samples
process generateDepthStats20x {
    publishDir "${params.output_dir}/stats", mode: 'symlink'
    
    input:
        path(depth_files)
    
    output:
        path "20x_depth_stats.txt"
        
    script:
    """
    # Initialize the statistics file
    echo "Sample_ID\tMean_Depth\tMedian_Depth\t${params.target_chr}" > 20x_depth_stats.txt
    
    # Process each depth file
    for depth_file in ${depth_files}; do
        # Extract sample ID from filename
        sample_id=\$(basename \$depth_file .20x.depth.txt)
        
        # Calculate mean depth
        mean_depth=\$(awk '{ sum += \$3 } END { print sum/NR }' \$depth_file)
        
        # Calculate median depth
        median_depth=\$(awk '{ depths[NR] = \$3 } END { 
            asort(depths); 
            if (NR % 2) { print depths[int(NR/2) + 1] } 
            else { print (depths[NR/2] + depths[NR/2 + 1]) / 2.0 } 
        }' \$depth_file)
        
        # Append to the statistics file
        echo -e "\$sample_id\t\$mean_depth\t\$median_depth" >> 20x_depth_stats.txt
    done
    
    # Calculate group statistics
    mean_group=\$(awk 'NR>1 {sum+=\$2} END {print sum/(NR-1)}' 20x_depth_stats.txt)
    median_group=\$(awk 'NR>1 {a[NR-1]=\$2} END {
        asort(a);
        n=length(a);
        if (n % 2) { print a[int(n/2) + 1] }
        else { print (a[n/2] + a[n/2 + 1]) / 2.0 }
    }' 20x_depth_stats.txt)
    
    # Append group statistics
    echo -e "\nGroup Statistics:" >> 20x_depth_stats.txt
    echo -e "Mean Depth: \$mean_group" >> 20x_depth_stats.txt
    echo -e "Median Depth: \$median_group" >> 20x_depth_stats.txt
    """
}

workflow {
    // Validate required parameters within the workflow
    if (params.cram_dir == null) {
        error "ERROR: Please specify the directory containing CRAM files using --cram_dir parameter"
    }

    if (params.reference_genome == null) {
        error "ERROR: Please specify a reference genome using --reference_genome parameter (required for CRAM processing)"
    }

    // Channel for CRAM files and their indexes
    ch_crams = Channel.fromFilePairs("${params.cram_dir}/*.{cram,cram.crai}", size: 2, flat: true)
        .map { file, cram, crai -> tuple(file, cram, crai) }
    
    // Reference genome channel
    ch_reference = Channel.value(file(params.reference_genome))
    
    // Downsample CRAM files
    ch_10x = downsampleOneThird(ch_crams, ch_reference)
    ch_20x = downsampleTwoThirds(ch_crams, ch_reference)
    
    // Calculate depth for downsampled files
    ch_depth_10x = calculateDepth10x(ch_10x.downsampled_cram, ch_reference)
    ch_depth_20x = calculateDepth20x(ch_20x.downsampled_cram, ch_reference)
    
    // Generate depth statistics summaries
    depth_files_10x = ch_depth_10x.map { _sample_id, depth_file -> depth_file }.collect()
    depth_files_20x = ch_depth_20x.map { _sample_id, depth_file -> depth_file }.collect()
    generateDepthStats10x(depth_files_10x)
    generateDepthStats20x(depth_files_20x)
}