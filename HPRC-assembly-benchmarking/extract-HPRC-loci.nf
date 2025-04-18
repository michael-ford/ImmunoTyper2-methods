#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Global parameters
params.download_grch38 = params.download_grch38 ?: false
params.grch38_url    = 'http://example.com/GRCh38_full_analysis_set_plus_decoy_hla.fa'
params.bed           = './loci-GRCh38-coordinates.bed'
params.assembly_dir  = './assembly-sequences'
params.use_index_cache = params.use_index_cache ?: false

// New parameters for digger
params.digger_reference_dir       = "allele-databases"


// Process 1: Download GRCh38_full_analysis_set_plus_decoy_hla.fa
process downloadGRCh38 {
    publishDir "../", mode: 'move'
    output:
        path "GRCh38_full_analysis_set_plus_decoy_hla.fa", emit: grch38_fa
    script:
    """
    wget -O GRCh38_full_analysis_set_plus_decoy_hla.fa ${params.grch38_url}
    """
}

// Process 2: Extract loci sequences from GRCh38 using bedtools.
// Each locus is saved as its own FASTA file named {LOCUS}.fasta in loci-GRCh38-sequences/
process extractRegions {
    publishDir "loci-GRCh38-sequences", mode: 'rellink'
    input:
        tuple path(grch38_fa), path(bed)
    output:
        path "*.fasta", emit: loci_fasta
    script:
    """
    # Extract all loci into one FASTA using the BED file (-name uses the 4th column as header)
    bedtools getfasta -fi ${grch38_fa} -bed ${bed} -name -fo all_loci.fasta

    # Split the combined FASTA into individual files.
    # Sanitize header names: replace ':' with '_' then split on "__" and use the first field as the locus name.
    awk '/^>/{ 
             if(f){ close(f) }
             header = substr(\$0, 2)
             gsub(/[:]/, "_", header)
             split(header, a, "__")
             f = a[1] ".fasta"
             print \$0 > f
             next 
         } 
         { print \$0 > f }' all_loci.fasta

    rm all_loci.fasta
    """
}

// Process 5: Create Minimap2 index (.mmi) and FASTA index (.fai) for each HPRC assembly FASTA.
process indexAssemblies {
    publishDir "indexed-assemblies", mode: 'copy'
    input:
        path assembly
    output:
        // Output a tuple: (assembly, minimap index, faidx index)
        tuple path(assembly), path("${assembly.baseName}.mmi"), path("${assembly}.fai")
    script:
    """
    minimap2 -d ${assembly.baseName}.mmi ${assembly}
    samtools faidx ${assembly}
    """
}

// Process 3: For each combination of a locus FASTA and an indexed HPRC assembly,
// run minimap2 mapping using the .mmi index.
// (The FASTA index is not needed here, so we pass it along unchanged.)
process minimapMapping {
    maxForks 50
    input:
        tuple path(locus_fasta), path(assembly), path(assembly_index), path(assembly_faidx)
    output:
        // Emit a tuple: (locus name, assembly, assembly base name, mapping SAM file, assembly_faidx)
        tuple val(locus_fasta.baseName), path(assembly), env('ass_base'), path("*.sam"), path(assembly_faidx)
    script:
    """
    ass_base=\$(basename ${assembly} .fa.gz)
    minimap2 -ax asm5 ${assembly_index} ${locus_fasta} > mapping_${locus_fasta.baseName}_\${ass_base}.sam
    """
}

// Process 4: Extract the mapped regions from the minimap2 SAM file.
// For each chromosome with a mapping, determine the minimum start and maximum end positions,
// then extract the corresponding region from the HPRC assembly using the provided FASTA index.
process extractMappedRegion {
    publishDir "assembly-loci-sequences", mode: 'rellink'
    input:
        // Now includes the faidx index as the fifth element.
        tuple val(locus), path(assembly), val(assembly_prefix), path(mapping_sam), path(assembly_faidx)
    output:
        // Emit a tuple: (locus, assembly_prefix, single extracted FASTA file with multiple entries)
        tuple val(locus), val(assembly_prefix), path("${assembly_prefix}/${locus}/${assembly_prefix}_${locus}.extracted.fasta"), emit: extracted_fasta
    script:
    """
    # Create directory for this locus inside assembly-loci-sequences.
    mkdir -p ${assembly_prefix}/${locus}
        
    # Compute the minimum start and maximum end positions for each chromosome with a mapping.
    samtools view -F 4 ${mapping_sam} | awk '{
      chr = \$3; start = \$4; end = \$4 + length(\$10);
      if (!(chr in min) || start < min[chr]) { min[chr] = start }
      if (!(chr in max) || end > max[chr]) { max[chr] = end }
    }
    END {
      for (c in min) {
         print c, min[c], max[c]
      }
    }' > region_bounds.txt

    # Define the output file (single FASTA with multiple entries)
    outfile=${assembly_prefix}/${locus}/${assembly_prefix}_${locus}.extracted.fasta
    > \$outfile   # Create/empty the file

    # Loop over each chromosome region and append the extracted region to the output file.
    while read chr start end; do
      echo "Extracting region for \$chr: \$start-\$end" 1>&2
      samtools faidx ${assembly} \$chr:\$start-\$end >> \$outfile
    done < region_bounds.txt

    rm region_bounds.txt
    """
}

process runDigger {
    publishDir "digger-results", mode: 'rellink'
    errorStrategy 'ignore'
    input:
        // Tuple from previous process: (locus, assembly_prefix, extracted_fasta)
        tuple val(locus), val(assembly_prefix), path(extracted_fasta)
        // Allele database inputs provided as files so they're staged in the work directory
        path(reference_dir)        
    output:
        path "${assembly_prefix}/${locus}/*.csv", emit: digger_csv
    script:
    """
    prefix=\$(basename ${extracted_fasta} .fasta)
    mkdir -p ${assembly_prefix}/${locus} 
    digger ${extracted_fasta} \\
        -v_ref ${reference_dir}/${locus}/${locus}V.fasta \\
        -v_ref_gapped ${reference_dir}/${locus}/${locus}V-aligned.fasta \\
        -ref imgt,${reference_dir}/${locus}/${locus}V.fasta \\
        -species human \\
        -locus ${locus} \\
        ${assembly_prefix}/${locus}/\${prefix}.csv
    """
}


workflow {

    // Determine GRCh38 channel.
    if( params.download_grch38 ) {
        ch_grch38 = downloadGRCh38()
    }
    else {
        ch_grch38 = Channel.fromPath("../GRCh38_full_analysis_set_plus_decoy_hla.fa")
    }

    // Process 2: Extract regions from GRCh38 using the BED file.
    grch38_extract_ch = ch_grch38.map { grch38 -> tuple(grch38, file(params.bed)) }
    ch_loci_fasta = extractRegions(grch38_extract_ch)

    // Channel for HPRC assembly files (all .fa.gz files in the given directory)
    ch_assemblies = Channel.fromPath("${params.assembly_dir}/*fa.gz")

    // Process 5: Create indexes for each assembly.
    // Conditionally create the 'ch_indexed' channel:
    if ( params.use_index_cache ) {
        ch_indexed = Channel.fromPath("indexed-assemblies/*fa.gz").map { assembly ->
            def mmiName = assembly.getName().replaceFirst(/\.fa\.gz$/, ".fa.mmi")
            def faiName = assembly.getName() + ".fai"
            return tuple( assembly, file("indexed-assemblies/${mmiName}"), file("indexed-assemblies/${faiName}") )
        }
    } else {
        ch_indexed = indexAssemblies(ch_assemblies)
    }

    // Create the Cartesian product of each locus FASTA with each indexed assembly.
    // Since ch_indexed emits tuples of [assembly, assembly_index, assembly_faidx],
    // combine() them with each locus FASTA and then create a tuple:
    // [locus_fasta, assembly, assembly_index, assembly_faidx].
    ch_mapping_inputs = ch_loci_fasta.flatten()
                         .combine(ch_indexed)

    // Process 3: Run minimap2 mapping for each locus-assembly combination.
    ch_mapping = minimapMapping(ch_mapping_inputs)

    // Process 4: Extract the mapped region from the minimap output.
    ch_extracted = extractMappedRegion(ch_mapping)

    ch_reference_dir = Channel.value(file(params.digger_reference_dir))

    runDigger(ch_extracted, ch_reference_dir)


    
}
