#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=norm
#SBATCH --time=10:00:00

nextflow run extract-HPRC-loci.nf -profile biowulf --use_index_cache -resume