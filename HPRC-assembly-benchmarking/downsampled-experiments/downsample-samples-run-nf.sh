#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=norm
#SBATCH --time=10:00:00

nextflow run downsample-samples.nf -profile biowulf -resume