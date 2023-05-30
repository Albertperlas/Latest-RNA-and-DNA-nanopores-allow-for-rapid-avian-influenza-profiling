#!/bin/bash

#slurm command

#SBATCH --ntasks=1
#SBATCH -p gpu_p
#SBATCH -q gpu
#SBATCH --gres=gpu:1

# Read command-line arguments
input_fast5_dir=$1
output_fastq_file=$2
reference_genome_file=$3
output_sam_file=$4

# Base calling with Guppy
./ont-guppy/bin/guppy_basecaller \
    -i "$input_fast5_dir" \
    -r -s "$output_fastq_file" \
    -c rna_r9.4.1_70bps_hac.cfg \
    -x "cuda:0"



