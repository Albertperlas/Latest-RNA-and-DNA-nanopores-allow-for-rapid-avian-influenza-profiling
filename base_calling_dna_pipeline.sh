#!/bin/bash

#slurm command

#SBATCH --ntasks=1
#SBATCH -p gpu_p
#SBATCH -q gpu
#SBATCH --gres=gpu:1

# Read command-line arguments
input_fast5_dir=$1
output_fastq_file=$2

# HAC Base calling with Guppy
./ont-guppy/bin/guppy_basecaller \
    -i "$input_fast5_dir" \
    -r -s "$output_fastq_file" \
    --detect_barcodes \
    -c dna_r10.4.1_e8.2_400bps_hac.cfg \
    -x "cuda:0"
