#!/bin/bash

# SLURM directives
#SBATCH -p gpu_p
#SBATCH -q gpu_short
#SBATCH --gres=gpu:1

# Read command-line arguments for input directory and output file
input_dir=$1
output_fastq_file=$2

# Dorado basecalling command 
/home/haicu/albert.perlas/dorado-0.4.3-linux-x64/bin/dorado basecaller \
    --emit-fastq \
    /home/haicu/albert.perlas/dorado-0.4.3-linux-x64/bin/rna004_130bps_hac@v3.0.1 \
    "$input_dir" \
    > "$output_fastq_file"
