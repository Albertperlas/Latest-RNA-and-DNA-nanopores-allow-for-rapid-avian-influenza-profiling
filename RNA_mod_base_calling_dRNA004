#!/bin/bash

# SLURM directives
#SBATCH --ntasks=1
#SBATCH -p gpu_p
#SBATCH -q gpu_short
#SBATCH --gres=gpu:1

# Read command-line arguments
input_dir=$1
output_bam_file=$2

# Dorado basecalling and call m6a modifications
/home/haicu/albert.perlas/dorado-0.4.3-linux-x64/bin/dorado basecaller \
    --modified-bases-models /home/haicu/albert.perlas/dorado-0.4.3-linux-x64/bin/rna004_130bps_sup@v3.0.1_m6A_DRACH@v1/ \
    --reference /home/haicu/albert.perlas/past/tfm/downsampled_fastq/plot_coverage/italy_ref_plot.fasta \
    /home/haicu/albert.perlas/dorado-0.4.3-linux-x64/bin/rna004_130bps_sup@v3.0.1/ \
    "$input_dir" \
    > "$output_bam_file"
