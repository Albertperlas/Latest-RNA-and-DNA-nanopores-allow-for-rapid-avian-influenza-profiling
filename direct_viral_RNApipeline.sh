#!/bin/bash

#slurm command

srun -p gpu_p -q gpu --gres=gpu:1 

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


minimap2 -ax splice -uf -k14 ref.fa direct-rna.fq > aln.sam


# Mapping with Minimap2 from direct RNA-seq data
minimap2 \
    -ax splice \
    -uf -k14 \
    "$reference_genome_file" \
    "$output_fastq_file" \
    -o "$output_sam_file"
