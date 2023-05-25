#!/bin/bash

# Read command-line arguments
input_fast5_dir=$1
output_fastq_file=$2
reference_genome_file=$3
output_sam_file=$4

# Base calling with Guppy
guppy_basecaller \
    -i "$input_fast5_dir" \
    -s "$output_fastq_file" \
    --flowcell FLO-MIN106 \
    --kit SQK-LSK109 \
    --barcode_kits EXP-PBC096 \
    --qscore_filtering

# Mapping with Minimap2
minimap2 \
    -ax map-ont \
    -t 8 \
    "$reference_genome_file" \
    "$output_fastq_file" \
    -o "$output_sam_file"
