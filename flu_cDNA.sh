#!/bin/bash
#SBATCH -p gpu_p
#SBATCH -q gpu
#SBATCH --gres=gpu:1

if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_file> <reference_file>"
    exit 1
fi

input_file="$1"
reference_file="$2"
output_file="filtered.fastq"
min_length=50

# Remove short reads
seqkit seq -m "$min_length" "$input_file" > "$output_file"

# Align to reference
minimap2 -ax map-ont "$reference_file" "$output_file" > aln.sam

# Convert to BAM, index, and sort
samtools view -bS aln.sam | samtools sort -o sorted.bam - && samtools index sorted.bam

# Generate a basic coverage report
samtools depth sorted.bam > coverage_report.txt
