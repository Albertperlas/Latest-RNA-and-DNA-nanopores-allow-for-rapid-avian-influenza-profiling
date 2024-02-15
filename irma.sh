#!/bin/bash
#SBATCH -p gpu_p
#SBATCH -q gpu_short
#SBATCH --gres=gpu:1


#/path/to/fastq_directory

input_dir="$1"

# Iterate over FASTQ files in the input directory
for fastq_file in "$input_dir"/*.fastq; do
    # Extract file name without extension
    file_name=$(basename "$fastq_file" .fastq)

    # Run IRMA FLU-minion command
    IRMA FLU-minion "$fastq_file" "${file_name}_MinION"
done

