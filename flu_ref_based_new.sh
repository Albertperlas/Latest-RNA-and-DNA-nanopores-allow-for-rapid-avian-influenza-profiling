#!/bin/bash
#SBATCH -p gpu_p
#SBATCH -q gpu_short
#SBATCH --gres=gpu:1

# Directory containing FASTQ files
fastq_dir="/path/to/fastq/files"
# Directory containing reference files for each segment
reference_dir="/path/to/reference/files"

min_length=50

# Loop over each FASTQ file
for fastq_file in "$fastq_dir"/*.fastq; do
    # Extract the name of the FASTQ file without extension
    fastq_name=$(basename "$fastq_file" .fastq)

    # Loop over each influenza segment reference
    for reference_file in "$reference_dir"/*.fasta; do
        # Extract the segment name from the reference file
        segment_name=$(basename "$reference_file" .fasta)

        # Define output file names
        output_file="${fastq_name}_${segment_name}_filtered.fastq"
        consensus_bcftools_file="${fastq_name}_${segment_name}_consensus_bcftools.fasta"
        consensus_ivar_file="${fastq_name}_${segment_name}_consensus_ivar.fa"

        # Remove short reads
        seqkit seq -m "$min_length" "$fastq_file" > "$output_file"

        # Initial alignment to reference database
        minimap2 -ax map-ont "$reference_file" "$output_file" > "${fastq_name}_${segment_name}_aln.sam"

        # Convert to BAM, index, and sort
        samtools view -bS "${fastq_name}_${segment_name}_aln.sam" | samtools sort -o "${fastq_name}_${segment_name}_sorted.bam" - && samtools index "${fastq_name}_${segment_name}_sorted.bam"

        # Select the best reference 
        samtools idxstats "${fastq_name}_${segment_name}_sorted.bam" > "${fastq_name}_${segment_name}_idxstats.txt"
        best_reference="$(awk '{if ($3 > max) {max=$3; ref=$1}} END {print ref}' "${fastq_name}_${segment_name}_idxstats.txt")"

        # Copy bam file and fasta file from the best reference 
        samtools view -b -o "${fastq_name}_${segment_name}_best_reference.bam" "${fastq_name}_${segment_name}_sorted.bam" "$best_reference"
        samtools faidx "$reference_file" "$best_reference" > "${fastq_name}_${segment_name}_best_reference.fasta" 

        # Align all reads again to the best reference
        minimap2 -ax map-ont "${fastq_name}_${segment_name}_best_reference.fasta" "$output_file" > "${fastq_name}_${segment_name}_aln_best_ref.sam"

        # Convert to BAM, index, and sort for best reference
        samtools view -bS "${fastq_name}_${segment_name}_aln_best_ref.sam" | samtools sort -o "${fastq_name}_${segment_name}_sorted_best_reference.bam" - && samtools index "${fastq_name}_${segment_name}_sorted_best_reference.bam"

        # Generate consensus with bcftools
        bcftools mpileup -Ou -f "${fastq_name}_${segment_name}_best_reference.fasta" "${fastq_name}_${segment_name}_sorted_best_reference.bam" | bcftools call -Oz -mv -o "${fastq_name}_${segment_name}_variants.vcf.gz"
        bcftools index "${fastq_name}_${segment_name}_variants.vcf.gz"
        bcftools consensus -f "${fastq_name}_${segment_name}_best_reference.fasta" "${fastq_name}_${segment_name}_variants.vcf.gz" > "$consensus_bcftools_file"

        # Generate consensus with ivar 
        samtools mpileup -aa -A -d 0 -Q 0 "${fastq_name}_${segment_name}_sorted_best_reference.bam" | ivar consensus -m 0 -q 0 -p "$consensus_ivar_file"

    done
done
