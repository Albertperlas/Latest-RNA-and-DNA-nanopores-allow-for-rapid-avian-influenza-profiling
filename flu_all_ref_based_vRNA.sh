#!/bin/bash
#SBATCH -p gpu_p
#SBATCH -q gpu
#SBATCH --gres=gpu:1

reference_file="$1"
fastq_directory="/home/haicu/albert.perlas/tfm/downsampled_fastq/final/fastq/RNA_files/"
gold_standard="$2"

min_length=50


# Get the base name of the reference file (without the extension)
reference_base=$(basename "$reference_file")
reference_base="${reference_base%.*}"


  # Loop over FASTQ files
  for fastq_file in "${fastq_directory}"/*.fastq; do
    # Get the base name of the FASTQ file (without the extension)
    fastq_base=$(basename "$fastq_file")
    fastq_base="${fastq_base%.*}"

    # Create a unique output directory for each combination of reference and FASTQ file
    output_directory="${fastq_base}_${reference_base}_output"
    mkdir -p "$output_directory"
    cd "$output_directory" || exit

    # Remove short reads
    seqkit seq -m "$min_length" "$fastq_file" > filtered.fastq

    # Align to reference
    minimap2 -ax splice -uf -k7 "$reference_file" filtered.fastq > "${fastq_base}_${reference_base}_aln.sam"

    # Convert to BAM, index, and sort
    samtools view -bS "${fastq_base}_${reference_base}_aln.sam" | samtools sort -o "${fastq_base}_${reference_base}_sorted.bam" - && samtools index "${fastq_base}_${reference_base}_sorted.bam"

    # Select the best reference
    samtools idxstats "${fastq_base}_${reference_base}_sorted.bam" > "${fastq_base}_${reference_base}_idxstats.txt"
    best_reference=$(awk '{if ($3 > max) {max=$3; ref=$1}} END {print ref}' "${fastq_base}_${reference_base}_idxstats.txt")

    # Copy BAM file and fasta file from the best reference
    samtools view -b -o "${fastq_base}_${reference_base}_best_reference.bam" "${fastq_base}_${reference_base}_sorted.bam" "$best_reference"
    samtools faidx "$reference_file" "$best_reference" > "${fastq_base}_${reference_base}_best_reference.fasta"

    # Sort and index BAM file from best reference
    samtools sort -o "${fastq_base}_${reference_base}_sorted_best_reference.bam" "${fastq_base}_${reference_base}_best_reference.bam" && samtools index "${fastq_base}_${reference_base}_sorted_best_reference.bam"

    # Generate consensus with bcftools
    bcftools mpileup -Ou -f "${fastq_base}_${reference_base}_best_reference.fasta" "${fastq_base}_${reference_base}_sorted_best_reference.bam" | bcftools call -Oz -mv -o "${fastq_base}_${reference_base}_variants.vcf.gz"
    bcftools index "${fastq_base}_${reference_base}_variants.vcf.gz"
    bcftools consensus -f "${fastq_base}_${reference_base}_best_reference.fasta" "${fastq_base}_${reference_base}_variants.vcf.gz" > "${fastq_base}_${reference_base}_consensus_bcftools.fasta"

    # Generate consensus with ivar
    samtools mpileup -aa -A -d 0 -Q 0 "${fastq_base}_${reference_base}_sorted_best_reference.bam" | ivar consensus -m 0 -q 0 -p "${fastq_base}_${reference_base}_consensus_ivar.fa"

    # Check results with blastn against gold standard
    blastn -query "${fastq_base}_${reference_base}_consensus_bcftools.fasta" -subject "$gold_standard" -out "${fastq_base}_${reference_base}_output_summary_bcftools.txt" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    blastn -query "${fastq_base}_${reference_base}_consensus_ivar.fa" -subject "$gold_standard" -out "${fastq_base}_${reference_base}_output_summary_ivar.txt" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

    # Remove intermediate files
    rm "${fastq_base}_${reference_base}_aln.sam"
    rm "${fastq_base}_${reference_base}_sorted.bam"
    rm "${fastq_base}_${reference_base}_sorted.bam.bai"
    rm "${fastq_base}_${reference_base}_best_reference.bam"
    rm "${fastq_base}_${reference_base}_best_reference.bam.bai"
    rm "${fastq_base}_${reference_base}_best_reference.fasta"
    rm "${fastq_base}_${reference_base}_best_reference.fasta.fai"
    rm "${fastq_base}_${reference_base}_consensus_ivar.qual.txt"
    rm "${fastq_base}_${reference_base}_idxstats.txt"
    rm "${fastq_base}_${reference_base}_sorted_best_reference.bam"
    rm "${fastq_base}_${reference_base}_sorted_best_reference.bam.bai"
    rm "filtered.fastq"
	
    cd ..
  done
done

