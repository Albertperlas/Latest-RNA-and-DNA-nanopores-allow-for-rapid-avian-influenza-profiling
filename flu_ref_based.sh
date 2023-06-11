#!/bin/bash
#SBATCH -p gpu_p
#SBATCH -q gpu
#SBATCH --gres=gpu:1


input_file="$1"
reference_file="$2"
output_file="filtered.fastq"
min_length=50
gold_standard="$3"

# Remove short reads
seqkit seq -m "$min_length" "$input_file" > "$output_file"

# Align to reference
minimap2 -ax map-ont "$reference_file" "$output_file" > aln.sam

# Convert to BAM, index, and sort
samtools view -bS aln.sam | samtools sort -o sorted.bam - && samtools index sorted.bam

# Select the best reference 

samtools idxstats sorted.bam > idxstats.txt

best_reference="$(awk '{if ($3 > max) {max=$3; ref=$1}} END {print ref}' idxstats.txt)"

# Copy bam file and fasta file from the best reference 

samtools view -b -o best_reference.bam sorted.bam "$best_reference"
samtools faidx "$reference_file" "$best_reference" > best_reference.fasta 

#sort and index bam file from best reference 

samtools sort -o sorted_best_reference.bam best_reference.bam && samtools index sorted_best_reference.bam

#generate consensus with bcftools

bcftools mpileup -Ou -f best_reference.fasta sorted_best_reference.bam | bcftools call -Oz -mv -o variants.vcf.gz
bcftools index variants.vcf.gz
bcftools consensus -f best_reference.fasta variants.vcf.gz > consensus_bcftools.fasta

#generate consensus with ivar 

samtools mpileup -aa -A -d 0 -Q 0  sorted_best_reference.bam | ivar consensus -m 0 -q 0 -p consensus_ivar.fa  



#check results with blastn againts gold standard
blastn -query consensus_bcftools.fasta -subject "$gold_standard" -out output_summary.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

blastn -query consensus_ivar.fa -subject "$gold_standard" -out output_summary.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

