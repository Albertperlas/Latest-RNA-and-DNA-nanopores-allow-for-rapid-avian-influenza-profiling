#!/bin/bash
#SBATCH -p gpu_p
#SBATCH -q gpu
#SBATCH --gres=gpu:1

#you can prepare your reference with makeblastdb 
#give reference_file without .fasta extension, blastn needs the indexed files

input_file="$1"
reference_file="$2"
output_file="filtered.fastq"
min_length=50

# Remove short reads
seqkit seq -m "$min_length" "$input_file" > "$output_file"


#we use flye for the assembly 

flye --nano-raw "$output_file" --out-dir . 

#polishing with minimap2 and racon

#minimap2 with assembly.fasta as input and fastq 

minimap2 -ax map-ont -t 8 assembly.fasta "$output_file" > aln_ass.sam

#racon consensus contruction

racon "$output_file" aln_ass.sam assembly.fasta > polyshed_consensus.fasta

# Output directory for blastn results
blastn_output_dir="blastn_results"

# Create output directory if it doesn't exist
mkdir -p "$blastn_output_dir"

# Step 1: Separate contigs of polyshed_consensus.fasta into individual FASTA files
awk '/^>/{s="polyshed_consensus_"++d".fasta"} {print > s}' polyshed_consensus.fasta

# Step 2: Perform blastn against the reference database for each contig
for fasta_file in polyshed_consensus_*.fasta; do
    # Get the contig name without the file extension
    contig_name=$(basename "$fasta_file" .fasta)
    
    # Perform blastn
    blastn -query "$fasta_file" -db "$reference_file" -out blastn_results/"$contig_name.txt" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads 4
done
    

