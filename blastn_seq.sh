#!/bin/bash

#blastn to all the segments and to all fasta files 


# Directory paths
dir1="./"
dir2="/home/haicu/albert.perlas/tfm/gold_standard"

# Sequence names
sequences=("PB2" "PB1" "PA" "NA" "MP" "HA" "NP" "NS")

# Iterate over sequences
for sequence in "${sequences[@]}"; do
    # Iterate over FASTA files in dir1
    for fasta_file1 in "$dir1"/*"${sequence}"*.fasta; do
        # Extract file name without extension and directory name
        file_name1=$(basename "$fasta_file1" .fasta)
        dir_name1=$(basename "$dir1")

        # Iterate over FASTA files in dir2
        for fasta_file2 in "$dir2"/*"${sequence}"*.fasta; do
            # Extract file name without extension and directory name
            file_name2=$(basename "$fasta_file2" .fasta)
            dir_name2=$(basename "$dir2")

            # Run blastn command
            blastn -query "$fasta_file1" -subject "$fasta_file2" -out "env_cDNA_${file_name1}.txt" -outfmt "6"
        done
    done
done

