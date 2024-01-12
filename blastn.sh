#!/bin/bash

# SLURM directives
#SBATCH -p gpu_p
#SBATCH -q gpu_short
#SBATCH --gres=gpu:1


# Directory containing consensus files
consensus_dir=$1
# Directory containing gold standard reference files
reference_dir=$2
# Output file
output_file="blastn_results.txt"

# Header for the output file
echo "ConsensusFile BLASTN_Output" > "$output_file"

# Loop over each consensus file
for consensus_file in "$consensus_dir"/*; do
    # Extract the name of the consensus file without the path
    consensus_name=$(basename "$consensus_file")

    # Determine the segment name from the consensus file name
    # Adjust the following line according to your file naming convention
    segment_name=$(echo "$consensus_name" | grep -oE '(_HA_|_MP_|_NA_|_NP_|_NS_|_PA_|_PB1_|_PB2_)')

    # Determine the corresponding gold standard reference file
    # Adjust this line to match your reference file naming pattern
    segment_name_no_underscore=$(echo "$segment_name" | sed 's/_//g')
    reference_file="${reference_dir}/${segment_name_no_underscore}.fasta" 

    # Run BLASTN
    # Adjust the BLASTN command according to your requirements
    blastn_output=$(blastn -query "$consensus_file" -subject "$reference_file" -outfmt 6)

    # Append the result to the output file
    echo "${consensus_name} ${blastn_output}" >> "$output_file"
done


