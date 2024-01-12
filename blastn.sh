
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
    segment_name=$(echo "$consensus_name" | grep -oE '(HA|MP|NA|NP|NS|PA|PB1|PB2)')

    # Determine the corresponding gold standard reference file
    # Adjust this line to match your reference file naming pattern
    reference_file="${reference_dir}/${segment_name}_gold_standard.fasta"

    # Run BLASTN
    # Adjust the BLASTN command according to your requirements
    blastn_output=$(blastn -query "$consensus_file" -subject "$reference_file" -outfmt 6)

    # Append the result to the output file
    echo "${consensus_name} ${blastn_output}" >> "$output_file"
done
