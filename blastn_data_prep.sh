#!/bin/bash

# Define paths to your input (blastn_results.txt) and output files
INPUT_FILE="/path/to/blastn_results.txt"
OUTPUT_FILE="/path/to/processed_blastn_results.txt"

# Run the Python script
python3 << EOF
import re

# Function to determine the downsampling value
def get_downsampling(name):
    if 'min' in name:
        return 'min'
    elif 'max' in name:
        return 'max'
    elif 'med' in name:
        return 'med'
    else:
        return 'unknown'

# Function to extract the segment
def get_segment(name):
    segments = ['_HA', '_NA', '_NP', '_NS', '_M', '_PA', '_PB1', '_PB2']
    for segment in segments:
        if segment in name:
            return segment.strip('_')
    return 'unknown'

# Function to determine the technique
def get_technique(name):
    if 'bcftools' in name:
        return 'bcftools'
    elif 'irma' in name:
        return 'irma'
    elif 'ivar' in name:
        return 'ivar'
    else:
        return 'unknown'

# Process the file and create the new data structure
def process_file(input_file_path, output_file_path):
    processed_data = []
    with open(input_file_path, 'r') as file:
        for line in file:
            parts = line.split()
            if len(parts) >= 13 and parts[0] != 'ConsensusFile':
                downsampling = get_downsampling(parts[0])
                nucleotide = 'cDNA'
                segment = get_segment(parts[0])
                technique = get_technique(parts[0])
                bit_score = parts[12]
                processed_data.append([downsampling, nucleotide, segment, technique, bit_score])

    # Create a new file with the processed data
    with open(output_file_path, 'w') as file:
        file.write("downsampling\tnucleotide\tsegment\ttechnique\tbit_score\n")
        for row in processed_data:
            file.write("\\t".join(row) + "\\n")

# Run the processing function
process_file('$INPUT_FILE', '$OUTPUT_FILE')
EOF

echo "Processing completed. Output file: $OUTPUT_FILE"
