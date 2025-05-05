#!/bin/bash

# Set the input directory (expand ~ to the full home directory path)
input_dir=""

# Set the output directory
output_dir=""

# Check if the input directory exists
if [[ ! -d "$input_dir" ]]; then
    echo "Error: Input directory $input_dir does not exist. Exiting."
    exit 1
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through all R1 files in the input directory and process them with their corresponding R2 files
for r1_file in "$input_dir"/*_1.fastq.gz; do
    # Get the base sample name (everything before _1.fastq.gz)
    base_name=$(basename "$r1_file" _1.fastq.gz)
    
    # Define the R2 file name
    r2_file="${input_dir}/${base_name}_2.fastq.gz"
    
    # Check if the R2 file exists
    if [[ ! -f "$r2_file" ]]; then
        echo "Warning: Paired file for $r1_file not found (expected $r2_file). Skipping."
        continue
    fi

    # Define output file names
    output_r1="${output_dir}/${base_name}_1.fastq.gz"
    output_r2="${output_dir}/${base_name}_2.fastq.gz"

    # Run cutadapt with overlap and error rate settings
    cutadapt -u 26 -U 26 --overlap 3 --error-rate 0.1 --cores 4 -o "$output_r1" -p "$output_r2" "$r1_file" "$r2_file"

    # Check the exit status of cutadapt
    if [[ $? -ne 0 ]]; then
        echo "Error trimming $r1_file and $r2_file. Check your files or parameters."
    else
        echo "Successfully trimmed $r1_file and $r2_file. Output saved to $output_dir."
    fi
done