#!/bin/bash

# Define paths
INPUT_DIR="./bed_files"
OUTPUT_DIR="./fasta_outputs"
REF_GENOME="/home/photon/cfg_project/hg38.fa"

# 1. Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# 2. Loop through the BED files
for file in "$INPUT_DIR"/*.bed; do
    # Check if the glob actually matched anything
    [ -e "$file" ] || continue

    # Extract the filename without the extension
    base=$(basename "$file" .bed)
    
    echo "Extracting sequences for: $base"
    
    # 3. Run bedtools and save to the new directory
    bedtools getfasta \
        -fi "$REF_GENOME" \
        -bed "$file" \
        -fo "$OUTPUT_DIR/${base}.fa"
done

echo "Done! Files are saved in $OUTPUT_DIR"