#!/bin/bash

# Check if input directory is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <path_to_vcf_directory>"
    echo "Example: add_chr_prefix_tovcf.sh /path/to/my/files"
    exit 1
fi

INPUT_DIR="$1"

# Check if directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Directory '$INPUT_DIR' does not exist."
    exit 1
fi

module load htslib

echo "Processing VCF files in: $INPUT_DIR"

# Loop through vcf files
for f in "$INPUT_DIR"/*.vcf.gz; do
    # Check if file exists
    [ -e "$f" ] || continue

    echo "Processing $f..."

    # Create temp file path
    temp_file="${f}.temp"

    # Add 'chr' prefix to chromosome names and update ID fields, then compress
    zcat "$f" | \
    sed -e '/^[^#]/ s/^/chr/' -e 's/ID=\([0-9XY]\)/ID=chr\1/' | \
    bgzip > "$temp_file"

    # Safety check: verify temp file was created and has size
    if [ -s "$temp_file" ]; then
        mv "$temp_file" "$f"
        echo "Successfully added 'chr' prefix to $f"
    else
        echo "Error: Failed to process $f. Original file kept."
        rm -f "$temp_file"
    fi
done