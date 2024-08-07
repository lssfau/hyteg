#!/bin/bash

# Define the source and destination directories
source_dir="./output_pgfplot"
destination_dir="test"

# Use find to locate files containing "None" and exclude directories
find "$source_dir" -type f -exec grep -l "None" {} + | while IFS= read -r file; do
    # Remove lines containing "None" and preserve the directory structure in the destination folder
    mkdir -p "$(dirname "$destination_dir/${file#$source_dir}")"
    grep -v "None" "$file" > "$destination_dir/${file#$source_dir}"
done
