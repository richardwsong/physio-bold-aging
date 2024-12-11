#!/bin/bash

# Description:
# Run FSL randomise for young vs. old baseline group comparisons for percent variance in BOLD explained by HRCO2, HR, and CO2 regressors.

# Load the FSL module for accessing FSL tools
module load FSL

# Define paths for mask, output directory, and input files
mask_path="/path/to/MNI152_T1_2mm_brain.nii"
output_dir="/path/to/output_dir"
randomise_output_dir="$output_dir/randomise"

# Create the randomise output directory if it doesn't exist
mkdir -p "$randomise_output_dir"

# Array of required files
declare -a required_files=(
    "$mask_path"
    "$output_dir/design.txt"
    "$output_dir/design.con"
    "$output_dir/hrco2_cov_young_old.nii.gz"
    "$output_dir/hr_cov_young_old.nii.gz"
    "$output_dir/co2_cov_young_old.nii.gz"
)

# Check if all required files exist
for file in "${required_files[@]}"; do
    if [ ! -e "$file" ]; then
        echo "Error: Required file not found: $file"
        exit 1
    fi
done

# Convert design and contrast text files to FSL-compatible .mat format
Text2Vest "$output_dir/design.txt" "$output_dir/design.mat"
Text2Vest "$output_dir/design.con" "$output_dir/design.con"

# Run randomise for HRCO2
echo "Running randomise for HRCO2..."
randomise -i "$output_dir/hrco2_cov_young_old.nii.gz" \
          -o "$randomise_output_dir/hrco2" \
          -m "$mask_path" \
          -d "$output_dir/design.mat" \
          -t "$output_dir/design.con" \
          -n 5000 -T

# Run randomise for HR
echo "Running randomise for HR..."
randomise -i "$output_dir/hr_cov_young_old.nii.gz" \
          -o "$randomise_output_dir/hr" \
          -m "$mask_path" \
          -d "$output_dir/design.mat" \
          -t "$output_dir/design.con" \
          -n 5000 -T

# Run randomise for CO2
echo "Running randomise for CO2..."
randomise -i "$output_dir/co2_cov_young_old.nii.gz" \
          -o "$randomise_output_dir/co2" \
          -m "$mask_path" \
          -d "$output_dir/design.mat" \
          -t "$output_dir/design.con" \
          -n 5000 -T

echo "All analyses completed successfully."
