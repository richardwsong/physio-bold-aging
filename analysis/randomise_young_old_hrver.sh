#!/bin/bash

# Description:
# Run FSL randomise for young vs. old baseline group comparisons for percent variance in BOLD explained by HRCO2, HR, and CO2 regressors.

# Load the FSL module for accessing FSL tools
module load FSL

# Define paths for mask, output directory, and input files
mask_path="data/masks/MNI152_T1_2mm_brain.nii"
output_dir="results/hrver_pve_results_gender_lf_hf_avg_hr" # Change based on the covariates used in the analysis
randomise_output_dir="$output_dir/randomise"

# Create the randomise output directory if it doesn't exist
mkdir -p "$randomise_output_dir"

# Convert design and contrast text files to FSL-compatible .mat format
Text2Vest "$output_dir/design_matrix.txt" "$output_dir/design.mat"
Text2Vest "$output_dir/contrast_matrix.txt" "$output_dir/design.con"

# Array of required files
declare -a required_files=(
    "$mask_path"
    "$output_dir/design.mat"
    "$output_dir/design.con"
    "$output_dir/hrco2_cov_young_old.nii.gz"
)

# Check if all required files exist
for file in "${required_files[@]}"; do
    if [ ! -e "$file" ]; then
        echo "Error: Required file not found: $file"
        exit 1
    fi
done

# Run randomise for HRCO2
echo "Running randomise for HRCO2..."
randomise -i "$output_dir/hrco2_cov_young_old.nii.gz" \
          -o "$randomise_output_dir/hrco2" \
          -m "$mask_path" \
          -d "$output_dir/design.mat" \
          -t "$output_dir/design.con" \
          -n 5000 -T

echo "All analyses completed successfully."
