#!/bin/bash

# Description: Run randomise analyses for comparing young vs. old percent variance in BOLD explained by ... 
# HRCO2, HR, and CO2 regressors post-pre for both Osc conditions

# Load the FSL module for accessing FSL tools
module load FSL

# Define paths for design files, mask, and output directories
mask_path="/path/to/masks/MNI152_T1_2mm_brain.nii"
output_dir="/path/to/output_dir"
randomise_output_dir="$output_dir/randomise"

# Create randomise output directory if it doesn't exist
mkdir -p "$randomise_output_dir"

# Check if required files exist
required_files=(
    "$design_con_file"
    "$mask_path"
    "$output_dir/hrco2_young_vs_old_osc_plus_design.txt"
    "$output_dir/hrco2_young_vs_old_osc_minus_design.txt"
    "$output_dir/hrco2_young_vs_old_osc_plus.nii.gz"
    "$output_dir/hrco2_young_vs_old_osc_minus.nii.gz"
    "$output_dir/hr_young_vs_old_osc_plus.nii.gz"
    "$output_dir/hr_young_vs_old_osc_minus.nii.gz"
    "$output_dir/co2_young_vs_old_osc_plus.nii.gz"
    "$output_dir/co2_young_vs_old_osc_minus.nii.gz"
)

for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "Error: Required file not found - $file"
        exit 1
    fi
done

# Convert text files to FSL-compatible .mat design matrices
Text2Vest "$output_dir/young_vs_old_osc_plus_design.txt" "$output_dir/young_vs_old_osc_plus_design.mat"
Text2Vest "$output_dir/young_vs_old_osc_minus_design.txt" "$output_dir/young_vs_old_osc_minus_design.mat"
Text2Vest "$output_dir/contrast.txt" "$output_dir/design.con"
design_con_file="$output_dir/design.con"

# Run randomise analyses
echo "Running randomise analyses..."

randomise -i "$output_dir/hrco2_young_vs_old_osc_minus.nii.gz" \
          -o "$randomise_output_dir/hrco2_young_vs_old_osc_minus" \
          -d "$output_dir/hrco2_young_vs_old_osc_minus_design.mat" \
          -t "$design_con_file" \
          -m "$mask_path" -T -n 5000

randomise -i "$output_dir/hrco2_young_vs_old_osc_plus.nii.gz" \
          -o "$randomise_output_dir/hrco2_young_vs_old_osc_plus" \
          -d "$output_dir/hrco2_young_vs_old_osc_plus_design.mat" \
          -t "$design_con_file" \
          -m "$mask_path" -T -n 5000

randomise -i "$output_dir/hr_young_vs_old_osc_minus.nii.gz" \
          -o "$randomise_output_dir/hr_young_vs_old_osc_minus" \
          -d "$output_dir/hrco2_young_vs_old_osc_minus_design.mat" \
          -t "$design_con_file" \
          -m "$mask_path" -T -n 5000

randomise -i "$output_dir/hr_young_vs_old_osc_plus.nii.gz" \
          -o "$randomise_output_dir/hr_young_vs_old_osc_plus" \
          -d "$output_dir/hrco2_young_vs_old_osc_plus_design.mat" \
          -t "$design_con_file" \
          -m "$mask_path" -T -n 5000

randomise -i "$output_dir/co2_young_vs_old_osc_minus.nii.gz" \
          -o "$randomise_output_dir/co2_young_vs_old_osc_minus" \
          -d "$output_dir/hrco2_young_vs_old_osc_minus_design.mat" \
          -t "$design_con_file" \
          -m "$mask_path" -T -n 5000

randomise -i "$output_dir/co2_young_vs_old_osc_plus.nii.gz" \
          -o "$randomise_output_dir/co2_young_vs_old_osc_plus" \
          -d "$output_dir/hrco2_young_vs_old_osc_plus_design.mat" \
          -t "$design_con_file" \
          -m "$mask_path" -T -n 5000

echo "Randomise analyses completed successfully."
