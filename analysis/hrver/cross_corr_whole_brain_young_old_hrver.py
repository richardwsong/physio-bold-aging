# Creates Whole Brain Cross Correlation Maps for each lag and group comparison (young vs. old)

import os
import numpy as np
import nibabel as nib
import pandas as pd
from scipy.io import loadmat
from scipy.signal import detrend
import concurrent.futures
import subprocess
import warnings
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from utils.file_paths_hrver import get_file_paths
warnings.filterwarnings("ignore")

def run_randomise(input_4d, output_prefix, design_mat, contrast_mat, n_permutations=5000):
    """
    Runs FSL's `randomise` command to perform non-parametric permutation testing.

    Parameters:
        input_4d (str): Path to the 4D NIfTI file containing combined group data.
        output_prefix (str): Prefix for the output files generated by `randomise`.
        design_mat (str): Path to the design matrix file.
        contrast_mat (str): Path to the contrast matrix file.
        n_permutations (int): Number of permutations to run (default: 5000).
    """
    cmd = [
        'randomise',
        '-i', input_4d,
        '-o', output_prefix,
        '-d', design_mat,
        '-t', contrast_mat,
        '-T',
        '-n', str(n_permutations),
        '-m', 'data/masks/MNI152_T1_2mm_brain.nii'
    ]
    subprocess.run(cmd, check=True)

def save_design_contrast_matrices(group1_size, group2_size):
    """
    Creates and saves design and contrast matrices for group comparisons.

    Parameters:
        group1_size (int): Number of participants in group 1 (e.g., young).
        group2_size (int): Number of participants in group 2 (e.g., old).

    Outputs:
        - Saves `design.txt` and `contrast.txt` in the `whole_brain_lag_analysis` directory.
    """
    design_mat = np.zeros((group1_size + group2_size, 2))
    design_mat[:group1_size, 0] = 1
    design_mat[group1_size:, 1] = 1

    contrast_mat = np.array([[1, -1], [-1, 1]])

    np.savetxt("results/hrver/cross_corr_whole_brain_young_old/design.txt", design_mat, fmt="%d")
    np.savetxt("results/hrver/cross_corr_whole_brain_young_old/contrast.txt", contrast_mat, fmt="%d")

    os.system("Text2Vest results/hrver/cross_corr_whole_brain_young_old/design.txt results/hrver/cross_corr_whole_brain_young_old/design.mat")
    os.system("Text2Vest results/hrver/cross_corr_whole_brain_young_old/contrast.txt results/hrver/cross_corr_whole_brain_young_old/design.con")

def concatenate_results(result_list, output_file):
    """
    Concatenates results along the fourth dimension to create a 4D NIfTI file.

    Parameters:
        result_list (list): List of paths to 3D NIfTI files for concatenation.
        output_file (str): Path to save the concatenated 4D NIfTI file.
    """
    concatenated_data = np.stack(result_list, axis=-1)
    reference_img = nib.load(result_list[0])
    concatenated_img = nib.Nifti1Image(concatenated_data, reference_img.affine)
    nib.save(concatenated_img, output_file)

def fsl_config():
    """
    Configures the FSL environment by setting necessary environment variables.
    This function must be run before using any FSL tools.
    """
    os.environ['FSLDIR'] = '/accre/common/easybuild/software/FSL/6.0.6.5'
    os.environ['FSLOUTPUTTYPE'] = 'NIFTI_GZ'
    os.environ['PATH'] = f"{os.environ['FSLDIR']}/bin:{os.environ['PATH']}"
    
    # Additional FSL environment variables that might be needed
    os.environ['FSLMULTIFILEQUIT'] = 'TRUE'
    os.environ['FSLLOCKDIR'] = ''
    os.environ['FSLMACHINELIST'] = ''
    os.environ['FSLREMOTECALL'] = ''

def process_single_lag(subs_id, ages, lag, start_lag, x, y, z):
    """
    Processes data for a specific lag across all participants and separates results 
    by age group (young vs. old).

    Parameters:
        subs_id (np.ndarray): Array of subject IDs.
        ages (np.ndarray): Array of ages corresponding to the subject IDs.
        lag (int): Current lag index.
        start_lag (int): Start lag value (e.g., -2).
        x, y, z (int): Dimensions of the brain volume.

    Returns:
        tuple: 4D arrays for heart rate (HR) and respiratory variation (RV) correlations 
               for young and old groups.
    """
    hr_young_data = []
    hr_old_data = []
    co2_young_data = []
    co2_old_data = []
    
    # Process participants in parallel for this lag
    num_cores = min(os.cpu_count()//2, 16)
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cores) as executor:
        futures = {executor.submit(process_participant_single_lag, sub_id, lag, start_lag): (sub_id, age) 
                  for sub_id, age in zip(subs_id, ages)}
        
        for future in concurrent.futures.as_completed(futures):
            sub_id, age = futures[future]
            try:
                hr_corr, co2_corr = future.result()
                if age < 50:
                    hr_young_data.append(hr_corr)
                    co2_young_data.append(co2_corr)
                else:
                    hr_old_data.append(hr_corr)
                    co2_old_data.append(co2_corr)
                print(f"Lag {lag}: Completed participant {sub_id}")
            except Exception as e:
                print(f"Error processing participant {sub_id} for lag {lag}: {e}")
    
    # Convert lists to arrays and create 4D files
    hr_young_4d = np.stack(hr_young_data, axis=3) if hr_young_data else np.zeros((x, y, z, 0))
    hr_old_4d = np.stack(hr_old_data, axis=3) if hr_old_data else np.zeros((x, y, z, 0))
    co2_young_4d = np.stack(co2_young_data, axis=3) if co2_young_data else np.zeros((x, y, z, 0))
    co2_old_4d = np.stack(co2_old_data, axis=3) if co2_old_data else np.zeros((x, y, z, 0))
    
    return hr_young_4d, hr_old_4d, co2_young_4d, co2_old_4d

def process_participant_single_lag(sub_id, lag, start_lag):
    """
    Processes data for a specific participant and lag, calculating correlations 
    between physiological signals (HR, RV) and fMRI BOLD signals.

    Parameters:
        sub_id (str): Subject ID.
        lag (int): Current lag index.
        start_lag (int): Start lag value.

    Returns:
        tuple: 3D arrays for HR and RV correlations.
    """
    data = loadmat(get_file_paths('pre', sub_id)['physio'])
    bold_path = get_file_paths('pre', sub_id)['bold']

    hr = data['REGS']['hr'][0, 0].flatten()
    co2 = data['RESP']['etco2'][0, 0].flatten()

    bold_img = nib.load(bold_path)
    bold_data = bold_img.get_fdata()
    mask = nib.load("data/masks/MNI152_T1_2mm_brain.nii").get_fdata().astype(bool)

    # Calculate correlations for just this lag
    bold_reshaped = bold_data.reshape(-1, bold_data.shape[-1])
    mask_indices = np.where(mask.flatten())[0]
    bold_reshaped = bold_reshaped[mask_indices, :]

    mean_bold = np.mean(bold_reshaped, axis=1, keepdims=True)
    std_bold = np.std(bold_reshaped, axis=1, keepdims=True)
    bold_normalized = (bold_reshaped - mean_bold) / std_bold

    hr = detrend(hr - np.mean(hr)) / np.std(hr)
    co2 = detrend(co2 - np.mean(co2)) / np.std(co2)

    x, y, z, _ = bold_data.shape
    bold_denom = np.sqrt(np.sum(bold_normalized ** 2, axis=1))

    # Calculate for specific lag
    current_lag = start_lag + lag
    lagged_hr = np.roll(hr, current_lag)

    corr_hr = np.sum(bold_normalized * lagged_hr, axis=1) / (bold_denom * np.linalg.norm(lagged_hr))
    corr_co2 = np.sum(bold_normalized * co2, axis=1) / (bold_denom * np.linalg.norm(co2))

    # Reshape back to 3D
    hr_3d = np.zeros((x, y, z))
    co2_3d = np.zeros((x, y, z))
    hr_3d.flat[mask_indices] = corr_hr
    co2_3d.flat[mask_indices] = corr_co2

    return hr_3d, co2_3d

def main(start_lag=-1, end_lag=9):
    """
    Main function to perform the entire analysis pipeline:
    1. Configure FSL environment.
    2. Load demographic and physiological data.
    3. Iterate over lags, calculating correlations for HR and RV.
    4. Run FSL's `randomise` for group comparisons.
    5. Save results as 4D NIfTI files.

    Parameters:
        start_lag (int): Starting lag for analysis (default: -2).
        end_lag (int): Ending lag for analysis (default: 20).

    Outputs:
        - Statistical maps for each lag and group comparison (e.g., `hr_young-gt-old_tstat_4d.nii.gz`).
        - Design and contrast matrices for reproducibility.
    """
    fsl_config()
    df = pd.read_csv("data/hrver_ses_pre_age_gender.csv")
    subs_id = df.iloc[:, 0].values
    ages = df.iloc[:, 2].values

    output_dir = "results/hrver/cross_corr_whole_brain_young_old"
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "temp"), exist_ok=True)

    x, y, z, _ = nib.load(get_file_paths('pre', subs_id[0])['bold']).shape
    num_lags = end_lag - start_lag + 1

    # Lists to store randomise results for both directions
    hr_tstat1_list = []  # Young > Old
    hr_tstat2_list = []  # Old > Young
    hr_corrp1_list = []  # TFCE p-values for Young > Old
    hr_corrp2_list = []  # TFCE p-values for Old > Young
    co2_tstat1_list = []
    co2_tstat2_list = []
    co2_corrp1_list = []
    co2_corrp2_list = []

    # Process one lag at a time
    for lag_idx in range(num_lags):
        print(f"\nProcessing lag {start_lag + lag_idx} ({lag_idx + 1}/{num_lags})")
        
        # Get correlations for this lag
        hr_young_4d, hr_old_4d, co2_young_4d, co2_old_4d = process_single_lag(
            subs_id, ages, lag_idx, start_lag, x, y, z
        )

        # Run randomise for this lag
        for measure, young_data, old_data, tstat1_list, tstat2_list, corrp1_list, corrp2_list in [
            ('hr', hr_young_4d, hr_old_4d, hr_tstat1_list, hr_tstat2_list, hr_corrp1_list, hr_corrp2_list),
            ('co2', co2_young_4d, co2_old_4d, co2_tstat1_list, co2_tstat2_list, co2_corrp1_list, co2_corrp2_list)
        ]:
            # Combine data and save temporary 4D file
            combined_data = np.concatenate([young_data, old_data], axis=3)
            temp_4d = os.path.join(output_dir, "temp", f"{measure}_lag_{lag_idx}_4D.nii.gz")
            nib.save(nib.Nifti1Image(combined_data, None), temp_4d)

            # Save design/contrast matrices if first lag
            if lag_idx == 0:
                save_design_contrast_matrices(
                    young_data.shape[3], 
                    old_data.shape[3]
                )

            # Run randomise
            pve_dir = "results/hrver/pve_results_gender_lf_hf_avg_hr_avg_co2" # If a PVE analysis was done, you can use design and contrast matrices from the PVE analysis
            run_randomise(
                temp_4d,
                os.path.join(output_dir, "temp", f"{measure}_lag_{lag_idx}"),
                os.path.join(pve_dir if os.path.exists(pve_dir) else output_dir, "design.mat"),
                os.path.join(pve_dir if os.path.exists(pve_dir) else output_dir, "design.con"),
                n_permutations=1000
            )

            # Store results for both directions
            tstat1_list.append(
                nib.load(os.path.join(output_dir, "temp", f"{measure}_lag_{lag_idx}_tstat1.nii.gz")).get_fdata()
            )
            tstat2_list.append(
                nib.load(os.path.join(output_dir, "temp", f"{measure}_lag_{lag_idx}_tstat2.nii.gz")).get_fdata()
            )
            corrp1_list.append(
                nib.load(os.path.join(output_dir, "temp", f"{measure}_lag_{lag_idx}_tfce_corrp_tstat1.nii.gz")).get_fdata()
            )
            corrp2_list.append(
                nib.load(os.path.join(output_dir, "temp", f"{measure}_lag_{lag_idx}_tfce_corrp_tstat2.nii.gz")).get_fdata()
            )

            # Clean up temporary files
            os.remove(temp_4d)

        # Clear any large arrays from memory
        del hr_young_4d, hr_old_4d, co2_young_4d, co2_old_4d, combined_data
    
    # Save final 4D results
    affine = nib.load(get_file_paths('pre', subs_id[0])['bold']).affine
    for data_list, filename in [
        (hr_tstat1_list, os.path.join(output_dir, "hr_young-gt-old_tstat_4d.nii.gz")),
        (hr_tstat2_list, os.path.join(output_dir, "hr_old-gt-young_tstat_4d.nii.gz")),
        (hr_corrp1_list, os.path.join(output_dir, "hr_young-gt-old_tfce_corrp_4d.nii.gz")),
        (hr_corrp2_list, os.path.join(output_dir, "hr_old-gt-young_tfce_corrp_4d.nii.gz")),
        (co2_tstat1_list, os.path.join(output_dir, "co2_young-gt-old_tstat_4d.nii.gz")),
        (co2_tstat2_list, os.path.join(output_dir, "co2_old-gt-young_tstat_4d.nii.gz")),
        (co2_corrp1_list, os.path.join(output_dir, "co2_young-gt-old_tfce_corrp_4d.nii.gz")),
        (co2_corrp2_list, os.path.join(output_dir, "co2_old-gt-young_tfce_corrp_4d.nii.gz"))
    ]:
        combined = np.stack(data_list, axis=3)
        nib.save(nib.Nifti1Image(combined, affine), filename)

    print("Analysis complete!")

if __name__ == "__main__":
    main()
