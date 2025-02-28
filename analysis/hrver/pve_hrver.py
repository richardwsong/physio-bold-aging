import os
import numpy as np
import nibabel as nib
from scipy.io import loadmat, savemat
from scipy.linalg import pinv
from scipy.signal import correlate, detrend, convolve
from scipy.signal import find_peaks, butter, filtfilt, welch
from scipy.interpolate import interp1d
import concurrent.futures
import pandas as pd
import matplotlib.pyplot as plt
import warnings 
import neurokit2 as nk
import sys
sys.path.append('/data1/neurdylab/songrw/scripts/physio-bold-aging')  # REPLACE WITH YOUR PATH  

from utils.file_paths_hrver import get_file_paths
warnings.filterwarnings("ignore")

dt = 2.4 # HRV-ER dataset has TR of 2.4 s, NKI uses 1.4 s TR 
t = np.arange(0, 40 + dt, dt)  # assume 40-s long total for each basis set 

# CRF basis
crf_p = 0.3 * t**2.7 * np.exp(-t / 1.6) - 1.05 * np.exp(-((t - 12)**2) / 18) # primary (Chang et al., 2009)
crf_t1 = 1.94 * t**1.7 * np.exp(-t / 1.6) - 0.45 * t**2.7 * np.exp(-t / 1.6) # time derivative 1 (Chen et al., 2020)
crf_t2 = 0.55 * (t - 12) * np.exp(-((t - 12)**2) / 18) # time derivative 2 (Chen et al., 2020)
crf_d1 = 0.056 * t**3.7 * np.exp(-t / 1.6) # dispersion 1 (Chen et al., 2020)
crf_d2 = 0.15 * (t - 12)**2 * np.exp(-((t - 12)**2) / 18) # dispersion 2 (Chen et al., 2020)

# Golestani CO2: Double gamma HRF with the parameterization: h(t) = a * t^b * exp(-t/c) + d * t^e * exp(-t/f)
def double_gamma_param(t, a=0.2, b=2.02, c=3.1, d=8.8*10**-12, e=23.65, f=0.31):
    """Double gamma HRF parameterized as two gamma functions."""
    gamma1 = a * t**b * np.exp(-t / c)
    gamma2 = d * t**e * np.exp(-t / f)
    return gamma1 + gamma2

# Temporal derivatives
def temporal_derivatives_param(t, a=0.2, b=2.02, c=3.1, d=8.8*10**-12, e=23.65, f=0.31):
    """First temporal derivatives for the two gamma components."""
    dgamma1_dt = a * (b * t**(b-1) - t**b / c) * np.exp(-t / c)
    dgamma2_dt = d * (e * t**(e-1) - t**e / f) * np.exp(-t / f)
    return dgamma1_dt + dgamma2_dt

# Normalize Function
def normalize_function(f):
    max_amplitude = np.max(np.abs(f))
    return f / max_amplitude

# Compute CO2 Response Function and its derivatives
co2rf_p = normalize_function(double_gamma_param(t))
co2rf_t1 = normalize_function(temporal_derivatives_param(t))

def create_physio_basis_HRV_ER(etco2, hr):
    """
    Create a physiological basis set for HRV and ER analysis by convolving 
    input physiological signals (etco2, heart rate, respiratory variation) 
    with their respective basis functions.

    This function generates a design matrix for modeling physiological effects 
    in fMRI data, specifically focusing on HRV (Heart Rate Variability) and 
    ER (End-tidal CO2 Responses). It uses convolution with predefined basis 
    functions for HRF (Hemodynamic Response Function), CRF (Cardiac Response 
    Function), and RRF (Respiratory Response Function).

    Parameters:
        etco2 (np.ndarray): End-tidal CO2 signal, a 1D array.
        hr (np.ndarray): Heart rate signal, a 1D array.

    Global Variables:
        hrf_p, hrf_t, hrf_d: HRF basis functions for CO2.
        crf_p, crf_t1, crf_t2, crf_d1, crf_d2: CRF basis functions for heart rate.
        rrf_p, rrf_t1, rrf_t2, rrf_d1, rrf_d2: RRF basis functions for respiratory variation.

    Returns:
        np.ndarray: A concatenated design matrix `X_physio`, with the following structure:
            - Column 1: Raw etco2 signal.
            - Columns 2-4: HRF convolved with etco2.
            - Columns 5-9: CRF convolved with heart rate.
            - Columns 10-14: RRF convolved with respiratory variation.
    """
    global crf_p, crf_t1, crf_t2, crf_d1, crf_d2, rrf_p, rrf_t1, rrf_t2, rrf_d1, rrf_d2, co2rf_p, co2rf_t1

    # Convolve etco2 with HRF basis
    CO2 = etco2
    CO2 = CO2.flatten()  # Ensure column vector

    CO2_conv = np.zeros((len(CO2), 3))
    CO2_conv[:, 0] = convolve(CO2, co2rf_p, mode='full')[:len(CO2)]
    CO2_conv[:, 1] = convolve(CO2, co2rf_t1, mode='full')[:len(CO2)]
    CO2_basis_regs = CO2_conv

    # Convolve HR with CRF basis
    HR = hr 
    HR = HR.flatten()  # Ensure column vector

    HR_conv = np.zeros((len(HR), 5))
    HR_conv[:, 0] = convolve(HR, crf_p, mode='full')[:len(HR)]
    HR_conv[:, 1] = convolve(HR, crf_t1, mode='full')[:len(HR)]
    HR_conv[:, 2] = convolve(HR, crf_t2, mode='full')[:len(HR)]
    HR_conv[:, 3] = convolve(HR, crf_d1, mode='full')[:len(HR)]
    HR_conv[:, 4] = convolve(HR, crf_d2, mode='full')[:len(HR)]
    HR_basis_regs = HR_conv

    # Concatenate etco2, CO2, HR, and RV basis functions
    X_physio = np.hstack([CO2_basis_regs, HR_basis_regs])

    return X_physio

def load_nii(file_path):
    """
    Load a NIfTI file and return the image data as a NumPy array.

    Parameters:
        file_path (str): Path to the NIfTI file.

    Returns:
        np.ndarray: The image data contained in the NIfTI file.
        np.ndarray: The affine transformation matrix from the NIfTI file.
    """
    temp = nib.load(file_path)
    return temp.get_fdata(), temp.affine


def load_pre(): 
    """
    Load and return the pre-session dataset.

    This function reads a CSV file containing information about subjects
    during a pre-session experiment, including age and gender.

    Returns:
        pd.DataFrame: A DataFrame containing the pre-session data.
    """
    infile_pre = "data/hrver_ses_pre_age_gender.csv"
    df_pre = pd.read_csv(infile_pre)
    return df_pre

def load_post():
    """
    Load and return the post-session dataset.

    This function reads a CSV file containing information about subjects
    during a post-session experiment, including age and gender.

    Returns:
        pd.DataFrame: A DataFrame containing the post-session data.
    """
    infile_post = "data/hrver_ses_post_age_gender.csv"
    df_post = pd.read_csv(infile_post)
    return df_post

def load_shared():
    """
    Load and return a subset of the pre-session dataset containing only
    subjects who participated in both the pre- and post-session experiments.

    This function identifies the common subject IDs between the pre-session
    and post-session datasets, and filters the pre-session dataset to include
    only these shared subjects.

    Returns:
        pd.DataFrame: A DataFrame containing data for shared subjects in the pre-session dataset.
    """
    # Load the pre-session data and extract subject IDs
    infile_pre = "data/hrver_ses_pre_age_gender.csv"
    df_pre = pd.read_csv(infile_pre)
    sub_pre = df_pre['subs_id'].values

    # Load the post-session data and extract subject IDs
    infile_post = "data/hrver_ses_post_age_gender.csv"
    df_post = pd.read_csv(infile_post)
    sub_post = df_post['subs_id'].values

    # Find the intersection of subject IDs between the two datasets
    shared_subs = np.intersect1d(sub_pre, sub_post)

    # Filter the pre-session dataset to include only the shared subjects
    df_shared = df_pre[df_pre['subs_id'].isin(shared_subs)]
    return df_shared

def calculate_breathing_rate_new(resp_signal, fs):
    # Process the respiration signal
    rsp_signals, _ = nk.rsp_process(resp_signal, sampling_rate=fs)

    # Extract the mean breathing rate
    return rsp_signals["RSP_Rate"].mean()

# Function to load and calculate HRV metrics and breathing rate
def calc_physio_features(name, session):
    data = loadmat(get_file_paths(session, name)['physio']) # Path to preprocessed physio data
    ibi = data['OUT_p']['IBI_clean'][0][0].flatten()
    ibi = ibi[~np.isnan(ibi)]
    ibi_ms = ibi * 1000  # Convert to milliseconds
    rmssd = np.log(np.sqrt(np.mean(np.diff(ibi_ms) ** 2)))

    hr = data['REGS']['hr'][0][0].flatten()
    co2 = data['RESP']['etco2'][0][0].flatten()
    resp = data['RESP']['cng_ds'][0][0].flatten()
    breathing_rate = calculate_breathing_rate_new(resp, fs=1000)
    breathing_rate = np.mean(breathing_rate)    

    avg_co2 = np.mean(co2)
    avg_hr = np.mean(hr)

    # Calculate power spectral density for LF and HF
    f, psd = welch(ibi_ms, fs=1 / np.mean(ibi), nperseg=len(ibi))
    lf_band = (f >= 0.04) & (f < 0.15)
    hf_band = (f >= 0.15) & (f < 0.4)

    lf = np.log(np.trapz(psd[lf_band], f[lf_band]))
    hf = np.log(np.trapz(psd[hf_band], f[hf_band]))
    hr = detrend(hr)
    co2 = detrend(co2)
    print(f"Processed Physio {name}")
    return hr, co2, lf, hf, avg_hr, avg_co2, breathing_rate

def calc_percent_variance_explained(name, mask_indices, session): 
    """
    Calculate the percent variance in fMRI data explained by various physiological signals.

    This function uses physiological signals (heart rate, end-tidal CO2, and respiratory variation)
    as regressors to estimate the variance they explain in fMRI voxel time series, 
    based on a mask of brain regions of interest. 

    The percent variance explained is calculated for:
    - Heart rate and CO2 combined
    - Heart rate only
    - CO2 only
    - Respiratory variation (RV) only
    - Heart rate and respiratory variation combined

    Parameters:
        name (str): Subject's name or identifier.
        session (str): Session identifier (e.g., "pre", "post").
        mask_indices (np.ndarray): Indices of voxels within the brain mask (based on an MNI template).

    Returns:
        tuple: Contains the percent variance explained for different regressors:
            - percent_variance_hrco2 (np.ndarray): Variance explained by heart rate and CO2 combined (voxels x 1).
            - percent_variance_hr (np.ndarray): Variance explained by heart rate only (voxels x 1).
            - percent_variance_co2 (np.ndarray): Variance explained by CO2 only (voxels x 1).
            - percent_variance_rv (np.ndarray): Variance explained by RV only (voxels x 1).
            - percent_variance_hrrv (np.ndarray): Variance explained by heart rate and RV combined (voxels x 1).
    """
    
    hr, co2, lf, hf, avg_hr, avg_co2, br = calc_physio_features(name, session)      
    nii_path = get_file_paths(session, name)['bold']

    nn, _ = load_nii(nii_path)
    dims = nn.shape
    V = nn.reshape(-1, dims[3]) # Reshape to 2D matrix (voxels x time)
    Y = V.T  # Transpose into (time x voxels) matrix

    # Normalize each column (voxel) to reflect percent signal change
    for k in range(Y.shape[1]):
        avg = np.mean(Y[:, k])
        if avg != 0:
            Y[:, k] = (Y[:, k] - avg) / avg

    Y_clean = Y.T # (voxel x time)
    Y_clean = Y_clean[mask_indices, :] # mask the brain regions inside the mask (MNI template)

    X = create_physio_basis_HRV_ER(co2, hr) # X is the design matrix for physiological regressors  

    # HRRV only
    X1 = np.hstack([np.ones((X.shape[0], 1)), X[:, :]])
    B1 = pinv(X1) @ Y_clean.T  # Least-squares regression coefficients
    percent_variance_hrco2 = np.var((X1 @ B1).T, axis=1, keepdims=True) / np.var(Y_clean, axis=1, keepdims=True)

    return percent_variance_hrco2, lf, hf, avg_hr, avg_co2, br 

# Function to process data for a single participant
def process_participant(i, subs_id, mask_indices, session):
    """
    Process data for a single participant.

    Parameters:
        i (int): Participant index.
        subs_id (list): List of subject IDs.
        mask_indices (np.ndarray): Indices of voxels within the brain mask.

    Returns:
        tuple: Participant index and arrays for percent variance explained 
               (HR, CO2, RV, HRRV, HR+CO2).
    """
    hrrv, lf, hf, avg_hr, sd_rv, br = calc_percent_variance_explained(subs_id[i], mask_indices, session)
    return i, hrrv, lf, hf, avg_hr, sd_rv, br

# Function to create a design matrix
def create_design_matrix(df, physio_stats, covariates=None):
    """
    Create a design matrix for group comparison (young vs. old) with flexible covariates.

    Parameters:
        df (pd.DataFrame): DataFrame containing participant information, including age.
        physio_stats (np.ndarray): Array of physiological statistics (lf, hf, avg_hr, avg_co2, br).
        covariates (list, optional): List of covariates to include. Options are:
            - 'gender': Include gender as covariate
            - 'lf': Include low frequency HRV
            - 'hf': Include high frequency HRV
            - 'avg_hr': Include average heart rate
            - 'avg_co2': Include average CO2
            - 'br': Include breathing rate
            If None, only includes young/old groups without covariates.

    Returns:
        np.ndarray: Design matrix with specified columns
    """
    # Base design matrix with young/old groups
    young = np.array(df['age'] < 50, dtype=float)
    old = np.array(df['age'] >= 50, dtype=float)
    design_matrix = np.concatenate((young[:, np.newaxis], old[:, np.newaxis]), axis=1)

    if covariates is None:
        return design_matrix

    # Dictionary mapping covariate names to their data
    covariate_map = {
        'gender': np.array(df['gender'] == 'MALE', dtype=float)[:, np.newaxis],
        'lf': physio_stats[:, 0:1],
        'hf': physio_stats[:, 1:2],
        'avg_hr': physio_stats[:, 2:3],
        'avg_co2': physio_stats[:, 3:4],
        'br': physio_stats[:, 4:5]
    }

    # Add requested covariates
    for cov in covariates:
        if cov in covariate_map:
            design_matrix = np.concatenate((design_matrix, covariate_map[cov]), axis=1)
        else:
            print(f"Warning: Covariate '{cov}' not recognized and will be skipped")

    return design_matrix

def create_contrast_matrix(n_covariates=0):
    """
    Create a contrast matrix for group comparison (young vs. old).

    Parameters:
        n_covariates (int): Number of covariates in the design matrix

    Returns:
        np.ndarray: Contrast matrix for statistical analysis.
    """
    # Create contrast vector with zeros for covariates
    contrast_vector = np.zeros(2 + n_covariates)
    contrast_vector[0:2] = [1, -1]  # young vs old contrast
    return np.vstack([contrast_vector, -contrast_vector])

def update_brain_map(data_update, mask_indices):
    """
    Map flattened data back into brain space.

    Parameters:
        data_update (np.ndarray): Flattened data.
        mask_indices (np.ndarray): Indices of brain mask voxels.

    Returns:
        np.ndarray: Reshaped data in brain space (91 x 109 x 91).
    """
    brain_map = np.zeros((91 * 109 * 91, 1))
    brain_map[mask_indices] = data_update
    return brain_map.reshape(91, 109, 91)

def save_nifti(data, affine, file_path):
    """
    Save data as a NIfTI file.

    Parameters:
        data (np.ndarray): Data to be saved.
        affine (np.ndarray): Affine transformation matrix.
        file_path (str): File path for saving the NIfTI file.
    """
    nii = nib.Nifti1Image(data, affine)
    nib.save(nii, file_path)

# Main function
def main():
    """
    Main function to perform the whole-brain analysis.
    - Loads demographic and brain mask data.
    - Initializes arrays for storing variance explained results.
    - Processes each participant's data in parallel.
    - Saves variance explained results as NIfTI files.
    - Generates and saves design and contrast matrices for further analysis.
    """
    # Load demographic data and brain mask
    session = 'pre'
    df = load_pre() 
    subs = df['sub_id'].values
    age = df['age'].values
    gender = df['gender'].values
    N = len(subs)

    mask, affine = load_nii("data/masks/MNI152_T1_2mm_brain.nii")
    mask = mask.astype(bool)
    mask_indices = np.where(mask.flatten())[0]

    # Initialize result arrays
    hrrv_cov = np.zeros((91, 109, 91, N))
    physio_stats = np.zeros((N, 5))

    # Parallel processing of participants
    num_cores = min(os.cpu_count()//2, 16)
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cores) as executor: 
        try:
            futures = [executor.submit(process_participant, i, subs, mask_indices, session) for i in range(N)]
            for future in concurrent.futures.as_completed(futures):
                index, hrrv_cov_update, lf, hf, avg_hr, avg_co2, br = future.result()

                # Reshape results back into brain space
                hrrv_cov[:, :, :, index] = update_brain_map(hrrv_cov_update, mask_indices)
                physio_stats[index, :] = lf, hf, avg_hr, avg_co2, br
                print(f"Processing complete for participant {subs[index]}")
        except Exception as e:
            print(f"Error: {e}")

    # Specify which covariates to include
    covariates = ['gender', 'lf', 'hf', 'avg_hr', 'avg_co2']  # Example: only control for gender

    # Create output directory name that includes all covariates
    output_dir = 'results/hrver/pve_results_' + '_'.join(covariates)
    os.makedirs(output_dir, exist_ok=True)

    # Save variance explained results as NIfTI files
    save_nifti(hrrv_cov, affine, os.path.join(output_dir, 'hrco2_cov_young_old.nii.gz'))

    # Generate and save design and contrast matrices
    design = create_design_matrix(df, physio_stats, covariates=covariates)
    np.savetxt(os.path.join(output_dir, 'design_matrix.txt'), design, fmt='%.6f')

    contrast = create_contrast_matrix(n_covariates=len(covariates))
    np.savetxt(os.path.join(output_dir, 'contrast_matrix.txt'), contrast, fmt='%d')
    
    print("Done!")
    #os.system('bash /path/to/randomise_young_old_baseline.sh') # Run randomise script for group comparison

if __name__ == "__main__":
    main()