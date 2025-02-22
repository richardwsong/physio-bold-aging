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
warnings.filterwarnings("ignore")

dt = 1.4 # HRV-ER dataset has TR of 2.4 s, NKI uses 1.4 s TR 
t = np.arange(0, 40 + dt, dt)  # assume 40-s long total for each basis set 

# CRF basis
crf_p = 0.3 * t**2.7 * np.exp(-t / 1.6) - 1.05 * np.exp(-((t - 12)**2) / 18) # primary (Chang et al., 2009)
crf_t1 = 1.94 * t**1.7 * np.exp(-t / 1.6) - 0.45 * t**2.7 * np.exp(-t / 1.6) # time derivative 1 (Chen et al., 2020)
crf_t2 = 0.55 * (t - 12) * np.exp(-((t - 12)**2) / 18) # time derivative 2 (Chen et al., 2020)
crf_d1 = 0.056 * t**3.7 * np.exp(-t / 1.6) # dispersion 1 (Chen et al., 2020)
crf_d2 = 0.15 * (t - 12)**2 * np.exp(-((t - 12)**2) / 18) # dispersion 2 (Chen et al., 2020)

# RRF basis
rrf_p = 0.6 * t**2.1 * np.exp(-t / 1.6) - 0.0023 * t**3.54 * np.exp(-t / 4.25) # primary (Birn et al., 2008)
rrf_t1 = -0.79 * t**2.1 * np.exp(-t / 1.6) + 2.66 * t**1.1 * np.exp(-t / 1.6) # time derivative 1 (Chen et al., 2020)
rrf_t2 = -0.069 * t**2.54 * np.exp(-t / 4.25) + 0.0046 * t**3.54 * np.exp(-t / 4.25) # time derivative 2 (Chen et al., 2020)
rrf_d1 = 0.16 * t**3.1 * np.exp(-t / 1.6) # dispersion 1 (Chen et al., 2020)
rrf_d2 = 0.00014 * t**4.54 * np.exp(-t / 4.25) # dispersion 2 (Chen et al., 2020)

def create_physio_basis(hr, rv):
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
        rv (np.ndarray): Respiratory variation signal, a 1D array.

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
    global crf_p, crf_t1, crf_t2, crf_d1, crf_d2, rrf_p, rrf_t1, rrf_t2, rrf_d1, rrf_d2

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

    # Convolve RV with RRF basis
    RV = rv
    RV = RV.flatten()  # Ensure column vector

    RV_conv = np.zeros((len(RV), 5))
    RV_conv[:, 0] = convolve(RV, rrf_p, mode='full')[:len(RV)]
    RV_conv[:, 1] = convolve(RV, rrf_t1, mode='full')[:len(RV)]
    RV_conv[:, 2] = convolve(RV, rrf_t2, mode='full')[:len(RV)]
    RV_conv[:, 3] = convolve(RV, rrf_d1, mode='full')[:len(RV)]
    RV_conv[:, 4] = convolve(RV, rrf_d2, mode='full')[:len(RV)]
    RV_basis_regs = RV_conv

    # Concatenate etco2, CO2, HR, and RV basis functions
    X_physio = np.hstack([HR_basis_regs, RV_basis_regs])

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

def calculate_breathing_rate_new(resp_signal, fs):
    # Process the respiration signal
    rsp_signals, _ = nk.rsp_process(resp_signal, sampling_rate=fs)

    # Extract the mean breathing rate
    return rsp_signals["RSP_Rate"].mean()

# Function to load and calculate HRV metrics and breathing rate
def calc_physio_features(subs_id):
    for i, name in enumerate(subs_id):
        if os.path.exists(f"/data1/neurdylab/datasets/nki_rockland/preproc_physio/{name}/{name}_ses-BAS1_task-rest_acq-1400_physio_physOUT.mat"):
            hrv_data = loadmat(f"/data1/neurdylab/datasets/nki_rockland/preproc_physio/{name}/{name}_ses-BAS1_task-rest_acq-1400_physio_physOUT.mat")
        else:
            hrv_data = loadmat(f"/data1/neurdylab/datasets/nki_rockland/preproc_physio_QA_bad/{name}/{name}_ses-BAS1_task-rest_acq-1400_physio_physOUT.mat")

        #hrv_data = loadmat(f"/data1/neurdylab/datasets/nki_rockland/preproc_physio/{name}/{name}_ses-BAS1_task-rest_acq-1400_physio_physOUT.mat")
        ibi = hrv_data['OUT_p']['IBI_clean'][0][0].flatten()
        ibi = ibi[~np.isnan(ibi)]
        ibi_ms = ibi * 1000  # Convert to milliseconds
        rmssd = np.log(np.sqrt(np.mean(np.diff(ibi_ms) ** 2)))

        data = loadmat(f"/data1/neurdylab/datasets/nki_rockland/preproc_physio2/{name}/{name}_ses-BAS1_task-rest_acq-1400_physio_physOUT.mat")

        hr = data['REGS']['hr'][0][0].flatten()
        rv = data['REGS']['rv'][0][0].flatten()
        resp = data['OUT_p']['resp'][0][0][0][0][1].flatten()
        breathing_rate = calculate_breathing_rate_new(resp, fs=62.5)

        sd_rv = np.std(rv)
        avg_hr = np.mean(hr)

        # Calculate power spectral density for LF and HF
        f, psd = welch(ibi_ms, fs=1 / np.mean(ibi), nperseg=len(ibi))
        lf_band = (f >= 0.04) & (f < 0.15)
        hf_band = (f >= 0.15) & (f < 0.4)

        lf = np.log(np.trapz(psd[lf_band], f[lf_band]))
        hf = np.log(np.trapz(psd[hf_band], f[hf_band]))
        hr = detrend(hr)
        rv = detrend(rv)
        print(f"Processed Physio {name}")
        return hr, rv, lf, hf, avg_hr, sd_rv, breathing_rate

def calc_percent_variance_explained(name, mask_indices): 
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
    
    hr, rv, lf, hf, avg_hr, sd_rv, br = calc_physio_features(name) 
    nii_path = f'/data1/neurdylab/datasets/nki_rockland/proc/{name}/ses-BAS1/func/ants_out/{name}_ses-BAS1_task-rest_acq-1400_bold_mo_EPI2MNI_sm_nr.nii.gz'

    nn, _ = load_nii(nii_path)

    # Detrend physiological signals to remove linear trends
    co2 = detrend(co2)
    hr = detrend(hr)
    rv = detrend(rv)

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

    X = create_physio_basis(hr, rv, 1.4) # X is the design matrix for physiological regressors  

    # HRRV only
    X5 = np.hstack([np.ones((X.shape[0], 1)), X[:, 4:]])
    B5 = pinv(X5) @ Y_clean.T  # Least-squares regression coefficients
    percent_variance_hrrv = np.var((X5 @ B5).T, axis=1, keepdims=True) / np.var(Y_clean, axis=1, keepdims=True)

    return percent_variance_hrrv, lf, hf, avg_hr, sd_rv, br 

# Function to process data for a single participant
def process_participant(i, subs_id, mask_indices):
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
    hrrv, lf, hf, avg_hr, sd_rv, br = calc_percent_variance_explained(subs_id[i], mask_indices)
    return i, hrrv, lf, hf, avg_hr, sd_rv, br

# Function to create a design matrix
def create_design_matrix(df):
    """
    Create a design matrix for group comparison (young vs. old).

    Parameters:
        df (pd.DataFrame): DataFrame containing participant information, including age.

    Returns:
        np.ndarray: Design matrix with one column for "young" and another for "old".
    """
    young = np.array(df['age'] < 50, dtype=int)
    old = np.array(df['age'] >= 50, dtype=int)
    gender = np.array(df['gender'] == 'M', dtype=int)
    return np.concatenate((young[:, np.newaxis], old[:, np.newaxis]), axis=1)

# Function to create a contrast matrix
def create_contrast_matrix():
    """
    Create a contrast matrix for group comparison (young vs. old).

    Returns:
        np.ndarray: Contrast matrix for statistical analysis.
    """
    return np.array([[1, -1], [-1, 1]])

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
    df = np.read_csv('../metadata/nki_age_gender_v2.csv')
    subs = df['subs_id'].values
    age = df['age'].values
    gender = df['gender'].values
    N = len(subs)

    mask, affine = load_nii("../metadata/nki_age_gender_v2.csv")
    mask = mask.astype(bool)
    mask_indices = np.where(mask.flatten())[0]

    # Initialize result arrays
    hrrv_cov = np.zeros((91, 109, 91, N))

    # Parallel processing of participants
    num_cores = 16
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cores) as executor:
        try: 
            futures = [executor.submit(process_participant, i, subs, mask_indices) for i in range(N)]
            for future in concurrent.futures.as_completed(futures):
                index, hrrv_cov_update, lf, hf, avg_hr, sd_rv, br = future.result()

                # Reshape results back into brain space
                hrrv_cov[:, :, :, index] = update_brain_map(hrrv_cov_update, mask_indices)

                print(f"Processing complete for participant {subs[index]}")
        except Exception as e:
            print(f"Error: {e}")

    output_dir = 'data/nki_pve_results'
    os.makedirs(output_dir, exist_ok=True)

    # Save variance explained results as NIfTI files
    save_nifti(hrrv_cov, affine, os.path.join(output_dir, 'hrrv_cov_young_old.nii.gz'))

    # Generate and save design and contrast matrices
    design = create_design_matrix(df)
    np.savetxt(os.path.join(output_dir, 'design_matrix.txt'), design, fmt='%d')

    contrast = create_contrast_matrix()
    np.savetxt(os.path.join(output_dir, 'contrast_matrix.txt'), contrast, fmt='%d')

    os.system('bash /path/to/randomise_young_old_baseline.sh') # Run randomise script for group comparison

if __name__ == "__main__":
    main()