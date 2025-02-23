import os
import numpy as np
import nibabel as nib
import pandas as pd
from scipy.io import loadmat
from scipy.signal import detrend
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import seaborn as sns
from concurrent.futures import ProcessPoolExecutor
import concurrent.futures
from scipy import stats
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from utils.file_paths_hrver import get_file_paths
import warnings
warnings.filterwarnings("ignore")

# Set random seed at the start and configure parameters of interest
np.random.seed(42)
tissue = "gray_matter"
age_group = "older"
osc_condition = "osc+"

def get_age_group(age):
    """Categorize participants into age groups."""
    if age <= 36:
        return "young"
    else:
        return "older"

def resample_signal(signal, old_fs, new_fs):
    """
    Resample signal from old sampling frequency to new sampling frequency.
    
    Parameters:
        signal (np.ndarray): Original signal
        old_fs (float): Original sampling frequency in Hz
        new_fs (float): Target sampling frequency in Hz
    
    Returns:
        np.ndarray: Resampled signal
    """
    # Create time vectors
    old_time = np.arange(len(signal)) / old_fs
    new_time = np.arange(0, old_time[-1], 1/new_fs)
    
    # Create interpolation function
    f = interp1d(old_time, signal, kind='cubic', fill_value='extrapolate')
    
    # Resample signal
    resampled_signal = f(new_time)
    
    return resampled_signal

def extract_gray_matter_timeseries(bold_data, schaefer_mask):
    """Extract average BOLD timeseries for all gray matter voxels."""
    gm_mask = schaefer_mask > 0
    gm_timeseries = np.mean(bold_data[gm_mask], axis=0)
    return gm_timeseries

def process_participant(sub_id, cond):
    """Process single participant's data with resampling."""
    # Load physiological data
    physio_data = loadmat(get_file_paths(cond, sub_id)['physio']) # Change based on the dataset
    hr = physio_data['REGS']['hr'][0][0].flatten()
    co2 = physio_data['RESP']['etco2'][0][0].flatten()    
    
    # Normalize physiological signals
    hr = detrend(hr)
    co2 = detrend(co2)
    
    # Load BOLD data
    bold_img = nib.load(get_file_paths(cond, sub_id)['bold']) # Change based on the dataset
    bold_data = bold_img.get_fdata()
    
    # Load Mask
    if tissue == "gray_matter":
        mask_dir = "data/masks/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_2mm.nii.gz" # Schaefer 100 parcel mask for gray matter
    elif tissue == "white_matter":
        mask_dir = "data/masks/avg152T1_white_bin90pct.nii" # White matter mask
    elif tissue == "ventricle":
        mask_dir = "data/masks/ventricle_mask_2mm.nii.gz" # Ventricle mask
    mask = nib.load(mask_dir).get_fdata()
    
    # Extract gray matter timeseries
    gm_timeseries = extract_gray_matter_timeseries(bold_data, mask)
    
    # Normalize timeseries
    gm_timeseries = (gm_timeseries - np.mean(gm_timeseries)) / np.std(gm_timeseries)
    
    # Resample all signals from 2.4s to 0.2s
    old_fs = 1/2.4  # ~0.417 Hz
    new_fs = 1/0.2  # 5 Hz
    
    gm_timeseries_resampled = resample_signal(gm_timeseries, old_fs, new_fs)
    hr_resampled = resample_signal(hr, old_fs, new_fs)
    co2_resampled = resample_signal(co2, old_fs, new_fs)
    
    # Normalize resampled signals
    gm_timeseries_resampled = (gm_timeseries_resampled - np.mean(gm_timeseries_resampled)) / np.std(gm_timeseries_resampled)
    hr_resampled = (hr_resampled - np.mean(hr_resampled)) / np.std(hr_resampled)
    co2_resampled = (co2_resampled - np.mean(co2_resampled)) / np.std(co2_resampled)
    
    return sub_id, cond, gm_timeseries_resampled, hr_resampled, co2_resampled

def compute_cross_correlation(signal1, signal2, max_lag):
    """Compute cross-correlation between two signals for a range of lags."""
    correlations = np.zeros(max_lag + 25)  # -24 to max_lag (2 seconds before at 0.2s resolution)
    
    for lag_idx, lag in enumerate(range(-24, max_lag + 1)):  # Start from -2s (10 samples at 0.2s)
        signal2_lagged = np.roll(signal2, lag)
        correlations[lag_idx] = np.corrcoef(signal1, signal2_lagged)[0, 1]
    
    return correlations

def main():
    # Load participant data
    df_pre = pd.read_csv("data/hrver_ses_pre_age_gender.csv")
    df_post = pd.read_csv("data/hrver_ses_post_age_gender.csv")

    # Filter participants based on age and oscillation status (CHANGE BASED ON ANALYSIS)
    if age_group == "older":
        df_pre = df_pre[(df_pre['age'] > 36) & (df_pre['osc'] == osc_condition)] # Older adults osc+ pre
        df_post = df_post[(df_post['age'] > 36) & (df_post['osc'] == osc_condition)] # Older adults osc+ post
    elif age_group == "younger":
        df_pre = df_pre[(df_pre['age'] <= 36) & (df_pre['osc'] == osc_condition)] # Younger adults osc+ pre
        df_post = df_post[(df_post['age'] <= 36) & (df_post['osc'] == osc_condition)] # Younger adults osc+ post

    df_pre = df_pre[df_pre['sub_id'].isin(df_post['sub_id'])]
    df_post = df_post[df_post['sub_id'].isin(df_pre['sub_id'])]
    
    results_dict = {'pre': {}, 'post': {}}
    with ProcessPoolExecutor(max_workers=16) as executor:
        # Create all futures
        futures = []
        for sub_id in df_pre['sub_id'].values:
            futures.append((executor.submit(process_participant, sub_id, 'pre'), sub_id, 'pre'))
            futures.append((executor.submit(process_participant, sub_id, 'post'), sub_id, 'post'))
        
        # Process results as they complete
        for future, sub_id, cond in futures:
            try:
                result_sub_id, result_cond, gm_ts, hr, co2 = future.result()
                assert result_sub_id == sub_id and result_cond == cond  # Sanity check
                results_dict[cond][sub_id] = (gm_ts, hr, co2)
                print(f"Completed processing participant {sub_id} condition {cond}")
            except Exception as e:
                print(f"Error processing participant {sub_id} condition {cond}: {e}")
    
    # Ensure consistent ordering of subjects
    subject_ids = sorted(set(results_dict['pre'].keys()) & set(results_dict['post'].keys()))
    
    # Organize results by condition with consistent subject ordering
    conditions = {
        'pre': [results_dict['pre'][sub_id] for sub_id in subject_ids],
        'post': [results_dict['post'][sub_id] for sub_id in subject_ids]
    }
    
    # Compute cross-correlations
    max_lag = 144  # 12 seconds at 0.2s resolution
    lags = np.arange(-24, max_lag + 1) * 0.2
    
    correlations = {
        signal_type: {
            cond: np.zeros((len(conditions[cond]), len(lags)))
            for cond in conditions
        }
        for signal_type in ['HR', 'CO2']
    }
    
    # Calculate correlations with consistent subject ordering
    for cond in conditions:
        for i, (gm_ts, hr, co2) in enumerate(conditions[cond]):
            correlations['HR'][cond][i] = compute_cross_correlation(
                gm_ts, hr, max_lag
            )
            correlations['CO2'][cond][i] = compute_cross_correlation(
                gm_ts, co2, max_lag
            )
    
    plot_correlations(correlations, lags)

def compute_ttest_at_lags(correlations, signal_type):
    """Compute t-test between pre and post conditions at each lag."""
    pre_data = correlations[signal_type]['pre']
    post_data = correlations[signal_type]['post']
    
    t_stats = np.zeros(pre_data.shape[1])
    p_values = np.zeros(pre_data.shape[1]) 
    
    for i in range(pre_data.shape[1]):
        t_stat, p_val = stats.ttest_rel(pre_data[:, i], post_data[:, i])
        t_stats[i] = t_stat
        p_values[i] = p_val
    
    return t_stats, p_values

def plot_correlations(correlations, lags):
    """Create plots of cross-correlations with significance testing."""
    signal_types = ['CO2', 'HR']
    conditions = ['pre', 'post']
    colors = {'pre': 'olive', 'post': 'purple'}
    
    plt.rcParams.update({
        'font.size': 18,
        'axes.titlesize': 16,
        'axes.labelsize': 20,
        'xtick.labelsize': 18,
        'ytick.labelsize': 18,
        'legend.fontsize': 18
    })
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    for idx, signal_type in enumerate(signal_types):
        ax = axes[idx]
        
        # Compute t-tests at each lag
        t_stats, p_values = compute_ttest_at_lags(correlations, signal_type)
        significant_bonferroni = p_values < 0.05 / (lags.shape[0]/12)  # Bonferroni corrected
        significant_uncorrected = p_values < 0.05  # Uncorrected threshold
        print(len(lags))
        print(f"{signal_type} p-values: {p_values}")
        
        # Plot correlations for each condition
        for cond in conditions:
            mean_corr = np.mean(correlations[signal_type][cond], axis=0)
            sem_corr = np.std(correlations[signal_type][cond], axis=0) / \
                    np.sqrt(correlations[signal_type][cond].shape[0])
            
            ax.plot(lags, mean_corr, color=colors[cond], 
                label=f'{cond.capitalize()} (n={len(correlations[signal_type][cond])})')
            ax.fill_between(lags, mean_corr - sem_corr, mean_corr + sem_corr,
                        color=colors[cond], alpha=0.2)
        
        # Add shading for significant differences
        ylims = ax.get_ylim()
        # First shade uncorrected significance (lighter gray)
        ax.fill_between(lags, ylims[0], ylims[1], 
                    where=significant_uncorrected,
                    color='gray', alpha=0.09)
        
        # Then shade Bonferroni-corrected significance (darker gray)
        ax.fill_between(lags, ylims[0], ylims[1], 
                    where=significant_bonferroni,
                    color='gray', alpha=0.18)
        
        ax.set_xlabel('Lag (seconds)')
        ax.set_ylabel('Correlation')
        ax.tick_params(axis='both', which='major', labelsize=18)
        ax.set_title(f'{signal_type}-BOLD Cross-Correlation', fontsize=20)
        ax.set_ylim(ylims)
        ax.grid(True)
        ax.legend()
        ax.axvline(x=0, color='black', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

    output_dir = f"results/hrver/cross_corr_avg_pre_post_{age_group}_{osc_condition}"
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(f"{output_dir}/{tissue}.png", dpi=300)

if __name__ == "__main__":
    main()