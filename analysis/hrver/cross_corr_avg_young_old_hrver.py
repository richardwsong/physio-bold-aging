# Old vs Young pre Cross Correlation per tissue (HRV-ER)
import os
import numpy as np
import nibabel as nib
import pandas as pd
from scipy.io import loadmat
from scipy.signal import detrend
from scipy import stats
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker 
import seaborn as sns
from concurrent.futures import ProcessPoolExecutor
import concurrent.futures 
import warnings
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from utils.file_paths_hrver import get_file_paths
warnings.filterwarnings("ignore")

tissue = "gray_matter"

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

def process_participant(sub_id, age):
    """Process single participant's data with resampling."""
    # Load physiological data
    physio_data = loadmat(get_file_paths('pre', sub_id)['physio']) # Change based on the dataset
    hr = physio_data['REGS']['hr'][0][0].flatten()
    co2 = physio_data['RESP']['etco2'][0][0].flatten()    
    
    # Normalize physiological signals
    hr = detrend(hr)
    co2 = detrend(co2)
    
    # Load BOLD data
    bold_img = nib.load(get_file_paths('pre', sub_id)['bold']) # Change based on the dataset
    bold_data = bold_img.get_fdata()
    
    # Load Schaefer atlas
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
    
    return get_age_group(age), gm_timeseries_resampled, hr_resampled, co2_resampled

def compute_cross_correlation(signal1, signal2, max_lag):
    """Compute cross-correlation between two signals for a range of lags."""
    correlations = np.zeros(max_lag + 25)  # -24 to max_lag (2 seconds before at 0.2s resolution)
    
    for lag_idx, lag in enumerate(range(-24, max_lag + 1)):  # Start from -2s (10 samples at 0.2s)
        signal2_lagged = np.roll(signal2, lag)
        correlations[lag_idx] = np.corrcoef(signal1, signal2_lagged)[0, 1]
    
    return correlations

def main():
    # Load participant data
    df = pd.read_csv("data/hrver_ses_pre_age_gender.csv")
    
    # Process all participants in parallel
    results = []
    num_cores = min(os.cpu_count()//2, 16)
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cores) as executor:
        future_to_subj = {executor.submit(process_participant, sub_id, age): (sub_id, age)
                         for sub_id, age in zip(df.iloc[:, 0].values, df.iloc[:, 2].values)}
        
        for future in concurrent.futures.as_completed(future_to_subj):
            sub_id, _ = future_to_subj[future]
            try:
                results.append(future.result())
                print(f"Completed processing participant {sub_id}")
            except Exception as e:
                print(f"Error processing participant {sub_id}: {e}")
    
    # Organize results by age group
    age_groups = {'young': [], 'older': []}
    for age_group, gm_ts, hr, co2 in results:
        age_groups[age_group].append((gm_ts, hr, co2))
    
    # Compute cross-correlations for each age group
    max_lag = 144  # 12 seconds at 0.2s resolution (was 12 TR at 2.4s)
    lags = np.arange(-24, max_lag + 1) * 0.2  # Convert lag indices to seconds
    
    # Initialize correlation arrays
    correlations = {
        signal_type: {
            age_group: np.zeros((len(age_groups[age_group]), len(lags)))
            for age_group in age_groups
        }
        for signal_type in ['HR', 'CO2']
    }
    
    # Calculate correlations
    for age_group in age_groups:
        for i, (gm_ts, hr, co2) in enumerate(age_groups[age_group]):
            correlations['HR'][age_group][i] = compute_cross_correlation(
                gm_ts, hr, max_lag
            )
            correlations['CO2'][age_group][i] = compute_cross_correlation(
                gm_ts, co2, max_lag
            )
    
    # Plot results
    plot_correlations(correlations, lags)

def compute_ttest_at_lags(correlations, signal_type):
    """Compute independent t-test between young and older groups at each lag."""
    young_data = correlations[signal_type]['young']
    older_data = correlations[signal_type]['older']
    
    t_stats = np.zeros(young_data.shape[1])
    p_values = np.zeros(young_data.shape[1])
    
    for i in range(young_data.shape[1]):
        t_stat, p_val = stats.ttest_ind(young_data[:, i], older_data[:, i])
        t_stats[i] = t_stat
        p_values[i] = p_val
    
    return t_stats, p_values

def plot_correlations(correlations, lags):
    """Create plots of cross-correlations with two levels of significance testing."""
    # Set global font sizes
    plt.rcParams.update({
        'font.size': 18,
        'axes.titlesize': 16,
        'axes.labelsize': 20,
        'xtick.labelsize': 18,
        'ytick.labelsize': 18,
        'legend.fontsize': 18
    })
    
    signal_types = ['CO2', 'HR']
    age_groups = ['young', 'older']
    colors = {'young': 'blue', 'middle': 'green', 'older': 'red'}
    
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    for idx, signal_type in enumerate(signal_types):
        ax = axes[idx]
        
        # Compute t-tests at each lag
        t_stats, p_values = compute_ttest_at_lags(correlations, signal_type)
        
        # Two levels of significance
        significant_bonferroni = p_values < 0.05/(lags.shape[0]/12)  # Bonferroni corrected
        significant_uncorrected = p_values < 0.05  # Uncorrected
        
        # Calculate y-axis limits from the data first
        all_means = []
        all_sems = []
        for age_group in age_groups:
            mean_corr = np.mean(correlations[signal_type][age_group], axis=0)
            sem_corr = np.std(correlations[signal_type][age_group], axis=0) / \
                      np.sqrt(correlations[signal_type][age_group].shape[0])
            all_means.extend([mean_corr + sem_corr, mean_corr - sem_corr])
        
        ymin, ymax = np.min(all_means), np.max(all_means)
        
        # Add shading for different significance levels
        # First, uncorrected significance (lighter gray)
        ax.fill_between(lags, ymin, ymax,
                       where=significant_uncorrected,
                       color='gray', alpha=0.09)
        
        # Then, Bonferroni-corrected significance (darker gray)
        ax.fill_between(lags, ymin, ymax,
                       where=significant_bonferroni,
                       color='gray', alpha=0.2)
        
        # Plot correlations for each age group
        for age_group in age_groups:
            mean_corr = np.mean(correlations[signal_type][age_group], axis=0)
            sem_corr = np.std(correlations[signal_type][age_group], axis=0) / \
                      np.sqrt(correlations[signal_type][age_group].shape[0])
            
            ax.plot(lags, mean_corr, color=colors[age_group], 
                   label=f'{age_group.capitalize()} (n={len(correlations[signal_type][age_group])})')
            ax.fill_between(lags, mean_corr - sem_corr, mean_corr + sem_corr,
                          color=colors[age_group], alpha=0.2)
        
        ax.set_xlabel('Lag (seconds)')
        ax.set_ylabel('Correlation')
        ax.set_title(f'{signal_type}-BOLD Cross-correlation', fontsize=20, pad=15)
        
        # Set exact y-axis limits
        ax.set_ylim(ymin, ymax)
        ax.set_xlim([-2.8, 30])  # Adjust x-axis limits to remove the gap
        ax.set_xticks(np.array([0, 10, 20, 30]))

        # Increase tick label sizes and make ticks more prominent
        ax.tick_params(axis='both', which='major', width=1.5, length=6)
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f'{y:.2f}'))
        
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right')
        ax.axvline(x=0, color='black', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    # Save the plot
    output_dir = "results/hrver/cross_corr_avg_young_old"
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(f'{output_dir}/{tissue}.png', dpi=300)

if __name__ == "__main__":
    main()