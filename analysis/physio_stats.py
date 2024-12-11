# Calculates and compares physio metrics between young and old participants (e.g., HRV, breathing rate) (Figure 2). 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.io import loadmat
from scipy.signal import find_peaks, butter, filtfilt, welch
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Define key metrics for heart rate variability (HRV) and other physiological measures
metrics = ['SDNN', 'RMSSD', 'LF', 'HF', 'Avg HR', 'Avg CO2', 'Breathing Rate']
titles = ['SDNN', 'RMSSD', 'LF Power', 'HF Power', 'Avg Heart Rate', 'Avg CO2', 'Breathing Rate']

def add_stat_annotation(ax, x1, x2, y, h, p):
    """
    Adds a statistical annotation (e.g., p-values and significance stars) to a plot.

    Parameters:
        ax (matplotlib.axes.Axes): The plot axis to annotate.
        x1, x2 (float): The x-coordinates of the groups being compared.
        y (float): The y-coordinate for placing the annotation line.
        h (float): Height of the annotation line.
        p (float): The p-value of the comparison.
    """
    y_padding = 0.15 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    y = y + y_padding
    ax.plot([x1, x2], [y + h, y + h], lw=1.5, c='black')
    stars = '*' if p < 0.05 else '**' if p < 0.01 else '***' if p < 0.001 else '****' if p < 0.0001 else 'n.s.'
    ax.text((x1 + x2) * .5, y + h, stars, ha='center', va='bottom', fontsize=15)

def load_data():
    """
    Loads and processes HRV data for both pre and post sessions.

    Returns:
        tuple: DataFrames for pre-session data, post-session data, and average data across sessions.
    """
    # Load pre- and post-session data from CSV files
    infile_pre = "/path/to/ses_pre_age_gender.csv"
    infile_post = "/path/to/ses_post_age_gender.csv"
    df_pre = pd.read_csv(infile_pre)
    df_post = pd.read_csv(infile_post)

    # Identify and filter for participants common to both sessions
    common_subs = set(df_pre.iloc[:, 0]).intersection(set(df_post.iloc[:, 0]))
    df_pre = df_pre[df_pre.iloc[:, 0].isin(common_subs)]
    df_post = df_post[df_post.iloc[:, 0].isin(common_subs)]
    subs_id = df_pre.iloc[:, 0].values
    n_subs = len(subs_id)

    # Categorize participants by age
    age = np.array([int(a) for a in df_pre.iloc[:, 2].values])
    age_cat = ['old' if a >= 50 else 'young' for a in age]

    # Initialize dictionaries for storing metric results
    results_pre = {metric: np.zeros(n_subs) for metric in metrics}
    results_post = {metric: np.zeros(n_subs) for metric in metrics}

    def calc_physio(subs_id, results, session):
        """
        Loads and calculates HRV metrics and breathing rate for each participant.

        Parameters:
            subs_id (list): List of subject IDs.
            results (dict): Dictionary to store calculated metrics.
            session (str): Session identifier ('pre' or 'post').
        """
        for i, name in enumerate(subs_id):
            # Load physiological data
            data = loadmat(f'/path/to/{name}_{session}/physio.mat')
            ibi = data['OUT_p']['IBI_clean'][0][0].flatten()
            ibi = ibi[~np.isnan(ibi)]
            ibi_ms = [i * 1000 for i in ibi]

            # Calculate HRV metrics
            results['SDNN'][i] = np.log(np.std(ibi_ms))
            results['RMSSD'][i] = np.log(np.sqrt(np.mean(np.diff(ibi_ms) ** 2)))

            # Extract additional metrics
            hr = data['REGS']['hr'][0][0].flatten()
            co2 = data['RESP']['etco2'][0][0].flatten()
            cng_ds = data['RESP']['cng_ds'][0][0].flatten()
            results['Avg HR'][i] = np.mean(hr)
            results['Avg CO2'][i] = np.mean(co2)

            # Power spectral density for LF and HF bands
            f, psd = welch(ibi_ms, fs=1 / np.mean(ibi), nperseg=len(ibi))
            lf_band = (f >= 0.04) & (f < 0.15)
            hf_band = (f >= 0.15) & (f < 0.4)
            results['LF'][i] = np.log(np.trapz(psd[lf_band], f[lf_band]))
            results['HF'][i] = np.log(np.trapz(psd[hf_band], f[hf_band]))

            # Calculate breathing rate from respiratory signal
            def lowpass_filter(signal, cutoff=0.5, fs=1000, order=4):
                """
                Applies a low-pass filter to a signal using a Butterworth filter to get rid of high-frequency noise.

                Parameters:
                    signal (ndarray): Input signal.
                    cutoff (float): Cutoff frequency in Hz.
                    fs (int): Sampling frequency in Hz.
                    order (int): Filter order.

                Returns:
                    ndarray: Filtered signal.
                """
                nyquist = 0.5 * fs # Assumes 30 breaths for minute as the maximum breathing rate
                normal_cutoff = cutoff / nyquist
                b, a = butter(order, normal_cutoff, btype='low', analog=False)
                return filtfilt(b, a, signal)

            # Filter and calculate breathing rate
            filtered_resp = lowpass_filter(cng_ds)
            peaks, _ = find_peaks(filtered_resp, distance=int(1000 * 2))
            peak_intervals = np.diff(peaks) / 1000
            results['Breathing Rate'][i] = 60 / np.mean(peak_intervals)

    # Calculate physio metrics for pre and post sessions
    calc_physio(subs_id, results_pre, 'pre')
    calc_physio(subs_id, results_post, 'post')

    # Create DataFrames for pre, post, and average sessions
    df_hrv_pre = pd.DataFrame({
        'Subject': subs_id,
        'OSC': df_pre['osc'].values,
        'Age Category': age_cat,
        **results_pre
    })

    df_hrv_post = pd.DataFrame({
        'Subject': subs_id,
        'OSC': df_post['osc'].values,
        'Age Category': age_cat,
        **results_post
    })

    df_hrv_average = pd.DataFrame({
        'Subject': subs_id,
        'Age Category': age_cat,
        **{metric: (results_pre[metric] + results_post[metric]) / 2 for metric in metrics}
    })

    return df_hrv_pre, df_hrv_post, df_hrv_average

def create_old_young_violins(df):
    """
    Creates violin plots comparing HRV metrics between young and old participants.

    Parameters:
        df (DataFrame): DataFrame containing HRV metrics and age categories.
    """
    sns.set_theme(style='white')
    fig, axes = plt.subplots(2, 3, figsize=(22, 14))
    metrics_hsv = ['RMSSD', 'LF', 'HF']
    metrics_other = ['Avg HR', 'Avg CO2', 'Breathing Rate']
    titles_hsv = ['ln RMSSD', 'ln LF Power', 'ln HF Power']
    titles_other = ['Avg Heart Rate', 'Avg CO2', 'Breathing Rate']

    for i, (metric, title) in enumerate(zip(metrics_hsv + metrics_other, titles_hsv + titles_other)):
        ax = axes[i // 3, i % 3]
        sns.violinplot(x='Age Category', y=metric, data=df, ax=ax, inner=None)
        sns.boxplot(x='Age Category', y=metric, data=df, width=0.3, showfliers=False, linewidth=2)
        sns.stripplot(x='Age Category', y=metric, data=df, jitter=True, ax=ax, size=5, alpha=0.5)

        # Add statistical annotations
        old_data = df[df['Age Category'] == 'old'][metric]
        young_data = df[df['Age Category'] == 'young'][metric]
        ttest = stats.ttest_ind(old_data, young_data)
        if ttest.pvalue < 0.05:
            add_stat_annotation(ax, 0, 1, max(old_data.max(), young_data.max()), 0.1, ttest.pvalue)

    plt.tight_layout()
    plt.show()

def create_summary_table(df_pre, df_post, metrics):
    """
    Creates a summary table of HRV metrics for pre and post sessions, 
    including within-group and between-group statistical comparisons.

    Parameters:
        df_pre (DataFrame): DataFrame containing pre-session data.
        df_post (DataFrame): DataFrame containing post-session data.
        metrics (list): List of HRV metrics to analyze.

    Returns:
        tuple:
            - DataFrame: Summary statistics for pre and post sessions across groups.
            - DataFrame: Within-group t-test statistics for pre vs. post comparisons.
            - DataFrame: Between-group t-test statistics for young vs. old comparisons.
    """
    # Define groups and sessions
    groups = ['Young OSC+', 'Young OSC-', 'Old OSC+', 'Old OSC-']
    sessions = ['Pre', 'Post']

    # Create a MultiIndex for the summary table columns
    columns = pd.MultiIndex.from_product([groups, sessions], names=['Group', 'Session'])

    # Determine sample sizes for each group and session
    sample_sizes = {}
    for group in ['young_osc+', 'young_osc-', 'old_osc+', 'old_osc-']:
        pre_count = ((df_pre['Age Category'] == group.split('_')[0]) & 
                     (df_pre['OSC'] == group.split('_')[1])).sum()
        post_count = ((df_post['Age Category'] == group.split('_')[0]) & 
                      (df_post['OSC'] == group.split('_')[1])).sum()
        sample_sizes[group] = (pre_count, post_count)

    # Initialize containers for summary data and statistical results
    data = []
    within_group_stats = []
    between_group_stats = []

    # Process each metric
    for metric in metrics:
        metric_data = []
        metric_within_stats = []
        metric_between_stats = {'Young vs Old OSC+': {}, 'Young vs Old OSC-': {}}

        for group in ['young_osc+', 'young_osc-', 'old_osc+', 'old_osc-']:
            age = group.split('_')[0]
            osc = group.split('_')[1]

            # Filter data for this group
            mask_pre = (df_pre['Age Category'] == age) & (df_pre['OSC'] == osc)
            mask_post = (df_post['Age Category'] == age) & (df_post['OSC'] == osc)

            pre_data = df_pre[mask_pre][metric]
            post_data = df_post[mask_post][metric]

            # Calculate mean and variance for pre and post sessions
            mean_pre = pre_data.mean()
            var_pre = pre_data.var()
            mean_post = post_data.mean()
            var_post = post_data.var()

            # Perform paired t-test for pre vs. post comparison
            t_stat, p_val = stats.ttest_rel(pre_data, post_data)

            # Append results to metric data
            metric_data.extend([
                f'{mean_pre:.2f} ± {var_pre:.2f}',
                f'{mean_post:.2f} ± {var_post:.2f}'
            ])

            # Store within-group statistics
            metric_within_stats.append({
                'group': group,
                't_stat': t_stat,
                'p_value': p_val
            })

            # Prepare for between-group comparisons
            if age == 'young':
                young_data = np.concatenate([pre_data])
                young_group = osc
            elif age == 'old':
                old_data = np.concatenate([pre_data])
                # Perform between-group t-test
                t_stat, p_val = stats.ttest_ind(young_data, old_data)
                metric_between_stats[f'Young vs Old {osc.upper()}'] = {
                    't_stat': t_stat,
                    'p_value': p_val
                }

        data.append(metric_data)
        within_group_stats.append(metric_within_stats)
        between_group_stats.append(metric_between_stats)

    # Create the summary DataFrame
    df_summary = pd.DataFrame(
        data,
        index=metrics,
        columns=columns
    )

    # Process within-group statistics
    within_group_results = []
    for metric_idx, metric in enumerate(metrics):
        for stat in within_group_stats[metric_idx]:
            within_group_results.append({
                'Metric': metric,
                'Group': stat['group'],
                't-statistic': f"{stat['t_stat']:.3f}",
                'p-value': f"{stat['p_value']:.3f}"
            })

    # Process between-group statistics
    between_group_results = []
    for metric_idx, metric in enumerate(metrics):
        for comparison, stats in between_group_stats[metric_idx].items():
            between_group_results.append({
                'Metric': metric,
                'Comparison': comparison,
                't-statistic': f"{stats['t_stat']:.3f}",
                'p-value': f"{stats['p_value']:.3f}"
            })

    df_within_stats = pd.DataFrame(within_group_results)
    df_between_stats = pd.DataFrame(between_group_results)

    # Flatten summary DataFrame for CSV export
    new_columns = []
    for group in groups:
        group_key = group.lower().replace(' ', '_')
        pre_n, post_n = sample_sizes[group_key]
        new_columns.extend([
            f'{group} (n={pre_n}/{post_n}) Pre', 
            f'{group} (n={pre_n}/{post_n}) Post'
        ])

    df_summary_flat = pd.DataFrame(
        data,
        index=metrics,
        columns=new_columns
    )

    # Save the summary and statistical results to CSV
    df_summary_flat.to_csv('hrv_summary_statistics.csv')
    df_within_stats.to_csv('hrv_within_group_statistics.csv')
    df_between_stats.to_csv('hrv_between_group_statistics.csv')

    return df_summary, df_within_stats, df_between_stats


# Load and process physio data
df_hrv_pre, df_hrv_post, df_hrv_average = load_data()
create_old_young_violins(df_hrv_pre)
