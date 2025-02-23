import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.io import loadmat
from scipy.signal import find_peaks, butter, filtfilt, welch
import neurokit2 as nk
import warnings
import sys 
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from utils.file_paths_hrver import get_file_paths
warnings.filterwarnings('ignore')

# Example usage:
metrics = ['SDNN', 'RMSSD', 'LF', 'HF', 'Avg HR', 'Avg CO2', 'Breathing Rate']
titles = ['SDNN', 'RMSSD', 'LF Power', 'HF Power', 'Avg Heart Rate', 'Avg CO2', 'Breathing Rate']

def add_stat_annotation(ax, x1, x2, y, h, p):
    y_padding = 0.15 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    y = y + y_padding
    ax.plot([x1, x2], [y+h, y+h], lw=1.5, c='black')
    stars = '*' if p < 0.05 else '**' if p < 0.01 else '***' if p < 0.001 else '****' if p < 0.0001 else 'n.s.'
    ax.text((x1+x2)*.5, y+h, stars, ha='center', va='bottom', fontsize=26, fontweight='bold')

def load_data():
    # Load participant list for ses-pre and ses-post
    infile_pre = "data/hrver_ses_pre_age_gender.csv"
    infile_post = "data/hrver_ses_post_age_gender.csv"
    df_pre = pd.read_csv(infile_pre)
    df_post = pd.read_csv(infile_post)

    # Identify common participants in both sessions
    common_subs = set(df_pre.iloc[:, 0]).intersection(set(df_post.iloc[:, 0]))
    common_subs = list(common_subs)

    # Filter for common subjects
    df_pre = df_pre[df_pre.iloc[:, 0].isin(common_subs)]
    df_post = df_post[df_post.iloc[:, 0].isin(common_subs)]
    subs_id = df_pre.iloc[:, 0].values
    n_subs = len(subs_id)

    # Extract age and categorize for common subjects
    age = np.array([int(a) for a in df_pre.iloc[:, 2].values])
    age_cat = ['old' if a >= 50 else 'young' for a in age]

    # Initialize arrays for metrics
    results_pre = {metric: np.zeros(n_subs) for metric in metrics}
    results_post = {metric: np.zeros(n_subs) for metric in metrics}

    # Function to load and calculate HRV metrics and breathing rate
    def load_hrv_data(subs_id, results, session):
        for i, name in enumerate(subs_id):
            data = loadmat(get_file_paths('pre', name)['physio'])
            ibi = data['OUT_p']['IBI_clean'][0][0].flatten()
            ibi = ibi[~np.isnan(ibi)]
            ibi_ms = [i*1000 for i in ibi]
            results['SDNN'][i] = np.log(np.std(ibi_ms))
            results['RMSSD'][i] = np.log(np.sqrt(np.mean(np.diff(ibi_ms) ** 2)))

            hr = data['REGS']['hr'][0][0].flatten()
            co2 = data['RESP']['etco2'][0][0].flatten()
            cng_ds = data['RESP']['cng_ds'][0][0].flatten()
            results['Avg HR'][i] = np.mean(hr)
            results['Avg CO2'][i] = np.mean(co2)

            # Calculate power spectral density for LF and HF
            f, psd = welch(ibi_ms, fs=1 / np.mean(ibi), nperseg=len(ibi))
            lf_band = (f >= 0.04) & (f < 0.15)
            hf_band = (f >= 0.15) & (f < 0.4)

            results['LF'][i] = np.log(np.trapz(psd[lf_band], f[lf_band]))
            results['HF'][i] = np.log(np.trapz(psd[hf_band], f[hf_band]))

            resp_signals, _ = nk.rsp_process(cng_ds, sampling_rate=1000)
            results['Breathing Rate'][i] = resp_signals['RSP_Rate'].mean()
            print(f"Subject {name} processed for {session} session")

    # Load HRV metrics and breathing rate for pre session
    load_hrv_data(subs_id, results_pre, 'pre')
    load_hrv_data(subs_id, results_post, 'post')

    # Create DataFrame for common subjects
    df_hrv_pre = pd.DataFrame({
        'Subject': subs_id,
        'OSC': df_pre['osc'].values,
        'Age Category': age_cat,
        'SDNN': results_pre['SDNN'],
        'RMSSD': results_pre['RMSSD'],
        'LF': results_pre['LF'],
        'HF': results_pre['HF'],
        'Avg CO2': results_pre['Avg CO2'],
        'Avg HR': results_pre['Avg HR'],
        'Breathing Rate': results_pre['Breathing Rate']
    })

    df_hrv_post = pd.DataFrame({
        'Subject': subs_id,
        'OSC': df_post['osc'].values,
        'Age Category': age_cat,
        'SDNN': results_post['SDNN'],
        'RMSSD': results_post['RMSSD'],
        'LF': results_post['LF'],
        'HF': results_post['HF'],
        'Avg CO2': results_post['Avg CO2'],
        'Avg HR': results_post['Avg HR'],
        'Breathing Rate': results_post['Breathing Rate']
    })

    df_hrv_average = pd.DataFrame({
        'Subject': subs_id,
        'Age Category': age_cat,
        'SDNN': (results_pre['SDNN'] + results_post['SDNN']) / 2,
        'RMSSD': (results_pre['RMSSD'] + results_post['RMSSD']) / 2,
        'LF': (results_pre['LF'] + results_post['LF']) / 2,
        'HF': (results_pre['HF'] + results_post['HF']) / 2,
        'Avg CO2': (results_pre['Avg CO2'] + results_post['Avg CO2']) / 2,
        'Avg HR': (results_pre['Avg HR'] + results_post['Avg HR']) / 2,
        'Breathing Rate': (results_pre['Breathing Rate'] + results_post['Breathing Rate']) / 2
    })

    return df_hrv_pre, df_hrv_post, df_hrv_average

def create_old_young_violins(df):
    sns.set_theme(style='white')
    fig, axes = plt.subplots(2, 3, figsize=(22, 14))  # Slightly larger figure size
    plt.rcParams['font.family'] = 'Arial' 
    plt.rcParams['font.sans-serif'] = ['Arial']

    palette = 'coolwarm'
    opacity = 0.5
    metrics_hsv = ['RMSSD', 'LF', 'HF']
    metrics_other = ['Avg HR', 'Avg CO2', 'Breathing Rate']
    titles_hsv = ['ln RMSSD', 'ln LF Power', 'ln HF Power']
    titles_other = ['Avg Heart Rate', 'Avg CO2', 'Breathing Rate']

    # Plot HRV metrics
    for i, (metric, title) in enumerate(zip(metrics_hsv + metrics_other, titles_hsv + titles_other)):
        ax = axes[i // 3, i % 3]
        sns.violinplot(x='Age Category', y=metric, data=df, palette=palette, ax=ax, inner=None)
        sns.boxplot(x='Age Category', y=metric, data=df, palette=palette, ax=ax, width=0.3, showfliers=False, linewidth=2)
        sns.stripplot(x='Age Category', y=metric, data=df, color='black', jitter=True, ax=ax, size=5, alpha=opacity)

        # Set plot title and calculate statistical significance
        ax.set_title(title)
        old_data = df[df['Age Category'] == 'old'][metric]
        young_data = df[df['Age Category'] == 'young'][metric]
        ttest = stats.ttest_ind(old_data, young_data)

        ymax = max(old_data.max(), young_data.max())
        if ttest.pvalue < 0.05:
            add_stat_annotation(ax, 0, 1, ymax, 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0]), ttest.pvalue)
        
        # increase the y-axis limit slightly for better visualization
        ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] + 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0]))

        # Customize y-axis ticks and labels
        ax.set_ylabel(title, fontsize=28, fontweight='bold', family='Arial')
        ax.tick_params(axis='y', labelsize=24)  # Increase y-axis tick size
        ax.tick_params(axis='x', labelsize=24)  # Increase x-axis tick size
        ax.set_xlabel("")
        ax.set_xticklabels(['Young', 'Old'], fontsize=28, fontweight='bold', family='Arial')
        ax.set_title("", fontsize=15)

    # Adjust spacing between subplots for clarity
    plt.tight_layout()
    fig.subplots_adjust(wspace=0.3, hspace=0.3)  # Increase spacing between plots
    plt.savefig('results/hrver/physio_stats_hrver/old_young_violins.png', dpi=300)

def create_summary_table(df_pre, df_post, metrics):
    from scipy import stats
    import numpy as np
    
    # Create MultiIndex columns for better organization
    groups = ['Young OSC+', 'Young OSC-', 'Old OSC+', 'Old OSC-']
    sessions = ['Pre', 'Post']
    
    # Create column MultiIndex
    columns = pd.MultiIndex.from_product([
        groups,
        sessions
    ], names=['Group', 'Session'])
    
    # Get sample sizes for each group
    sample_sizes = {}
    for group in ['young_osc+', 'young_osc-', 'old_osc+', 'old_osc-']:
        pre_count = ((df_pre['Age Category'] == group.split('_')[0]) & 
                    (df_pre['OSC'] == group.split('_')[1])).sum()
        post_count = ((df_post['Age Category'] == group.split('_')[0]) & 
                     (df_post['OSC'] == group.split('_')[1])).sum()
        sample_sizes[group] = (pre_count, post_count)
    
    # Initialize data dictionary and statistical results
    data = []
    within_group_stats = []
    between_group_stats = []
    
    # Fill data and calculate statistics for each metric
    for metric in metrics:
        metric_data = []
        metric_within_stats = []
        metric_between_stats = {'Young vs Old OSC+': {}, 'Young vs Old OSC-': {}}
        
        for group in ['young_osc+', 'young_osc-', 'old_osc+', 'old_osc-']:
            age = group.split('_')[0]
            osc = group.split('_')[1]
            
            # Get pre and post data for this group
            mask_pre = (df_pre['Age Category'] == age) & (df_pre['OSC'] == osc)
            mask_post = (df_post['Age Category'] == age) & (df_post['OSC'] == osc)
            
            pre_data = df_pre[mask_pre][metric]
            post_data = df_post[mask_post][metric]
            
            # Calculate means and variances
            mean_pre = pre_data.mean()
            var_pre = pre_data.var()
            mean_post = post_data.mean()
            var_post = post_data.var()
            
            # Perform paired t-test for pre vs post
            t_stat, p_val = stats.ttest_rel(pre_data, post_data)
            
            metric_data.extend([
                f'{mean_pre:.2f} ± {var_pre:.2f}',  # Changed from SEM to variance
                f'{mean_post:.2f} ± {var_post:.2f}'  # Changed from SEM to variance
            ])
            
            # Store within-group statistical results
            metric_within_stats.append({
                'group': group,
                't_stat': t_stat,
                'p_value': p_val
            })
            
            # Prepare data for between-group comparisons
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
    
    # Create main summary DataFrame
    df_summary = pd.DataFrame(
        data,
        index=metrics,
        columns=columns
    )
    
    # Create statistical results DataFrames
    within_group_results = []
    for metric_idx, metric in enumerate(metrics):
        for stat in within_group_stats[metric_idx]:
            within_group_results.append({
                'Metric': metric,
                'Group': stat['group'],
                't-statistic': f"{stat['t_stat']:.3f}",
                'p-value': f"{stat['p_value']:.3f}"
            })
    
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
    
    # Create flat version for CSV export
    new_columns = []
    for group in groups:
        group_key = group.lower().replace(' ', '_')
        pre_n, post_n = sample_sizes[group_key]
        new_columns.extend([f'{group} (n={pre_n}/{post_n}) Pre', 
                          f'{group} (n={pre_n}/{post_n}) Post'])
    
    df_summary_flat = pd.DataFrame(
        data,
        index=metrics,
        columns=new_columns
    )
    
    # Save results to CSV files
    df_summary_flat.to_csv('results/hrver/physio_stats_hrver/hrv_summary_statistics.csv')
    df_within_stats.to_csv('results/hrver/physio_stats_hrver/hrv_within_group_statistics.csv')
    df_between_stats.to_csv('results/hrver/physio_stats_hrver/hrv_between_group_statistics.csv')
    
    return df_summary, df_within_stats, df_between_stats

# Create and display the summary table
os.makedirs('results/hrver/physio_stats_hrver', exist_ok=True)  
metrics = ['SDNN', 'RMSSD', 'LF', 'HF', 'Avg HR', 'Avg CO2', 'Breathing Rate']
df_hrv_pre, df_hrv_post, df_hrv_average = load_data()
summary_table, within_stats, between_stats = create_summary_table(df_hrv_pre, df_hrv_post, metrics)

# View the statistical results
print("Within-group statistics (Pre vs Post):")
print(within_stats)
print("\nBetween-group statistics (Young vs Old):")
print(between_stats)

create_old_young_violins(df_hrv_pre)
