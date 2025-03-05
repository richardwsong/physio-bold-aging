import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.io import loadmat
from scipy.signal import find_peaks, butter, filtfilt, welch
import neurokit2 as nk
import warnings
import concurrent.futures
import sys
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from utils.file_paths_hrver import get_file_paths

warnings.filterwarnings('ignore')

def add_stat_annotation(ax, x1, x2, y, h, p):
    y_padding = 0.15 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    y = y + y_padding
    ax.plot([x1, x2], [y+h, y+h], lw=1.5, c='black')
    stars = '*' if p < 0.05 else '**' if p < 0.01 else '***' if p < 0.001 else '****' if p < 0.0001 else 'n.s.'
    ax.text((x1+x2)*.5, y+h, stars, ha='center', va='bottom', fontsize=26, fontweight='bold')
    
def process_subject(name, session):
    try:
        data = loadmat(get_file_paths(session, name)['physio']) # Path to preprocessed physio data
        # Process IBI data
        ibi = data['OUT_p']['IBI_clean'][0][0].flatten()
        ibi = ibi[~np.isnan(ibi)]
        ibi_ms = [i*1000 for i in ibi]
        
        # Calculate HRV metrics
        sdnn = np.log(np.std(ibi_ms))
        rmssd = np.log(np.sqrt(np.mean(np.diff(ibi_ms) ** 2)))
        
        # Extract other physiological signals
        hr = data['REGS']['hr'][0][0].flatten()
        co2 = data['RESP']['etco2'][0][0].flatten()
        cng_ds = data['RESP']['cng_ds'][0][0].flatten()
        avg_hr = np.mean(hr)
        avg_co2 = np.mean(co2)
        
        # Calculate power spectral density for LF and HF
        f, psd = welch(ibi_ms, fs=1 / np.mean(ibi), nperseg=len(ibi))
        lf_band = (f >= 0.04) & (f < 0.15)
        hf_band = (f >= 0.15) & (f < 0.4)
        
        lf = np.log(np.trapz(psd[lf_band], f[lf_band]))
        hf = np.log(np.trapz(psd[hf_band], f[hf_band]))
        
        # Calculate breathing rate
        resp_signals, _ = nk.rsp_process(cng_ds, sampling_rate=1000)
        breathing_rate = resp_signals['RSP_Rate'].mean()
        
        return {
            'Subject': name,
            'SDNN': sdnn,
            'RMSSD': rmssd,
            'LF': lf,
            'HF': hf,
            'Avg HR': avg_hr,
            'Avg CO2': avg_co2,
            'Breathing Rate': breathing_rate
        }
    except Exception as e:
        print(f"Error processing subject {name} for session {session}: {e}")
        return None

def load_data():
    # Load participant list for ses-pre and ses-post
    infile_pre = "/data1/neurdylab/songrw/derivates/hrv_er/ses_pre_age_gender.csv"
    infile_post = "/data1/neurdylab/songrw/derivates/hrv_er/ses_post_age_gender.csv"
    df_pre = pd.read_csv(infile_pre)
    df_post = pd.read_csv(infile_post)

    # Identify common participants in both sessions
    common_subs = set(df_pre.iloc[:, 0]).intersection(set(df_post.iloc[:, 0]))
    common_subs = list(common_subs)

    # Filter for common subjects
    df_pre_shared = df_pre[df_pre.iloc[:, 0].isin(common_subs)]
    df_post = df_post[df_post.iloc[:, 0].isin(common_subs)]
    subs_id_shared = df_post.iloc[:, 0].values
    
    # Use all pre subjects for pre_only
    subs_id_pre_only = df_pre.iloc[:, 0].values
    
    # Process all pre-session data in parallel
    print(f"Processing pre-session data for {len(subs_id_pre_only)} subjects using ProcessPoolExecutor...")
    pre_results = []
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=16) as executor:
        # Submit all tasks and create a future-to-subject mapping
        future_to_subject = {
            executor.submit(process_subject, name, 'pre'): name 
            for name in subs_id_pre_only
        }
        
        # Process results as they complete
        for future in concurrent.futures.as_completed(future_to_subject):
            subject = future_to_subject[future]
            try:
                result = future.result()
                if result is not None:
                    pre_results.append(result)
                    print(f"Completed processing pre-session subject: {subject}")
                else:
                    print(f"Failed to process pre-session subject: {subject}")
            except Exception as exc:
                print(f"Pre-session subject {subject} generated an exception: {exc}")
    
    print(f"Completed processing {len(pre_results)} of {len(subs_id_pre_only)} pre-session subjects")
    
    # Create pre_only DataFrame from all pre results
    df_hrv_pre_only = pd.DataFrame(pre_results)
    
    # Join with demographic info - using the correct column name from df_pre
    # First, print the columns to verify correct names
    print("df_pre columns:", df_pre.columns.tolist())
    
    # Merge using the correct ID column (first column)
    id_column = df_pre.columns[0]  # Get the actual name of the first column
    df_hrv_pre_only = df_hrv_pre_only.merge(
        df_pre[[id_column, 'osc', 'age']], 
        left_on='Subject', 
        right_on=id_column, 
        how='left'
    )
    
    # Add age category
    df_hrv_pre_only['Age Category'] = df_hrv_pre_only['age'].apply(lambda a: 'old' if int(a) >= 50 else 'young')
    df_hrv_pre_only.rename(columns={'osc': 'OSC'}, inplace=True)
    
    # Drop the redundant ID column if it's not the same as 'Subject'
    if id_column != 'Subject':
        df_hrv_pre_only.drop(columns=[id_column], inplace=True)
    df_hrv_pre_only.drop(columns=['age'], inplace=True)
    
    # Process post-session data only for shared subjects
    print(f"Processing post-session data for {len(subs_id_shared)} subjects using ProcessPoolExecutor...")
    post_results = []
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=16) as executor:
        # Submit all tasks and create a future-to-subject mapping
        future_to_subject = {
            executor.submit(process_subject, name, 'post'): name 
            for name in subs_id_shared
        }
        
        # Process results as they complete
        for future in concurrent.futures.as_completed(future_to_subject):
            subject = future_to_subject[future]
            try:
                result = future.result()
                if result is not None:
                    post_results.append(result)
                    print(f"Completed processing post-session subject: {subject}")
                else:
                    print(f"Failed to process post-session subject: {subject}")
            except Exception as exc:
                print(f"Post-session subject {subject} generated an exception: {exc}")
    
    print(f"Completed processing {len(post_results)} of {len(subs_id_shared)} post-session subjects")
    # Create post DataFrame
    df_hrv_post = pd.DataFrame(post_results)
    
    # Join with demographic info - using the correct column name
    id_column = df_post.columns[0]  # Get the actual name of the first column
    df_hrv_post = df_hrv_post.merge(
        df_post[[id_column, 'osc', 'age']], 
        left_on='Subject', 
        right_on=id_column, 
        how='left'
    )
    
    # Add age category
    df_hrv_post['Age Category'] = df_hrv_post['age'].apply(lambda a: 'old' if int(a) >= 50 else 'young')
    df_hrv_post.rename(columns={'osc': 'OSC'}, inplace=True)
    
    # Drop the redundant ID column if it's not the same as 'Subject'
    if id_column != 'Subject':
        df_hrv_post.drop(columns=[id_column], inplace=True)
    df_hrv_post.drop(columns=['age'], inplace=True)
    
    # Create pre DataFrame using only the subjects also in post (shared)
    df_hrv_pre = df_hrv_pre_only[df_hrv_pre_only['Subject'].isin(df_hrv_post['Subject'])]
    
    return df_hrv_pre, df_hrv_post, df_hrv_pre_only

# Modified function to handle outliers by metric rather than by participant
def handle_metric_outliers(df, columns=None, multiplier=1.5):
    """
    Mark outliers in specific columns using the IQR method, by setting them to NaN
    instead of removing entire rows.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        The DataFrame containing the data
    columns : list or None
        List of column names to check for outliers. If None, all numeric columns are used.
    multiplier : float
        The multiplier for the IQR to determine the outlier threshold (default: 1.5)
        
    Returns:
    --------
    pandas.DataFrame
        A copy of the DataFrame with outliers set to NaN in specific columns
    """
    if columns is None:
        # Use all numeric columns
        columns = df.select_dtypes(include=['number']).columns.tolist()
    
    df_clean = df.copy()
    outlier_counts = {}
    
    for col in columns:
        if col in df.columns and df[col].dtype.kind in 'fc':  # Only process numeric columns
            Q1 = df[col].quantile(0.25)
            Q3 = df[col].quantile(0.75)
            IQR = Q3 - Q1
            
            lower_bound = Q1 - multiplier * IQR
            upper_bound = Q3 + multiplier * IQR
            
            # Find indices of outliers in this column
            outliers = (df[col] < lower_bound) | (df[col] > upper_bound)
            outlier_indices = df.index[outliers]
            outlier_counts[col] = len(outlier_indices)
            
            # Set outliers to NaN in this specific column only
            df_clean.loc[outlier_indices, col] = np.nan
    
    # Print summary of removed outliers by column
    total_cells = sum(outlier_counts.values())
    print(f"Marked {total_cells} outlier values as NaN across {len(columns)} metrics:")
    for col, count in outlier_counts.items():
        if count > 0:
            print(f"  - {col}: {count} outliers ({count/len(df)*100:.1f}%)")
    
    return df_clean

# Updated create_summary_table function to handle NaN values in metrics
# Updated create_summary_table function to ensure proper subject matching for paired t-tests
def create_summary_table(df_pre, df_post, df_pre_only, metrics):
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
    
    # Get sample sizes for each group and metric
    sample_sizes = {}
    for group in ['young_osc+', 'young_osc-', 'old_osc+', 'old_osc-']:
        age, osc = group.split('_')
        
        # Calculate base counts
        pre_mask = (df_pre['Age Category'] == age) & (df_pre['OSC'] == osc)
        post_mask = (df_post['Age Category'] == age) & (df_post['OSC'] == osc)
        pre_only_mask = (df_pre_only['Age Category'] == age) & (df_pre_only['OSC'] == osc)
        
        # Count total participants in each group
        pre_total = pre_mask.sum()
        post_total = post_mask.sum()
        pre_only_total = pre_only_mask.sum()
        
        # Initialize metric-specific sample sizes
        metric_sample_sizes = {}
        for metric in metrics:
            # Count non-NA values for each metric
            pre_count = df_pre.loc[pre_mask, metric].notna().sum()
            post_count = df_post.loc[post_mask, metric].notna().sum()
            pre_only_count = df_pre_only.loc[pre_only_mask, metric].notna().sum()
            metric_sample_sizes[metric] = (pre_count, post_count, pre_only_count)
        
        sample_sizes[group] = {
            'total': (pre_total, post_total, pre_only_total),
            'by_metric': metric_sample_sizes
        }
    
    # Initialize data dictionary and statistical results
    data = []
    within_group_stats = []
    between_group_stats = []
    
    # Fill data and calculate statistics for each metric
    for metric in metrics:
        metric_data = []
        metric_within_stats = []
        metric_between_stats = {'Young vs Old': {}}
        
        for group in ['young_osc+', 'young_osc-', 'old_osc+', 'old_osc-']:
            age, osc = group.split('_')
            
            # Get pre and post data for this group (for within-group)
            mask_pre = (df_pre['Age Category'] == age) & (df_pre['OSC'] == osc)
            mask_post = (df_post['Age Category'] == age) & (df_post['OSC'] == osc)
            
            pre_data_all = df_pre[mask_pre].copy()
            post_data_all = df_post[mask_post].copy()
            
            # Only keep subjects that are in both pre and post
            common_subjects = set(pre_data_all['Subject']).intersection(set(post_data_all['Subject']))
            
            # FIX: Ensure consistent ordering by subject ID before statistical tests
            # Filter to common subjects and sort by Subject ID for consistency
            pre_data = pre_data_all[pre_data_all['Subject'].isin(common_subjects)].sort_values('Subject')
            post_data = post_data_all[post_data_all['Subject'].isin(common_subjects)].sort_values('Subject')
            
            # Verify subjects are in the same order in both dataframes
            assert all(pre_data['Subject'].values == post_data['Subject'].values), "Subject mismatch in pre and post data"
            
            # FIX: Create matched pairs by ensuring subject alignment and then dropping NaNs in either dataset
            paired_data = pd.DataFrame({
                'Subject': pre_data['Subject'].values,
                'pre': pre_data[metric].values,
                'post': post_data[metric].values
            })
            
            # Now we can safely drop NaN values knowing subjects are aligned
            paired_data = paired_data.dropna()
            
            # Calculate means and variances
            if len(paired_data) > 0:
                mean_pre = paired_data['pre'].mean()
                var_pre = paired_data['pre'].var()
                mean_post = paired_data['post'].mean()
                var_post = paired_data['post'].var()
                n_paired = len(paired_data)
            else:
                mean_pre, var_pre = np.nan, np.nan
                mean_post, var_post = np.nan, np.nan
                n_paired = 0
            
            # Perform paired t-test for pre vs post on the valid pairs
            if n_paired > 1:  # Need at least 2 pairs for a paired t-test
                t_stat, p_val = stats.ttest_rel(paired_data['pre'], paired_data['post'])
            else:
                t_stat, p_val = np.nan, np.nan
            
            # Store results using standard deviation (not variance)
            metric_data.extend([
                f'{mean_pre:.2f} ({var_pre:.2f})' if not np.isnan(mean_pre) else 'N/A',
                f'{mean_post:.2f} ({var_post:.2f})' if not np.isnan(mean_post) else 'N/A'
            ])
            
            # Store within-group statistical results
            metric_within_stats.append({
                'group': group,
                't_stat': t_stat,
                'p_value': p_val,
                'n_samples': n_paired
            })
        
        # Get young and old data for between-group analysis from pre_only
        young_data = df_pre_only[df_pre_only['Age Category'] == 'young'][metric].dropna()
        old_data = df_pre_only[df_pre_only['Age Category'] == 'old'][metric].dropna()
        
        # Perform independent t-test
        if len(young_data) > 1 and len(old_data) > 1:  # Need at least 2 samples per group
            t_stat, p_val = stats.ttest_ind(young_data, old_data, equal_var=False)  # Using Welch's t-test
        else:
            t_stat, p_val = np.nan, np.nan
            
        metric_between_stats['Young vs Old'] = {
            't_stat': t_stat,
            'p_value': p_val,
            'n_young': len(young_data),
            'n_old': len(old_data)
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
                't-statistic': f"{stat['t_stat']:.3f}" if not np.isnan(stat['t_stat']) else "N/A",
                'p-value': f"{stat['p_value']:.4f}" if not np.isnan(stat['p_value']) else "N/A",
                'n': stat['n_samples']
            })
    
    between_group_results = []
    for metric_idx, metric in enumerate(metrics):
        for comparison, stats in between_group_stats[metric_idx].items():
            between_group_results.append({
                'Metric': metric,
                'Comparison': comparison,
                't-statistic': f"{stats['t_stat']:.3f}" if not np.isnan(stats['t_stat']) else "N/A",
                'p-value': f"{stats['p_value']:.4f}" if not np.isnan(stats['p_value']) else "N/A",
                'n_young': stats['n_young'],
                'n_old': stats['n_old']
            })
    
    df_within_stats = pd.DataFrame(within_group_results)
    df_between_stats = pd.DataFrame(between_group_results)
    
    # Create flat version for CSV export with metric-specific sample sizes
    flat_columns = []
    for group in groups:
        group_key = group.lower().replace(' ', '_')
        for session_idx, session in enumerate(['Pre', 'Post']):
            metric_specific_counts = []
            for metric in metrics:
                counts = sample_sizes[group_key]['by_metric'][metric]
                metric_specific_counts.append(f"{counts[session_idx]}")
            
            # Add the column with sample size range if needed
            unique_counts = set(metric_specific_counts)
            if len(unique_counts) == 1:
                sample_info = f"(n={next(iter(unique_counts))})"
            else:
                min_count = min(int(c) for c in unique_counts)
                max_count = max(int(c) for c in unique_counts)
                sample_info = f"(n={min_count}-{max_count})"
            
            flat_columns.append(f'{group} {sample_info} {session}')
    
    df_summary_flat = pd.DataFrame(
        data,
        index=metrics,
        columns=flat_columns
    )
    
    return df_summary, df_within_stats, df_between_stats

# Updated violin plot function to handle NaN values
def create_old_young_violins(df):
    sns.set_theme(style='white')
    fig, axes = plt.subplots(2, 3, figsize=(22, 14))
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
        
        # Create a temporary dataframe without NaN values for this metric
        temp_df = df.dropna(subset=[metric]).copy()
        
        # Get actual sample sizes for this metric
        young_count = temp_df[temp_df['Age Category'] == 'young'].shape[0]
        old_count = temp_df[temp_df['Age Category'] == 'old'].shape[0]
        
        # Update the Age Category labels with sample sizes
        temp_df['Age Label'] = temp_df['Age Category'].apply(
            lambda x: f"Young" if x == 'young' else f"Old"
        )
        
        # Plot with the updated labels
        sns.violinplot(x='Age Label', y=metric, data=temp_df, palette=palette, ax=ax, inner=None)
        sns.boxplot(x='Age Label', y=metric, data=temp_df, palette=palette, ax=ax, width=0.3, 
                    showfliers=False, linewidth=2)
        sns.stripplot(x='Age Label', y=metric, data=temp_df, color='black', jitter=True, 
                    ax=ax, size=5, alpha=opacity)
        
        # Statistical test
        old_data = temp_df[temp_df['Age Category'] == 'old'][metric]
        young_data = temp_df[temp_df['Age Category'] == 'young'][metric]
        
        if len(old_data) > 1 and len(young_data) > 1:  # Need at least 2 samples per group
            ttest = stats.ttest_ind(old_data, young_data, equal_var=False)  # Using Welch's t-test
            
            ymax = max(old_data.max(), young_data.max())
            
            if ttest.pvalue < 0.05:
                add_stat_annotation(ax, 0, 1, ymax, 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0]), ttest.pvalue)
        
        # increase the y-axis limit slightly for better visualization
        ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] + 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0]))

        # Customize y-axis ticks and labels
        ax.set_ylabel(title, fontsize=28, fontweight='bold', family='Arial')
        ax.tick_params(axis='y', labelsize=24)
        ax.tick_params(axis='x', labelsize=24)
        ax.set_xlabel("")
        ax.set_title("", fontsize=15)

        # bold the x-axis labels
        ax.set_xticklabels(ax.get_xticklabels(), fontweight='bold')

    plt.tight_layout()
    fig.subplots_adjust(wspace=0.3, hspace=0.3)

    # save the figure
    if not os.path.exists('results/hrver/physio_stats_hrver'):
        os.makedirs('results/hrver/physio_stats_hrver')
    fig.savefig('results/hrver/physio_stats_hrver/hrver_physio_violins.png', dpi=300)

# Modify the main analysis execution to use the updated functions
def run_updated_analysis():
    # Load the data as before
    df_hrv_pre, df_hrv_post, df_hrv_pre_only = load_data()
    
    # Apply the new outlier handling approach - mark outliers as NaN instead of removing rows
    physio_metrics = ['SDNN', 'RMSSD', 'LF', 'HF', 'Avg HR', 'Avg CO2', 'Breathing Rate']
    
    print("\nHandling outliers in Pre-only dataset:")
    df_hrv_pre_only_clean = handle_metric_outliers(df_hrv_pre_only, columns=physio_metrics)
    
    print("\nHandling outliers in Pre dataset:")
    df_hrv_pre_clean = handle_metric_outliers(df_hrv_pre, columns=physio_metrics)
    
    print("\nHandling outliers in Post dataset:")
    df_hrv_post_clean = handle_metric_outliers(df_hrv_post, columns=physio_metrics)
    
    # Report on remaining data after outlier handling
    print("\nData remaining after outlier handling:")
    for metric in physio_metrics:
        pre_count = df_hrv_pre_clean[metric].notna().sum()
        pre_total = len(df_hrv_pre_clean)
        pre_pct = pre_count / pre_total * 100
        
        post_count = df_hrv_post_clean[metric].notna().sum()
        post_total = len(df_hrv_post_clean)
        post_pct = post_count / post_total * 100
        
        pre_only_count = df_hrv_pre_only_clean[metric].notna().sum()
        pre_only_total = len(df_hrv_pre_only_clean)
        pre_only_pct = pre_only_count / pre_only_total * 100
        
        print(f"{metric}: pre={pre_count}/{pre_total} ({pre_pct:.1f}%), "
              f"post={post_count}/{post_total} ({post_pct:.1f}%), "
              f"pre_only={pre_only_count}/{pre_only_total} ({pre_only_pct:.1f}%)")
    
    # Create and display violin plots (using pre_only)
    create_old_young_violins(df_hrv_pre_only_clean)
    
    # Create summary tables with within and between group statistics
    summary_table, within_stats, between_stats = create_summary_table(
        df_hrv_pre_clean, df_hrv_post_clean, df_hrv_pre_only_clean, physio_metrics
    )
    
    # View the statistical results
    print("\nWithin-group statistics (Pre vs Post):")
    print(within_stats)
    print("\nBetween-group statistics (Young vs Old using pre_only data):")
    print(between_stats)
    
    return df_hrv_pre_clean, df_hrv_post_clean, df_hrv_pre_only_clean, summary_table, within_stats, between_stats

# Replace the main code execution
metrics = ['SDNN', 'RMSSD', 'LF', 'HF', 'Avg HR', 'Avg CO2', 'Breathing Rate']
df_hrv_pre_clean, df_hrv_post_clean, df_hrv_pre_only_clean, summary_table, within_stats, between_stats = run_updated_analysis()

# save summary table, within_stats, between_stats to csv files
if not os.path.exists('results/hrver/physio_stats_hrver'):
    os.makedirs('results/hrver/physio_stats_hrver')

summary_table.to_csv('results/hrver/physio_stats_hrver/hrver_summary_statistics_new.csv')
within_stats.to_csv('results/hrver/physio_stats_hrver/hrver_within_group_statistics_new.csv')
between_stats.to_csv('results/hrver/physio_stats_hrver/hrver_between_group_statistics_new.csv')
