# physio-bold-aging

Source code for "Physiological Component of the BOLD Signal: Influence of Age and Heart Rate Variability Biofeedback Training"

## Overview

This repository contains the analysis code for investigating the relationship between physiological signals and the BOLD fMRI response across different age groups, and how this relationship is modulated by heart rate variability (HRV) biofeedback training. The analysis examines:

- Age-related differences in physiological signals 
- Aged-related differences in percent variance explained of BOLD signal by physiological signals 
- Age-related differences in cross-correlation between physiological measures and BOLD signal
- Effects of HRV biofeedback training on cross-correlation between physiological measures and BOLD signal

## Repository Structure
```
.
├── analysis/
│   ├── pve_nki.py      # NKI dataset analysis
│   ├── pve_hrver.py    # HRV-ER dataset analysis
│   ├── cross_corr.py   # Cross-correlation analysis between physio and BOLD
│   ├── physio_stats.py # Statistical analysis of physiological measures
│   ├── randomise_nki_young_old.sh    # FSL randomise script for NKI analysis
│   └── randomise_hrver_young_old.sh  # FSL randomise script for HRV-ER analysis
├── data/
│   ├── nki_age_gender_v2.csv          # NKI participant information
│   ├── hrver_ses_pre_age_gender.csv   # HRV-ER pre-session data
│   ├── hrver_ses_post_age_gender.csv  # HRV-ER post-session data
│   └── MNI152_T1_2mm_brain.nii        # MNI template for masking
│   └── masks/         # Brain masks and templates
├── results/           # Analysis outputs (created during runtime)
│   ├── nki_pve_results_*/
│   └── hrver_pve_results_*/
└── logs/              # Analysis logs with timestamps (you can create this)
```

## Dependencies

Python version:
- Python 3.8.8

Additional requirements:
- FSL 6.0.1
- AFNI 21.1.03
- ANTs 2.3.0
- OpenMPI 3.1.4
- MATLAB 2020b

## Setup

1. Clone the repository:
```bash
git clone https://github.com/neurdylab/physio-bold-aging.git
cd physio-bold-aging
```

2. Create Python virtual environment and install dependencies:
```bash
mkdir -p logs
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

3. **Change File Paths**:
Before running the analysis, ensure that the file paths in the utility scripts match your directory structure.
Open the `utils/file_paths_hrver.py` file and verify the following path definitions:
   - The physiological data path is defined as:
       ```python
       physio_path = os.path.join(f'/data1/neurdylab/songrw/derivates/hrv_er/preproc_physio_ses-{session}', sub_id, f'{sub_id}_ses-{session}_task-rest_physio_physOUT.mat')
       ```
   - The BOLD data path is defined as:
       ```python
       bold_path = os.path.join(f'/data1/neurdylab/datasets/HRV-ER/HRV-ER_proc/{sub_id}/ses-{session}/func/ants_out', f'{sub_id}_ses-{session}_task-rest_bold_mo_EPI2MNI_sm_nr.nii.gz')
       ```
   Ensure that these paths reflect the actual locations of your physiological and BOLD data files.

   Before running the preprocessing, ensure that the following file paths in `preprocessing/preproc_imaging.sh` match your directory structure:
   - `maindir_raw` - The main directory containing the raw fMRI data:
     ```bash
     maindir_raw=/data1/neurdylab/datasets/HRV-ER/HRV-ER_raw            #Directory of raw data 
     maindir_proc=/data1/neurdylab/datasets/HRV-ER/HRV-ER_proc          #Directory to save preprocessed data 
     scripts_path="/data1/neurdylab/scripts/vu_meica_pipeline"          #Parent directory of AFNI cmake
     afni_init="singularity exec --bind /data1:/data1 ${scripts_path}/afni_cmake_build_AFNI_21.1.03.sif"     #Full path to AFNI cmake 
     ```

   You also need to ensure the following file paths in `preprocessing/preproc_physio.m` match your directory structure: 
     ```matlab 
      dir = '/path/to/output/';              % Output directory for processed data
      D = "/path/to/physio_data/";           % Source directory containing raw physiological data
      session = 'pre';                       % Session identifier (e.g., 'pre', 'post')
      data = readtable(append('/hrver_ses_', session, '_age_gender.csv'));  % Read the CSV file
      id = string(data.sub_id);            % Array of subject IDs
     ```

## Preprocessing Instructions

Before running the analysis, you'll need to preprocess both the physiological and fMRI data. Note that you'll need to modify the file paths in these scripts to match your directory structure (specifically for the locations of the physio and imaging data)

### Physiological Data
```bash
# Edit paths in preproc_physio.m first
matlab preprocessing/preproc_physio.m
```

### fMRI Data
```bash
# Edit paths in preproc_imaging.sh first
bash preprocessing/preproc_imaging.sh
```

## Running the Analysis

### NKI Dataset Analysis

The `pve_nki.py` script analyzes physiological variance in the NKI dataset.

Available covariates:
- gender
- lf (low frequency HRV)
- hf (high frequency HRV)
- avg_hr (average heart rate)
- sd_rv (standard deviation of respiratory variation)
- br (breathing rate)

To run the analysis and subsequent FSL randomise:
```bash
# Run PVE analysis
nohup python -u analysis/nki/pve_nki.py > logs/pve_nki_$(date +%Y%m%d_%H%M%S).log 2>&1 &

# Run FSL randomise on the results
nohup bash analysis/nki/randomise_nki_young_old.sh > logs/randomise_nki_$(date +%Y%m%d_%H%M%S).log 2>&1 &
```

### HRV-ER Dataset Analysis

#### 1. PVE Analysis (Young vs Old)
The `pve_hrver.py` script analyzes the PVE of BOLD explained by HR and CO2 in the HRV-ER dataset.

Available covariates:
- gender
- lf (low frequency HRV)
- hf (high frequency HRV)
- avg_hr (average heart rate)
- avg_co2 (average end-tidal CO2)
- br (breathing rate)

To run the analysis and subsequent FSL randomise:
```bash
# Run PVE analysis
nohup python -u analysis/hrver/pve_hrver.py > logs/pve_hrver_$(date +%Y%m%d_%H%M%S).log 2>&1 &

# Run FSL randomise on the results
nohup bash analysis/hrver/randomise_hrver_young_old.sh > logs/randomise_hrver_$(date +%Y%m%d_%H%M%S).log 2>&1 &
```
The analysis creates directories named according to the covariates used. For example:
- `results/nki/pve_results_gender/` - results controlling for gender only
- `results/hrver/pve_results_gender_avg_hr/` - results controlling for gender and average heart rate

Each directory contains:
- `*_cov_young_old.nii.gz` - Variance explained maps
- `design_matrix.txt` - Design matrix for FSL randomise
- `contrast_matrix.txt` - Contrast matrix for FSL randomise

#### 2. Cross-Correlation (tissue specific)

To run the cross-correlation analysis for young vs. old participants:

- Open the `analysis/hrver/cross_corr_avg_young_old_hrver.py` file.
- Modify the `tissue` variable to specify the tissue type you want to analyze (e.g., "gray_matter", "white_matter", or "ventricle").
- Run the script:
  ```bash
  nohup python -u analysis/hrver/cross_corr_avg_young_old_hrver.py > logs/cross_corr_young_old_$(date +%Y%m%d_%H%M%S).log 2>&1 &
  ```

To run the cross-correlation analysis for pre vs. post conditions:

- Open the `analysis/hrver/cross_corr_avg_pre_post_hrver.py` file.
- Modify the `age_group` variable to specify the age group you want to analyze (e.g., "young" or "older").
- Modify the `osc_condition` variable to specify the oscillation condition (e.g., "osc+" or "osc-").
- Modify the `tissue` variable to specify the tissue type you want to analyze (e.g., "gray_matter", "white_matter", or "ventricle").
- Run the script:
  ```bash
  nohup python -u analysis/hrver/cross_corr_avg_pre_post_hrver.py > logs/cross_corr_pre_post_$(date +%Y%m%d_%H%M%S).log 2>&1 &
  ```

The analysis creates directories named according to the type of cross correlation analysis performed. For example: 
- `results/hrver/cross_corr_avg_pre_post_older_osc+` - results comparing cross correlation before and after osc+ in older adults 
- `results/hrver/cross_corr_avg_young_old` - results comparing cross correlation between older and younger adults 

Each directory contains contains the cross correlation analysis performed on the specified tissue type. For example: 
- `gray_matter.png` - HR-BOLD and CO2-BOLD cross correlation analysis done on BOLD signal averagd over gray matter, which significant regions shaded.


#### 3. Cross-Correlation (whole brain): Young vs. Old Analysis

To run the cross-correlation analysis for young vs. old participants:

- Open the `analysis/hrver/cross_corr_whole_brain_young_old_hrver.py` file.
- To expedite the analysis, you can modify the lags that the analysis is performed on as well as the number of tfce permutations: 
```python
# Lines 259 - 265
run_randomise(
      temp_4d,
      os.path.join(output_dir, "temp", f"{measure}_lag_{lag_idx}"),
      os.path.join(output_dir, "design.mat"),
      os.path.join(output_dir, "design.con"),
      n_permutations=1000 # 5000 used in the paper, but 1000 gives similar results 
)
```
```python
# Line 193
def main(start_lag=-1, end_lag=9): # These are the lags (in TRs) used in the paper 
```
- Run the script:
  ```bash
  nohup python -u analysis/hrver/cross_corr_whole_brain_young_old_hrver.py > logs/cross_corr_whole_brain_young_old_$(date +%Y%m%d_%H%M%S).log 2>&1 &
  ```

The whole brain cross-correlation analysis results are stored in the `results/hrver/cross_corr_whole_brain_young_old` directory. The structure of this directory is as follows:

```
└── cross_corr_whole_brain_young_old/
    ├── temp/                                   - Temporary files generated during the analysis, including intermediate results.
    ├── hr_young-gt-old_tstat_4d.nii.gz        - Statistical map showing the t-statistics for HR correlations between young and old participants.
    ├── co2_young-gt-old_tstat_4d.nii.gz       - Statistical map showing the t-statistics for CO2 correlations between young and old participants.
    ├── design.mat                              - Design matrix used for the FSL randomise analysis.
    ├── design.con                              - Contrast matrix used for the FSL randomise analysis.
    ├── hr_young-gt-old_corrp_4d.nii.gz        - TFCE-corrected p-values for HR correlations (young > old).
    ├── co2_young-gt-old_corrp_4d.nii.gz       - TFCE-corrected p-values for CO2 correlations (young > old).
    ├── hr_old-gt-young_corrp_4d.nii.gz        - TFCE-corrected p-values for HR correlations (old > young).
    └── co2_old-gt-young_corrp_4d.nii.gz       - TFCE-corrected p-values for CO2 correlations (old > young).
```
Each statistical map is saved in NIfTI format and can be visualized using neuroimaging software such as AFNI.

#### 4. Physio stats 
To obtain stats for the physiological metrics (rmssd, lf hrv, hf hrv, breathing rate, average heart rate, average end tidal co2), run: 
 ```bash
  nohup python -u analysis/hrver/physio_stats_hrver.py > logs/physio_stats_hrver_$(date +%Y%m%d_%H%M%S).log 2>&1 &
 ```
The outputs are stored under the directory `results/hrver/physio_stats_hrver`, which will contain: 
- `hrv_summary_statistics.csv` - mean and variances of the physiological metrics for each group (old/young, osc +/-, pre/post)
- `hrv_between_group_statistics.csv` - old vs young t-stat and p-values for each physiological metric 
- `hrv_within_group_statistics.csv` - pre vs post (osc +/- and old/young) t-stat and p-values for each physiological metric 

## Data Availability 
The HRV-ER dataset can be found online on [OpenNeuro][openneuro-link].
Preprocessed Physiological Data for HRV-ER dataset is available on [Box][box-link].

[openneuro-link]: https://openneuro.org/datasets/ds003823/versions/1.2.0
[box-link]: https://vanderbilt.app.box.com/s/2v0qwfitb07crqjgtorlg5wtgs2y3mhk

## Contact 
If you have any questions, please email me at richard.w.song@vanderbilt.edu
