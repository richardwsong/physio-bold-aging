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

## Data Preprocessing Pipeline

### Physiological Data Processing
- Processes respiration, CO2, and cardiac measures from fMRI sessions
- Generates physiological regressors for fMRI analysis
- Performs outlier detection and cleaning of heart rate data
- Creates visualization plots for QA

### fMRI Preprocessing
1. Motion correction using AFNI's 3dvolreg
2. Slice timing correction
3. Multi-echo ICA denoising using tedana
4. Registration to anatomical T1 and MNI space
5. Spatial smoothing and nuisance regression

## Analysis Pipeline

1. **Cross-correlation Analysis**
   - Whole-brain cross-correlation maps for physiological signals
   - Group comparisons between young and old participants
   - Statistical testing using FSL's randomise

2. **Variance Explained Analysis**
   - Calculates percent variance in BOLD explained by:
     - Heart rate
     - Respiratory variation
     - End-tidal CO2
     - Combined physiological measures
   - Age group comparisons
   - Pre vs post HRV Biofeedback training comparisons

3. **Physiological Statistics**
   - HRV metrics calculation (SDNN, RMSSD, LF/HF)
   - Breathing rate analysis
   - Group comparisons and statistical testing

## Dependencies

Python version:
- Python 3.8.8

Python packages (specified versions in requirements.txt):
- matplotlib
- nibabel
- numpy
- pandas
- scipy
- seaborn

Additional requirements:
- FSL 6.0.1
- AFNI
- ANTs
- MATLAB (for physiological preprocessing)

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

## Preprocessing Instructions

Before running the analysis, you'll need to preprocess both the physiological and fMRI data. Note that you'll need to modify the file paths in these scripts to match your directory structure (specifically for the locations of the physio and imaging data)

### Physiological Data
```bash
matlab preprocessing/preproc_physio.m
```

### fMRI Data
```bash
# Edit paths in preproc_imaging.sh first
# Then run the preprocessing pipeline
bash preprocessing/preproc_imaging.sh

# The script performs:
# 1. Motion correction
# 2. Slice timing correction
# 3. Registration to MNI space
# 4. Spatial smoothing
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

The cross-correlation analysis can be performed to compare physiological signals (HR and CO2) with BOLD signals across different conditions (young vs. old and pre vs. post).

#### 2. Cross-Correlation (tissue specific): Young vs. Old Analysis

To run the cross-correlation analysis for young vs. old participants:

- Open the `analysis/hrver/cross_corr_avg_young_old_hrver.py` file.
- Modify the `tissue` variable to specify the tissue type you want to analyze (e.g., "gray_matter", "white_matter", or "ventricle").
- Run the script:
  ```bash
  nohup python -u analysis/hrver/cross_corr_avg_young_old_hrver.py > logs/cross_corr_young_old_$(date +%Y%m%d_%H%M%S).log 2>&1 &
  ```

#### 3. Cross-Correlation (tissue specific): Pre vs. Post Analysis

To run the cross-correlation analysis for pre vs. post conditions:

- Open the `analysis/hrver/cross_corr_avg_pre_post_hrver.py` file.
- Modify the `age_group` variable to specify the age group you want to analyze (e.g., "young" or "older").
- Modify the `osc_condition` variable to specify the oscillation condition (e.g., "osc+" or "osc-").
- Modify the `tissue` variable to specify the tissue type you want to analyze (e.g., "gray_matter", "white_matter", or "ventricle").
- Run the script:
  ```bash
  nohup python -u analysis/hrver/cross_corr_avg_pre_post_hrver.py > logs/cross_corr_pre_post_$(date +%Y%m%d_%H%M%S).log 2>&1 &
  ```

#### 4. Cross-Correlation (whole brain): Young vs. Old Analysis

To run the cross-correlation analysis for young vs. old participants:

- Open the `analysis/hrver/cross_corr_whole_brain_young_old_hrver.py` file.
- To expedite the analysis, you can modify the lags that the analysis is performed on as well as the number of tfce permutations. 
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

## Output

### PVE Analysis
The analysis creates directories named according to the covariates used. For example:
- `results/nki/pve_results_gender/` - results controlling for gender only
- `results/hrver/pve_results_gender_avg_hr/` - results controlling for gender and average heart rate

Each directory contains:
- `*_cov_young_old.nii.gz` - Variance explained maps
- `design_matrix.txt` - Design matrix for FSL randomise
- `contrast_matrix.txt` - Contrast matrix for FSL randomise

### Cross Correlation Analysis (tissue specific)
The analysis creates directories named according to the type of cross correlation analysis performed. For example: 
- `results/hrver/cross_corr_avg_pre_post_older_osc+` - results comparing cross correlation before and after osc+ in older adults 
- `results/hrver/cross_corr_avg_young_old` - results comparing cross correlation between older and younger adults 

Each directory contains contains the cross correlation analysis performed on the specified tissue type. For example: 
- `gray_matter.png` - HR-BOLD and CO2-BOLD cross correlation analysis done on BOLD signal averagd over gray matter, which significant regions shaded.

## Logs

All script outputs are stored in the `logs/` directory with timestamps:
- `pve_nki_YYYYMMDD_HHMMSS.log` - Output from NKI analysis
- `pve_hrver_YYYYMMDD_HHMMSS.log` - Output from HRV-ER analysis
- `randomise_*_YYYYMMDD_HHMMSS.log` - Output from FSL randomise

## Data Availability 
The HRV-ER dataset can be found online on [OpenNeuro][openneuro-link].
Preprocessed Physiological Data for HRV-ER dataset is available on [Box][box-link].

## Notes

- The analysis uses parallel processing with 16 cores by default
- Age groups are defined as:
  - Young: < 50 years
  - Old: ≥ 50 years
- Design matrices are formatted for FSL compatibility
- All physiological measures are detrended before analysis

[openneuro-link]: https://openneuro.org/datasets/ds003823/versions/1.2.0
[box-link]: https://vanderbilt.app.box.com/s/2v0qwfitb07crqjgtorlg5wtgs2y3mhk

