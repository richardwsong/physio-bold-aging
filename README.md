# physio-bold-aging

Source code for "Physiological Component of the BOLD Signal: Influence of Age and Heart Rate Variability Biofeedback Training"

## Overview

This repository contains the analysis code for investigating the relationship between physiological signals and the BOLD fMRI response across different age groups, and how this relationship is modulated by heart rate variability (HRV) biofeedback training. The analysis examines:

- Age-related differences in physiological contributions to BOLD signal
- Effects of HRV biofeedback training on BOLD signal physiology
- Cross-correlation between physiological measures and BOLD signal
- Variance explained by different physiological components

## Repository Structure
```
.
├── analysis/
│   ├── pve_nki.py      # NKI dataset analysis
│   ├── pve_hrver.py    # HRV-ER dataset analysis
│   ├── cross_corr.py   # Cross-correlation analysis between physio and BOLD
│   └── physio_stats.py # Statistical analysis of physiological measures
├── metadata/
│   ├── nki_age_gender_v2.csv          # NKI participant information
│   ├── hrver_ses_pre_age_gender.csv   # HRV-ER pre-session data
│   ├── hrver_ses_post_age_gender.csv  # HRV-ER post-session data
│   └── MNI152_T1_2mm_brain.nii        # MNI template for masking
├── data/               # Analysis outputs (created during runtime)
│   ├── nki_pve_results_*/
│   └── hrver_pve_results_*/
└── logs/              # Analysis logs with timestamps
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
nohup python -u analysis/pve_nki.py > logs/pve_nki_$(date +%Y%m%d_%H%M%S).log 2>&1 &

# Run FSL randomise on the results
nohup bash analysis/randomise_nki_young_old.sh > logs/randomise_nki_$(date +%Y%m%d_%H%M%S).log 2>&1 &
```

### HRV-ER Dataset Analysis

The `pve_hrver.py` script analyzes the HRV-ER dataset.

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
nohup python -u analysis/pve_hrver.py > logs/pve_hrver_$(date +%Y%m%d_%H%M%S).log 2>&1 &

# Run FSL randomise on the results
nohup bash analysis/randomise_hrver_young_old.sh > logs/randomise_hrver_$(date +%Y%m%d_%H%M%S).log 2>&1 &
```

## Output

The analysis creates directories named according to the covariates used. For example:
- `data/nki_pve_results_gender/` - results controlling for gender only
- `data/hrver_pve_results_gender_avg_hr/` - results controlling for gender and average heart rate

Each directory contains:
- `*_cov_young_old.nii.gz` - Variance explained maps
- `design_matrix.txt` - Design matrix for FSL randomise
- `contrast_matrix.txt` - Contrast matrix for FSL randomise

## Logs

All script outputs are stored in the `logs/` directory with timestamps:
- `pve_nki_YYYYMMDD_HHMMSS.log` - Output from NKI analysis
- `pve_hrver_YYYYMMDD_HHMMSS.log` - Output from HRV-ER analysis
- `randomise_*_YYYYMMDD_HHMMSS.log` - Output from FSL randomise

## Data Availability 
The HRV-ER dataset can be found online on OpenNeuro: https://openneuro.org/datasets/ds003823/versions/1.2.0 
Preprocessed Physiological Data for HRV-ER dataset: https://vanderbilt.app.box.com/s/2v0qwfitb07crqjgtorlg5wtgs2y3mhk 

## Notes

- The analysis uses parallel processing with 16 cores by default
- Age groups are defined as:
  - Young: < 50 years
  - Old: ≥ 50 years
- Design matrices are formatted for FSL compatibility
- All physiological measures are detrended before analysis
