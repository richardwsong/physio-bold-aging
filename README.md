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
├── preprocessing/
│   ├── preproc_physio.m    # Physiological data preprocessing
│   └── preproc_imaging.sh  # fMRI preprocessing pipeline
├── analysis/
│   ├── cross_corr.py       # Cross-correlation analysis
│   ├── percent_var_explained.ipynb  # Percent Variance Explained analysis
│   └── physio_stats.py     # Statistical analysis of physiological measures
├── metadata/               # Subject information and experiment metadata
└── requirements.txt        # Python dependencies
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

## Installation and Setup

1. Clone the repository:
   ```bash
   git clone [repository-url]
   ```

2. Install Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Ensure FSL, AFNI, and ANTs are installed and properly configured in your environment.

## Usage

1. **Physiological Data Preprocessing**
   ```matlab
   % In MATLAB
   run preprocessing/preproc_physio.m
   ```

2. **fMRI Preprocessing**
   ```bash
   bash preprocessing/preproc_imaging.sh
   ```

3. **Analysis**
   - Run cross-correlation analysis:
     ```python
     python analysis/cross_corr.py
     ```
   - Run variance explained analysis using Jupyter notebook:
     ```bash
     jupyter notebook analysis/percent_var_explained.ipynb
     ```
   - Run physiological statistics:
     ```python
     python analysis/physio_stats.py
     ```

## Data Availability

The preprocessed physiological data used in this analysis is available at:
https://vanderbilt.app.box.com/s/2v0qwfitb07crqjgtorlg5wtgs2y3mhk

## Contact

Email: richard.w.song@vanderbilt.edu 
