import os
import numpy as np
import nibabel as nib
from scipy.io import loadmat, savemat
from scipy.linalg import pinv
from scipy.signal import correlate, detrend, convolve
from scipy.interpolate import interp1d
import concurrent.futures
import pandas as pd
import matplotlib.pyplot as plt
import warnings 
warnings.filterwarnings("ignore")

dt = 2.4 # HRV-ER dataset has TR of 2.4 s, NKI uses 1.4 s TR 
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