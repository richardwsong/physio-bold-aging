#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 17:50:17 2017

Apply concat rigid->affine->syn transforms from pre-computed Ants reg.
This is usually a good way to do nonlinear warping from T1->standard MNI space

@author: catie

13 June 2017: modifier to allow 'dimens' input, in case of 4d
nifti. Also changed inputs to antsRegFunctions.py.
11 March 2020: modified paths for merlion
"""

import argparse
import sys, os

sys.path.insert(0, '/data1/neurdylab/eegfmri_vu_pipeline/scripts/ants_python/')
from antsRegFunctions import ants_transformation

parser = argparse.ArgumentParser(description='Apply Ants Alignment')
parser.add_argument('--in_nii', dest='in_nii',
                   default='', 
                   help='input image to be registered')
parser.add_argument('--out_nii', dest='out_nii',
                   default='', 
                   help='registered output')
parser.add_argument('--ref_nii', dest='ref_nii',
                   default='',
                   help='reference (usually brain-extracted MNI), perhaps defines the space?')
parser.add_argument('--tform_dir', dest='tform_dir',
                   default='', 
                   help='output dir for writing xforms')
args = parser.parse_args()

my_in_nii = args.in_nii
my_ref_nii = args.ref_nii
my_out_nii = args.out_nii
tform_dir = args.tform_dir

if my_out_nii=='':
    prefix,ext = os.path.splitext(my_in_nii)
    my_out_nii = prefix + '_warpAnts.nii'

# build the transform and call the function
ants_transformation(
    in_nii = my_in_nii,
    ref_nii = my_ref_nii,
    out_nii = my_out_nii,
    in_tform = ' '.join([tform_dir + 'tform3_1Warp.nii.gz',
                        tform_dir  + 'tform3_0GenericAffine.mat',
                        tform_dir  + 'tform2_0GenericAffine.mat',
                        tform_dir  + 'tform1_0GenericAffine.mat']),
    )
