#!/usr/bin/python
'''
 calls the functions in antsRegFunctions.py to align a T1 to MNI image
 but could be more general! just align in_nii to ref_nii.
 author: catie
    version: 25 May 2017

11 March 2020: modified paths for merlion
'''

import argparse
import sys, os

sys.path.insert(0, '/data1/neurdylab/eegfmri_vu_pipeline/scripts/ants_python/')
from antsRegFunctions import ants_registration


parser = argparse.ArgumentParser(description='Spatially align in_nii to ref_nii')
parser.add_argument('--in_nii', dest='in_nii',
                   default='', 
                   help='the input T1 image')
parser.add_argument('--ref_nii', dest='ref_nii',
                   default='',
                   help='reference (usually MNI) image to which we are aligning the T1')
parser.add_argument('--out_dir', dest='out_dir',
                   default='', 
                   help='output dir for writing xforms')
args = parser.parse_args()

my_in_nii = args.in_nii
my_ref_nii = args.ref_nii
out_dir = args.out_dir

# default output directory will be same as input filename
if out_dir=='':
    out_dir = os.path.dirname(my_in_nii)
    
#%%
# for testing: T1 to MNI
#my_in_nii = '/Users/catie/brain/sleepy/test/test_align/anat_bet1.nii'
#my_ref_nii = '/Users/catie/brain/sleepy/test/test_align/MNI152_T1_1mm_brain.nii'
#out_dir = '/Users/catie/brain/sleepy/test/test_align/';
# for testing: EPI to T1
#my_in_nii = '/Users/catie/brain/sleepy/test/test_align_func/oneVol_homcor_restore.nii'
#my_ref_nii = '/Users/catie/brain/sleepy/test/test_align/anat_bet1.nii'
#out_dir = '/Users/catie/brain/sleepy/test/test_align_func/'


#%%
# a series of rigid, affine, syn transforms here
ants_registration(
    in_nii = my_in_nii,
    ref_nii = my_ref_nii,
    out_nii =   out_dir +  'reg_tmp_r.nii',
    out_tform = out_dir +  'tform1_',
    tip = 'r'
)
ants_registration(
    in_nii = out_dir + 'reg_tmp_r.nii',
    ref_nii = my_ref_nii,
    out_nii =   out_dir +     'reg_tmp_a.nii',
    out_tform = out_dir +     'tform2_',
    tip = 'a'
)
ants_registration(
    in_nii = out_dir + 'reg_tmp_a.nii',
    ref_nii = my_ref_nii,
    out_nii = out_dir +     'reg_tmp_s.nii',
    out_tform = out_dir +   'tform3_',
    tip = 's'
)
