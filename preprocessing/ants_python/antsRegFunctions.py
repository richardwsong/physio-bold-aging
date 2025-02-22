##% module containing functions for ants alignment
# note, paths are currently set for running on amri raid
# catie
import os

def ants_registration(in_nii, ref_nii, out_nii, out_tform, tip, restrict=None, exe=True):
    '''ants registration'''

    lin_tform_params = ' '.join([
        '--metric MI[' + ref_nii + ',' + in_nii + ',1,32,Regular,0.25]',
        '--convergence [1000x500x250x125]',
        '--shrink-factors 12x8x4x2',
        '--smoothing-sigmas 4x3x2x1vox',
    ])
    syn_tform_params = ' '.join([
        '--metric CC[' + ref_nii + ',' + in_nii + ',1,4]',
        '--convergence [100x100x70x50x20]',
        '--shrink-factors 10x6x4x2x1',
        '--smoothing-sigmas 5x3x2x1x0vox',
    ])
    antsRegistration_call = ' '.join([
        'antsRegistration',
        '--initial-moving-transform [' + ref_nii + ',' + in_nii + ',1]',
        '--output [' + out_tform + ',' + out_nii + ']',
        '--dimensionality 3',
        '--float 1',
        '--interpolation Linear',
        '--winsorize-image-intensities [0.005,0.995]',
        '--use-histogram-matching 0',
        ('--restrict-deformation ' + restrict if restrict else ''),
        ('--transform Translation[0.1] '      + lin_tform_params if 't' in tip else ''),
        ('--transform Rigid[0.1] '            + lin_tform_params if 'r' in tip else ''),
        ('--transform Similarity[0.1] '       + lin_tform_params if 'i' in tip else ''),
        ('--transform Affine[0.1] '           + lin_tform_params if 'a' in tip else ''),
        ('--transform SyN[0.1,3,0] '          + syn_tform_params if 's' in tip else ''),
        ('--transform BSplineSyN[0.1,26,0,3]' + syn_tform_params if 'b' in tip else '')
    ])
    if exe:
        os.system(antsRegistration_call)
    else:
        return(antsRegistration_call)


def ants_transformation(in_nii, ref_nii, out_nii, in_tform, interpolation='Linear'):
    '''application of ants transform'''
    print("hi")
    antsTransformation_call = ' '.join([
        'antsApplyTransforms',
        '--dimensionality 3',
        '--input', in_nii,
        '--reference-image', ref_nii,
        '--output', out_nii,
        '--interpolation', interpolation,
        '--transform', in_tform,
        '--float',
        '--verbose', "1",
        '-e',"3",
    ])
    os.system(antsTransformation_call)
