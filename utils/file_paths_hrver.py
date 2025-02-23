import os

def get_file_paths(session, sub_id):
    """
    Generate file paths for physiological and BOLD data for the HRV-ER dataset.

    Parameters:
        session (str): The session identifier ('pre' or 'post').
        sub_id (str): The subject ID.

    Returns:
        dict: A dictionary containing the file paths for physiological and BOLD data.
    """
    base_path = '/data1/neurdylab/songrw/derivates/hrv_er'
    
    # Define paths for physiological data
    physio_path = os.path.join(base_path, f'preproc_physio_ses-{session}', sub_id, f'{sub_id}_ses-{session}_task-rest_physio_physOUT.mat')
    
    # Define paths for BOLD data
    bold_path = os.path.join('/data1/neurdylab/datasets/HRV-ER/HRV-ER_proc', sub_id, f'ses-{session}', 'func/ants_out', f'{sub_id}_ses-{session}_task-rest_bold_mo_EPI2MNI_sm_nr.nii.gz')
    
    return {
        'physio': physio_path,
        'bold': bold_path
    } 