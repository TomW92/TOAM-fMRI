# In[]
import nitransforms as nt


# In[]
itk_file = nt.linear.load(
    '/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/hipp/bids/derivatives/fmriprep/sub-60601/ses-1/anat/sub-60601_ses-1_acq-t1csmp2ragesag06mmUNIDEN_from-orig_to-T1w_mode-image_xfm.txt',
    fmt='itk',
    moving='/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/hipp/bids/sub-60601/ses-1/anat/sub-60601_ses-1_acq-t1csmp2ragesag06mmUNIDEN_T1w.nii.gz',
    reference='/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/hipp/bids/derivatives/fmriprep/sub-60601/anat/sub-60601_acq-t1csmp2ragesag06mmUNIDEN_desc-preproc_T1w.nii.gz'
)
itk_file.to_filename('transform.lta',fmt='fs')


# In[] %%
#compare affines!
import nibabel as nb
files={
    'ses-1anat':'/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/hipp/bids/derivatives/spmMB/sub-60601/preproc/ses-1/anat/sub-60601_ses-1_acq-t1csmp2ragesag06mmUNIDEN_T1w.nii',
    'also-ses-1anat':'/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/hipp/bids/sub-60601/ses-1/anat/sub-60601_ses-1_acq-t1csmp2ragesag06mmUNIDEN_T1w.nii.gz',
    'robust_template_t1w':'/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/hipp/bids/derivatives/fmriprep/sub-60601/anat/sub-60601_acq-t1csmp2ragesag06mmUNIDEN_desc-preproc_T1w.nii.gz',
    'prob_wrongCA1Mask':'/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/hipp/bids/derivatives/fmriprep/sub-60601/ses-1/CA3_mask_prob_wrong_space.nii',
    'prob_rightCA1Mask':'/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/hipp/bids/derivatives/fmriprep/sub-60601/ses-1/CA3_mask_ses-1_space.nii',
    'initial_maskCA1Mask':'/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/hipp/bids/derivatives/freesurfer/sub-60601/mri/CA3_mask.nii',
    'rawavg':'/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/hipp/bids/derivatives/freesurfer/sub-60601/mri/rawavg.mgz' 
}

affines = []
names = []
for key,value in files.items():
    nii = nb.load(value)
    hdr = nii.header.copy()
    aff = nii.affine.copy()
    names.append(key)
    affines.append(aff)
    print(key + 'has the following affine:')
    print(aff)
# %%


import numpy
from nilearn.image import resample_img, load_img
def readITKtransform(transform_file):
    '''
    '''

    # read the transform
    transform = None
    with open(transform_file, 'r') as f:
        for line in f:

            # check for Parameters:
            if line.startswith('Parameters:'):
                values = line.split(': ')[1].split(' ')

                # filter empty spaces and line breaks
                values = [float(e) for e in values if (e != '' and e != '\n')]
                # create the upper left of the matrix
                transform_upper_left = numpy.reshape(values[0:9], (3, 3))
                # grab the translation as well
                translation = values[9:]

            # check for FixedParameters:
            if line.startswith('FixedParameters:'):
                values = line.split(': ')[1].split(' ')

                # filter empty spaces and line breaks
                values = [float(e) for e in values if (e != '' and e != '\n')]
                # setup the center
                center = values

    # compute the offset
    offset = numpy.ones(4)
    for i in range(0, 3):
        offset[i] = translation[i] + center[i];
        for j in range(0, 3):
            offset[i] -= transform_upper_left[i][j] * center[i]

    # add the [0, 0, 0] line
    transform = numpy.vstack((transform_upper_left, [0, 0, 0]))
    # and the [offset, 1] column
    transform = numpy.hstack((transform, numpy.reshape(offset, (4, 1))))

    return transform



config=setup_configNL(60601)
config['ses-2_Affine']  = '/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/hipp/bids/derivatives/fmriprep/sub-60601/ses-2/anat/sub-60601_ses-2_acq-t1csmp2ragesag06mmUNIDEN_from-orig_to-T1w_mode-image_xfm.txt'

transform = readITKtransform(config['ses-2_Affine'])


print(transform)

ses_2_native_t1     = load_img("/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/hipp/bids/derivatives/spmMB/sub-60601/preproc/ses-2/anat/sub-60601_ses-2_acq-t1csmp2ragesag06mmUNIDEN_T1w.nii")
ses_2_NATIVE_meanimg = load_img("/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/hipp/bids/derivatives/spm/sub-60601/preproc/func/ses-2/meansub-60601_ses-2_task-forgetexp_acq-epfid2d1200_bold.nii")
ses_2_t1w_meanimg = resample_img(ses_2_native_t1,transform)
ses_2_t1w_meanimg.to_filename("/Users/tomwillems/Desktop/testtransformANAT.nii")


## possibly, the transform is named wrong , see https://neurostars.org/t/fmriprep-applying-transforms-from-orig-to-t1w-mode-image-query/6729
## so let's try the other wayround.
tw1space_mask = config['fusiform_mask']
ses2_native_mask = resample_img(tw1space_mask,transform,interpolation="nearest")
ses2_native_mask.to_filename("/Users/tomwillems/Desktop/otherwayroundtransformANAT.nii")
