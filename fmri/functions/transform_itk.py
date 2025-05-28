import nibabel as nib
import numpy as np
from nipype.interfaces.c3 import C3dAffineTool  
from pathlib import Path
import subprocess

# this code is based on a thread in the neurostars forum: https://neurostars.org/t/what-do-itk-snap-registrations-compute-exactly/24292/14
# it requires you to produce a r.mat file first with itksnap.


basePath = Path('/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/hipp/bids/derivatives/fmriprep/')

sub = 'sub-68227'
subjPath = basePath / sub
# argsSession2 = f"c3d_affine_tool -itk {subjPath}ses-2/anat/sub-68227_ses-2_acq-t1csmp2ragesag06mmUNIDEN_from-orig_to-T1w_mode-image_xfm.txt -o ses-2_to_t1w_PYthon.mat"
# c3 = C3dAffineTool()
# c3.inputs.reference_file = basePath / subjPath / 'ses-1/anat/sub-68227_ses-1_acq-t1csmp2ragesag06mmUNIDEN_from-orig_to-T1w_mode-image_xfm.txt'
# c3.inputs.fsl2ras = False
# c3.cmdline = argsSession2
# c3.run()
matfile = subjPath / f'ses-1/anat/{sub}_ses-1_to_t1w.mat'
itkfile = subjPath / f"ses-1/anat/{sub}_ses-1_acq-t1csmp2ragesag06mmUNIDEN_from-orig_to-T1w_mode-image_xfm.txt"
# R2 = np.loadtxt("ses-2_to_t1w_NoInfo.mat", delimiter=" ")
R = np.loadtxt(matfile, delimiter=" ")
# Rfunc = np.loadtxt("ses-2_to_t1wfunc.mat", delimiter=" ")
# img2 = nib.load("sub-68227_ses-2_acq-t1csmp2ragesag06mmUNIDEN_T1w.nii.gz")
img2 = nib.load("meanasub-68227_ses-1_task-retrieval_acq-epfid2d1200_bold.nii")
# Rinversed = np.loadtxt("ses-2_to_t1w_inv.mat", delimiter=" ")


img2_registered_on_img1 = nib.nifti1.Nifti1Image(
    img2.get_fdata(),
    np.linalg.inv(R) @ img2.affine
)

img2_registered_on_img1.to_filename(
    "ses1mean-registered-on_t1w.nii.gz"
)
