# %%
import os, sys
from pathlib import Path

home = Path(os.environ["HOME"])
sys.path.append(str(home))
from TOAM.fmri.mvpa.mvpa_utils.get_datasets import *
from nilearn.glm.contrasts import compute_contrast
from nilearn import plotting, datasets, surface, image, datasets, masking
from nilearn.plotting import plot_stat_map, plot_img_on_surf, plot_anat
from nilearn.glm import threshold_stats_img
from pathlib import Path
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
from matplotlib.transforms import Affine2D
from scipy.stats import norm, gaussian_kde, pearsonr, zscore

seq = "hipp"
basePath = Path("/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/")
dataPath = basePath / "data/fMRI/"
spmMB = Path(
    "/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/"
    + seq
    + "/bids/derivatives/spmMB"
)
Workspace = Path(
    "/storage/workspaces/psy_memory_wfg_psy/hpc_henke_wfg/s2019_twillems_fMRI_silent_engram/data/fMRI/"
)
fmriPrep = Path(
    "/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/"
    + seq
    + "/bids/derivatives/fmriprep/"
)
deriv = Path(
    "/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/"
    + seq
    + "/bids/derivatives/"
)
RFX = spmMB / "RFX_TContrasts_fl_explmask"
corRet3 = spmMB / "RFX_old_Correlations" / "RFX_ret3unconsc_RFXCorr_s6wra" / "ses-1"
bg_image = dataPath / "utilities" / "mniStudyTemplate.nii"
fname = (
    corRet3
    / "CORREL_faceRecog_correct_unconsc 1+2_faceRecog_incorrect_unconsc 1+2_X_accuracy"
    / "con_0001.nii"
)
con1 = (
    RFX
    / "faceRecog_correct_unconsc 1+2_faceRecog_incorrect_unconsc 1+2"
    / "con_0001.nii"
)
con2 = (
    corRet3
    / "CORREL_faceRecog_consc 1+2_faceRecog_unconsc 1+2_X_accuracy"
    / "con_0001.nii"
)
con3 = corRet3 / "CORREL_faceRecog_correct_unconsc 2_number_X_accuracy" / "con_0001.nii"
con4 = (
    corRet3
    / "CORREL_faceRecog_correct_unconsc 1+2_faceRecog_incorrect_unconsc 1+2_X_accuracy"
    / "con_0003.nii"
)  # actual correlation
home = Path("/storage/homefs/tw18a205")
correlation_cluster = (
    spmMB
    / "2nd_level"
    / "regressions"
    / "con_images"
    / "faceRecog_correct_unconsc_1_vs_faceRecog_incorrect_unconsc_1_x_unconsciousRetrievalAccuracy_Day1"
    / "spmT_0001.nii"
)
correlation_cluster_day2 = (
    spmMB
    / "2nd_level"
    / "regressions"
    / "con_images"
    / "faceRecog_correct_unconsc_2_vs_faceRecog_incorrect_unconsc_2_x_unconsciousRetrievalAccuracy_Day2"
    / "spmT_0001.nii"
)
behav_correlationdt = (
    spmMB
    / "2nd_level"
    / "regressions"
    / "con_images"
    / "faceRecog_correct_unconsc_1_vs_faceRecog_incorrect_unconsc_1_x_unconsciousRetrievalAccuracy_Day1"
    / "corfile.csv"
)
bg_img = datasets.load_mni152_template(resolution=0.5)

seq = "hipp"
basePath = Path("/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/")
dataPath = basePath / "data/fMRI/"
spmMB = Path(
    "/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/"
    + seq
    + "/bids/derivatives/spmMB"
)
ppiconRet2 = (
    spmMB
    / "2nd_level"
    / "PPI_noDur"
    / "ConnProject_hipp"
    / "results"
    / "secondlevel"
    / "gPPI_face"
    / "diff_d2d1_Accd2"
    / "rHippHead.rHippBody.rHippTail"
    / "spmT_0001.nii"
)
Py = (
    spmMB
    / "/sub-60601/preproc/ses-1/anat/sub-60601_ses-1_acq-t1csmp2ragesag06mmUNIDEN_T1w.nii"
)
Pnative = (
    Workspace
    / "hipp/searchlight/sub-60601/nativeunconscious_object_categoryRetrieval_sub-60601_idx0.nii"
)

df = pd.read_csv(deriv / "behav/df_hipp_raw.csv")
mask = "/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/wb/bids/derivatives/fmriprep/sub-62426/anat/mask_GreyAndWhiteMattersub-62426.nii"

motor_images = datasets.fetch_neurovault_motor_task()
stat_img = motor_images.images[0]
fsaverage = datasets.fetch_surf_fsaverage()

# Read in your nifti (.nii/.nii.gz) file you want as the background
# Or load a common template is are included with nilearn

# bg_img = bg_image
# Read in your nifti file (same dimentions as bg_img) that want as the overlay.
img = image.load_img(fname)

# Or load a overlay this is are included with nilearn
img = datasets.load_mni152_template(1)

mydf = pd.read_csv(behav_correlationdt)


seq = "wb"
basePath = Path("/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/")
dataPath = basePath / "data/fMRI/"
spmMB = Path(
    "/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/"
    + seq
    + "/bids/derivatives/spmMB"
)
Workspace = Path(
    "/storage/workspaces/psy_memory_wfg_psy/hpc_henke_wfg/s2019_twillems_fMRI_silent_engram/data/fMRI/"
)
fmriPrep = Path(
    "/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/"
    + seq
    + "/bids/derivatives/fmriprep/"
)

spmtfile = (
    spmMB
    / "RFX_Correlations"
    / "RFX_l1_mni_delret_subsequent_forgetexp_s6wr_ret4_maxTrials"
    / "CORREL_later_consc_stays_unconsc_X_accuracy"
    / "spmT_0001.nii"
)
spmtfilefMRIPREP = '/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/wb/bids/derivatives/spmMB/RFX_Correlations/RFX_l1_mni_delret_subsequent_forgetexp_s6wrfmriPrepAROMA_s6ret4fmriPrepAROMAmaxTrials/CORREL_later_consc_stays_unconsc_X_accuracy/spmT_0001.nii'
fsaverage = datasets.fetch_surf_fsaverage()

# %%
# plot 1
sns.set_theme(style="white", palette="binary", font_scale=1.5)  # font_scale=1.5)


fig, (ax0, ax1) = plt.subplots(
    nrows=1, ncols=2, figsize=(14, 4)
)  # sharex=True,gridspec_kw={'width_ratios': [1, 1]}


statimg, thresh = threshold_stats_img(spmtfilefMRIPREP, alpha=0.001, cluster_threshold=10)  # , height_control="fdr")
print(thresh)
subseq_recog = plot_stat_map(
    statimg,
    bg_img=bg_img,
    threshold=thresh,  # 1?
    display_mode="x",  #'yx',
    #title="24h retrieval: inaccessible > unavailable",  # names[i],
    colorbar=True,
    draw_cross=False,
    dim=-0.25,
    cut_coords=[2],#cut_coords=[-28, 2, 15],
    axes=ax0,
    # cmap = 'inferno',
    black_bg=False,  # True,
)


seq = "hipp"
basePath = Path("/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/")
dataPath = basePath / "data/fMRI/"
spmMB = Path(
    "/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/"
    + seq
    + "/bids/derivatives/spmMB"
)

alt_contrast_file_recog = (
    spmMB / "2nd_level" / "contrasts" / "fullFactorial_Face" / "spmT_0197.nii"
)
contrast_file_recog = spmMB / '2nd_level' / 'contrasts' / 'FullFactorial_Recognition' / 'spmT_0009.nii'
contrast_file1, thresh2 = threshold_stats_img(
    contrast_file_recog, alpha=0.005, cluster_threshold=10
)
recog_plot = plot_stat_map(
    contrast_file_recog,
    bg_img=bg_img,
    threshold=thresh2,
    display_mode="x",
    #title="recognition: correct guesses > correct sure",
    colorbar=True,
    draw_cross=False,
    dim=-0.25,
    cut_coords=[-28],#cut_coords=[-32, -28, 37],
    axes=ax1,
    black_bg=False,
)

ax0.text(
    -0.001,
    1,
    "A",
    transform=ax0.transAxes,
    fontsize=16,
    fontweight="bold",
    va="top",
    ha="right",
)
ax1.text(
    -0.001,
    1,
    "B",
    transform=ax1.transAxes,
    fontsize=16,
    fontweight="bold",
    va="top",
    ha="right",
)


# Write "t-values" in the bottom right corner of the second row, second column (ax3)
ax0.text(
    1,
    -0.05,
    "t-values",
    ha="right",
    va="bottom",
    transform=ax0.transAxes,
    fontsize=12,
)
# Write "t-values" in the bottom right corner of the second row, second column (ax3)
ax1.text(
    1,
    -0.05,
    "t-values",
    ha="right",
    va="bottom",
    transform=ax1.transAxes,
    fontsize=12,
)

width, height = 6 / 10, 1.8 / 10  # Fraction of the image dimensions
x, y = 0.0425, 0.37  # Bottom-left corner (normalized coordinates)
corner_radius = 0.05  # Rounded corners (in normalized units)

recog_plot.axes[-28.0].ax.set_xlim(-70, 30)
recog_plot.axes[-28.0].ax.set_ylim(-50, 30)

# %% Figure alone:


fig, ax0New = plt.subplots(
    nrows=1, ncols=1, figsize=(7, 4)
)  # sharex=True,gridspec_kw={'width_ratios': [1, 1]}
subseq_recog = plot_stat_map(
    statimg,
    bg_img=bg_img,
    threshold=thresh,  # 1?
    display_mode="x",  #'yx',
    # title="24h retrieval: inaccessible > unavailable",  # names[i],
    colorbar=True,
    draw_cross=False,
    dim=-0.25,
    axes=ax0New,
    cut_coords=[2],  # cut_coords=[-28, 2, 15],
    # cmap = 'inferno',
    black_bg=False,  # True,
)
ax0New.text(
    1,
    -0.05,
    "t-values",
    ha="right",
    va="bottom",
    transform=ax0New.transAxes,
    fontsize=12,
)

outdir = home / "TOAM" / "paper" / "src" / "figures" / "FigureS1_inaccessible_subsequent.png"
plt.savefig(outdir, dpi=300, bbox_inches="tight")
# %%
