# %% imports
import os, sys

sys.path.append(os.environ["HOME"])
from pathlib import Path
from TOAM.fmri.mvpa.mvpa_utils.get_datasets import *
from nilearn import plotting, datasets
from nilearn.plotting import plot_stat_map
from nilearn.glm import threshold_stats_img
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from scipy.stats import pearsonr, zscore
from matplotlib import cm
from matplotlib.transforms import Affine2D
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset


sns.set_context("notebook", rc={"lines.linewidth": 2}, font_scale=1.5)
threshz = 3

# %% paths
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
    / "spmT_0001.nii"  # / "con_0001.nii"
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
home = Path(os.environ["HOME"])
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

threshold = "99%"
contrast_images = [fname, correlation_cluster, con1, con2, con3, con4]
### PPI

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
ppi = fmriPrep / "2nd_level" / "PPI_noDurations" / "new"
ppiconEnc = (
    ppi
    / "conn_project_supraliminal_wb_fmriprep_maxTrials_noDurations"
    / "results"
    / "secondlevel"
    / "gPPIFaceHipp"
    / "EncUncAccBase"
    / "rightHeadHipp"
    / "spmT_0001.nii"
)
ppiconRet = (
    ppi
    / "conn_project_supraliminal_wb_fmriprep_maxTrials_noDurations"
    / "results"
    / "secondlevel"
    / "gPPIFaceHipp"
    / "day1uncAccBase"
    / "rightHeadHipp.rightBodyHipp.rightTailHipp"
    / "spmT_0001.nii"
)
seed = ppi / "ROIs" / "wb" / "ses-1" / "sub-60100" / "w232_mask_rh.nii"
seed_ret1 = ppi / "ROIs" / "wb" / "ses-1" / "sub-60100" / "w226_mask_rh.nii"
seed_ret = dataPath / "utilities" / "ROI_masks" / "SPM_Hippocampus_mask_bilat.nii"

ppidfEnc = pd.read_csv(
    f'{os.environ["HOME"]}/TOAM/data/Figure4_Encoding_Connectivity_anterior_rightHipp.csv',
    sep=";",
    header=0,
)
ppidfEnc["ret2_z"] = zscore(ppidfEnc["unconsc Accuracy Day 1"])
ppidfEnc["betas_z"] = zscore(ppidfEnc["ACC and mPFC"])
ppidfEnc.drop(np.where(np.abs(ppidfEnc["ret2_z"]) > threshz)[0], inplace=True)
ppidfEnc.drop(np.where(np.abs(ppidfEnc["betas_z"]) > threshz)[0], inplace=True)

ppidfRet = pd.read_csv(
    f'{os.environ["HOME"]}/TOAM/data/Figure4_Retrieval1_Connectivity_whole_rightHipp.csv',
    sep=";",
    header=0,
)
ppidfRet["ret2_z"] = zscore(ppidfRet["unconsc Accuracy Day 1"])
ppidfRet["betas_z"] = zscore(ppidfRet["right ACC and mPFC"])
ppidfRet.drop(np.where(np.abs(ppidfRet["ret2_z"]) > threshz)[0], inplace=True)
ppidfRet.drop(np.where(np.abs(ppidfRet["betas_z"]) > threshz)[0], inplace=True)

names = [
    "Unconscious Trials: Correct vs. Incorrect",
    "Guessed Trials: Correct vs. Incorrect",
    "correct_unconsc_vs_incorrect_unconsc_2",
    "face_consc_vs_unconsc",
    "face-unconsc2-vs-number",
    "correct_unconsc_vs_incorrect_unconsc_CORR_unconscRet3",
]
thresholds = [6.0]
hemisphere = ["right", "right", "left"]

chords = [
    [30, -28.4, -9.6],
    [26.8, -7.6, -24.8],
    [32, -38, 3],
    [-25, -37.3, -3.96],
    [-25, -37.3, -3.96],
    [-25, -37.3, -3.96],
]
i = 1
# for i in range(0,len(contrast_images[:2])):
contrast_image = contrast_images[i]


# %% Settings
sns.set_theme(style="white", palette="summer", font_scale=1.25)
# plt.tight_layout()
plt.subplots_adjust(
    left=0.05, right=0.95, top=0.95, bottom=0.05, wspace=0.2, hspace=0.2
)
viridis = cm.get_cmap("viridis")
max_color = viridis(0)


contrast_image, thresh1 = threshold_stats_img(
    contrast_image, alpha=0.005, cluster_threshold=20
)
ppiEncTresh, encThresh2 = threshold_stats_img(
    ppiconEnc, alpha=0.001, cluster_threshold=30
)
ppistatRet, retThresh_ppi = threshold_stats_img(
    ppiconRet, alpha=0.001, cluster_threshold=30
)


# %% Figure alternative with different
fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(
    nrows=2, ncols=2, gridspec_kw={"width_ratios": [1, 1]}, figsize=(12, 10)
)  # sharex=True,

plt.subplots_adjust(
    hspace=0.45,
)  # top=0.99, bottom=0.01,wspace=0.4

# Create a grey overlay covering the entire plot
width, height = 6.25 / 10, 1.85 / 10  # Fraction of the image dimensions
x, y = -0.0225, 0.375  # Bottom-left corner (normalized coordinates)
corner_radius = 0.05  # Rounded corners (in normalized units)


ax0.text(
    -0.15,
    1.15,
    "A",
    transform=ax0.transAxes,
    fontsize=14,
    fontweight="bold",
    va="top",
    ha="right",
)
ax2.text(
    -1.35,
    -0.25,
    "B",
    transform=ax1.transAxes,
    fontsize=14,
    fontweight="bold",
    va="top",
    ha="right",
)


sns.regplot(
    x="ret2_z",
    y="betas_z",
    data=ppidfEnc,
    fit_reg=True,
    ax=ax0,
    line_kws={
        "color": max_color,
    },  # Line color
    scatter_kws={"color": max_color, "s": 110, "alpha": 0.5},
)

ax0.set_xlabel("30min retrieval guessing accuracy (z-scored)")
ax0.set_ylabel("Seed-based conn. (z-scored)")
ax0.set_ylim(-3, 3)
x = ppidfEnc["ret2_z"]
y = ppidfEnc["betas_z"]
t = pearsonr(x, y)
ax0.annotate(
    f"R = {t[0]:.2f}",#, p= {t[1]:.5f}",
    xy=(0.05, 0.9),
    xycoords="axes fraction",
    fontsize=12,
)
plot_stat_map(
    ppiEncTresh,
    bg_img=bg_img,
    threshold=encThresh2,  # 1?
    display_mode="x",  #'yx',
    # title="Encodig", #names[i],
    colorbar=True,
    draw_cross=False,
    figure=fig,
    axes=ax1,
    dim=-0.25,
    cut_coords=[12],  # cut_coords=[-47,-51,-2],
    # cmap = 'inferno',
    black_bg=False,  # True,
)

# arr_img = plt.imread('rightTailHipp.png')
arr_img = plt.imread(f'{os.environ["HOME"]}/TOAM/data/rightAntHippSeed.png')

imagebox = OffsetImage(arr_img, zoom=0.35)
imagebox.image.axes = ax1

xy = [0.23, 0.75]
ab3 = AnnotationBbox(
    imagebox,
    xy,
    xybox=(100, -150),  # xybox= (.64, 0.25),#(.3, xy[1]),
    xycoords="data",
    boxcoords="offset points",
)

ax1.add_artist(ab3)
ab3.set_zorder(0)
destinations = [(0.58, 0.44)]  # Ç (0.23, 0.75),(0.23,0.45)]
styles = ["arc3,rad=0.3", "arc3,rad=-0.1", "arc3,rad=-0.3"]

circle1 = Circle((0.57, 0.48), 0.06, lw=3, color=max_color, fill=False, zorder=0)
ax1.add_patch(
    circle1,
)

# Create arrows to different destinations

source = [0.57, 0.025]
arrows = []
for ind in range(0, len(destinations)):
    dest = destinations[ind]
    style = styles[ind]
    dx = dest[0] - source[0]
    dy = dest[1] - source[1]
    arrow = FancyArrowPatch(
        source,
        dest,
        connectionstyle=style,
        lw=4,
        arrowstyle="->",
        color="red",
        mutation_scale=15,
        alpha=0.35,
    )  # "->"
    arrows.append(arrow)

# Add the arrows to the plot
for arrow in arrows:
    arrow.set_zorder(10)
    ax1.add_patch(arrow)

fig.text(
    0.1,
    0.90,
    "Encoding connectivity for later correct > incorrect guess • 30min retrieval guessing accuracy",
    ha="left",
    va="center",
    # pad=10,
    fontweight="bold",
    fontsize=12,
)

sns.regplot(
    x="ret2_z",
    y="betas_z",
    data=ppidfRet,
    fit_reg=True,
    ax=ax2,
    line_kws={
        "color": max_color,
    },  # Line color
    scatter_kws={"color": max_color, "s": 110, "alpha": 0.5},
)
ax2.set_xlabel("30min retrieval guessing accuracy (z-scored)")
ax2.set_ylabel("Seed-based conn. (z-scored)")
ax2.set_yticks(range(-3, 3))
x = ppidfRet["unconsc Accuracy Day 1"]
y = ppidfRet["right ACC and mPFC"]
t = pearsonr(x, y)
ax2.annotate(
    f"R = {t[0]:.2f}",#, p= {t[1]:.5f}",
    xy=(0.05, 0.9),
    xycoords="axes fraction",
    fontsize=12,
)

plot_stat_map(
    ppistatRet,
    bg_img=bg_img,
    threshold=retThresh_ppi,  # 1?
    # display_mode="tiled",  #'yx',
    display_mode="x",
    # title="Retrieval Day 1", #names[i],
    colorbar=True,
    draw_cross=False,
    figure=fig,
    axes=ax3,
    dim=-0.25,
    # cut_coords=[8, 40, 20],
    cut_coords=[6.1],
    # cmap = 'inferno',
    black_bg=False,  # True,
)
arr_img = plt.imread(f'{os.environ["HOME"]}/TOAM/data/rightHippSeed.png')
imagebox = OffsetImage(arr_img, zoom=0.35)
imagebox.image.axes = ax3


xy = [0.23, 0.75]  # (0.66, 0.25)  # xy =

ab2 = AnnotationBbox(
    imagebox,
    xy,
    xybox=(100, -150),
    xycoords="data",
    boxcoords=("offset points"),
    pad=0.5,
)
# annotation_clip = False,
# bboxprops=dict(alpha=0.5))
ax3.add_artist(ab2)


circle1 = Circle((0.55, 0.56), 0.07, lw=3, color=max_color, fill=False, zorder=0)

ax3.add_patch(
    circle1,
)

source = [0.57, 0.025]
arrows = []
destinations = [(0.54, 0.50)]  # Ç (0.23, 0.75),(0.23,0.45)]
for ind in range(0, len(destinations)):
    dest = destinations[ind]
    style = styles[ind]
    dx = dest[0] - source[0]
    dy = dest[1] - source[1]
    arrow = FancyArrowPatch(
        source,
        dest,
        connectionstyle=style,
        lw=4,
        arrowstyle="->",
        color="red",
        mutation_scale=15,
        alpha=0.35,
    )  # "->"
    arrows.append(arrow)

# Add the arrows to the plot
for arrow in arrows:
    arrow.set_zorder(10)
    ax3.add_patch(arrow)

fig.text(
    0.1,
    0.45,
    "30min retrieval connectivity for correct > incorrect guess • 30min retrieval guessing accuracy",
    ha="left",
    va="center",
    # pad=10,
    fontweight="bold",
    fontsize=12,
)
ax1.text(
    1, -0.05, "t-values", ha="right", va="bottom", transform=ax1.transAxes, fontsize=12
)

# Write "t-values" in the bottom right corner of the second row, second column (ax3)
ax3.text(
    1,
    -0.05,
    "t-values",
    ha="right",
    va="bottom",
    transform=ax3.transAxes,
    fontsize=12,
)


plt.savefig(
    home / "TOAM" / "paper" / "src" / "figures" / "Figure4_gPPI_Hipp_prefrontal.png",
    bbox_inches="tight",
    dpi=333,
)

# %%
