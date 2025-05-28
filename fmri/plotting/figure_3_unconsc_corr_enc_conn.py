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
# df_corr = pd.read_csv(fmriPrep /'/sub-69632/ses-1/func/sub-69632_ses-1_task-forgetexp_events_corrected.tsv',sep='\t')
# df_maxT = pd.read_csv(fmriPrep /'/sub-69632/ses-1/func/sub-69632_ses-1_task-forgetexp_events_maxTrials.tsv',sep='\t')
# df_orig = pd.read_csv(fmriPrep /'/sub-69632/ses-1/func/sub-69632_ses-1_task-forgetexp_events.tsv',sep='\t')
df = pd.read_csv(deriv / "behav/df_hipp_raw.csv")
mask = "/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/wb/bids/derivatives/fmriprep/sub-62426/anat/mask_GreyAndWhiteMattersub-62426.nii"


mydf = pd.read_csv(behav_correlationdt)

mydf["ret2_z"] = zscore(mydf["ret2"])
mydf["betas_z"] = zscore(mydf["betas_hipp_ret1"])
mydf.drop(np.where(np.abs(mydf["ret2_z"]) > threshz)[0], inplace=True)
mydf.drop(np.where(np.abs(mydf["betas_z"]) > threshz)[0], inplace=True)


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
    / "rightTailHipp"
    / "spmT_0001.nii"
)
seed = ppi / "ROIs" / "wb" / "ses-1" / "sub-60100" / "w232_mask_rh.nii"
seed_ret1 = ppi / "ROIs" / "wb" / "ses-1" / "sub-60100" / "w226_mask_rh.nii"
seed_ret = dataPath / "utilities" / "ROI_masks" / "SPM_Hippocampus_mask_bilat.nii"


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
    ppiconEnc, alpha=0.005, cluster_threshold=20
)
ppistatRet, retThresh_ppi = threshold_stats_img(
    ppiconRet, alpha=0.005, cluster_threshold=20
)


# %% FIGURE 2 ##############################
fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(
    nrows=2, ncols=2, gridspec_kw={"width_ratios": [1, 1]}, figsize=(12, 12)
)  # sharex=True,

plt.subplots_adjust(
    hspace=0.3,
)  # top=0.99, bottom=0.01,wspace=0.4

ax0.text(
    -0.15,
    1.1,
    "C",
    transform=ax0.transAxes,
    fontsize=14,
    fontweight="bold",
    va="top",
    ha="right",
)

# Add shared titles for each row
fig.text(
    0.1,
    0.90,
    "30min retrieval correct > incorrect guess • 30min retrieval guessing accuracy",
    ha="left",
    va="center",
    # pad=10,
    fontweight="bold",
    fontsize=12,
)
fig.text(
    0.1,
    0.46,
    "24h retrieval correct > incorrect guess • 24h retrieval guessing accuracy",
    ha="left",
    va="center",
    # pad=10,
    fontweight="bold",
    fontsize=12,
)

sns.regplot(
    x="ret2_z",
    y="betas_z",
    data=mydf,
    fit_reg=True,
    ax=ax0,
    line_kws={
        "color": max_color,
    },  # Line color
    scatter_kws={"color": max_color, "s": 110, "alpha": 0.5},
)
ax0.set_xlabel("30min retrieval guessing accuracy (z-scored)")
ax0.set_ylabel("Estimated responses (betas, z-scored)")
x = mydf["ret2_z"]
y = mydf["betas_z"]
t = pearsonr(x, y)
ax0.annotate(
    f"R = {t[0]:.2f}",# p= {t[1]:.6f}",
    xy=(0.05, 0.9),
    xycoords="axes fraction",
    fontsize=12,
)


circleOne = Circle((0.37, 0.39), 0.06, lw=3, color=max_color, fill=False, zorder=0)
ax1.add_patch(
    circleOne,
)


ax0.patch.set_edgecolor("black")
ax0.patch.set_linewidth(1)

# brainplot
day1 = plot_stat_map(
    contrast_image,
    bg_img=bg_img,
    threshold=thresh1,  # 1?
    display_mode="x",
    colorbar=True,
    draw_cross=False,
    figure=fig,
    axes=ax1,
    # title="Retrieval 1",
    cut_coords=[27],  # chords[i][0]
    # cmap = 'inferno',
    dim=-0.25,
    black_bg=False,  # True,
    resampling_interpolation="continuous",
)
day1.axes[27.0].ax.set_xlim(-60, 30)
day1.axes[27.0].ax.set_ylim(-40, 20)

mydf_ret2 = pd.read_csv(behav_correlationdt)

ret2betas = pd.read_csv(
    f'{os.environ["HOME"]}/TOAM/data/Figure3_retrieval2_rightanthipp_cluster_params.csv',
    sep=";",
    header=0,
)

mydf_ret2["betas_hipp_ret2"] = ret2betas["Var2"]
mydf_ret2["ret3_z"] = zscore(mydf_ret2["ret3"])
mydf_ret2["betas_z"] = zscore(mydf_ret2["betas_hipp_ret2"])
mydf_ret2.drop(np.where(np.abs(mydf_ret2["ret3_z"]) > threshz)[0], inplace=True)
mydf_ret2.drop(np.where(np.abs(mydf_ret2["betas_z"]) > threshz)[0], inplace=True)

cluster_img, trsh = threshold_stats_img(
    correlation_cluster_day2, alpha=0.005, cluster_threshold=10
)
outname = (
    home / "TOAM" / "paper" / "src" / "figures" / "fMRI_UniFuncCorrDay2_corrShown.png"
)
sns.regplot(
    x="ret3_z",
    y="betas_z",
    data=mydf_ret2,
    fit_reg=True,
    ax=ax2,
    line_kws={
        "color": max_color,
    },  # Line color
    scatter_kws={"color": max_color, "s": 110, "alpha": 0.5},
)
ax2.set(xlabel=None)
ax2.set_ylabel("Estimated responses (betas, z-scored)")
ax2.set_xlabel("24h retrieval guessing accuracy (z-scored)")
x = mydf_ret2["ret3_z"]
y = mydf_ret2["betas_z"]
t = pearsonr(x, y)
ax2.annotate(
    f"R = {t[0]:.2f}",#, p= {t[1]:.4f}",
    xy=(0.05, 0.9),
    xycoords="axes fraction",
    fontsize=12,
)

circle1 = Circle((0.37, 0.38), 0.05, lw=3, color=max_color, fill=False, zorder=0)
ax3.add_patch(
    circle1,
)

day2 = plot_stat_map(
    cluster_img,
    bg_img=bg_img,
    threshold=trsh,
    display_mode="x",
    colorbar=True,
    draw_cross=False,
    cut_coords=[27],
    dim=-0.25,
    figure=fig,
    axes=ax3,
    black_bg=False,
)
day2.axes[27.0].ax.set_xlim(-60, 30)
day2.axes[27.0].ax.set_ylim(-40, 20)

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
    home
    / "TOAM"
    / "paper"
    / "src"
    / "figures"
    / "Figure3c_correlationHipp30min_24h.png",
    bbox_inches="tight",
    dpi=333,
)

# %%
