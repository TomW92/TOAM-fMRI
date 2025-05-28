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
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
from scipy.stats import norm, gaussian_kde, pearsonr, zscore
from matplotlib import cm
from matplotlib.transforms import Affine2D

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
bg_img = datasets.load_mni152_template(resolution=.5)

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
    / "rHippHead" #/ "rHippHead.rHippBody.rHippTail"
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

viridis = cm.get_cmap("viridis")
max_color = viridis(0)

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
fsaverage = datasets.fetch_surf_fsaverage()

# %%
# plot 1
sns.set_theme(style="white", palette="binary", font_scale=1.2)  # font_scale=1.5)


spmTfile = (
    spmMB / "2nd_level" / "contrasts" / "fullFactorial_Face" / "spmT_0056.nii"
)  #'spmT_0028.nii'
spmTfile = (
    fmriPrep
    / "2nd_level"
    / "contrasts_noDurations"
    / "fullFactorial_Face"
    / "spmT_0056.nii"
)
spmTfile2 = spmMB / "2nd_level" / "contrasts" / "fullFactorial_Face" / "spmT_0029.nii"
spmTfile2 = (
    fmriPrep
    / "2nd_level"
    / "contrasts_noDurations"
    / "Doppelkontraste"
    / "spmT_0015.nii"
)


spm1con, spm1thresh = threshold_stats_img(spmTfile, alpha=0.001, cluster_threshold=10)
spm2con, spm2thresh = threshold_stats_img(spmTfile2, alpha=0.005, cluster_threshold=10)


ppidf = pd.read_csv(
    f"{home}/TOAM/data/Figure5_Retrieval2_vs_Retrieval1_HippBody.csv", sep=";", header=0
)
ppidf["ret2_z"] = zscore(ppidf["unconscious Accuracy Day 2"])
ppidf["betas_z"] = zscore(ppidf["left Hipp"])
ppidf.drop(np.where(np.abs(ppidf["ret2_z"]) > 3)[0], inplace=True)
ppidf.drop(np.where(np.abs(ppidf["betas_z"]) > 3)[0], inplace=True)


# %% plot
from nilearn.plotting import plot_surf_roi, plot_surf, plot_surf_roi

hipp_path = "/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/utilities/ROI_masks/SPM_Hippocampus_mask_bilat.nii"
hipp_mask = image.load_img(hipp_path)

# %% ppi alone
fig, (axx1, axx2) = plt.subplots(
    nrows=1, ncols=2, gridspec_kw={"width_ratios": [1, 1]}, figsize=(12, 6)
)  # sharex=True,


sns.regplot(
    x="ret2_z",
    y="betas_z",
    data=ppidf,
    fit_reg=True,
    ax=axx1,
    line_kws={
        "color": max_color,
    },  # Line color
    scatter_kws={"color": max_color, "s": 110, "alpha": 0.5},
)
axx1.set_xlabel("24h retrieval guessing accuracy (z-scored)")
axx1.set_ylabel("Seed-based conn. (z-scored)")
x = ppidf["ret2_z"]
y = ppidf["betas_z"]
t = pearsonr(x, y)
axx1.annotate(
    f"R = {t[0]:.2f}",#, p= {t[1]:.4f}",
    xy=(0.05, 0.9),
    xycoords="axes fraction",
    fontsize=12,
)

ppiconRet2_stat, threshppricon2 = threshold_stats_img(
    ppiconRet2, alpha=0.001, cluster_threshold=10
)

plot_stat_map(
    ppiconRet2_stat,
    bg_img=bg_img,
    threshold=threshppricon2,  # 1?
    display_mode="y",  #'yx',
    # title="Retrieval Day 1", #names[i],
    colorbar=True,
    draw_cross=False,
    figure=fig,
    axes=axx2,
    dim=-0.25,
    cut_coords=[-28],
    # cmap = 'inferno',
    black_bg=False,  # True,
)
arr_img = plt.imread(f"{home}/TOAM/data/rightAntHippSeed.png")
imagebox = OffsetImage(arr_img, zoom=0.5)
imagebox.image.axes = axx1
xy = [0.5, 0.5]
ab3 = AnnotationBbox(
    imagebox,
    xy,
    xybox=(95, -112.5),  # (.3, xy[1]),
    xycoords="data",
    # box_alignment=(0, 3),
    boxcoords=("offset points"),  #
    pad=0.5,
)

axx1.add_artist(ab3)

axx1.spines["right"].set_visible(False)
axx1.spines["top"].set_visible(False)
cluster = (0.23, 0.45)
circle1 = Circle(cluster, 0.05, lw=3, color=max_color, fill=False, zorder=0)
axx2.add_patch(
    circle1,
)

destinations = [(0.63, 0.46)]#, (0.72, 0.44)]  # Ç (0.23, 0.75),(0.23,0.45)]
styles = ["arc3,rad=0.3", "arc3,rad=0.1", "arc3,rad=-0.3"]
# Create arrows to different destinations

source = [0.45, 0.24]
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
    # axx1.add_patch(arrow)
    fig.add_artist(arrow)

width, height = 6 / 10, 1.8 / 10  # Fraction of the image dimensions
x, y = 0.0425, 0.37  # Bottom-left corner (normalized coordinates)
corner_radius = 0.05  # Rounded corners (in normalized units)

# Create a rectangle with rounded edges
rectangle = FancyBboxPatch(
    (x, y),
    width,
    height,
    boxstyle=f"round,pad=0.02,rounding_size={corner_radius}",
    edgecolor="black",
    linewidth=3,
    facecolor="none",
    alpha=1,
    zorder=10,
)

textx = x + width + 0.05
texty = y + height / 2

axx2.annotate(
    "FOV",
    xy=(textx, texty),  # Center of the rectangle
    ha="left",
    va="center",
    fontsize=20,
    color="black",
    transform=1,
)

axx2.add_patch(rectangle)
# Write "t-values" in the bottom right corner of the second row, second column (ax3)
axx2.text(
    1,
    -0.05,
    "t-values",
    ha="right",
    va="bottom",
    transform=axx2.transAxes,
    fontsize=12,
)

plt.savefig(
    home
    / "TOAM"
    / "paper"
    / "src"
    / "figures"
    / "Figure_5_gPPI_Hipp_bilat.png",
    dpi=333,
    bbox_inches="tight",
)
plt.show()

# %%
