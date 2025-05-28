# %%
import os, sys
from pathlib import Path

home = Path(os.environ["HOME"])
sys.path.append(str(home))
from TOAM.fmri.mvpa.mvpa_utils.get_datasets import *
from nilearn import datasets, image, datasets
from nilearn.plotting import plot_stat_map
from nilearn.glm import threshold_stats_img
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

seq = "hipp"
basePath = Path("/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/")
dataPath = basePath / "data/fMRI/"
spmMB = Path(
    "/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/"
    + seq
    + "/bids/derivatives/spmMB"
)

bg_img = datasets.load_mni152_template(resolution=0.5)

conj_contrast_1 = (
    spmMB
    / "2nd_level"
    / "contrasts"
    / "fullFactorial_Face"
    / "conjunction_day_1_spm.nii"  # / "conjunction conscious unconscious day 1.nii"
)
conj_contrast_2 = (
    spmMB
    / "2nd_level"
    / "contrasts"
    / "fullFactorial_Face"
    / "conjunction_day_2_spm.nii"  # / "conjunction conscious unconscious day 2.nii"
)


# %%
sns.set_theme(
    style="white",
    palette="summer",
    # font_scale=1.5,
)
fig, (ax0, ax1) = plt.subplots(
    nrows=1, ncols=2, figsize=(12, 3)
)  # sharex=True,gridspec_kw={'width_ratios': [1, 1]}

conj_contrast_1, thresh1 = threshold_stats_img(
    conj_contrast_1, alpha=0.001, cluster_threshold=10
)
conj_contrast_2, thresh2 = threshold_stats_img(
    conj_contrast_2, alpha=0.005, cluster_threshold=10
)

ax0.text(
    -0.1,
    1,
    "A",
    transform=ax0.transAxes,
    fontsize=14,
    fontweight="bold",
    va="top",
    ha="right",
)
ax1.text(
    -0.1,
    1,
    "B",
    transform=ax1.transAxes,
    fontsize=14,
    fontweight="bold",
    va="top",
    ha="right",
)

width, height = 6.25 / 10, 1.85 / 10  # Fraction of the image dimensions
x, y = -0.0225, 0.375  # Bottom-left corner (normalized coordinates)
corner_radius = 0.05  # Rounded corners (in normalized units)

conj1 = plot_stat_map(
    conj_contrast_1,
    bg_img=bg_img,
    threshold=thresh1,
    display_mode="x",
    title="30min retrieval conjunction",  #: sure correct > guess incorrect ∧ guess correct > guess incorrect",  # names[i],
    colorbar=True,
    draw_cross=False,
    figure=fig,
    axes=ax0,
    cut_coords=[29],
    dim=-0.25,
    black_bg=False,
    resampling_interpolation="continuous",
)

# below: hiess ehemals con1
conj2 = plot_stat_map(
    conj_contrast_2,
    bg_img=bg_img,
    threshold=thresh2,
    display_mode="x",
    title="24h retrieval conjunction (p<0.005)",  #: sure correct > guess incorrect ∧ guess correct > guess incorrect",  # names[i],
    colorbar=True,
    draw_cross=False,
    figure=fig,
    axes=ax1,
    cut_coords=[26],
    dim=-0.25,
    black_bg=False,
    resampling_interpolation="continuous",
)


conj1.axes[29.0].ax.set_xlim(-60, 30)
conj1.axes[29.0].ax.set_ylim(-40, 20)
conj2.axes[26.0].ax.set_xlim(-60, 30)
conj2.axes[26.0].ax.set_ylim(-40, 30)


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


color = "limegreen"

plt.savefig(
    home / "TOAM" / "paper" / "src" / "figures" / "Figure3ab_conjunction.png",
    dpi=333,
    bbox_inches="tight",
)
plt.show()

# %%
