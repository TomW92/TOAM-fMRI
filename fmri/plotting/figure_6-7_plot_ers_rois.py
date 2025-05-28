# %% imports: compute_erst_results.py
import sys, os
from pathlib import Path
from datetime import datetime
import pandas as pd
from nilearn.plotting import plot_anat, show, view_img

sys.path.append(str(Path.home()))
from TOAM.fmri.mvpa.mvpa_utils.get_datasets import setup_configNL  # type: ignore
import time
import matplotlib.pyplot as plt
import nilearn as nl
import numpy as np
import seaborn as sns
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import zscore, ks_2samp, pearsonr, ttest_1samp
import matplotlib.colors as mcolors

temp = nl.datasets.load_mni152_template()

freesurfer_dir = Path(
    f"/storage/workspaces/psy_memory_wfg_psy/hpc_henke_wfg/s2019_twillems_fMRI_silent_engram/data/fMRI/hipp/freesurfer/ses-1/sub-60601/mri"
)

# %% plot 1
masks = [
    "hipp_head_rh",
    "hipp_body_rh",
    "hipp_tail_rh",
    "parahipp_rh_no_er_mask",
    "entorhinal_rh_mask",
]

colors = [
    (230 / 255, 159 / 255, 0 / 255),  # Orange
    (86 / 255, 180 / 255, 233 / 255),  # Sky Blue
    (0 / 255, 158 / 255, 115 / 255),  # Bluish Green
    (213 / 255, 94 / 255, 0 / 255),  # Vermilion (Red-Orange)
    (240 / 255, 228 / 255, 66 / 255),  # Yellow
]

config = setup_configNL(60601)
anatimg = "/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/hipp/bids/derivatives/spmMB/sub-60601/preproc/ses-1/anat/wskullStrippedsub-60601_ses-1_acq-t1csmp2ragesag06mmUNIDEN_T1w.nii"

display = plot_anat(
    config["spmBiasCorr"], #config['normAnatFiles'][6],#$freesurfer_dir / "T1.mgz",
    display_mode="x",
    cut_coords=[25],
    draw_cross=False,
    dim=-.5,
)

for mask, color in zip(masks, colors):
    display.add_contours(
        freesurfer_dir / config["freesurfer_ses-1"][mask],
        colors=[color],
        filled=True,
    )

# %% plot 2
masks = [
    "hipp_head_rh",
    "hipp_body_rh",
    "hipp_tail_rh",
    "parahipp_mask",
]

colors = [
    (230 / 255, 159 / 255, 0 / 255),  # Orange
    (86 / 255, 180 / 255, 233 / 255),  # Sky Blue
    (0 / 255, 158 / 255, 115 / 255),  # Bluish Green
    (213 / 255, 94 / 255, 0 / 255),  # Vermilion (Red-Orange)
]

# %% plot 2
display2 = plot_anat(
    config["spmBiasCorr"],  
    display_mode="x",
    title="T1",
    cut_coords=[25],
    draw_cross=False,
    dim=-0.5,
)

for mask, color in zip(masks, colors):
    display2.add_contours(
        freesurfer_dir / config["freesurfer_ses-1"][mask],
        colors=[color],
        filled=True,
    )

# %% Save
home = os.environ["HOME"]
display.savefig(f'{home}/TOAM/data/roi_ers_plot.png', dpi=333)
display2.savefig(f"{home}/TOAM/data/roi_ers_plot2.png", dpi=333)
# %%
