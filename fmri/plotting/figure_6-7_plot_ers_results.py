# %% imports: compute_erst_results.py
import sys, os
from pathlib import Path
from datetime import datetime

import pandas as pd
import matplotlib as mpl

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
from pingouin import bayesfactor_ttest, pairwise_corr

# from TOAM.fmri.mvpa.rsa.ers.encoding_retrieval_sim import *
sns.set_context("paper", rc={"lines.linewidth": 2}, font_scale=2)

WB_MAPPING = {
    # "rh_hipp": "R-H",
    "parahipp_rh_mask": "R-PH",
    # "lh_hipp": "L-H",
    "parahipp_lh_mask": "L-PH",
    "inferior_temp_sulc_mask": "ITS",
    "mid_temp_gyr_mask": "MTG",
    "inf_temp_gyr_mask": "ITG",
    "cuneus_mask": "CUN",
    "precuneus_mask": "PCUN",
    "acc_mask": "ACC",
    "front_sup_mask": "FS",
}
WB_MAPPING_DET = {
    # "rh_hipp": "R. Hipp.",
    "parahipp_rh_mask": "R. parahipp. g.",
    # "lh_hipp": "L. Hipp.",
    "parahipp_lh_mask": "L. parahipp g.",
    "inferior_temp_sulc_mask": "B. inf. temp. s.",
    "mid_temp_gyr_mask": "B. mid. temp. g.",
    "inf_temp_gyr_mask": "B. inf. temp. g.",
    "cuneus_mask": "B. cuneus",
    "precuneus_mask": "B. precuneus",
    "acc_mask": "B. ant. cing. g.",
    "front_sup_mask": "B. sup. front. g.",
}
WB_RECOG_MAPPING = {
    # "rh_hipp": "R-H",
    # "lh_hipp": "L-H",
    "parahipp_lh_no_er_mask": "L-PH",
    "parahipp_rh_no_er_mask": "R-PH",
    "entorhinal_lh_mask": "L-EC",
    "entorhinal_rh_mask": "R-EC",
    "inferior_temp_sulc_mask": "ITS",
    "mid_temp_gyr_mask": "MTG",
    "inf_temp_gyr_mask": "ITG",
    "cuneus_mask": "CUN",
    "precuneus_mask": "PCUN",
    "acc_mask": "ACC",
    "front_sup_mask": "FS",
}
WB_RECOG_MAPPING_DET = {
    # "rh_hipp": "R. Hipp.",
    # "lh_hipp": "L. Hipp.",
    "parahipp_lh_no_er_mask": "L. parahipp. g.",
    "parahipp_rh_no_er_mask": "R. parahipp. g.",
    "entorhinal_lh_mask": "L. entorh. c.",
    "entorhinal_rh_mask": "R. entorh. c.",
    "inferior_temp_sulc_mask": "B. inf. temp. s.",
    "mid_temp_gyr_mask": "B. mid. temp. g.",
    "inf_temp_gyr_mask": "B. inf. temp. g.",
    "cuneus_mask": "B. cuneus",
    "precuneus_mask": "B. precuneus",
    "acc_mask": "B. ant. cing. g.",
    "front_sup_mask": "B. sup. front. g.",
}
HIPP_MAPPING = {
    "hipp_head_rh": "R-HH",
    "hipp_body_rh": "R-HB",
    "hipp_tail_rh": "R-HT",
    "rh_hipp": "R-H",
    "parahipp_rh_mask": "R-PH",
    "hipp_head_lh": "L-HH",
    "hipp_body_lh": "L-HB",
    "hipp_tail_lh": "L-HT",
    "lh_hipp": "L-H",
    "parahipp_lh_mask": "L-PH",
    # "inferior_temp_sulc_mask": "ITS",
    # "mid_temp_gyr_mask": "MTG",
    # "inf_temp_gyr_mask": "ITG",
    # "cuneus_mask": "CUN",
    # "precuneus_mask": "PCUN",
}
HIPP_MAPPING_DET = {
    "hipp_head_rh": "R. hipp. head",
    "hipp_body_rh": "R. hipp. body",
    "hipp_tail_rh": "R. hipp. tail",
    "rh_hipp": "R. hipp.",
    "parahipp_rh_mask": "R. parahipp. g.",
    "hipp_head_lh": "L. hipp. head",
    "hipp_body_lh": "L. hipp. body",
    "hipp_tail_lh": "L. hipp. tail",
    "lh_hipp": "L. hipp.",
    "parahipp_lh_mask": "L. parahipp. g.",
    # "inferior_temp_sulc_mask": "B. inf. temp. s.",
    # "mid_temp_gyr_mask": "B. mid. temp. g.",
    # "inf_temp_gyr_mask": "B. inf. temp. g.",
    # "cuneus_mask": "B. cuneus",
    # "precuneus_mask": "B. precuneus",
}
HIPP_RECOG_MAPPING = {
    "hipp_head_rh": "R-HH",
    "hipp_body_rh": "R-HB",
    "hipp_tail_rh": "R-HT",
    "rh_hipp": "R-H",
    "hipp_head_lh": "L-HH",
    "hipp_body_lh": "L-HB",
    "hipp_tail_lh": "L-HT",
    "lh_hipp": "L-H",
    "parahipp_lh_no_er_mask": "L-PH",
    "parahipp_rh_no_er_mask": "R-PH",
    "entorhinal_lh_mask": "L-EC",
    "entorhinal_rh_mask": "R-EC",
    # "mid_temp_gyr_mask": "MTG",
    # "cuneus_mask": "CUN",
    # "precuneus_mask": "PCUN",
}
HIPP_RECOG_MAPPING_DET = {
    "hipp_head_rh": "R. hipp. head",
    "hipp_body_rh": "R. hipp. body",
    "hipp_tail_rh": "R. hipp. tail",
    "rh_hipp": "R. hipp.",
    "parahipp_rh_no_er_mask": "R. parahipp. g.",
    "entorhinal_rh_mask": "R. entorh. c.",
    "hipp_head_lh": "L. hipp. head",
    "hipp_body_lh": "L. hipp. body",
    "hipp_tail_lh": "L. hipp. tail",
    "lh_hipp": "L. hipp.",
    "parahipp_lh_no_er_mask": "L. parahipp. g.",
    "entorhinal_lh_mask": "L. entorh. c.",
    # "mid_temp_gyr_mask": "B. mid. temp. g.",
    # "cuneus_mask": "B. cuneus",
    # "precuneus_mask": "B. precuneus",
}


# Function to rename columns based on the mapping
def rename_column_values(df, column, mapping):
    df1 = df.copy()
    df1[column] = df[column].replace(mapping)
    return df1


def get_limited_palette(cmap, start=0.2, end=0.8, n_colors=0):
    return sns.color_palette(
        cmap(mcolors.Normalize(vmin=start, vmax=end)(range(n_colors)))
    )


## %% define functions
def plot_ers_results(
    data=None,
    datapath=None,
    exp_id="test_exp_id",
    remove_outliers=False,
    corrFile=None,
    plot_outliers=True,
    plot_correlations=False,
    **kwargs,
):
    """
    plots the ERS analysis and performs a significance test
    takes the dataframe corrected above
    """
    printnotagain = False
    ax = "failed to create a subplot?"
    title = ""
    # paper_fig_names = [
    #     "enc_ret1_perm_consc.png",
    #     "__enc_ret2_perm.png",
    #     "enc_ret1_perm.png",
    #     "__enc_ret2_perm.png",
    #     "__enc_ret2_perm_consc.png",
    #     "__ret1_ret2_perm_consc.png",
    #     "__ret1_ret2_perm.png",
    #     "enc_recog_perm_consc.png",
    #     "wb-all_memory_n_enc_recog_perm.png",
    # ]
    paper_fig_names = ["unconsc_incorr.png", "unconsc_incorr_"]

    # load data
    if isinstance(datapath, list):
        df = pd.concat([pd.read_csv(str(f) + ".txt", sep="\t") for f in datapath])
        datapath = str(datapath[0]) + "mixed"
    elif data is not None:
        df = data
        printnotagain = True
    elif datapath is not None:
        df = pd.read_csv(str(datapath) + ".txt", sep="\t")  # ,header=None)
        if "enc" in datapath.split("_") and "ret1" in datapath.split("_"):
            title = "Encoding – 30min retrieval similarity for"
        elif "enc" in datapath.split("_") and "ret2" in datapath.split("_"):
            title = "Encoding - 24h retrieval similarity for"
        elif "ret1" in datapath.split("_") and "ret2" in datapath.split("_"):
            title = "30min – 24h retrieval similarity for"
        elif "recog" in datapath.split("_"):
            title = "Encoding - recognition similarity for"
    else:
        raise ValueError("Either 'dataframe' or 'datapath' must be provided")

    if datapath == "facetGrid_plot":
        pass
    else:
        if str(datapath).startswith("/"):
            if "enc" in datapath.split("_") and "ret1" in datapath.split("_"):
                title = "Encoding - 30min retrieval similarity for"
            elif "enc" in datapath.split("_") and "ret2" in datapath.split("_"):
                title = "Encoding - 24h retrieval similarity for"
            elif "ret1" in datapath.split("_") and "ret2" in datapath.split("_"):
                title = "30min – 24h retrieval similarity for"
            elif "recog" in datapath.split("_"):
                title = "Encoding - recognition similarity for"

            if "consc_consc" in str(datapath):
                title = title + " sure responses"
            if "consc_unconsc" in str(datapath):
                title = title + " guess responses"
        else:
            title = str(datapath).split("-")[-1]
            title = str(title).replace("_", " ").capitalize()

    # outlier detection
    if remove_outliers == True:
        if "consc" in str(datapath):
            df = df[
                (np.abs(df["unconscious"]) < 3.5)
                & (np.abs(df["neither"]) < 3.5)
                & (np.abs(df["conscious"]) < 3.5)
            ]
        else:
            df = df[(np.abs(df["simVal"]) < 3.5)]

    
    correct_results = ['consc_consc-unconsc_incorr','consc_unconsc_corr-unconsc_incorr']
    study = []
    order = []
    t_tests = []
    if "hipp" in str(datapath) or df["subject"].str.contains("sub-60601").any():
        study = "hipp"
        mapping = HIPP_MAPPING_DET
        if (
            "recog" in str(datapath)
            or "L-EC" in df.roi.unique()
            or "L. entorh. c." in df.roi.unique()
        ):
            mapping = HIPP_RECOG_MAPPING_DET
        order = list(mapping.values())
    elif "wb" in str(datapath) or df["subject"].str.contains("sub-60100").any():
        study = "wb"
        mapping = WB_MAPPING_DET
        if (
            "recog" in str(datapath)
            or "L-EC" in df.roi.unique()
            or "L. entorh. c." in df.roi.unique()
        ):
            mapping = WB_RECOG_MAPPING_DET
        order = list(mapping.values())
    else:
        raise ValueError("Unknown study")
    datapath = str(datapath).split("/")[-1]
    OutputDir = Path(
        f"/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/{study}/bids/derivatives/nilearn/ERS/results/{exp_id}"
    )
    OutputDir.mkdir(parents=True, exist_ok=True)
    skip_rest = False

    rois_to_keep = list(mapping.keys())
    df = df[df["roi"].isin(rois_to_keep)]
    df = rename_column_values(df, "roi", mapping)

    if str(datapath).endswith("_consc"):
        numcols = round(len(df["roi"].unique()))
        df = df.dropna()
        dfCopy = df.copy()
        df = df.melt(
            id_vars=list(df.columns[[0, 1]]),
            value_vars=list(df.columns[[2, 3, 4]]),
            var_name="consciousnessLevel",
            value_name="simVal",
        )

        # EINKLAMMERN FÜR OVERALL CONSC RESULTS!
        # fig,ax = plt.subplots(figsize=(numcols,8))
        # if plot_outliers == True:
        #     sns.boxplot(x='roi',y='simVal',hue='consciousnessLevel',data=df,ax=ax)
        # else:
        #     sns.boxplot(x='roi',y='simVal',hue='consciousnessLevel',data=df,ax=ax,showfliers=False)
        # ax.axhline(0, color='r',linestyle='--')
        # ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
        # ax.set_title('ERS Results across ROIs and Consc Levels')
        # plt.savefig(datapath+'.png')

        t_tests = df.groupby(["roi", "consciousnessLevel"], observed=True)[
            "simVal"
        ].apply(lambda x: ttest_1samp(x.dropna(), 0))
        dfCopy.sort_index(
            axis=1, inplace=True
        )  #  so it is always the same order of conscious, unconscious neither, etc. that I have to consider
        columnsOfInterest = [
            col for col in dfCopy.columns if col not in ["subject", "roi"]
        ]
        new_row = pd.DataFrame(
            {
                "subject": [df["subject"].sort_values().iloc[0]],
                "roi": ["no_roi"],
                "simVal": [0],
                "consciousnessLevel": ["plot_roi"],
            }
        )
        # df = pd.concat([df,new_row],axis=0, ignore_index=True)
        tempdf = rename_column_values(
            df,
            "consciousnessLevel",
            {
                "unconscious_correct": "Correct guess",
                "unconscious_incorrect": "Incorrect guess",
                "conscious_correct": "Correct sure",
            },
        )

        g = sns.FacetGrid(
            tempdf,
            col="consciousnessLevel",
            col_wrap=3,
            height=7,
            aspect=1.25,
            despine=True,
            margin_titles=False,
        )

        g.map_dataframe(
            plot_ers_results,
            datapath="facetGrid_plot",
            exp_id=exp_id,
            remove_outliers=False,
            plot_outliers=True,
            plot_correlations=False,
            order=order,
        )
        # g.add_legend(title="ROI", bbox_to_anchor=(1, 0), loc="lower right")
        g.set_titles("{col_name}", size=16)
        g.fig.subplots_adjust(hspace=0.4)
        # g.fig.suptitle(f"{title} pairwise RSA across consciousness levels")
        g.fig.suptitle(
            f"{title} pairwise RSA across consciousness levels", y=1.05, fontsize=16
        )

        for axes in g.axes.flat:
            _ = axes.set_xticks(
                axes.get_xticks(),
                axes.get_xticklabels(),
                rotation=0,
                ha="center",
                rotation_mode="anchor",
            )

        if any(pattern in (datapath + ".png") for pattern in paper_fig_names):
            alt_OutputDir = (
                os.environ["HOME"] + f"/TOAM/paper/src/figures/{datapath}_FACET.png"
            )
            g.savefig(alt_OutputDir, bbox_inches="tight", dpi=333)

        outname = OutputDir / f"{datapath}_FACET.png"
        g.savefig(outname, bbox_inches="tight", dpi=333)
        if corrFile is not None:
            for columnName, corrF in zip(columnsOfInterest, corrFile):
                df = dfCopy[["subject", "roi", columnName]]
                df = df.rename(columns={columnName: "simVal"})
                columnName = datapath + "-" + columnName + "_trials"
                plot_ers_results(
                    data=df,
                    datapath=columnName,
                    exp_id=exp_id,
                    remove_outliers=False,
                    corrFile=corrF,
                    plot_outliers=True,
                    plot_correlations=False,
                )

    else:
        ax = plt.gca()

        if ("consciousnessLevel" in df.columns) and (
            "plot_roi" in df["consciousnessLevel"].unique()
        ):
            skip_rest = True
            roi_paths = [config["freesurfer_ses-1"][roi] for roi in order]
            from nilearn.plotting import plot_roi, plot_prob_atlas

            plot_prob_atlas(
                roi_paths,
                title="ROIs",
                bg_img=config["spmCoregAnat"],
                display_mode="ortho",
                draw_cross=False,
                annotate=False,
                axes=ax,
            )
        else:
            palette_dict = None
            #   if datapath == "facetGrid_plot":
            # df["roi"] = pd.Categorical(
            #    df["roi"], categories=order, ordered=True, observed=True
            # )
            checkme = df.copy()

            if "hipp" in str(dataPath):
                filter_values = ["R. hipp.", "L. hipp."]
                df = df[~df.isin(filter_values).any(axis=1)]
                order = [x for x in order if x not in filter_values]
            elif "wb" in str(dataPath):
                filter_values = ["R. entorh. c.", "L. entorh. c."]
                df = df[~df.isin(filter_values).any(axis=1)]
                order = [x for x in order if x not in filter_values]
            df = df.assign(
                roi=pd.Categorical(df["roi"], categories=order, ordered=True)
            )
            df = df.sort_values("roi")

            blues = sns.color_palette("Blues", as_cmap=True)
            greens = sns.color_palette("Greens", as_cmap=True)
            rest_pallete = sns.color_palette("Oranges", as_cmap=True)

            # Get limited palettes
            limited_blues = get_limited_palette(blues, start=-1, end=6, n_colors=6)
            limited_greens = get_limited_palette(greens, start=-1, end=6, n_colors=6)
            limited_rest = get_limited_palette(
                rest_pallete, start=-1, end=6, n_colors=7
            )

            if "hipp" in str(dataPath):
                palette = limited_blues[:5] + limited_greens[:5] + limited_rest[:5]
                palette = [
                    (230 / 255, 159 / 255, 0 / 255),  # Orange
                    (86 / 255, 180 / 255, 233 / 255),  # Sky Blue
                    (0 / 255, 158 / 255, 115 / 255),  # Bluish Green
                    (213 / 255, 94 / 255, 0 / 255),  # Vermilion (Red-Orange)
                    (230 / 255, 159 / 255, 0 / 255),  # Orange
                    (86 / 255, 180 / 255, 233 / 255),  # Sky Blue
                    (0 / 255, 158 / 255, 115 / 255),  # Bluish Green
                    (213 / 255, 94 / 255, 0 / 255),  # Vermilion (Red-Orange)
                ]

                if "R. entorh. c." in order:
                    limited_rest = [limited_rest[i] for i in [1, 3, 4]]
                    palette = limited_blues + limited_greens + limited_rest
                    palette = [
                        (230 / 255, 159 / 255, 0 / 255),  # Orange
                        (86 / 255, 180 / 255, 233 / 255),  # Sky Blue
                        (0 / 255, 158 / 255, 115 / 255),  # Bluish Green
                        (213 / 255, 94 / 255, 0 / 255),  # Vermilion (Red-Orange)
                        (240 / 255, 228 / 255, 66 / 255),  # Yellow
                        (230 / 255, 159 / 255, 0 / 255),  # Orange
                        (86 / 255, 180 / 255, 233 / 255),  # Sky Blue
                        (0 / 255, 158 / 255, 115 / 255),  # Bluish Green
                        (213 / 255, 94 / 255, 0 / 255),  # Vermilion (Red-Orange)
                        (240 / 255, 228 / 255, 66 / 255),  # Yellow
                    ]

            if "wb" in str(dataPath):
                palette = (
                    limited_blues[-1:] + limited_greens[-1:] + limited_rest
                )  # vorher -3:-1
                palette = [
                    (0.4, 0.5, 0.8),  # right ph
                    (0.7, 0.3, 0.0),  # left ph
                    (0.6, 0.6, 0.6),  # ITS
                    (0.8, 0.2, 0.1),  # MTG
                    (0.0, 0.45, 0.7),  # ITG
                    (0.3, 0.7, 0.9),  # CUN
                    (0.8, 0.6, 0.7),  # PRECUN
                    (0.9, 0.6, 0.0),  # ACC
                    (0.0, 0.6, 0.5),
                ]

                if "R. entorh. c." in order:
                    palette = (
                        limited_blues[-3:-1] + limited_greens[-3:-1] + limited_rest
                    )
                    palette = [
                        (0.4, 0.5, 0.8),  # right ph
                        (0.7, 0.3, 0.0),  # left ph
                        (0.6, 0.2, 0.0),  # left ec
                        (0.7, 0.5, 0.6),  # right ec
                        (0.6, 0.6, 0.6),  # ITS
                        (0.8, 0.2, 0.1),  # MTG
                        (0.0, 0.45, 0.7),  # ITG
                        (0.3, 0.7, 0.9),  # CUN
                        (0.8, 0.6, 0.7),  # PRECUN
                        (0.9, 0.6, 0.0),  # ACC
                        (0.0, 0.6, 0.5),  # FS
                    ]

            palette_dict = {category: color for category, color in zip(order, palette)}
            import matplotlib.patches as mpatches

            handles = [
                mpatches.Patch(color=palette_dict[category], label=category)
                for category in order
            ]
            non_cat_df = df.copy()
            non_cat_df["roi"] = non_cat_df["roi"].astype(str)
            if plot_outliers == True:
                sns.boxplot(
                    y="roi",
                    x="simVal",
                    data=non_cat_df,
                    ax=ax,
                    fliersize=1.5,
                    order=order,
                    palette=palette_dict,
                )  # $showfliers=False)
            else:
                sns.boxplot(
                    y="roi", x="simVal", data=df, ax=ax, showfliers=False, order=order
                )

            # plt.legend(handles=handles, title="ROIs", bbox_to_anchor=(1.05, 1), loc="upper left")
            ax.axvline(0, color="r", linestyle="--")
            # ax.set_xticklabels(ax.get_xticklabels(), rotation=75, ha='right')
            ax.set_xticks(
                ax.get_xticks(),
                ax.get_xticklabels(),
                rotation=0,
                ha="center",
                rotation_mode="anchor",
            )

            # ax.set_xlabel("ROI",fontsize=12)
            ax.set_xlabel("Reinstatement (z-scored similarity)")
            ax.set_ylabel("Regions of interest")
            if not datapath == "facetGrid_plot":
                ax.set_title(f"{title}", size=16)
            ymin, ymax = ax.get_xlim()
            position = ymax - (ymax - ymin) * 0.06
            plt.show()

            # plot significance stars
            rois = df["roi"].unique()
            t_tests = df.groupby("roi", observed=True)["simVal"].apply(
                lambda x: ttest_1samp(x.dropna(), 0)
            )
            t_tests_dict = t_tests.to_dict()
            t_tests_list = [result.pvalue for result in t_tests]
            _, pvals_corrected = fdrcorrection(
                t_tests_list, alpha=0.05, method="indep", is_sorted=False
            )
            pvals_corrected_dict = dict(zip(t_tests_dict.keys(), pvals_corrected))
            sim_tests_df = pd.DataFrame(
                zip(
                    list(t_tests_dict.keys()),
                    list(t_tests_dict.values()),
                    pvals_corrected,
                ),
                columns=["roi", "test", "p"],
            )

            bf_list = []
            for i, roi in enumerate(rois):
                # Calculate Bayes Factor
                bf_list.append(
                    bayesfactor_ttest(
                        t_tests_dict[roi].statistic, t_tests_dict[roi].df
                    )
                )
                if t_tests_dict[roi].pvalue < 0.05:
                    if pvals_corrected_dict[roi] < 0.05:
                        ax.annotate(
                            "**", (position, i), va="center", fontsize=16, weight="bold"
                        )
                    else:
                        ax.annotate(
                            "*", (position, i), va="center", fontsize=16, weight="bold"
                        )
                        
            #add bf to output file
            sim_tests_df["bf"] = bf_list

            # Add legend to the plot
            if datapath == "facetGrid_plot":
                current_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
                datapath = datapath + current_time
                plt.legend(
                    title="ROIs",
                    handles=handles,
                    bbox_to_anchor=(1.25, 1),
                    loc="upper left",
                )
            else:
                outname = OutputDir / f"{datapath}.png"
                sim_tests_df.to_csv(
                    OutputDir / f"{datapath}_sim_t_tests.csv", index=False
                )
                plt.savefig(outname, bbox_inches="tight", dpi=333)

            if any(pattern in (datapath + ".png") for pattern in paper_fig_names):
                alt_OutputDir = (
                    os.environ["HOME"] + f"/TOAM/paper/src/figures/{datapath}.png"
                )
                plt.savefig(alt_OutputDir, bbox_inches="tight", dpi=333)

            plt.close()
            plt.cla()
            plt.clf()

    if (
        corrFile is not None
        and not isinstance(corrFile, list)
        and plot_correlations is True
    ):
        dfCorrTest = df.dropna()
        # corrFile["Acc"] = zscore(corrFile["Acc"])
        dfCorrTest = dfCorrTest.merge(corrFile, on="subject", how="outer")
        dfCorrTest = dfCorrTest.dropna().groupby("roi", observed=True)

        t_tests_correlation = dfCorrTest.apply(
            lambda x: list(pearsonr(x["simVal"], x["Acc"])), include_groups=False
        )
        corr_t_test_dict = t_tests_correlation.to_dict()
        corr_t_tests_list = [result[1] for result in corr_t_test_dict.values()]
        _, corr_p_vals = fdrcorrection(
            corr_t_tests_list, alpha=0.05, method="indep", is_sorted=False
        )
        corr_p_vals_dict = dict(zip(corr_t_test_dict.keys(), corr_p_vals))

        for roi in df["roi"].unique():
            dfroi = df[df["roi"] == roi]
            dfAll = dfroi.merge(corrFile, on="subject", how="outer")
            dfAll = dfAll.dropna()
            dfAll["AccZ"] = zscore(dfAll["Acc"])
            dfAll.drop(np.where(np.abs(dfAll["AccZ"]) > 3)[0], inplace=True)
            
            bf_pearson = pairwise_corr(
                    dfAll, ['simVal', 'Acc'], method='pearson',
                )

            x = dfAll["simVal"]
            y = dfAll["AccZ"]

            # save table
            file_path = OutputDir / f"{datapath}_correlationResultsTable.csv"
            # Überprüfen Sie, ob die Datei existiert
            if not os.path.exists(file_path):
                # Wenn die Datei nicht existiert, erstellen Sie sie
                with open(file_path, "w") as f:
                    f.write(f"roi,corrected_p_value,r-value,raw_p_value,bayes_pearson\n")
                    f.write(
                        f"{roi},{corr_p_vals_dict[roi]},{corr_t_test_dict[roi][0]:.2f},{corr_t_test_dict[roi][1],bf_pearson}\n"
                    )
            else:
                # Wenn die Datei existiert, öffnen Sie sie und fügen Sie eine neue Zeile hinzu
                with open(file_path, "a") as f:
                    f.write(
                        f"{roi},{corr_p_vals_dict[roi]},{corr_t_test_dict[roi][0]:.2f},{corr_t_test_dict[roi][1],bf_pearson}\n"
                    )

            if corr_t_test_dict[roi][1] < 0.05:
                viridis = mpl.colormaps["viridis"]
                if "unconsc_corr" in str(datapath):
                    colorchoice = 0
                    myylabel = "Guessing accuracy (z-scored)"
                    if "30min – 24h retrieval similarity for" in str(title):
                        myylabel = "Guessing accuracy delta (z-scored)"
                        colorchoice = 0
                else:
                    colorchoice = 0.5
                    myylabel = "Sure accuracy (z-scored)"
                    if "30min – 24h retrieval similarity for" in str(title):
                        myylabel = "Sure accuracy delta (z-scored)"
                        colorchoice = 0.5

                max_color = viridis(colorchoice)
                fig, ax = plt.subplots(figsize=(6, 6))
                sns.regplot(
                    x=x,
                    y=y,
                    fit_reg=True,
                    ax=ax,
                    line_kws={
                        "color": max_color,
                    },  # Line color
                    scatter_kws={"color": max_color, "s": 110, "alpha": 0.5},
                )
                ax.set_xlabel(f"{roi} (z-scored similarity)", size=16)
                ax.set_ylabel(myylabel, size=16)
                ax.set_title(f"{title}", size=16)

                ax.annotate(
                    f"R= {corr_t_test_dict[roi][0]:.2f}, p= {corr_t_test_dict[roi][1]:.3f}",  # , p_corrected= {corr_p_vals_dict[roi]:.2f}",
                    xy=(0.05, 0.9),
                    xycoords="axes fraction",
                    fontsize=12,
                )

                if datapath.endswith("unconsc_incorr"):
                    outname1 = (
                        os.environ["HOME"]
                        + f"/TOAM/paper/src/figures/{datapath}_{roi}_ERS_correl.png"
                    )
                else:
                    outname1 = OutputDir / f"{datapath}_{roi}_ERS_corr.png"
                if correct_results[0] in str(datapath) or correct_results[1] in str(datapath):
                    fig.savefig(outname1, bbox_inches="tight", dpi=333)
                plt.cla()
                plt.clf()
                plt.close()

    # Correction for multiple testing
    # FDR, Benjamini/Hochberg
    # https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html
    if not skip_rest:
        t_tests_list = []
        for tTest in t_tests:
            t_tests_list.append(tTest.pvalue)

        corrected_t_tests = fdrcorrection(
            t_tests_list, alpha=0.05, method="indep", is_sorted=False
        )
        # print(f"\nresults for {datapath}")
        for i in range(0, len(corrected_t_tests[0])):
            if (
                corrected_t_tests[1][i] < 0.05 and not printnotagain
            ):  # or t_tests.index[i] == 'rh_hipp':
                ind_list = t_tests.index.to_list()
                print(
                    f"sig: ROIS: {ind_list[i]}"  # orig pvalue: {t_tests.iloc[i].pvalue}, corrected p-value: {round(corrected_t_tests[1][i],4)}, statistic: {t_tests.iloc[i].statistic}"
                )

    return


def plot_multiple_ers(*args, **kwargs):
    ax = plot_ers_results(*args, **kwargs)
    print(ax)
    return ax


def perform_test_on_grouped(group):
    """
    performs a smolgorov-smirnoff on the group
    needs to have two distributions you're comparing
    """
    return list(ks_2samp(group["unconscious_correct"], group["unconscious_incorrect"]))


def plot_consc_unconsc_and_double_contrast(
    datapath,
    exp_id,
    corrFile,
    consc_level,
    data=None,
    remove_outliers=False,
    plot_outliers=True,
    plot_correlations=False,
):
    """
    plots the consc-separated similarity results and performs a significance test
    """
    # load data
    assert str(datapath).endswith("consc")
    if isinstance(datapath, list):
        df = pd.concat([pd.read_csv(str(f) + ".txt", sep="\t") for f in datapath])
        datapath = str(datapath[0]) + "mixed"
    elif data is not None:
        df = data
    else:
        df = pd.read_csv(str(datapath) + ".txt", sep="\t")  # ,header=None)

    if df["subject"].isin(["sub-60100"]).any():
        study = "wb"
    elif df["subject"].isin(["sub-60601"]).any():
        study = "hipp"
    else:
        raise ValueError("Unknown study")
    csv_OutputDir = Path(
        f"/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/{study}/bids/derivatives/nilearn/ERS/results/{exp_id}"
    )

    # numcols = round(len(df["roi"].unique()))
    if consc_level == "consc":
        delta = "consc-unconsc_incorr"
        minuend = "conscious_correct"
    elif consc_level == "unconsc":
        delta = "unconsc_corr-unconsc_incorr"
        minuend = "unconscious_correct"

    dfCopy = df.copy()
    dfCopy.loc[:, delta] = dfCopy[minuend] - dfCopy["unconscious_incorrect"]

    # do stat test on delta
    grouped = dfCopy.groupby("roi", observed=True)
    unconsc_diff_results = grouped.apply(perform_test_on_grouped, include_groups=False)
    out_file_path = csv_OutputDir / f"{datapath.split('/')[-1]}_{delta}.csv"
    # unconsc_diff_results.to_csv(out_file_path)

    df = dfCopy.loc[:, ["subject", "roi", delta]]
    df = df.rename(columns={delta: "simVal"})

    outstring = delta
    # do plotting on delta
    columnName = datapath + f"_{outstring}"
    plot_ers_results(
        datapath=columnName,
        exp_id=exp_id,
        data=df,
        remove_outliers=False,
        corrFile=corrFile,
        plot_outliers=True,
        plot_correlations=True,
    )


def convert_pvalue_to_asterisks(pvalue):
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
    return "ns"


### %% loop for actual plotting
experiment_id = time.strftime("%Y%m%d_%H%M%S")
sequences = ["hipp", "wb"]
for seq in sequences:

    fmri = Path(
        "/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/"
    )
    deriv = fmri / seq / "bids" / "derivatives"
    df = pd.read_csv(deriv / f"behav/df_{seq}_raw.csv")
    accs = pd.read_csv(fmri / "behav" / "accuracies_per_sub_retall_maxTrials.csv")
    behavData = fmri / "behav" / "accuracies_per_sub_ret_conscall_maxTrials.csv"
    behavData = pd.read_csv(behavData)
    if seq == "hipp":
        config = setup_configNL(60601)
    elif seq == "wb":
        config = setup_configNL(60100)
    else:
        raise ValueError("Unknown scanning sequence")
    accsHipp = accs[
        (accs["id"].isin(config["allSubsNum"])) & (accs["retrieval"] == 2)
    ].reset_index()
    # accsHipp = accsHipp[& accs['retrieval']=="2
    accsHipp["subject"] = accsHipp["id"].apply(lambda x: "sub-{:03d}".format(x))
    accsHipp = accsHipp[["subject", "Acc"]]

    dataPath = deriv / "nilearn" / "ERS" / "ers_data"

    config["allSubsStr"]
    numConscAnswers = (
        behavData[
            (behavData["id"].isin(config["allSubsNum"]))
            & (behavData["retrieval"] == 2)
            & (behavData["consc"] == 4)
        ]
        .loc[:, ["id", "n"]]
        .reset_index(drop=True)
        .rename(columns={"id": "subject", "n": "Acc"})
        .groupby("subject")
        .mean()
        .reset_index()
    )
    unconscAccuracy = (
        behavData[
            (behavData["id"].isin(config["allSubsNum"]))
            & (behavData["retrieval"] == 2)
            & (behavData["consc"] == 2)
        ]
        .loc[:, ["id", "Acc"]]
        .reset_index(drop=True)
        .rename(columns={"id": "subject", "n": "Acc"})
        .groupby("subject")
        .mean()
        .reset_index()
    )
    conscAccuracy = (
        behavData[
            (behavData["id"].isin(config["allSubsNum"])) & (behavData["retrieval"] == 2)
        ]
        .loc[:, ["id", "Acc"]]
        .reset_index(drop=True)
        .rename(columns={"id": "subject", "n": "Acc"})
        .groupby("subject")
        .mean()
        .reset_index()
    )

    numConscAnswersRet2 = (
        behavData[
            (behavData["id"].isin(config["allSubsNum"]))
            & (behavData["retrieval"] == 3)
            & (behavData["consc"] == 4)
        ]
        .loc[:, ["id", "n"]]
        .reset_index(drop=True)
        .rename(columns={"id": "subject", "n": "Acc"})
        .groupby("subject")
        .mean()
        .reset_index()
    )
    unconscAccuracyRet2 = (
        behavData[
            (behavData["id"].isin(config["allSubsNum"]))
            & (behavData["retrieval"] == 3)
            & (behavData["consc"] == 2)
        ]
        .loc[:, ["id", "Acc"]]
        .reset_index(drop=True)
        .rename(columns={"id": "subject", "n": "Acc"})
        .groupby("subject")
        .mean()
        .reset_index()
    )
    conscAccuracyRet2 = (
        behavData[
            (behavData["id"].isin(config["allSubsNum"])) & (behavData["retrieval"] == 3)
        ]
        .loc[:, ["id", "Acc"]]
        .reset_index(drop=True)
        .rename(columns={"id": "subject", "n": "Acc"})
        .groupby("subject")
        .mean()
        .reset_index()
    )

    conscAccuracyRet3 = (
        behavData[
            (behavData["id"].isin(config["allSubsNum"])) & (behavData["retrieval"] == 4)
        ]
        .loc[:, ["id", "Acc"]]
        .reset_index(drop=True)
        .rename(columns={"id": "subject", "n": "Acc"})
        .groupby("subject")
        .mean()
        .reset_index()
    )
    numConscAnswersRet3 = (
        behavData[
            (behavData["id"].isin(config["allSubsNum"]))
            & (behavData["retrieval"] == 4)
            & (behavData["consc"] == 4)
        ]
        .loc[:, ["id", "n"]]
        .reset_index(drop=True)
        .rename(columns={"id": "subject", "n": "Acc"})
        .groupby("subject")
        .mean()
        .reset_index()
    )
    unconscAccuracyRet3 = (
        behavData[
            (behavData["id"].isin(config["allSubsNum"]))
            & (behavData["retrieval"] == 4)
            & (behavData["consc"] == 2)
        ]
        .loc[:, ["id", "Acc"]]
        .reset_index(drop=True)
        .rename(columns={"id": "subject", "n": "Acc"})
        .groupby("subject")
        .mean()
        .reset_index()
    )

    numConscAnswers["subject"] = numConscAnswers["subject"].apply(
        lambda x: "sub-{:03d}".format(x)
    )
    unconscAccuracy["subject"] = unconscAccuracy["subject"].apply(
        lambda x: "sub-{:03d}".format(x)
    )
    conscAccuracy["subject"] = conscAccuracy["subject"].apply(
        lambda x: "sub-{:03d}".format(x)
    )
    numConscAnswersRet2["subject"] = numConscAnswersRet2["subject"].apply(
        lambda x: "sub-{:03d}".format(x)
    )
    unconscAccuracyRet2["subject"] = unconscAccuracyRet2["subject"].apply(
        lambda x: "sub-{:03d}".format(x)
    )
    conscAccuracyRet2["subject"] = conscAccuracyRet2["subject"].apply(
        lambda x: "sub-{:03d}".format(x)
    )
    numConscAnswersRet3["subject"] = numConscAnswersRet3["subject"].apply(
        lambda x: "sub-{:03d}".format(x)
    )
    unconscAccuracyRet3["subject"] = unconscAccuracyRet3["subject"].apply(
        lambda x: "sub-{:03d}".format(x)
    )
    conscAccuracyRet3["subject"] = conscAccuracyRet3["subject"].apply(
        lambda x: "sub-{:03d}".format(x)
    )
    ConscAnswersDelta = numConscAnswers.merge(
        numConscAnswersRet2, on="subject", suffixes=("_ret1", "_ret2")
    )
    ConscAnswersDelta["Acc"] = (
        ConscAnswersDelta["Acc_ret2"] - ConscAnswersDelta["Acc_ret1"]
    )
    ConscAnswersDelta.drop(columns=["Acc_ret1", "Acc_ret2"], inplace=True)

    unconscAccDelta = unconscAccuracy.merge(
        unconscAccuracyRet2, on="subject", suffixes=("_ret1", "_ret2")
    )
    unconscAccDelta["Acc"] = unconscAccDelta["Acc_ret2"] - unconscAccDelta["Acc_ret1"]
    unconscAccDelta.drop(columns=["Acc_ret1", "Acc_ret2"], inplace=True)

    accDelta = conscAccuracy.merge(
        conscAccuracyRet2, on="subject", suffixes=("_ret1", "_ret2")
    )
    accDelta["Acc"] = accDelta["Acc_ret2"] - accDelta["Acc_ret1"]
    accDelta.drop(columns=["Acc_ret1", "Acc_ret2"], inplace=True)

    ## %% hipp
    if seq == "hipp":
        dataList = [
            "20240530_hipp-array_n_enc_ret1_perm",
            "20240530_hipp-array_n_enc_ret1_perm_consc",
            "20240530_hipp-array_n_enc_ret2_perm",
            "20240530_hipp-array_n_enc_ret2_perm_consc",
            "20240530_hipp-array_n_ret1_ret2_perm",
            "20240530_hipp-array_n_ret1_ret2_perm_consc",
            "20240710_hipp-all_memory_n_enc_recog_perm",
            "20240710_hipp-all_memory_n_enc_recog_perm_consc",
        ]
    elif seq == "wb":
        dataList = [
            "20240530_wb-array_n__enc_ret1_perm",
            "20240530_wb-array_n__enc_ret1_perm_consc",
            "20240530_wb-array_n__enc_ret2_perm",
            "20240530_wb-array_n__enc_ret2_perm_consc",
            "20240530_wb-array_n__ret1_ret2_perm",
            "20240530_wb-array_n__ret1_ret2_perm_consc",
            "20240710_wb-all_memory_n_enc_recog_perm",
            "20240710_wb-all_memory_n_enc_recog_perm_consc",
        ]
    else:
        raise ValueError("Unknown scanning sequence")

    corrList = [
        conscAccuracy,
        [numConscAnswers, conscAccuracy, unconscAccuracy],
        conscAccuracyRet2,
        [numConscAnswersRet2, conscAccuracyRet2, unconscAccuracyRet2],
        accDelta,
        [ConscAnswersDelta, accDelta, unconscAccDelta],
        conscAccuracyRet3,  # recog
        [numConscAnswersRet3, conscAccuracyRet3, unconscAccuracyRet3],  # recog
    ]

    for i, (my_dataset, corr_file) in enumerate(zip(dataList, corrList)):
        ers_data = dataPath / my_dataset
        df = pd.read_csv(dataPath / f"{my_dataset}.txt", sep="\t")
        duplicates = df.duplicated()
        # print(sum(duplicates))
        print(f"working on dataset: {my_dataset}")
        plot_ers_results(
            datapath=str(ers_data),
            exp_id=experiment_id,
            remove_outliers=False,
            corrFile=corr_file,
            plot_correlations=True,
        )
        if "consc" in str(ers_data):
            consclevels = ["consc", "unconsc"]
            corrfiles = [corr_file[0], corr_file[-1]]
            for consclevel, correct_corrfile in zip(consclevels, corrfiles):
                plot_consc_unconsc_and_double_contrast(
                    str(ers_data),
                    exp_id=experiment_id,
                    consc_level=consclevel,
                    remove_outliers=False,
                    corrFile=correct_corrfile,  # just take the unconscious accuracy
                    plot_correlations=True,
                )

    print(f"putting output to experiment folder: {experiment_id}")

# %%
