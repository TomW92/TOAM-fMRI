# load libraries
import sys

sys.path.append("/storage/homefs/tw18a205/TOAM")
import numpy as np
from nilearn.masking import apply_mask
from nilearn.image import resample_to_img, resample_img, load_img
from nibabel.nifti1 import load
from nibabel.nifti1 import Nifti1Image
import time

# from sys import platform
from fmri.mvpa.mvpa_utils.get_datasets import *
from scipy import stats
from nilearn.plotting import plot_stat_map, plot_img_on_surf, plot_anat, view_img
from sklearn.metrics import pairwise_distances
import os
import sys
import sklearn as skl
from pathlib import Path
import warnings
from joblib import Parallel, delayed


def pairwise_rsa(
    subject,
    experiment_id,
    roi,
    run1,
    run2,
    force=True,
    perm="doit",
):
    """
    Get the data for a given subject and ROI and create an RSA dataset
    Makes use of the function create_mv_data

    subject = 'int for subject number, or string for subject name'
    roi = 'string for ROI name' valid inputs: everythingfound under config['freesurfer_ses-1'] and config['freesurfer_ses-2']
    outputfilename = 'string for the output filename'
    run2 = 'string for the retrieval session' valid inputs: 'ret1' or 'ret2'
    """
    warnings.resetwarnings()
    print("experiment_id: ", experiment_id)
    # get subject info
    config = setup_configNL(subject)
    print(f"running RSA for {roi} {run1} {run2}")
    outputdir = Path(
        f'/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/{config["study"]}/bids/derivatives/nilearn/ERS/'
    )
    # save Outputs in nilearn Path
    rsa_data_outputfilename = f"{experiment_id}_{run1}_{run2}"
    rsa_outputdir = outputdir / "ers_data"
    if not os.path.exists(rsa_outputdir):
        os.makedirs(rsa_outputdir)

    encTmaps = []
    enc_labels = []
    if run1 == "enc":
        enc_labels, enc_betas, _ = create_mv_data(
            "Encoding", "face_id", config, pf="s1r"
        )  #'s1r'=smoothed with 1mm kernel
        # enc_labels, enc_betas, _ = get_lss_tMaps('enc',config)
        encTmaps = enc_betas.replace("beta", "spmT", regex=True)
    else: # run1 == "ret1":
        enc_labels, ret_betas, _ = create_mv_data(
            "Retrieval", "face_id", config, pf="s1r"
        )  #'s1r'=smoothed with 1mm kernel
        encTmaps = ret_betas.replace("beta", "spmT", regex=True)

    retTmaps = []
    ret_labels = []
    # get the retrieval data, if retrieval 2, register, if not just resample.
    if (run2 == "ret2") or (run2 == "recog"):
        if run2 == "ret2":
            ret_labels, ret_betas_ses2Space, _ = create_mv_data(
                "DelayedRetrieval", "face_id", config, pf="s1r"
            )
            accLabels, _, _ = create_mv_data(
                "DelayedRetrieval", "accuracy", config, pf="s1r"
            )  #'s1r'=smoothed with 1mm kernel
            conscLabels, _, _ = create_mv_data(
                "DelayedRetrieval", "consciousAll", config, pf="s1r"
            )  #'s1r'=smoothed with 1mm kernel
        elif run2 == "recog":
            ret_labels, ret_betas_ses2Space, _ = create_mv_data(
                "Recognition", "face_id", config, pf="s1r"
            )
            accLabels, _, _ = create_mv_data(
                "Recognition", "accuracy", config, pf="s1r"
            )
            conscLabels, _, _ = create_mv_data(
                "Recognition", "consciousAll", config, pf="s1r"
            )

        # replace beta with tmap, as done in hebscher 2021
        retTmaps_Ses2Space = ret_betas_ses2Space.replace("beta", "spmT", regex=True)
        RSes2 = np.loadtxt(config["affineMatrix-ses2ToT1w"], delimiter=" ")
        for tMap in retTmaps_Ses2Space:
            tMap = load(tMap)
            tMap_t1w_space = Nifti1Image(
                tMap.get_fdata(), np.linalg.inv(RSes2) @ tMap.affine
            )
            retTmaps.append(tMap_t1w_space)

    else: # run2 == "ret1":
        ret_labels, ret_betas, _ = create_mv_data(
            "Retrieval", "face_id", config, pf="s1r"
        )  #'s1r'=smoothed with 1mm kernel
        # ret_labels, ret_betas, _ = get_lss_tMaps('ret1',config)
        accLabels, _, _ = create_mv_data(
            "Retrieval", "accuracy", config, pf="s1r"
        )  #'s1r'=smoothed with 1mm kernel
        conscLabels, _, _ = create_mv_data(
            "Retrieval", "consciousAll", config, pf="s1r"
        )  #'s1r'=smoothed with 1mm kernel
        retTmaps = ret_betas.replace("beta", "spmT", regex=True)

    # get mask, resample to image, and then extract data
    mask = config["freesurfer_ses-1"][
        roi
    ]  # use session one because effectively, the T1w is the ses-1 space. the unbiased template was created based on ses-1 space.
    # anat = config['spmCoregAnat']
    biasCorrAnat = load(config["spmBiasCorr"])

    # resample the already registered images to the anatomical sample that is in the same
    resampledFolder = rsa_outputdir / f"resampledNiftis"
    os.makedirs(resampledFolder, exist_ok=True)
    resampled_1 = os.path.join(resampledFolder, f"{subject}_{roi}_{run1}.npy")
    resampled_2 = os.path.join(resampledFolder, f"{subject}_{roi}_{run2}.npy")

    if os.path.exists(resampled_1) and not force:
        maskedDataEnc = np.load(resampled_1, allow_pickle=True)
        print(f"loaded resampled data: {run1}")
    else:
        start = time.time()
        # resampledEncImages = [resample_to_img(image_path, biasCorrAnat, interpolation='continuous',copy=False) for image_path in encTmaps]
        resampledEncImages = Parallel(n_jobs=-1, verbose=1, backend="multiprocessing")(
           delayed(resample_to_img)(
               image_path, biasCorrAnat, interpolation="continuous", copy=False
           )
           for image_path in encTmaps
        )
        stop = time.time()
        print(f"run1 resampling done, it took {stop-start} seconds")
        # mask Enc and RET Data
        start = time.time()
        rsmpldMask = resample_to_img(
            mask, biasCorrAnat, interpolation="nearest"
        )  # mask
        # maskedDataEnc = apply_mask(resampledEncImages, rsmpldMask)
        maskedDataEnc = Parallel(n_jobs=-1, verbose=1, backend="multiprocessing")(
           delayed(apply_mask)(image, rsmpldMask) for image in resampledEncImages
        )

        np.save(resampled_1, maskedDataEnc)
        stop = time.time()
        print(
            f"Masking (extracting ROI Data -> Nilearn Masking) done, it took {stop-start} seconds"
        )
    if os.path.exists(resampled_2) and not force:
        maskedDataRet = np.load(resampled_2, allow_pickle=True)
        print(f"loaded resampled data: {run2}")
    else:
        start = time.time()
        resampledRetImages = Parallel(n_jobs=-1, verbose=1, backend="multiprocessing")(
            delayed(resample_to_img)(
                image_path, biasCorrAnat, interpolation="continuous", copy=False
            )
            for image_path in retTmaps
        )
        # resampledRetImages = [resample_to_img(image_path, biasCorrAnat, interpolation='continuous',copy=False) for image_path in retTmaps]
        rsmpldMask = resample_to_img(
            mask, biasCorrAnat, interpolation="nearest"
        )  # mask image is coregistered, but not in the same space, not the same dimensions. that is fixed here.
        stop = time.time()
        print(f"run 2 resampling done, it took {stop-start} seconds")
        
        # mask Enc and RET Data
        start = time.time()
        # maskedDataRet = apply_mask(resampledEncImages, rsmpldMask)
        maskedDataRet = Parallel(n_jobs=-1, verbose=1, backend="multiprocessing")(
            delayed(apply_mask)(image, rsmpldMask) for image in resampledRetImages
        )
        np.save(resampled_2, maskedDataRet)
        stop = time.time()
        print(
            f"Masking (extracting ROI Data -> Nilearn Masking) done, it took {stop-start} seconds"
        )

    df1 = pd.DataFrame(
        list(zip(enc_labels, maskedDataEnc )),
        columns=["Label", roi + "EncData"],
    )
    df2 = pd.DataFrame(
        list(zip(ret_labels, maskedDataRet, accLabels, conscLabels)),
        columns=["Label", roi + "RetData", "Accuracy", "conscLabels"],
    )

    # merge the dataframes
    df = df1.merge(
        df2, on="Label", how="inner"
    )  # changed from outer to inner, because of LSS computations, which are independent
    df["conscLabels"] = df["conscLabels"].astype(str) + "_" + df["Accuracy"].astype(str)
    # calculrate Encoding - Retrieval Similarity; do it pairwise for encoding retrieval pairs that are on the same row in deffernt columns
    similarities = pairwise_distances(
        np.stack(np.array(df[roi + "EncData"])),
        np.stack(np.array(df[roi + "RetData"])),
        metric="correlation",
    )
    print("max simvalue: ", np.max(similarities))
    print("min simvalue: ", np.min(similarities))

    # fisher-z-transform the correlation metrix
    zTransSimilarities = np.arctanh(1 - similarities)

    # calculate the average similarity for pairs and non-pairs (that is, diagonal and the off-diagonal)
    arrMask = np.ones(zTransSimilarities.shape, dtype=bool)
    np.fill_diagonal(arrMask, 0)
    pairSimAverage = np.mean(np.diagonal(zTransSimilarities))
    nonPairSimAverage = np.mean(zTransSimilarities[arrMask])
    simVal = pairSimAverage - nonPairSimAverage

    # for zvalues
    filename = rsa_data_outputfilename + ".txt"

    # save the data
    file_path = rsa_outputdir / filename
    # Überprüfen Sie, ob die Datei existiert
    if not os.path.exists(file_path):
        # Wenn die Datei nicht existiert, erstellen Sie sie
        with open(file_path, "w") as f:
            f.write(f"subject\troi\tsimVal\n")
            f.write(f"{subject}\t{roi}\t{simVal}\n")
    else:
        # Wenn die Datei existiert, öffnen Sie sie und fügen Sie eine neue Zeile hinzu
        with open(file_path, "a") as f:
            f.write(f"{subject}\t{roi}\t{simVal}\n")

    print("doing permutation...")
    nPerm = 1000
    perm_stats = run_perms(df,roi,nPerm)

    # create a z-value as simVal based on the permutation analysis
    permMean = np.mean(perm_stats)
    print("permMean: ", permMean)
    permStdDev = np.std(perm_stats)
    perm_simVal = (simVal - permMean) / permStdDev

    # for zvalues
    filename = rsa_data_outputfilename + "_perm.txt"
    print("performed permutation test")

    # save the data
    file_path = rsa_outputdir / filename
    # Überprüfen Sie, ob die Datei existiert
    if not os.path.exists(file_path):
        # Wenn die Datei nicht existiert, erstellen Sie sie
        with open(file_path, "w") as f:
            f.write(f"subject\troi\tsimVal\n")
            f.write(f"{subject}\t{roi}\t{perm_simVal}\n")
    else:
        # Wenn die Datei existiert, öffnen Sie sie und fügen Sie eine neue Zeile hinzu
        with open(file_path, "a") as f:
            f.write(f"{subject}\t{roi}\t{perm_simVal}\n")

    # create simvals based on conscLabel
    labelsOfInterest = [
        "unconscious_correct",
        "unconscious_incorrect",
        "conscious_correct",
    ]
    dfConsc = df[df["conscLabels"].isin(labelsOfInterest)]
    conscLabels = dfConsc["conscLabels"].unique()
    sorted_indices = np.argsort(conscLabels)
    conscLabels = conscLabels[sorted_indices] #get the same and right order for every participants!
    conscSimVals = []
    conscPermSimVals = []
    for label in conscLabels:
        dfConscLabel = dfConsc[dfConsc["conscLabels"] == label]
        similarities = pairwise_distances(
            np.stack(np.array(dfConscLabel[roi + "EncData"])),
            np.stack(np.array(dfConscLabel[roi + "RetData"])),
            metric="correlation",
        )
        zTransSimilarities = np.arctanh(1 - similarities)
        arrMask = np.ones(zTransSimilarities.shape, dtype=bool)
        np.fill_diagonal(arrMask, 0)
        pairSimAverage = np.mean(np.diagonal(zTransSimilarities))
        nonPairSimAverage = np.mean(zTransSimilarities[arrMask])
        consc_simVal = pairSimAverage - nonPairSimAverage

        print(f"doing permutation for {label}...")
        nPerm = 1000
        perm_stats = run_perms(dfConscLabel,roi,nPerm)

        # create a z-value as simVal based on the permutation analysis
        permMean = np.mean(perm_stats)
        print(f"permMean: {permMean}")
        permStdDev = np.std(perm_stats)
        perm_simVal = (consc_simVal - permMean) / permStdDev

        conscPermSimVals.append(perm_simVal)
        conscSimVals.append(consc_simVal)

    filename = rsa_data_outputfilename + "_consc.txt"
    file_path = rsa_outputdir / filename
    # Überprüfen Sie, ob die Datei existiert
    # If the file does not exist, create it
    if not os.path.exists(file_path):
        with open(file_path, "w") as f:
            f.write("subject\troi\t" + "\t".join(conscLabels) + "\n")
            f.write(f"{subject}\t{roi}\t" + "\t".join(map(str, conscSimVals)) + "\n")
    else:
        # If the file exists, open it and add a new line
        with open(file_path, "a") as f:
            print(conscLabels)
            f.write(f"{subject}\t{roi}\t" + "\t".join(map(str, conscSimVals)) + "\n")

    filename = rsa_data_outputfilename + "_perm_consc.txt"
    file_path = rsa_outputdir / filename
    # Überprüfen Sie, ob die Datei existiert
    # If the file does not exist, create it
    if not os.path.exists(file_path):
        with open(file_path, "w") as f:
            f.write("subject\troi\t" + "\t".join(conscLabels) + "\n")
            f.write(
                f"{subject}\t{roi}\t" + "\t".join(map(str, conscPermSimVals)) + "\n"
            )
    else:
        # If the file exists, open it and add a new line
        with open(file_path, "a") as f:
            print(conscLabels)
            f.write(
                f"{subject}\t{roi}\t" + "\t".join(map(str, conscPermSimVals)) + "\n"
            )

    print("calculated RSA")
    return

def run_one_perm(df,roi,i):
    # shuffle the encoding data, that should be enough to create mismatched pairs.
    df.reset_index(drop=True, inplace=True)
    df[roi + "EncData"] = (
        df[roi + "EncData"].sample(frac=1).reset_index(drop=True)
    )
    # calculrate Encoding - Retrieval Similarity
    similarities = pairwise_distances(
        np.stack(np.array(df[roi + "EncData"])),
        np.stack(np.array(df[roi + "RetData"])),
        metric="correlation",
    )  # metric = 'correlation' was the default for now
    # fisher-z-transform the correlation metrix
    zTransSimilarities = np.arctanh(1 - similarities)

    # calculate the average similarity for pairs and non-pairs (that is, diagonal and off-diagonal)
    arrMask = np.ones(zTransSimilarities.shape, dtype=bool)
    np.fill_diagonal(arrMask, 0)
    pairSimAverage = np.mean(np.diagonal(zTransSimilarities))
    nonPairSimAverage = np.mean(zTransSimilarities[arrMask])
    permVal = pairSimAverage - nonPairSimAverage
    return permVal

def run_perms(df,roi,n_perms):
    permStats = Parallel(n_jobs=-1,verbose=2,backend='multiprocessing')(
            delayed(run_one_perm)(df,roi,i)
            for i in range(n_perms)
        )
    # for i in range(n_perms):
    #     permStats = []
    #     permStats.append(run_one_perm(df,roi,i))
    return permStats

def main():
    combinations = [
        ["enc",
        "recog"]
    ]  # [["enc", "ret1"], ["enc", "ret2"], ["ret1", "ret2"]]
    config = setup_configNL(60601)
    subjects = list(config["allSubsStr"])
    # subject = subjects[subject_id]
    rois = []
    filestring = []
    if config['study'] == 'hipp':
        filestring = "hipp"
        rois = [
            "cuneus_mask",
            "precuneus_mask",
            #"inf_temp_gyr_mask",
            "mid_temp_gyr_mask",
            #"inferior_temp_sulc_mask",
            #"parahipp_lh_mask",
            #"parahipp_rh_mask",
            "entorhinal_lh_mask",
            "entorhinal_rh_mask",
            "parahipp_lh_no_er_mask",
            "parahipp_rh_no_er_mask",
            "lh_hipp",
            "rh_hipp",
            "hipp_tail_lh",
            "hipp_tail_rh",
            "hipp_body_lh",
            "hipp_body_rh",
            "hipp_head_lh",
            "hipp_head_rh",
        ]
    elif config['study'] == 'wb':
        filestring = "wb"
        rois = [
            "cuneus_mask",
            "acc_mask",
            "front_sup_mask",
            "precuneus_mask",
            "inf_temp_gyr_mask",
            "mid_temp_gyr_mask",
            "inferior_temp_sulc_mask",
            # "parahipp_lh_mask",
            # "parahipp_rh_mask",
            "entorhinal_lh_mask",
            "entorhinal_rh_mask",
            "parahipp_lh_no_er_mask",
            "parahipp_rh_no_er_mask",
            "lh_hipp",
            "rh_hipp",
        ]
    else:
        ValueError("study not recognized")

    indexed_list = [(sub, r) for sub in subjects for r in rois]
    print(len(indexed_list))
    # print("rois: ", rois)
    # print("subbjects: ", subjects)
    exp_id = time.strftime("%Y%m%d_") + f"{filestring}-all_failed_memory_n"   
    array_ID = int(os.environ["SLURM_ARRAY_TASK_ID"])
    print("array_ID: ", array_ID)
    subject = indexed_list[array_ID][0]
    roi = indexed_list[array_ID][1]

    print("working on subject: ", subject)
    # Parallel(n_jobs=-1, verbose=2, max_nbytes=None, backend="multiprocessing")(
    #     delayed(pairwise_rsa)(subject, exp_id, r, c[0], c[1])  # , backend="multiprocessing"
    #     for r in roi
    for c in combinations:
        pairwise_rsa(subject, exp_id, roi, c[0], c[1])
        print("working on combinations", c)
        print(f"done with {subject} and {roi}")

def main_interactive(subject,id):
    combinations = [["enc", "recog"]]#[["enc", "ret1"], ["enc", "ret2"], ["ret1", "ret2"]]
    #c = combinations[0]
    #config = setup_configNL(60601)
    #subjects = config["allSubsStr"]
    #subject_id = int(os.environ["SLURM_ARRAY_TASK_ID"])
    #subject = subjects[subject_id]
    roi = [
        "cuneus_mask",
        "precuneus_mask",
        "inf_temp_gyr_mask",
        "mid_temp_gyr_mask",
        "inferior_temp_sulc_mask",
        "parahipp_lh_mask",
        "parahipp_rh_mask",
        "lh_hipp",
        "rh_hipp",
        "hipp_tail_lh",
        "hipp_tail_rh",
        "hipp_body_lh",
        "hipp_body_rh",
        "hipp_head_lh",
        "hipp_head_rh",
    ]
    exp_id = id#time.strftime("%Y%m%d_%H")
    print("working on subject: ", subject)
    #print("working on combinations", c)
    Parallel(n_jobs=-1, verbose=2, max_nbytes=None, backend="multiprocessing")(
        delayed(pairwise_rsa)(subject, exp_id, r, c[0], c[1])  # , backend="multiprocessing"
        for r in roi
        for c in combinations
    )


if __name__ == "__main__":
    main()
    # globals()[args[1]](*args[2:])
