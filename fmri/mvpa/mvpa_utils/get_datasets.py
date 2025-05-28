import numpy as np
import pandas as pd
import os
import fnmatch
import matplotlib.colors as colors
from sys import platform
import nipype.interfaces.spm as spm
import glob
from pathlib import Path
from scipy.io import loadmat


def RDMcolormapObject(direction=1):
    """
    Returns a matplotlib color map object for RSA and brain plotting
    """
    if direction == 0:
        cs = ["yellow", "red", "gray", "turquoise", "blue"]
    elif direction == 1:
        cs = ["blue", "turquoise", "gray", "red", "yellow"]
    else:
        raise ValueError("Direction needs to be 0 or 1")
    
    cmap = colors.LinearSegmentedColormap.from_list("", cs)
    return cmap


def readITKtransform(transform_file):
    """ """

    # read the transform
    transform = None
    with open(transform_file, "r") as f:
        for line in f:

            # check for Parameters:
            if line.startswith("Parameters:"):
                values = line.split(": ")[1].split(" ")

                # filter empty spaces and line breaks
                values = [float(e) for e in values if (e != "" and e != "\n")]
                # create the upper left of the matrix
                transform_upper_left = np.reshape(values[0:9], (3, 3))
                # grab the translation as well
                translation = values[9:]

            # check for FixedParameters:
            if line.startswith("FixedParameters:"):
                values = line.split(": ")[1].split(" ")

                # filter empty spaces and line breaks
                values = [float(e) for e in values if (e != "" and e != "\n")]
                # setup the center
                center = values

    # compute the offset
    offset = np.ones(4)
    for i in range(0, 3):
        offset[i] = translation[i] + center[i]
        for j in range(0, 3):
            offset[i] -= transform_upper_left[i][j] * center[i]

    # add the [0, 0, 0] line
    transform = np.vstack((transform_upper_left, [0, 0, 0]))
    # and the [offset, 1] column
    transform = np.hstack((transform, np.reshape(offset, (4, 1))))

    return transform


def infer_button_press(row, session):
    if (session == "Retrieval") or (session == "DelayedRetrieval"):
        if (row["ret"] == "correct" and row["obj_cat"] == "organic") or (
            row["ret"] == "incorrect" and row["obj_cat"] == "inorganic"
        ):
            return pd.Series(["middle"])
        else:
            return pd.Series(["index"])
    elif session == "ImmRecall":
        if (row["Recall"] == "correct" and row["obj_cat"] == "organic") or (
            row["Recall"] == "incorrect" and row["obj_cat"] == "inorganic"
        ):
            return pd.Series(["middle"])
        else:
            return pd.Series(["index"])


def create_model_RDM(label_vector):
    # input is label vector, with 0 and 1 for a certain category like remembered, not remembered
    # Author: Tom Willems
    num_labels = len(label_vector)
    model_rdm = np.ones((num_labels, num_labels))

    for row_ind in range(0, num_labels):
        for column_ind in range(0, num_labels):
            if label_vector[row_ind] == 1 & label_vector[column_ind] == 1:
                model_rdm[row_ind, column_ind] = 0
    return model_rdm


def normalize_SPM_segmented_data(files, session, config):
    if session == "day1":
        deformation_field = config["defSess1"]
    else:
        deformation_field = config["defSess2"]
    norm12 = spm.Normalize12()
    norm12.inputs.deformation_file = deformation_field  # deformation field for this session, to spatially normalize all files for this subject
    norm12.inputs.apply_to_files = files.to_list()  # betas_train.to_list()
    norm12.inputs.write_voxel_sizes = [1, 1, 1]
    norm12.inputs.jobtype = "write"
    norm12.run()  # run normalization


def get_lss_tMaps(session, config):
    """
    get_lss_tMaps - get the t-maps from the LSS model
    Syntax: [tMaps] = get_lss_tMaps(session,config)

    Inputs:
    session - the session you want to get the t-maps from, enc or ret1 at the moment
    config - the config file for this particular subject

    Outputs:
    tMaps - the t-maps for the LSS model
    """

    folder_names = []
    spmT_files = []
    tMapPath = Path(
        config["genOut"] + "/spmMB/" + config["subStr"] + f"/LSS_GLMs_all/{session}_LSS"
    )
    for folder in os.listdir(tMapPath):
        if os.path.isdir(os.path.join(tMapPath, folder)):
            folder_names.append(folder)
            for file in os.listdir(os.path.join(tMapPath, folder)):
                if "spmT_0001" in file:
                    spmT_files.append(os.path.join(tMapPath, folder, file))
                    break  # stop searching after finding the first 'spmTfile' in the folder
    labels = ["_".join(folderName.split("_")[5:7]) for folderName in folder_names]
    labels = pd.Series(labels)
    t_maps = pd.Series(spmT_files)
    sl_mask = config["STE_native_s1r"]["STE_enc_mask"]

    return labels, t_maps, sl_mask


def create_mv_data(session, category, config, space="native", pf="s1r"):
    """
    creates multivariate data for decoding/RSA analyses.
    input:
        session (enc, ret30min, ret24h, recognition),
        category (object,sex,consciousness_level, accuracy) i.e., what is to be decoded
        config: the config file for this particular subject.
        space: the space in which images are loaded (options: "mni", "native")
        prefix: either "s1r" or "r" (default). selector for data that is either smoothed with a 1mm gauss. kernel, or not.
    """
    labels = pd.Series()

    if (
        (session == "Encoding" and category == "unconscious_accuracy")
        or (session == "Encoding" and category == "unconscious_object_category")
        or (session == "Encoding" and category == "buttonpress")
    ):
        raise Exception(
            "Invalid combination of category and session, combination was: {} and {}".format(
                session, category
            )
        )

    if session == "ImmRecall":
        behav = pd.read_csv(config["eventsEncoding"], sep="\t")  # get encoding behav
        behav["item"] = pd.Categorical(
            behav["item"], categories=behav["item"].unique(), ordered=True
        )  # otherwise pd.pivot messes up the directions. depreceated?
        behav = pd.pivot(behav, index="item", columns="trial_type", values="accuracy")

        betas = config["STE_" + space + "_" + pf][
            "STE_" + session
        ]  # add betas from config from this subject
        spm = loadmat(
            betas[0].replace("beta_0001.nii", "SPM.mat")
        )  # get spm file that created the beta nifti's.
        beta_ident = spm["SPM"][0][0]["xX"]["name"][0][0][
            0
        ]  # in order to prevent mismatched labels
        beta_names = [
            str(item.flat[0][6:19]) for item in beta_ident[0 : len(betas)]
        ]  # by merging the betas
        beta_files_dict = {
            "item": beta_names,
            "betas": betas,
        }  # with the behavioral file
        beta_files = pd.DataFrame(
            beta_files_dict
        )  # on the item (filename of the presented face, i.e. fol_060253.png):
        df = behav.merge(beta_files, how="inner", on="item")
        df.reset_index(inplace=True)
        sl_mask = config["STE_" + space + "_" + pf]["STE_ImmRecall_mask"]

    elif session == "Encoding":
        behav = pd.read_csv(config["events" + session], sep="\t")  # get encoding behav
        behav = behav[
            behav["trial_type"] == "Enc"
        ]  # its in long format, only take the enc trials, so you have the same length
        betas = config["STE_" + space + "_" + pf][
            "STE_enc"
        ]  # add betas from config from this subject
        spm = loadmat(
            betas[0].replace("beta_0001.nii", "SPM.mat")
        )  # get spm file that created the beta nifti's.
        beta_ident = spm["SPM"][0][0]["xX"]["name"][0][0][
            0
        ]  # in order to prevent mismatched labels
        beta_names = [
            str(item.flat[0][6:19]) for item in beta_ident[0 : len(betas)]
        ]  # by merging the betas
        beta_files_dict = {
            "item": beta_names,
            "betas": betas,
        }  # with the behavioral file
        beta_files = pd.DataFrame(
            beta_files_dict
        )  # on the item (filename of the presented face, i.e. fol_060253.png):
        df = behav.merge(beta_files, how="inner", on="item")
        df.reset_index(inplace=True)

        # merge again with retrieval data to get info about encoding sussex.
        ret_behav = pd.read_csv(
            config["eventsRetrieval"], sep="\t"
        )  # get retrieval behav data for consc- level of enc
        ret_behav = ret_behav[
            ret_behav["trial_type"] == "retRatingFeedback"
        ]  # its in long format, only take the ret trials with consc level, so you have same length
        df = df.merge(
            ret_behav, how="inner", on="item", suffixes=("_enc", "_ret")
        )  # merge on trial item (filename of the presented face)
        df.reset_index(inplace=True)

        sl_mask = config["STE_" + space + "_" + pf]["STE_enc_mask"]

    elif session == "Retrieval":
        # just take the relevant trials, onsets, etc.
        # behav = behav[behav['trial_type'] == 'retRatingFeedback']  # reduce to the correct number of trials
        # instead of selecting by retRatingFeedback, now pivot wider to keep your options for different labels later on.
        # filtering for retRatingFeedback only keeps labels for consciousness ratings, and you would have to define more than 4 sessions, which is not logical
        behav = pd.read_csv(config["events" + session], sep="\t")
        behav["item"] = pd.Categorical(
            behav["item"], categories=behav["item"].unique(), ordered=True
        )  # otherwise pd.pivot messes up the directions. depreceated?
        behav = pd.pivot(behav, index="item", columns="trial_type", values="accuracy")

        betas = config["STE_" + space + "_" + pf]["STE_ret"]
        spm = loadmat(
            betas[0].replace("beta_0001.nii", "SPM.mat")
        )  # get spm file that created the beta nifti's.
        beta_ident = spm["SPM"][0][0]["xX"]["name"][0][0][
            0
        ]  # in order to prevent mismatched labels
        beta_names = [
            str(item.flat[0][6:19]) for item in beta_ident[0 : len(betas)]
        ]  # by merging the betas
        beta_files_dict = {
            "item": beta_names,
            "betas": betas,
        }  # with the behavioral file
        beta_files = pd.DataFrame(
            beta_files_dict
        )  # on the item (filename of the presented face, i.e. fol_060253.png):
        df = behav.merge(beta_files, how="inner", on="item")
        df.reset_index(inplace=True)

        sl_mask = config["STE_" + space + "_" + pf]["STE_ret_mask"]

    elif session == "DelayedRetrieval":
        behav = pd.read_csv(config["events" + session], sep="\t")
        behav["item"] = pd.Categorical(
            behav["item"], categories=behav["item"].unique(), ordered=True
        )  # otherwise pd.pivot messes up the directions. depreceated?
        behav = pd.pivot(behav, index="item", columns="trial_type", values="accuracy")

        betas = config["STE_" + space + "_" + pf]["STE_delret"]  # load beta files
        spm = loadmat(betas[0].replace("beta_0001.nii", "SPM.mat"))  # get spm file
        beta_ident = spm["SPM"][0][0]["xX"]["name"][0][0][
            0
        ]  # in order to prevent mismatched labels
        beta_names = [
            str(item.flat[0][6:19]) for item in beta_ident[0 : len(betas)]
        ]  # by merging the betas
        beta_files_dict = {
            "item": beta_names,
            "betas": betas,
        }  # with the behavioral file
        beta_files = pd.DataFrame(
            beta_files_dict
        )  # on the item (filename of the presented face, i.e. fol_060253.png)
        df = behav.merge(beta_files, how="inner", on="item")
        df.reset_index(inplace=True)

        sl_mask = config["STE_" + space + "_" + pf]["STE_delret_mask"]

    elif session == "Recognition":
        behav = pd.read_csv(config["eventsDelayedRetrieval"], sep="\t")
        behav["item"] = pd.Categorical(
            behav["item"], categories=behav["item"].unique(), ordered=True
        )  # otherwise pd.pivot messes up the directions. depreceated?
        behav = pd.pivot(
            behav, index="item", columns="trial_type", values="accuracy"
        )  # change to wide format (as many rows as betas.)

        betas = config["STE_" + space + "_" + pf]["STE_recog"]  # load beta files
        spm = loadmat(
            betas[0].replace("beta_0001.nii", "SPM.mat")
        )  # get spm file that created the beta nifti's.
        beta_ident = spm["SPM"][0][0]["xX"]["name"][0][0][
            0
        ]  # in order to prevent mismatched labels
        beta_names = [
            str(item.flat[0][6:19]) for item in beta_ident[0 : len(betas)]
        ]  # by merging the betas
        beta_files_dict = {
            "item": beta_names,
            "betas": betas,
        }  # with the behavioral file
        beta_files = pd.DataFrame(
            beta_files_dict
        )  # on the item (filename of the presented face, i.e. fol_060253.png):
        df = behav.merge(beta_files, how="inner", on="item")
        df.reset_index(inplace=True)

        sl_mask = config["STE_" + space + "_" + pf]["STE_recog_mask"]

    elif session == "Inaccessible":
        behav = pd.read_csv(config["eventsDelayedRetrieval"], sep="\t")

    #### Category: Make final choices for labels and betas depending on the categories you want to look at
    # for example, only choose the unconscious trials if you'd like to decode the category.
    if category == "sex":
        labels = df.item.apply(lambda x: "male" if x[0] == "m" else "female")
        # mask_list = ['fusiform_mask', 'STS_mask', ]  # took iog out, because its not part of func image.
    elif category == "object_category":
        labels = df.item.apply(lambda x: "organic" if x[1] == "o" else "inorganic")
        # mask_list = []  # add brain masks for object category sensitive regions.
    elif category == "face_id":
        labels = df.item
    elif category == "consciousAll":
        if session == "Retrieval" or session == "DelayedRetrieval":
            labels = df["retRatingFeedback"]
        elif session == "Recognition":
            labels = df["recogRatingFeedback"]
        elif session == "Encoding":
            labels = df["accuracy_ret"]
    elif category == "conscious":
        if session == "Retrieval" or session == "DelayedRetrieval":
            df = df[df["retRatingFeedback"] != "neither"]
            labels = df["retRatingFeedback"]
        elif session == "Recognition":
            df = df[
                df["recogRatingFeedback"] != "neither"
            ]  # almost no unconscious trials, we could also add this to raised expections and don't run it at all.
            labels = df["recogRatingFeedback"]
        elif session == "Encoding":
            df = df[df["accuracy_ret"] != "neither"]
            labels = df["accuracy_ret"]
        # mask_list = ['CA1_lh','CA1_rh','CA1_bilat', 'CA3_lh','CA3_rh','CA3_bilat', 'CA4_lh','CA4_rh','CA4_bilat','fusiform_mask', 'STS_mask', ]  # add the hippocampal masks here.
    elif category == "unconscious_accuracy":
        if session == "Retrieval" or session == "DelayedRetrieval":
            df = df[df["retRatingFeedback"] == "unconscious"]
            labels = df["ret"]
        elif session == "Recognition":
            df = df[df["recogRatingFeedback"] == "unconscious"]
            labels = df["recog"]
    elif category == "accuracy":
        if session == "Retrieval" or session == "DelayedRetrieval":
            labels = df["ret"]
        elif session == "Recognition":
            labels = df['recog']
    elif category == "unconscious_object_category":
        if session == "Retrieval" or session == "DelayedRetrieval":
            df = df[df["retRatingFeedback"] == "unconscious"]
        elif session == "Recognition":
            df = df[df["recogRatingFeedback"] == "unconscious"]
        labels = df.item.apply(lambda x: "organic" if x[1] == "o" else "inorganic")
    elif category == "conscious_object_category":
        if session == "Retrieval" or session == "DelayedRetrieval":
            df = df[df["retRatingFeedback"] != "unconscious"]
        elif session == "Recognition":
            df = df[df["recogRatingFeedback"] != "unconscious"]
        labels = df.item.apply(lambda x: "organic" if x[1] == "o" else "inorganic")
    elif category == "buttonpress":
        if session == "Recognition":
            labels = df["recogFeedback"]
        elif session == "ImmRecall":
            df["obj_cat"] = df.item.apply(
                lambda x: "organic" if x[1] == "o" else "inorganic"
            )
            labels = df.apply(lambda x: infer_button_press(x, session), axis=1)
    else:
        ValueError("Category not recognized")
        
    betas = df["betas"]
    labels = labels.reset_index(drop=True)
    betas = betas.reset_index(drop=True)

    return labels, betas, sl_mask


def setup_configNL(subj):
    """
    setup_config - Prepares most other scripts with setup.
    sets up Datapath for the analysis
    if running on mac (assuming that this mac has connection to the
    research storage and able to login with smb connection)
    Syntax: [config] = setup_config(id,qnap)

    Inputs:
    id - participant number

    Outputs:
    config - general inforamtion of subject path setup.

    Initially written for Matlab: Mirko Bristle
    Adapted for python by: Tom Willems

    to work with dataset in the BIDS format
    should now work with other BIDS format datasets.
    """
    # if platform == "darwin":
    #    volume = "/Volumes/"
    # elif platform == "linux":
    #    volume = "/storage/research/"
    # else:
    #    volume = "Z:\\"
    volume = "/storage/research/"

    participantFiles = {}  # create empty dict
    genDerivPath = "/bids/derivatives/"
    derivPath = "/bids/derivatives/fmriprep/"  # maybe not always FMRIPREP_0808 und sowas dabei haben?

    if type(subj) != str:
        subj = str(subj)
    if len(subj) < 6:
        subj = "sub-" + subj

    subjStr = subj
    basePath = volume + "psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/"

    hippParticipants = pd.read_csv(
        basePath + "hipp/bids/" + "participants.tsv", sep="\t"
    )
    wbParticipants = pd.read_csv(basePath + "wb/bids/" + "participants.tsv", sep="\t")

    if hippParticipants["participant_id"].str.contains(subjStr).any():
        study = "hipp"
        participants = hippParticipants
        print("subj from group HIPP")
    elif wbParticipants["participant_id"].str.contains(subjStr).any():
        study = "wb"
        participants = wbParticipants
        print("subj from group WB")
    else:
        raise Exception("OH BOY, this subject is not part of the study!")

    htmlPath = basePath + study + derivPath
    genOutPath = basePath + study + genDerivPath
    participantFiles["workspace"] = (
        "/storage/workspaces/psy_memory_wfg_psy/hpc_henke_wfg/s2019_twillems_fMRI_silent_engram/data/fMRI/"
    )
    participantFiles["genOut"] = genOutPath
    participantFiles["study"] = study

    participantFiles["SLout"] = participantFiles["workspace"] + study
    participantFiles["subStr"] = subjStr
    participantFiles["allSubsStr"] = participants["participant_id"]
    participantFiles["allSubsNum"] = participantFiles["allSubsStr"].apply(
        lambda x: int(x[4:])
    )
    fmriprep_types = [
        "sub*UNIDEN_desc-preproc_T1w.nii*",
        "*MNI*desc-preproc_T1w.nii*",
        "*GM_probse*nii*",
        "*ses-1_task-forgetexp_acq*T1w*preproc*nii*",
        "*ses-1_task-consolidation*T1w*preproc*nii*",
        "*ses-1_task-retrieval*T1w*preproc*nii*",
        "*ses-2_task-forgetexp*T1w*preproc*nii*",
        "*ses-1_task-forgetexp*MNI*preproc*nii*",
        "*ses-1_task-consolidation*MNI*preproc*nii*",
        "*ses-1_task-retrieval*MNI*preproc*nii*",
        "*ses-2_task-forgetexp*MNI*preproc*nii*",
        "*ses-1_task-forgetexp*desc-confounds_timeseries*tsv*",
        "*ses-1_task-consolidation*desc-confounds_timeseries*tsv*",
        "*ses-1_task-retrieval*desc-confounds_timeseries*tsv*",
        "*ses-2_task-forgetexp*desc-confounds_timeseries*tsv*",
        "sub*ses-1_task-forgetexp*events_maxTrials.tsv",
        "sub*ses-1_task-consolidation*events_maxTrials.tsv",
        "sub*ses-1_task-retrieval*events_maxTrials.tsv",
        "sub*ses-2_task-forgetexp*events_maxTrials.tsv",
        "sub*ses-2_to_t1w.mat",
        "sub*ses-1_to_t1w.mat",
        "sub*ses-2_*from-orig_to-T1w_mode-image_xfm.txt",
        "sub*ses-1_*from-orig_to-T1w_mode-image_xfm.txt",
    ]

    fmriprep_names = [
        "nativeT1",
        "mniT1",
        "nativeGMMask",
        "nativeEncoding",
        "nativeConsolidation",
        "nativeRetrieval",
        "nativeDelayedRetrieval",
        "mniEncoding",
        "mniConsolidation",
        "mniRetrieval",
        "mniDelayedRetrieval",
        "confoundsEncoding",
        "confoundsConsolidation",
        "confoundsRetrieval",
        "confoundsDelayedRetrieval",
        "eventsEncoding",
        "eventsConsolidation",
        "eventsRetrieval",
        "eventsDelayedRetrieval",
        "affineMatrix-ses2ToT1w",
        "affineMatrix-ses1ToT1w",
        "affineMatrix-ses2ToT1w_ANTS",
        "affineMatrix-ses1ToT1w_ANTS",
    ]

    spmDeriv_types = [
        "sub*ses-1*UNIDEN_T1w.nii",
        "msub*ses-1*UNIDEN_T1w.nii",
        "rsub*ses-1*UNIDEN_T1w.nii",
        "meansub*ses-1_task-forgetexp*",
    ]
    spmDeriv_names = ["spmCoregAnat", "spmBiasCorr", "spmReslicedAnat", "meanEncFunc"]

    spmMB_types = [
        "rs*T1w.nii",
        "meana*.nii",
        "wra*consolidation*.nii",
        "s6wra*consolidation*.nii",
        "rp_*consolidation*.txt",
        "y_*ses-1*.nii",
        "y_*ses-2*.nii",
    ]
    spmMB_names = [
        "reslicedAnat",
        "meanCoregFunc",
        "consolidationNormFunc",
        "consolidationSmoothNormFunc",
        "spmMotionParamsConsolidation",
        "defSess1",
        "defSess2",
    ]

    nilearn_types = [subjStr + "*RDM_brain.nii", subjStr + "*evalScore.npy"]
    nilearn_names = ["searchlightRDMBrain", "searlightRDMdata"]

    freesurfer_types = [
        "206_mask_lh.nii",
        "206_mask_rh.nii",
        "206_mask_bilat.nii",
        "208_mask_lh.nii",
        "208_mask_rh.nii",
        "208_mask_bilat.nii",
        "209_mask_lh.nii",
        "209_mask_rh.nii",
        "209_mask_bilat.nii",
        "21_bilat.nii*",
        "74_bilat.nii*",
        "02_bilat.nii*",
        "23_bilat.nii*",
        "23_lh.nii*",
        "23_rh.nii*",
        "25_bilat.nii*",
        "26_bilat.nii*",
        "27_bilat.nii*",
        "206_mask_lh_ses-2.nii",
        "206_mask_rh_ses-2.nii",
        "206_mask_bilat_ses-2.nii",
        "208_mask_lh_ses-2.nii",
        "208_mask_rh_ses-2.nii",
        "208_mask_bilat_ses-2.nii",
        "209_mask_lh_ses-2.nii",
        "209_mask_rh_ses-2.nii",
        "209_mask_bilat_ses-2.nii",
        "21_bilat_ses-2.nii*",
        "74_bilat_ses-2.nii*",
        "02_bilat_ses-2.nii*",
        "23_bilat_ses-2.nii*",
        "25_bilat_ses-2.nii*",
        "26_bilat_ses-2.nii*",
        "27_bilat_ses-2.nii*",
        "lh_hipp.nii",
        "rh_hipp.nii",
        "226_mask_lh*",
        "226_mask_rh*",
        "231_mask_lh*",
        "231_mask_rh*",
        "232_mask_lh*",
        "232_mask_rh*",
        "11_bilat.nii*",
        "06_bilat.nii*",
        "16_bilat.nii*",
        "30_bilat.nii*",
        "37_bilat.nii*",
        "38_bilat.nii*",
        "73_bilat.nii*",
        "06_DKT_lh.nii*",
        "06_DKT_rh.nii*",
        "16_DKT_lh.nii*",
        "16_DKT_rh.nii*",
    ]
    freesurfer_names = [
        "CA1_lh",
        "CA1_rh",
        "CA1_bilat",
        "CA3_lh",
        "CA3_rh",
        "CA3_bilat",
        "CA4_lh",
        "CA4_rh",
        "CA4_bilat",
        "fusiform_mask",
        "STS_mask",
        "IOG_mask",
        "parahipp_mask",
        "parahipp_lh_mask",
        "parahipp_rh_mask",
        "angularis_mask",
        "supramarg_mask",
        "supParGyr_mask",
        "CA1_lh_ses-2",
        "CA1_rh_ses-2",
        "CA1_bilat_ses-2",
        "CA3_lh_ses-2",
        "CA3_rh_ses-2",
        "CA3_bilat_ses-2",
        "CA4_lh_ses-2",
        "CA4_rh_ses-2",
        "CA4_bilat_ses-2",
        "fusiform_mask_ses-2",
        "STS_mask_ses-2",
        "IOG_mask_ses-2",
        "Parahipp_mask_ses-2",
        "angularis_mask_ses-2",
        "supramarg_mask_ses-2",
        "supParGyr_mask_ses-2",
        "lh_hipp",
        "rh_hipp",
        "hipp_tail_lh",
        "hipp_tail_rh",
        "hipp_body_lh",
        "hipp_body_rh",
        "hipp_head_lh",
        "hipp_head_rh",
        "cuneus_mask",
        "acc_mask",
        "front_sup_mask",
        "precuneus_mask",
        "inf_temp_gyr_mask",
        "mid_temp_gyr_mask",
        "inferior_temp_sulc_mask",
        "entorhinal_lh_mask",
        "entorhinal_rh_mask",
        "parahipp_lh_no_er_mask",
        "parahipp_rh_no_er_mask", 
    ]
    participantFiles["id"] = subj  # first entry is identifier

    fmriprepPath = htmlPath + subjStr  # subject specific path
    for ind in range(
        0, len(fmriprep_types)
    ):  # loop by ind to also get the prettier name from the names list
        scan_regex = fmriprep_types[ind]  # create variable for the right naming pattern
        scan_name = fmriprep_names[ind]  # create variable for a pretty name
        for root, subdirs, files in os.walk(
            fmriprepPath
        ):  # walk through directories and subdirectories
            for filename in files:
                if fnmatch.fnmatch(
                    filename, scan_regex
                ):  # check for filetypes defined in "types"
                    fullFilePath = os.path.join(
                        root, filename
                    )  # join the fullpath together
                    participantFiles[scan_name] = (
                        fullFilePath  # create new dictionary input with pretty name and the full path for the file
                    )

    spmPath = basePath + study + genDerivPath + "spm/"
    spmPath2 = basePath + study + genDerivPath + "spmMB/"
    spmDerivPath = spmPath2 + subjStr + "/preproc/"
    for ind in range(
        0, len(spmDeriv_types)
    ):  # loop by ind to also get the prettier name from the names list
        scan_regex = spmDeriv_types[ind]  # create variable for the right naming pattern
        scan_name = spmDeriv_names[ind]  # create variable for a pretty name
        for root, subdirs, files in os.walk(
            spmDerivPath
        ):  # walk through directories and subdirectories
            for filename in files:
                if fnmatch.fnmatch(
                    filename, scan_regex
                ):  # check for filetypes defined in "types"
                    fullFilePath = os.path.join(
                        root, filename
                    )  # join the fullpath together
                    participantFiles[scan_name] = (
                        fullFilePath  # create new dictionary input with pretty name and the full path for the file
                    )

    spmMBDerivPath = spmPath2 + subjStr + "/preproc/"

    # participantFiles['normFiles'] = glob.glob(spmMBDerivPath+'**/w*.nii',recursive=True)
    # participantFiles['normFilesSPMold'] = glob.glob(spmPath+subjStr+'/preproc/'+'**/w*.nii',recursive=True)
    # participantFiles['normFilesSPMMBold'] = glob.glob(basePath+study+genDerivPath+'spmMB_old/'+subjStr+'/preproc/'+'**/w*.nii',recursive=True)

    participantFiles["normFuncFiles"] = glob.glob(spmMBDerivPath + "**/func/s3w*.nii")
    participantFiles["normFuncFilesNoSmooth"] = glob.glob(spmMBDerivPath + "**/func/w*.nii")
    participantFiles["normFuncFiles1s"] = glob.glob(spmDerivPath + "**/func/s1w*.nii")
    participantFiles['normAnatFiles'] = glob.glob(spmMBDerivPath + "**/anat/w*.nii")
    participantFiles["normFilesSPMold"] = glob.glob(
        spmPath + subjStr + "/preproc/" + "**/**/w*.nii"
    )
    participantFiles["normFilesSPMMBold"] = glob.glob(
        basePath
        + study
        + genDerivPath
        + "spmMB_old/"
        + subjStr
        + "/preproc/"
        + "**/**/w*.nii"
    )

    for ind in range(
        0, len(spmMB_types)
    ):  # loop by ind to also get the prettier name from the names list
        scan_regex = spmMB_types[ind]  # create variable for the right naming pattern
        scan_name = spmMB_names[ind]  # create variable for a pretty name
        for root, subdirs, files in os.walk(
            spmMBDerivPath
        ):  # walk through directories and subdirectories
            for filename in files:
                if fnmatch.fnmatch(
                    filename, scan_regex
                ):  # check for filetypes defined in "types"
                    fullFilePath = os.path.join(
                        root, filename
                    )  # join the fullpath together
                    participantFiles[scan_name] = (
                        fullFilePath  # create new dictionary input with pretty name and the full path for the file
                    )

    nilearnPath = basePath + study + genDerivPath + "nilearn/"
    for ind in range(
        0, len(nilearn_types)
    ):  # loop by ind to also get the prettier name from the names list
        scan_regex = nilearn_types[ind]  # create variable for the right naming pattern
        scan_name = nilearn_names[ind]  # create variable for a pretty name
        for root, subdirs, files in os.walk(
            nilearnPath
        ):  # walk through directories and subdirectories
            for filename in files:
                if fnmatch.fnmatch(
                    filename, scan_regex
                ):  # check for filetypes defined in "types"
                    fullFilePath = os.path.join(
                        root, filename
                    )  # join the fullpath together
                    participantFiles[scan_name] = (
                        fullFilePath  # create new dictionary input with pretty name and the full path for the file
                    )

    # freesurferPath = basePath + study + genDerivPath + 'freesurfer/' + subjStr
    for sess in ["ses-1", "ses-2"]:
        freesurfer_dict = dict()
        freesurferPath = os.path.join(
            participantFiles["workspace"], study, "freesurfer/", sess, subjStr
        )
        for ind in range(0, len(freesurfer_types)):
            scan_regex = freesurfer_types[ind]
            scan_name = freesurfer_names[ind]
            for root, subdirs, files in os.walk(freesurferPath):
                for filename in files:
                    if fnmatch.fnmatch(
                        filename, scan_regex
                    ):  # check for filetypes defined in "types"
                        fullFilePath = os.path.join(
                            root, filename
                        )  # join the fullpath together
                        freesurfer_dict[scan_name] = (
                            fullFilePath  # create new dictionary input with pretty name and the full path for the file
                        )
        participantFiles["freesurfer_" + sess] = freesurfer_dict

    participantFiles["SingleTrialEncEventsBetasFolder"] = (
        spmPath + subjStr + "/stats/single_trial_events_enc_new"
    )

    STEPath = os.path.join(
        "/storage/workspaces/psy_memory_wfg_psy/hpc_henke_wfg/s2019_twillems_fMRI_silent_engram/data/fMRI/",
        study,
        "1st_level/",
    )
    STE_names = ["STE_enc", "STE_ret", "STE_delret", "STE_recog"]
    spaces = ["native", "mni"]
    for space in spaces:
        if space == "native":
            spaceStr = "_native"
            if study == "hipp":
                prefix = ["s1r", "r"]
            else:
                prefix = ["s1r", "r"]
        elif space == "mni":
            spaceStr = space
            prefix = ["wr"]
        for pf in prefix:
            STE_folders = [
                "/l1_" + spaceStr + "_enc_single_forgetexp_" + pf,
                "/l1_" + spaceStr + "_ret_single_retrieval_" + pf,
                "/l1_" + spaceStr + "_ret_single_forgetexp_" + pf,
                "/l1_" + spaceStr + "_recog_single_forgetexp_" + pf,
            ]
            if space == "native" and study == "wb":
                STE_names.append("STE_ImmRecall")
                STE_folders = [folder + "PPI" for folder in STE_folders]
                STE_folders.append(
                    "/l1__native_enc_single_forgetexp_" + pf + "MaskedbuttonPress"
                )

            space_dict = dict()
            for ind in range(0, len(STE_folders)):
                scan_folder = STEPath + subjStr + STE_folders[ind]
                # print('scan_folder:',scan_folder)
                files = glob.glob(scan_folder + "/beta*.nii")
                files.sort()
                files = files[:-7]
                # print('files:',files[0])
                mask_img = glob.glob(scan_folder + "/mask.nii")
                if pf == "s1r":
                    residuals = glob.glob(scan_folder + "/Res_0*.nii")
                    residualMS = glob.glob(scan_folder + "/ResMS.nii")
                    residuals.sort()
                    # residuals = residuals[:-7]
                    # beta_df  = pd.DataFrame(files,residuals, columns = ['betas','residuals'])
                    tMaps = glob.glob(scan_folder + "/spmT_*")
                    tMaps.sort()
                    space_dict[STE_names[ind] + "_residuals"] = residuals
                    space_dict[STE_names[ind] + "_tMaps"] = tMaps
                    space_dict[STE_names[ind] + "_ResMS"] = residualMS

                space_dict[STE_names[ind]] = files
                space_dict[STE_names[ind] + "_mask"] = mask_img
            participantFiles["STE_" + space + "_" + pf] = space_dict

    if study == "hipp":
        FL_names = ["Encoding", "Retrieval", "DelayedRetrieval"]
        FL_folders = [
            "l1__native_enc_subsequent_forgetexp_s1r_RSAready",
            "l1__native_ret_retrieval_s1r_RSAready",
            "l1__native_delret_forgetexp_s1r_RSAready",
        ]
        participantFiles["nonSTE_RSAready"] = dict()

        for ind in range(0, len(FL_folders)):
            sessionDict = dict()
            scan_folder = STEPath + subjStr + "/" + FL_folders[ind]
            # print('scan_folder:',scan_folder)
            files = glob.glob(scan_folder + "/beta*.nii")
            files.sort()
            files = files[:-7]
            # print('files:',files[0])
            mask_img = glob.glob(scan_folder + "/mask.nii")
            residuals = glob.glob(scan_folder + "/Res_0*.nii")
            residualMS = glob.glob(scan_folder + "/ResMS.nii")
            residuals.sort()
            # residuals = residuals[:-7]
            # beta_df  = pd.DataFrame(files,residuals, columns = ['betas','residuals'])
            tMaps = glob.glob(scan_folder + "/spmT_*")
            tMaps.sort()

            spm = loadmat(scan_folder + "/SPM.mat")
            poplist = list()
            for t_ind in range(0, len(spm["SPM"]["xCon"][0][0][0])):
                tMap_fromSPM = spm["SPM"]["xCon"][0][0][0][t_ind]["Vspm"]["descrip"][0][
                    0
                ][0][
                    spm["SPM"]["xCon"][0][0][0][t_ind]["Vspm"]["descrip"][0][0][0].find(
                        ":"
                    )
                    + 2 :
                ]
                beta_ident = spm["SPM"][0][0]["xX"]["name"][0][0][0]
                beta_names = [
                    str(item.flat[0][item.flat[0].find(")") + 2 : -6])
                    for item in beta_ident[0 : len(beta_ident[:-7])]
                ]  # by merging the betas
                if tMap_fromSPM not in beta_names:
                    poplist.append(t_ind)

            tMaps = [v for i, v in enumerate(tMaps) if i not in frozenset((poplist))]

            sessionDict[FL_names[ind] + "_residuals"] = residuals
            sessionDict[FL_names[ind] + "_tMaps"] = tMaps
            sessionDict[FL_names[ind] + "_ResMS"] = residualMS
            sessionDict[FL_names[ind] + "_files"] = files
            sessionDict[FL_names[ind] + "_mask"] = mask_img
            participantFiles["nonSTE_RSAready"][FL_names[ind]] = sessionDict

    eventsFilePath = basePath + "/utilities/eventFiles/"
    encodingSingleEventsFilesPath = eventsFilePath + "/single_trial_events_encoding/"

    participantFiles["SingleTrialEncEvents"] = (
        encodingSingleEventsFilesPath + subjStr + "_subsequent_consc.csv"
    )
    print("got all files for :" + subjStr)
    return participantFiles

    # behav_ret1 = pd.read_csv(config['eventsRetrieval'], sep="\t")
    # behav_ret2 = pd.read_csv(config['eventsDelayedRetrieval'], sep="\t")
    # acc_ret1 = behav_ret1[behav_ret1['trial_type'] == 'ret']['accuracy'].reset_index(drop=True)
    # acc_ret2 = behav_ret2[behav_ret2['trial_type'] == 'ret']['accuracy'].reset_index(drop=True)
    # behav_ret1 = behav_ret1[behav_ret1['trial_type'] == 'retRatingFeedback'].reset_index(drop=True)
    # behav_ret2 = behav_ret2[behav_ret2['trial_type'] == 'retRatingFeedback'].reset_index(drop=True)
    # behav_ret1.loc[:, 'betas'] = config['STE_ret']  # add beta pathnames
    # behav_ret2.loc[:, 'betas'] = config['STE_delret']  # add beta pathnames
    # behav_ret1.loc[:, 'label'] = acc_ret1
    # behav_ret2.loc[:, 'label'] = acc_ret2
    # behav_ret1 = behav_ret1[behav_ret1['accuracy'] == 'unconscious']
    # behav_ret2 = behav_ret2[behav_ret2['accuracy'] == 'unconscious']
    # if multiSession == "multi":
    # df = pd.concat([behav_ret1, behav_ret2], ignore_index=True, axis=0)
    # okay but what to do with masks. probably make masks part of this script as well... and then return it as part of df...
    # especially if you run over multiple sessions, you want different masks.... okay yes do this, shouldn't be too hard
