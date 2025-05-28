# %% load libraries
import sys,os
sys.path.append(os.environ['HOME'])
from TOAM.fmri.mvpa.mvpa_utils.get_datasets import *
import numpy as np
from nilearn.glm.first_level import FirstLevelModel


def get_recognition_betas(label, config, smooth, force):
    """description of function

    Args:
        label (string): which label you want to get back, i.e. whether recognition was conscious or not. put 'recogRatingFeedback' for conscious or not
        config (object): just the config object with all the necessary information and paths for this subject, including the sequence
        smooth (integer): mm by which you want the first level to smooth the data
        force (bool): wether you want to force the new calculation of the first level
    """

    preproc_folder = Path(f'{config["genOut"]}/{config["id"]}/preproc/ses-2/func')
    preproc_files = 's1rasub-60601_ses-2_task-forgetexp_acq-epfid2d1200_bold.nii'
    df = pd.read_csv(config["eventsDelayedRetrieval"], sep="\t")
    labels = df[df['trial_type'] == label]['accuracy']
    items = df[df["trial_type"] == "recog"]["item"]
    func_img = preproc_folder / preproc_files
    out= (f'{config["workspace"]}/{config["study"]}/1st_level/{config["id"]}/l1_recog_single_s{smooth}r')

    ret0l,ret0b, _= create_mv_data('Retrieval','face_id',config)
    retl, retb, _ = create_mv_data('DelayedRetrieval','face_id',config)
    recl, recb, _ = create_mv_data('Recognition','face_id',config)

    retl2, _, _   = create_mv_data("DelayedRetrieval","accuracy",config)
    recl2, _, _   = create_mv_data("Recognition","accuracy",config)

    retl3, _, _   = create_mv_data("DelayedRetrieval","consciousAll",config)
    retl4, _, _   = create_mv_data("Recognition","consciousAll",config)


    beta_filenames='betas_recog_'+label+'_smooth_'+str(smooth)+'.npy'
    betafiles= outdir / beta_filenames
    if betafiles.exists() and not force:
        betas = np.load(betafiles)
    else:
        # create first level object
        fl_glm = FirstLevelModel(
            t_r=config.tr,
            noise_model="ar1",
            smoothing_fwhm=smooth,
            standardize=False,
            hrf_model="spm",
            drift_model="cosine",
            high_pass=0.01,
        )

        fl_glm.fit(func_img, design_matrices=recog_events)
        betas = fl_glm.get_params()

    return labels, betas

if __name__ == '__main__':
    get_recognition_betas('recogRatingFeedback', setup_configNL(60100), 6, False)
