from TOAM.fmri.mvpa.mvpa_utils.get_datasets import *

config = setup_configNL(60601)
subjects = config["allSubsStr"]
for subj in subjects:
    config = setup_configNL(subj)
    combinations = [
        {"object_category": "Encoding"},
        {"object_category": "Retrieval"},
        {"object_category": "DelayedRetrieval"},
    ]
    for combination in combinations:
        ((category, session),) = combination.items()
        # get all paths for this subject

        labels, betas, sl_mask = create_mv_data(
            session, category, config, "native", pf="s1r"
        )
        if session == "DelayedRetrieval":
            session_def = "day2"
        else:
            session_def = "day1"

        normalize_SPM_segmented_data(betas, session_def, config)
