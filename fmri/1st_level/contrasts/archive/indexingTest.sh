STORAGE=/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI

INPUT=$1
subjectNumber=$2
STUDY="$STORAGE/$INPUT"
BIDS_DIR="$STUDY/bids"

# get all directorie names in bids dir
subjects=($(ls -d ${BIDS_DIR}/sub-*))

# get only the subject number
subject=${subjects[$subjectNumber]##*-}

echo "Running subject $subjectNumber: $subject"