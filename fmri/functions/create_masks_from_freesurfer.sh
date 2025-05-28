#!/bin/bash
#SBATCH --mail-user=tom.willems@unibe.ch
#SBATCH --mail-type=fail,end
#SBATCH --job-name=XX_FS_Masks
#
#SBATCH --time=00:40:01
#SBATCH --mem-per-cpu=16G
#SBATCH --tmp=2048
#SBATCH --qos=job_cpu
#SBATCH --partition=epyc2,bdw

# parallel jobs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -o /storage/homefs/tw18a205/log/%x-%A-%a.out

echo "SLURM_JOBID"="$SLURM_JOBID"
echo "SLURM_JOB_NODELIST"="$SLURM_JOB_NODELIST"
echo "SLURM_NNODES"="$SLURM_NNODES"
echo "SLURMTMPDIR"="$SLURMTMPDIR"
echo "SLURM_ARRAYID"="$SLURM_ARRAYID"
echo "SLURM_ARRAY_JOB_ID"="$SLURM_ARRAY_JOB_ID"
echo "SLURM_ARRAY_TASK_ID"="$SLURM_ARRAY_TASK_ID"
echo "working directory = ""$SLURM_SUBMIT_DIR"


git log --decorate=short --oneline -1
HPC_WORKSPACE=hpc_henke_wfg module load Workspace

#take one input, which is which study to run this in (hipp or wb)
INPUT=$1
STUDY="$STORAGE/$INPUT"
BIDS_DIR="$STUDY/bids"
shortStudy=$(basename "$STUDY")
subject=$( sed -n -E "$((SLURM_ARRAY_TASK_ID + 1))s/sub-(\S*)\>.*/\1/gp" "${BIDS_DIR}"/participants.tsv )
subject=sub-${subject}
echo "$subject"

#loop over sessions
for ses in {ses-1,ses-2}
do
  export SUBJECTS_DIR=${WORKSPACE}/s2019_twillems_fMRI_silent_engram/data/fMRI/${shortStudy}/freesurfer/${ses}
  echo "$SUBJECTS_DIR"

  # shellcheck disable=SC2164
  cd "$HOME" || exit

  # check if recon-all failed. If so, exit
  # the line below is only present in the logfile if recon-all exited with errors
  EXPECTED_PATTERN="To report a problem"
  LOG_FILE="${SUBJECTS_DIR}/${subject}/scripts/recon-all.log"
  ACTUAL_LINE=$(tail -n 1 "$LOG_FILE")
  echo "$ACTUAL_LINE"
  if [[ $ACTUAL_LINE == "$EXPECTED_PATTERN"* ]]; then
    echo "Error: "Last line of "$LOG_FILE" starts with "$EXPECTED_PATTERN"
    echo "Freesurfer recon all failed, moving on to next session"
    echo -e "\n"
    continue
  else
    echo "Freesurfer recon all ran successfully"
    # continue with the rest of your script here
  fi

  if [ -f "${SUBJECTS_DIR}"/"${subject}"/scripts/IsRunningHPsubT1.lh+rh ]; then
    echo "Is runningFile exists. deleting attempt."
    rm "${SUBJECTS_DIR}"/"${subject}"/scripts/IsRunningHPsubT1.lh+rh
  else
    echo "Is runningFile doesn't exist. continuing."
  fi 

  if [ ! -f "${SUBJECTS_DIR}"/"${subject}"/mri/lh.hippoAmygLabels-T1.v21.CA.mgz ]; then
    echo "segmentHA_T1.sh didn't run yet, executing..."
    #segment images
    segmentHA_T1.sh "$subject"
  else
    echo "segmentHA_T1.sh already ran, continuing."
  fi
  echo -e "\n\n"
  # t1w space from fmriprep is an unbiased template created with freesurfers _mri_robust_template_ 
  # and thus niether session one anat nor session 2 anat is necessarily aligned to it (the input images)
  # the surface reconstructions are also put into that space and not the space of any of the input images. 
  # so, I need to rewrite the hipp labels to the input in order to extract the right voxels.
  # turn labels to volume (from: https://neurobren.com/freesurfer-masks/ )

  cd "$SUBJECTS_DIR" || exit
  cd "${subject}"/mri || exit

  #create CA masks
  # for label in {206,208,209} #CA1,CA3,CA4;
  # do
  #     #right hemisphere
  #     mri_label2vol \
  #     --seg rh.hippoAmygLabels-T1.v21.CA.FSvoxelSpace.mgz \
  #     --temp rawavg.mgz \
  #     --regheader rh.hippoAmygLabels-T1.v21.CA.FSvoxelSpace.mgz \
  #     --o rh_hippAmygLabels-CA-T1.mgz

  #     mri_binarize --i rh_hippAmygLabels-CA-T1.mgz --o "${label}"_mask_rh.mgz --match "$label"

  #     mri_convert --in_type mgz --out_type nii "${label}"_mask_rh.mgz "${label}"_mask_rh.nii

  #     #left hemisphere
  #     mri_label2vol \
  #     --seg lh.hippoAmygLabels-T1.v21.CA.FSvoxelSpace.mgz \
  #     --temp rawavg.mgz \
  #     --regheader lh.hippoAmygLabels-T1.v21.CA.FSvoxelSpace.mgz \
  #     --o lh_hippAmygLabels-CA-T1.mgz

  #     mri_binarize --i lh_hippAmygLabels-CA-T1.mgz --o "${label}"_mask_lh.mgz --match "$label"

  #     mri_convert --in_type mgz --out_type nii "${label}"_mask_lh.mgz "${label}"_mask_lh.nii

  #     #combine hemispheres
  #     mri_concat --i "${label}"_mask_lh.nii --i "${label}"_mask_rh.nii --o "${label}"_mask_bilat.nii --combine
  # done

  # # create Head/Body/Tail masks
  # for label in {226,231,232} #tail body head
  # do
  #   #right hemisphere
  #   mri_label2vol \
  #   --seg rh.hippoAmygLabels-T1.v21.HBT.FSvoxelSpace.mgz \
  #   --temp rawavg.mgz \
  #   --regheader rh.hippoAmygLabels-T1.v21.HBT.FSvoxelSpace.mgz \
  #   --o rh_hippAmygLabels-HBT-T1.mgz

  #   mri_binarize --i rh_hippAmygLabels-HBT-T1.mgz --o "${label}"_mask_rh.mgz --match "$label"

  #   mri_convert --in_type mgz --out_type nii "${label}"_mask_rh.mgz "${label}"_mask_rh.nii

  #   #left hemisphere  
  #   mri_label2vol \
  #   --seg lh.hippoAmygLabels-T1.v21.HBT.FSvoxelSpace.mgz \
  #   --temp rawavg.mgz \
  #   --regheader lh.hippoAmygLabels-T1.v21.HBT.FSvoxelSpace.mgz \
  #   --o lh_hippAmygLabels-HBT-T1.mgz

  #   mri_binarize --i lh_hippAmygLabels-HBT-T1.mgz --o "${label}"_mask_lh.mgz --match "$label"

  #   mri_convert --in_type mgz --out_type nii "${label}"_mask_lh.mgz "${label}"_mask_lh.nii

  #   #combine hemispheres
  #   mri_concat --i "${label}"_mask_lh.nii --i "${label}"_mask_rh.nii --o "${label}"_mask_bilat.nii --combine
  # done

  # #create cortical masks 
  # # inferior occ- gyrus (02), ACC (06), cuneus (11), frontalsupgyrus (16), fusiform gyrus (21), parahipp (23), angularis(25), supramarginalis(26), superior parietal gyrus(27), 
  # # precuneus (30), inferior temporal gyrus (37), middle temporal gyrus (38),inferior temporal SULCUS (73),superior temporal SULCUS (74),
  # for label in {02,06,11,16,21,23,25,26,27,30,37,38,73,74} 
  # do
  #   mri_label2vol \
  #       --seg aparc.a2009s+aseg.mgz  \
  #       --temp rawavg.mgz \
  #       --regheader aparc.a2009s+aseg.mgz  \
  #       --o aparc.a2009s+aseg_t1space.mgz

  #   mri_binarize \
  #       --i aparc.a2009s+aseg_t1space.mgz \
  #       --o "${label}"_lh.mgz \
  #       --match 111"${label}"

  #   mri_convert --in_type mgz --out_type nii "${label}"_lh.mgz "${label}"_lh.nii.gz

  #   mri_binarize \
  #       --i aparc.a2009s+aseg_t1space.mgz \
  #       --o "${label}"_rh.mgz \
  #       --match 121"${label}"

  #   mri_convert --in_type mgz --out_type nii "${label}"_rh.mgz "${label}"_rh.nii.gz

  #   mri_binarize \
  #       --i aparc.a2009s+aseg_t1space.mgz \
  #       --o "${label}"_bilat.mgz \
  #       --match 111"${label}" 121"${label}"

  #   mri_convert --in_type mgz --out_type nii "${label}"_bilat.mgz "${label}"_bilat.nii.gz
  # done
  
  for label in {06,16} #06 entorhinal, 16 parahipp
  do
    mri_label2vol \
        --seg aparc.DKTatlas+aseg.mgz  \
        --temp rawavg.mgz \
        --regheader aparc.DKTatlas+aseg.mgz  \
        --o aparc.DKTatlas+aseg_t1space.mgz

    mri_binarize \
        --i aparc.DKTatlas+aseg_t1space.mgz \
        --o "${label}"_DKT_lh.mgz \
        --match 10"${label}"

    mri_convert --in_type mgz --out_type nii "${label}"_DKT_lh.mgz "${label}"_DKT_lh.nii.gz

    mri_binarize \
        --i aparc.DKTatlas+aseg_t1space.mgz \
        --o "${label}"_DKT_rh.mgz \
        --match 20"${label}"

    mri_convert --in_type mgz --out_type nii "${label}"_DKT_rh.mgz "${label}"_DKT_rh.nii.gz

    mri_binarize \
        --i aparc.DKTatlas+aseg_t1space.mgz \
        --o "${label}"_DKT_bilat.mgz \
        --match 10"${label}" 20"${label}"

    mri_convert --in_type mgz --out_type nii "${label}"_DKT_bilat.mgz "${label}"_DKT_bilat.nii.gz
  done

  # # create hippocampus masks
  # mri_binarize \
  #   --i aparc.a2009s+aseg_t1space.mgz \
  #   --o lh_hipp.mgz \
  #   --match 17

  # mri_convert --in_type mgz --out_type nii lh_hipp.mgz lh_hipp.nii.gz
  
  # mri_binarize \
  #   --i aparc.a2009s+aseg_t1space.mgz \
  #   --o rh_hipp.mgz \
  #   --match 53  

  # #convert to nifti
  # mri_convert --in_type mgz --out_type nii rh_hipp.mgz rh_hipp.nii.gz
  
  # #create bilateral mask
  # mri_concat --i lh_hipp.nii.gz --i rh_hipp.nii.gz --o hipp_bilat.nii --combine

done 

echo "done"
