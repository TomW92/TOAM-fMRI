#!/bin/bash
#SBATCH -J fmriprep_slurm_file
#SBATCH --time=24:00:00
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
# Outputs ----------------------------------
#SBATCH -o /storage/homefs/tw18a205/log/%x-%A-%a.out
#SBATCH -e /storage/homefs/tw18a205/log/%x-%A-%a.err

#SBATCH --mail-user=tom.willems@unibe.ch
#SBATCH --mail-type=fail,end,ARRAY_TASKS
# ------------------------------------------

#load workspace to make sure $SCRATCH is defined
echo "TMPDIR="$TMPDIR
module load Workspace
export TMPDIR=$SCRATCH

#/storage/research/psy_wfg_henke/s2023_kzervas_twillems_fMRI_supraliminal/data/fmri/pilot/raw/sub-002
BIDS_DIR="/storage/research/psy_wfg_henke/s2023_kzervas_twillems_fMRI_supraliminal/data/fmri/memoryTrace_remoteRecall/bids"
DERIVS_DIR="derivatives/fMRIprepNoSbref" #fMRIprepNoSbref
FREESURFER_DIR="derivatives/freesurfer"

## TRY OUT CLEAN WORKDIR INSTEAD IF NECESSARY and check SEFF for 1 and 2 fmriprepNoSBREF

# Prepare some writeable bind-mount points.
TEMPLATEFLOW_HOST_HOME=$WORKSPACE/.cache/templateflow #Changed this from $HOME to $WORKSPACE
TEMPLATEFLOW_HOME=$WORKSPACE/.templateflow #ADDED this and did the changes around to account for missing access to cache from other lab members
FMRIPREP_HOST_CACHE=$WORKSPACE/.cache/fmriprep #Changed this from $HOME to $WORKSPACE
mkdir -p ${TEMPLATEFLOW_HOST_HOME}
mkdir -p ${FMRIPREP_HOST_CACHE}

# Prepare derivatives folder
mkdir -p ${BIDS_DIR}/${DERIVS_DIR}
mkdir -p ${BIDS_DIR}/${FREESURFER_DIR}

LOCAL_FREESURFER_DIR="${BIDS_DIR}/${FREESURFER_DIR}"

# Make sure FS_LICENSE is defined in the container.
#
export FREESURFER=$HOME/toolboxes/freesurfer
export SINGULARITYENV_FS_LICENSE=$HOME/toolboxes/freesurfer/fs_license.txt
echo $HOME
# Designate a templateflow bind-mount point
#export SINGULARITYENV_TEMPLATEFLOW_HOME="/templateflow"
SINGULARITY_CMD="singularity run --cleanenv \
  -B $BIDS_DIR:/bids \
  -B $SCRATCH:/work -B ${LOCAL_FREESURFER_DIR}:/fsdir \
  $HOME/my_images/fmriprep-23-2-0.simg" #$HOME/my_images/fmriprep-20-2-3.simg"

  #-B ${TEMPLATEFLOW_HOST_HOME}:${SINGULARITYENV_TEMPLATEFLOW_HOME} \

# Parse the participants.tsv file and extract one subject ID
# from the line corresponding to this SLURM task.
subject=$( sed -n -E "$((${SLURM_ARRAY_TASK_ID} + \
  1))s/sub-(\S*)\>.*/\1/gp" ${BIDS_DIR}/participants.tsv ) 


# Remove IsRunning files from FreeSurfer
find ${LOCAL_FREESURFER_DIR}/sub-$subject/ -name "*IsRunning*" -type f -delete


# Compose the command line
cmd="${SINGULARITY_CMD} /bids /bids/${DERIVS_DIR} participant --participant-label \
  $subject -w /work/ -vv --omp-nthreads 12 --nthreads 16 --mem_mb 60000 \
  --fs-subjects-dir /fsdir \
  --bold2t1w-init register  \
  --output-spaces MNI152NLin2009cAsym anat \
  --ignore {sbref,t2w}"
  
# --ignore sbref
# --use-aroma deprecated
# --output-spaces MNI152NLin2009cAsym:res-2 anat fsnative fsaverage5 \
# --bold2t1w-init header --use-plugin $HOME/fmri/plugin.yml
# Either register the default to initialize volumes at center or header to use the header information when coregistering BOLD to T1w images.

# Setup done, run the command
echo Running task ${SLURM_ARRAY_TASK_ID}
echo Commandline: $cmd
eval $cmd
exitcode=$?

# Output results to a table file 
echo "sub-$subject   ${SLURM_ARRAY_TASK_ID}    $exitcode" \
      >> ${SLURM_JOB_NAME}.${SLURM_ARRAY_JOB_ID}.tsv
echo Finished tasks ${SLURM_ARRAY_TASK_ID} with exit code $exitcode
exit $exitcode
