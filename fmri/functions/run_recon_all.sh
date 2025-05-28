#!/bin/bash
#SBATCH --mail-user=tom.willems@unibe.ch
#SBATCH --mail-type=fail,end,ARRAY_TASKS
#SBATCH --job-name=recon-all-clean
#
#SBATCH --time=24:00:01
#SBATCH --mem-per-cpu=16G
#SBATCH --tmp=2048
#
# Partition
#### SBATCH --partition=all
# parallel jobs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20 #openmp
#SBATCH --time=48:00:00

#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH -o /storage/homefs/tw18a205/log/%x-%A-%a.out

#paths and preparation
git log --decorate=short --oneline -1
INPUT=$1
STUDY="$STORAGE/$INPUT"
BIDS_DIR="$STUDY/bids"
ses=ses-1
shortStudy=$(basename "$STUDY")
subject=$( sed -n -E "$((${SLURM_ARRAY_TASK_ID} + 1))s/sub-(\S*)\>.*/\1/gp" ${BIDS_DIR}/participants.tsv )
subject=sub-${subject}
export SUBJECTS_DIR=${WORKSPACE}/s2019_twillems_fMRI_silent_engram/data/fMRI/${shortStudy}/freesurfer/${ses}
echo $SUBJECTS_DIR

bias_corrected=${BIDS_DIR}/derivatives/spmMB/${subject}/preproc/${ses}/anat/m${subject}_${ses}_acq-t1csmp2ragesag06mmUNIDEN_T1w.nii
echo $bias_corrected
#execute the actual command
#recon-all -subject ${subject} -i ${bias_corrected} -all -hires -openmp 20 -debug

#try fix for sub 17
recon-all -subject ${subject} -i ${bias_corrected} -all -hires -openmp 20 -clean -wsless

#recon-all -subject ${subject} -all -hires -openmp 20 -debug

# earlier options
# -autorecon1 -noskullstrip -mris_inflate -n 50
# -autorecon2 -volonly
# -autorecon3 -nosphere
#             -nosurfreg
#             -nojacobian_white
#             -noavgcurv
#             -nocortparc
#             -nopial
#             -noparcstats
#             -nocortparc2
#             -noparcstats2
#             -nocortparc3
#             -noparcstats3
#             -nopctsurfcon
#             -nocortribbon
#             -nobalabels



