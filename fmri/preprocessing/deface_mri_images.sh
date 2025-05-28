#!/bin/bash
#SBATCH --mail-user=tom.willems@unibe.ch
#SBATCH --mail-type=fail,end #,ARRAY_TASKS
#SBATCH --job-name=deface
#
#SBATCH --time=03:30:01
#SBATCH --mem-per-cpu=12G
#SBATCH --tmp=2048
#
# parallel jobs

#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH -o /storage/homefs/tw18a205/log/%x-%A-%a.out

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "SLURM_ARRAYID="$SLURM_ARRAYID
echo "SLURM_ARRAY_JOB_ID"=$SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID"=$SLURM_ARRAY_TASK_ID
echo "working directory = "$SLURM_SUBMIT_DIR
echo "slurm submit dir= "$SLURM_SUBMIT_DIR
git log --decorate=short --oneline -1

INPUT=$1
STUDY="$STORAGE/$INPUT"
BIDS_DIR="$STUDY/bids"

subject=$( sed -n -E "$((${SLURM_ARRAY_TASK_ID} + 1))s/sub-(\S*)\>.*/\1/gp" ${BIDS_DIR}/participants.tsv )
echo $subject

module load Anaconda3
eval "$(conda shell.bash hook)"
conda init 
conda activate toam
module load Workspace
module load fsl


export TMPDIR=$SCRATCH/tmp_${SLURM_JOB_ID}
mkdir -p $TMPDIR
export NIPYPE_WORK_DIR=$TMPDIR

cd /storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/${INPUT}/bids/sub-${subject}/ses-1/anat/
pydeface --force --nocleanup sub-${subject}_ses-1_acq-t1csmp2ragesag06mmUNIDEN_T1w.nii.gz 
    
cd /storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/${INPUT}/bids/sub-${subject}/ses-2/anat/
pydeface --force --nocleanup sub-${subject}_ses-2_acq-t1csmp2ragesag06mmUNIDEN_T1w.nii.gz 