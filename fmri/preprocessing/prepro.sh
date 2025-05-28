#!/bin/bash
#SBATCH --mail-user=tom.willems@unibe.ch
#SBATCH --mail-type=fail,end #,ARRAY_TASKS
#SBATCH --job-name=prepro
#
#
#SBATCH --time=03:30:01
#SBATCH --mem-per-cpu=12G
#SBATCH --tmp=2048
#
#  Partition
#### SBATCH --partition=bdw

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

module load MATLAB/2021b
shortStudy=$(basename "$STUDY")
BIDS_DIR="$STUDY/bids"
subject=$( sed -n -E "$((${SLURM_ARRAY_TASK_ID} + 1))s/sub-(\S*)\>.*/\1/gp" ${BIDS_DIR}/participants.tsv )
echo $subject

#srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nosplash -nodisplay -singleCompThread -r "settings; run('/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/utilities/gunziploopHIPP($subject, \"$shortStudy\")'); exit;"
#wait
srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nosplash -nodisplay -singleCompThread -r "settings; run('/storage/homefs/tw18a205/TOAM/fmri/preprocessing/smallFOV_preprocessing_1($subject)'); exit;"
srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nosplash -nodisplay -singleCompThread -r "settings; run('/storage/homefs/tw18a205/TOAM/fmri/preprocessing/smallFOV_preprocessing_1s($subject)'); exit;"
wait
srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nosplash -nodisplay -singleCompThread -r "settings; run('/storage/homefs/tw18a205/TOAM/fmri/preprocessing/smallFOV_preprocessing_2($subject, 1)'); exit;"
