#!/bin/bash
#SBATCH --mail-user=tom.willems@unibe.ch
#SBATCH --mail-type=fail,end
#SBATCH --job-name=psth
#
#
#SBATCH --time=04:59:01
#SBATCH --mem-per-cpu=8G
#SBATCH --tmp=6048
#
#  Partition 
#SBATCH --partition=bdw

# parallel jobs
#SBATCH --ntasks=20
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

module load MATLAB

INPUT=$1
STUDY="$STORAGE/$INPUT"
BIDS_DIR="$STUDY/bids"
subject=$( sed -n -E "$((SLURM_ARRAY_TASK_ID + 1))s/sub-(\S*)\>.*/\1/gp" "${BIDS_DIR}"/participants.tsv )
srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nodisplay -singleCompThread -batch "cd('/storage/homefs/tw18a205/TOAM/fmri/functions'); normalisation_spm($subject); exit;"
#srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nodisplay -singleCompThread -batch "cd('/storage/homefs/tw18a205/TOAM/fmri/functions'); normalize_native_decoding_results($subject); settings; exit;"
#srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nodisplay -singleCompThread -batch "cd('/storage/homefs/tw18a205/TOAM/fmri/functions'); extract_PSTH_version_kosta; settings; exit;"


