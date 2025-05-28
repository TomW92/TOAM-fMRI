#!/bin/bash
#SBATCH --mail-user=tom.willems@unibe.ch
#SBATCH --mail-type=start,fail,end
#SBATCH --job-name=xx_lsa_first_level
#SBATCH --time=4:59:01
#SBATCH --mem-per-cpu=12G
#SBATCH --tmp=2048
#
### SBATCH --partition=epyc2

# parallel jobs
#SBATCH --ntasks=1
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

git log --decorate=short --oneline -1

module load MATLAB/2023b

STORAGE=/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI

INPUT=$1
STUDY="$STORAGE/$INPUT"
BIDS_DIR="$STUDY/bids"
subject=$( sed -n -E "$((${SLURM_ARRAY_TASK_ID} + 1))s/sub-(\S*)\>.*/\1/gp" ${BIDS_DIR}/participants.tsv )


#univariate MNI multisession FL
###### srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nodisplay -singleCompThread -batch "cd('/storage/homefs/tw18a205/TOAM/fmri/1st_level'); multi_session_first_level($subject,'expl_mask'); exit;"
###### srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nodisplay -singleCompThread -batch "cd('/storage/homefs/tw18a205/TOAM/fmri/1st_level'); multi_session_first_level($subject,'expl_mask'); exit;"

#single trial estimates for RSA, multivariate analyses
#srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nodisplay -singleCompThread -batch "cd('/storage/homefs/tw18a205/TOAM/fmri/1st_level'); first_level($subject,1,'forgetexp',0,'enc_single','enc_bids',1,0); settings; exit;"
#srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nodisplay -singleCompThread -batch "cd('/storage/homefs/tw18a205/TOAM/fmri/1st_level'); first_level($subject,1,'retrieval',0,'ret_single','ret_bids',1,0); settings; exit;"
#srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nodisplay -singleCompThread -batch "cd('/storage/homefs/tw18a205/TOAM/fmri/1st_level'); first_level($subject,2,'forgetexp',0,'ret_single','delret_bids',1,0); settings; exit;"
#srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nodisplay -singleCompThread -batch "cd('/storage/homefs/tw18a205/TOAM/fmri/1st_level'); first_level($subject,2,'forgetexp',0,'recog_single','delret_bids',1,0); settings; exit;"

# subsequent recognition trials!
#srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nodisplay -singleCompThread -batch "cd('/storage/homefs/tw18a205/TOAM/fmri/1st_level'); first_level($subject,2,'forgetexp',1,'delret_subsequent','delret_subseq',1,1); settings; exit;"

