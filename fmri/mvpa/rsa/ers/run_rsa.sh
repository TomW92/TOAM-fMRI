#!/bin/bash
#SBATCH --mail-user=tom.willems@unibe.ch
#SBATCH --mail-type=fail,end
#SBATCH --job-name=aa_hipp_failed_recog_rsa
#SBATCH --partition=epyc2,bdw
#SBATCH --time=2:30:01
#SBATCH --array=101,103,116,118,11,131,133,13,146,148,161,163,176,178,191,193,206,208,221,223,236,238,251,253,266,268,26,281,283,28,296,298,41,43,56,58,71,73,86,88
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=32G
#SBATCH --tmp=4048

#SBATCH -o /storage/homefs/tw18a205/log/%x-%A-%a.out

echo "TMPDIR="$TMPDIR
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURM_ARRAYID="$SLURM_ARRAYID
echo "SLURM_ARRAY_JOB_ID"=$SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID"=$SLURM_ARRAY_TASK_ID
echo "working directory = "$SLURM_SUBMIT_DIR

git log --decorate=short --oneline -1
HPC_WORKSPACE=psy_wfg_henke module load Workspace
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate toam
exp_id=$(date "+%Y%m%d_%H")

#take one input, which is which study to run this in (hipp or wb)
INPUT=$1
export STUDY="$STORAGE/$INPUT"
BIDS_DIR="$STUDY/bids"
subject=$( sed -n -E "$((${SLURM_ARRAY_TASK_ID} + 1))s/sub-(\S*)\>.*/\1/gp" ${BIDS_DIR}/participants.tsv )

srun python /storage/homefs/tw18a205/TOAM/fmri/mvpa/rsa/ers/encoding_retrieval_sim.py 

echo "finished enc_recog_rsa"
