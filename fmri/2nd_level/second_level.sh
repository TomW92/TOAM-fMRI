#!/bin/bash
#SBATCH --mail-user=tom.willems@unibe.ch
#SBATCH --mail-type=start,fail,end,ARRAY_TASKS
#SBATCH --job-name=2ndlevel_deco_unconsc
#
#
#SBATCH --time=02:10:01
#SBATCH --mem-per-cpu=8G
#SBATCH --tmp=2048
#
#  Partition
#SBATCH --partition=bdw

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

module load MATLAB/2021b

INPUT=$1
STUDY="$STORAGE/$INPUT"
BIDS_DIR="$STUDY/bids"

#how to submit this script:
# sbatch --array=1 second_level.sh hipp/wb


subject=$( sed -n -E "$((${SLURM_ARRAY_TASK_ID} + 1))s/sub-(\S*)\>.*/\1/gp" ${BIDS_DIR}/participants.tsv )

#
#srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nodisplay -singleCompThread -batch "cd('/storage/homefs/tw18a205/TOAM/fmri/2nd_level'); second_level_multiSession($subject,1,'retrieval','ret'); exit;"


# Define an array of input strings
#input_strings=("Retrieval" "DelayedRetrieval" "Recognition")

# Loop over the input strings and run each one
# for session in "${input_strings[@]}"
# do
#     # Run MATLAB with the current input string
#     #matlab -nodisplay -nosplash -nodesktop -r "second_level_decodingMaps('$input'); exit"
#     srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nodisplay -singleCompThread -batch "cd('/storage/homefs/tw18a205/TOAM/fmri/2nd_level'); second_level_decodingMaps($subject,$session); exit;"

#     # Wait for 5 seconds before running the next job
#     sleep 5
# done

for session in {"Retrieval","DelayedRetrieval"}
do 
   srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nodisplay -singleCompThread -batch "cd('/storage/homefs/tw18a205/TOAM/fmri/2nd_level'); second_level_decodingMapsModels($subject,'$session'); exit;"
done

# for correl_name in {"wpt","ret4","ret4unconsc","ret4unconscAboveChance","numConscAnswersRet4","ret3","ret3unconsc","ret3unconscAboveChance","numConscAnswersRet3"}
# do 
#    srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nodisplay -singleCompThread -batch "cd('/storage/homefs/tw18a205/TOAM/fmri/2nd_level'); second_level_correlation($subject,'$correl_name','maxTrials'); exit;"
# done

#srun -N 1 -n 1 -c 1 --cpu_bind=cores matlab -nodisplay -singleCompThread -batch "cd('/storage/homefs/tw18a205/TOAM/fmri/2nd_level'); second_level_delret_subseq($subject,'ret4','maxTrials'); exit;"

