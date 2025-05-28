#!/bin/bash

# create a folder called my_images or something similar in your homedirectory
# then run this script to download the fmriprep image
# don't try to hand it in as a job, via s_run or sbatch, just run it in the terminal via: bash get_fmriprep.sh

# possible problems: singularity doesn't find directories to store temporary files.
# is workspace loaded? what does "echo $WORKSPACE" return?
# if it returns nothing, load workspace via "module load Workspace"
# another way is to export the environment variable of the tmp directory to the basic scratch folder via:
export SINGULARITY_TMPDIR=/scratch/local
singularity build ~/my_images/fmriprep-newest.simg docker://nipreps/fmriprep:23.2.0
singularity build ~/my_images/mriqc_latest.sif docker://nipreps/mriqc:latest 
