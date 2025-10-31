# TOAM - Trajectory of a Memory Trace

## Contents
This repository contains code and data files required to reproduce the results of the study entitled "Neural Traces of Forgotten Memories Persist in Humans and are Behaviorally Relevant" by Tom Willems, Konstantinos Zervas, Luzius Brogli, Finn Rabe, Andrea Federspiel, and Katharina Henke. A preprint of the study can be found on [bioRxiv](https://doi.org/10.1101/2025.06.02.656652). The fMRI data are available on OpenNeuro:
- [small FOV fMRI group](https://doi.org/10.18112/openneuro.ds006265.v1.0.0)
- [whole-brain fMRI group](https://doi.org/10.18112/openneuro.ds006266.v1.0.0)

___

## Behavior

The analysis scripts can be found in the subfolder `/behavior`.
The behavioral data can be found in the subfolder `/data/behavioral_data`.

___

## fMRI
### Preprocessing
The preprocessing of the fMRI data was done with fMRIprep (whole-brain fMRI sequence) and SPM12 (small FOV fMRI sequence). The preprocessing scripts can be found in the subfolder `/fMRI/preprocessing`.

### First Level
The first-level analysis can be executed with the `first_level.m` function for a specific session of a specific subject.
It dependents on a successful preprocessing of the fMRI and behavioral data. The contrast must be specified in a `fMRI/1st_level/contrasts/contrast_CONTRAST-NAME.m` file. This file is called based on the contrast name entered to the function.
The first-level analysis scripts can be found in the subfolder `/fMRI/1st_level`.

### Second Level
There are 3 second-level analyses: `second_level_multiSession.m`, which computes the second-level analysis for all sessions across all subjects, `second_level_delret_subseq`, which computes the second-level analysis for a specific contrast used to analyse subsequent recovery of forgotten memories, and `second_level_correlation`, which computes second-level brainbehavior correlations for a specific contrast.
The second-level analysis scripts can be found in the subfolder `/fMRI/2nd_level`.

### RSA
The RSA analysis is found in the subfolder `/fMRI/mvpa/rsa/ers`. The main analysis script is `encoding_retrieval_rsa.m`, which computes the RSA for all encoding and retrieval conditions.

### Plotting
The scripts that create all the main and supplementary figures can be found in `/fMRI/plotting`.

### Functions
Utility functions used in the fMRI analysis can be found in the subfolder `/fMRI/functions`.

## Paper
The supplementary material and all Figures can be found in `/paper/src`.

## License
Licensed under Creative Commons Attribution 4.0 International (CC BY 4.0)
https://creativecommons.org/licenses/by/4.0/

## Acknowledgements
This work was supported by Sitem-Insel Support Funds SISF 2019 to K. Henke and by the SNFS Advanced Grant TMAG-1_209374 to K. Henke. The authors thank Mirco Bristle for his support during data analysis. Calculations were performed on UBELIX (https://www.id.unibe.ch/hpc), the HPC cluster at the University of Bern.

## Contact
tom.willems@unibe.ch
