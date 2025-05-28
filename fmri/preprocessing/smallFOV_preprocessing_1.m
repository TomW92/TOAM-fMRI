function preprocessing_1(subj)
% preprocessing_1 -  Does first step of preprocessing:
% DOES: slice-time, realign&unwarp for both mri sequences of the TOAM study
% Inputs:
% subj - string of subject name
% Other m-files required: setup_config.m
% Subfunctions: none
% MAT-files required: none
%
% Author: Mirko Bristle, Tom Willems

% IF no argument is supplied, the subject is asked
if nargin<1
    subj = 60601;
end
%%
p=pwd;
idcs=strfind(p,'/');
addpath(fullfile(p(1:idcs(end)-1),'functions'))

%Setup the scripts according to setup_config in the same script folder
config = setup_config(subj);
nbatch=0;

%% Slice Timing Correction
%--------------------------------------------------------------------------
%config.session contains 4 cell arrays with all the nifti functional nifti files for each of the 4 sessions
%config.session_signs contains just the information about wether the
%session was recorded on day one or two, in this study we had 3
%runs/sessions on day 1 and 1 run/session on day two.
for i=1:length(config.session_signs) 
    config.current_session = config.session{i}; 
    nbatch=nbatch+1;
    matlabbatch{nbatch}.spm = slice_time_correction(config);
end

%% Realign & Unwarp
%--------------------------------------------------------------------------
for i=1:length(config.session_signs)
    config.data.deriv.a_files={cellstr(spm_file(config.session{i},'prefix','a'))};
    nbatch=nbatch+1;
    matlabbatch{nbatch}.spm = realign(config.data.deriv.a_files);
end

%% Save the script used
save(fullfile(config.data.bold,'1_slicerealign'),'matlabbatch')
%Run 
spm_jobman('run',matlabbatch);

disp(['... pre-processing part 1 done on subj ', num2str(subj)])
return
