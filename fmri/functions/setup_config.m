function [config] = setup_config(id,rstorage)
%setup_config - Prepares most other scripts with setup. 
%  sets up Datapath for the analysis
%  if running on mac (assuming that this mac has connection to the
%  research storage and able to login with smb connection)
% Syntax: [config] = setup_config(id,qnap)
%
% Inputs:
% id - participant number 
% rstorage - boolean if connect to research storage,
%
% Outputs:
% config - general inforamtion of subject path setup.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Mirko Bristle
% Adapted by: Tom Willems
%
% to work with dataset in the BIDS format
% should now work with other BIDS format datasets.
assert(nargin>=1, 'Wrong number of input arguments: no Subject')
if nargin == 1 && ~exist('qnap','var')
   rstorage=false; 
end

%Print out which Subjects will be
config = struct();
config.sessions  = cellstr(['ses-1'; 'ses-2']); %why does this throw an error if strings aren't the same length
% cast subj to be numeric
if ~isnumeric(id)
    config.id = str2double(id);
else
    config.id = id;
end
config.strID = ['sub-',num2str(id)];

%Get who's calling me (script)
st = dbstack('-completenames');

%Useful for UBELIX debugging
disp(['Script: ',st(end).file])


%search and find sequence for subject that was provided
shouldbeOne = 0;
isUbelix = 0;   
search_counter = 0;
sequences= {'hipp','wb'};
% system depending variables
while (shouldbeOne < 1 && search_counter<length(sequences))
    for seq = sequences
        if ismac 
            addpath('/Users/tomwillems/MATLAB/spm12')
            if config.id==6601 && ~rstorage
                config.data.path ='/Users/tomwillems/Downloads/fmri/sub-60601';
            else
                fmriPath = '/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/';
                config.data.workspace=fullfile('/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/',seq{1},'bids/derivatives/spmMB/');
                data_path = fullfile(fmriPath, seq{1});
                if ~isfolder(data_path)
                    system('open smb://resstore.unibe.ch/psy_wfg_henke/ -gj','-echo');
                end
                config.data.path = fullfile(data_path,'bids');
            end
        elseif ispc
            addpath('C:\Users\mk24i760\Documents\MATLAB\spm12')
            fmriPath = 'X:\s2019_twillems_fMRI_silent_engram\data\fMRI';
            data_path = fullfile(fmriPath, seq{1});
            config.data.path = fullfile(data_path,'bids');

        else %assumes is on UBELIX 
            isUbelix = 1;
            addpath('/storage/homefs/tw18a205/toolboxes/spm12')
            fmriPath = '/storage/research/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/';
            basePath=fullfile(fmriPath, seq{1});
            config.data.path=fullfile(basePath,'bids');
            config.data.workspace = fullfile('/storage/workspaces/psy_memory_wfg_psy/hpc_henke_wfg/s2019_twillems_fMRI_silent_engram/data/fMRI/',seq{1});
        end
        config.seq.participants = tdfread(fullfile(config.data.path,'participants.tsv'));
        config.seq.participants.participant_id = cellstr(config.seq.participants.participant_id);
        
        % check if given ID is actually part of this list, i.e. if you have given a
        % subject name from the correct sequence.
        for i=1:length(config.seq.participants.participant_id)
            if strcmp(config.seq.participants.participant_id{i},config.strID)
                shouldbeOne = 1;                
            end
        end
        try 
            assert(shouldbeOne==1)
            break
        catch
            disp(['Given subject ID is not a part of group: ',seq{1}])
            search_counter = search_counter + 1;
        end
            
    end
end


disp(['Participant: ',num2str(config.id)])
assert(shouldbeOne==1,'Given subject is not part of any group')

if strcmp(seq,'hipp')
    spm1 = 'spm';
    spm2 = 'spmMB';
    
elseif (strcmp(seq,'wb'))
    spm1 = 'spm_TB';
    spm2 = 'spmMB';
end

%get infos from json sidecar from this participant
config.seq.scans = tdfread(char(fullfile(config.data.path,config.strID,config.sessions(1),[char(config.strID),'_',char(config.sessions(1)),'_scans.tsv'])));
config.seq.scans.filename = cellstr(config.seq.scans.filename);
config.seq.seq = seq;

%Initalize spm
global initspm;
if isempty(initspm)
    spm('Defaults','fMRI');
    spm_jobman('initcfg');
    initspm=1;
end

%add files session-wise
for i=1:size(config.sessions,1)
        %anat path
        config.data.ses(i).bidsAnat = fullfile(config.data.path,config.strID,config.sessions(i),'anat');
        config.data.ses(i).bidsAnatFiles = cellstr(spm_select('FPListRec',config.data.ses(i).bidsAnat,'^s.*\.nii.gz'));


        %func_path
        config.data.ses(i).bidsFunc = fullfile(config.data.path,config.strID,config.sessions(i),'func');
        config.data.ses(i).bidsFuncFiles = cellstr(spm_select('FPListRec',config.data.ses(i).bidsFunc,'^s.*\.nii*'));

        %func_sidecar
        config.data.ses(i).bidsFuncJSON   = regexprep(config.data.ses(i).bidsFuncFiles{1},'nii*.*','json');
        config.data.ses(i).bidsFuncParams = jsondecode(fileread(config.data.ses(i).bidsFuncJSON));

        %derivs
        config.data.deriv.folder = fullfile(config.data.path,'derivatives');
 
        %depreceated spm derivatives
        config.data.deriv.spm.ses(i).anat =  fullfile(config.data.deriv.folder,spm1,config.strID,'preproc/anat',config.sessions(i));
        config.data.deriv.spm.ses(i).anatFiles = cellstr(spm_select('FPListRec',config.data.deriv.spm.ses(i).anat,'\.nii'));
        config.data.deriv.spm.ses(i).func =  fullfile(config.data.deriv.folder,spm1,config.strID,'preproc/func',config.sessions(i));
        config.data.deriv.spm.ses(i).funcFiles = cellstr(spm_select('FPListRec',config.data.deriv.spm.ses(i).func,'\.nii'));
        
        %current spm derivatives!
        config.data.bold = fullfile(config.data.deriv.folder,spm2);
        config.data.deriv.spmMB.ses(i).anat = fullfile(config.data.deriv.folder,spm2,config.strID,'preproc',config.sessions(i), 'anat');
        config.data.deriv.spmMB.ses(i).anatFiles = cellstr(spm_select('FPListRec',config.data.deriv.spmMB.ses(i).anat,'\.nii'));
        config.data.deriv.spmMB.ses(i).func = fullfile(config.data.deriv.folder,spm2,config.strID,'preproc',config.sessions(i), 'func');
        config.data.deriv.spmMB.ses(i).funcFiles = cellstr(spm_select('FPListRec',config.data.deriv.spmMB.ses(i).func,'\.nii'));   
        config.data.deriv.spmMB.ses(i).FLMask = cellstr(spm_select('FPListRec',fullfile(config.data.bold,config.strID,'preproc'),'^crossSessionWhiteandGreyMatterMask'));
        config.data.deriv.spmMB.ses(i).nativeFLMask = cellstr(spm_select('FPListRec',config.data.deriv.spmMB.ses(i).anat,'^rfbrainMask')); 

        %fmriprep
        if strcmp(seq,'wb')
            config.data.deriv.fmriprep.ses(i).func = fullfile(config.data.deriv.folder,'fmriprep',config.strID, config.sessions(i),'func');
        end

        %behavioral files
        onset_types = {'maxTrials','corrected'};
        for onset_type = onset_types
            config.data.behav.(onset_type{1}).enc_subsequent = fullfile(fmriPath,['utilities/eventFiles/subs_memory_all/',config.strID,'_subsequent_consc_',onset_type{1},'.csv']);
            config.data.behav.(onset_type{1}).consolidation  = fullfile(fmriPath,['utilities/eventFiles/consolidation/',config.strID,'_consolidataion.csv']);
            config.data.behav.(onset_type{1}).enc_bids       = fullfile(config.data.deriv.folder,['fmriprep/',config.strID,'/ses-1/func/',config.strID,'_ses-1_task-forgetexp_events_',onset_type{1},'.tsv']);
            config.data.behav.(onset_type{1}).ret            = fullfile(fmriPath,['utilities/eventFiles/retrieval_all/',config.strID,'_retrieval_all_events_',onset_type{1},'.csv']);
            config.data.behav.(onset_type{1}).ret_bids       = fullfile(config.data.deriv.folder,['fmriprep/',config.strID,'/ses-1/func/',config.strID,'_ses-1_task-retrieval_events_',onset_type{1},'.tsv']);
            config.data.behav.(onset_type{1}).delret         = fullfile(fmriPath,['utilities/eventFiles/delayed_retrieval_all/',config.strID,'_delayed_retrieval_all_events_',onset_type{1},'.csv']);
            config.data.behav.(onset_type{1}).delret_bids    = fullfile(config.data.deriv.folder,['fmriprep/',config.strID,'/ses-2/func/',config.strID,'_ses-2_task-forgetexp_events_',onset_type{1},'.tsv']);
            config.data.behav.(onset_type{1}).wideSubs       = fullfile(fmriPath,['behav/fMRI_behav_allgroups_wide_',onset_type{1},'.csv']);
            config.data.behav.(onset_type{1}).longSubs       = fullfile(fmriPath,['behav/fMRI_behav_allgroups_long_',onset_type{1},'.csv']);
            config.data.behav.(onset_type{1}).wpt            = fullfile(fmriPath,'behav/WPT/WPT_final_results.csv');
            config.data.behav.(onset_type{1}).unconscRet     = fullfile(config.data.deriv.folder,['behav/mean_unconsc_accuravy_wideFormat_',seq{1},'_',onset_type{1},'.csv']);
            config.data.behav.(onset_type{1}).allRetsAcc     = fullfile(fmriPath,['behav/accuracies_per_sub_retall_',onset_type{1},'.csv']);
            config.data.behav.(onset_type{1}).allRetsAccConsc= fullfile(fmriPath,['behav/accuracies_per_sub_ret_conscall_',onset_type{1},'.csv']);
            config.data.behav.(onset_type{1}).delret_subseq  = fullfile(fmriPath,['utilities/eventFiles/subs_memory_all/',config.strID,'_delret_subsequent_consc_',onset_type{1},'.csv']);
        end
        
        if isUbelix
            config.data.deriv.searchlight.folder = fullfile(config.data.workspace,'/searchlight/',config.strID);
            config.data.deriv.freesurfer.ses(i).mri  = fullfile(config.data.workspace,'/freesurfer/',config.sessions(i), config.strID,'/mri');
            config.data.deriv.freesurfer.ses(i).masks = cellstr(spm_select('FPListRec',config.data.deriv.freesurfer.ses(i).mri,'\.nii'));
        end


end


%assumes that params are the same for all sessions, thus just takes info
%from session 1 JSON sidecar.
if strcmp(seq,'hipp')
    config.params.nslices = length(config.data.ses(1).bidsFuncParams.SliceTiming);
    config.params.TR = config.data.ses(1).bidsFuncParams.RepetitionTime; % repetition time  - dicom_tag (0018,0080)
    config.params.TA = double(config.params.TR-(config.params.TR/config.params.nslices));
    config.params.slice_order = 1000*config.data.ses(1).bidsFuncParams.SliceTiming;
    config.params.ref_slice = config.params.slice_order(config.params.nslices/2);      %reference TIME OF middle slice because above is the real time
    config.params.res = 0.8;
elseif strcmp(seq,'wb')
    config.params.nslices = length(config.data.ses(1).bidsFuncParams.SliceTiming);
    config.params.TR = config.data.ses(1).bidsFuncParams.RepetitionTime; % repetition time  - dicom_tag (0018,0080)
    config.params.TA = config.params.TR-(config.params.TR/config.params.nslices);
    config.params.slice_order = 1000*config.data.ses(1).bidsFuncParams.SliceTiming;
    config.params.ref_slice = config.params.slice_order(config.params.nslices/2);    %reference to middle slice  #before INT8 here and reference slice, but WHY
    config.params.res = 1.1;
end

%to loop through sessions efficiently:
config.session_names = {'forgetexp','consolidation','retrieval','forgetexp'};
config.session_signs = {1,1,1,2};

% get all files func files for efficient preprocessing
for i=1:length(config.session_signs)
    config.session{i}=cellstr(spm_select('ExtFPList', config.data.deriv.spmMB.ses(config.session_signs{i}).func, ['^sub.*' config.session_names{i} '.*\.nii$']));
end
return