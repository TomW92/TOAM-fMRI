function second_level_delret_subseq(subj, myCorr, onset_type )
% second_level_acc_cat_cond - calculates teh second level analysis 
%  - GLM and behavioral correlations
%
% Inputs:
% subj - string of subject number
%
% Other m-files required: run2ndLvlTcon.m, setup_config.m, 
% Subfunctions: none
% MAT-files required: SPM.mat of first level
%
% See also: first_level.m, 
% Author: Mirko Bristle, Tom Willems

if nargin <1
   subj=60601;   %input('Subject: ');
   myCorr= 'ret';
   onset_type = 'maxTrials';
end
if nargin <3
   onset_type = 'maxTrials';
end


clear matlabbatch;

%% Setup
analysis_name=[myCorr,'fmriPrepAROMA',onset_type];%,
%L1_p='l1_mni_delret_subsequent_forgetexp_s6wrwithfMRIPREPfiles';
L1_p ='l1_mni_delret_subsequent_forgetexp_s6wrfmriPrepAROMA';

p=pwd;
idcs=strfind(p,'/');
addpath(fullfile(p(1:idcs(end)-1),'functions'))

%Setup the scripts according to setup_config in the same script folder
config = setup_config(subj);
disp("currently working with correlation: ")
disp(myCorr)

subjList=config.seq.participants.participant_id; %unique(behavior.id);


%get correlations
if strcmp(myCorr,'wpt')
    behavior = readtable(fullfile(config.data.behav.(onset_type).wpt));
    subjIdx  = cellfun(@(x) ismember(str2double(x(5:end)),behavior.id),subjList)>0;
    subjList=subjList(subjIdx);
    behavior.accuracy = behavior.mean_all;
elseif strcmp(myCorr,'ret4')
    behavior = readtable(fullfile(config.data.behav.(onset_type).allRetsAcc),"VariableNamingRule","preserve");
    behavior = behavior(behavior.retrieval==4,:);
    behavior.accuracy = behavior.Acc;
elseif strcmp(myCorr,'ret4unconsc')
    behavior = readtable(fullfile(config.data.behav.(onset_type).unconscRet),"VariableNamingRule","preserve");
    behavior.accuracy = behavior.ret4;
elseif strcmp(myCorr,'ret4unconscAboveChance')
    behavior = readtable(fullfile(config.data.behav.(onset_type).unconscRet),"VariableNamingRule","preserve");
    behavior.accuracy = behavior.ret4;
    behavior = behavior(behavior.accuracy>0.5,:);
elseif strcmp(myCorr,'numConscAnswersRet4')
    behavior = readtable(fullfile(config.data.behav.(onset_type).allRetsAccConsc),"VariableNamingRule","preserve");
    behavior = behavior(behavior.retrieval==4,:);
    behavior = behavior(behavior.consc==4,:);
    behavior.accuracy = behavior.freq;
elseif strcmp(myCorr,'numConscAnswersRet3')
    behavior = readtable(fullfile(config.data.behav.(onset_type).allRetsAccConsc),"VariableNamingRule","preserve");
    behavior = behavior(behavior.retrieval==3,:);
    behavior = behavior(behavior.consc==4,:);
    behavior.accuracy = behavior.freq;
elseif strcmp(myCorr,'ret3unconsc')
    behavior = readtable(fullfile(config.data.behav.(onset_type).unconscRet),"VariableNamingRule","preserve");
    behavior.accuracy = behavior.ret3;
elseif strcmp(myCorr,'ret3unconscAboveChance')
    behavior = readtable(fullfile(config.data.behav.(onset_type).unconscRet),"VariableNamingRule","preserve");
    behavior.accuracy = behavior.ret3;
    behavior = behavior(behavior.accuracy>0.5,:);
end

%get subjList with L1 analysis
%idcs=strfind(config.data.path,filesep); SMART CODE; NEED THIS LATER!

subjIdx=arrayfun(@(x) isfile(fullfile(config.data.bold,x,L1_p,'con_0001.nii')),subjList)>0;
subjList=subjList(subjIdx);
subjPathList = strcat(fullfile(config.data.bold),filesep,subjList);
nsubj=size(subjList,1);

subjListNUM = cellfun(@(x) str2double(x(1,5:end)),subjList);
behavior=behavior(ismember(behavior.id,subjListNUM),:);

spm_p=fullfile(subjPathList(1,:),L1_p,'SPM.mat');
load(spm_p{1})
assert(nsubj>0,'No Subject available for second-level analysis. Check session')

% Creation of directory for the rfx
rfxDir = fullfile(config.data.bold,'RFX_Correlations',['RFX_',L1_p,'_s6',analysis_name]);
% Creation of directory for the rfx
if exist(rfxDir,'dir')
    rmdir(rfxDir,'s')
end
mkdir(rfxDir)
cd(rfxDir)

%% Save the used batch in the folder
% Con Images: {'intact_correct';'intact_incorrect';'broken_correct';'broken_incorrect';'baseline';'intact > broken';'intact > baseline';'broken > intact';'broken > baseline';'baseline > intact';'baseline > broken';'correct > incorrect';'correct > baseline';'incorrect > correct';'incorrect > baseline';'baseline > correct';'baseline > incorrect'}

cov_names={'accuracy'};
analysis= 1:6; % whats with the 233rd conttrast??
for icon=analysis
    % Perform for each Cov seperatly
    for cov_name=cov_names
        disp([SPM.xCon(icon).name,cov_name])

        %Create the folder of that contrast
        analysis_name = ['CORREL_',strrep(SPM.xCon(icon).name,' > ','_'),...
            '_X_',cov_name{1}];
        mkdir(analysis_name);
        %contrast = SPM.xCon(i_con).Vcon.fname;
        
        files = strcat(subjPathList,filesep,fullfile(L1_p,sprintf('con_%04d.nii',icon)),',1');
        
        mbatch.conName    =SPM.xCon(icon).name;
        mbatch.outDir     =cellstr(fullfile(rfxDir,analysis_name));
        mbatch.scans      =cellstr(files);
        mbatch.covariates = {char(cov_name)};
        mbatch.covarvecs  = {behavior.(cov_name{1})};
        mbatch.contrasts  = { ['Contrast ',SPM.xCon(icon).name],...
                            ['Contrast inverted ',SPM.xCon(icon).name],...
                            ['Correl: ',SPM.xCon(icon).name],...
                            ['AntiCorrel: ',SPM.xCon(icon).name]};
        mbatch.convecs    = { [1 0],...
                            [-1 0],...
                            [0 1],...
                            [0 -1]};
        run2ndLvlTcon(mbatch);
    end
end
end
