function multi_session_first_level( subj, flName )
% first_level - Calculates the first level analysis for each subject
%
% INPUT:    - subject: string of numbers or number of subject (either 1-n or
%               participant name
%           - flName: Name of the first level folder. 
%
% Other m-files required: contrast_multi_session, setup_config
% Subfunctions: spm_jobman
%
% See also: SPM
% Authors: Mirko Bristle, Tom Willems

if nargin <2
   subj=60601;
   flName='pleaserenamethisfirstlevelfolder';
end

%% SETUP
% add functions and contrasts
p=pwd;
idcs=strfind(p,'/');
addpath(fullfile(p(1:idcs(end)-1),'functions'))
addpath('contrasts')

%Setup the scripts according to setup_config in the same script folder
config = setup_config(subj);
myPrefix = 's6wr';
myMovementPrefix = 'rp_';

%config.c_session=config.sessions{sess};
config.l1.outdirs = {'firstLevel'};
config.l1.contrasts={};
config.l1.conds={};
config.l1.convecs={};
config.l1.onsets={};

nbatch=0;
clear matlabbatch;

%% Model declaration (SPM.mat)
nbatch=nbatch+1;
outdir=fullfile(config.data.bold,config.strID,['l1_', flName ,'_',myPrefix]);
matlabbatch{nbatch}.spm.stats.fmri_spec.dir = cellstr(outdir);

%Load struct from multi-session contrast file.
multiSessionContrasts = contrast_multi_session(config);

for sess_ind=1:length(config.session_signs)
    %session specific names now:
    behaviorName = config.session_names{sess_ind};
    
    %% matlabbatch first level
    matlabbatch{nbatch}.spm.stats.fmri_spec.sess(sess_ind).scans = multiSessionContrasts.ses(sess_ind).files;
    % write trialonset for each subj, session, cond
    % Conditions --------------------------------------------------------------
    for k = 1:length(multiSessionContrasts.ses(sess_ind).conds)
        %set parameters
        matlabbatch{nbatch}.spm.stats.fmri_spec.sess(sess_ind).cond(k).name = multiSessionContrasts.ses(sess_ind).conds{k};
        matlabbatch{nbatch}.spm.stats.fmri_spec.sess(sess_ind).cond(k).onset = multiSessionContrasts.ses(sess_ind).onsets{k};
        matlabbatch{nbatch}.spm.stats.fmri_spec.sess(sess_ind).cond(k).duration = 0; %multiSessionContrasts.ses(sess_ind).durations{k}; %set 0 since event related design; OR NOT??? TB says you can still add seconds, AF says no way, implicitly differnt under the hood because no longer assumes eventRelated design
        matlabbatch{nbatch}.spm.stats.fmri_spec.sess(sess_ind).cond(k).tmod = 0;
        matlabbatch{nbatch}.spm.stats.fmri_spec.sess(sess_ind).cond(k).pmod = struct('name', {}, 'param', {}, 'poly', {}); % another option to add parametric modulators such as RT to conditions
    end
    
    matlabbatch{nbatch}.spm.stats.fmri_spec.sess(sess_ind).multi = {''}; %here you can put input a mat file with conditions.
    matlabbatch{nbatch}.spm.stats.fmri_spec.sess(sess_ind).regress = struct('name', {}, 'val', {});
    matlabbatch{nbatch}.spm.stats.fmri_spec.sess(sess_ind).multi_reg = cellstr(spm_select('FPList', config.data.deriv.spmMB.ses(config.session_signs{sess_ind}).func,[myMovementPrefix,'.*',behaviorName,'.*.txt']));
    matlabbatch{nbatch}.spm.stats.fmri_spec.sess(sess_ind).hpf = 128; %okay, esp. for eventrelated designs, might be a bit to high for longer, complex WM, block designs.
end

% Parameters --------------------------------------------------------------
matlabbatch{nbatch}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{nbatch}.spm.stats.fmri_spec.timing.RT = config.params.TR;   % repetition time.
matlabbatch{nbatch}.spm.stats.fmri_spec.timing.fmri_t = config.params.nslices; % this must be the same as defined in preprocessing: slice-time correction, nslices
matlabbatch{nbatch}.spm.stats.fmri_spec.timing.fmri_t0 = config.params.nslices/2; % this must be the same as defined in preprocessing: slice-time correction
matlabbatch{nbatch}.spm.stats.fmri_spec.mthresh =  -Inf; %-Inf; %0.05; %0.8
matlabbatch{nbatch}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{nbatch}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{nbatch}.spm.stats.fmri_spec.volt = 1;
matlabbatch{nbatch}.spm.stats.fmri_spec.global = 'None';

assert(~isempty(config.data.deriv.spmMB.ses(1).FLMask),'NO mask image found')
matlabbatch{nbatch}.spm.stats.fmri_spec.mask = config.data.deriv.spmMB.ses(1).FLMask; %session indicator doesnt matter because its MNI space.
matlabbatch{nbatch}.spm.stats.fmri_spec.cvi = 'AR(1)';

%% Estimate
% spm12 patched by: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=spm;2c0ae193.1611
nbatch=nbatch+1;
spmfile=cellstr(fullfile(outdir,'SPM.mat'));
matlabbatch{nbatch}.spm.stats.fmri_est.spmmat = spmfile;
matlabbatch{nbatch}.spm.stats.fmri_est.method.Classical = 1;

%% Make contrasts
nbatch=nbatch+1;
matlabbatch{nbatch}.spm.stats.con.spmmat = spmfile;
matlabbatch{nbatch}.spm.stats.con.delete = 1;

%for each contrast
for k = 1:length(multiSessionContrasts.fullContrasts)
    matlabbatch{nbatch}.spm.stats.con.consess{k}.tcon.name = multiSessionContrasts.fullContrasts{k};
    matlabbatch{nbatch}.spm.stats.con.consess{k}.tcon.convec = multiSessionContrasts.fullConvecs{k};
    matlabbatch{nbatch}.spm.stats.con.consess{k}.tcon.sessrep = 'none';
end

%% save and run
if isfolder(outdir)
    rmdir(outdir,'s')
end
mkdir(outdir)
save(fullfile(outdir,'/1st_level'),'matlabbatch')
disp(outdir)

%Run (only now will the actual calculations happen)
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

disp("model estimated and contrasts computed!")
end
