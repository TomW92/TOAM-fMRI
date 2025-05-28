function first_level(subj, sess, sessionName, mniSpace ,myContrast, myBehavior, smooth, test_cases)
% first_level - Calculates the first level analysis for each subject
%
% Syntax: first_level( subj, sess, behavior_name, myClassic,myContrast )
%
% INPUT:    - subject: string of numbers or number of subject (either 1-n or
%               participant name
%           - sess: number or string of number, session number 1-3
%           - myBehavior: corresponding bahvior name
%           - myContrast: name of the Contrast used 
%
% Other m-files required: contrast_*.m, setup_config
% Subfunctions: spm_jobman
%
% See also: SPM
% Authors: Mirko Bristle, Tom Willems

if nargin <1
   subj=60601;   %input('Subject: ');
   sess=1;
   sessionName='forgetexp';
   mniSpace = '1';
   myContrast= 'enc_subsequent';
   myBehavior= 'enc_subsequent';
end
if nargin <8
    test_cases=0;
end
if ~isnumeric(sess)
    sess = str2double(sess);
end
if ~isnumeric(mniSpace)
	mniSpace = str2double(mniSpace);
end

%% SETUP
myMix='';
myMovementPrefix = 'rp_';

if mniSpace
    if smooth
        myPrefix = [myMix,'s6wr'];
    else
        myPrefix = [myMix,'wr'];
    end
    mniString = 'mni_';
    myMask = 'FLMask';
else
    if smooth
        myPrefix = [myMix,'s1r'];
    else
        myPrefix = [myMix,'r'];
    end
    mniString = '_native_';
    myMask = 'nativeFLMask';
end

% add functions
p=pwd;
idcs=strfind(p,filesep);
addpath(fullfile(p(1:idcs(end)-1),'functions'))
addpath('contrasts')

%define contrast files.
if class(myContrast) == 'cell'
    myLSSTrial = myContrast{2};
    myContrast = myContrast{1};
else 
    myLSSTrial = cell(1);
end
specify_contrasts = str2func(['contrasts_',myContrast]);

%Setup the scripts according to setup_config in the same script folder
config = setup_config(subj);
config.c_session=config.sessions{sess};

test_string='PPI';
if test_cases == 1
   %myMovementPrefix = 'PCA_60601';
   test_string = 'fmriPrepAROMA';
   %config.data.workspace = fullfile('/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/',config.seq.seq{1});
end
config.l1.outdirs = {'firstLevel'};

nbatch=0;
clear matlabbatch;

%read behavior (output from R script)
behavior = readtable([config.data.behav.maxTrials.(myBehavior)],'FileType','delimitedtext'); % the sessions data is stored)
assert(any(cell2mat(arrayfun(@(x) x==config.id,unique(behavior.id),'UniformOutput',false))),'subject was discarded in behavior analysis')

%read scans
dataPre=myPrefix((length(myMix)+1):end); % Substring if mixedAnalysis
files=cellstr(spm_select('ExtFPList', fullfile(config.data.deriv.spmMB.ses(sess).func), ['^',dataPre,'.*',sessionName,'.*\.nii$']));
%files = cellstr(spm_select('ExtFPList', fullfile(config.data.deriv.fmriprep.ses(sess).func), ['^sub.*smoothAROMAnonaggr_bold.nii$'])); %for subsequentdelret with fmriprep files

disp("herecom the files")
disp(files{1})

if myLSSTrial{1}
    behavior = [{behavior},{myLSSTrial}];
end

%get contrast depended variables
[conds, contrasts, convecs, onsets, durations] = specify_contrasts(behavior);
disp("conds specified: ")
disp(conds)
config.l1.conds = conds;
disp("contrasts")
disp(contrasts)
config.l1.contrasts=contrasts;
disp("convecs")
disp(convecs)   
config.l1.convecs=convecs;
disp("onsets")
disp(onsets)
config.l1.onsets=onsets;
disp("durations")
disp(durations)
config.l1.duration = durations;
%config.l1.pmod = struct('name', {}, 'param', {}, 'poly', {});
%config.l1.pmod = pmods;

%% Model declaration (SPM.mat)

% nbatch = nbatch + 1;
%matlabbatch{nbatch}.spm.spatial.smooth.data   = files;
%matlabbatch{nbatch}.spm.spatial.smooth.fwhm   = [6,6,6];
%matlabbatch{nbatch}.spm.spatial.smooth.dtype  = 0;
%matlabbatch{nbatch}.spm.spatial.smooth.im     = 0;
%matlabbatch{nbatch}.spm.spatial.smooth.prefix = ['s6'];
%files = cellstr(spm_file(files,'prefix','s6'));

if contains(myContrast,'_LSS')
    outdir=fullfile(config.data.bold,config.strID,'LSS_GLMs_all',myContrast,['l1_',mniString,myContrast,'_',myLSSTrial{1}(1:9),'_', sessionName,'_',myPrefix,test_string]);
else
    outdir=fullfile(config.data.workspace,'1st_level',config.strID,['l1_',mniString,myContrast,'_', sessionName,'_',myPrefix,test_string,]);
end

if isfolder(outdir)
    disp("current_dir:")
    disp(outdir)
    disp("deleting prior attempt")
    rmdir(outdir,'s')
end
mkdir(outdir)

nbatch = nbatch + 1;
matlabbatch{nbatch}.spm.stats.fmri_spec.dir = cellstr(outdir);
matlabbatch{nbatch}.spm.stats.fmri_spec.sess.scans = files;
matlabbatch{nbatch}.spm.stats.fmri_spec.sess.multi_reg = cellstr(spm_select('FPList', config.data.deriv.spmMB.ses(sess).func,['^',myMovementPrefix,'.*',sessionName,'.*.txt']));

% Conditions --------------------------------------------------------------
% write trialonset for each subj, session, cond
for k = 1:length(config.l1.conds)

    %set parameters
    matlabbatch{nbatch}.spm.stats.fmri_spec.sess.cond(k).name = config.l1.conds{k};
    matlabbatch{nbatch}.spm.stats.fmri_spec.sess.cond(k).onset = config.l1.onsets{k};
    matlabbatch{nbatch}.spm.stats.fmri_spec.sess.cond(k).duration = config.l1.duration{k}; %set 0 since event related design, now not annymore, TB said its jjust more information added even if it is just 2.5 seconds or sth comparably short
    matlabbatch{nbatch}.spm.stats.fmri_spec.sess.cond(k).tmod = 0;
    matlabbatch{nbatch}.spm.stats.fmri_spec.sess.cond(k).pmod = struct('name', {}, 'param', {}, 'poly', {}); %config.l1.pmod{k};
    matlabbatch{nbatch}.spm.stats.fmri_spec.sess.cond(k).orth = 0;

end

% Parameters --------------------------------------------------------------
config.l1.nslices =  config.params.nslices; % this must be the same as defined in preprocessing: slice-time correction
config.l1.TR = config.params.TR;   % repetition time.
matlabbatch{nbatch}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{nbatch}.spm.stats.fmri_spec.timing.RT = config.params.TR;   % repetition time.
matlabbatch{nbatch}.spm.stats.fmri_spec.timing.fmri_t = config.l1.nslices; % this must be the same as defined in preprocessing: slice-time correction, nslices
matlabbatch{nbatch}.spm.stats.fmri_spec.timing.fmri_t0 = config.l1.nslices/2; % this must be the same as defined in preprocessing: slice-time correction
matlabbatch{nbatch}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{nbatch}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{nbatch}.spm.stats.fmri_spec.mthresh = -Inf; %0.8; %0.05;
matlabbatch{nbatch}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{nbatch}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});

%matlabbatch{nbatch}.spm.stats.fmri_spec.bases.fir.length = 12.26;
%matlabbatch{nbatch}.spm.stats.fmri_spec.bases.fir.order = 5 ;
matlabbatch{nbatch}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];

matlabbatch{nbatch}.spm.stats.fmri_spec.volt = 1;
matlabbatch{nbatch}.spm.stats.fmri_spec.global = 'None';

%Use an explicit mask
explicitMask = config.data.deriv.spmMB.ses(sess).(myMask);
assert(~isempty(explicitMask),'NO mask image found')
matlabbatch{nbatch}.spm.stats.fmri_spec.mask = explicitMask; %session indicator doesnt matter because its MNI space.
matlabbatch{nbatch}.spm.stats.fmri_spec.cvi = 'AR(1)';

%% Estimate
% spm12 patched by: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=spm;2c0ae193.1611
nbatch=nbatch+1;
spmfile=cellstr(fullfile(outdir,'SPM.mat'));
matlabbatch{nbatch}.spm.stats.fmri_est.spmmat = spmfile;
matlabbatch{nbatch}.spm.stats.fmri_est.write_residuals = 0; 
matlabbatch{nbatch}.spm.stats.fmri_est.method.Classical = 1;

%% Make contrasts
if length(config.l1.contrasts)>1
    disp("you should again see this")
    nbatch=nbatch+1;
    
    matlabbatch{nbatch}.spm.stats.con.spmmat = spmfile;
    matlabbatch{nbatch}.spm.stats.con.delete = 0;
    
    %for each contrast
    for k = 1:length(config.l1.contrasts)
        matlabbatch{nbatch}.spm.stats.con.consess{k}.tcon.name = config.l1.contrasts{k};
        matlabbatch{nbatch}.spm.stats.con.consess{k}.tcon.convec = config.l1.convecs{k};
        matlabbatch{nbatch}.spm.stats.con.consess{k}.tcon.sessrep = 'none';    
    end
end

%% save and run
if contains(myContrast,'_LSS')  
    save(fullfile(outdir,['1st_level',config.strID,'_l1_',mniString,myContrast,'_',myLSSTrial{1}(1:9),'_',sessionName,'_',myPrefix,test_string]),'matlabbatch')
else
    save(fullfile(outdir,['1st_level',config.strID,'_l1_',mniString,myContrast,'_', sessionName,'_',myPrefix,test_string]),'matlabbatch')
end 

%Run (only now will the actual calculations happen)
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

disp("model estimated and contrasts computed!")
end
