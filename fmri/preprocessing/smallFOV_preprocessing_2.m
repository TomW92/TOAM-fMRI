function preprocessing_2(subj, fwhm)
% preprocessing_1s -  Does second step of preprocessing:
% DOES: coregister, normalisation functional and t1, smoothing
%
% Syntax: [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
% subj - string of subject name
% classic - deprecated boolean if run with or without fieldmap/dartel
%
% Other m-files required: setup_config.m
% Subfunctions: none
% MAT-files required: none
%
% Author: Mirko Bristle, Tom Willems

% IF no argument is supplied, the subject is asked
if nargin<1
    subj = 60601; 
    fwhm = 1;
end
p=pwd;
idcs=strfind(p,'/');
addpath(fullfile(p(1:idcs(end)-1),'functions'))
if ~isnumeric(subj)
    subj = str2double(subj);
end

%get subject files.
config = setup_config(subj);

%select functional files
config.data.ra_files=arrayfun(@(x){cellstr(spm_select('ExtFPList',...
    fullfile(config.data.deriv.spmMB.ses(config.session_signs{x}).func),['^ras.*' config.session_names{x} '.*.nii$']))},1:length(config.session_signs),'UniformOutput', 1);
mean_images = arrayfun(@(x){cellstr(spm_select('ExtFPList',...
    fullfile(config.data.deriv.spmMB.ses(config.session_signs{x}).func),['^mean.*' config.session_names{x} '.*.nii$']))},1:length(config.session_signs),'UniformOutput', 1);
assert(~isempty(config.data.ra_files));

%select anatomical files
biasCorrected   = arrayfun(@(x){cellstr(spm_select('FPList', config.data.deriv.spmMB.ses(x).anat, '^m.*T1w.nii$'))},1:length(config.sessions),'UniformOutput',1);
skullStrippedt1 = arrayfun(@(x){cellstr(spm_select('FPList', config.data.deriv.spmMB.ses(x).anat, '^skullStripped.*T1w.nii$'))},1:length(config.sessions),'UniformOutput',1);
brainMasks =  arrayfun(@(x){cellstr(spm_select('FPList', config.data.deriv.spmMB.ses(x).anat, '^brainMask.*T1w.nii$'))},1:length(config.sessions),'UniformOutput',1);
c1 = arrayfun(@(x){cellstr(spm_select('FPList', config.data.deriv.spmMB.ses(x).anat, '^c1s.*T1w.nii$'))},1:length(config.sessions),'UniformOutput',1);

clear matlabbatch
nbatch=0;

%% Coregister the mean func image of each session to the biasCorrected image of that very session,
%% apply to other functional files of that session NOT THE SKULLSTRIPPED
% ONE; THAT FAILED with coregistration for the hippocampal slices. 
for i=1:length(config.session_signs)
    ind=config.session_signs{i};
    nbatch = nbatch + 1;

    matlabbatch{nbatch}.spm.spatial.coreg.estimate.ref = cellstr(biasCorrected{ind});
    matlabbatch{nbatch}.spm.spatial.coreg.estimate.source = cellstr(mean_images{i}); % config.data.ra_files{i}(1); <- mit dem ersten func, war schlechter bei mir.
    % meanabold file is created by preprocessing_1.m
    assert(~isempty(matlabbatch{nbatch}.spm.spatial.coreg.estimate.ref),...
        'Can''t find meana.*.nii, was preprocessing_1 succsessfull?')

    matlabbatch{nbatch}.spm.spatial.coreg.estimate.other             = config.data.ra_files{i}; %%config.data.ra_files{i}(2:end); <- wenn erstes func verwendet wird
    matlabbatch{nbatch}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{nbatch}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
    matlabbatch{nbatch}.spm.spatial.coreg.estimate.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{nbatch}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];
end

%reslice anatomical masks to the mean image that was just coregistered into
%their space in the code above. those masks are then also normalized and
%used for the estimation of the first level 
for i=1:length(config.sessions)
    nbatch = nbatch + 1;
    if i == 2
        meanimg = mean_images{4};
    else
        meanimg = mean_images{i};
    end
    matlabbatch{nbatch}.spm.spatial.coreg.write.ref = meanimg;
    matlabbatch{nbatch}.spm.spatial.coreg.write.source = [brainMasks{i} c1{i}]';
    matlabbatch{nbatch}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{nbatch}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{nbatch}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{nbatch}.spm.spatial.coreg.write.roptions.prefix = 'rf';

    %% FOR FIRST LEVEL THAT MASKING THAT IS NOT MNI SPACE!!! Binarize brain mask (grey+white matter) FOR FIRST LEVEL
    %--------------------------------------------------------------------------
    nbatch = nbatch + 1;
    %Copied from the default batch script
    % Grey matter, white matter
    % but leave csf out!
    matlabbatch{nbatch}.spm.util.imcalc.input          = spm_file(brainMasks{i},'prefix','rf');...
    outfile=spm_file(brainMasks{i},'prefix','rf');
    matlabbatch{nbatch}.spm.util.imcalc.output         = outfile{1};
    matlabbatch{nbatch}.spm.util.imcalc.outdir         = cellstr(fullfile(config.data.bold,config.strID,'preproc',config.sessions{i},'anat'));
    matlabbatch{nbatch}.spm.util.imcalc.expression     = 'i1 > 0.2'; % 0.2 from a reference from jonathan peele: http://jpeelle.net/mri/misc/creating_explicit_mask.html
    matlabbatch{nbatch}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
    matlabbatch{nbatch}.spm.util.imcalc.options.dmtx   = 0;
    matlabbatch{nbatch}.spm.util.imcalc.options.mask   = 0;
    matlabbatch{nbatch}.spm.util.imcalc.options.interp = 1;
    matlabbatch{nbatch}.spm.util.imcalc.options.dtype  = 2;

end   


%% Normalise:brain masks
%-------------------------------------------------------------------------
% normalise anatomical files
for i = 1:length(config.sessions)
    nbatch=nbatch+1;
    matlabbatch{nbatch}.spm.spatial.normalise.write.subj.def = cellstr(spm_select('FPList',config.data.deriv.spmMB.ses(i).anat ,'^y_sub'));
    matlabbatch{nbatch}.spm.spatial.normalise.write.subj.resample = [skullStrippedt1{i} cellstr(spm_file(c1{i},'prefix','rf')) cellstr(spm_file(brainMasks{i},'prefix','rf'))]';
    matlabbatch{nbatch}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72 
                                                                    90 90 108];
    matlabbatch{nbatch}.spm.spatial.normalise.write.woptions.vox = [0.6 0.6 0.6];
    matlabbatch{nbatch}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{nbatch}.spm.spatial.normalise.write.woptions.prefix = 'w';

    %% binarize brain masks again 
    % because coregistration and normalization includes some tranformation and interpolation which creates bad masks for fl explicit masking
    % binarize brain mask (grey matter)
    %--------------------------------------------------------------------------
    nbatch = nbatch + 1;
    %Copied from the default batch script
    % Grey matter, white matter
    % but leave csf out!
    matlabbatch{nbatch}.spm.util.imcalc.input          = spm_file(c1{i},'prefix','wrf');...
    outfile=spm_file(c1{i},'prefix','wrf');
    matlabbatch{nbatch}.spm.util.imcalc.output         = outfile{1};
    matlabbatch{nbatch}.spm.util.imcalc.outdir         = cellstr(fullfile(config.data.bold,config.strID,'preproc',config.sessions{i},'anat'));
    matlabbatch{nbatch}.spm.util.imcalc.expression     = 'i1 > 0.2'; 
    matlabbatch{nbatch}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
    matlabbatch{nbatch}.spm.util.imcalc.options.dmtx   = 0;
    matlabbatch{nbatch}.spm.util.imcalc.options.mask   = 0;
    matlabbatch{nbatch}.spm.util.imcalc.options.interp = 1;
    matlabbatch{nbatch}.spm.util.imcalc.options.dtype  = 2;
    
    % Binarize brain mask (grey+white matter)
    %--------------------------------------------------------------------------
    nbatch = nbatch + 1;
    %Copied from the default batch script
    % Grey matter, white matter
    % but leave csf out!
    matlabbatch{nbatch}.spm.util.imcalc.input          = spm_file(brainMasks{i},'prefix','wrf');...
    outfile=spm_file(brainMasks{i},'prefix','wrf');
    matlabbatch{nbatch}.spm.util.imcalc.output         = outfile{1};
    matlabbatch{nbatch}.spm.util.imcalc.outdir         = cellstr(fullfile(config.data.bold,config.strID,'preproc',config.sessions{i},'anat'));
    matlabbatch{nbatch}.spm.util.imcalc.expression     = 'i1 > 0.2'; % 0.2 from a reference from jonathan peele: http://jpeelle.net/mri/misc/creating_explicit_mask.html
    matlabbatch{nbatch}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
    matlabbatch{nbatch}.spm.util.imcalc.options.dmtx   = 0;
    matlabbatch{nbatch}.spm.util.imcalc.options.mask   = 0;
    matlabbatch{nbatch}.spm.util.imcalc.options.interp = 1;
    matlabbatch{nbatch}.spm.util.imcalc.options.dtype  = 2;
end

%% Create cross-session grey+white matter mask.
%--------------------------------------------------------------------------
nbatch = nbatch + 1;
% Copied from the default batch script
% Grey matter, white matter
% but leave csf out!
matlabbatch{nbatch}.spm.util.imcalc.input          = [spm_file(brainMasks{1},'prefix','wrf'),spm_file(brainMasks{2},'prefix','wrf')]';
matlabbatch{nbatch}.spm.util.imcalc.output         = 'crossSessionWhiteandGreyMatterMaskMNI';
matlabbatch{nbatch}.spm.util.imcalc.outdir         = cellstr(fullfile(config.data.bold,config.strID,'preproc'));
matlabbatch{nbatch}.spm.util.imcalc.expression     = '(i1 + i2) >= 1'; % 0.2 from a reference from jonathan peele: http://jpeelle.net/mri/misc/creating_explicit_mask.html
matlabbatch{nbatch}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
matlabbatch{nbatch}.spm.util.imcalc.options.dmtx   = 0;
matlabbatch{nbatch}.spm.util.imcalc.options.mask   = 0;
matlabbatch{nbatch}.spm.util.imcalc.options.interp = 1;
matlabbatch{nbatch}.spm.util.imcalc.options.dtype  = 2;



%% normalise: functional files
for i=1:length(config.session_signs)
    ind=config.session_signs{i};
    nbatch=nbatch+1;
    matlabbatch{nbatch}.spm.spatial.normalise.write.subj.def = cellstr(spm_select('FPList',config.data.deriv.spmMB.ses(ind).anat ,'^y_sub'));
    matlabbatch{nbatch}.spm.spatial.normalise.write.subj.resample = config.data.ra_files{i};
    matlabbatch{nbatch}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72 
                                                                    90 90 108];
    matlabbatch{nbatch}.spm.spatial.normalise.write.woptions.vox = [config.params.res,  config.params.res,  config.params.res];
    matlabbatch{nbatch}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{nbatch}.spm.spatial.normalise.write.woptions.prefix = 'w';
    
    %% SMOOTH: 
    % functional native space
    nbatch = nbatch + 1;
    matlabbatch{nbatch}.spm.spatial.smooth.data   = config.data.ra_files{i};
    matlabbatch{nbatch}.spm.spatial.smooth.fwhm   = [fwhm,fwhm,fwhm];
    matlabbatch{nbatch}.spm.spatial.smooth.dtype  = 0;
    matlabbatch{nbatch}.spm.spatial.smooth.im     = 0;
    matlabbatch{nbatch}.spm.spatial.smooth.prefix = ['s',num2str(fwhm)];

    % SMOOTH: functional normalized
    nbatch = nbatch + 1;
    matlabbatch{nbatch}.spm.spatial.smooth.data   = cellstr(spm_file(config.data.ra_files{i},'prefix','w')) ;
    matlabbatch{nbatch}.spm.spatial.smooth.fwhm   = [6,6,6];
    matlabbatch{nbatch}.spm.spatial.smooth.dtype  = 0;
    matlabbatch{nbatch}.spm.spatial.smooth.im     = 0;
    matlabbatch{nbatch}.spm.spatial.smooth.prefix = 's6';
    
    % SMOOTH: functional normalized
    nbatch = nbatch + 1;
    matlabbatch{nbatch}.spm.spatial.smooth.data   = cellstr(spm_file(config.data.ra_files{i},'prefix','w')) ;
    matlabbatch{nbatch}.spm.spatial.smooth.fwhm   = [3,3,3];
    matlabbatch{nbatch}.spm.spatial.smooth.dtype  = 0;
    matlabbatch{nbatch}.spm.spatial.smooth.im     = 0;
    matlabbatch{nbatch}.spm.spatial.smooth.prefix = 's3';

end
%% Save batch and run it
save(fullfile(config.data.bold,'2_restofpreprocessing'),'matlabbatch')
spm_jobman('run',matlabbatch);

disp(['... pre-processing part 2 done on subj ', num2str(subj)])
return
