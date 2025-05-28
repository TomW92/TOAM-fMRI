function preprocessing_1s(subj)
% preprocessing_1s -  Does second step of preprocessing:
% DOES: segment, brain extraction, create grey+white matter mask
% Syntax: [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
% subj - string of subject name
%
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

if ~isnumeric(subj)
    subj = str2double(subj);
end

%Setup the scripts according to setup_config in the same script folder
config = setup_config(subj);

for sess_ind = 1:length(config.sessions)

    %make sure no other batch is loaded
    clear matlabbatch
    
    %get t1
    t1 = spm_select('ExtFPList', config.data.deriv.spmMB.ses(sess_ind).anat, '^sub.*T1w.nii$');
    assert(~isempty(t1),'NO structural image found')

    %counter
    nbatch=0;

    
    %% Segment and implicitly noramlize (via TPM, creates forward and inverse Deformation fields)
    %--------------------------------------------------------------------------
    nbatch = nbatch+1;
    matlabbatch{nbatch}.spm.spatial.preproc.channel.vols  = {t1};
    matlabbatch{nbatch}.spm.spatial.preproc.channel.write = [1 1];
    
    %Save only c1, c2, c3 and rc1, rc2
    matlabbatch{nbatch}.spm.spatial.preproc.tissue        = struct('native',{[1 0],[1 0],[1 0],[0 0],[0 0],[0 0]},...
        'ngaus',{1,1,2,3,4,2},...
        'warped',[0 0],...
        'tpm',cellfun(@(x) {horzcat(fullfile(spm('dir'),'tpm','TPM.nii,'),x)},cellstr(num2str([1:6]')),'UniformOutput',false)');
    matlabbatch{nbatch}.spm.spatial.preproc.warp.write = [1 1];
    
    
    
    %% Create brain image
    %--------------------------------------------------------------------------
    nbatch = nbatch + 1;
    %Copied from the default batch script
    % Grey matter, white matter and csf
    matlabbatch{nbatch}.spm.util.imcalc.input          = {spm_file(t1,'prefix','c1');...
        spm_file(t1,'prefix','c2');...
        spm_file(t1,'prefix','c3');...
        t1};
    outfile=spm_file(t1,'prefix','skullStripped');
    matlabbatch{nbatch}.spm.util.imcalc.output         = outfile(1:end-2);
    matlabbatch{nbatch}.spm.util.imcalc.outdir         = cellstr(fullfile(config.data.bold,config.strID,'preproc',config.sessions{sess_ind},'anat'));
    matlabbatch{nbatch}.spm.util.imcalc.expression     = '(i1 + i2 + i3) .* i4';
    matlabbatch{nbatch}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
    matlabbatch{nbatch}.spm.util.imcalc.options.dmtx   = 0;
    matlabbatch{nbatch}.spm.util.imcalc.options.mask   = 0;
    matlabbatch{nbatch}.spm.util.imcalc.options.interp = 1;
    matlabbatch{nbatch}.spm.util.imcalc.options.dtype  = 4;
    
    
    %% Create brain mask (grey+white matter)
    %--------------------------------------------------------------------------
    nbatch = nbatch + 1;
    %Copied from the default batch script
    % Grey matter, white matter
    % but leave csf out!
    matlabbatch{nbatch}.spm.util.imcalc.input          = {spm_file(t1,'prefix','c1');...
        spm_file(t1,'prefix','c2')};
    outfile=spm_file(t1,'prefix','brainMask_');
    matlabbatch{nbatch}.spm.util.imcalc.output         = outfile(1:end-2);
    matlabbatch{nbatch}.spm.util.imcalc.outdir         = cellstr(fullfile(config.data.bold,config.strID,'preproc',config.sessions{sess_ind},'anat'));
    matlabbatch{nbatch}.spm.util.imcalc.expression     = '(i1 + i2) > 0.2'; % 0.2 from a reference from jonathan peele: http://jpeelle.net/mri/misc/creating_explicit_mask.html
    matlabbatch{nbatch}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
    matlabbatch{nbatch}.spm.util.imcalc.options.dmtx   = 0;
    matlabbatch{nbatch}.spm.util.imcalc.options.mask   = 0;
    matlabbatch{nbatch}.spm.util.imcalc.options.interp = 1;
    matlabbatch{nbatch}.spm.util.imcalc.options.dtype  = 2;

    %% Save the script used
    save(fullfile(config.data.bold,config.strID,'preproc',config.sessions{sess_ind},'/anat/','structbatch.mat'),'matlabbatch')

    %Run
    spm_jobman('run',matlabbatch);
    clear matlabbatch
    disp(['... structural pre-processing done for session ', num2str(sess_ind)])
end

disp(['... structural pre-processing done on subj ', num2str(subj)])

