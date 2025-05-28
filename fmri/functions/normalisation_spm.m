function normalisation_spm(subj)
%%%
%
% Normalisation of the freesurfer masks with the deformation field for that
% very session. Freesurfer - recon all - was run with the biascorrected image and the
% deformation field from this segmentation (SPM) is now used to create the
% mni registered masks.
%
%
%
% IF no argument is supplied, the subject is asked
if nargin<1
    subj = 60601;
    
end
%% 
p=pwd;
idcs=strfind(p,'/');
%addpath(fullfile(p(1:idcs(end)-1),'functions'))
disp(pwd)
disp("pwd")
if ~isnumeric(subj)
    subj = str2double(subj);
end

addpath('/storage/homefs/tw18a205/toolboxes/spm12');
spm('Defaults','fMRI');
spm_jobman('initcfg');


config = setup_config(subj);

nbatch = 0;
for i = 1:length(config.sessions)

    %load mask_images
    mask_images = config.data.deriv.freesurfer.ses(i).masks;

    %because spm only accepts unzipped niftis, unzip the gzipped masks
    %only unzip if file doesn't exist already.
    for mask_ind = 1:length(mask_images)
        if endsWith(mask_images{mask_ind},'.gz')
            %[parent,file,~] = fileparts(mask_images{mask_ind});
            %if ~isfile(fullfile(parent,file))
            gunzip(mask_images{mask_ind});
            %end
        end
    end

    %reload just the unzipped niftimasks
    unzipped_mask_images = cellstr(spm_select('FPListRec',config.data.deriv.freesurfer.ses(i).mri,'^(?!w).*^(?!.*mni).*(?!.*ses-2).*(?:mask|bilat|lh|rh|hipp).*\.nii$'));

    %complete batch
    nbatch=nbatch+1;
    matlabbatch{nbatch}.spm.spatial.normalise.write.subj.def = cellstr(spm_select('FPList',config.data.deriv.spmMB.ses(i).anat ,'^y_sub'));
    matlabbatch{nbatch}.spm.spatial.normalise.write.subj.resample = unzipped_mask_images;
    matlabbatch{nbatch}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72 
                                                                    90 90 108];
    matlabbatch{nbatch}.spm.spatial.normalise.write.woptions.vox = [config.params.res config.params.res config.params.res];
    matlabbatch{nbatch}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{nbatch}.spm.spatial.normalise.write.woptions.prefix = 'w';

    for mask_ind = 1:length(unzipped_mask_images)
        nbatch = nbatch + 1;
        %Copied from the default batch script
        % Grey matter, white matter
        % but leave csf out!
        matlabbatch{nbatch}.spm.util.imcalc.input          = cellstr(spm_file(unzipped_mask_images{mask_ind},'prefix','w'));...
        outfile=spm_file(unzipped_mask_images{mask_ind},'prefix','w');
        matlabbatch{nbatch}.spm.util.imcalc.output         = outfile;
        matlabbatch{nbatch}.spm.util.imcalc.outdir         = cellstr(fullfile(config.data.deriv.freesurfer.ses(i).mri));
        matlabbatch{nbatch}.spm.util.imcalc.expression     = 'i1 > 0.2'; % 0.2 from a reference from jonathan peele: http://jpeelle.net/mri/misc/creating_explicit_mask.html
        matlabbatch{nbatch}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
        matlabbatch{nbatch}.spm.util.imcalc.options.dmtx   = 0;
        matlabbatch{nbatch}.spm.util.imcalc.options.mask   = 0;
        matlabbatch{nbatch}.spm.util.imcalc.options.interp = 1;
        matlabbatch{nbatch}.spm.util.imcalc.options.dtype  = 2;
    end
end

save(fullfile(config.data.workspace,'/freesurfer/','normalise_freesurferMasks'),'matlabbatch')
spm_jobman('run',matlabbatch);

disp(['normalized subject-specific mask images of the hippocampus and the hippocampal subfields: ', num2str(subj)])

end
