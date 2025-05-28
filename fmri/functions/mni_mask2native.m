function mni_mask2native(subj)
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
addpath(fullfile(p(1:idcs(end)-1),'functions'))

if ~isnumeric(subj)
    subj = str2double(subj);
end

%addpath('/storage/homefs/tw18a205/toolboxes/spm12');
addpath('/storage/homefs/fr22c605/matlab/spm12');
spm('Defaults','fMRI');
spm_jobman('initcfg');


config = setup_config(subj);

ses1 = {'l1__native_enc_single_forgetexp_r','l1__native_ret_single_retrieval_r'};
ses2 = {'l1__native_ret_single_forgetexp_r','l1__native_recog_single_forgetexp_r'};
sess = {ses1,ses2};

odir= '/storage/workspaces/psy_memory_wfg_psy/hpc_henke_wfg/s2019_twillems_fMRI_silent_engram/data/fMRI/wb/1st_level/';

for i = 1:length(sess)
    for tp = 1:length(sess{i})
        nbatch = 0;
        unzipped_mask_images = cellstr('/storage/homefs/fr22c605/toolboxes/fsl/data/standard/MNI152_T1_1mm_brain_mask.nii');
        %unzipped_mask_images = cellstr('/storage/homefs/fr22c605/matlab/spm12/canonical/avg152T1.nii');
        %dl=[];
        %dl = [dl; '/storage/homefs/fr22c605/toolboxes/fsl/data/standard/MNI152_T1_1mm_brain_mask.nii'];
        %dl = [dl; '/storage/homefs/fr22c605/toolboxes/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii'];
        %unzipped_mask_images = cellstr(dl)
        
        %display(cellstr(spm_select('FPList',config.data.deriv.spmMB.ses(i).anat ,'^iy_sub')))
        %complete batch
        nbatch=nbatch+1;
        matlabbatch{nbatch}.spm.spatial.normalise.write.subj.def = cellstr(spm_select('FPList',config.data.deriv.spmMB.ses(i).anat ,'^iy_sub'));
        matlabbatch{nbatch}.spm.spatial.normalise.write.subj.resample = {'/storage/homefs/fr22c605/toolboxes/fsl/data/standard/MNI152_T1_1mm_brain_mask.nii,1'};
        matlabbatch{nbatch}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72 
                                                                        90 90 108];
        %matlabbatch{nbatch}.spm.spatial.normalise.write.woptions.vox = [0.8 0.8 0.8];
        matlabbatch{nbatch}.spm.spatial.normalise.write.woptions.vox = [1.104762 1.104762 1.1];
        matlabbatch{nbatch}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{nbatch}.spm.spatial.normalise.write.woptions.prefix = 'native_';

        spm_jobman('run',matlabbatch);
        outdir = fullfile(odir,['sub-' num2str(subj)], sess{i}{tp})
        movefile ('/storage/homefs/fr22c605/matlab/spm12/canonical/native_avg152T1.nii, outdir)
        clear matlabbatch
        display('batch done')
end
end
