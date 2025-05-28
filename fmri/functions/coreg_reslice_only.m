function coreg_reslice_only(ref,moving,prefix)
%%
% input: 
% ref: reference image, to which the moving image is resliced, without
% estimation of coregistration. 
%
% moving: moving image, which is resliced to image grid of the reference
% image
%
%
%
%%
    matlabbatch{1}.spm.spatial.coreg.write.ref = ref;
    matlabbatch{1}.spm.spatial.coreg.write.source = moving;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = prefix;
    
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
end




