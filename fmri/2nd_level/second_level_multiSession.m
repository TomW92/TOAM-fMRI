function second_level_multiSession(subj, sess, behavior_name, myContrast)
% second_level_acc_cat_cond - calculates the second level analysis 
%  - GLM and behavioral correlation
%
% Other m-files required: run2ndLvlTcon.m, setup_config.m, 
% Subfunctions: none
% MAT-files required: SPM.mat of first level
%
% See also: first_level.m, 
% Author: Mirko Bristle, Tom Willems

if nargin <1
   subj=60100;   %input('Subject: ');
   sess=1;
   behavior_name='r';
   myClassic = '1';
   myContrast= 'wpt';
   myMix='';
end

myMix='';
myClassic=1;

if ~isnumeric(sess)
    sess = str2double(sess);
end
if ~isnumeric(myClassic)
	myClassic = str2double(myClassic);
end

%% Setup

if myClassic
    myPrefix = [myMix,'s6wra'];
else
    myPrefix = [myMix,'sr'];
end

analysis_name=['TContrasts_fl_explmask'];
L1_p='l1_expl_mask_s6wr';

nbatch=0;
clear matlabbatch;

p=pwd;
idcs=strfind(p,'/');
addpath(fullfile(p(1:idcs(end)-1),'functions'))

%Setup the scripts according to setup_config in the same script folder
config = setup_config(subj);
session=config.sessions{sess};
disp(session)


subjList=config.seq.participants.participant_id; %unique(behavior.id);
subjIdx=arrayfun(@(x) isfile(fullfile(config.data.bold,x,L1_p,'con_0001.nii')),subjList)>0;
subjList=subjList(subjIdx);

subjPathList = strcat(config.data.bold,filesep,subjList);
nsubj=size(subjList,1);

spm_p=fullfile(subjPathList(1,:),L1_p,'SPM.mat');
load(spm_p{1})

assert(nsubj>0,['No Subject available for second-level analysis. Check session: ',session])

% Creation of directory for the rfx
rfxDir = fullfile(config.data.bold,['RFX_',analysis_name]);
mkdir(rfxDir)
cd(rfxDir)


%% Analysing each t contrast
% load the names of the contrasts from the first subject
% it use the contrasts_name matrix created in the contrasts.m script

arr={SPM.xCon.name};

repeatedContrastCounter = 1;
for icon = 1:length(arr) %f_cons+1:(nb_con+f_cons)
    nbatch=nbatch+1;
    % loop over contrasts

    disp(['Processing contrast : ', SPM.xCon(icon).name]);

    %Create the folder of that contrast
    condDir = strrep(SPM.xCon(icon).name,' > ','_');
    
    if isfolder(condDir)
        if strcmp(SPM.xCon(icon).name,'faceRecog_correct_consc')
            repeatedContrastCounter=repeatedContrastCounter+1;
        end
        condDir = [condDir, '_', num2str(repeatedContrastCounter)];
        mkdir(condDir);
        
    else 
        mkdir(condDir);
    end
 

    %Put the files in there
    FPcondDir=cellstr(fullfile(rfxDir,condDir));
    matlabbatch{nbatch}.spm.stats.factorial_design.dir = FPcondDir;

    %Expect all the contrast structure to be identical
    files = strcat(subjPathList,filesep,fullfile(L1_p,sprintf('con_%04d.nii',icon)),',1');
    matlabbatch{nbatch}.spm.stats.factorial_design.des.t1.scans = cellstr(files);

    %Rest was left by default
    matlabbatch{nbatch}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{nbatch}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{nbatch}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{nbatch}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{nbatch}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{nbatch}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{nbatch}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{nbatch}.spm.stats.factorial_design.globalm.glonorm = 1;

    %Estimate it directly
    nbatch=nbatch+1;
    matlabbatch{nbatch}.spm.stats.fmri_est.spmmat(1) =cellstr(fullfile(FPcondDir,'SPM.mat'));
    matlabbatch{nbatch}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{nbatch}.spm.stats.fmri_est.method.Classical = 1;

    %Create a t-test contrast
    nbatch=nbatch+1;
    matlabbatch{nbatch}.spm.stats.con.spmmat(1) = cellstr(fullfile(FPcondDir,'SPM.mat'));
    matlabbatch{nbatch}.spm.stats.con.consess{1}.tcon.name = SPM.xCon(icon).name;
    matlabbatch{nbatch}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{nbatch}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{nbatch}.spm.stats.con.delete = 0;

    nbatch=nbatch+1;
    matlabbatch{nbatch}.spm.stats.results.spmmat(1) = cellstr(fullfile(FPcondDir,'SPM.mat'));
    matlabbatch{nbatch}.spm.stats.results.conspec.titlestr = SPM.xCon(icon).name;
    matlabbatch{nbatch}.spm.stats.results.conspec.contrasts = 1;
    matlabbatch{nbatch}.spm.stats.results.conspec.threshdesc = 'none';
    matlabbatch{nbatch}.spm.stats.results.conspec.thresh = 0.001;
    matlabbatch{nbatch}.spm.stats.results.conspec.extent = 0;
    matlabbatch{nbatch}.spm.stats.results.conspec.conjunction = 1;
    matlabbatch{nbatch}.spm.stats.results.conspec.mask.none = 1;
    matlabbatch{nbatch}.spm.stats.results.units = 1;
    matlabbatch{nbatch}.spm.stats.results.export{1}.ps = true;
    matlabbatch{nbatch}.spm.stats.results.export{1}.pdf = true;


end %end of condition loop

%% Save the used batch in the folder
save(fullfile(rfxDir,filesep,'batch.mat'), 'matlabbatch')

%Do the actual calculation
spm_jobman('run',matlabbatch)
clear matlabbatch

end





