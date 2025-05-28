function run1stLevelLSS(subject,session)
%
% this is looping over all trials per session per subject
%
if nargin<1
    subject = 60601;
    session = 'Enc';
end
%
%
% add functions
p=pwd;
idcs=strfind(p,filesep);
addpath(fullfile(p(1:idcs(end)-1),'functions'))
addpath('contrasts')

if strcmp(session,'Enc')
    sess = 1;
    sessionName = 'forgetexp';
    myContrast = [{'enc_LSS'},{''}];
    myBehavior = 'enc_bids';
elseif strcmp(session,'Ret1')
    sess = 1;
    sessionName = 'retrieval';
    myContrast = [{'ret1_LSS'},{''}];
    myBehavior = 'ret_bids';
elseif strcmp(session,'Ret2')
    sess = 1;
    sessionName = 'retrieval2';
    myContrast = [{'ret2_LSS'},{''}];
    myBehavior = 'ret2_bids';
end

mniSpace = 0;
smooth = 1;
test_cases = 0;

config=setup_config(subject);
behavior = readtable([config.data.behav.maxTrials.(myBehavior)],'FileType','delimitedtext'); 
trials = unique(behavior.item)';

for trial = trials 
    % disp(trial)
    myContrast{2} = trial;
 first_level(subject, sess, sessionName, mniSpace ,myContrast, myBehavior, smooth, test_cases)
end