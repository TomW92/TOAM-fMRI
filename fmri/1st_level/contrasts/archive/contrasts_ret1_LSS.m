function [conds, contrast, convecs, onsets, durations] = contrasts_ret1_LSS(behavior,condition,model_duration)
% contrast_acc_cond - Calculates the T and F contrast of the first level model
%
% Syntax: [conds, contrast, convecs, onsets] = contrast_acc_cond(config,behavior)
%
% Inputs:
% config - config struct. holds general information of the analysis
%               (session, id, sessionOnsetCorrection)
% behavior - input table of the behavior with onset times, conditions, 
%               subjects, session, etc. 
%
% Outputs:
% conds - names of conditions of analysis
% contrast - contrast names of analysis
% convecs - convecs for spm for corresponding to contrasts
% onsets - onsets of contrasts
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: setup_config.m first_level.m
% Author: Mirko Bristle, Tom Willems

if nargin<3
    model_duration=0;
    condition='Ret1';
end

trialOfInterest = behavior{2};
behavior = behavior{1};

index_ret1_trials = strcmp(behavior.trial_type, 'faceRecog');
ret1_trials = behavior(index_ret1_trials, :);

% CREATE THE CONDITION FOR TRIAL OF INTEREST

conds = {trialOfInterest{1}, 'rest'};

durations = repelem(num2cell(0),length(conds));

% GET THE ONSETS OF TRIAL OF INTEREST AND ALL OTHERS

% Create a logical index that is true for each cell in the 'item' column that equals the trialOfInterest and use it to get the corresponding 'onset' values
logicalIndex2 = strcmp(trialOfInterest, ret1_trials.item);
onsetTrialOfInterest = ret1_trials.onset(logicalIndex2);

logicalIndex3 = ~strcmp(trialOfInterest, ret1_trials.item);
onsetsOthers = ret1_trials.onset(logicalIndex3);

onsets = {onsetTrialOfInterest, onsetsOthers}; % should be cell array with one row and two columns: onset of trial of interest vs. onsets of all other trials

% CREATE A CONTRAST AND THE CONTRAST vector
contrast = {trialOfInterest{1},'restOfTrials'};
convecs = {[1, 0],[0,1]};

% CHECK OUTPUT
disp(['contrast: ', contrast])
disp(['condition', conds])
disp(['contrast vector', convecs])
disp(['onsets of interest', onsets])




% LSS: just one contrast: {first} item vs. all other items
% contrast = conds{1}
% convec: contains just the first contrast vector: first item - all other items 
% convec = cell(length(conds),1);
% convec_first = double(convecs{1, 0})
% convec {1, -1} = 1x1 cell array with the first column containing contrast for first item (1) and the second column containing contrast for all other items (-1)

% old version: filteredItems = ttypeEvents.item(logicalIndex1, column);
% old version: conds = {ttypeEvents.item(i), ttypeEvents.item(filtereditems)}