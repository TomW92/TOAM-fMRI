function [conds, contrast, convecs, onsets, durations] = contrasts_ret_FIR(behavior)
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
% Author: Mirko Bristle

contrast = 0;
convecs  = 0;
conds = {'faceRecog';'ret'};
onsets=cell(1,length(conds));
durations=cell(1,length(conds));

%loop over trial types
    for ttypes_ind = 1:length(conds)
        ttype = char(conds(ttypes_ind));
        ttypeRows = strcmp(behavior.trial_type,ttype);
        ttypeEvents = behavior(ttypeRows,:);

        %Create the cell arrays
        onsets{ttypes_ind} = ttypeEvents.onset;
        durations{ttypes_ind} = ttypeEvents.duration;
    end

    
%% check output
    assert(all( cellfun(@(x) ischar(x),conds)))
    %assert(all( cellfun(@(x) ischar(x),contrast)))
    %assert(all( cellfun(@(x) isnumeric(x),convecs))) 
end