function [conds, contrast, convecs, onsets, durations] = contrasts_ret_single(behavior,condition,model_duration)
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

if nargin<3
    model_duration=0;
    condition='faceRecog';
end

contrast = 0;
convecs  = 0;
trial_types = unique(behavior.trial_type);

%loop over trial types
    for ttypes_ind = 1:length(trial_types)
        ttype = char(trial_types(ttypes_ind));
        ttypeRows = strcmp(behavior.trial_type,ttype);
        ttypeEvents = behavior(ttypeRows,:);


        if strcmp(ttype,condition)
    
            %Create the cell arrays
            onsets = num2cell(ttypeEvents.onset)';
            conds = ttypeEvents.item';         
            if model_duration
                durations = repelem(num2cell(1),length(conds));
            else
                durations = repelem(num2cell(0),length(conds));
            end
            
        end
        %pmod(ttypes_ind).name=strcat(ttype,"_SubsequentConsc");
        %pmod(ttypes_ind).param=num2cell(ttypeEvents.accuracy,1);
        %pmod(ttypes_ind).poly={1};
    end
    
    k=1;
    contrast = cell(length(conds),1);
    convecs  = cell(length(conds),1);
    for i = 1:length(conds)
       contrast{k}=conds{i};
       convecs{k}=double(cellfun(@(x,y) strcmp(x,conds{i}),conds));
       k=k+1;
    end

%% check output
    assert(all( cellfun(@(x) ischar(x),conds)))
    assert(all( cellfun(@(x) ischar(x),contrast)))
    assert(all( cellfun(@(x) isnumeric(x),convecs))) 
end