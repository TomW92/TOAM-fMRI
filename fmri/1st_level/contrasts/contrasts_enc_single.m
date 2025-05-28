function [conds, contrast, convecs, onsets, durations] = contrasts_enc_single(behavior,condition,model_duration)
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
    condition='Enc';
end


trial_types = unique(behavior.trial_type); %list of trial types: Face - Enc - Recall - Feedback
    
%loop over trial types
    for ttypes_ind = 1:length(trial_types)
        ttype = char(trial_types(ttypes_ind));
        ttypeRows = strcmp(behavior.trial_type,ttype);
        ttypeEvents = behavior(ttypeRows,:); % selects all rows with the trial type Enc --> ttypeEvents


        if strcmp(ttype,condition) %Returns true if ttype and condition are the same string --> loop is only executed for the condition Enc
    
            %Create the cell arrays
            onsets = num2cell(ttypeEvents.onset)';
            conds = ttypeEvents.item';         %{[char(acc_Events.trial_type(1)), '_', conscLabel '_' char(current_accLevel)]}; %ja, char!
            if model_duration
                durations = repelem(num2cell(2.75),length(conds)); %repeats elements of cell array for length of conds
            else
                durations = repelem(num2cell(0),length(conds));
            end
              
            %actually save them. no we do no longer do this.
            %filename = [events_files(1).folder, filesep, events_files(1).name(1:9),'_EncEvent_5s_SINGLE_TRIAL_SECONDS_ONSET.mat'];
            %save(filename,'names','onsets','durations')%,'pmod')
        end
        %pmod(ttypes_ind).name=strcat(ttype,"_SubsequentConsc");
        %pmod(ttypes_ind).param=num2cell(ttypeEvents.accuracy,1);
        %pmod(ttypes_ind).poly={1};
    end
    
    k=1;
    contrast = cell(length(conds),1); %cell array with one column
    convecs  = cell(length(conds),1);
    for i = 1:length(conds) % Creates a contrast vector for each condition = item
       contrast{k}=conds{i};
       convecs{k}=double(cellfun(@(x,y) strcmp(x,conds{i}),conds)); %creates a vector with 1 for the condition and 0 for all other conditions
       k=k+1;
    end

%% check output
    assert(all( cellfun(@(x) ischar(x),conds)))
    assert(all( cellfun(@(x) ischar(x),contrast)))
    assert(all( cellfun(@(x) isnumeric(x),convecs))) 
end