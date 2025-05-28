function [conds, contrast, convecs, onsets, durations] = contrasts_delret_subsequent(behavior,model_duration)
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

if nargin<2
    model_duration=0;
end

custom_contrast = 0;
custom_convecs  = 0;
trial_types = unique(behavior.trial_type);
trial_types = trial_types(~strcmp(trial_types,'retFeedback'));
trial_types = trial_types(~strcmp(trial_types,'retRatingFeedback'));
trial_types = trial_types(~strcmp(trial_types,'recogRatingFeedback'));


%loop over trial types
onsets=cell(1,2);
durations=cell(1,2);
k=0;
for ttypes_ind = 1:length(trial_types)
    ttype = char(trial_types(ttypes_ind));
    ttypeRows = strcmp(behavior.trial_type,ttype);
    ttypeEvents = behavior(ttypeRows,:);

    ttype_cat = unique(ttypeEvents.accuracy);
    for cat_ind = 1:length(ttype_cat)     
        k=k+1;
        cat = char(ttype_cat(cat_ind));
        catRows= strcmp(ttypeEvents.accuracy,cat);
        catEvents = ttypeEvents(catRows,:);

        onsets{k} = catEvents.onset;
        conds{k}  = strcat(ttype,'_',cat);         
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
contrast = cell(2,1);
convecs  = cell(2,1);
for i = 1:2%length(conds)
   contrast{k}=conds{i};
   convecs{k}=double(cellfun(@(x,y) strcmp(x,conds{i}),conds));
   k=k+1;
end

custom_contrast = {'later_consc > stays_unconsc','stays_unconsc > later_consc',...
                    'faceRecog_consc > later_consc',...
                    'faceRecog_consc > stays_unconsc'};
custom_convecs = cell(length(custom_contrast),1);
for con_num = 1:length(custom_contrast)
    con_names = strsplit(custom_contrast{con_num},{'>',' '});
    level = con_names{1};
    negLevel = con_names{2};
    custom_convecs{con_num}=cellfun(@(x,y) ~isempty(regexp(x,['^.*',level,'(_|\>)'],'match')),conds)...  % '(_|\>)' in regexp means to look for underscore (_) , but because it should work at the end of the string to, make logical OR (|) and look for the end of the word (\>)
                    +(-1)*cellfun(@(x,y) ~isempty(regexp(x,['^.*',negLevel,'(_|\>)'],'match')),conds);
            
end

contrast = vertcat(contrast,custom_contrast');
convecs  = vertcat(convecs,custom_convecs);

%% check output
assert(all( cellfun(@(x) ischar(x),conds)))
assert(all( cellfun(@(x) ischar(x),contrast)))
assert(all( cellfun(@(x) isnumeric(x),custom_convecs))) 
end