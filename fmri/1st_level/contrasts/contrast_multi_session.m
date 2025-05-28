function multiSessionContrasts = contrast_multi_session(config)
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

if nargin<1
    config=setup_config(60100);
end

    %% setup Contrasts -------------------------------------------------------
behaviorNames = {'forgetexp','consolidation','retrieval','forgetexp'};
possibleContrasts   = {{'enc_subsequent','consolidation','ret_twoEvents','delret'},...
                    {'another_enc_contrast','enc_bids_forinstance?','add_later','just to make it more flexible'}};
onsetFilesNames = {'enc_subsequent','consolidation','ret','delret'};

%this takes onset files with max amount of trials, minimal scrubbing of
%incomplete cases etc. 
onset_type = 'maxTrials';
myContrasts = possibleContrasts{1};
myPrefix = 's6wra';

multiSessionContrasts = struct;
multiSessionContrasts.fullConvecs = {};
multiSessionContrasts.fullContrasts = {};
multiSessionContrasts.util.fullCondsWithMPs = {};

mvmParams = repelem({'mvm_params'},6);
condition_counter = 0;

for sess_ind = 1:length(myContrasts)
    myContrast = myContrasts{sess_ind};
    behaviorName = behaviorNames{sess_ind};
    onsets = onsetFilesNames{sess_ind};

    %read behavior (output from R script)
    disp("creating onsets with this file:")
    disp(config.data.behav.(onset_type).(onsets))
    behavior = readtable([config.data.behav.(onset_type).(onsets)]);
    assert(any(cell2mat(arrayfun(@(x) x==config.id,unique(behavior.id),'UniformOutput',false))),'subject was discarded in behavior analysis')
    specify_contrasts = str2func(['contrasts_',myContrast]);
    
    [conds, contrasts, convecs, onsets, durations] = specify_contrasts(behavior);
    
    %the following loop takes the conditions and searches in the earlier
    %done conditions, adds a 2 to the condition in case of repetition
    %only start doing so after the first session is collected, for obvious
    %reasons.
    if sess_ind>1
        for cond_ind = 1:length(conds)
            for second_sess_ind = 1:length(multiSessionContrasts.ses)
                conds_struct2search = multiSessionContrasts.ses(second_sess_ind).conds;
                condition = conds{cond_ind};
                if any(strcmp(conds_struct2search, condition))
                    conds{cond_ind} = [condition,' ', num2str(config.session_signs{sess_ind})];
                end
            end
        end
    end

    multiSessionContrasts.ses(sess_ind).durations = durations;
    multiSessionContrasts.ses(sess_ind).conds = conds;
    multiSessionContrasts.ses(sess_ind).onsets = onsets;

    %for later use in contrast
    convecs = cellfun(@(x) [zeros(1,condition_counter),x],convecs,'UniformOutput',false);
    condition_counter = condition_counter + length(conds);
    condition_counter = condition_counter + 6; % movement parameters
    
    multiSessionContrasts.fullConvecs = vertcat(multiSessionContrasts.fullConvecs,convecs);
    multiSessionContrasts.fullContrasts= vertcat(multiSessionContrasts.fullContrasts,contrasts);

    condsMVM = horzcat(conds, mvmParams);
    if sess_ind == 3
        condsMVM = strcat(condsMVM,' 1');
    elseif sess_ind == 4
        condsMVM = strcat(condsMVM,' 2');
    end
    multiSessionContrasts.util.fullCondsWithMPs = horzcat(multiSessionContrasts.util.fullCondsWithMPs,condsMVM);
    
    dataPre=myPrefix; % Substring if mixedAnalysis
    multiSessionContrasts.ses(sess_ind).files =cellstr(spm_select('ExtFPList', fullfile(config.data.deriv.spmMB.ses(config.session_signs{sess_ind}).func), ['^',dataPre,'.*',behaviorName,'.*\.nii$']));
end

%always put in contrast names in this way
customContrasts = {'faceRecog_consc 1+2 > faceRecog_unconsc 1+2',...
                   'faceRecog_correct_unconsc 1+2 > faceRecog_incorrect_unconsc 1+2',...
                   'faceRecog_correct_unconsc 1 > faceRecog_incorrect_unconsc 1',...
                   'faceRecog_correct_unconsc 2 > faceRecog_incorrect_unconsc 2',...
                   'faceRecog_incorrect_unconsc 1+2 > faceRecog_correct_unconsc 1+2',...
                   'ret_correct_unconsc 1+2 > ret_incorrect_unconsc 1+2',...
                   'faceRecog_consc 1 > Face',...
                   'faceRecog_consc 2 > Face',...
                   'faceRecog_consc 2 > faceRecog_consc 1',...
                   'faceRecog_unconsc 2 > faceRecog_unconsc 1',...
                   'recog_correct_unconsc 2 > recog_incorrect_unconsc 2',...
                   'faceRecog_correct_unconsc 1+2 > number',...
                   'faceRecog_correct_unconsc 2 > number'
                    };

new_custom_convecs = cell(length(customContrasts),1);
for con_num = 1:length(customContrasts)
    con_names = strsplit(customContrasts{con_num},{'>',' ','_'});
    if  length(con_names)==4
        levels = con_names(1:2);
        whichSession = con_names(3);
        negLevel = con_names(4);
        new_custom_convecs{con_num}=...
            cellfun(@(x,y) length([regexp(x,['^(\w*_|)',levels{1},'(_|\>)'],'match'),regexp(x,['^(\w*_|)',levels{2},'(_|\>)'],'match')])==2 && ~isempty(regexp(x,whichSession{1},'match')),multiSessionContrasts.util.fullCondsWithMPs)...  % '(_|\>)' in regexp means to look for underscore (_) , but because it should work at the end of the string to, make logical OR (|) and look for the end of the word (\>)
             +(-1)*cellfun(@(x,y) ~isempty(regexp(x,['^(\w*_|)',negLevel{1},'(_|\>)'],'match')),multiSessionContrasts.util.fullCondsWithMPs);

    elseif length(con_names)==6
        levels = con_names(1:2);
        whichSession1 = con_names{3};
        negLevels = con_names(4:5);
        whichSession2 = con_names{6};

        if (strcmp(whichSession1,'1+2') && strcmp(whichSession2,'1+2'))
            new_custom_convecs{con_num}=...
                cellfun(@(x,y) length([regexp(x,['^(\w*_|)',levels{1},'(_|\>)'],'match'),regexp(x,['^(\w*_|)',levels{2},'(_|\>)'],'match')])==2,multiSessionContrasts.util.fullCondsWithMPs)...  % '(_|\>)' in regexp means to look for underscore (_) , but because it should work at the end of the string to, make logical OR (|) and look for the end of the word (\>)
                 +(-1)*cellfun(@(x,y) length([regexp(x,['^(\w*_|)',negLevels{1},'(_|\>)'],'match'),regexp(x,['^(\w*_|)',negLevels{2},'(_|\>)'],'match')])==2,multiSessionContrasts.util.fullCondsWithMPs);
        else
            new_custom_convecs{con_num}=...
                cellfun(@(x,y) length([regexp(x,['^(\w*_|)',levels{1},'(_|\>)'],'match'),regexp(x,['^(\w*_|)',levels{2},'(_|\>)'],'match')])==2 && ~isempty(regexp(x,whichSession1,'match')),multiSessionContrasts.util.fullCondsWithMPs)...  % '(_|\>)' in regexp means to look for underscore (_) , but because it should work at the end of the string to, make logical OR (|) and look for the end of the word (\>)
                 +(-1)*cellfun(@(x,y) length([regexp(x,['^(\w*_|)',negLevels{1},'(_|\>)'],'match'),regexp(x,['^(\w*_|)',negLevels{2},'(_|\>)'],'match')])==2 && ~isempty(regexp(x,whichSession2,'match')),multiSessionContrasts.util.fullCondsWithMPs);
        end 

    else 
        con_names_simple = strsplit(customContrasts{con_num},{' ','>'});

        if ~strcmp(con_names_simple{2},'1+2')
            con_names_simple{1} = strcat([con_names_simple{1},' ',con_names_simple{2}]);
        end
        if length(con_names_simple)>3 && ~strcmp(con_names_simple{4},'1+2')
            con_names_simple{3} = strcat([con_names_simple{3},' ',con_names_simple{4}]);
        end
        level = strfind(multiSessionContrasts.util.fullCondsWithMPs,con_names_simple{1});        
        negLevel = strfind(multiSessionContrasts.util.fullCondsWithMPs,con_names_simple{3});
        pos_con = ~cellfun(@isempty,level);
        neg_con = ~cellfun(@isempty,negLevel);   

        new_custom_convecs{con_num}=zeros(length(multiSessionContrasts.util.fullCondsWithMPs),1);
        new_custom_convecs{con_num}(pos_con)=1;
        new_custom_convecs{con_num}(neg_con)=-1;
    end
end

for myind = 1:length(new_custom_convecs)
    convecRow = new_custom_convecs{myind};
    if ~(sum(convecRow) == 0)
        %disp(contrast(len+myind))
        numOnes = find(convecRow==1);
        numMinusOnes = find(convecRow==-1);
        convecRow(numOnes) = 1/length(numOnes);
        convecRow(numMinusOnes) = -1/length(numMinusOnes);
        new_custom_convecs(myind) = {convecRow};
    end
end

multiSessionContrasts.fullConvecs   = vertcat(multiSessionContrasts.fullConvecs,new_custom_convecs);
multiSessionContrasts.fullContrasts = vertcat(multiSessionContrasts.fullContrasts,customContrasts');

