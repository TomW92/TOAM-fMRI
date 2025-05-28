function [conds, contrast, convecs, onsets, durations] = contrasts_enc_subsequent(behavior)
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

    %% setup Contrasts -------------------------------------------------------
    NUM_OF_CONTRAST=4*3+4*2+3*2+2*4+2*3;
     
    %hard coded conditions!
    conditions.levels{1}={'Face','Enc','Recall','Feedback'};
    conditions.levels{2}={'correct','incorrect'};
    conditions.levels{3}={'subseqconsc','subsequnconsc','subseqbetween' };
   
    k=1;
    for i=1:length(conditions.levels{1})
        for ii=1:length(conditions.levels{2})
            for iii=1:length(conditions.levels{3})
                conditions.names(k)=strcat(conditions.levels{1}(i),'_',...
                    conditions.levels{2}(ii),'_',...
                    conditions.levels{3}(iii));
                k=k+1;
            end
        end
    end
    %conditions.names(k)=cellstr(conditions.baseline);
    
    %delete consc_incorrect (not having it!)
    S = ~contains(conditions.names,'incorrect_subseqconsc');
    conditions.names = conditions.names(S);

    %concat levels of all factors
    conds = conditions.names;
    len=length(conds);

    k=1;

    contrast=cell(NUM_OF_CONTRAST,1);
    convecs=cell(NUM_OF_CONTRAST,1);

    %% F contrasts
    for i = 1:length(conditions.names)
       contrast{k}=conditions.names{i};
       convecs{k}=double(cellfun(@(x,y) strcmp(x,conditions.names{i}),conditions.names));
       k=k+1;
    end
    %% T contrasts
    for i = 1:length(conditions.levels)
        for l = 1:length(conditions.levels{i})
            levels=[conditions.levels{i}]; % ,conditions.baseline]
            level=levels{l};
            for negl = find(~strcmp(levels,level))
                negLevel=levels{negl};
                contrast{k}=[level,' > ',negLevel];
                convecs{k}=...
                    cellfun(@(x,y) ~isempty(regexp(x,['^(\w*_|)',level,'.*'],'match')),conds)...
                    +(-1)*cellfun(@(x,y) ~isempty(regexp(x,['^(\w*_|)',negLevel,'.*'],'match')),conds);
                k=k+1;
            end
        end
    end


    for myind = 1:length(convecs(len+1:end))
        convecRow = convecs{len+myind};
        if ~(sum(convecRow) == 0)
            %disp(contrast(len+myind))
            numOnes = find(convecRow==1);
            numMinusOnes = find(convecRow==-1);
            convecRow(numOnes) = 1/length(numOnes);
            convecRow(numMinusOnes) = -1/length(numMinusOnes);
            convecs(len+myind) = {convecRow};
        end
    end


    %convecs=vertcat(convecs(1:len), cellfun(@(x) [ x(1:end-1),abs(sum(x(1:end-1))).*x(end)],convecs((len+1):end),'UniformOutput',false)); %whats that though?
    assert(all(cellfun(@(x) sum(x),convecs((len+1):end))==0))

    %% Calc onsets
    onsets=cell(1,length(conds));
    durations=cell(1,length(conds));
    for i = 1:length(conds)
        %disp(conds{i})
        %calc onsets
        %c_session = strcmp(behavior.session,config.c_session);
        %c_levels=regexp(conds{i},'(?<cond>[a-z]*)_(?<acc>[a-z]*)','names','ignorecase');
        c_levels=regexp(conds{i},'(?<cond>[a-z]*)_(?<acc>[a-z]*)_(?<conscLevel>\w*)','names','ignorecase');

%         if ~isempty(c_levels.baseline)
%             c_levels.cond=c_levels.baseline;
%             c_levels.acc=c_levels.baseline;
%         end

        if isnumeric(behavior.rating)
            if strcmp(c_levels.conscLevel,'subsequnconsc')
                c_levels.conscLevel = 2;
            elseif strcmp(c_levels.conscLevel,'subseqbetween')
                c_levels.conscLevel = 3;
            elseif strcmp(c_levels.conscLevel,'subseqconsc')
                c_levels.conscLevel = 4;
            end
        end


        %enc_trial_type = strcmp(behavior.trial_type,'Enc');
        c_cond = behavior.rating == c_levels.conscLevel;
        c_acc = strcmp(behavior.accuracy,c_levels.acc);
        c_type = strcmp(behavior.trial_type,c_levels.cond);
        onsets{i} = behavior.onset(c_type...
            &c_cond...
            &c_acc);%-config.sessionOnsetCorrection;
        durations{i} = behavior.duration(c_type...
            &c_cond...
            &c_acc);
        length(onsets{i});
        %assert(~isempty(onsets{i}),['no onsets found for this condition: ', conds{i}])
        if (isempty(onsets{i}))
            warning(['no onsets found for this condition: ', conds{i}])
        end
    end
    
%% check output
    assert(all( cellfun(@(x) ischar(x),conds)))
    assert(all( cellfun(@(x) ischar(x),contrast)))
    assert(all( cellfun(@(x) isnumeric(x),convecs))) 
end