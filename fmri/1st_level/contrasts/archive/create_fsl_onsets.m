function [] = create_fsl_onsets(subject,special_case)

p=pwd;
idcs=strfind(p,'/');
addpath(fullfile(p(1:idcs(end)-1),'functions'))
addpath('contrasts')


base_config = setup_config(subject);
all_subs = base_config.seq.participants.participant_id;
for ind0 = 1:length(all_subs)
    subj = str2double(all_subs{ind0}(5:end));
    config = setup_config(subj);
    outnames ={'encoding','consolidation','retrieval30min','retrieval24h'};
    sessions ={'enc_subsequent','consolidation','ret','delret'};
    behavNames = {'enc_subsequent','consolidation','ret','delret'};
    for ind = 1:length(sessions)
        behaviorName = behavNames{ind};
        session = sessions{ind};
        outname = outnames{ind};
        behavior = readtable([config.data.behav.(behaviorName)],'FileType','delimitedtext'); % the sessions data is stored)
    
        outdir=fullfile('/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/utilities/fsl_onset_files/',outname,'/consc_onsets_twoEvents/',['sub-',num2str(subj)]);
        if ~isfolder(outdir)
            mkdir(outdir)
        end
    
        specify_contrasts = str2func(['contrasts_',session]);
        [conds, ~ , ~ , onsets, durations] = specify_contrasts(behavior);
        disp(length(conds))
        disp(length(onsets))
        disp(behaviorName)
        for ind2 = 1:length(conds)
            cond = conds{ind2};
            if special_case == 1
                special_cases = {'faceRecog','retRating_','recog_','recogRating_'};
                if contains(cond,special_cases{1}) || contains(cond,special_cases{2}) || contains(cond,special_cases{3}) || contains(cond,special_cases{4})
                    fname = [num2str(subj),'_',cond,'_','fsl_onset.txt'];
                    outpath = fullfile(outdir,fname);
                    if contains(cond,special_cases{1})
                        durations{ind2} = repelem(4.75,length(durations{ind2}))';
                    end
                    T = table(onsets{ind2},durations{ind2});
                    T.parammodul= repelem(1,height(T))';
                    writetable(T,outpath,'Delimiter',' ',"WriteVariableNames",0)
                else
                % Do Something else
                end
                
            end
            fname = [num2str(subj),'_',cond,'_','fsl_onset.txt'];
            outpath = fullfile(outdir,fname);
            T = table(onsets{ind2},durations{ind2});
            T.parammodul= repelem(1,height(T))';
            %writetable(T,outpath,'Delimiter',' ',"WriteVariableNames",0)  
        end
    end
end