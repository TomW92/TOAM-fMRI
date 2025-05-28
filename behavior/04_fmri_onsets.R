require(pacman)
library(viridisLite)
library(reshape2)
library(jsonlite)
library(tidyverse)
library(ggpubr)
library(hrbrthemes)
library(alluvial)
library(patchwork)
library(reshape2)
library(wesanderson)
library(factoextra)
source('~/GitHub/TOAM/behavior/archive/util.R')

# Script to extract and calculate the behavioral fMRI onsets

#check OS ----
if (.Platform$OS.type == 'windows'){
  volume <- 'Z:'
} else if (.Platform$OS.type == 'unix') {
  volume <- '/Volumes/psy_wfg_henke'
}

## create paths ----
basepath <- paste0(volume, "/s2019_twillems_fMRI_silent_engram/data/fMRI/")
raw <- "raw/"
bids <- "bids/"
behav <- "/behav/"
output <- 'derivatives/onsetfiles/'
ses1 <- "ses-1"
func <- "func"
subject_tsv <- "participants.tsv"


createOnsets <- function(sequence,group){
  if (grepl( 'hipp', sequence, fixed = TRUE) == TRUE){
    TR <- 2.12
  } else {
    TR <- 2.70
  }
  
  #load dataframe to delete onsets that were cut before due to different reasons such as to fast RT etc. multiipls options
  df_corrected <- read.csv(paste0(basepath,behav,"df_bothgroups",group,".csv"))

  #added this code so we delete onsets that we have no fMRI data for, because the acquisition ended prematurely
  wb_ons_del   <- read.csv(paste0(basepath,'utilities/wb_onsets_to_delete.csv'))
  hipp_ons_del   <- read.csv(paste0(basepath,'utilities/hipp_onsets_to_delete.csv'))
  ons_del <- rbind(wb_ons_del,hipp_ons_del)
  ons_del <- ons_del %>% select(c(id,item)) %>% distinct()
  for(ind in seq_len(nrow(ons_del))) {
    del_item <- ons_del[ind,]$item
    del_id <- ons_del[ind,]$id
    df_corrected <- df_corrected %>% filter(!((item %in% del_item + id %in% del_id ) == 2))
  }

  # set wd so output files are saved in the right directory with the right name According to the BIDS standard
  outputdir <- paste0(basepath, sequence, bids)
  setwd(outputdir)
  
  # Read the participant.TSV file which is an output of the bits coiner in order to have to participate names for output and file-naming
  subjects <- paste0(basepath, sequence, bids, subject_tsv)
  if (grepl( 'wb', sequence, fixed = TRUE) == TRUE){
    subjects <- paste0(basepath, sequence, bids, 'participants.tsv')
  }
  subjects <- read_tsv(subjects, col_names = TRUE)
  subjects <- subjects$participant_id
  
  #Create  empty data frames
  df_detailed <- data.frame()
  df_obj <- data.frame()
  df_recogbutton <- data.frame()
  
  #loop through each subject within one of the sequences
  for (subj in subjects){
    fullpath <- paste0(basepath, sequence, raw, subj, behav)
    files <- list.files(path=fullpath, no.. = TRUE)
    
   
    # Through all behavioural output files for this subject
    for (file in files){
      print(paste("reading: ",file))
      
      
      #get encoding files ----
      if (grepl("-encoding_run1_logfile.txt",file, fixed = TRUE) == TRUE){
        full_file <- paste0(fullpath, file)
        enc <- read.delim(full_file,)
        print(paste("wrangling: ",full_file))
        
        if (subj == "sub-65399"){
          print("stopheree")
        }
        
        #This saves the time of the first pulse the mri systems sends out. the very first pulse, even before the 5 dummy scans
        firstPulseTime <- enc$First_pulse_time
        
        #now take the encoding file and select the relevant columns
        enc_behav <- enc %>%
          select(ID,
                 Presented_Person,
                 FaceBase_OnsetTime,
                 #add duration
                 FaceBase_fMRI_Pulse,
                 FaceBase_fMRI_PulseTime, #sample to be delivered in pulses, or better the exact pulses as calculated with marc in the matlab script
                 #add trialtype
                 #add dummy response time
                 #add HED
                 Assoc_Encoding_OnsetTime,
                 #Add duration
                 Assoc_Encoding_fMRI_Pulse,
                 Assoc_Encoding_fMRI_PulseTime,
                 #add Trial_type
                 #Add dummy response_time (NA)
                 #add HED
                 Immediate_Recall_OnsetTime,
                 #add duration
                 Immediate_Recall_fMRI_Pulse,
                 Immediate_Recall_fMRI_PulseTime,
                 #trial_type
                 Immediate_Recall_RT,
                 Immediate_Recall_Accuracy,
                 
                 Immediate_Recall_Feedback_OnsetTime,
                 Immediate_Recall_Feedback_fMRI_Pulse,
                 Immediate_Recall_Feedback_fMRI_PulseTime) %>%
          
          #rename the relevant columns for reshaping
          rename(id=ID,item=Presented_Person,
                 
                 onset.Face = FaceBase_OnsetTime,
                 sample.Face = FaceBase_fMRI_Pulse,
                 PulseTime.Face = FaceBase_fMRI_PulseTime,
                 
                 onset.Enc = Assoc_Encoding_OnsetTime,
                 sample.Enc = Assoc_Encoding_fMRI_Pulse,
                 PulseTime.Enc = Assoc_Encoding_fMRI_PulseTime,
                 
                 onset.Recall = Immediate_Recall_OnsetTime,
                 sample.Recall = Immediate_Recall_fMRI_Pulse,
                 PulseTime.Recall = Immediate_Recall_fMRI_PulseTime,
                 response_time.Recall = Immediate_Recall_RT, 
                 accuracy.Recall = Immediate_Recall_Accuracy,
                 
                 onset.Feedback = Immediate_Recall_Feedback_OnsetTime,
                 sample.Feedback = Immediate_Recall_Feedback_fMRI_Pulse,
                 PulseTime.Feedback = Immediate_Recall_Feedback_fMRI_PulseTime) %>% 
          
          #mutate the columns to conform to bids format 
          mutate(duration.Face = 1,
                 response_time.Face = NA,
                 accuracy.Face = NA,
                 
                 duration.Enc = 2.75,
                 response_time.Enc = NA,
                 accuracy.Enc = NA,
                 
                 response_time.Recall = response_time.Recall/1000,
                 duration.Recall = 0,
                 
                 duration.Feedback = 0, 
                 response_time.Feedback = NA,
                 accuracy.Feedback = NA,
          )


        enc_info <- df_corrected %>% filter(id==substr(subj,5,9)) %>% select(c("item","first.retrieval","first.RT"))

        enc_behav <- enc_behav %>%
          #then delete trials that are not relevant anymore
          filter(item %in% unique(enc_info$item))  %>% right_join(y=enc_info,by="item") %>%
          #now, replace the answers, that were corrected in post and are not misses anymore.
          mutate(accuracy.Recall=first.retrieval,
                 response_time.Recall = first.RT) %>%
            select(-c("first.retrieval","first.RT")) %>%

          #change orientation
          pivot_longer(
            cols = 3:26,
            names_to=c(".value","trial_type"),
            names_sep = ('[.]')) %>% 
          
          mutate(
            onset = onset/1000,  #turn into seconds
            PulseTime = PulseTime/1000, #turn into seconds
            round_sample = sample,
            sample = sample + (onset - PulseTime)/TR, #calculate the "exact" pulse as onsets.
            onsets_sec = onset, #export onsets in seconds
            onset = onsets_sec - enc[1,]$First_pulse_time/1000
          ) %>% 
          
          select(-c(PulseTime,onsets_sec))

        #save ouput ----
        fname <- paste0(subj, "/", "ses-1/func/", subj, "_ses-1", '_task-forgetexp_events',group,'.tsv')
        print(paste("reading: ", as.character(file)))
        #write.table(enc_behav, file = fname, row.names=FALSE, sep="\t")
        print(paste("written: ",as.character(subj)))

        buttonpress <- enc_behav %>% filter(trial_type=="Recall")
        #write.table(buttonpress, file=paste0("/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/utilities/eventFiles/enc_buttonpress/",subj,'_enc_immrecall_onsets',group,'.txt'))
        
        #for fsl
        onset_types <- c("STEs","trial_types")
        for (onset_type in onset_types){
          fsl_outdir <- paste0(basepath,'utilities/fsl_onset_files/','encoding/',onset_type,'/',subj)
          dir.create(file.path(fsl_outdir), showWarnings = FALSE,recursive = TRUE)
          create_fsl_onsets(enc_behav,'enc',onset_type,fsl_outdir,group)
        }

        #for spm
        face_onsets <- enc_behav[enc_behav$trial_type=="Face",]$sample
        enc_onsets <- enc_behav[enc_behav$trial_type=="Enc",]$sample
        recall_onsets <- enc_behav[enc_behav$trial_type=="Recall",]$sample
        feedb_onsets <- enc_behav[enc_behav$trial_type=="Feedback",]$sample
        spm_onsets <- tibble(face_onsets, enc_onsets, recall_onsets, feedb_onsets)
        ##write.table(spm_onsets, file = paste0("/Users/tomwillems/Downloads/fmri/spm/",subj,'_task-forgetexp_events_SPM.txt'))

        
        memory <- enc_behav[enc_behav$trial_type=="Recall",]$accuracy
        single_trial_estimates_enc <- enc_behav %>% filter(trial_type == "Enc") %>% mutate(
          trial_type = item,
          accuracy = memory
        )
        ##write.table(single_trial_estimates_enc, file = paste0("/Users/tomwillems/Downloads/fmri/single_trial_events/",subj,'_task-forgetexp_events_SINGLETRIAL.txt'))
      }
      
      
      
      #get retrieval ----
      if (grepl("-retrieval_run1_logfile.txt",file, fixed = TRUE) == TRUE){
        full_file <- paste0(fullpath, file)
        print(paste("wrangling: ",full_file))
        ret <- read.delim(full_file,)
        fifth_pulse <- ret$num_pulse[1]
        fifth_pulseTime <- ret$fifth_pulse_time[1]
        ret_behav <- ret %>%
          select(ID,
                 Presented_Person,
                 FaceRecog_OnsetTime,
                 #add duration
                 FaceRecog_fMRI_Pulse,
                 FaceRecog_fMRI_PulseTime, #sample to be delivered in pulses, or better the exact pulses as calculated with marc in the matlab script
                 #add trialtype
                 #add dummy response time
                 #add HED
                 Assoc_Retrieval_OnsetTime,
                 #Add duration
                 Assoc_Retrieval_fMRI_Pulse,
                 Assoc_Retrieval_fMRI_PulseTime,
                 #add Trial_type
                 Assoc_Retrieval_RT,
                 Assoc_Retrieval_Accuracy,
                 #add HED
                 Assoc_Retrieval_Feedback_OnsetTime,
                 #add duration
                 Assoc_Retrieval_Feedback_fMRI_Pulse,
                 Assoc_Retrieval_Feedback_fMRI_PulseTime,
                 #trial_type
                 Retrieval_consc_Rating_OnsetTime,
                 Retrieval_consc_Rating_fMRI_Pulse,
                 Retrieval_consc_Rating_fMRI_PulseTime,
                 Retrieval_consc_Rating_RT,
                 Retrieval_consc_Rating_Pressed_Button,
                 
                 Retrieval_consc_Rating_Feedback_OnsetTime,
                 Retrieval_consc_Rating_Feedback_fMRI_Pulse,
                 Retrieval_consc_Rating_Feedback_fMRI_PulseTime) %>%
          
          rename(id=ID,
                 item=Presented_Person,
                 
                 onset.faceRecog = FaceRecog_OnsetTime,
                 sample.faceRecog = FaceRecog_fMRI_Pulse,
                 PulseTime.faceRecog = FaceRecog_fMRI_PulseTime,
                 
                 onset.ret = Assoc_Retrieval_OnsetTime,
                 sample.ret = Assoc_Retrieval_fMRI_Pulse,
                 PulseTime.ret = Assoc_Retrieval_fMRI_PulseTime,
                 response_time.ret = Assoc_Retrieval_RT,
                 accuracy.ret = Assoc_Retrieval_Accuracy,
                 
                 onset.retFeedback = Assoc_Retrieval_Feedback_OnsetTime,
                 sample.retFeedback = Assoc_Retrieval_Feedback_fMRI_Pulse,
                 PulseTime.retFeedback = Assoc_Retrieval_Feedback_fMRI_PulseTime,
                 
                 onset.retRating = Retrieval_consc_Rating_OnsetTime,
                 sample.retRating = Retrieval_consc_Rating_fMRI_Pulse,
                 PulseTime.retRating = Retrieval_consc_Rating_fMRI_PulseTime,
                 response_time.retRating = Retrieval_consc_Rating_RT, 
                 accuracy.retRating = Retrieval_consc_Rating_Pressed_Button,
                 
                 onset.retRatingFeedback = Retrieval_consc_Rating_Feedback_OnsetTime,
                 sample.retRatingFeedback = Retrieval_consc_Rating_Feedback_fMRI_Pulse,
                 PulseTime.retRatingFeedback = Retrieval_consc_Rating_Feedback_fMRI_PulseTime) %>% 
          
          mutate(duration.faceRecog = 1,
                 response_time.faceRecog = NA,
                 accuracy.faceRecog = NA,
                 
                 accuracy.ret = case_when(accuracy.ret == "correct" ~ "correct",
                                          accuracy.ret == "incorrect" ~ "incorrect",
                                          TRUE ~ "miss"),  #here and below with the delayed retrieval, I delete "MISS" answers
                 duration.ret = 2.75,
                 response_time.ret = response_time.ret/1000,
                 
                 duration.retFeedback = 0,
                 response_time.retFeedback = NA,
                 accuracy.retFeedback = NA,
                 
                 duration.retRating = 2.75,
                 accuracy.retRating = as.character(accuracy.retRating),
                 response_time.retRating = response_time.retRating/1000,
                 
                 duration.retRatingFeedback = 0, 
                 response_time.retRatingFeedback = NA)

          #filter out missed trials, as they would lead to more noisy data.
          #first select the important columns.
          ret_info <- df_corrected %>% filter(id==substr(subj,5,9)) %>% select(c("item","second.retrieval","second.rating","second.RT"))

          ret_behav <- ret_behav %>%
          #then delete trials that are not relevant anymore
          filter(item %in% unique(ret_info$item))  %>% right_join(y=ret_info,by="item") %>%
          #now, replace the answers, that were corrected in post and are not misses anymore.
          mutate(accuracy.ret=second.retrieval,
                 accuracy.retRating=as.character(second.rating),
                 response_time.ret = second.RT,
                 accuracy.retRatingFeedback = case_when(as.integer(accuracy.retRating) == 4 ~ "conscious",  #accuracy.ret == "correct" &  was here before, but for the sake of just doing an consc vs unconsc we leave the accuracy
                                                        as.integer(accuracy.retRating) == 3 ~ "neither",
                                                        as.integer(accuracy.retRating) == 2 ~ "unconscious",
                                                          TRUE ~ "miss") ) %>%
            select(-c("second.retrieval","second.rating","second.RT")) %>%

          pivot_longer(
            cols = 3:32,
            names_to=c(".value","trial_type"),
            names_sep = ('[.]')) %>% 
          
          # ENTSCHEIDENDER PUNKT HIER, BERECHNUNG DER SAMPLES/ONSETS/SCANS/SEKUNDEN ####
          mutate(
            onset = onset/1000,  #turn into seconds
            PulseTime = PulseTime/1000, #turn into seconds
            round_sample = sample - ( fifth_pulse -5 ) , #pulse counted forward and I saved the fifth pulse, then substract 5 again
            sample = round_sample + (onset - PulseTime)/TR, #calculate the "exact" pulse as onsets.
            onsets_orig = onset, #export onsets in seconds
            onset = onsets_orig - ((ret$fifth_pulse_time[1]/1000)-TR*4), #s only subtract 4 times, at the moment of the 5th pulse, only 4*TR has happened, not 5
            onset = round(onset, 2)
          ) %>% 
          filter_at(vars(onset), all_vars(!is.na(.)))  %>% # pc crashed once, missing last 10 trials, filter those out here (sub-65399)
          select(-c(PulseTime,onsets_orig))

        ratings <- ret_behav %>%
          filter(trial_type == "retRatingFeedback") %>%
          select(accuracy) %>%
          mutate(consc_level = accuracy) %>%
          select(-accuracy)

        consc_scores <- ret_behav %>%
          filter(trial_type=="ret") %>%
          mutate(consc_level = ratings,
                 pmod = case_when(
             accuracy == "incorrect" & consc_level == "unconscious" ~ "_guess_incorrect" ,
             accuracy == "correct" & consc_level == "unconscious" ~ "_guess_correct",
             accuracy == "incorrect" & consc_level == "neither" ~ "_unsure_incorrect",
             accuracy == "correct" & consc_level == "neither" ~ "_unsure_correct",
             accuracy == "incorrect" & consc_level == "conscious" ~ "_sure_incorrect",
             accuracy == "correct" & consc_level == "conscious" ~ "_sure_correct"))
        
        ret_behav <- ret_behav %>% mutate(trialtype_new = rep(consc_scores$pmod, each=length(unique(trial_type))))

        bids_ret = ret_behav %>% 
          mutate(
            trial_type = paste0(trial_type,trialtype_new),
            duration = 0,) %>% 
          select(-c(sample,id,accuracy,round_sample,trialtype_new,response_time)) %>% 
          subset(select=c(3,4,2,1))
        
        fname <- paste0(subj, "/", "ses-1/func/", subj, "_ses-1", '_task-retrieval_events',group,'.tsv')
        if (group == "_maxTrials"){
          fname <- paste0(subj, "/", "ses-1/func/", subj, "_ses-1", '_task-retrieval_events.tsv')
          write.table(bids_ret, file = fname, row.names=FALSE, sep="\t", quote = FALSE)
        }

        # save ouput ----
        fname <- paste0(subj, "/", "ses-1/func/", subj, "_ses-1", '_task-retrieval_events',group,'.tsv')
        print(paste("reading: ", as.character(file)))
        #write.table(ret_behav, file = fname, row.names=FALSE, sep="\t")
        #write.csv(ret_behav, file = paste0("/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/utilities/eventFiles/retrieval_all/",subj,'_retrieval_all_events',group,'.csv'))
        print(paste("written: ", as.character(subj), "to: ", fname))

        #ret_behav <- ret_behav %>% select(-pmod)

        #onsets for picture viewing vs no pic
        target <- c("faceRecog", "ret", "retFeedback")
        pictureOnsets <- ret_behav %>%  filter(trial_type %in% target) %>% mutate(trial_type = "pictureOnScreen")
        #write.table(pictureOnsets, file=paste0("/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/utilities/eventFiles/pictureOnsets/",subj,'_pictureOnsets',group,'.txt'))
        
        #onsets for scripted nilearn 
        conscinfo <- ret_behav %>% filter(trial_type == "retRatingFeedback") %>% select(accuracy)
        retConscEvents <- ret_behav %>% filter(trial_type == "ret") %>% mutate(trial_type = conscinfo$accuracy)
        #write.table(retConscEvents, file=paste0("/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/utilities/eventFiles/consc_unconsc_neither_retCOMBINED/",subj,'_ret_combined',group,'.txt'))

        #for fsl
        onset_types <- c("STEs","trial_types")
        for (onset_type in onset_types){
          fsl_outdir <- paste0(basepath,'utilities/fsl_onset_files/','retrieval30min/',onset_type,'/',subj)
          dir.create(file.path(fsl_outdir), showWarnings = FALSE,recursive = TRUE)
          create_fsl_onsets(ret_behav,'faceRecog',onset_type,fsl_outdir,group)
        }

        #create spm onsets
        face_onsets_ret <- ret_behav[ret_behav$trial_type=="faceRecog",]$sample
        ret_onsets_ret <- ret_behav[ret_behav$trial_type=="ret",]$sample
        retfeed_onsets_ret <- ret_behav[ret_behav$trial_type=="retFeedback",]$sample
        rating_onsets_ret <- ret_behav[ret_behav$trial_type=="retRating",]$sample
        ratingfeed_onsets_ret <- ret_behav[ret_behav$trial_type=="retRatingFeedback",]$sample
        spm_onsets_ret <- tibble(face_onsets_ret, ret_onsets_ret, retfeed_onsets_ret, rating_onsets_ret, ratingfeed_onsets_ret)
        ##write.table(spm_onsets_ret, file = paste0("/Users/tomwillems/Downloads/fmri/spm/",subj,'_task-retrieval_events_SPM.txt'))
        
        #differnt spm_onsets
        unconsc_items <- ret_behav %>% filter(accuracy=="unconscious") %>% select(item)
        consc_items <- ret_behav %>% filter(accuracy=="conscious") %>% select(item)
        neither_items <- ret_behav %>% filter(accuracy=="neither") %>% select(item)
        unconsc_ret_trials <- ret_behav %>% filter(item %in% unconsc_items$item) %>% filter(trial_type=="ret")
        consc_ret_trials <- ret_behav %>% filter(item %in% consc_items$item) %>% filter(trial_type=="ret")
        neither_ret_trials <- ret_behav %>% filter(item %in% neither_items$item) %>% filter(trial_type=="ret")
        #write.table(unconsc_ret_trials, file = paste0("/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/utilities/eventFiles/consc_unconsc_neither_ret/",subj,'_ret_unconsc_ret',group,'.txt'))
        #write.table(consc_ret_trials,   file = paste0("/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/utilities/eventFiles/consc_unconsc_neither_ret/",subj,'_ret_consc_ret',group,'.txt'))
        #write.table(neither_ret_trials, file = paste0("/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/utilities/eventFiles/consc_unconsc_neither_ret/",subj,'_ret_neither_ret',group,'.txt'))
        
        memory <- ret_behav[ret_behav$trial_type=="retRatingFeedback",]$accuracy
        single_trial_estimates_ret <- ret_behav %>% filter(trial_type == "ret") %>% mutate(
          trial_type = item,
          accuracy = memory
        )
        ##write.table(single_trial_estimates_ret, file = paste0("/Users/tomwillems/Downloads/fmri/single_trial_events/",subj,'_task-retrieval_events_SINGLETRIAL.txt'))
        
      } 
      
      
      #get delayed ret ----
      if (grepl("delayed_retrieval_run1_logfile.txt",file, fixed = TRUE) == TRUE){
        full_file <- paste0(fullpath, file)
        print(paste("wrangling DELRET: ",full_file))
        delret <- read.delim(full_file)
        fifth_pulse_delret <- delret$num_pulse[1]
        fifth_pulseTime_delret <- delret$fifth_pulse_time[1]
        delret_behav <- delret %>%
          select(ID,
                 Presented_Person,
                 FaceRecog_OnsetTime,
                 #add duration
                 FaceRecog_fMRI_Pulse,
                 FaceRecog_fMRI_PulseTime, #sample to be delivered in pulses, or better the exact pulses as calculated with marc in the matlab script
                 #add trialtype
                 #add dummy response time
                 #add HED
                 Assoc_Retrieval_OnsetTime,
                 #Add duration
                 Assoc_Retrieval_fMRI_Pulse,
                 Assoc_Retrieval_fMRI_PulseTime,
                 #add Trial_type
                 Assoc_Retrieval_RT,
                 Assoc_Retrieval_Accuracy,
                 #add HED
                 Assoc_Retrieval_Feedback_OnsetTime,
                 #add duration
                 Assoc_Retrieval_Feedback_fMRI_Pulse,
                 Assoc_Retrieval_Feedback_fMRI_PulseTime,
                 #trial_type
                 Retrieval_consc_Rating_OnsetTime,
                 Retrieval_consc_Rating_fMRI_Pulse,
                 Retrieval_consc_Rating_fMRI_PulseTime,
                 Retrieval_consc_Rating_RT,
                 Retrieval_consc_Rating_Pressed_Button,
                 
                 Retrieval_consc_Rating_Feedback_OnsetTime,
                 Retrieval_consc_Rating_Feedback_fMRI_Pulse,
                 Retrieval_consc_Rating_Feedback_fMRI_PulseTime,
                 
                 ForcedChoice_OnsetTime,
                 #add duration
                 ForcedChoice_fMRI_Pulse,
                 ForcedChoice_fMRI_PulseTime,
                 ForcedChoice_Pressed_Button,
                 #trial_type
                 ForcedChoice_Accuracy,
                 ForcedChoice_RT,
                 
                 ForcedChoice_Feedback_OnsetTime,
                 ForcedChoice_Feedback_fMRI_Pulse,
                 ForcedChoice_Feedback_fMRI_PulseTime,
                 
                 Recognition_consc_Rating_OnsetTime,
                 Recognition_consc_Rating_fMRI_Pulse,
                 Recognition_consc_Rating_fMRI_PulseTime,
                 Recognition_consc_Rating_RT,
                 Recognition_consc_Rating_Pressed_Button,
                 
                 Recognition_consc_Rating_Feedback_OnsetTime,
                 Recognition_consc_Rating_Feedback_fMRI_Pulse,
                 Recognition_consc_Rating_Feedback_fMRI_PulseTime) %>% 
          
          rename(id=ID,
                 item=Presented_Person,
                 
                 onset.faceRecog = FaceRecog_OnsetTime,
                 sample.faceRecog = FaceRecog_fMRI_Pulse,
                 PulseTime.faceRecog = FaceRecog_fMRI_PulseTime,
                 
                 onset.ret = Assoc_Retrieval_OnsetTime,
                 sample.ret = Assoc_Retrieval_fMRI_Pulse,
                 PulseTime.ret = Assoc_Retrieval_fMRI_PulseTime,
                 response_time.ret = Assoc_Retrieval_RT,
                 accuracy.ret = Assoc_Retrieval_Accuracy,
                 
                 onset.retFeedback = Assoc_Retrieval_Feedback_OnsetTime,
                 sample.retFeedback = Assoc_Retrieval_Feedback_fMRI_Pulse,
                 PulseTime.retFeedback = Assoc_Retrieval_Feedback_fMRI_PulseTime,
                 
                 onset.retRating = Retrieval_consc_Rating_OnsetTime,
                 sample.retRating = Retrieval_consc_Rating_fMRI_Pulse,
                 PulseTime.retRating = Retrieval_consc_Rating_fMRI_PulseTime,
                 response_time.retRating = Retrieval_consc_Rating_RT,
                 accuracy.retRating = Retrieval_consc_Rating_Pressed_Button,
                 
                 onset.retRatingFeedback = Retrieval_consc_Rating_Feedback_OnsetTime,
                 sample.retRatingFeedback = Retrieval_consc_Rating_Feedback_fMRI_Pulse,
                 PulseTime.retRatingFeedback = Retrieval_consc_Rating_Feedback_fMRI_PulseTime,
                 
                 onset.recog = ForcedChoice_OnsetTime,
                 sample.recog = ForcedChoice_fMRI_Pulse,
                 PulseTime.recog = ForcedChoice_fMRI_PulseTime,
                 response_time.recog = ForcedChoice_RT,
                 accuracy.recog = ForcedChoice_Accuracy,

                 accuracy.recogFeedback = ForcedChoice_Pressed_Button, #in order to have the button presses set
                 onset.recogFeedback = ForcedChoice_Feedback_OnsetTime,
                 sample.recogFeedback = ForcedChoice_Feedback_fMRI_Pulse,
                 PulseTime.recogFeedback = ForcedChoice_Feedback_fMRI_PulseTime,

                 onset.recogRating = Recognition_consc_Rating_OnsetTime,
                 sample.recogRating = Recognition_consc_Rating_fMRI_Pulse,
                 PulseTime.recogRating = Recognition_consc_Rating_fMRI_PulseTime,
                 response_time.recogRating = Recognition_consc_Rating_RT,
                 accuracy.recogRating = Recognition_consc_Rating_Pressed_Button,
                 
                 onset.recogRatingFeedback = Recognition_consc_Rating_Feedback_OnsetTime,
                 sample.recogRatingFeedback = Recognition_consc_Rating_Feedback_fMRI_Pulse,
                 PulseTime.recogRatingFeedback = Recognition_consc_Rating_Feedback_fMRI_PulseTime) %>% 
          
          mutate(duration.faceRecog = 1,
                 response_time.faceRecog = NA,
                 accuracy.faceRecog = NA,
                 
                 accuracy.ret = case_when(accuracy.ret == "correct" ~ "correct",
                                          accuracy.ret == "incorrect" ~ "incorrect",
                                          TRUE ~ "incorrect"),
                 duration.ret = 2.75, #response_time.ret/1000,
                 response_time.ret = response_time.ret/1000,
                 
                 
                 duration.retFeedback = 0,
                 response_time.retFeedback = NA,
                 accuracy.retFeedback = NA,
                 
                 duration.retRating = 2.75, #response_time.retRating/1000,
                 accuracy.retRating = as.character(accuracy.retRating),
                 response_time.retRating = response_time.retRating/1000,
                 
                 duration.retRatingFeedback = 0, 
                 response_time.retRatingFeedback = NA,
                 accuracy.retRatingFeedback = case_when(as.integer(accuracy.retRating) == 4 ~ "conscious",  #accuracy.ret == "correct" &  was here before, but for the sake of just doing an consc vs unconsc we leave the accuracy
                                                        as.integer(accuracy.retRating) == 3 ~ "neither",  ##if_else((accuracy.ret == "correct" & as.integer(accuracy.retRating) > 2),"remembered","forgotten"),
                                                        as.integer(accuracy.retRating) == 2 ~ "unconscious",
                                                        TRUE ~ "unconscious"),
                 

                 
                 duration.recog = 3, #response_time.recog/1000,
                 response_time.recog = response_time.recog/1000,
                 accuracy.recog = case_when(accuracy.recog == "correct" ~ "correct",
                                          accuracy.recog == "incorrect" ~ "correct",
                                          TRUE ~ "incorrect"),
                 
                 duration.recogFeedback = 0,
                 response_time.recogFeedback = NA,
                 accuracy.recogFeedback = if_else(accuracy.recogFeedback==2,'index','middle'),
                 
                 duration.recogRating = 2.75, #response_time.retRating/1000,
                 accuracy.recogRating = as.character(accuracy.recogRating),
                 response_time.recogRating = response_time.retRating/1000,
                 
                 duration.recogRatingFeedback = 0, #duration.recogRating - (response_time.recogRating), 
                 response_time.recogRatingFeedback = NA) 


          #filter out missed trials, as they would lead to more noisy data.
          #first select the important columns.
          delret_info <- df_corrected %>% filter(id==substr(subj,5,9)) %>% select(c("item","third.retrieval","third.rating","third.RT","fourth.retrieval","fourth.rating","fourth.RT"))

          delret_behav <- delret_behav %>%
          #then delete trials that are not relevant anymore
          filter(item %in% unique(delret_info$item))  %>% right_join(y=delret_info,by="item") %>%
          #now, replace the answers, that were corrected in post and are not misses anymore.
          mutate(accuracy.ret=third.retrieval,
                 accuracy.retRating=as.character(third.rating),
                 response_time.ret = third.RT,
                 accuracy.retRatingFeedback = case_when(as.integer(accuracy.retRating) == 4 ~ "conscious",  #accuracy.ret == "correct" &  was here before, but for the sake of just doing an consc vs unconsc we leave the accuracy
                                                        as.integer(accuracy.retRating) == 3 ~ "neither",
                                                        as.integer(accuracy.retRating) == 2 ~ "unconscious",
                                                          TRUE ~ "miss"),
                 accuracy.recog = fourth.retrieval,
                 accuracy.recogRating = as.character(fourth.rating),
                 response_time.recog = fourth.RT,
                 accuracy.recogRatingFeedback = case_when(as.integer(accuracy.recogRating) == 4 ~ "conscious",  #accuracy.ret == "correct" &  was here before, but for the sake of just doing an consc vs unconsc we leave the accuracy
                                                          as.integer(accuracy.recogRating) == 3 ~ "neither",  ##if_else((accuracy.ret == "correct" & as.integer(accuracy.retRating) > 2),"remembered","forgotten"),
                                                          as.integer(accuracy.recogRating) == 2 ~ "unconscious",
                                                          TRUE ~ "unconscious")) %>%
            
          select(-c("third.retrieval","third.rating","third.RT","fourth.retrieval","fourth.rating","fourth.RT")) %>%

          pivot_longer(cols = 3:56,
                       names_to=c(".value","trial_type"),
                       names_sep = ('[.]')) %>%

          # ENTSCHEIDENDER PUNKT HIER, BERECHNUNG DER SAMPLES/ONSETS/SCANS/SEKUNDEN ####
          mutate(onset = onset/1000,  #turn into seconds
                 PulseTime = PulseTime/1000, #turn into seconds
                 round_sample = sample,
                 sample = sample + (onset - PulseTime)/TR, #calculate the "exact" pulse as onsets.
                 onsets_sec = onset, #export onsets in seconds
                 onset = onsets_sec - delret$First_Pulse_Time[1]/1000,
                 #onset = sample + (onset - PulseTime)/TR, #calculate the "exact" pulse as onsets.
                 onset = round(onset, 2)
          ) %>%
          select(-c(PulseTime,onsets_sec))

        ratings <- delret_behav %>%
          filter(trial_type == "retRatingFeedback") %>%
          select(accuracy) %>%
          mutate(consc_level = accuracy) %>%
          select(-accuracy)
        
        recogRatings = delret_behav %>%
          filter(trial_type == "recogRatingFeedback") %>%
          select(accuracy) %>%
          mutate(consc_level = accuracy) %>%
          select(-accuracy)

        consc_scores <- delret_behav %>%
          filter(trial_type=="ret") %>%
          mutate(consc_level = ratings,
                 pmod = case_when(
                   accuracy == "incorrect" & consc_level == "unconscious" ~ "_guess_incorrect" ,
                   accuracy == "correct" & consc_level == "unconscious" ~ "_guess_correct",
                   accuracy == "incorrect" & consc_level == "neither" ~ "_unsure_incorrect",
                   accuracy == "correct" & consc_level == "neither" ~ "_unsure_correct",
                   accuracy == "incorrect" & consc_level == "conscious" ~ "_sure_incorrect",
                   accuracy == "correct" & consc_level == "conscious" ~ "_sure_correct"))
       
        delret_only = delret_behav %>% filter(trial_type == "faceRecog" | trial_type == "ret" | trial_type == "retFeedback" | trial_type == "retRating" | trial_type == "retRatingFeedback")
        delret_only_behav <- delret_only %>% mutate(trialtype_new = rep(consc_scores$pmod, each=length(unique(trial_type))))

        consc_scores <- delret_behav %>%
          filter(trial_type=="recog") %>%
          mutate(consc_level = recogRatings,
                 pmod = case_when(
                   accuracy == "incorrect" & consc_level == "unconscious" ~ "_guess_incorrect" ,
                   accuracy == "correct" & consc_level == "unconscious" ~ "_guess_correct",
                   accuracy == "incorrect" & consc_level == "neither" ~ "_unsure_incorrect",
                   accuracy == "correct" & consc_level == "neither" ~ "_unsure_correct",
                   accuracy == "incorrect" & consc_level == "conscious" ~ "_sure_incorrect",
                   accuracy == "correct" & consc_level == "conscious" ~ "_sure_correct"))
        
        recog_only = delret_behav %>% filter(trial_type == "recog" | trial_type == "recogFeedback" | trial_type == "recogRating" | trial_type == "recogRatingFeedback" )
        recog_only_behav <- recog_only %>% mutate(trialtype_new = rep(consc_scores$pmod, each=length(unique(trial_type))))
        
        bids_delret = rbind(delret_only_behav,recog_only_behav)
        bids_delret = bids_delret[order(bids_delret$onset), ]
        bids_delret = bids_delret %>% 
          mutate(
            trial_type = paste0(trial_type,trialtype_new),
            duration = 0,) %>% 
          select(-c(sample,id,accuracy,round_sample,trialtype_new,response_time)) %>% 
          subset(select=c(3,4,2,1))
        
        if (group == "_maxTrials"){
          fname <- paste0(subj, "/", "ses-2/func/", subj, "_ses-2", '_task-forgetexp_events.tsv')
          write.table(bids_delret, file = fname, row.names=FALSE, sep="\t", quote = FALSE)
        }

        #save ouput ----
        fname <- paste0(subj, "/", "ses-2/func/", subj, "_ses-2", '_task-forgetexp_events',group,'.tsv')
        print(paste("reading: ", as.character(file)))
        #write.table(delret_behav, file = fname, row.names=FALSE, sep="\t")
        #write.csv(delret_behav, file = paste0("/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/utilities/eventFiles/delayed_retrieval_all/",subj,'_delayed_retrieval_all_events',group,'.csv'))
        print(paste("written: ",as.character(subj), "to: ", fname))

        #delret_behav <- delret_behav %>% select(-pmod)

        #for fsl
        onset_types <- c("STEs","trial_types")
        STE_types <- c("faceRecog","recog")
        for (onset_type in onset_types){
          for (STE_type in STE_types){
            if (onset_type == "STEs"){
              fsl_outdir <- paste0(basepath,'utilities/fsl_onset_files/retrieval24h/',onset_type,'/',STE_type,'/',subj)
            } else {
              fsl_outdir <- paste0(basepath,'utilities/fsl_onset_files/retrieval24h/',onset_type,'/',subj)
            }
            dir.create(file.path(fsl_outdir), showWarnings = FALSE,recursive = TRUE)
            create_fsl_onsets(delret_behav,STE_type,onset_type,fsl_outdir,group)
          }
        }

        
        #SPM onsets
        face_onsets_delret <- delret_behav[delret_behav$trial_type=="faceRecog",]$sample
        ret_onsets_delret <- delret_behav[delret_behav$trial_type=="ret",]$sample
        retfeed_onsets_delret <- delret_behav[delret_behav$trial_type=="retFeedback",]$sample
        rating_onsets_delret <- delret_behav[delret_behav$trial_type=="retRating",]$sample
        ratingfeed_onsets_delret <- delret_behav[delret_behav$trial_type=="retRatingFeedback",]$sample
        recog_onset_delret <- delret_behav[delret_behav$trial_type=="recog",]$sample
        recogfeed_onset_delret <- delret_behav[delret_behav$trial_type=="recogFeedback",]$sample
        recograting_onset_delret <- delret_behav[delret_behav$trial_type=="recogRating",]$sample
        recogratingfeed_onset_delret <- delret_behav[delret_behav$trial_type=="recogRatingFeedback",]$sample
        spm_onsets_delret <- tibble(face_onsets_delret, ret_onsets_delret, retfeed_onsets_delret, rating_onsets_delret, ratingfeed_onsets_delret, recog_onset_delret,
                                        recogfeed_onset_delret, recograting_onset_delret, recogratingfeed_onset_delret)
        ##write.table(spm_onsets_delret, file = paste0("/Users/tomwillems/Downloads/",subj,'_task-delayed_retrieval_events_SPM.txt'))
        
        delret_subseq <- delret_behav %>% pivot_wider(id_cols = c(id,item),values_from = c(onset,sample,accuracy),names_from = trial_type)
        delret_subseq <- delret_subseq %>% mutate(accuracy_faceRecog = case_when(
              accuracy_retRating == "4" ~ 'consc',
              accuracy_retRating == "3" ~ 'neither',
              accuracy_retRating == "2" & accuracy_recogRating == "4" ~ 'later_consc',
              accuracy_retRating == "2" & accuracy_recogRating != "4" ~ 'stays_unconsc',
              TRUE ~  'undefined'
            ))
        delret_subseq <- delret_subseq %>% pivot_longer(cols=3:29,names_to=c(".value",'trial_type'),names_sep=('[_]'))
        #write.csv(delret_subseq, file = paste0("/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/utilities/eventFiles/subs_memory_all/",subj,'_delret_subsequent_consc',group,'.csv'))




        #delret consc vs unconsc
        unconsc_items <- delret_behav %>% filter(trial_type != "recogRatingFeedback") %>%  filter(accuracy=="forgotten") %>% select(item)
        consc_items <- delret_behav %>% filter(trial_type != "recogRatingFeedback") %>%  filter(accuracy=="remembered") %>% select(item)
        unconsc_delret_trials <- delret_behav %>% filter(item %in% unconsc_items$item) %>% filter(trial_type=="faceRecog")
        consc_delret_trials <- delret_behav %>% filter(item %in% consc_items$item) %>% filter(trial_type=="faceRecog")
        #write.table(unconsc_delret_trials, file = paste0("/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/utilities/eventFiles/consc_unconsc_del_ret/",subj,'_delret_unconsc_faceBase',group,'.txt'))
        #write.table(consc_delret_trials,   file = paste0("/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/utilities/eventFiles/consc_unconsc_del_ret/",subj,'_delret_consc_faceBase',group,'.txt'))
        
        #del recog consc vs unconsc
        unconsc_items <- delret_behav %>% filter(trial_type != "retRatingFeedback") %>%  filter(accuracy=="forgotten") %>% select(item)
        consc_items <- delret_behav %>% filter(trial_type != "retRatingFeedback") %>%  filter(accuracy=="remembered") %>% select(item)
        unconsc_delret_trials <- delret_behav %>% filter(item %in% unconsc_items$item) %>% filter(trial_type=="faceRecog")
        consc_delret_trials <- delret_behav %>% filter(item %in% consc_items$item) %>% filter(trial_type=="faceRecog")
        #write.table(unconsc_delret_trials, file = paste0("/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/utilities/eventFiles/consc_unconsc_del_recog/",subj,'_delret_unconsc_faceBase',group,'.txt'))
        #write.table(consc_delret_trials,   file = paste0("/Volumes/psy_wfg_henke/s2019_twillems_fMRI_silent_engram/data/fMRI/utilities/eventFiles/consc_unconsc_del_recog/",subj,'_delret_consc_faceBase',group,'.txt'))
      }
    } #FILES
    
    
    #Subsequenct Memory onsets.
    subsec_enc <- enc_behav %>% select(id, item, trial_type, onset, sample, response_time) #%>% filter(trial_type == "Enc")
    
    subsec_ret <- ret_behav %>% filter(trial_type == "ret") %>% select(id, item, accuracy) #%>% mutate(conscious = if_else(accuracy > 3, TRUE, FALSE))
    subsec_recog <- ret_behav %>%  filter(trial_type == "retRating" ) %>%  select(id, item, accuracy) %>% mutate(rating = accuracy) %>% select(-accuracy)
    subsec_ret <- right_join(subsec_ret, subsec_recog, by = c("id", "item"))
    subsequent_mem_onsets <- right_join(subsec_enc, subsec_ret, by=c("id", "item"))
    subsequent_mem_onsets <- subsequent_mem_onsets %>% mutate(
      duration = case_when(
        trial_type == "Face" ~ 1,
        trial_type == "Enc" ~ 2.5,
        trial_type == "Recall" ~ 2.75,
        trial_type == "RecallFeedback" ~ 0,
        TRUE ~ 0,
      ))
    
    #bids_output:
    bids_subsequent_enc = subsequent_mem_onsets %>% 
      mutate(
        trial_type = case_when(
          accuracy=="correct" & rating=="4" ~ paste(trial_type,"_sure_correct",sep=""),
          accuracy=="correct" & rating=="3" ~  paste(trial_type,"_unsure_correct",sep=""),
          accuracy=="correct" & rating=="2" ~  paste(trial_type,"_guess_correct",sep=""),
          accuracy=="incorrect" & rating=="4" ~ paste(trial_type,"_sure_incorrect",sep=""),
          accuracy=="incorrect" & rating=="3" ~ paste(trial_type,"_unsure_incorrect",sep=""),
          accuracy=="incorrect" & rating=="2" ~ paste(trial_type,"_guess_incorrect",sep="")),
        duration = 0,) %>% 
      select(-c(sample,id,accuracy,rating,response_time)) %>% 
      mutate(onset = round(onset, 2)) %>% 
      subset(select=c(3,4,2,1))
    
    fname <- paste0(subj, "/", "ses-1/func/", subj, "_ses-1", '_task-retrieval_events',group,'.tsv')
    if (group == "_maxTrials"){
      fname <- paste0(subj, "/", "ses-1/func/", subj, "_ses-1", '_task-forgetexp_events.tsv')
      write.table(bids_subsequent_enc, file = fname, row.names=FALSE, sep="\t", quote = FALSE)
      
      fname <- paste0(subj, "/", "ses-1/func/", subj, "_ses-1", '_task-forgetexp_events',group,'.tsv')
      print(paste("reading: ", as.character(file)))
      #write.table(enc_behav, file = fname, row.names=FALSE, sep="\t")
    }
    print("encoding subsequent")
  } #SUBJECTS
print("next Subject")
} #FUNCTION

sequences <- list("hipp/", "wb/")
groups <- list("_corrected","_maxTrials")
# main loop over sequences -----
for (sequence in sequences){
  for (group in groups){
      createOnsets(sequence,group)
  }
}
