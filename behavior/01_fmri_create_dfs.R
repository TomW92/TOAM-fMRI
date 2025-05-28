require(pacman)
library(jsonlite)
library(tidyverse)
library(patchwork)
library(reshape2)

#check OS ----
if (Sys.info()['sysname']=='Linux'){
  volume <- '/storage/research/psy_wfg_henke'
} else if (Sys.info()['sysname']=='Darwin') {
  volume <- '/Volumes/psy_wfg_henke'
} else {
  volume <- 'Z:' #or the letter your windows mounted the research storage under
}

## create paths ----
basepath <- paste0(volume, "/s2019_twillems_fMRI_silent_engram/data/fMRI/")
sequences <- list("hipp/", "wb/")
raw <- "raw/"
bids <- "bids/"
behav <- "/behav/"
output <- 'derivatives/behav/'
ses1 <- "ses-1"
func <- "func"
subject_tsv <- "participants.tsv"
specpath <- "sub-"

#create empty dataframes for behavioral data from both acquisitions
df_detailed_hipp <- data.frame()
df_detailed_wb <- data.frame()
df_detailed_all <- data.frame()

#now go through the files for one sequences
for (sequence in sequences){
  print(sequence)
  ## get participant names with the participants file from bids ----
  subjects <- paste0(basepath, sequence, bids, subject_tsv)
  subjects <- read_tsv(subjects, col_names = TRUE, col_types = cols())
  subjects <- subjects$participant_id
  df_detailed <- data.frame()
  
  for (subj in subjects){
    fullpath <- paste0(basepath, sequence, raw, subj, behav)
    files <- list.files(path=fullpath, no.. = TRUE)


    recogCorrect <- read.delim(paste0(fullpath,files[grepl("delayed_retrieval_object_pos_run1",files)]),header=FALSE)

    for (file in files){
      #get encoding
      if (grepl("-encoding_run1_logfile.txt",file, fixed = TRUE) == TRUE){
        full_file <- paste0(fullpath, file)
        enc <- read.delim(full_file,)
        enc_behav <- enc %>% select(ID, Presented_Person, Immediate_Recall_Pressed_Button, Immediate_Recall_Accuracy, Immediate_Recall_RT) %>%
          rename(id=ID,item=Presented_Person,first.button=Immediate_Recall_Pressed_Button,first.retrieval=Immediate_Recall_Accuracy,first.RT=Immediate_Recall_RT) %>% 
          mutate(first.rating=4)
        print(paste("reading: ", as.character(file)))
     
      } #parse encoding
      
      
      #get retrieval
      if (grepl("-retrieval_run",file, fixed = TRUE) == TRUE){
        full_file <- paste0(fullpath, file)
        ret <- read.delim(full_file)
        ret_behav <- ret %>% select(ID, Presented_Person, Assoc_Retrieval_Pressed_Button, Assoc_Retrieval_Accuracy, Assoc_Retrieval_RT, Retrieval_consc_Rating_Pressed_Button) %>%
          rename(id=ID,item=Presented_Person,second.button = Assoc_Retrieval_Pressed_Button, second.retrieval = Assoc_Retrieval_Accuracy, second.RT=Assoc_Retrieval_RT,second.rating=Retrieval_consc_Rating_Pressed_Button)
        print(paste("reading: ", as.character(file)))
        
      } 
      
      #get delayed ret
      if (grepl("delayed_retrieval_run",file, fixed = TRUE) == TRUE){
        full_file <- paste0(fullpath, file)
        delret <- read.delim(full_file)
        delret_behav <- delret %>% select(ID, Presented_Person, Assoc_Retrieval_Pressed_Button, Assoc_Retrieval_Accuracy, Assoc_Retrieval_RT, Retrieval_consc_Rating_Pressed_Button,
                                          ForcedChoice_Pressed_Button, ForcedChoice_Accuracy, ForcedChoice_RT, Recognition_consc_Rating_Pressed_Button) %>%
          rename(id=ID,item=Presented_Person,third.button=Assoc_Retrieval_Pressed_Button, third.retrieval=Assoc_Retrieval_Accuracy, third.RT = Assoc_Retrieval_RT,
                 third.rating = Retrieval_consc_Rating_Pressed_Button,fourth.button = ForcedChoice_Pressed_Button, fourth.retrieval = ForcedChoice_Accuracy,
                 fourth.RT = ForcedChoice_RT, fourth.rating = Recognition_consc_Rating_Pressed_Button) %>%
          mutate(TwoCorrectAnswerIfValueFour = recogCorrect$V1,
                 fourth.retrieval = case_when(
                   (fourth.button == 2 & TwoCorrectAnswerIfValueFour == 4) | (fourth.button == 3 & TwoCorrectAnswerIfValueFour == 6) ~ "correct",
                   (fourth.button == 3 & TwoCorrectAnswerIfValueFour == 4) | (fourth.button == 2 & TwoCorrectAnswerIfValueFour == 6) ~ "incorrect",
                   TRUE ~ "miss"
                 )) %>%
          select(-TwoCorrectAnswerIfValueFour)
        print(paste("reading: ", as.character(file)))
        
      }
    }
    df_behav <- right_join(enc_behav, ret_behav, by=c("id", "item"))
    df_behav <- right_join(df_behav, delret_behav, by=c("id", "item"))
    df_detailed <- rbind(df_detailed, df_behav)
  }
  
  if (sequence == 'hipp/') {
    df_detailed$group <- "hipp"
    df_detailed <- df_detailed %>% group_by(id) %>% mutate(rownums = row_number()) %>% mutate(upperhalf = if_else(rownums < 49, "upper","lower")) %>% select(-rownums)
    df_detailed_hipp <- rbind(df_detailed_hipp, df_detailed)
  } else if (sequence == 'wb/'){
    df_detailed$group <- "wb"
    df_detailed <- df_detailed %>% group_by(id) %>% mutate(rownums = row_number()) %>% mutate(upperhalf = if_else(rownums < 49, "upper","lower")) %>% select(-rownums)
    df_detailed_wb <- rbind(df_detailed_wb, df_detailed)
  }
  df_detailed_all <- rbind(df_detailed_all, df_detailed)
}  

write.csv(df_detailed_hipp,paste0(basepath,sequences[1],bids,output,"df_hipp_raw.csv"),row.names = FALSE)
write.csv(df_detailed_wb,paste0(basepath,sequences[2],bids,output,"df_wb_raw.csv"),row.names = FALSE)
          