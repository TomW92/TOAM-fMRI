library(pacman)
library(jsonlite)
library(tidyverse)
source('~/GitHub/TOAM/behavior/r_utils.R')

#check OS ----
if (.Platform$OS.type == 'windows'){
  filesep <- .Platform$file.sep
  volume <- 'Z:'  #or the letter your windows mounted the research storage unter
} else if (Sys.info()["sysname"] == "Darwin") { #check if mac/OSX
  volume <- '/Volumes/psy_wfg_henke'
  filesep <- .Platform$file.sep
} else if (Sys.info()["sysname"] == "Linux") {
  volume <- 'storage/research/psy_wf_henke'
  filesep <- .Platform$file.sep
}
sequences <- c("hipp", "wb")
basepath <- paste0(volume, "/s2019_twillems_fMRI_silent_engram/data/fMRI/")

for (sequence in sequences){
  dfMaxTrials <- findLateResponses(sequence)
  dfCorrectedClean <- cleanUpBehavData(dfMaxTrials)

  #for "Recog-Trial-by-exclusion" check
  dfCorrectedClean <- dfCorrectedClean %>% select(-rowNum)
  dfMaxTrials <- dfMaxTrials %>%  select(-rowNum)

  write.csv(dfCorrectedClean,paste0(basepath,sequence,"/bids/derivatives/behav/df_",sequence,"_corrected.csv"),row.names=FALSE)
  write.csv(dfMaxTrials,paste0(basepath,sequence,"/bids/derivatives/behav/df_",sequence,"_maxTrials.csv"),row.names=FALSE)
}