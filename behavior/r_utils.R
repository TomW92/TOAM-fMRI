library(pacman)
library(jsonlite)
library(tidyverse)

findLateResponses <- function(sequence){
  wideDataFramePath <- paste0(basepath, sequence, "/bids/derivatives/behav/df_", sequence, "_raw.csv")
  df <- read.csv(wideDataFramePath)
  df$rowNum <- as.integer(rownames(df))
  dfZero <- df  %>%  filter_all(any_vars(. %in% 0), .preserve = TRUE) # get df entries with misses
  
  if (nrow(dfZero) > 0){ #if there are no misses, move on
    print(paste("looking for potential misses, in this number of Trials:",as.character(nrow(dfZero))))
    for (row in seq_len(nrow(dfZero))){
      subject <- dfZero[row,]$id
      fullpath <- paste0(basepath, sequence, "/raw/sub-", subject, "/behav/")
      logFilePath <- paste0(fullpath, list.files(path=fullpath, pattern="-forgetting_experiment"))
      logFile <- read.table(logFilePath, sep = "\t", header=TRUE, skip=3, skipNul = TRUE, fill=TRUE)
      delLogFilePath <- paste0(fullpath, list.files(path=fullpath, pattern="-second_retrieval"))
      delLogFile <- read.table(delLogFilePath, sep = "\t", header=TRUE, skip=3, skipNul = TRUE, fill=TRUE)
      
      currentTrial <- dfZero[row,]
      trialName <- currentTrial$item
      for (colName in names(currentTrial[,currentTrial==0])) {
        if (colName == "first.button"){
          #encoding  #below: check where this trial is in the original logfile which also contains the late response misses
          rowNum <- which(logFile$Code %in% paste0("immediate_recall_", trialName))
          missingTrialTime <- logFile[rowNum,]$Time
          logFileCode <- "start" #arbitrary value for logfilecode
          while (grepl("Face_Baseline",logFileCode)==FALSE){ #as long es we're not at the next trial
            rowNum <- rowNum + 1
            logFileCode <- logFile[rowNum,]$Code
            if ((logFile[rowNum,]$Event.Type == "Response") & (logFile[rowNum,]$Code %in% c("2","3"))){ #if you find a response 
              ResponseTrial <- logFile[rowNum,]
              missResponse <- ResponseTrial$Code
              missRT <- (ResponseTrial$Time-missingTrialTime)/10
              df[currentTrial$rowNum,]$first.button <- missResponse
              df[currentTrial$rowNum,]$first.RT <- missRT
              if ((substr(trialName,2,2) == "a" & missResponse == "2") | (substr(trialName,2,2) == "o" & missResponse == "3")) {
                df[currentTrial$rowNum,]$first.retrieval <- "correct"
              } else {
                df[currentTrial$rowNum,]$first.retrieval <- "incorrect"
              }
            }  # if you find a response
          } # while next trial is not reached
          
        } else if (colName == "second.button"){
          #below: check where this trial is in the original logfile which also contains the late reponse misses
          rowNum <- as.integer(row.names(logFile[startsWith(logFile$Code, trialName),]))
          missingTrialTime <- logFile[rowNum,]$Time
          logFileCode <- "start" #arbitrary value for logfilecode
          while (grepl("conscious_rating_event",logFileCode)==FALSE){ #as long es we're not at the next trial
            rowNum <- rowNum + 1
            logFileCode <- logFile[rowNum,]$Code
            if ((logFile[rowNum,]$Event.Type == "Response") & (logFile[rowNum,]$Code %in% c("2","3"))){ #if you find a response that is one of the options
              ResponseTrial <- logFile[rowNum,]
              missResponse <- ResponseTrial$Code
              missRT <- (ResponseTrial$Time-missingTrialTime)/10
              df[currentTrial$rowNum,]$second.button <- missResponse
              df[currentTrial$rowNum,]$second.RT <- missRT
              if ((substr(trialName,2,2) == "a" & missResponse == "2") | (substr(trialName,2,2) == "o" & missResponse == "3")) {
                df[currentTrial$rowNum,]$second.retrieval <- "correct"
              } else {
                df[currentTrial$rowNum,]$second.retrieval <- "incorrect"
              }
            } # if you find a response
          } # while next trial
          #first retrieval
          
        } else if (colName == "second.rating"){
          #first rating
          #below: check where this trial is in the original logfile which also contains the late reponse misses
          rowNum <- as.integer(row.names(logFile[startsWith(logFile$Code, trialName),]))
          logFileCode <- "start"
          while (grepl("conscious_rating_event",logFileCode)==FALSE){ #lets look for the next rating-event,but do nothing
            rowNum <- rowNum + 1
            logFileCode <- logFile[rowNum,]$Code
          } 
          while (logFileCode != "face recognition"){ #as long as next trial doesnt start
            rowNum <- rowNum + 1
            logFileCode <- logFile[rowNum,]$Code
            if (logFile[rowNum,]$Event.Type == "Response"){
              ResponseTrial <- logFile[rowNum,]
              missResponse <- ResponseTrial$Code
              df[currentTrial$rowNum,]$second.rating <- missResponse
            }    # put in right response
          }
          
        } else if (colName == "third.button"){
          #second retrieval
          #below: check where this trial is in the original logfile which also contains the late reponse misses
          rowNum <- as.integer(row.names(delLogFile[startsWith(delLogFile$Code, paste0("Retrieval of: ", trialName)),]))
          missingTrialTime <- delLogFile[rowNum,]$Time
          logFileCode <- "start" #arbitrary value for logfilecode
          while (grepl("conscious_rating_event",logFileCode)==FALSE){ #as long es we're not at the next trial
            rowNum <- rowNum + 1
            logFileCode <- delLogFile[rowNum,]$Code
            if ((delLogFile[rowNum,]$Event.Type == "Response") & (delLogFile[rowNum,]$Code %in% c("2","3"))){ #if you find a response that is one of the options
              ResponseTrial <- delLogFile[rowNum,]
              missResponse <- ResponseTrial$Code
              missRT <- (ResponseTrial$Time-missingTrialTime)/10
              df[currentTrial$rowNum,]$third.button <- missResponse
              df[currentTrial$rowNum,]$third.RT <- missRT
              if ((substr(trialName,2,2) == "a" & missResponse == "2") | (substr(trialName,2,2) == "o" & missResponse == "3")) {
                df[currentTrial$rowNum,]$third.retrieval <- "correct"
              } else {
                df[currentTrial$rowNum,]$third.retrieval <- "incorrect"
              }
            } # if you find a response
          } # while next trial is not reached
          
        } else if (colName == "third.rating"){
          #second rating
          #below: check where this trial is in the original delLogFile which also contains the late reponse misses
          rowNum <- as.integer(row.names(delLogFile[startsWith(delLogFile$Code, paste0("Retrieval of: ", trialName)),]))
          logFileCode <- "start"
          while (grepl("conscious_rating_event",logFileCode)==FALSE){ #lets look for the next rating-event,but do nothing
            rowNum <- rowNum + 1
            logFileCode <- delLogFile[rowNum,]$Code
          } 
          while (!(startsWith(logFileCode,"recog: "))){ #as long as next trial doesnt start
            rowNum <- rowNum + 1
            logFileCode <- delLogFile[rowNum,]$Code
            if (delLogFile[rowNum,]$Event.Type == "Response"){
              ResponseTrial <- delLogFile[rowNum,]
              missResponse <- ResponseTrial$Code
              df[currentTrial$rowNum,]$third.rating <- missResponse
            }    # put in right response
          } # while loop
          
        } else if (colName == "fourth.button"){
          #recog
          rowNum <- as.integer(row.names(delLogFile[startsWith(delLogFile$Code, paste0("recog: ", trialName)),]))
          missingTrialTime <- delLogFile[rowNum,]$Time
          logFileCode <- "start" #arbitrary value for logfilecode
          while (grepl("conscious_rating_event",logFileCode)==FALSE){ #as long es we're not at the next trial
            rowNum <- rowNum + 1
            logFileCode <- delLogFile[rowNum,]$Code
            if ((delLogFile[rowNum,]$Event.Type == "Response") & (delLogFile[rowNum,]$Code %in% c("2","3"))){ #if you find a response that is one of the options
              ResponseTrial <- delLogFile[rowNum,]
              missResponse <- ResponseTrial$Code
              missRT <- (ResponseTrial$Time-missingTrialTime)/10
              df[currentTrial$rowNum,]$fourth.button <- missResponse
              df[currentTrial$rowNum,]$fourth.RT <- missRT
              if ((substr(trialName,2,2) == "a" & missResponse == "2") | (substr(trialName,2,2) == "o" & missResponse == "3")) {
                df[currentTrial$rowNum,]$fourth.retrieval <- "correct"
              } else {
                df[currentTrial$rowNum,]$fourth.retrieval <- "incorrect"
              }
            } # if you find a response
          } 
          
        } else if (colName == "fourth.rating"){
          #recog rating
          #below: check where this trial is in the original delLogFile which also contains the late reponse misses
          rowNum <- as.integer(row.names(delLogFile[startsWith(delLogFile$Code, paste0("recog: ", trialName)),]))
          logFileCode <- "start"
          while (grepl("conscious_rating_event",logFileCode)==FALSE){ #lets look for the next rating-event,but do nothing
            rowNum <- rowNum + 1
            logFileCode <- delLogFile[rowNum,]$Code
          } 
          while (logFileCode != "face recognition"){ #as long as next trial doesn't start
            rowNum <- rowNum + 1
            logFileCode <- delLogFile[rowNum,]$Code
            if (delLogFile[rowNum,]$Event.Type == "Response"){
              ResponseTrial <- delLogFile[rowNum,]
              missResponse <- ResponseTrial$Code
              df[currentTrial$rowNum,]$fourth.rating <- missResponse
            }    # put in right response
          } # while loop
        }
      }
    }
  }
  stillNotCleanRows <- df %>% filter_all(any_vars(. %in% 0), .preserve = TRUE) %>% nrow() # count misses after correction
  print(paste("corrected Number of trials:", nrow(dfZero)-stillNotCleanRows))
  print(paste("This number of rows could not be corrected with the help of the original log-file: ", as.character(stillNotCleanRows)))
  
  
  return(df)
}

cleanUpBehavData <- function(df){
  # also delete unattended or not understood trials.
  rowsBefore <- nrow(df)
  df <- df %>% filter(fourth.RT > 100) %>%
    filter(first.retrieval == "correct")
  rowsAfter <- nrow(df)
  print(paste0("deleted because of not paying attention during imm. recall or too short reaction time in recog:",rowsBefore-rowsAfter))
  
  # delete false Memories.
  rowsBefore <- nrow(df)
  df <- df %>% filter(
    !(second.retrieval == "incorrect" & second.rating == "4") &
      !(third.retrieval == "incorrect" & third.rating == "4") &
      !(fourth.retrieval == "incorrect" & fourth.rating == "4"))
  rowsAfter <- nrow(df)
  print(paste0("deleted false memory trials (high confidence and wrong answers):",rowsBefore-rowsAfter))
  
  #only allow valid answers
  rowsBefore <- nrow(df)
  df <- df %>% 
    filter(if_any(first.button, ~  (.x %in% c("2", "3")))) %>%
    filter(if_any(second.button, ~ (.x %in% c("2", "3")))) %>%
    filter(if_any(third.button, ~  (.x %in% c("2", "3")))) %>%
    filter(if_any(fourth.button, ~ (.x %in% c("2", "3")))) %>%
    filter(if_any(first.rating, ~  (.x %in% c(2, 3,4)))) %>%
    filter(if_any(second.rating, ~ (.x %in% c(2, 3,4)))) %>%
    filter(if_any(third.rating, ~  (.x %in% c(2, 3,4)))) %>%
    filter(if_any(fourth.rating, ~ (.x %in% c(2, 3,4))))
  rowsAfter <- nrow(df)
  print(paste0("deleted because of out of bounds answer (pressed button 1 instead of 2 or 3): ",rowsBefore-rowsAfter))

  return(df)
}